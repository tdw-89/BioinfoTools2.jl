module Paralogs

using DataFrames
using Graphs
using SparseArrays
using StructArrays

using ..Reference

function _match_relation(colname)
    key = lowercase(String(colname))
    return key == "dn" ? :dN :
           key == "ds" ? :dS :
           key == "id_subject_query" ? :id_subject_query :
           key == "id_query_subject" ? :id_query_subject : nothing
end

# Marker distinguishing the raw-fields constructor (below) from the public
# `GeneFamily(genome, pairs)` one; kept unexported so ordinary callers can't
# reach it and reintroduce arbitrary construction.
struct _RawFields end

struct GeneFamily
    # The gene coordinates keyed by scaffold name/id
    intervals::Dict{String,StructArray{Reference.IntervalSimple}}

    # The linear index ranges (inclusive) spanned by each scaffold
    scaffold_ranges::Dict{String,UnitRange{Int}}

    # Reverse lookup from gene ID to its linear index, for string-based `getindex`
    id_to_index::Dict{String,Int}

    topology::SparseMatrixCSC{Bool,UInt32}

    dN::Union{Nothing,SparseMatrixCSC{Float64,UInt32}}
    dS::Union{Nothing,SparseMatrixCSC{Float64,UInt32}}
    id_subject_query::Union{Nothing,SparseMatrixCSC{Float64,UInt32}}
    id_query_subject::Union{Nothing,SparseMatrixCSC{Float64,UInt32}}

    """
        GeneFamily(genome::Reference.Genome, pairs::DataFrame)

    Build a `GeneFamily` from a reference `genome` and a paralog-pair table.

    `pairs` has 2-6 columns; the first two are required and taken by position,
    the rest optional and matched to the relation matrices *by name*
    (case-insensitive):

    1. Query ID (`String`)
    2. Subject ID (`String`)
    3. `dN`               → pairwise dN
    4. `dS`               → pairwise dS
    5. `id_subject_query` → % identity subject → query
    6. `id_query_subject` → % identity query → subject

    With only the two ID columns just `topology` is populated; every recognised
    value column adds its matrix and the rest stay `nothing`.

    Each gene is placed by looking its ID up in `genome`, inheriting the
    feature's 1-based bounds and full 64-bit metadata code (strand, SO term,
    metadata index). Genes are grouped by scaffold and handed contiguous linear
    indices — `scaffold_ranges[name]` is that block. IDs absent from `genome`,
    and any pair row referencing one, are skipped.

    Every matrix is `n_genes x n_genes` over those indices, oriented so the axis
    you look up by is the (column-major) fast axis:

    - `topology`         — symmetric `Bool` adjacency; a pair sets both `[a, b]`
      and `[b, a]`, so the relation is undirected.
    - `dN`, `dS`         — symmetric by definition, so stored **once** in the
      upper triangle: the value for genes `a`, `b` is at `[min(a, b), max(a, b)]`.
    - `id_subject_query` — %ID subject → query with **subjects on the columns**
      (`[query, subject]`), so `M[:, s]` gathers subject `s`'s scores in O(nnz).
    - `id_query_subject` — %ID query → subject with **queries on the columns**
      (`[subject, query]`), so `M[:, q]` gathers query `q`'s scores in O(nnz).

    Repeated or reciprocal entries for a cell collapse to one value.

    This is the only public constructor, so `intervals`, `scaffold_ranges` and
    `id_to_index` are always concordant by construction. Query/subject ID
    columns must hold `AbstractString`s and every matched value column must
    hold `Real`s; anything else raises an informative `ArgumentError`.

    See `getindex` for per-gene (`gf[i]`/`gf["id"]`) and sub-family
    (`gf[indices]`/`gf[ids]`) indexing.
    """
    function GeneFamily(genome::Reference.Genome, pairs::DataFrame)
        ncol(pairs) >= 2 ||
            throw(ArgumentError("`pairs` needs at least a query and a subject column"))

        eltype(pairs[:, 1]) <: AbstractString || throw(
            ArgumentError(
                "Query ID column (column 1, \"$(names(pairs)[1])\") must contain strings, got element type $(eltype(pairs[:, 1]))",
            ),
        )
        eltype(pairs[:, 2]) <: AbstractString || throw(
            ArgumentError(
                "Subject ID column (column 2, \"$(names(pairs)[2])\") must contain strings, got element type $(eltype(pairs[:, 2]))",
            ),
        )

        query_ids = string.(pairs[:, 1])
        subject_ids = string.(pairs[:, 2])

        # Resolve every ID to its feature.
        records = genome[unique(vcat(query_ids, subject_ids))]

        # Group the found features by scaffold, discarding repeated IDs.
        scaffold_intervals = Dict{String,Vector{Reference.IntervalSimple}}()
        scaffold_id_lists = Dict{String,Vector{String}}()
        seen = Set{String}()
        for record in records
            record.id in seen && continue
            push!(seen, record.id)
            push!(
                get!(scaffold_intervals, record.chromosome, Reference.IntervalSimple[]),
                Reference.IntervalSimple(record.start_pos, record.end_pos, record.code),
            )
            push!(get!(scaffold_id_lists, record.chromosome, String[]), record.id)
        end

        # Lay the scaffolds out end to end, ordering the genes within each by
        # position, and remember where every gene lands in the linear indexing.
        intervals = Dict{String,StructArray{Reference.IntervalSimple}}()
        scaffold_ranges = Dict{String,UnitRange{Int}}()
        id_to_index = Dict{String,Int}()
        cursor = 1
        for name in sort!(collect(keys(scaffold_intervals)))
            ivs = scaffold_intervals[name]
            ids = scaffold_id_lists[name]
            order = collect(1:length(ivs))
            sort!(order; by = k -> (ivs[k].start_pos, ivs[k].end_pos, ids[k]))
            ivs, ids = ivs[order], ids[order]

            n = length(ivs)
            scaffold_ranges[name] = cursor:(cursor+n-1)
            intervals[name] = StructArray(ivs)
            for (offset, id) in enumerate(ids)
                id_to_index[id] = cursor + offset - 1
            end
            cursor += n
        end
        n_genes = cursor - 1

        # Match the optional value columns to their target matrices by name.
        value_columns = Dict{Symbol,Int}()
        for c = 3:ncol(pairs)
            field = _match_relation(names(pairs)[c])
            isnothing(field) || (value_columns[field] = c)
        end

        for (field, c) in value_columns
            eltype(pairs[:, c]) <: Real || throw(
                ArgumentError(
                    "Column \"$(names(pairs)[c])\" (matched to `$field`) must contain numeric values, got element type $(eltype(pairs[:, c]))",
                ),
            )
        end

        has(field) = haskey(value_columns, field)

        # Per-matrix coordinates/values, each oriented so the axis you look up
        # by is the (fast) column axis; dN/dS are symmetric and fold into a
        # single upper triangle. See the docstring for the exact conventions.
        topo_i, topo_j = UInt32[], UInt32[]
        tri_i, tri_j = UInt32[], UInt32[]          # shared upper-triangle coords
        dN_v, dS_v = Float64[], Float64[]
        sq_i, sq_j, sq_v = UInt32[], UInt32[], Float64[]
        qs_i, qs_j, qs_v = UInt32[], UInt32[], Float64[]

        want_tri = has(:dN) || has(:dS)
        for r = 1:nrow(pairs)
            q, s = query_ids[r], subject_ids[r]
            (haskey(id_to_index, q) && haskey(id_to_index, s)) || continue
            iq, is = id_to_index[q], id_to_index[s]

            push!(topo_i, iq)
            push!(topo_j, is)

            if want_tri
                lo, hi = minmax(iq, is)
                push!(tri_i, lo)
                push!(tri_j, hi)
                has(:dN) && push!(dN_v, Float64(pairs[r, value_columns[:dN]]))
                has(:dS) && push!(dS_v, Float64(pairs[r, value_columns[:dS]]))
            end

            if has(:id_subject_query)
                push!(sq_i, iq)               # subjects on the column axis
                push!(sq_j, is)
                push!(sq_v, Float64(pairs[r, value_columns[:id_subject_query]]))
            end
            if has(:id_query_subject)
                push!(qs_i, is)               # queries on the column axis
                push!(qs_j, iq)
                push!(qs_v, Float64(pairs[r, value_columns[:id_query_subject]]))
            end
        end

        # `max` collapses repeated/reciprocal writes to one value rather than
        # summing them; topology is symmetrised (both directions) and ORed.
        topology = sparse(
            vcat(topo_i, topo_j),
            vcat(topo_j, topo_i),
            trues(2 * length(topo_i)),
            n_genes,
            n_genes,
            |,
        )
        dN = has(:dN) ? sparse(tri_i, tri_j, dN_v, n_genes, n_genes, max) : nothing
        dS = has(:dS) ? sparse(tri_i, tri_j, dS_v, n_genes, n_genes, max) : nothing
        id_subject_query =
            has(:id_subject_query) ? sparse(sq_i, sq_j, sq_v, n_genes, n_genes, max) :
            nothing
        id_query_subject =
            has(:id_query_subject) ? sparse(qs_i, qs_j, qs_v, n_genes, n_genes, max) :
            nothing

        return new(
            intervals,
            scaffold_ranges,
            id_to_index,
            topology,
            dN,
            dS,
            id_subject_query,
            id_query_subject,
        )
    end

    # Internal-only: build directly from already-resolved fields, used by the
    # sub-family `getindex` methods below. The `_RawFields` marker keeps this
    # unreachable from `GeneFamily(...)` calls outside this module.
    function GeneFamily(
        ::_RawFields,
        intervals::Dict{String,StructArray{Reference.IntervalSimple}},
        scaffold_ranges::Dict{String,UnitRange{Int}},
        id_to_index::Dict{String,Int},
        topology::SparseMatrixCSC{Bool,UInt32},
        dN,
        dS,
        id_subject_query,
        id_query_subject,
    )
        return new(
            intervals,
            scaffold_ranges,
            id_to_index,
            topology,
            dN,
            dS,
            id_subject_query,
            id_query_subject,
        )
    end
end

_n_genes(gf::GeneFamily) = size(gf.topology, 1)

function Base.show(io::IO, gf::GeneFamily)
    n = _n_genes(gf)
    nscaff = length(gf.scaffold_ranges)
    npairs = nnz(gf.topology) ÷ 2
    relations = [
        String(f) for f in (:dN, :dS, :id_subject_query, :id_query_subject) if
        !isnothing(getfield(gf, f))
    ]
    rel_str = isempty(relations) ? "none" : join(relations, ", ")
    print(
        io,
        "GeneFamily($n gene$(n == 1 ? "" : "s"), $nscaff scaffold$(nscaff == 1 ? "" : "s"), $npairs pair$(npairs == 1 ? "" : "s"); relations: $rel_str)",
    )
end

"""
    gf[i::Integer]
    gf[id::AbstractString]

Return column `i` (or the column for gene `id`) from every relation matrix
present on `gf`, as a `Dict{Symbol, Vector}` keyed by `:topology` and whichever
of `:dN`, `:dS`, `:id_subject_query`, `:id_query_subject` are non-`nothing`.

`dN`/`dS` store only their upper triangle (see the `GeneFamily` docstring), so
column `i` alone only surfaces partners with a smaller index; combine with row
`i` for the complete symmetric relation.
"""
function Base.getindex(gf::GeneFamily, i::Integer)
    n = _n_genes(gf)
    1 <= i <= n || throw(ArgumentError("gene index $i out of bounds (1:$n)"))
    cols = Dict{Symbol,Vector}(:topology => Vector(gf.topology[:, i]))
    for field in (:dN, :dS, :id_subject_query, :id_query_subject)
        m = getfield(gf, field)
        isnothing(m) || (cols[field] = Vector(m[:, i]))
    end
    return cols
end

function Base.getindex(gf::GeneFamily, id::AbstractString)
    haskey(gf.id_to_index, id) || throw(ArgumentError("Unknown gene ID: \"$id\""))
    return gf[gf.id_to_index[id]]
end

"""
    gf[indices::AbstractVector{<:Integer}]
    gf[ids::AbstractVector{<:AbstractString}]

Return a new `GeneFamily` restricted to the given genes (a "sub-family"): every
relation matrix present on `gf` is sliced to just those genes, and
`intervals`/`scaffold_ranges` are rebuilt to stay concordant with the new,
contiguous linear indexing.
"""
function Base.getindex(gf::GeneFamily, indices::AbstractVector{<:Integer})
    isempty(indices) && throw(ArgumentError("`indices` must be non-empty"))
    n = _n_genes(gf)
    all(i -> 1 <= i <= n, indices) ||
        throw(ArgumentError("gene index out of bounds (1:$n)"))

    order = sort!(unique(collect(Int, indices)))
    order_set = Set(order)
    index_to_id = Dict(v => k for (k, v) in gf.id_to_index)

    new_intervals = Dict{String,StructArray{Reference.IntervalSimple}}()
    new_ranges = Dict{String,UnitRange{Int}}()
    new_id_to_index = Dict{String,Int}()
    cursor = 1
    for name in sort!(collect(keys(gf.scaffold_ranges)))
        rng = gf.scaffold_ranges[name]
        kept = filter(gi -> gi in order_set, rng)
        isempty(kept) && continue

        ivs = [gf.intervals[name][gi-first(rng)+1] for gi in kept]
        m = length(kept)
        new_intervals[name] = StructArray(ivs)
        new_ranges[name] = cursor:(cursor+m-1)
        for (offset, gi) in enumerate(kept)
            new_id_to_index[index_to_id[gi]] = cursor + offset - 1
        end
        cursor += m
    end

    relation(field) = (m = getfield(gf, field); isnothing(m) ? nothing : m[order, order])

    return GeneFamily(
        _RawFields(),
        new_intervals,
        new_ranges,
        new_id_to_index,
        gf.topology[order, order],
        relation(:dN),
        relation(:dS),
        relation(:id_subject_query),
        relation(:id_query_subject),
    )
end

function Base.getindex(gf::GeneFamily, ids::AbstractVector{<:AbstractString})
    isempty(ids) && throw(ArgumentError("`ids` must be non-empty"))
    indices = map(ids) do id
        haskey(gf.id_to_index, id) || throw(ArgumentError("Unknown gene ID: \"$id\""))
        gf.id_to_index[id]
    end
    return gf[indices]
end

function rbh_ds(paralog_df::DataFrame)
    @assert typeof(paralog_df[1, 1]) <: AbstractString
    @assert typeof(paralog_df[1, 2]) <: AbstractString
    @assert typeof(paralog_df[1, 3]) <: AbstractFloat

    unique_ids = unique(vcat(paralog_df[:, 1], paralog_df[:, 2]))
    ids_to_ind_dict = Dict(unique_ids[i] => i for i in eachindex(unique_ids))
    ind_to_ids_dict = Dict(i => unique_ids[i] for i in eachindex(unique_ids))

    # Create a matrix of zeros
    rbh_matrix = zeros(Float64, length(unique_ids), length(unique_ids))
    orig_mat = zeros(Float64, length(unique_ids), length(unique_ids))

    # Fill in the matrix, treating gene 'i' as the query and gene 'j' as the subject
    for row in eachrow(paralog_df)

        i = ids_to_ind_dict[row[1]]
        j = ids_to_ind_dict[row[2]]

        orig_mat[i, j] = row[3]
        orig_mat[j, i] = row[3]
        rbh_matrix[i, j] = row[3]
        rbh_matrix[j, i] = row[3]
    end

    rbh_gene, rbh_paralog = String[], String[]
    matched_inds = Int[]
    score_i_origs, score_j_origs, min_scores, mean_scores =
        Float64[], Float64[], Float64[], Float64[]
    for i = 1:size(rbh_matrix)[1]

        _, min_i = findmin(rbh_matrix[i, :])[1:2]
        _, min_j = findmin(rbh_matrix[:, min_i])[1:2]
        score_i_orig = orig_mat[i, min_i]
        score_j_orig = orig_mat[min_i, i]
        min_score = min(score_i_orig, score_j_orig)
        mean_score = (score_i_orig + score_j_orig) / 2

        if i == min_j && min_i ∉ matched_inds && min_j ∉ matched_inds

            id_i = ind_to_ids_dict[min_i]
            id_j = ind_to_ids_dict[min_j]

            push!(rbh_gene, id_i)
            push!(rbh_paralog, id_j)
            push!(matched_inds, min_i)
            push!(matched_inds, min_j)
            push!(score_i_origs, score_i_orig)
            push!(score_j_origs, score_j_orig)
            push!(min_scores, min_score)
            push!(mean_scores, mean_score)
        end
    end

    return DataFrame(
        "GeneID" => rbh_gene,
        "ParalogID" => rbh_paralog,
        "ds" => score_i_origs,
        "min_ds" => min_scores,
    )
end

"""
    rbh(paralog_df; scoring="max")

Identify reciprocal best hits (RBH) between paralogs based on similarity scores.

# Arguments
- `paralog_df::DataFrame`: DataFrame with at least 4 columns:
  1. GeneID (String)
  2. ParalogID (String)
  3. Percent identity from gene to paralog (Float), or dS value if using `scoring="ds"`
  4. Percent identity from paralog to gene (Float), or ignored if using `scoring="ds"`
- `scoring::String="max"`: Scoring method for determining best hits
  - `"max"` or `"maximum"`: Use maximum of bidirectional scores
  - `"mean"`, `"avg"`, or `"average"`: Use mean of bidirectional scores
  - `"double_max"`: Use original bidirectional scores
  - `"ds"`: Use dS values (column 3 should contain dS values)

# Returns
- `DataFrame`: Contains reciprocal best hit pairs with columns:
  - `GeneID`: First gene in RBH pair
  - `ParalogID`: Second gene in RBH pair
  - `perc_1`: Original percent identity (gene → paralog)
  - `perc_2`: Original percent identity (paralog → gene)
  - `max_perc`: Maximum of the two scores
  - `mean_perc`: Mean of the two scores

# Examples
```julia
df = DataFrame(
    GeneID = ["A", "B", "C"],
    ParalogID = ["B", "A", "D"],
    Perc1 = [95.0, 94.0, 80.0],
    Perc2 = [94.0, 95.0, 85.0]
)
rbh_pairs = rbh(df; scoring="max")
```
"""
function rbh(paralog_df::DataFrame; scoring::String = "max")

    scoring = lowercase(scoring)

    if !(scoring in ["max", "maximum", "double_max", "mean", "average", "avg", "ds"])

        error(
            "Invalid scoring method. Must be 'ds', 'max', 'maximum', 'double_max', 'mean', 'avg', or 'average'.",
        )
    end

    scoring =
        scoring in ["max", "maximum"] ? "max" :
        scoring in ["mean", "avg", "average"] ? "mean" : scoring

    if scoring == "ds"
        return rbh_ds(paralog_df)
    end

    @assert typeof(paralog_df[1, 1]) <: AbstractString
    @assert typeof(paralog_df[1, 2]) <: AbstractString
    @assert typeof(paralog_df[1, 3]) <: AbstractFloat
    @assert typeof(paralog_df[1, 4]) <: AbstractFloat

    unique_ids = unique(vcat(paralog_df[:, 1], paralog_df[:, 2]))
    ids_to_ind_dict = Dict(unique_ids[i] => i for i in eachindex(unique_ids))
    ind_to_ids_dict = Dict(i => unique_ids[i] for i in eachindex(unique_ids))

    # Create a matrix of zeros
    rbh_matrix = zeros(Float64, length(unique_ids), length(unique_ids))
    orig_mat = zeros(Float64, length(unique_ids), length(unique_ids))

    # Fill in the matrix, treating gene 'i' as the query and gene 'j' as the subject
    for row in eachrow(paralog_df)

        i = ids_to_ind_dict[row[1]]
        j = ids_to_ind_dict[row[2]]

        orig_mat[i, j] = row[3]
        orig_mat[j, i] = row[4]

        if scoring == "max"

            max_perc_temp = max(row[3], row[4])
            rbh_matrix[i, j] = max_perc_temp
            rbh_matrix[j, i] = max_perc_temp
        elseif scoring == "mean"

            mean_perc_temp = (row[3] + row[4]) / 2
            rbh_matrix[i, j] = mean_perc_temp
            rbh_matrix[j, i] = mean_perc_temp
        else

            rbh_matrix[i, j] = orig_mat[i, j]
            rbh_matrix[j, i] = orig_mat[j, i]
        end
    end

    rbh_gene, rbh_paralog = String[], String[]
    matched_inds = Int[]
    perc_i_origs, perc_j_origs, max_percs, mean_percs =
        Float64[], Float64[], Float64[], Float64[]

    for i = 1:size(rbh_matrix)[1]

        _, max_i = findmax(rbh_matrix[i, :])[1:2]
        _, max_j = findmax(rbh_matrix[:, max_i])[1:2]
        perc_i_orig = orig_mat[i, max_i]
        perc_j_orig = orig_mat[max_i, i]
        max_perc = max(perc_i_orig, perc_j_orig)
        mean_perc = (perc_i_orig + perc_j_orig) / 2

        if i == max_j && max_i ∉ matched_inds && max_j ∉ matched_inds

            id_i = ind_to_ids_dict[max_i]
            id_j = ind_to_ids_dict[max_j]

            push!(rbh_gene, id_i)
            push!(rbh_paralog, id_j)
            push!(matched_inds, max_i)
            push!(matched_inds, max_j)
            push!(perc_i_origs, perc_i_orig)
            push!(perc_j_origs, perc_j_orig)
            push!(max_percs, max_perc)
            push!(mean_percs, mean_perc)
        end
    end

    return DataFrame(
        "GeneID" => rbh_gene,
        "ParalogID" => rbh_paralog,
        "perc_1" => perc_i_origs,
        "perc_2" => perc_j_origs,
        "max_perc" => max_percs,
        "mean_perc" => mean_percs,
    )
end

export GeneFamily, rbh, rbh_ds

end

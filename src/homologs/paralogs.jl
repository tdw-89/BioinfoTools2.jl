"""
Homology *within* a genome: gene families of paralogous pairs, the relation
matrices over them, and reciprocal-best-hit detection.

See [`Orthologs`](@ref) for the between-genome counterpart.
"""
module Paralogs

using DataFrames
using Graphs
using SimpleWeightedGraphs
using SparseArrays
using StructArrays

# Nested one level deeper than the other modules, so `Reference` is three dots up.
using ...Reference

"""Relation matrix an optional `pairs` column name maps to, or `nothing`."""
function _match_relation(colname)
    key = lowercase(String(colname))
    key == "dn" && return :dN
    key == "ds" && return :dS
    key == "id_subject_query" && return :id_subject_query
    key == "id_query_subject" && return :id_query_subject
    return nothing
end

"""The optional relation fields of a [`GeneFamily`](@ref), in declaration order."""
const RELATIONS = (:dN, :dS, :id_subject_query, :id_query_subject)

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

    1. Query ID (`String`; mandatory)
    2. Subject ID (`String`; mandatory)
    3. `dN`               → pairwise dN (Float; optional)
    4. `dS`               → pairwise dS (Float; optional)
    5. `id_subject_query` → % identity subject → query (Float; optional)
    6. `id_query_subject` → % identity query → subject (Float; optional)

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

    Repeated or reciprocal entries for a cell collapse to one value. Rows
    carrying a `NaN` in any matched value column are dropped up front, with a
    warning: `NaN != NaN` would break the symmetry the `*_graph` views check.

    See `getindex` for per-gene (`gf[i]`/`gf["id"]`) and sub-family
    (`gf[indices]`/`gf[ids]`) indexing, and `topology_graph`/`dN_graph`/
    `dS_graph`/`id_subject_query_graph`/`id_query_subject_graph` for `Graphs`/
    `SimpleWeightedGraphs` views of the relation matrices.
    """
    function GeneFamily(genome::Reference.Genome, pairs::DataFrame)
        ncol(pairs) >= 2 ||
            throw(ArgumentError("`pairs` needs at least a query and a subject column"))

        for (column, role) in ((1, "Query"), (2, "Subject"))
            eltype(pairs[:, column]) <: AbstractString || throw(
                ArgumentError(
                    "$role ID column (column $column, \"$(names(pairs)[column])\") must contain strings, got element type $(eltype(pairs[:, column]))",
                ),
            )
        end

        # Match the optional value columns to their target matrices by name.
        value_columns = Dict{Symbol,Int}()
        for column = 3:ncol(pairs)
            field = _match_relation(names(pairs)[column])
            isnothing(field) || (value_columns[field] = column)
        end

        for (field, column) in value_columns
            eltype(pairs[:, column]) <: Real || throw(
                ArgumentError(
                    "Column \"$(names(pairs)[column])\" (matched to `$field`) must contain numeric values, got element type $(eltype(pairs[:, column]))",
                ),
            )
        end

        if !isempty(value_columns)
            keep = trues(nrow(pairs))
            for column in values(value_columns)
                keep .&= .!isnan.(pairs[:, column])
            end
            n_dropped = count(!, keep)
            if n_dropped > 0
                @warn "GeneFamily: dropped $n_dropped row$(n_dropped == 1 ? "" : "s") with a NaN relation value"
                pairs = pairs[keep, :]
            end
        end

        query_ids = string.(pairs[:, 1])
        subject_ids = string.(pairs[:, 2])
        # Materialised up front: indexing a DataFrame cell per row is both
        # type-unstable and far slower than indexing a `Vector{Float64}`.
        values_by_field = Dict{Symbol,Vector{Float64}}(
            field => Float64.(pairs[:, column]) for (field, column) in value_columns
        )

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
            scaffold_ivs = scaffold_intervals[name]
            ids = scaffold_id_lists[name]
            order = sortperm(
                eachindex(scaffold_ivs);
                by = k -> (scaffold_ivs[k].start_pos, scaffold_ivs[k].end_pos, ids[k]),
            )

            n_scaffold_genes = length(order)
            scaffold_ranges[name] = cursor:(cursor+n_scaffold_genes-1)
            intervals[name] = StructArray(scaffold_ivs[order])
            for (offset, k) in enumerate(order)
                id_to_index[ids[k]] = cursor + offset - 1
            end
            cursor += n_scaffold_genes
        end
        n_genes = cursor - 1

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
        for row = 1:nrow(pairs)
            query, subject = query_ids[row], subject_ids[row]
            (haskey(id_to_index, query) && haskey(id_to_index, subject)) || continue
            query_index, subject_index = id_to_index[query], id_to_index[subject]

            push!(topo_i, query_index)
            push!(topo_j, subject_index)

            if want_tri
                lo, hi = minmax(query_index, subject_index)
                push!(tri_i, lo)
                push!(tri_j, hi)
                has(:dN) && push!(dN_v, values_by_field[:dN][row])
                has(:dS) && push!(dS_v, values_by_field[:dS][row])
            end

            if has(:id_subject_query)
                push!(sq_i, query_index)          # subjects on the column axis
                push!(sq_j, subject_index)
                push!(sq_v, values_by_field[:id_subject_query][row])
            end
            if has(:id_query_subject)
                push!(qs_i, subject_index)        # queries on the column axis
                push!(qs_j, query_index)
                push!(qs_v, values_by_field[:id_query_subject][row])
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
        relation(rows, cols, vals) = sparse(rows, cols, vals, n_genes, n_genes, max)

        return new(
            intervals,
            scaffold_ranges,
            id_to_index,
            topology,
            has(:dN) ? relation(tri_i, tri_j, dN_v) : nothing,
            has(:dS) ? relation(tri_i, tri_j, dS_v) : nothing,
            has(:id_subject_query) ? relation(sq_i, sq_j, sq_v) : nothing,
            has(:id_query_subject) ? relation(qs_i, qs_j, qs_v) : nothing,
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

"""Number of genes a family spans."""
_n_genes(gf::GeneFamily) = size(gf.topology, 1)

"""Reverse of `gf.id_to_index`: linear gene index to gene ID."""
_index_to_id(gf::GeneFamily) = Dict(index => id for (id, index) in gf.id_to_index)

function Base.show(io::IO, gf::GeneFamily)
    n = _n_genes(gf)
    n_scaffolds = length(gf.scaffold_ranges)
    n_pairs = nnz(gf.topology) ÷ 2
    present = [String(field) for field in RELATIONS if !isnothing(getfield(gf, field))]
    print(
        io,
        "GeneFamily($n gene$(n == 1 ? "" : "s"), $n_scaffolds scaffold$(n_scaffolds == 1 ? "" : "s"), $n_pairs pair$(n_pairs == 1 ? "" : "s"); relations: $(isempty(present) ? "none" : join(present, ", ")))",
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
    columns = Dict{Symbol,Vector}(:topology => Vector(gf.topology[:, i]))
    for field in RELATIONS
        matrix = getfield(gf, field)
        isnothing(matrix) || (columns[field] = Vector(matrix[:, i]))
    end
    return columns
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
    index_to_id = _index_to_id(gf)

    new_intervals = Dict{String,StructArray{Reference.IntervalSimple}}()
    new_ranges = Dict{String,UnitRange{Int}}()
    new_id_to_index = Dict{String,Int}()
    cursor = 1
    for name in sort!(collect(keys(gf.scaffold_ranges)))
        range = gf.scaffold_ranges[name]
        kept = filter(in(order_set), range)
        isempty(kept) && continue

        new_intervals[name] =
            StructArray([gf.intervals[name][gene-first(range)+1] for gene in kept])
        new_ranges[name] = cursor:(cursor+length(kept)-1)
        for (offset, gene) in enumerate(kept)
            new_id_to_index[index_to_id[gene]] = cursor + offset - 1
        end
        cursor += length(kept)
    end

    function relation(field)
        matrix = getfield(gf, field)
        return isnothing(matrix) ? nothing : matrix[order, order]
    end

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

#= Graph views =#

"""
    topology_graph(gf::GeneFamily) -> SimpleGraph

Unweighted `Graphs.SimpleGraph` view of `gf.topology`: vertex `v` is gene `v`
(see `gf.id_to_index` for the name), and edge `u-v` marks a topology pair.
"""
topology_graph(gf::GeneFamily) = SimpleGraph(gf.topology)

_weighted_undirected_graph(::Nothing) = nothing
# dN/dS store only their upper triangle; SimpleWeightedGraph requires a fully
# symmetric adjacency matrix, so mirror it across the diagonal first.
_weighted_undirected_graph(m::SparseMatrixCSC) = SimpleWeightedGraph(m + permutedims(m))

"""
    dN_graph(gf::GeneFamily) -> Union{Nothing, SimpleWeightedGraph}

Undirected `SimpleWeightedGraphs.SimpleWeightedGraph` view of `gf.dN`, or
`nothing` if `dN` wasn't supplied.
"""
dN_graph(gf::GeneFamily) = _weighted_undirected_graph(gf.dN)

"""
    dS_graph(gf::GeneFamily) -> Union{Nothing, SimpleWeightedGraph}

Same as [`dN_graph`](@ref), for `gf.dS`.
"""
dS_graph(gf::GeneFamily) = _weighted_undirected_graph(gf.dS)

_weighted_directed_graph(::Nothing) = nothing
_weighted_directed_graph(m::SparseMatrixCSC) = SimpleWeightedDiGraph(m)

"""
    id_subject_query_graph(gf::GeneFamily) -> Union{Nothing, SimpleWeightedDiGraph}

Directed `SimpleWeightedGraphs.SimpleWeightedDiGraph` view of
`gf.id_subject_query`, or `nothing` if it wasn't supplied. Edge `i -> j` carries
weight `M[i, j]`, following the stored `[query, subject]` axis convention (see
the `GeneFamily` docstring) — no symmetrisation, since the relation is directed.
"""
id_subject_query_graph(gf::GeneFamily) = _weighted_directed_graph(gf.id_subject_query)

"""
    id_query_subject_graph(gf::GeneFamily) -> Union{Nothing, SimpleWeightedDiGraph}

Same as [`id_subject_query_graph`](@ref), for `gf.id_query_subject` (stored
`[subject, query]`).
"""
id_query_subject_graph(gf::GeneFamily) = _weighted_directed_graph(gf.id_query_subject)

#= Reciprocal best hits =#

"""
Canonical name for a scoring alias (`"maximum"` → `"max"`, `"avg"` → `"mean"`,
…), or `nothing` when the alias is unrecognised. Which names a given `rbh`
method accepts is up to that method.
"""
function _canonical_scoring(scoring::AbstractString)
    key = lowercase(scoring)
    key in ("max", "maximum") && return "max"
    key in ("min", "minimum") && return "min"
    key in ("mean", "avg", "average") && return "mean"
    key in ("double_max", "ds") && return key
    return nothing
end

# For edge {a, b}, resolve which end was the original query/subject using
# whichever id_* relation is present (their storage axis records it exactly);
# falls back to the lower index as query when neither is present.
function _edge_query_subject(gf::GeneFamily, a::Integer, b::Integer)
    if !isnothing(gf.id_subject_query)
        return gf.id_subject_query[a, b] != 0 ? (a, b) : (b, a)
    elseif !isnothing(gf.id_query_subject)
        return gf.id_query_subject[b, a] != 0 ? (a, b) : (b, a)
    end
    return minmax(a, b)
end

"""
    edge_identity(gf, query, subject, scoring)

The %ID for edge `(query, subject)`, combining whichever of
`id_subject_query`/`id_query_subject` are present via `scoring` (`"mean"`,
`"min"` or `"max"`). Returns `nothing` when neither relation is present.
"""
function edge_identity(gf::GeneFamily, query::Integer, subject::Integer, scoring::String)
    subject_query =
        isnothing(gf.id_subject_query) ? nothing : gf.id_subject_query[query, subject]
    query_subject =
        isnothing(gf.id_query_subject) ? nothing : gf.id_query_subject[subject, query]
    isnothing(subject_query) && return query_subject
    isnothing(query_subject) && return subject_query
    scoring == "mean" && return (subject_query + query_subject) / 2
    scoring == "min" && return min(subject_query, query_subject)
    return max(subject_query, query_subject)
end

# Group `graph`'s edges by connected component, returning `(components, edges)`
# as parallel vectors. One pass over the edges, rather than one pass per
# component.
function _component_edges(graph::SimpleGraph)
    components = connected_components(graph)
    component_of = zeros(Int, nv(graph))
    for (index, component) in enumerate(components), vertex in component
        component_of[vertex] = index
    end

    grouped = [eltype(edges(graph))[] for _ in components]
    for edge in edges(graph)
        push!(grouped[component_of[src(edge)]], edge)
    end
    return components, grouped
end

"""
    rbh(gf::GeneFamily; scoring::String = "mean") -> DataFrame

Reciprocal best hit per weakly-connected component of `gf.topology`. Within
each component, edges are ranked ascending by `dS`, then `dN`, then %ID
(combining `id_subject_query`/`id_query_subject` via `scoring` — `"mean"`
(default), `"min"`, or `"max"`), skipping whichever levels `gf` lacks. The
top-ranked edge of each multi-gene component is kept; singleton (edge-less)
components are skipped. A tie at the top rank is broken by the lowest gene
indices, and warned about (naming the tied component's gene indices).

Returns a `DataFrame` shaped like `GeneFamily`'s constructor input — `query`,
`subject`, plus whichever of `dN`/`dS`/`id_subject_query`/`id_query_subject`
are present on `gf` — containing only the chosen pairs.
"""
function rbh(gf::GeneFamily; scoring::String = "mean")
    scoring = something(_canonical_scoring(scoring), "")
    scoring in ("mean", "min", "max") ||
        throw(ArgumentError("scoring must be \"mean\", \"min\", or \"max\""))

    any(!isnothing, (gf.dN, gf.dS, gf.id_subject_query, gf.id_query_subject)) ||
        throw(ArgumentError("gf has no dN, dS, or %ID relation to rank edges by"))

    # Ascending ranking key: dS, then dN, then -%ID. A level `gf` lacks
    # contributes a constant, so it never discriminates; `primary` records which
    # slot is the first real one, for the tie warning below. `(lo, hi)` is a
    # final, deterministic tie-break (lowest indices win).
    primary = !isnothing(gf.dS) ? 1 : !isnothing(gf.dN) ? 2 : 3
    function rank_key(a, b)
        lo, hi = minmax(a, b)
        pid = edge_identity(gf, _edge_query_subject(gf, a, b)..., scoring)
        return (
            (
                isnothing(gf.dS) ? 0.0 : gf.dS[lo, hi],
                isnothing(gf.dN) ? 0.0 : gf.dN[lo, hi],
                isnothing(pid) ? 0.0 : -pid,
            ),
            (lo, hi),
        )
    end

    index_to_id = _index_to_id(gf)
    components, component_edges = _component_edges(topology_graph(gf))

    query, subject = String[], String[]
    dN_out, dS_out, sq_out, qs_out = Float64[], Float64[], Float64[], Float64[]
    tied_components = Vector{Int}[]

    for (component, edge_list) in zip(components, component_edges)
        isempty(edge_list) && continue   # singleton gene, no partner to hit

        ranked = sort([(rank_key(src(e), dst(e)), src(e), dst(e)) for e in edge_list])
        best_metric = ranked[1][1][primary]
        count(r -> r[1][primary] == best_metric, ranked) > 1 &&
            push!(tied_components, sort(component))

        _, a, b = ranked[1]
        gene, paralog = _edge_query_subject(gf, a, b)
        push!(query, index_to_id[gene])
        push!(subject, index_to_id[paralog])
        isnothing(gf.dN) || push!(dN_out, gf.dN[minmax(gene, paralog)...])
        isnothing(gf.dS) || push!(dS_out, gf.dS[minmax(gene, paralog)...])
        isnothing(gf.id_subject_query) || push!(sq_out, gf.id_subject_query[gene, paralog])
        isnothing(gf.id_query_subject) || push!(qs_out, gf.id_query_subject[paralog, gene])
    end

    if !isempty(tied_components)
        groups = join(("[$(join(c, ", "))]" for c in tied_components), ", ")
        @warn "rbh: tied best-hit distance in $(length(tied_components)) group(s): $groups"
    end

    columns = Pair{String,Vector}["query"=>query, "subject"=>subject]
    isnothing(gf.dN) || push!(columns, "dN" => dN_out)
    isnothing(gf.dS) || push!(columns, "dS" => dS_out)
    isnothing(gf.id_subject_query) || push!(columns, "id_subject_query" => sq_out)
    isnothing(gf.id_query_subject) || push!(columns, "id_query_subject" => qs_out)
    return DataFrame(columns...)
end

"""
Row and column best hits of a sparse score matrix held as
`(row, column) => score`, over `n_genes` genes. Returns
`(best_in_row, best_in_column)`, each a `Vector{Int}` of gene indices, `0` where
a gene has no scored partner. `better(a, b)` decides whether `a` beats `b`
(`>` for a similarity, `<` for a distance); ties go to the lowest index, so the
result does not depend on iteration order.
"""
function _best_hits(scores::Dict{Tuple{Int,Int},Float64}, n_genes::Int, better::F) where {F}
    best_in_row = zeros(Int, n_genes)
    best_in_column = zeros(Int, n_genes)
    row_best = Vector{Float64}(undef, n_genes)
    column_best = Vector{Float64}(undef, n_genes)

    for ((row, column), score) in scores
        if best_in_row[row] == 0 ||
           better(score, row_best[row]) ||
           (score == row_best[row] && column < best_in_row[row])
            best_in_row[row] = column
            row_best[row] = score
        end
        if best_in_column[column] == 0 ||
           better(score, column_best[column]) ||
           (score == column_best[column] && row < best_in_column[column])
            best_in_column[column] = row
            column_best[column] = score
        end
    end
    return best_in_row, best_in_column
end

# Score matrices for a flat paralog table: `original` holds the two directional
# values as given, `ranked` the values best hits are chosen by (identical to
# `original` under "double_max" and "ds"). Later rows overwrite earlier ones for
# a repeated pair, and both are keyed by the gene indices of `id_to_index`.
function _pair_matrices(
    paralog_df::DataFrame,
    id_to_index::Dict{<:Any,Int},
    scoring::String,
)
    original = Dict{Tuple{Int,Int},Float64}()
    ranked = Dict{Tuple{Int,Int},Float64}()

    for row in eachrow(paralog_df)
        gene = id_to_index[string(row[1])]
        paralog = id_to_index[string(row[2])]
        forward = Float64(row[3])
        backward = scoring == "ds" ? forward : Float64(row[4])

        original[(gene, paralog)] = forward
        original[(paralog, gene)] = backward
        combined =
            scoring == "max" ? max(forward, backward) :
            scoring == "mean" ? (forward + backward) / 2 : nothing
        ranked[(gene, paralog)] = isnothing(combined) ? forward : combined
        ranked[(paralog, gene)] = isnothing(combined) ? backward : combined
    end

    return original, ranked
end

# Genes that are each other's best hit and not already paired, as
# `(gene, best_partner)` pairs in ascending gene order.
function _reciprocal_pairs(
    ranked::Dict{Tuple{Int,Int},Float64},
    n_genes::Int,
    better::F,
) where {F}
    best_in_row, best_in_column = _best_hits(ranked, n_genes, better)
    matched = Set{Int}()
    pairs = Tuple{Int,Int}[]

    for gene = 1:n_genes
        partner = best_in_row[gene]
        (partner == 0 || best_in_column[partner] != gene) && continue
        (gene in matched || partner in matched) && continue
        push!(matched, gene)
        push!(matched, partner)
        push!(pairs, (gene, partner))
    end
    return pairs
end

"""
    rbh(paralog_df::DataFrame; scoring = "max") -> DataFrame

Identify reciprocal best hits (RBH) between paralogs listed one pair per row.

`paralog_df` needs at least three columns (four unless `scoring = "ds"`):

1. `GeneID` (`String`)
2. `ParalogID` (`String`)
3. percent identity gene → paralog, or the pair's dS when `scoring = "ds"`
4. percent identity paralog → gene (ignored when `scoring = "ds"`)

`scoring` picks the value pairs are ranked by:

- `"max"`/`"maximum"` (default): the larger of the two directions
- `"mean"`/`"avg"`/`"average"`: their mean
- `"double_max"`: each direction on its own
- `"ds"`: column 3 read as a dS value and ranked **ascending** (see
  [`rbh_ds`](@ref))

A pair is kept when each gene is the other's best hit and neither has already
been paired. Returns a `DataFrame` of `GeneID`, `ParalogID`, `perc_1`, `perc_2`
(the two original scores), `max_perc` and `mean_perc`; `scoring = "ds"` returns
[`rbh_ds`](@ref)'s columns instead. An empty input gives an empty, correctly
shaped result.

**NOTE:** only scored pairs are ranked. A gene absent from a pair has no best
hit there, rather than being treated as scoring zero against it.
"""
function rbh(paralog_df::DataFrame; scoring::String = "max")
    canonical = _canonical_scoring(scoring)
    canonical in ("max", "mean", "double_max", "ds") || throw(
        ArgumentError(
            "Invalid scoring method. Must be 'ds', 'max', 'maximum', 'double_max', 'mean', 'avg', or 'average'.",
        ),
    )
    canonical == "ds" && return rbh_ds(paralog_df)

    _check_pair_columns(paralog_df, 4)
    ids, id_to_index = _index_pair_ids(paralog_df)
    original, ranked = _pair_matrices(paralog_df, id_to_index, canonical)

    gene_ids, paralog_ids = String[], String[]
    forwards, reverses, maxima, means = Float64[], Float64[], Float64[], Float64[]
    for (gene, partner) in _reciprocal_pairs(ranked, length(ids), >)
        forward = get(original, (gene, partner), 0.0)
        backward = get(original, (partner, gene), 0.0)
        push!(gene_ids, ids[partner])
        push!(paralog_ids, ids[gene])
        push!(forwards, forward)
        push!(reverses, backward)
        push!(maxima, max(forward, backward))
        push!(means, (forward + backward) / 2)
    end

    return DataFrame(
        "GeneID" => gene_ids,
        "ParalogID" => paralog_ids,
        "perc_1" => forwards,
        "perc_2" => reverses,
        "max_perc" => maxima,
        "mean_perc" => means,
    )
end

"""
    rbh_ds(paralog_df::DataFrame) -> DataFrame

Reciprocal best hits ranked by dS — the lowest value wins, dS being a distance
rather than a similarity. Column 3 of `paralog_df` holds the pair's dS and is
taken to apply in both directions; see [`rbh`](@ref) for the pairing rule.

Returns a `DataFrame` of `GeneID`, `ParalogID`, `ds` and `min_ds`.
"""
function rbh_ds(paralog_df::DataFrame)
    _check_pair_columns(paralog_df, 3)
    ids, id_to_index = _index_pair_ids(paralog_df)
    original, ranked = _pair_matrices(paralog_df, id_to_index, "ds")

    gene_ids, paralog_ids = String[], String[]
    ds_values, minima = Float64[], Float64[]
    for (gene, partner) in _reciprocal_pairs(ranked, length(ids), <)
        forward = get(original, (gene, partner), 0.0)
        backward = get(original, (partner, gene), 0.0)
        push!(gene_ids, ids[partner])
        push!(paralog_ids, ids[gene])
        push!(ds_values, forward)
        push!(minima, min(forward, backward))
    end

    return DataFrame(
        "GeneID" => gene_ids,
        "ParalogID" => paralog_ids,
        "ds" => ds_values,
        "min_ds" => minima,
    )
end

# Validate a flat paralog table: two ID columns followed by `n_columns - 2`
# numeric ones.
function _check_pair_columns(paralog_df::DataFrame, n_columns::Int)
    ncol(paralog_df) >= n_columns || throw(
        ArgumentError(
            "`paralog_df` needs at least $n_columns columns, got $(ncol(paralog_df))",
        ),
    )
    for column = 1:2
        eltype(paralog_df[:, column]) <: AbstractString || throw(
            ArgumentError(
                "Column $column (\"$(names(paralog_df)[column])\") must contain gene IDs (strings), got element type $(eltype(paralog_df[:, column]))",
            ),
        )
    end
    for column = 3:n_columns
        eltype(paralog_df[:, column]) <: Real || throw(
            ArgumentError(
                "Column $column (\"$(names(paralog_df)[column])\") must contain scores (numbers), got element type $(eltype(paralog_df[:, column]))",
            ),
        )
    end
    return nothing
end

# Gene IDs in first-appearance order, and the reverse lookup to their indices.
function _index_pair_ids(paralog_df::DataFrame)
    ids = unique(vcat(string.(paralog_df[:, 1]), string.(paralog_df[:, 2])))
    return ids, Dict(id => index for (index, id) in enumerate(ids))
end

export GeneFamily,
    dN_graph,
    dS_graph,
    edge_identity,
    id_query_subject_graph,
    id_subject_query_graph,
    rbh,
    rbh_ds,
    topology_graph

end

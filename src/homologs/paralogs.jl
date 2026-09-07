"""
Homology *within* a genome: paralog groups, the relation matrices over them,
and reciprocal-best-hit detection.

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
    key == "lca_depth" && return :lca_depth
    key == "lca" && return :lca
    return nothing
end

"""The numeric relation fields of a [`ParalogGroup`](@ref), in declaration order."""
const RELATIONS = (:dN, :dS, :id_subject_query, :id_query_subject, :lca_depth)

"""
The %ID relations, whose value decides whether a pair is connected at all: a
row scoring 0 in every one of them is dropped (see the `ParalogGroup`
constructor).
"""
const IDENTITY_RELATIONS = (:id_subject_query, :id_query_subject)

# Marker distinguishing the raw-fields constructor (below) from the public
# `ParalogGroup(genome, pairs)` one; kept unexported so ordinary callers can't
# reach it and reintroduce arbitrary construction.
struct _RawFields end

struct ParalogGroup
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

    # How deep the pair's duplication sits in the hierarchy `lca` names.
    lca_depth::Union{Nothing,SparseMatrixCSC{Float64,UInt32}}

    # Last common ancestor of a pair, as a 1-based index into `lca_labels`;
    # a stored 0 is impossible, so an unstored cell means "unknown".
    lca::Union{Nothing,SparseMatrixCSC{UInt32,UInt32}}

    # Interned LCA names, in first-appearance order; empty when `lca` is nothing.
    lca_labels::Vector{String}

    """
        ParalogGroup(genome::Reference.Genome, pairs::DataFrame)

    Build a `ParalogGroup` from a reference `genome` and a paralog-pair table.

    `pairs` has 2-8 columns; the first two are required and taken by position,
    the rest optional and matched to the relation matrices *by name*
    (case-insensitive):

    1. Query ID (`String`; mandatory)
    2. Subject ID (`String`; mandatory)
    3. `dN`               → pairwise dN (Float; optional)
    4. `dS`               → pairwise dS (Float; optional)
    5. `id_subject_query` → % identity subject → query (Float; optional)
    6. `id_query_subject` → % identity query → subject (Float; optional)
    7. `lca`              → last common ancestor of the pair (`String`; optional)
    8. `lca_depth`        → depth of that ancestor (Float; optional)

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
    - `lca_depth`        — symmetric and upper-triangular like dN/dS: the depth
      of the speciation event preceding the pair's duplication, higher being
      more recent. A paralog group lives in one species, whose ancestors form a
      single lineage, so every duplication maps onto that one chain and depths
      **are** comparable across pairs — equal depth is the same ancestral
      lineage. Number the root 1: 0 is reserved for "not recorded", and
      [`rbh`](@ref) refuses to rank by it.
    - `lca`              — symmetric, so stored in the upper triangle like
      dN/dS, but holding a 1-based index into `lca_labels` rather than a value;
      read it back with [`lca_label`](@ref). An empty label is taken as unknown
      and left unstored, and a repeated pair keeps the first label seen.

    Only genes that end up **weakly connected** join the group: a row is kept
    when both IDs resolve against `genome` and the pair shows some similarity,
    and a gene is indexed only if it appears in a kept row. Two rules drop rows
    up front, each warning about how many it took:

    - a `NaN` in any matched numeric column — `NaN != NaN` would break the
      symmetry the `*_graph` views check;
    - a 0 in *every* matched %ID column, i.e. no alignment in either direction.
      A pair with no %ID column at all is taken on trust.

    Repeated or reciprocal entries for a cell collapse to one value.

    See `getindex` for per-gene (`pg[i]`/`pg["id"]`) and sub-group
    (`pg[indices]`/`pg[ids]`) indexing, and `topology_graph`/`dN_graph`/
    `dS_graph`/`id_subject_query_graph`/`id_query_subject_graph` for `Graphs`/
    `SimpleWeightedGraphs` views of the relation matrices.
    """
    function ParalogGroup(genome::Reference.Genome, pairs::DataFrame)
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
        # `lca` is categorical, so it is kept apart from the numeric relations.
        value_columns = Dict{Symbol,Int}()
        lca_column = nothing
        for column = 3:ncol(pairs)
            field = _match_relation(names(pairs)[column])
            isnothing(field) && continue
            field === :lca ? (lca_column = column) : (value_columns[field] = column)
        end

        for (field, column) in value_columns
            eltype(pairs[:, column]) <: Real || throw(
                ArgumentError(
                    "Column \"$(names(pairs)[column])\" (matched to `$field`) must contain numeric values, got element type $(eltype(pairs[:, column]))",
                ),
            )
        end

        if !isnothing(lca_column)
            eltype(pairs[:, lca_column]) <: AbstractString || throw(
                ArgumentError(
                    "Column \"$(names(pairs)[lca_column])\" (matched to `lca`) must contain ancestor names (strings), got element type $(eltype(pairs[:, lca_column]))",
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
                @warn "ParalogGroup: dropped $n_dropped row$(n_dropped == 1 ? "" : "s") with a NaN relation value"
                pairs = pairs[keep, :]
            end
        end

        # A pair scoring 0 in every %ID direction it reports is no pair at all.
        identity_columns =
            [value_columns[f] for f in IDENTITY_RELATIONS if haskey(value_columns, f)]
        if !isempty(identity_columns)
            keep = falses(nrow(pairs))
            for column in identity_columns
                keep .|= .!iszero.(pairs[:, column])
            end
            n_dropped = count(!, keep)
            if n_dropped > 0
                @warn "ParalogGroup: dropped $n_dropped row$(n_dropped == 1 ? "" : "s") with 0% identity in every direction"
                pairs = pairs[keep, :]
            end
        end

        query_ids = string.(pairs[:, 1])
        subject_ids = string.(pairs[:, 2])
        lca_names = isnothing(lca_column) ? String[] : string.(pairs[:, lca_column])
        # Materialised up front: indexing a DataFrame cell per row is both
        # type-unstable and far slower than indexing a `Vector{Float64}`.
        values_by_field = Dict{Symbol,Vector{Float64}}(
            field => Float64.(pairs[:, column]) for (field, column) in value_columns
        )

        # Resolve every ID to its feature.
        records = genome[unique(vcat(query_ids, subject_ids))]

        # A gene joins the group only through a pair whose *other* end also
        # resolved, so no isolated vertex can reach the matrices below.
        resolved = Set(record.id for record in records)
        connected = Set{String}()
        for row in eachindex(query_ids)
            query, subject = query_ids[row], subject_ids[row]
            (query in resolved && subject in resolved) || continue
            push!(connected, query)
            push!(connected, subject)
        end

        # Group the found features by scaffold, discarding repeated IDs.
        scaffold_intervals = Dict{String,Vector{Reference.IntervalSimple}}()
        scaffold_id_lists = Dict{String,Vector{String}}()
        seen = Set{String}()
        for record in records
            (record.id in connected && !(record.id in seen)) || continue
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
        dN_v, dS_v, depth_v = Float64[], Float64[], Float64[]
        sq_i, sq_j, sq_v = UInt32[], UInt32[], Float64[]
        qs_i, qs_j, qs_v = UInt32[], UInt32[], Float64[]
        lca_i, lca_j, lca_v = UInt32[], UInt32[], UInt32[]
        lca_labels = String[]
        lca_codes = Dict{String,UInt32}()

        want_tri = has(:dN) || has(:dS) || has(:lca_depth)
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
                has(:lca_depth) && push!(depth_v, values_by_field[:lca_depth][row])
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
            if !isnothing(lca_column) && !isempty(lca_names[row])
                code = get!(lca_codes, lca_names[row]) do
                    push!(lca_labels, lca_names[row])
                    UInt32(length(lca_labels))
                end
                lo, hi = minmax(query_index, subject_index)
                push!(lca_i, lo)
                push!(lca_j, hi)
                push!(lca_v, code)
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
            has(:lca_depth) ? relation(tri_i, tri_j, depth_v) : nothing,
            # `min` keeps the first label seen, codes being assigned in order.
            isnothing(lca_column) ? nothing :
            sparse(lca_i, lca_j, lca_v, n_genes, n_genes, min),
            lca_labels,
        )
    end

    # Internal-only: build directly from already-resolved fields, used by the
    # sub-group `getindex` methods below. The `_RawFields` marker keeps this
    # unreachable from `ParalogGroup(...)` calls outside this module.
    function ParalogGroup(
        ::_RawFields,
        intervals::Dict{String,StructArray{Reference.IntervalSimple}},
        scaffold_ranges::Dict{String,UnitRange{Int}},
        id_to_index::Dict{String,Int},
        topology::SparseMatrixCSC{Bool,UInt32},
        dN,
        dS,
        id_subject_query,
        id_query_subject,
        lca_depth,
        lca,
        lca_labels::Vector{String},
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
            lca_depth,
            lca,
            lca_labels,
        )
    end
end

"""Number of genes a group spans."""
_n_genes(pg::ParalogGroup) = size(pg.topology, 1)

"""Reverse of `pg.id_to_index`: linear gene index to gene ID."""
_index_to_id(pg::ParalogGroup) = Dict(index => id for (id, index) in pg.id_to_index)

function Base.show(io::IO, pg::ParalogGroup)
    n = _n_genes(pg)
    n_scaffolds = length(pg.scaffold_ranges)
    n_pairs = nnz(pg.topology) ÷ 2
    present = [String(field) for field in RELATIONS if !isnothing(getfield(pg, field))]
    isnothing(pg.lca) || push!(present, "lca")
    print(
        io,
        "ParalogGroup($n gene$(n == 1 ? "" : "s"), $n_scaffolds scaffold$(n_scaffolds == 1 ? "" : "s"), $n_pairs pair$(n_pairs == 1 ? "" : "s"); relations: $(isempty(present) ? "none" : join(present, ", ")))",
    )
end

"""
    pg[i::Integer]
    pg[id::AbstractString]

Return column `i` (or the column for gene `id`) from every relation matrix
present on `pg`, as a `Dict{Symbol, Vector}` keyed by `:topology` and whichever
of `:dN`, `:dS`, `:id_subject_query`, `:id_query_subject`, `:lca_depth`, `:lca`
are non-`nothing`. The `:lca` entry holds the ancestor *names*, `""` where the
pair has none.

`dN`/`dS`/`lca_depth`/`lca` store only their upper triangle (see the
`ParalogGroup` docstring), so column `i` alone only surfaces partners with a
smaller index; combine with row `i` for the complete symmetric relation.

A relation stores an explicit zero, so a recorded 0 is a real value here; a 0
for a pair that was never recorded looks identical, and `:topology` in the same
dict is the authority on which pairs exist at all.
"""
function Base.getindex(pg::ParalogGroup, i::Integer)
    n = _n_genes(pg)
    1 <= i <= n || throw(ArgumentError("gene index $i out of bounds (1:$n)"))
    columns = Dict{Symbol,Vector}(:topology => Vector(pg.topology[:, i]))
    for field in RELATIONS
        matrix = getfield(pg, field)
        isnothing(matrix) || (columns[field] = Vector(matrix[:, i]))
    end
    isnothing(pg.lca) ||
        (columns[:lca] = [iszero(code) ? "" : pg.lca_labels[code] for code in pg.lca[:, i]])
    return columns
end

function Base.getindex(pg::ParalogGroup, id::AbstractString)
    haskey(pg.id_to_index, id) || throw(ArgumentError("Unknown gene ID: \"$id\""))
    return pg[pg.id_to_index[id]]
end

"""
    pg[indices::AbstractVector{<:Integer}]
    pg[ids::AbstractVector{<:AbstractString}]

Return a new `ParalogGroup` restricted to the given genes (a "sub-group"): every
relation matrix present on `pg` is sliced to just those genes, and
`intervals`/`scaffold_ranges` are rebuilt to stay concordant with the new,
contiguous linear indexing. `lca_labels` is carried over whole, so a label the
sub-group no longer uses simply goes unreferenced.
"""
function Base.getindex(pg::ParalogGroup, indices::AbstractVector{<:Integer})
    isempty(indices) && throw(ArgumentError("`indices` must be non-empty"))
    n = _n_genes(pg)
    all(i -> 1 <= i <= n, indices) ||
        throw(ArgumentError("gene index out of bounds (1:$n)"))

    order = sort!(unique(collect(Int, indices)))
    order_set = Set(order)
    index_to_id = _index_to_id(pg)

    new_intervals = Dict{String,StructArray{Reference.IntervalSimple}}()
    new_ranges = Dict{String,UnitRange{Int}}()
    new_id_to_index = Dict{String,Int}()
    cursor = 1
    for name in sort!(collect(keys(pg.scaffold_ranges)))
        range = pg.scaffold_ranges[name]
        kept = filter(in(order_set), range)
        isempty(kept) && continue

        new_intervals[name] =
            StructArray([pg.intervals[name][gene-first(range)+1] for gene in kept])
        new_ranges[name] = cursor:(cursor+length(kept)-1)
        for (offset, gene) in enumerate(kept)
            new_id_to_index[index_to_id[gene]] = cursor + offset - 1
        end
        cursor += length(kept)
    end

    function relation(field)
        matrix = getfield(pg, field)
        return isnothing(matrix) ? nothing : matrix[order, order]
    end

    return ParalogGroup(
        _RawFields(),
        new_intervals,
        new_ranges,
        new_id_to_index,
        pg.topology[order, order],
        relation(:dN),
        relation(:dS),
        relation(:id_subject_query),
        relation(:id_query_subject),
        relation(:lca_depth),
        relation(:lca),
        pg.lca_labels,
    )
end

function Base.getindex(pg::ParalogGroup, ids::AbstractVector{<:AbstractString})
    isempty(ids) && throw(ArgumentError("`ids` must be non-empty"))
    indices = map(ids) do id
        haskey(pg.id_to_index, id) || throw(ArgumentError("Unknown gene ID: \"$id\""))
        pg.id_to_index[id]
    end
    return pg[indices]
end

"""
    lca_label(pg::ParalogGroup, a, b) -> Union{Nothing, String}

Name of the last common ancestor recorded for the pair of genes `a` and `b`
(linear indices or gene IDs), or `nothing` when `pg` has no `lca` relation or
the pair has no ancestor recorded. The relation is symmetric, so argument order
does not matter.
"""
function lca_label(pg::ParalogGroup, a::Integer, b::Integer)
    isnothing(pg.lca) && return nothing
    n = _n_genes(pg)
    all(i -> 1 <= i <= n, (a, b)) || throw(ArgumentError("gene index out of bounds (1:$n)"))
    code = pg.lca[minmax(a, b)...]
    return iszero(code) ? nothing : pg.lca_labels[code]
end

function lca_label(pg::ParalogGroup, a::AbstractString, b::AbstractString)
    for id in (a, b)
        haskey(pg.id_to_index, id) || throw(ArgumentError("Unknown gene ID: \"$id\""))
    end
    return lca_label(pg, pg.id_to_index[a], pg.id_to_index[b])
end

#= Graph views =#

"""
    topology_graph(pg::ParalogGroup) -> SimpleGraph

Unweighted `Graphs.SimpleGraph` view of `pg.topology`: vertex `v` is gene `v`
(see `pg.id_to_index` for the name), and edge `u-v` marks a topology pair.
"""
topology_graph(pg::ParalogGroup) = SimpleGraph(pg.topology)

_weighted_undirected_graph(::Nothing) = nothing
# dN/dS store only their upper triangle; SimpleWeightedGraph requires a fully
# symmetric adjacency matrix, so mirror it across the diagonal first.
_weighted_undirected_graph(m::SparseMatrixCSC) = SimpleWeightedGraph(m + permutedims(m))

"""
    dN_graph(pg::ParalogGroup) -> Union{Nothing, SimpleWeightedGraph}

Undirected `SimpleWeightedGraphs.SimpleWeightedGraph` view of `pg.dN`, or
`nothing` if `dN` wasn't supplied.

A weighted graph represents "no edge" as weight 0, and the symmetrisation above
prunes zeros besides, so a pair whose value is 0 — a dN of 0 between recent
duplicates, say — has no edge here. Use `pg.topology` to ask which pairs exist.
"""
dN_graph(pg::ParalogGroup) = _weighted_undirected_graph(pg.dN)

"""
    dS_graph(pg::ParalogGroup) -> Union{Nothing, SimpleWeightedGraph}

Same as [`dN_graph`](@ref), for `pg.dS`.
"""
dS_graph(pg::ParalogGroup) = _weighted_undirected_graph(pg.dS)

_weighted_directed_graph(::Nothing) = nothing
_weighted_directed_graph(m::SparseMatrixCSC) = SimpleWeightedDiGraph(m)

"""
    id_subject_query_graph(pg::ParalogGroup) -> Union{Nothing, SimpleWeightedDiGraph}

Directed `SimpleWeightedGraphs.SimpleWeightedDiGraph` view of
`pg.id_subject_query`, or `nothing` if it wasn't supplied. Edge `i -> j` carries
weight `M[i, j]`, following the stored `[query, subject]` axis convention (see
the `ParalogGroup` docstring) — no symmetrisation, since the relation is directed.
"""
id_subject_query_graph(pg::ParalogGroup) = _weighted_directed_graph(pg.id_subject_query)

"""
    id_query_subject_graph(pg::ParalogGroup) -> Union{Nothing, SimpleWeightedDiGraph}

Same as [`id_subject_query_graph`](@ref), for `pg.id_query_subject` (stored
`[subject, query]`).
"""
id_query_subject_graph(pg::ParalogGroup) = _weighted_directed_graph(pg.id_query_subject)

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
function _edge_query_subject(pg::ParalogGroup, a::Integer, b::Integer)
    if !isnothing(pg.id_subject_query)
        return pg.id_subject_query[a, b] != 0 ? (a, b) : (b, a)
    elseif !isnothing(pg.id_query_subject)
        return pg.id_query_subject[b, a] != 0 ? (a, b) : (b, a)
    end
    return minmax(a, b)
end

"""
    edge_identity(pg, query, subject, scoring)

The %ID for edge `(query, subject)`, combining whichever of
`id_subject_query`/`id_query_subject` are present via `scoring` (`"mean"`,
`"min"` or `"max"`). Returns `nothing` when neither relation is present.
"""
function edge_identity(pg::ParalogGroup, query::Integer, subject::Integer, scoring::String)
    subject_query =
        isnothing(pg.id_subject_query) ? nothing : pg.id_subject_query[query, subject]
    query_subject =
        isnothing(pg.id_query_subject) ? nothing : pg.id_query_subject[subject, query]
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
The relations [`rbh`](@ref) can rank edges by. `:identity` covers the
`id_subject_query`/`id_query_subject` pair, combined via `rbh`'s `scoring`.
"""
const RANK_LEVELS = (:dS, :dN, :identity, :lca_depth)

"""Whether `pg` carries the relation a [`RANK_LEVELS`](@ref) entry needs."""
_has_level(pg::ParalogGroup, level::Symbol) =
    level === :identity ?
    (!isnothing(pg.id_subject_query) || !isnothing(pg.id_query_subject)) :
    !isnothing(getfield(pg, level))

# `rbh`'s ranking levels: everything `pg` carries bar `:lca_depth`, which is
# opt-in, or the caller's own list, which must name only relations `pg` has.
function _rank_levels(pg::ParalogGroup, levels::Nothing)
    available = filter(level -> _has_level(pg, level), [:dS, :dN, :identity])
    isempty(available) && throw(
        ArgumentError(
            "pg has no dN, dS, or %ID relation to rank edges by; pass `levels` to rank by another relation",
        ),
    )
    return available
end

function _rank_levels(pg::ParalogGroup, levels)
    requested = collect(Symbol, levels)
    isempty(requested) && throw(ArgumentError("`levels` must name at least one relation"))
    allunique(requested) ||
        throw(ArgumentError("`levels` names the same relation more than once"))
    for level in requested
        level in RANK_LEVELS || throw(
            ArgumentError(
                "Unknown ranking level `$level`; expected one of $(join(RANK_LEVELS, ", "))",
            ),
        )
        _has_level(pg, level) ||
            throw(ArgumentError("pg carries no `$level` relation to rank edges by"))
    end
    return requested
end

# One level's contribution to the ascending rank key for edge {a, b}. Similarity
# levels are negated so that, throughout, smaller means better.
function _level_value(pg::ParalogGroup, level::Symbol, a::Integer, b::Integer, scoring)
    lo, hi = minmax(a, b)
    level === :dS && return pg.dS[lo, hi]
    level === :dN && return pg.dN[lo, hi]
    level === :lca_depth && return -pg.lca_depth[lo, hi]
    return -edge_identity(pg, _edge_query_subject(pg, a, b)..., scoring)
end

"""
    rbh(pg::ParalogGroup; scoring::String = "mean", levels = nothing) -> DataFrame

Reciprocal best hit per weakly-connected component of `pg.topology`. Within each
component, edges are ranked by `levels` — the first entry primary, each
subsequent entry breaking ties in the one before it — and the top-ranked edge of
each multi-gene component is kept; singleton (edge-less) components are skipped.
An edge tying on *every* level is resolved by the lowest gene indices, and
warned about (naming the tied component's gene indices).

`levels` names the relations to rank by, most significant first:

- `:dS`, `:dN`  — ascending: the smaller distance is the closer paralog.
- `:identity`   — descending, combining `id_subject_query`/`id_query_subject`
  via `scoring` (`"mean"` (default), `"min"`, or `"max"`).
- `:lca_depth`  — descending: the deeper duplication is the more recent one.

`levels = nothing` (the default) uses whichever of `:dS`, `:dN`, `:identity`
`pg` carries, in that order. `:lca_depth` is never included implicitly — ask for
it, e.g. `levels = [:lca_depth, :dS]` to let the phylogeny outrank the molecular
clock. A missing rank is never silently defaulted: naming a level `pg` lacks is
an error, as is ranking by `:lca_depth` when any edge has no depth recorded.

Returns a `DataFrame` shaped like `ParalogGroup`'s constructor input — `query`,
`subject`, plus whichever of `dN`/`dS`/`id_subject_query`/`id_query_subject`/
`lca_depth`/`lca` are present on `pg` — containing only the chosen pairs.
"""
function rbh(pg::ParalogGroup; scoring::String = "mean", levels = nothing)
    scoring = something(_canonical_scoring(scoring), "")
    scoring in ("mean", "min", "max") ||
        throw(ArgumentError("scoring must be \"mean\", \"min\", or \"max\""))

    ranking = _rank_levels(pg, levels)

    index_to_id = _index_to_id(pg)
    graph = topology_graph(pg)

    # Depths start at 1, so a 0 (stored or not) means "not recorded".
    if :lca_depth in ranking
        unranked = [
            (index_to_id[src(e)], index_to_id[dst(e)]) for
            e in edges(graph) if iszero(pg.lca_depth[minmax(src(e), dst(e))...])
        ]
        isempty(unranked) || throw(
            ArgumentError(
                "rbh: ranking by `lca_depth`, but $(length(unranked)) pair(s) have none recorded, e.g. $(first(unranked))",
            ),
        )
    end

    # Ascending throughout (see `_level_value`), with `(lo, hi)` as a final,
    # deterministic tie-break: lowest indices win.
    rank_key(a, b) =
        ([_level_value(pg, level, a, b, scoring) for level in ranking], minmax(a, b))

    components, component_edges = _component_edges(graph)

    query, subject = String[], String[]
    dN_out, dS_out, sq_out, qs_out = Float64[], Float64[], Float64[], Float64[]
    depth_out = Float64[]
    lca_out = String[]
    tied_components = Vector{Int}[]

    for (component, edge_list) in zip(components, component_edges)
        isempty(edge_list) && continue   # singleton gene, no partner to hit

        ranked = sort([(rank_key(src(e), dst(e)), src(e), dst(e)) for e in edge_list])
        # Only a tie on every level is arbitrary; a lower level resolving the
        # top one is the ranking working as asked.
        best_key = ranked[1][1][1]
        count(r -> r[1][1] == best_key, ranked) > 1 &&
            push!(tied_components, sort(component))

        _, a, b = ranked[1]
        gene, paralog = _edge_query_subject(pg, a, b)
        push!(query, index_to_id[gene])
        push!(subject, index_to_id[paralog])
        isnothing(pg.dN) || push!(dN_out, pg.dN[minmax(gene, paralog)...])
        isnothing(pg.dS) || push!(dS_out, pg.dS[minmax(gene, paralog)...])
        isnothing(pg.id_subject_query) || push!(sq_out, pg.id_subject_query[gene, paralog])
        isnothing(pg.id_query_subject) || push!(qs_out, pg.id_query_subject[paralog, gene])
        isnothing(pg.lca_depth) || push!(depth_out, pg.lca_depth[minmax(gene, paralog)...])
        isnothing(pg.lca) || push!(lca_out, something(lca_label(pg, gene, paralog), ""))
    end

    if !isempty(tied_components)
        groups = join(("[$(join(c, ", "))]" for c in tied_components), ", ")
        @warn "rbh: tied best-hit distance in $(length(tied_components)) group(s): $groups"
    end

    columns = Pair{String,Vector}["query"=>query, "subject"=>subject]
    isnothing(pg.dN) || push!(columns, "dN" => dN_out)
    isnothing(pg.dS) || push!(columns, "dS" => dS_out)
    isnothing(pg.id_subject_query) || push!(columns, "id_subject_query" => sq_out)
    isnothing(pg.id_query_subject) || push!(columns, "id_query_subject" => qs_out)
    isnothing(pg.lca_depth) || push!(columns, "lca_depth" => depth_out)
    isnothing(pg.lca) || push!(columns, "lca" => lca_out)
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

export ParalogGroup,
    dN_graph,
    dS_graph,
    edge_identity,
    id_query_subject_graph,
    id_subject_query_graph,
    lca_label,
    rbh,
    rbh_ds,
    topology_graph

end

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

#= Shared helpers =#

"""
Given a `count` and a singular `noun`, phrase the two together. Returns e.g.
`"1 pair"` or `"3 pairs"`.
"""
_plural(count::Integer, noun::AbstractString) = "$count $noun$(count == 1 ? "" : "s")"

"""
Given a table `column` required to hold `T` values (described by `expected`),
check its element type. Returns `nothing`, or throws an `ArgumentError` naming
the column after `role`.
"""
function _require_eltype(
    df::DataFrame,
    column::Integer,
    ::Type{T},
    expected::AbstractString,
    role::AbstractString = "Column",
) where {T}
    eltype(df[:, column]) <: T && return nothing
    throw(
        ArgumentError(
            "$role $column (\"$(names(df)[column])\") must contain $expected, got element type $(eltype(df[:, column]))",
        ),
    )
end

"""
Given an interning table, look `name` up in it, appending it when new. Returns
its 1-based code, so that a 0 is never a valid code.
"""
function _intern!(labels::Vector{String}, codes::Dict{String,UInt32}, name::AbstractString)
    return get!(codes, String(name)) do
        push!(labels, String(name))
        UInt32(length(labels))
    end
end

"""Given a `pairs` column name, return the relation it maps to, or `nothing`."""
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

#= Relation matrix construction =#

"""
Given coordinates and values over `n_genes` genes, build a relation matrix.
Returns a `SparseMatrixCSC` in which `combine` settles a repeated cell and an
explicit zero survives, so a meaningful dN or dS of 0 is not pruned.
"""
_relation(rows, cols, vals, n_genes::Integer, combine = max) =
    sparse(UInt32.(rows), UInt32.(cols), vals, n_genes, n_genes, combine)

"""
Given one end of each edge in `edge_i`/`edge_j`, build the adjacency over
`n_genes` genes. Returns a symmetric `SparseMatrixCSC{Bool}` setting both
`[a, b]` and `[b, a]`.
"""
_topology(edge_i, edge_j, n_genes::Integer) = sparse(
    UInt32.(vcat(edge_i, edge_j)),
    UInt32.(vcat(edge_j, edge_i)),
    trues(2 * length(edge_i)),
    n_genes,
    n_genes,
    |,
)

"""
Given two gene indices, read the single cell a symmetric relation stores their
pair in. Returns `matrix[min(a, b), max(a, b)]`.
"""
@inline _tri(matrix::SparseMatrixCSC, a::Integer, b::Integer) = matrix[minmax(a, b)...]

#= ParalogGroup =#

"""
Marker distinguishing the raw-fields `ParalogGroup` constructor from the public
`ParalogGroup(genome, pairs)` one. Kept unexported, so ordinary callers cannot
reach it and reintroduce arbitrary construction.
"""
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
        ParalogGroup(genome, pairs::DataFrame)

    Given a reference `genome` and a table of paralog pairs, index the genes
    those pairs connect and build the relation matrices over them. Returns a
    `ParalogGroup`. `docs/paralog_group.md` lists the column names matched out
    of columns 3-8 and the layout each relation matrix follows.

    Rows carrying a `NaN`, or 0% identity in every direction, are dropped with a
    warning, and a gene is indexed only once a pair leaves it with a partner.
    """
    function ParalogGroup(genome::Reference.Genome, pairs::DataFrame)
        ncol(pairs) >= 2 ||
            throw(ArgumentError("`pairs` needs at least a query and a subject column"))
        _require_eltype(pairs, 1, AbstractString, "gene IDs (strings)", "Query ID column")
        _require_eltype(pairs, 2, AbstractString, "gene IDs (strings)", "Subject ID column")

        value_columns, lca_column = _pair_value_columns(pairs)
        pairs = _filter_pair_rows(pairs, value_columns)

        query_ids = string.(pairs[:, 1])
        subject_ids = string.(pairs[:, 2])
        # Materialised up front: indexing a DataFrame cell per row is both
        # type-unstable and far slower than indexing a `Vector{Float64}`.
        values_by_field = Dict{Symbol,Vector{Float64}}(
            field => Float64.(pairs[:, column]) for (field, column) in value_columns
        )

        records = _connected_records(genome, query_ids, subject_ids)
        intervals, scaffold_ranges, id_to_index = _layout_genes(
            (
                record.chromosome,
                Reference.IntervalSimple(record.start_pos, record.end_pos, record.code),
                record.id,
            ) for record in records
        )
        n_genes = length(id_to_index)

        queries, subjects, rows = _pair_edges(query_ids, subject_ids, id_to_index)
        lower, upper = _upper_triangle(queries, subjects)

        # Symmetric relations fold into the shared upper triangle; the directed
        # ones are oriented so the axis you look up by is the fast (column) one.
        symmetric(field) =
            haskey(value_columns, field) ?
            _relation(lower, upper, values_by_field[field][rows], n_genes) : nothing
        directed(field, matrix_rows, matrix_cols) =
            haskey(value_columns, field) ?
            _relation(matrix_rows, matrix_cols, values_by_field[field][rows], n_genes) :
            nothing

        lca, lca_labels =
            isnothing(lca_column) ? (nothing, String[]) :
            _lca_relation(lower, upper, string.(pairs[:, lca_column]), rows, n_genes)

        return new(
            intervals,
            scaffold_ranges,
            id_to_index,
            _topology(queries, subjects, n_genes),
            symmetric(:dN),
            symmetric(:dS),
            directed(:id_subject_query, queries, subjects),
            directed(:id_query_subject, subjects, queries),
            symmetric(:lca_depth),
            lca,
            lca_labels,
        )
    end

    """
    Given already-resolved fields, build a group directly. Returns it; the
    [`_RawFields`](@ref) marker keeps this unreachable from `ParalogGroup(...)`
    calls outside the module, so [`_rebuild`](@ref) is the only way in.
    """
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

"""
Given a per-relation `relation(field)` transform, build a group carrying `pg`'s
identity. Returns a new `ParalogGroup`, taking every field not passed as a
keyword from `pg` unchanged.
"""
function _rebuild(
    pg::ParalogGroup,
    relation::F;
    intervals = pg.intervals,
    scaffold_ranges = pg.scaffold_ranges,
    id_to_index = pg.id_to_index,
    topology = pg.topology,
    lca_labels = pg.lca_labels,
) where {F}
    return ParalogGroup(
        _RawFields(),
        intervals,
        scaffold_ranges,
        id_to_index,
        topology,
        relation(:dN),
        relation(:dS),
        relation(:id_subject_query),
        relation(:id_query_subject),
        relation(:lca_depth),
        relation(:lca),
        lca_labels,
    )
end

#= Constructor stages =#

"""
Given a `pairs` table, match its optional columns to the relations they carry
and check their element types. Returns `(value_columns, lca_column)`, a
`Symbol => column` map of the numeric relations plus the categorical `lca`
column, which is `nothing` when absent.
"""
function _pair_value_columns(pairs::DataFrame)
    value_columns = Dict{Symbol,Int}()
    lca_column = nothing
    for column = 3:ncol(pairs)
        field = _match_relation(names(pairs)[column])
        isnothing(field) && continue
        field === :lca ? (lca_column = column) : (value_columns[field] = column)
    end

    for column in values(value_columns)
        _require_eltype(pairs, column, Real, "numbers", "Relation column")
    end
    isnothing(lca_column) || _require_eltype(
        pairs,
        lca_column,
        AbstractString,
        "ancestor names (strings)",
        "Ancestor column",
    )
    return value_columns, lca_column
end

"""
Given a row mask and the `reason` the rejected rows failed it, restrict `pairs`
to the rows `keep` admits. Returns the survivors, warning with a count whenever
anything is dropped.
"""
function _drop_rows(pairs::DataFrame, keep::AbstractVector{Bool}, reason::AbstractString)
    n_dropped = count(!, keep)
    n_dropped == 0 && return pairs
    @warn "ParalogGroup: dropped $(_plural(n_dropped, "row")) $reason"
    return pairs[keep, :]
end

"""
Given the matched numeric columns of `pairs`, drop the rows that describe no
pair at all. Returns the surviving rows.

A `NaN` would break the symmetry the `*_graph` views check, and a 0 in every
matched %ID column is no alignment in either direction; a table carrying no %ID
column is taken on trust.
"""
function _filter_pair_rows(pairs::DataFrame, value_columns::Dict{Symbol,Int})
    if !isempty(value_columns)
        keep = trues(nrow(pairs))
        for column in values(value_columns)
            keep .&= .!isnan.(pairs[:, column])
        end
        pairs = _drop_rows(pairs, keep, "with a NaN relation value")
    end

    identity_columns =
        [value_columns[f] for f in IDENTITY_RELATIONS if haskey(value_columns, f)]
    if !isempty(identity_columns)
        keep = falses(nrow(pairs))
        for column in identity_columns
            keep .|= .!iszero.(pairs[:, column])
        end
        pairs = _drop_rows(pairs, keep, "with 0% identity in every direction")
    end
    return pairs
end

"""
Given every row's query and subject ID, resolve them against `genome` and keep
the features paired with another resolved one. Returns those `FeatureRecord`s,
each at most once, so no isolated vertex can reach the relation matrices.
"""
function _connected_records(genome::Reference.Genome, query_ids, subject_ids)
    records = genome[unique(vcat(query_ids, subject_ids))]
    resolved = Set(record.id for record in records)

    connected = Set{String}()
    for row in eachindex(query_ids)
        query, subject = query_ids[row], subject_ids[row]
        (query in resolved && subject in resolved) || continue
        push!(connected, query)
        push!(connected, subject)
    end

    seen = Set{String}()
    return filter(records) do record
        (record.id in connected && !(record.id in seen)) || return false
        push!(seen, record.id)
        return true
    end
end

"""
Given `(scaffold, interval, id)` triples, lay the scaffolds out end to end and
order the genes within each by position. Returns
`(intervals, scaffold_ranges, id_to_index)`, a group's linear gene indexing.
"""
function _layout_genes(genes)
    Entry = Tuple{Reference.IntervalSimple,String}
    by_scaffold = Dict{String,Vector{Entry}}()
    for (scaffold, interval, id) in genes
        push!(get!(by_scaffold, scaffold, Entry[]), (interval, id))
    end

    intervals = Dict{String,StructArray{Reference.IntervalSimple}}()
    scaffold_ranges = Dict{String,UnitRange{Int}}()
    id_to_index = Dict{String,Int}()
    cursor = 1
    for name in sort!(collect(keys(by_scaffold)))
        on_scaffold = sort!(
            by_scaffold[name];
            by = entry -> (entry[1].start_pos, entry[1].end_pos, entry[2]),
        )
        scaffold_ranges[name] = cursor:(cursor+length(on_scaffold)-1)
        intervals[name] =
            StructArray(Reference.IntervalSimple[interval for (interval, _) in on_scaffold])
        for (offset, (_, id)) in enumerate(on_scaffold)
            id_to_index[id] = cursor + offset - 1
        end
        cursor += length(on_scaffold)
    end
    return intervals, scaffold_ranges, id_to_index
end

"""
Given every row's query and subject ID, resolve the rows whose ends are both
indexed. Returns `(query_indices, subject_indices, rows)` as parallel vectors.
"""
function _pair_edges(query_ids, subject_ids, id_to_index::Dict{String,Int})
    queries, subjects, rows = Int[], Int[], Int[]
    for row in eachindex(query_ids)
        query = get(id_to_index, query_ids[row], 0)
        subject = get(id_to_index, subject_ids[row], 0)
        (query == 0 || subject == 0) && continue
        push!(queries, query)
        push!(subjects, subject)
        push!(rows, row)
    end
    return queries, subjects, rows
end

"""
Given the two ends of each pair, order them. Returns `(lower, upper)`, the
upper-triangle cell a symmetric relation stores each pair in.
"""
_upper_triangle(queries::Vector{Int}, subjects::Vector{Int}) = (
    Int[min(query, subject) for (query, subject) in zip(queries, subjects)],
    Int[max(query, subject) for (query, subject) in zip(queries, subjects)],
)

"""
Given each pair's ancestor name, build the `lca` relation over `n_genes` genes.
Returns `(matrix, labels)`, the matrix holding 1-based codes into `labels`.

An empty name is taken as unknown and left unstored, and a repeated pair keeps
the first name seen.
"""
function _lca_relation(lower, upper, names::Vector{String}, rows::Vector{Int}, n_genes)
    labels = String[]
    codes = Dict{String,UInt32}()
    cell_i, cell_j, cell_codes = Int[], Int[], UInt32[]
    for (entry, row) in enumerate(rows)
        isempty(names[row]) && continue
        push!(cell_i, lower[entry])
        push!(cell_j, upper[entry])
        push!(cell_codes, _intern!(labels, codes, names[row]))
    end
    # `min` keeps the first name seen, codes being assigned in order.
    return _relation(cell_i, cell_j, cell_codes, n_genes, min), labels
end

#= Gene lookup =#

"""Given a gene index, check it lies within `pg`. Returns it as an `Int`."""
function _gene_index(pg::ParalogGroup, index::Integer)
    n = _n_genes(pg)
    1 <= index <= n || throw(ArgumentError("gene index $index out of bounds (1:$n)"))
    return Int(index)
end

"""Given a gene ID, look up its linear index in `pg`. Returns it as an `Int`."""
function _gene_index(pg::ParalogGroup, id::AbstractString)
    haskey(pg.id_to_index, id) || throw(ArgumentError("Unknown gene ID: \"$id\""))
    return pg.id_to_index[id]
end

function Base.show(io::IO, pg::ParalogGroup)
    present = [String(field) for field in RELATIONS if !isnothing(getfield(pg, field))]
    isnothing(pg.lca) || push!(present, "lca")
    print(
        io,
        "ParalogGroup($(_plural(_n_genes(pg), "gene")), $(_plural(length(pg.scaffold_ranges), "scaffold")), $(_plural(nnz(pg.topology) ÷ 2, "pair")); relations: $(isempty(present) ? "none" : join(present, ", ")))",
    )
end

"""
    pg[i::Integer]
    pg[id::AbstractString]

Given a gene, gather its column from every relation matrix present on `pg`.
Returns a `Dict{Symbol, Vector}` keyed by `:topology` and the relations `pg`
holds, `:lca` surfacing ancestor names rather than codes.

**NOTE:** the symmetric relations store only their upper triangle, so a column
alone surfaces just the partners with a smaller index; and since a stored 0 is a
real value, `:topology` is the authority on which pairs exist.
"""
function Base.getindex(pg::ParalogGroup, i::Integer)
    gene = _gene_index(pg, i)
    columns = Dict{Symbol,Vector}(:topology => Vector(pg.topology[:, gene]))
    for field in RELATIONS
        matrix = getfield(pg, field)
        isnothing(matrix) || (columns[field] = Vector(matrix[:, gene]))
    end
    isnothing(pg.lca) || (
        columns[:lca] =
            [iszero(code) ? "" : pg.lca_labels[code] for code in pg.lca[:, gene]]
    )
    return columns
end

Base.getindex(pg::ParalogGroup, id::AbstractString) = pg[_gene_index(pg, id)]

"""
    pg[indices::AbstractVector{<:Integer}]
    pg[ids::AbstractVector{<:AbstractString}]

Given a set of genes, slice `pg` down to them. Returns a new `ParalogGroup` whose
relation matrices, `intervals` and `scaffold_ranges` are all concordant with the
new, contiguous linear indexing.

`lca_labels` is carried over whole, so a label the sub-group no longer uses
simply goes unreferenced.
"""
function Base.getindex(pg::ParalogGroup, indices::AbstractVector{<:Integer})
    isempty(indices) && throw(ArgumentError("`indices` must be non-empty"))
    order = sort!(unique(Int[_gene_index(pg, index) for index in indices]))
    order_set = Set(order)
    index_to_id = _index_to_id(pg)

    intervals, scaffold_ranges, id_to_index = _layout_genes(
        (name, pg.intervals[name][gene-first(range)+1], index_to_id[gene]) for
        (name, range) in pg.scaffold_ranges for gene in range if gene in order_set
    )

    return _rebuild(
        pg,
        function (field)
            matrix = getfield(pg, field)
            return isnothing(matrix) ? nothing : matrix[order, order]
        end;
        intervals = intervals,
        scaffold_ranges = scaffold_ranges,
        id_to_index = id_to_index,
        topology = pg.topology[order, order],
    )
end

function Base.getindex(pg::ParalogGroup, ids::AbstractVector{<:AbstractString})
    isempty(ids) && throw(ArgumentError("`ids` must be non-empty"))
    return pg[Int[_gene_index(pg, id) for id in ids]]
end

"""
    lca_label(pg::ParalogGroup, a, b) -> Union{Nothing, String}

Given two genes (linear indices or gene IDs), look up the last common ancestor
recorded for the pair. Returns its name, or `nothing` when `pg` has no `lca`
relation or the pair has none recorded. Argument order does not matter.
"""
function lca_label(pg::ParalogGroup, a::Integer, b::Integer)
    isnothing(pg.lca) && return nothing
    code = _tri(pg.lca, _gene_index(pg, a), _gene_index(pg, b))
    return iszero(code) ? nothing : pg.lca_labels[code]
end

lca_label(pg::ParalogGroup, a::AbstractString, b::AbstractString) =
    lca_label(pg, _gene_index(pg, a), _gene_index(pg, b))

#= Graph views =#

"""
    topology_graph(pg::ParalogGroup) -> SimpleGraph

Given a group, view its topology as an unweighted `Graphs.SimpleGraph`. Returns
that graph, in which vertex `v` is gene `v` (see `pg.id_to_index` for the name)
and edge `u-v` marks a topology pair.
"""
topology_graph(pg::ParalogGroup) = SimpleGraph(pg.topology)

"""
Given a symmetric relation, view it as a weighted graph. Returns that graph, or
`nothing` for an absent relation; `SimpleWeightedGraph` needs a fully symmetric
adjacency matrix, so the stored upper triangle is mirrored first.
"""
_weighted_undirected_graph(::Nothing) = nothing
_weighted_undirected_graph(m::SparseMatrixCSC) = SimpleWeightedGraph(m + permutedims(m))

"""
    dN_graph(pg::ParalogGroup) -> Union{Nothing, SimpleWeightedGraph}

Given a group, view its `dN` relation as an undirected
`SimpleWeightedGraphs.SimpleWeightedGraph`. Returns that graph, or `nothing`
when `dN` wasn't supplied.

**NOTE:** a weighted graph represents "no edge" as weight 0, and symmetrising
prunes zeros besides, so a pair whose dN is 0 has no edge here. Ask
`pg.topology` which pairs exist.
"""
dN_graph(pg::ParalogGroup) = _weighted_undirected_graph(pg.dN)

"""
    dS_graph(pg::ParalogGroup) -> Union{Nothing, SimpleWeightedGraph}

Same as [`dN_graph`](@ref), for `pg.dS`.
"""
dS_graph(pg::ParalogGroup) = _weighted_undirected_graph(pg.dS)

"""
Given a directed relation, view it as a weighted digraph. Returns that graph, or
`nothing` for an absent relation; the relation is already directed, so nothing
is symmetrised.
"""
_weighted_directed_graph(::Nothing) = nothing
_weighted_directed_graph(m::SparseMatrixCSC) = SimpleWeightedDiGraph(m)

"""
    id_subject_query_graph(pg::ParalogGroup) -> Union{Nothing, SimpleWeightedDiGraph}

Given a group, view its `id_subject_query` relation as a
`SimpleWeightedGraphs.SimpleWeightedDiGraph`. Returns that graph, or `nothing`
when the relation wasn't supplied. Edge `i -> j` carries weight `M[i, j]`, with
no symmetrisation, the relation being directed.
"""
id_subject_query_graph(pg::ParalogGroup) = _weighted_directed_graph(pg.id_subject_query)

"""
    id_query_subject_graph(pg::ParalogGroup) -> Union{Nothing, SimpleWeightedDiGraph}

Same as [`id_subject_query_graph`](@ref), for `pg.id_query_subject` (stored
`[subject, query]`).
"""
id_query_subject_graph(pg::ParalogGroup) = _weighted_directed_graph(pg.id_query_subject)

#= OrthoFinder duplication ingestion =#

"""
Given a Newick tree, read the depth of every labelled node, the root at 1.
Returns a `label => depth` map.

A label's depth is its parenthesis nesting level plus one, so an internal label
(`N3`) and a terminal one (a species name) are read the same way and no tree
needs building. Branch lengths are skipped.
"""
function _newick_depths(newick::AbstractString)
    depths = Dict{String,Int}()
    level = 0
    label_start = 0
    label_level = 0
    reading_length = false

    for (index, char) in pairs(newick)
        if char in ('(', ')', ',', ';', ':')
            if label_start > 0
                label = strip(newick[label_start:prevind(newick, index)])
                isempty(label) || (depths[String(label)] = label_level + 1)
                label_start = 0
            end
            reading_length = char == ':'
            char == '(' && (level += 1)
            char == ')' && (level -= 1)
        elseif !reading_length && label_start == 0 && !isspace(char)
            label_start = index
            label_level = level
        end
    end
    return depths
end

"""Species tree OrthoFinder writes alongside a `Duplications.tsv`."""
_species_tree_path(duplications::AbstractString) = joinpath(
    dirname(dirname(duplications)),
    "Species_Tree",
    "SpeciesTree_rooted_node_labels.txt",
)

"""
Given a `Genes 1`/`Genes 2` field, pick out the genes of one species. Returns
those carrying `prefix`, with it stripped and `transform` applied.

**NOTE:** the prefix cannot be found by splitting on the first `_`, since a
species name may hold one itself (`pongo_abelii_XP_054409609.1`).
"""
function _species_genes(
    field::AbstractString,
    prefix::AbstractString,
    transform::F,
) where {F}
    genes = String[]
    for name in eachsplit(field, ',')
        name = strip(name)
        startswith(name, prefix) || continue
        push!(genes, String(transform(chop(name; head = length(prefix), tail = 0))))
    end
    return genes
end

"""
Given a cell predicate, restrict `m`'s stored entries to the cells it admits.
Returns the restricted matrix, or `nothing` for a relation `pg` never held.

Done over `findnz` rather than by multiplying with a mask, which would prune the
explicit zeros a meaningful dN/dS of 0 relies on.
"""
function _mask_relation(m::SparseMatrixCSC, keep::F, n_genes::Integer) where {F}
    rows, cols, vals = findnz(m)
    selected = findall(entry -> keep(rows[entry], cols[entry]), eachindex(rows))
    return _relation(rows[selected], cols[selected], vals[selected], n_genes)
end

_mask_relation(::Nothing, keep, n_genes::Integer) = nothing

"""
Given a `Duplications.tsv` header, locate the columns this reader needs.
Returns them as a `NamedTuple` of indices.
"""
function _duplication_columns(header::Vector{<:AbstractString}, path::AbstractString)
    function column(name)
        index = findfirst(==(name), header)
        isnothing(index) && throw(
            ArgumentError("\"$path\" has no \"$name\" column; is it a Duplications.tsv?"),
        )
        return index
    end
    return (
        node = column("Species Tree Node"),
        support = column("Support"),
        first_genes = column("Genes 1"),
        second_genes = column("Genes 2"),
    )
end

"""
Given the two gene sets of one duplication, record it against every pair of `pg`
they span. Returns how many of those pairs were already recorded at a different
depth, the shallower of the two — the one that actually separates the pair —
being kept.
"""
function _record_duplication!(
    annotations::Dict{Tuple{Int,Int},Tuple{String,Int}},
    pg::ParalogGroup,
    left::Vector{Int},
    right::Vector{Int},
    node::String,
    depth::Int,
)
    n_conflicts = 0
    for a in left, b in right
        (a == b || !pg.topology[a, b]) && continue
        pair = minmax(a, b)
        recorded = get(annotations, pair, nothing)
        if isnothing(recorded)
            annotations[pair] = (node, depth)
        elseif recorded[2] != depth
            n_conflicts += 1
            recorded[2] > depth && (annotations[pair] = (node, depth))
        end
    end
    return n_conflicts
end

"""
Given an open `Duplications.tsv`, collect the duplication each of `pg`'s pairs
descends from. Returns `(annotations, n_conflicts)`, the annotations mapping an
upper-triangle cell to its `(node, depth)`.

Rows below `min_support`, or whose support cannot be read, are skipped, as are
nodes absent from `depths`.
"""
function _read_duplications(
    io::IO,
    path::AbstractString,
    pg::ParalogGroup,
    prefix::AbstractString,
    transform::F,
    depths::Dict{String,Int},
    min_support::Real,
) where {F}
    columns = _duplication_columns(split(readline(io), '\t'), path)
    n_columns = maximum(columns)
    annotations = Dict{Tuple{Int,Int},Tuple{String,Int}}()
    n_conflicts = 0

    indices(field) = [
        pg.id_to_index[gene] for
        gene in _species_genes(field, prefix, transform) if haskey(pg.id_to_index, gene)
    ]

    for line in eachline(io)
        fields = split(line, '\t')
        length(fields) >= n_columns || continue

        support = tryparse(Float64, strip(fields[columns.support]))
        (isnothing(support) || support < min_support) && continue

        node = String(strip(fields[columns.node]))
        depth = get(depths, node, 0)
        depth == 0 && continue

        left = indices(fields[columns.first_genes])
        isempty(left) && continue
        right = indices(fields[columns.second_genes])
        isempty(right) && continue

        n_conflicts += _record_duplication!(annotations, pg, left, right, node, depth)
    end
    return annotations, n_conflicts
end

"""
Given the collected duplications, check their nodes lie on one lineage. Returns
`nothing`, or throws when two nodes share a depth — which would break the
comparability `lca_depth` rests on.
"""
function _check_one_lineage(annotations, tree_path::AbstractString)
    depth_of_node = Dict{String,Int}()
    for (node, depth) in values(annotations)
        depth_of_node[node] = depth
    end
    allunique(values(depth_of_node)) || throw(
        ArgumentError(
            "Nodes $(join(sort!(collect(keys(depth_of_node))), ", ")) do not lie on one lineage of \"$tree_path\"",
        ),
    )
    return nothing
end

"""
Given the collected duplications, turn them into relation coordinates. Returns
`(edge_i, edge_j, lca_codes, depths, labels)`, ordered by pair so that the
interned label order does not follow `Dict` hash order.
"""
function _annotation_entries(annotations::Dict{Tuple{Int,Int},Tuple{String,Int}})
    edge_i, edge_j = Int[], Int[]
    lca_codes, depth_values = UInt32[], Float64[]
    labels = String[]
    codes = Dict{String,UInt32}()

    for pair in sort!(collect(keys(annotations)))
        node, depth = annotations[pair]
        push!(edge_i, pair[1])
        push!(edge_j, pair[2])
        push!(depth_values, Float64(depth))
        push!(lca_codes, _intern!(labels, codes, node))
    end
    return edge_i, edge_j, lca_codes, depth_values, labels
end

"""
Given what [`add_duplications_of`](@ref) dropped, report it. Returns `nothing`,
warning once for the pairs left unannotated and once for the pairs reported at
more than one node.
"""
function _warn_duplication_losses(n_pairs::Int, n_genes::Int, n_conflicts::Int)
    if n_pairs > 0
        orphans =
            n_genes > 0 ? ", leaving $(_plural(n_genes, "gene")) without a partner" : ""
        @warn "add_duplications_of: dropped $(_plural(n_pairs, "pair")) with no duplication recorded$orphans"
    end
    n_conflicts > 0 &&
        @warn "add_duplications_of: $(_plural(n_conflicts, "pair")) reported at more than one node; kept the shallowest"
    return nothing
end

"""
    add_duplications_of(pg, path, species; species_tree = nothing,
                        min_support = 0.0, transform = identity) -> ParalogGroup

Given an OrthoFinder `Duplications.tsv` at `path`, annotate `pg`'s pairs with the
duplication each descends from. Returns a **new** group; only `species`' own
duplications count, which is what puts every depth on one lineage. Depths come
from `species_tree`, by default the `SpeciesTree_rooted_node_labels.txt` beside
`path`; `transform` maps a gene name, species prefix stripped, to `pg`'s own ID.

**NOTE:** `lca`/`lca_depth` are overwritten, and a pair with no duplication
recorded is dropped, so the result is always rankable by `:lca_depth`.
"""
function add_duplications_of(
    pg::ParalogGroup,
    path::AbstractString,
    species::AbstractString;
    species_tree::Union{Nothing,AbstractString} = nothing,
    min_support::Real = 0.0,
    transform = identity,
)
    tree_path = something(species_tree, _species_tree_path(path))
    isfile(tree_path) ||
        throw(ArgumentError("No species tree at \"$tree_path\"; pass `species_tree`"))
    depths = _newick_depths(read(tree_path, String))
    haskey(depths, species) ||
        throw(ArgumentError("\"$species\" is not a node of \"$tree_path\""))

    annotations, n_conflicts = open(path) do io
        _read_duplications(io, path, pg, species * "_", transform, depths, min_support)
    end
    isempty(annotations) && throw(
        ArgumentError(
            "No pair of `pg` matched a duplication of \"$species\" in \"$path\"; check `species` and `transform`",
        ),
    )
    _check_one_lineage(annotations, tree_path)

    n_genes = _n_genes(pg)
    edge_i, edge_j, lca_codes, depth_values, labels = _annotation_entries(annotations)
    kept_cells = Set(zip(edge_i, edge_j))
    kept_genes = sort!(unique(vcat(edge_i, edge_j)))
    _warn_duplication_losses(
        nnz(pg.topology) ÷ 2 - length(edge_i),
        n_genes - length(kept_genes),
        n_conflicts,
    )

    lca_depth = _relation(edge_i, edge_j, depth_values, n_genes)
    lca = _relation(edge_i, edge_j, lca_codes, n_genes, min)
    annotated = _rebuild(
        pg,
        field ->
            field === :lca_depth ? lca_depth :
            field === :lca ? lca :
            _mask_relation(
                getfield(pg, field),
                (row, col) -> minmax(Int(row), Int(col)) in kept_cells,
                n_genes,
            );
        topology = _topology(edge_i, edge_j, n_genes),
        lca_labels = labels,
    )
    return annotated[kept_genes]
end

#= Reciprocal best hits =#

"""
Given a scoring alias, canonicalise it (`"maximum"` → `"max"`, `"avg"` →
`"mean"`, …). Returns the canonical name, or `nothing` when unrecognised; which
names a given `rbh` method accepts is up to that method.
"""
function _canonical_scoring(scoring::AbstractString)
    key = lowercase(scoring)
    key in ("max", "maximum") && return "max"
    key in ("min", "minimum") && return "min"
    key in ("mean", "avg", "average") && return "mean"
    key in ("double_max", "ds") && return key
    return nothing
end

"""
Given an undirected edge, work out which end was the original query. Returns
`(query, subject)`, read off whichever `id_*` relation is present — their
storage axis records it exactly — and falling back to the lower index as query.
"""
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

Given an edge, combine whichever of `id_subject_query`/`id_query_subject` are
present via `scoring` (`"mean"`, `"min"` or `"max"`). Returns the combined %ID,
or `nothing` when neither relation is present.
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

"""
Given a graph, group its edges by connected component. Returns
`(components, edges)` as parallel vectors, built in one pass over the edges
rather than one pass per component.
"""
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

"""
Given no explicit `levels`, pick the ranking levels implicitly. Returns
everything `pg` carries bar `:lca_depth`, which is opt-in.
"""
function _rank_levels(pg::ParalogGroup, levels::Nothing)
    available = filter(level -> _has_level(pg, level), [:dS, :dN, :identity])
    isempty(available) && throw(
        ArgumentError(
            "pg has no dN, dS, or %ID relation to rank edges by; pass `levels` to rank by another relation",
        ),
    )
    return available
end

"""
Given the caller's own `levels`, check them against `pg`. Returns them
unchanged; a level `pg` lacks is an error, never a silently skipped rank.
"""
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

"""
Given one ranking level, score edge `{a, b}` by it. Returns the score with
similarity levels negated, so that smaller means better throughout.
"""
function _level_value(pg::ParalogGroup, level::Symbol, a::Integer, b::Integer, scoring)
    level === :dS && return _tri(pg.dS, a, b)
    level === :dN && return _tri(pg.dN, a, b)
    level === :lca_depth && return -_tri(pg.lca_depth, a, b)
    return -edge_identity(pg, _edge_query_subject(pg, a, b)..., scoring)
end

"""First few of `indices`, so a warning naming a large component stays readable."""
_truncated(indices, limit = 4) =
    join(Iterators.take(indices, limit), ", ") * (length(indices) > limit ? ", …" : "")

"""
Given a group about to be ranked by `:lca_depth`, check every edge has one
recorded. Returns `nothing`; depths start at 1, so a 0 means "not recorded" and
is an error rather than a silent default.
"""
function _check_depths(pg::ParalogGroup, graph::SimpleGraph)
    index_to_id = _index_to_id(pg)
    unranked = [
        (index_to_id[src(edge)], index_to_id[dst(edge)]) for
        edge in edges(graph) if iszero(_tri(pg.lca_depth, src(edge), dst(edge)))
    ]
    isempty(unranked) || throw(
        ArgumentError(
            "rbh: ranking by `lca_depth`, but $(_plural(length(unranked), "pair")) have none recorded, e.g. $(first(unranked))",
        ),
    )
    return nothing
end

"""
Given one component's edges, pick the best-ranking one. Returns
`(query, subject, tied)`, `tied` marking a tie on *every* level — which the
lowest gene indices then settle.
"""
function _best_edge(pg::ParalogGroup, edge_list, ranking::Vector{Symbol}, scoring::String)
    # Ascending throughout (see `_level_value`), with `(lo, hi)` as a final,
    # deterministic tie-break: lowest indices win.
    rank_key(a, b) =
        ([_level_value(pg, level, a, b, scoring) for level in ranking], minmax(a, b))
    ranked = sort([(rank_key(src(e), dst(e)), src(e), dst(e)) for e in edge_list])

    best_key = ranked[1][1][1]
    tied = count(entry -> entry[1][1] == best_key, ranked) > 1
    _, a, b = ranked[1]
    return (_edge_query_subject(pg, a, b)..., tied)
end

"""
Given the pairs [`rbh`](@ref) chose, assemble its output table. Returns a
`DataFrame` shaped like `ParalogGroup`'s constructor input, carrying whichever
relations `pg` holds.
"""
function _pair_table(pg::ParalogGroup, chosen::Vector{Tuple{Int,Int}})
    index_to_id = _index_to_id(pg)
    columns = Pair{String,Vector}[
        "query"=>String[index_to_id[query] for (query, _) in chosen],
        "subject"=>String[index_to_id[subject] for (_, subject) in chosen],
    ]
    add(name, values) = push!(columns, name => values)

    isnothing(pg.dN) ||
        add("dN", Float64[_tri(pg.dN, query, subject) for (query, subject) in chosen])
    isnothing(pg.dS) ||
        add("dS", Float64[_tri(pg.dS, query, subject) for (query, subject) in chosen])
    isnothing(pg.id_subject_query) || add(
        "id_subject_query",
        Float64[pg.id_subject_query[query, subject] for (query, subject) in chosen],
    )
    isnothing(pg.id_query_subject) || add(
        "id_query_subject",
        Float64[pg.id_query_subject[subject, query] for (query, subject) in chosen],
    )
    isnothing(pg.lca_depth) || add(
        "lca_depth",
        Float64[_tri(pg.lca_depth, query, subject) for (query, subject) in chosen],
    )
    isnothing(pg.lca) || add(
        "lca",
        String[something(lca_label(pg, query, subject), "") for (query, subject) in chosen],
    )
    return DataFrame(columns...)
end

"""
    rbh(pg::ParalogGroup; scoring = "mean", levels = nothing) -> DataFrame

Given a group, keep the best-ranking edge of each weakly-connected component of
`pg.topology`. Returns those pairs as a `DataFrame` shaped like the constructor's
input; singleton (edge-less) components are skipped.

`levels` ranks edges, most significant first: `:dS`/`:dN` ascending, `:identity`
(combined via `scoring`) and `:lca_depth` descending. It defaults to whichever of
`:dS`, `:dN`, `:identity` `pg` carries — `:lca_depth` is never implicit, and a
level `pg` lacks is an error rather than a skipped rank. A tie on *every* level
goes to the lowest gene indices, and warns.
"""
function rbh(pg::ParalogGroup; scoring::String = "mean", levels = nothing)
    scoring = something(_canonical_scoring(scoring), "")
    scoring in ("mean", "min", "max") ||
        throw(ArgumentError("scoring must be \"mean\", \"min\", or \"max\""))
    ranking = _rank_levels(pg, levels)

    graph = topology_graph(pg)
    :lca_depth in ranking && _check_depths(pg, graph)

    chosen = Tuple{Int,Int}[]
    tied_components = Vector{Int}[]
    for (component, edge_list) in zip(_component_edges(graph)...)
        isempty(edge_list) && continue   # singleton gene, no partner to hit
        query, subject, tied = _best_edge(pg, edge_list, ranking, scoring)
        push!(chosen, (query, subject))
        tied && push!(tied_components, sort(component))
    end

    if !isempty(tied_components)
        groups = join(("[$(_truncated(c))]" for c in first(tied_components, 3)), ", ")
        @warn "rbh: tied best-hit rank in $(_plural(length(tied_components), "group")): $groups$(length(tied_components) > 3 ? ", …" : "")"
    end
    return _pair_table(pg, chosen)
end

#= Reciprocal best hits over a flat pair table =#

"""
Given a flat paralog table, validate it: two ID columns followed by
`n_columns - 2` numeric ones. Returns `nothing`.
"""
function _check_pair_columns(paralog_df::DataFrame, n_columns::Int)
    ncol(paralog_df) >= n_columns || throw(
        ArgumentError(
            "`paralog_df` needs at least $n_columns columns, got $(ncol(paralog_df))",
        ),
    )
    for column = 1:2
        _require_eltype(paralog_df, column, AbstractString, "gene IDs (strings)")
    end
    for column = 3:n_columns
        _require_eltype(paralog_df, column, Real, "scores (numbers)")
    end
    return nothing
end

"""
Given a flat paralog table, collect its gene IDs. Returns `(ids, id_to_index)`,
the IDs in first-appearance order and the reverse lookup to their positions.
"""
function _index_pair_ids(paralog_df::DataFrame)
    ids = unique(vcat(string.(paralog_df[:, 1]), string.(paralog_df[:, 2])))
    return ids, Dict(id => index for (index, id) in enumerate(ids))
end

"""
Given a flat paralog table, score its pairs. Returns `(original, ranked)`, both
`(gene, partner) => value` dictionaries: `original` holds the two directional
values as given, `ranked` the values best hits are chosen by.

A later row overwrites an earlier one for a repeated pair.
"""
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

"""
Given a sparse score matrix held as `(row, column) => score`, find each gene's
best partner. Returns `(best_in_row, best_in_column)`, gene indices with `0`
where a gene has no scored partner.

`better(a, b)` decides whether `a` beats `b` (`>` for a similarity, `<` for a
distance); ties go to the lowest index, so the result is iteration-order
independent.
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

"""
Given the ranked scores, find the genes that are each other's best hit. Returns
those `(gene, best_partner)` pairs in ascending gene order, each gene appearing
at most once.
"""
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
Given a flat pair table's score dictionaries, resolve its reciprocal best hits
back to names and scores. Returns `(gene_ids, paralog_ids, forwards, backwards)`
as parallel vectors, the scores being the two directions as originally given.
"""
function _reciprocal_scores(
    original::Dict{Tuple{Int,Int},Float64},
    ranked::Dict{Tuple{Int,Int},Float64},
    ids::Vector{String},
    better::F,
) where {F}
    gene_ids, paralog_ids = String[], String[]
    forwards, backwards = Float64[], Float64[]

    for (gene, partner) in _reciprocal_pairs(ranked, length(ids), better)
        push!(gene_ids, ids[partner])
        push!(paralog_ids, ids[gene])
        push!(forwards, get(original, (gene, partner), 0.0))
        push!(backwards, get(original, (partner, gene), 0.0))
    end
    return gene_ids, paralog_ids, forwards, backwards
end

"""
    rbh(paralog_df::DataFrame; scoring = "max") -> DataFrame

Given paralogs one pair per row — two ID columns, then percent identity
gene → paralog and paralog → gene — keep the pairs that are each other's best
hit. Returns a `DataFrame` of `GeneID`, `ParalogID`, `perc_1`, `perc_2`,
`max_perc` and `mean_perc`. `scoring` picks what pairs are ranked by: `"max"`
(the default), `"mean"`, `"double_max"` (each direction on its own), or `"ds"`,
which reads column 3 as a dS and defers to [`rbh_ds`](@ref).

**NOTE:** only scored pairs are ranked. A gene absent from a pair has no best hit
there, rather than scoring zero against it.
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
    gene_ids, paralog_ids, forwards, backwards =
        _reciprocal_scores(original, ranked, ids, >)

    return DataFrame(
        "GeneID" => gene_ids,
        "ParalogID" => paralog_ids,
        "perc_1" => forwards,
        "perc_2" => backwards,
        "max_perc" => max.(forwards, backwards),
        "mean_perc" => (forwards .+ backwards) ./ 2,
    )
end

"""
    rbh_ds(paralog_df::DataFrame) -> DataFrame

Given paralogs listed one pair per row, with column 3 holding the pair's dS,
keep the pairs that are each other's best hit. Returns a `DataFrame` of
`GeneID`, `ParalogID`, `ds` and `min_ds`.

dS is a distance, so the **lowest** value wins; see [`rbh`](@ref) for the
pairing rule.
"""
function rbh_ds(paralog_df::DataFrame)
    _check_pair_columns(paralog_df, 3)
    ids, id_to_index = _index_pair_ids(paralog_df)
    original, ranked = _pair_matrices(paralog_df, id_to_index, "ds")
    gene_ids, paralog_ids, forwards, backwards =
        _reciprocal_scores(original, ranked, ids, <)

    return DataFrame(
        "GeneID" => gene_ids,
        "ParalogID" => paralog_ids,
        "ds" => forwards,
        "min_ds" => min.(forwards, backwards),
    )
end

export ParalogGroup,
    add_duplications_of,
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

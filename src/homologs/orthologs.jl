"""
Homology *between* genomes: a [`Clade`](@ref) of species, its species tree and
every OrthoFinder gene tree, laid out for fast ortholog and common-ancestor
queries. `docs/clade.md` is the reference for the layout and the rules.

See [`Paralogs`](@ref) for the within-genome counterpart.
"""
module Orthologs

using DataFrames

# Nested one level deeper than the other modules, so `Reference` is three dots up.
using ...Reference
using ..Paralogs: _newick_tokens

"""The most species a [`Clade`](@ref) holds: one bit each in a `UInt64` mask."""
const MAX_SPECIES = 64

"""A tree read from `Resolved_Gene_Trees.txt`."""
const TREE_RESOLVED = 0x01

"""A tree-less orthogroup of three or more genes, read as a star (one polytomy)."""
const TREE_STAR = 0x02

"""
A tree-less orthogroup of two genes: a star whose genes are never called
orthologs, as OrthoFinder calls none (see `docs/clade.md`).
"""
const TREE_PAIR = 0x03

#= Types =#

"""
A species tree in preorder: `parent` (0 at the root), `subtree_end` (a node's
last descendant), `branch_length` (`NaN` when absent), `labels`, `depth` (root 1,
as `lca_depth` counts) and a `species` bitmask per node, bit `s - 1` standing for
the `s`-th leaf. Internal nodes without a name, or named by a support value, are
named `N0`, `N1`, … in preorder, as OrthoFinder names them.
"""
struct SpeciesTree
    labels::Vector{String}
    parent::Vector{Int32}
    subtree_end::Vector{Int32}
    branch_length::Vector{Float32}
    depth::Vector{Int32}
    species::Vector{UInt64}
end

"""
Every gene tree of a [`Clade`](@ref), concatenated in preorder: per node its
`parent` (0 at a root), `subtree_end`, `species` bitmask, `branch_length`, and,
for a leaf, the metadata index of its `gene` (0 when unplaced) and its `protein`.
Per tree, `tree_start` offsets, sorted `orthogroups` names and a `source`
(`TREE_RESOLVED`/`TREE_STAR`/`TREE_PAIR`). Layout in `docs/clade.md`.
"""
struct GeneForest
    parent::Vector{Int32}
    subtree_end::Vector{Int32}
    species::Vector{UInt64}
    branch_length::Vector{Float32}
    gene::Vector{UInt32}
    protein::Vector{String}
    tree_start::Vector{Int32}
    orthogroups::Vector{String}
    source::Vector{UInt8}
end

"""Create an empty `GeneForest`, ready for trees to be appended."""
GeneForest() = GeneForest(
    Int32[],
    Int32[],
    UInt64[],
    Float32[],
    UInt32[],
    String[],
    Int32[1],
    String[],
    UInt8[],
)

"""
One species' lookup from an ID to its gene-tree leaf: sorted `keys`, with the
leaf `nodes` alongside. Binary searched rather than hashed, the way `Reference`
resolves IDs.
"""
struct GeneIndex
    keys::Vector{String}
    nodes::Vector{Int32}
end

"""
Several [`Species`](@ref), ordered as the leaves of their `species_tree`, with
every gene tree relating their genes (`gene_trees`) and, per species, a gene-ID
and a protein-accession [`GeneIndex`](@ref). A species with no genome loaded is
reachable by protein accession only. Build with [`load_clade`](@ref).
"""
struct Clade
    species::Vector{Species}
    species_tree::SpeciesTree
    gene_trees::GeneForest
    genes::Vector{GeneIndex}
    proteins::Vector{GeneIndex}
end

"""
One orthogroup's tree inside a [`Clade`](@ref), as a view: `orthogroup` is its
position in the forest. Nodes are the forest's preorder indices
([`tree_nodes`](@ref)), which every accessor takes and returns.
"""
struct GeneTree
    clade::Clade
    orthogroup::Int32
end

#= Newick =#

"""
Given a Newick tree, lay it out in preorder. Returns `(parent, subtree_end,
branch_length, labels)`, parents local to this tree and 0 at the root; an
unnamed internal node's label is `""`. Throws on unbalanced parentheses or an
unlabelled leaf.
"""
function _preorder_layout(newick::AbstractString)
    parent, subtree_end, branch_length, labels = Int32[], Int32[], Float32[], String[]
    open_nodes = Int32[]
    current = Int32(0)
    previous = nothing
    for token in _newick_tokens(newick; lengths = true)
        if token === '(' || (token isa String && previous !== ')')
            token isa String &&
                isempty(open_nodes) &&
                throw(ArgumentError("Newick leaf `$token` outside any clade"))
            push!(parent, isempty(open_nodes) ? Int32(0) : last(open_nodes))
            push!(subtree_end, Int32(length(parent)))
            push!(branch_length, NaN32)
            push!(labels, token === '(' ? "" : token)
            current = Int32(length(parent))
            token === '(' && push!(open_nodes, current)
        elseif token === ')'
            isempty(open_nodes) && throw(ArgumentError("unbalanced Newick: `$newick`"))
            current = pop!(open_nodes)
            subtree_end[current] = Int32(length(parent))
        elseif token isa Float64
            previous isa String ||
                previous === ')' ||
                throw(
                    ArgumentError("Newick branch length on an unlabelled leaf: `$newick`"),
                )
            branch_length[current] = token
        else
            labels[current] = token
        end
        previous = token
    end
    isempty(open_nodes) || throw(ArgumentError("unbalanced Newick: `$newick`"))
    return (; parent, subtree_end, branch_length, labels)
end

"""Given a node label, tell whether it is empty or a support value rather than a name."""
_unnamed(label::AbstractString) = isempty(label) || !isnothing(tryparse(Float64, label))

"""
Given a Newick species tree, lay it out as a [`SpeciesTree`](@ref). Returns it;
throws on an unnamed leaf, a repeated name, or more than [`MAX_SPECIES`](@ref) leaves.
"""
function SpeciesTree(newick::AbstractString)
    layout = _preorder_layout(newick)
    n_nodes = length(layout.parent)
    labels, depth, species =
        copy(layout.labels), zeros(Int32, n_nodes), zeros(UInt64, n_nodes)
    n_internal, n_leaves = 0, 0
    for node = 1:n_nodes
        parent = layout.parent[node]
        depth[node] = parent == 0 ? 1 : depth[parent] + 1
        if layout.subtree_end[node] > node
            _unnamed(labels[node]) && (labels[node] = "N$n_internal")
            n_internal += 1
        else
            isempty(labels[node]) &&
                throw(ArgumentError("species tree has an unnamed leaf"))
            n_leaves += 1
            n_leaves <= MAX_SPECIES ||
                throw(ArgumentError("a Clade holds at most $MAX_SPECIES species"))
            species[node] = one(UInt64) << (n_leaves - 1)
        end
    end
    allunique(labels) || throw(ArgumentError("species tree repeats a node name"))
    for node = n_nodes:-1:1
        parent = layout.parent[node]
        parent == 0 || (species[parent] |= species[node])
    end
    return SpeciesTree(
        labels,
        layout.parent,
        layout.subtree_end,
        layout.branch_length,
        depth,
        species,
    )
end

"""Given a species tree, its leaf names in preorder: the species, in `Clade` order."""
species_names(tree::SpeciesTree) =
    [tree.labels[node] for node in eachindex(tree.labels) if tree.subtree_end[node] == node]

"""
Given a species tree and a set of species, find the deepest node holding them
all, as OrthoFinder places a gene-tree node. Returns its index.
"""
function _species_node(tree::SpeciesTree, mask::UInt64)
    best = 1
    for node in eachindex(tree.species)
        tree.species[node] & mask == mask &&
            tree.depth[node] > tree.depth[best] &&
            (best = node)
    end
    return best
end

#= Building the forest =#

"""
Given species names, the prefixes OrthoFinder gives their leaves, longest first
so that `pongo_abelii_` is tried before a `pongo_`. Returns `prefix => index` pairs.
"""
_leaf_prefixes(names::Vector{String}) = sort!(
    [name * "_" => index for (index, name) in enumerate(names)];
    by = pair -> -length(first(pair)),
)

"""
Given an OrthoFinder leaf label `<species>_<protein>`, split it on the first of
`prefixes` it starts with. Returns `(species_index, protein)`; throws when no
species prefixes it.
"""
function _split_leaf(label::AbstractString, prefixes::Vector{Pair{String,Int}})
    for (prefix, index) in prefixes
        startswith(label, prefix) || continue
        return index, String(chop(label; head = length(prefix), tail = 0))
    end
    throw(
        ArgumentError("gene-tree leaf `$label` starts with no species of the species tree"),
    )
end

"""Given a forest, append one node. Returns its index."""
function _push_node!(
    forest::GeneForest,
    parent::Int32,
    subtree_end::Int32,
    mask::UInt64,
    branch_length::Float32,
    protein::String,
)
    push!(forest.parent, parent)
    push!(forest.subtree_end, subtree_end)
    push!(forest.species, mask)
    push!(forest.branch_length, branch_length)
    push!(forest.gene, zero(UInt32))
    push!(forest.protein, protein)
    return Int32(length(forest.parent))
end

"""Given a forest whose last tree is complete, record that tree. Returns `nothing`."""
function _close_tree!(forest::GeneForest, orthogroup::String, source::UInt8)
    push!(forest.tree_start, Int32(length(forest.parent) + 1))
    push!(forest.orthogroups, orthogroup)
    push!(forest.source, source)
    return nothing
end

"""
Given a resolved gene tree's Newick, append it to `forest`. Returns `nothing`;
throws on a leaf no species prefixes, or an internal label other than
OrthoFinder's preorder `n0`, `n1`, … (support values are tolerated).
"""
function _append_resolved!(forest, newick, prefixes, orthogroup::String)
    layout = _preorder_layout(newick)
    offset = Int32(length(forest.parent))
    n_internal = 0
    for node in eachindex(layout.parent)
        label = layout.labels[node]
        if layout.subtree_end[node] > node
            _unnamed(label) ||
                label == "n$n_internal" ||
                throw(
                    ArgumentError(
                        "$orthogroup: internal node `$label` is not OrthoFinder's `n$n_internal`",
                    ),
                )
            n_internal += 1
            mask, protein = zero(UInt64), ""
        else
            species, protein = _split_leaf(label, prefixes)
            mask = one(UInt64) << (species - 1)
        end
        local_parent = layout.parent[node]
        _push_node!(
            forest,
            local_parent == 0 ? Int32(0) : local_parent + offset,
            layout.subtree_end[node] + offset,
            mask,
            layout.branch_length[node],
            protein,
        )
    end
    return _close_tree!(forest, orthogroup, TREE_RESOLVED)
end

"""
Given a tree-less orthogroup's `(species, protein)` members, append it to
`forest` as a star: one root over every gene. Returns `nothing`.
"""
function _append_star!(
    forest::GeneForest,
    members::Vector{Tuple{Int,String}},
    orthogroup::String,
)
    root = Int32(length(forest.parent) + 1)
    _push_node!(forest, Int32(0), root + Int32(length(members)), zero(UInt64), NaN32, "")
    for (species, protein) in members
        leaf = Int32(length(forest.parent) + 1)
        _push_node!(forest, root, leaf, one(UInt64) << (species - 1), NaN32, protein)
    end
    return _close_tree!(forest, orthogroup, length(members) == 2 ? TREE_PAIR : TREE_STAR)
end

"""
Given a forest, fill every internal node's species mask from its leaves.
Returns `forest`.
"""
function _fill_masks!(forest::GeneForest)
    # Preorder puts every parent before its children, so one reverse sweep suffices.
    for node in reverse(eachindex(forest.parent))
        parent = forest.parent[node]
        parent == 0 || (forest.species[parent] |= forest.species[node])
    end
    return forest
end

"""
Given `Resolved_Gene_Trees.txt` (lines `<orthogroup>: <newick>`), read every
tree. Returns an `orthogroup => newick` `Dict`; throws on a malformed line.
"""
function _read_gene_trees(path::AbstractString)
    trees = Dict{String,String}()
    for (line_number, line) in enumerate(eachline(path))
        isempty(strip(line)) && continue
        parts = split(line, ": "; limit = 2)
        length(parts) == 2 ||
            throw(ArgumentError("$path:$line_number is not `<orthogroup>: <newick>`"))
        trees[String(parts[1])] = String(parts[2])
    end
    return trees
end

"""
Given `Orthogroups.tsv` (one column of `, `-separated proteins per species),
read the orthogroups not in `skip`. Returns `orthogroup => [(species, protein)]`,
species indexed into `names`; throws on a column naming no species of the tree.
"""
function _read_orthogroups(path::AbstractString, names::Vector{String}, skip)
    members = Dict{String,Vector{Tuple{Int,String}}}()
    open(path) do io
        header = split(readline(io), '\t')
        columns = map(header[2:end]) do name
            index = findfirst(==(name), names)
            isnothing(index) && throw(
                ArgumentError("$path: column `$name` is no species of the species tree"),
            )
            index
        end
        for line in eachline(io)
            fields = split(line, '\t')
            orthogroup = String(fields[1])
            orthogroup in skip && continue
            genes = Tuple{Int,String}[]
            for (column, species) in enumerate(columns),
                protein in eachsplit(get(fields, column + 1, ""), ',')

                protein = strip(protein)
                isempty(protein) || push!(genes, (species, String(protein)))
            end
            isempty(genes) || (members[orthogroup] = genes)
        end
    end
    return members
end

"""
Given the resolved trees and the tree-less orthogroups' members, lay every
orthogroup out in one forest, sorted by name. Returns the [`GeneForest`](@ref).
"""
function _build_forest(trees::Dict{String,String}, members, prefixes)
    forest = GeneForest()
    for orthogroup in sort!(collect(union(keys(trees), keys(members))))
        newick = get(trees, orthogroup, nothing)
        if isnothing(newick)
            _append_star!(forest, members[orthogroup], orthogroup)
        else
            _append_resolved!(forest, newick, prefixes, orthogroup)
        end
    end
    length(forest.parent) < typemax(Int32) ||
        throw(ArgumentError("too many gene-tree nodes"))
    return _fill_masks!(forest)
end

#= Placing leaves on genes =#

"""Given a species position, its bit in a species mask."""
_bit(position::Integer) = one(UInt64) << (position - 1)

"""Given a forest node, tell whether it is a leaf."""
_is_leaf(forest::GeneForest, node::Integer) = forest.subtree_end[node] == node

"""Given a forest leaf, its species' position in the clade."""
_leaf_position(forest::GeneForest, node::Integer) = trailing_zeros(forest.species[node]) + 1

"""
Given a forest and its species count, partition its leaves by species position,
in preorder. Returns a vector of leaf vectors, one per species.
"""
function _species_leaves(forest::GeneForest, n_species::Int)
    by_species = [Int32[] for _ = 1:n_species]
    for node in eachindex(forest.parent)
        _is_leaf(forest, node) || continue
        push!(by_species[_leaf_position(forest, node)], Int32(node))
    end
    return by_species
end

"""
Given a genome and gene IDs, resolve each ID's metadata index in one walk,
preferring a `gene` feature where several share the ID. Returns `ID => index`.
"""
function _gene_indices(genome::Genome, ids::Vector{String})
    records = genome[ids]
    indices = Dict{String,UInt32}()
    for record in records
        record.feature_type == :gene || continue
        haskey(indices, record.id) || (indices[record.id] = parse_index(record.code))
    end
    for record in records
        haskey(indices, record.id) || (indices[record.id] = parse_index(record.code))
    end
    return indices
end

"""
Given a forest and one species' `leaves`, point each at the gene its protein
codes for, through `proteins` (a [`protein_gene_map`](@ref)). Returns `nothing`;
leaves it cannot place keep gene 0, and are counted in a warning.
"""
function _place_genes!(
    forest::GeneForest,
    leaves::Vector{Int32},
    member::Species,
    proteins::AbstractDict,
)
    gene_ids = [get(proteins, forest.protein[leaf], nothing) for leaf in leaves]
    indices = _gene_indices(
        member.genome,
        unique!(String[id for id in gene_ids if !isnothing(id)]),
    )
    n_unmapped, n_unloaded = 0, 0
    for (leaf, gene_id) in zip(leaves, gene_ids)
        index = isnothing(gene_id) ? nothing : get(indices, gene_id, nothing)
        if isnothing(index)
            isnothing(gene_id) ? (n_unmapped += 1) : (n_unloaded += 1)
        else
            forest.gene[leaf] = index
        end
    end
    n_unmapped + n_unloaded == 0 || @warn string(
        "$(member.name): $(n_unmapped + n_unloaded) of $(length(leaves)) gene-tree proteins ",
        "have no gene ($n_unmapped not in its protein map, $n_unloaded whose gene is not loaded)",
    )
    return nothing
end

"""
Given one species' `leaves`, index them by `key_of(leaf)`, skipping a leaf whose
key is `nothing`. Returns a [`GeneIndex`](@ref); a repeated key warns, its first
leaf (in orthogroup order) answering lookups.
"""
function _build_index(leaves::Vector{Int32}, key_of, what::String, name::String)
    keys, nodes = String[], Int32[]
    for leaf in leaves
        key = key_of(leaf)
        isnothing(key) && continue
        push!(keys, key)
        push!(nodes, leaf)
    end
    order = sortperm(keys)   # stable, so repeats keep their forest order
    keys, nodes = keys[order], nodes[order]
    n_repeats = count(index -> keys[index] == keys[index-1], 2:length(keys); init = 0)
    n_repeats == 0 ||
        @warn "$name: $n_repeats $what IDs label several gene-tree leaves; the first answers lookups"
    return GeneIndex(keys, nodes)
end

"""
Given the species tree's names and the loaded `species`, line them up in tree
order, a species without one getting an empty genome. Returns that vector;
throws on a `Species` the tree lacks, or two sharing a name.
"""
function _clade_members(species::AbstractVector{Species}, names::Vector{String})
    by_name = Dict(member.name => member for member in species)
    length(by_name) == length(species) || throw(ArgumentError("two Species share a name"))
    unknown = sort!(collect(setdiff(keys(by_name), names)))
    isempty(unknown) || throw(
        ArgumentError(
            "not in the species tree: $(join(unknown, ", ")); it holds $(join(names, ", "))",
        ),
    )
    return [get(() -> Species(name), by_name, name) for name in names]
end

"""Given an OrthoFinder results directory, the files a `Clade` reads from it by default."""
_orthofinder_files(directory::AbstractString) = (
    species_tree = joinpath(
        directory,
        "Species_Tree",
        "SpeciesTree_rooted_node_labels.txt",
    ),
    gene_trees = joinpath(directory, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt"),
    orthogroups = joinpath(directory, "Orthogroups", "Orthogroups.tsv"),
)

"""
    Clade(species, orthofinder_dir; protein_maps, species_tree, gene_trees, orthogroups)

Given loaded `species` and an OrthoFinder results directory, build their
[`Clade`](@ref). Returns it. `protein_maps` holds a `species name =>`
[`protein_gene_map`](@ref) per genome; the three paths default to OrthoFinder's
layout, and `orthogroups = nothing` skips the tree-less orthogroups.
"""
function Clade(
    species::AbstractVector{Species},
    orthofinder_dir::AbstractString;
    protein_maps::AbstractDict = Dict{String,Dict{String,String}}(),
    species_tree::AbstractString = _orthofinder_files(orthofinder_dir).species_tree,
    gene_trees::AbstractString = _orthofinder_files(orthofinder_dir).gene_trees,
    orthogroups::Union{Nothing,AbstractString} = _orthofinder_files(orthofinder_dir).orthogroups,
)
    tree = SpeciesTree(read(species_tree, String))
    names = species_names(tree)
    members = _clade_members(species, names)
    trees = _read_gene_trees(gene_trees)
    tree_less =
        isnothing(orthogroups) ? Dict{String,Vector{Tuple{Int,String}}}() :
        _read_orthogroups(orthogroups, names, keys(trees))
    forest = _build_forest(trees, tree_less, _leaf_prefixes(names))

    leaves_by_species = _species_leaves(forest, length(members))
    for (leaves, member) in zip(leaves_by_species, members)
        proteins = get(protein_maps, member.name, nothing)
        if !isnothing(proteins)
            _place_genes!(forest, leaves, member, proteins)
        elseif !isempty(member.genome.scaffolds)
            @warn "$(member.name) has a genome but no protein map; its leaves have no gene"
        end
    end
    genes = [
        _build_index(
            leaves,
            leaf -> _leaf_gene_id(forest, member, leaf),
            "gene",
            member.name,
        ) for (leaves, member) in zip(leaves_by_species, members)
    ]
    proteins = [
        _build_index(leaves, leaf -> forest.protein[leaf], "protein", member.name) for
        (leaves, member) in zip(leaves_by_species, members)
    ]
    return Clade(members, tree, forest, genes, proteins)
end

"""Given a leaf and its species, its gene ID. Returns it, or `nothing` when unplaced."""
function _leaf_gene_id(forest::GeneForest, member::Species, leaf::Integer)
    index = forest.gene[leaf]
    return index == 0 ? nothing : get_metadata_id(member.genome, index)
end

"""
    load_clade(orthofinder_dir, gffs; load_options = (feature_types = :gene,), paths...)

Given an OrthoFinder results directory and a `species name => GFF3 path` map,
load each genome ([`add_features!`](@ref) with `load_options`) and its
[`protein_gene_map`](@ref), species concurrently, then build the [`Clade`](@ref).
Returns it; `paths` pass to `Clade`. A species of the tree missing from `gffs`
gets an empty genome; a name the tree lacks throws before anything loads.
"""
function load_clade(
    orthofinder_dir::AbstractString,
    gffs::AbstractDict{<:AbstractString,<:AbstractString};
    load_options::NamedTuple = (feature_types = :gene,),
    paths...,
)
    tree_path = get(paths, :species_tree, _orthofinder_files(orthofinder_dir).species_tree)
    _clade_members(
        [Species(String(name)) for name in keys(gffs)],
        species_names(SpeciesTree(read(tree_path, String))),
    )

    sanitize_ids = get(load_options, :sanitize_ids, true)
    tasks = map(collect(gffs)) do (name, path)
        Threads.@spawn begin
            gff_path = String(path)
            proteins_task = Threads.@spawn protein_gene_map(gff_path; sanitize_ids)
            member = Species(String(name))
            add_features!(gff_path, member; load_options...)
            member => fetch(proteins_task)
        end
    end
    species, protein_maps = Species[], Dict{String,Dict{String,String}}()
    for task in tasks
        member, proteins = fetch(task)
        push!(species, member)
        protein_maps[member.name] = proteins
    end
    return Clade(species, orthofinder_dir; protein_maps, paths...)
end

#= Queries =#

"""Given a clade, its species names, in species-tree order."""
species_names(clade::Clade) = [member.name for member in clade.species]

"""
Given a clade and a species name, find its position. Returns that index; throws
an `ArgumentError` listing the clade's species when it holds no such species.
"""
function _species_position(clade::Clade, name::AbstractString)
    position = findfirst(member -> member.name == name, clade.species)
    isnothing(position) && throw(
        ArgumentError(
            "no species `$name` in the clade; it holds $(join(species_names(clade), ", "))",
        ),
    )
    return position
end

"""Given an index and a key, binary search it. Returns the key's leaf, or 0."""
function _lookup(index::GeneIndex, key::AbstractString)
    position = searchsortedfirst(index.keys, key)
    found = position <= length(index.keys) && index.keys[position] == key
    return found ? index.nodes[position] : Int32(0)
end

"""
Given a clade, a species position and an ID — a gene ID, else a protein
accession — find its leaf. Returns the node, or 0 when it is in no orthogroup.
"""
function _leaf(clade::Clade, position::Int, id::AbstractString)
    node = _lookup(clade.genes[position], id)
    return node == 0 ? _lookup(clade.proteins[position], id) : node
end

"""
The children of one forest node, iterated without allocating: `first` and then
each next sibling, `subtree_end[child] + 1`, up to the node's `subtree_end`.
"""
struct _Children
    subtree_end::Vector{Int32}
    node::Int32
    first::Int32
end

_Children(subtree_end::Vector{Int32}, node::Integer) =
    _Children(subtree_end, Int32(node), Int32(node) + Int32(1))

Base.iterate(children::_Children, child::Int32 = children.first) =
    child > children.subtree_end[children.node] ? nothing :
    (child, children.subtree_end[child] + Int32(1))
Base.IteratorSize(::Type{_Children}) = Base.SizeUnknown()
Base.eltype(::Type{_Children}) = Int32

"""Given a forest node, the tree it belongs to. Returns that tree's index."""
_tree_of(forest::GeneForest, node::Integer) = searchsortedlast(forest.tree_start, node)

"""Given a forest tree, its preorder node range."""
_tree_nodes(forest::GeneForest, tree::Integer) =
    forest.tree_start[tree]:(forest.tree_start[tree+1]-Int32(1))

"""Given a forest node, its leaves of the species in `mask`, lazily."""
_leaves_in(forest::GeneForest, node::Integer, mask::UInt64) = (
    leaf for leaf = node:forest.subtree_end[node] if
    forest.subtree_end[leaf] == leaf && forest.species[leaf] & mask != 0
)

"""
Given a leaf, collect its orthologs among the species in `wanted`: walking up,
every sibling sharing no species with the side holding the leaf. Returns those
leaves, nearest first; none in a [`TREE_PAIR`](@ref) orthogroup.
"""
function _ortholog_leaves(forest::GeneForest, leaf::Int32, wanted::UInt64)
    found = Int32[]
    forest.source[_tree_of(forest, leaf)] == TREE_PAIR && return found
    side, node = leaf, forest.parent[leaf]
    while node != 0
        if forest.species[node] & wanted != 0
            for child in _Children(forest.subtree_end, node)
                child == side && continue
                (
                    forest.species[child] & wanted != 0 &&
                    forest.species[child] & forest.species[side] == 0
                ) || continue
                append!(found, _leaves_in(forest, child, wanted))
            end
        end
        side, node = node, forest.parent[node]
    end
    return found
end

"""Given species names (one, several, or `nothing` for all), their mask."""
_species_mask(clade::Clade, ::Nothing) = typemax(UInt64)
_species_mask(clade::Clade, name::AbstractString) = _bit(_species_position(clade, name))
_species_mask(clade::Clade, names) =
    reduce(|, (_species_mask(clade, name) for name in names); init = zero(UInt64))

"""Given a clade's leaf, its gene ID, or `missing` when it is placed on none."""
function _gene_id(clade::Clade, leaf::Integer)
    forest = clade.gene_trees
    member = clade.species[_leaf_position(forest, leaf)]
    return something(_leaf_gene_id(forest, member, leaf), missing)
end

"""
    orthologs(clade, species, id; species_wanted = nothing)

Given a gene (its gene ID or, for a genome-less species, protein accession),
find its orthologs, optionally only in `species_wanted` (a name or several).
Returns a `DataFrame` (`species`, `gene`, `protein`, `one_to_one`), empty for a
gene in no orthogroup; `one_to_one` marks a pair each the other's only ortholog
in that species. `gene` is `missing` for a genome-less species.
"""
function orthologs(
    clade::Clade,
    species::AbstractString,
    id::AbstractString;
    species_wanted = nothing,
)
    position = _species_position(clade, species)
    wanted = _species_mask(clade, species_wanted)
    leaf = _leaf(clade, position, id)
    forest = clade.gene_trees
    found = leaf == 0 ? Int32[] : _ortholog_leaves(forest, leaf, wanted)

    positions = [_leaf_position(forest, node) for node in found]
    per_species = Dict{Int,Int}()
    for found_position in positions
        per_species[found_position] = get(per_species, found_position, 0) + 1
    end
    one_to_one = Bool[
        per_species[found_position] == 1 &&
            length(_ortholog_leaves(forest, node, _bit(position))) == 1 for
        (node, found_position) in zip(found, positions)
    ]
    return DataFrame(
        species = species_names(clade)[positions],
        gene = Union{Missing,String}[_gene_id(clade, node) for node in found],
        protein = forest.protein[found],
        one_to_one = one_to_one,
    )
end

"""
Given two sibling subtrees `from` and `to`, append every `(x, y)` with `x` a
`bit_a` leaf of `from` and `y` a `bit_b` leaf of `to` to `left`/`right`.
Returns `nothing`.
"""
function _append_leaf_pairs!(
    left::Vector{Int32},
    right::Vector{Int32},
    forest::GeneForest,
    from::Int32,
    to::Int32,
    bit_a::UInt64,
    bit_b::UInt64,
)
    (forest.species[from] & bit_a != 0 && forest.species[to] & bit_b != 0) || return nothing
    for x in _leaves_in(forest, from, bit_a), y in _leaves_in(forest, to, bit_b)
        push!(left, x)
        push!(right, y)
    end
    return nothing
end

"""
Given a forest node, append every ortholog pair it decides — a leaf of `bit_a`
under one child, a leaf of `bit_b` under another sharing no species — to
`left`/`right`. Returns `nothing`.
"""
function _node_pairs!(
    left,
    right,
    forest::GeneForest,
    node::Int32,
    bit_a::UInt64,
    bit_b::UInt64,
)
    both = bit_a | bit_b
    for first_child in _Children(forest.subtree_end, node)
        first_mask = forest.species[first_child]
        first_mask & both == 0 && continue
        later =
            _Children(forest.subtree_end, node, forest.subtree_end[first_child] + Int32(1))
        for second_child in later
            second_mask = forest.species[second_child]
            (second_mask & both != 0 && first_mask & second_mask == 0) || continue
            _append_leaf_pairs!(
                left,
                right,
                forest,
                first_child,
                second_child,
                bit_a,
                bit_b,
            )
            _append_leaf_pairs!(
                left,
                right,
                forest,
                second_child,
                first_child,
                bit_a,
                bit_b,
            )
        end
    end
    return nothing
end

"""
    ortholog_pairs(clade, species_a, species_b; codeml = nothing)

Given two species of `clade`, list every ortholog pair between them, as
OrthoFinder's `Orthologues/<a>__v__<b>.tsv` does but one pair per row. Returns a
`DataFrame` (`orthogroup`, `gene_a`, `gene_b`, `protein_a`, `protein_b`,
`one_to_one`), `one_to_one` marking a pair whose genes have no other ortholog.
With `codeml` ([`CodemlTrees`](@ref)), a `dS` column adds each pair's `patristic`
dS on its orthogroup's divergence tree, `NaN` where that tree lacks either gene.
"""
function ortholog_pairs(
    clade::Clade,
    species_a::AbstractString,
    species_b::AbstractString;
    codeml = nothing,
)
    bit_a = _bit(_species_position(clade, species_a))
    bit_b = _bit(_species_position(clade, species_b))
    bit_a == bit_b && throw(ArgumentError("ortholog_pairs needs two different species"))
    forest = clade.gene_trees
    left, right, trees = Int32[], Int32[], Int32[]
    for tree in eachindex(forest.orthogroups)
        forest.source[tree] == TREE_PAIR && continue
        for node in _tree_nodes(forest, tree)
            mask = forest.species[node]
            (_is_leaf(forest, node) || mask & bit_a == 0 || mask & bit_b == 0) && continue
            _node_pairs!(left, right, forest, node, bit_a, bit_b)
            append!(trees, Iterators.repeated(Int32(tree), length(left) - length(trees)))
        end
    end
    table = _pair_table(clade, left, right, trees)
    isnothing(codeml) || (table.dS = _pair_dS(clade, codeml, left, right, trees))
    return table
end

"""
Given ortholog pairs as forest leaf vectors and their trees, read each pair's dS
off its orthogroup's divergence tree, built once per orthogroup. Returns a
`Vector{Float64}`, `NaN` without a codeml tree or where it lacks either gene.
"""
function _pair_dS(clade::Clade, codeml, left, right, trees)
    dS = fill(NaN, length(left))
    built = Dict{Int32,Any}()
    for index in eachindex(left)
        divergence =
            get!(() -> _divergence_of(clade, codeml, trees[index]), built, trees[index])
        isnothing(divergence) && continue
        tree, by_forest = divergence
        first_leaf, second_leaf =
            get(by_forest, left[index], 0), get(by_forest, right[index], 0)
        (first_leaf == 0 || second_leaf == 0) && continue
        dS[index] = patristic(tree, first_leaf, second_leaf)
    end
    return dS
end

"""Given ortholog pairs as leaf vectors, tabulate them for [`ortholog_pairs`](@ref)."""
function _pair_table(
    clade::Clade,
    left::Vector{Int32},
    right::Vector{Int32},
    trees::Vector{Int32},
)
    forest = clade.gene_trees
    counts = Dict{Int32,Int}()
    foreach(
        leaf -> counts[leaf] = get(counts, leaf, 0) + 1,
        Iterators.flatten((left, right)),
    )
    return DataFrame(
        orthogroup = forest.orthogroups[trees],
        gene_a = Union{Missing,String}[_gene_id(clade, leaf) for leaf in left],
        gene_b = Union{Missing,String}[_gene_id(clade, leaf) for leaf in right],
        protein_a = forest.protein[left],
        protein_b = forest.protein[right],
        one_to_one = Bool[counts[x] == 1 && counts[y] == 1 for (x, y) in zip(left, right)],
    )
end

"""Given a forest node and a node below it, the child of `node` holding `below`."""
function _child_holding(forest::GeneForest, node::Int32, below::Integer)
    for child in _Children(forest.subtree_end, node)
        child <= below <= forest.subtree_end[child] && return child
    end
    throw(ArgumentError("node $below is not below node $node"))
end

"""
    common_ancestor(clade, species_a => id_a, species_b => id_b)

Given two genes, find the gene-tree node where their lineages meet. Returns
`(node, event, species_node, depth)`, or `nothing` when they share no orthogroup.
`event` is `:speciation` or `:duplication`, decided by the two children holding
them (so per pair at a polytomy), or `:unresolved` in a two-gene orthogroup;
`species_node` is where `node` maps in the species tree, `depth` as `lca_depth`.
"""
function common_ancestor(
    clade::Clade,
    first_gene::Pair{<:AbstractString,<:AbstractString},
    second_gene::Pair{<:AbstractString,<:AbstractString},
)
    forest = clade.gene_trees
    x = _leaf(clade, _species_position(clade, first_gene.first), first_gene.second)
    y = _leaf(clade, _species_position(clade, second_gene.first), second_gene.second)
    (x == 0 || y == 0) && return nothing
    x == y && throw(ArgumentError("both sides name the same gene"))
    tree = _tree_of(forest, x)
    tree == _tree_of(forest, y) || return nothing

    side, node = x, forest.parent[x]
    while !(node <= y <= forest.subtree_end[node])
        side, node = node, forest.parent[node]
    end
    other = _child_holding(forest, node, y)
    event = if forest.source[tree] == TREE_PAIR
        :unresolved
    else
        forest.species[side] & forest.species[other] == 0 ? :speciation : :duplication
    end
    mapped = _species_node(clade.species_tree, forest.species[node])
    return (
        node = node,
        event = event,
        species_node = clade.species_tree.labels[mapped],
        depth = Int(clade.species_tree.depth[mapped]),
    )
end

"""
    orthogroup(clade, species, id)

Given a gene (its gene ID, or protein accession), find its orthogroup. Returns
the orthogroup's name, or `nothing` when OrthoFinder assigned it none.
"""
function orthogroup(clade::Clade, species::AbstractString, id::AbstractString)
    leaf = _leaf(clade, _species_position(clade, species), id)
    return leaf == 0 ? nothing :
           clade.gene_trees.orthogroups[_tree_of(clade.gene_trees, leaf)]
end

"""Given a clade, every orthogroup's name, sorted."""
orthogroups(clade::Clade) = copy(clade.gene_trees.orthogroups)

#= Gene-tree access =#

"""
Given a clade and an orthogroup name, view its tree. Returns a [`GeneTree`](@ref);
throws an `ArgumentError` when the clade has no such orthogroup.
"""
function gene_tree(clade::Clade, orthogroup::AbstractString)
    names = clade.gene_trees.orthogroups
    position = searchsortedfirst(names, orthogroup)
    (position <= length(names) && names[position] == orthogroup) ||
        throw(ArgumentError("no orthogroup `$orthogroup` in the clade"))
    return GeneTree(clade, Int32(position))
end

"""Given a gene tree, its nodes: a range of forest indices, root first (preorder)."""
tree_nodes(tree::GeneTree) = _tree_nodes(tree.clade.gene_trees, tree.orthogroup)

"""Given a gene tree, its root node."""
root(tree::GeneTree) = first(tree_nodes(tree))

"""
Given a gene tree and a node, check the node is the tree's. Returns it as an
`Int32`; throws an `ArgumentError` when it lies in another tree.
"""
function _checked(tree::GeneTree, node::Integer)
    node in tree_nodes(tree) || throw(ArgumentError("node $node is not in $(tree)"))
    return Int32(node)
end

"""Given a gene tree node, tell whether it is a leaf."""
is_leaf(tree::GeneTree, node::Integer) =
    _is_leaf(tree.clade.gene_trees, _checked(tree, node))

"""Given a gene tree node, its children, in Newick order."""
children(tree::GeneTree, node::Integer) =
    collect(_Children(tree.clade.gene_trees.subtree_end, _checked(tree, node)))

"""Given a gene tree node, its parent, or `nothing` at the root."""
function Base.parent(tree::GeneTree, node::Integer)
    parent = tree.clade.gene_trees.parent[_checked(tree, node)]
    return parent == 0 ? nothing : parent
end

"""Given a gene tree node (the root by default), the leaves under it, in Newick order."""
leaves(tree::GeneTree, node::Integer = root(tree)) =
    collect(_leaves_in(tree.clade.gene_trees, _checked(tree, node), typemax(UInt64)))

"""Given a gene tree node, the length of the branch above it (`NaN` when absent)."""
branch_length(tree::GeneTree, node::Integer) =
    tree.clade.gene_trees.branch_length[_checked(tree, node)]

"""Given a gene tree leaf, its species' name."""
leaf_species(tree::GeneTree, leaf::Integer) =
    tree.clade.species[_leaf_position(tree.clade.gene_trees, _checked(tree, leaf))].name

"""Given a gene tree leaf, its protein accession."""
leaf_protein(tree::GeneTree, leaf::Integer) =
    tree.clade.gene_trees.protein[_checked(tree, leaf)]

"""Given a gene tree leaf, its gene ID, or `missing` when it is placed on none."""
leaf_gene(tree::GeneTree, leaf::Integer) = _gene_id(tree.clade, _checked(tree, leaf))

"""
Given a gene tree node, its OrthoFinder label: `<species>_<protein>` for a leaf,
`n<k>` for the `k`-th internal node in preorder (the root being `n0`).
"""
function node_label(tree::GeneTree, node::Integer)
    forest = tree.clade.gene_trees
    node = _checked(tree, node)
    _is_leaf(forest, node) &&
        return string(leaf_species(tree, node), '_', forest.protein[node])
    return "n$(count(earlier -> !_is_leaf(forest, earlier), root(tree):(node-Int32(1))))"
end

"""
Given a gene tree, write it as Newick, with OrthoFinder's labels and every
stored branch length. Returns the `String`, ending in `;`.
"""
function newick(tree::GeneTree)
    io = IOBuffer()
    _write_newick!(io, tree, root(tree), Ref(0))
    print(io, ';')
    return String(take!(io))
end

"""Given a node, write its subtree as Newick, numbering internal nodes from `n_internal`."""
function _write_newick!(io::IO, tree::GeneTree, node::Int32, n_internal::Base.RefValue{Int})
    forest = tree.clade.gene_trees
    if _is_leaf(forest, node)
        print(io, node_label(tree, node))
    else
        label = "n$(n_internal[])"
        n_internal[] += 1
        print(io, '(')
        for (position, child) in enumerate(_Children(forest.subtree_end, node))
            position > 1 && print(io, ',')
            _write_newick!(io, tree, child, n_internal)
        end
        print(io, ')', label)
    end
    branch = forest.branch_length[node]
    isnan(branch) || print(io, ':', branch)
    return nothing
end

include("divergence_trees.jl")

#= Base.show overloads =#

Base.show(io::IO, tree::SpeciesTree) =
    print(io, "SpeciesTree(", join(species_names(tree), ", "), ")")

Base.show(io::IO, forest::GeneForest) = print(
    io,
    "GeneForest(",
    length(forest.orthogroups),
    " trees, ",
    length(forest.parent),
    " nodes)",
)

Base.show(io::IO, clade::Clade) = print(
    io,
    "Clade(",
    join(species_names(clade), ", "),
    "; ",
    length(clade.gene_trees.orthogroups),
    " orthogroups)",
)

Base.show(io::IO, tree::GeneTree) = print(
    io,
    "GeneTree(\"",
    tree.clade.gene_trees.orthogroups[tree.orthogroup],
    "\", ",
    length(leaves(tree)),
    " leaves)",
)

export Clade,
    CodemlTrees,
    DivergenceTree,
    GeneForest,
    GeneIndex,
    GeneTree,
    SpeciesTree,
    branch_length,
    children,
    common_ancestor,
    divergence_tree,
    gene_tree,
    is_leaf,
    leaf_gene,
    leaf_protein,
    leaf_species,
    leaves,
    load_clade,
    nearest_orthologs,
    newick,
    node_ages,
    node_event,
    node_label,
    ortholog_pairs,
    orthogroup,
    orthogroups,
    orthologs,
    patristic,
    prune,
    root,
    root_distances,
    species_names,
    species_node,
    tree_nodes

end

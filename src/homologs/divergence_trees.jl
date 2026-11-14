#= Divergence trees: codeml's per-branch dS/dN, rooted on OrthoFinder's gene trees.
   Included into `Orthologs`; `docs/divergence_tree.md` is the reference. =#

"""
A directory of codeml trees, `<orthogroup>_dS.nwk` (and optionally `_dN.nwk`) per
orthogroup. `orthogroups` lists, sorted, those with a dS tree. Build with
`CodemlTrees(directory)`; read one with [`divergence_tree`](@ref).
"""
struct CodemlTrees
    directory::String
    orthogroups::Vector{String}
end

"""
Given a directory of codeml trees, index it. Returns a [`CodemlTrees`](@ref);
throws when the directory holds no `<orthogroup>_dS.nwk`.
"""
function CodemlTrees(directory::AbstractString)
    found = String[]
    for name in readdir(directory)
        matched = match(r"^(.+)_dS\.nwk$", name)
        isnothing(matched) || push!(found, String(matched.captures[1]))
    end
    isempty(found) && throw(ArgumentError("no `<orthogroup>_dS.nwk` in $directory"))
    return CodemlTrees(String(directory), sort!(found))
end

"""Given codeml trees, their orthogroups (those with a dS tree), sorted."""
orthogroups(trees::CodemlTrees) = copy(trees.orthogroups)

"""Given codeml trees and an orthogroup name, tell whether it has a dS tree."""
function Base.in(orthogroup::AbstractString, trees::CodemlTrees)
    position = searchsortedfirst(trees.orthogroups, orthogroup)
    return position <= length(trees.orthogroups) &&
           trees.orthogroups[position] == orthogroup
end

"""
One orthogroup's codeml tree, rooted where OrthoFinder roots its gene tree, in
preorder with local node indices (root 1). Per node: `parent` (0 at the root),
`subtree_end`, the `dS`/`dN` of the branch above it, its age in each
(`dS_age`/`dN_age`, see [`node_ages`](@ref)), `forest_node` (a leaf's
forest leaf; else the forest LCA of its leaves), `exact` (OrthoFinder has this
clade), `species` (its lineage's mask, dropped genes included) and `event`.
`dropped` holds the orthogroup's forest leaves codeml lacks. Rules in
`docs/divergence_tree.md`.
"""
struct DivergenceTree
    clade::Clade
    orthogroup::Int32
    parent::Vector{Int32}
    subtree_end::Vector{Int32}
    dS::Vector{Float32}
    dN::Vector{Float32}
    dS_age::Vector{Float64}
    dN_age::Vector{Float64}
    forest_node::Vector{Int32}
    exact::Vector{Bool}
    species::Vector{UInt64}
    event::Vector{Symbol}
    dropped::Vector{Int32}
end

"""One direction of an undirected codeml branch: the node it leads `to`, and its lengths."""
struct _Arc
    to::Int32
    dS::Float32
    dN::Float32
end

#= Reading =#

"""
Given a codeml dS Newick and optionally its dN twin, lay both out. Returns
`(layout, dN)`, `dN` the dN tree's branch lengths (`NaN` without one); throws
when the two differ in shape or labels.
"""
function _codeml_layout(dS_newick::AbstractString, dN_newick)
    layout = _preorder_layout(dS_newick)
    isnothing(dN_newick) && return layout, fill(NaN32, length(layout.parent))
    dN_layout = _preorder_layout(dN_newick)
    (dN_layout.parent == layout.parent && dN_layout.labels == layout.labels) ||
        throw(ArgumentError("the dN tree differs from the dS tree in shape or labels"))
    return layout, dN_layout.branch_length
end

"""
Given a forest tree and a codeml layout, find each codeml leaf's forest leaf by
protein accession. Returns a vector over codeml nodes, 0 at internal ones;
throws on an accession the orthogroup lacks or one leaf repeated.
"""
function _codeml_leaves(forest::GeneForest, tree::Integer, layout)
    by_protein = Dict(
        forest.protein[node] => node for
        node in _tree_nodes(forest, tree) if _is_leaf(forest, node)
    )
    forest_leaf = zeros(Int32, length(layout.parent))
    for node in eachindex(layout.parent)
        layout.subtree_end[node] == node || continue
        leaf = get(by_protein, layout.labels[node], nothing)
        isnothing(leaf) && throw(
            ArgumentError(
                "codeml leaf `$(layout.labels[node])` is not in $(forest.orthogroups[tree])",
            ),
        )
        forest_leaf[node] = leaf
    end
    leaves = filter(!iszero, forest_leaf)
    allunique(leaves) || throw(ArgumentError("a codeml tree repeats a leaf"))
    length(leaves) >= 2 || throw(ArgumentError("a codeml tree needs two leaves"))
    return forest_leaf
end

"""
Given a codeml layout, its dN lengths and leaf map, unroot it: per node, the
branches to its neighbours. Returns that adjacency, a node of degree 2 (codeml's
root, when bifurcating) dissolved into one branch with summed lengths.
"""
function _undirected(layout, dN::Vector{Float32}, forest_leaf::Vector{Int32})
    adjacency = [_Arc[] for _ in layout.parent]
    for (node, parent) in enumerate(layout.parent)
        parent == 0 && continue
        dS = layout.branch_length[node]
        push!(adjacency[node], _Arc(parent, dS, dN[node]))
        push!(adjacency[parent], _Arc(node, dS, dN[node]))
    end
    for node in eachindex(adjacency)
        (length(adjacency[node]) == 2 && forest_leaf[node] == 0) || continue
        first_arc, second_arc = adjacency[node]
        dS, dN_sum = first_arc.dS + second_arc.dS, first_arc.dN + second_arc.dN
        _redirect!(adjacency[first_arc.to], node, _Arc(second_arc.to, dS, dN_sum))
        _redirect!(adjacency[second_arc.to], node, _Arc(first_arc.to, dS, dN_sum))
        empty!(adjacency[node])
    end
    return adjacency
end

"""Given a node's arcs, replace the one leading to `from` with `arc`. Returns `nothing`."""
function _redirect!(arcs::Vector{_Arc}, from::Integer, arc::_Arc)
    arcs[findfirst(existing -> existing.to == from, arcs)] = arc
    return nothing
end

#= Rooting =#

"""Given two forest nodes, their lowest common ancestor."""
function _forest_lca(forest::GeneForest, first_node::Integer, second_node::Integer)
    low, high = minmax(first_node, second_node)
    node = Int32(low)
    while forest.subtree_end[node] < high
        node = forest.parent[node]
    end
    return node
end

"""Given forest leaves, their lowest common ancestor."""
_forest_lca(forest::GeneForest, leaves::AbstractVector{<:Integer}) =
    _forest_lca(forest, minimum(leaves), maximum(leaves))

"""
Given the codeml-to-forest leaf map, group codeml's leaves by the child of their
forest LCA holding them, as OrthoFinder's root splits them. Returns `(lca,
group, sizes)`: `group` per codeml node (0 at internal ones), `sizes` per group.
"""
function _root_groups(forest::GeneForest, forest_leaf::Vector{Int32})
    lca = _forest_lca(forest, filter(!iszero, forest_leaf))
    tops = collect(_Children(forest.subtree_end, lca))
    group = zeros(Int32, length(forest_leaf))
    for (node, leaf) in enumerate(forest_leaf)
        leaf == 0 && continue
        group[node] = findfirst(top -> top <= leaf <= forest.subtree_end[top], tops)
    end
    sizes = [count(==(index), group) for index in eachindex(tops)]
    keep = findall(!iszero, sizes)
    renumber = zeros(Int32, length(tops))
    renumber[keep] .= 1:length(keep)
    return lca, [index == 0 ? Int32(0) : renumber[index] for index in group], sizes[keep]
end

"""
Given an adjacency and a start node, walk it depth first. Returns `(order,
from, last)`: nodes in preorder, each node's predecessor (0 at `start`), and the
`order` position of each node's last descendant.
"""
function _walk(adjacency::Vector{Vector{_Arc}}, start::Integer)
    order, from, last =
        Int32[], zeros(Int32, length(adjacency)), zeros(Int32, length(adjacency))
    stack = [(Int32(start), Int32(0), false)]
    while !isempty(stack)
        node, predecessor, done = pop!(stack)
        if done
            last[node] = length(order)
            continue
        end
        push!(order, node)
        from[node] = predecessor
        push!(stack, (node, predecessor, true))
        for arc in Iterators.reverse(adjacency[node])
            arc.to == predecessor || push!(stack, (arc.to, node, false))
        end
    end
    return order, from, last
end

"""
Given a walk and each leaf's root group, find the branch splitting off each
group. Returns, per group, `(node, inside)`: the branch above `node` in the
walk, `inside` true when the group lies below it. Throws when one is missing.
"""
function _group_branches(order, from, last, group::Vector{Int32}, sizes::Vector{Int})
    position = zeros(Int32, length(from))
    position[order] .= 1:length(order)
    counts = zeros(Int32, length(sizes), length(from))
    n_below = zeros(Int32, length(from))
    for node in Iterators.reverse(order)
        group[node] == 0 || (counts[group[node], node] += 1; n_below[node] += 1)
        from[node] == 0 && continue
        counts[:, from[node]] .+= view(counts, :, node)
        n_below[from[node]] += n_below[node]
    end
    total = sum(sizes)
    return map(eachindex(sizes)) do index
        size = sizes[index]
        for node in order
            from[node] == 0 && continue
            inside = counts[index, node] == size == n_below[node]
            outside = counts[index, node] == 0 && n_below[node] == total - size
            (inside || outside) && return (node, inside)
        end
        throw(
            ArgumentError(
                "the codeml tree does not split off OrthoFinder's root clade $index",
            ),
        )
    end
end

"""
Given a forest node and an ancestor, sum the forest branch lengths between them.
Returns the sum, `NaN` when any is absent.
"""
function _forest_path_length(forest::GeneForest, node::Integer, ancestor::Integer)
    total = 0.0
    while node != ancestor
        total += forest.branch_length[node]
        node = forest.parent[node]
    end
    return total
end

"""
Given OrthoFinder's root and its two codeml leaf groups, find where the root
falls along codeml's branch between them. Returns the fraction from group 1's
side: forest path lengths below the root in proportion, else 0.5.
"""
function _root_fraction(forest::GeneForest, lca::Integer, forest_leaf, group)
    lengths = map(1:2) do index
        leaves = forest_leaf[group .== index]
        _forest_path_length(forest, _forest_lca(forest, leaves), lca)
    end
    total = sum(lengths)
    return isfinite(total) && total > 0 ? lengths[1] / total : 0.5
end

"""
Given an adjacency, a start node and the nodes it may cross (`core`), sum dS
outward from `start`. Returns a `node => distance` `Dict` over `core`.
"""
function _core_distances(adjacency::Vector{Vector{_Arc}}, core::Set{Int32}, start::Int32)
    distances, stack = Dict(start => 0.0), [start]
    while !isempty(stack)
        node = pop!(stack)
        for arc in adjacency[node]
            (arc.to in core && !haskey(distances, arc.to)) || continue
            distances[arc.to] = distances[node] + arc.dS
            push!(stack, arc.to)
        end
    end
    return distances
end

"""
Given each root group's branch, find the codeml nodes outside every group (where
OrthoFinder's root polytomy falls) and pick the one whose farthest group is
nearest by dS. Returns that node; throws when no node lies outside them all.
"""
function _core_root(adjacency, order, from, last, branches)
    position = zeros(Int32, length(from))
    position[order] .= 1:length(order)
    covered, attached = falses(length(order)), Int32[]
    for (node, inside) in branches
        if inside
            covered[position[node]:last[node]] .= true
            push!(attached, from[node])
        else
            covered[1:(position[node]-1)] .= true
            covered[(last[node]+1):end] .= true
            push!(attached, node)
        end
    end
    core = Set(order[.!covered])
    isempty(core) &&
        throw(ArgumentError("the codeml tree conflicts with OrthoFinder's root"))
    best, best_reach = Int32(0), Inf
    for candidate in sort!(collect(core))
        distances = _core_distances(adjacency, core, candidate)
        reach = maximum(distances[node] for node in attached)
        reach < best_reach && ((best, best_reach) = (candidate, reach))
    end
    return best
end

"""
Given the branch between `near` and `far`, place a new root on it, `fraction` of
its lengths from `near`. Returns the new root's node, appended to `adjacency`
(and a 0 to `forest_leaf`).
"""
function _split_branch!(adjacency, forest_leaf, near::Int32, far::Int32, fraction::Real)
    arc = adjacency[near][findfirst(existing -> existing.to == far, adjacency[near])]
    root = Int32(length(adjacency) + 1)
    near_dS, near_dN = Float32(fraction * arc.dS), Float32(fraction * arc.dN)
    far_dS, far_dN = arc.dS - near_dS, arc.dN - near_dN
    _redirect!(adjacency[near], far, _Arc(root, near_dS, near_dN))
    _redirect!(adjacency[far], near, _Arc(root, far_dS, far_dN))
    push!(adjacency, [_Arc(near, near_dS, near_dN), _Arc(far, far_dS, far_dN)])
    push!(forest_leaf, 0)
    return root
end

"""
Given the unrooted codeml tree and its forest leaves, root it where OrthoFinder
roots the orthogroup. Returns the root node: a new one on the branch between
OrthoFinder's two root clades, else the centre of its root polytomy.
"""
function _place_root!(adjacency, forest_leaf::Vector{Int32}, forest::GeneForest)
    lca, group, sizes = _root_groups(forest, forest_leaf)
    order, from, last = _walk(adjacency, findfirst(!isempty, adjacency))
    branches = _group_branches(order, from, last, group, sizes)
    length(sizes) > 2 && return _core_root(adjacency, order, from, last, branches)
    node, inside = first(branches)
    near, far = inside ? (node, from[node]) : (from[node], node)
    fraction = _root_fraction(forest, lca, forest_leaf, group)
    return _split_branch!(adjacency, forest_leaf, near, far, fraction)
end

"""
Given the rooted adjacency, lay it out in preorder, children ordered by their
first forest leaf so leaves follow OrthoFinder's order. Returns `(parent,
subtree_end, dS, dN, leaf)`, `leaf` the forest leaf (0 at internal nodes).
"""
function _rooted_layout(adjacency, root::Integer, forest_leaf::Vector{Int32})
    order, from, _ = _walk(adjacency, root)
    key = [leaf == 0 ? typemax(Int32) : leaf for leaf in forest_leaf]
    for node in Iterators.reverse(order)
        from[node] == 0 || (key[from[node]] = min(key[from[node]], key[node]))
    end
    parent, subtree_end, dS, dN, leaf = Int32[], Int32[], Float32[], Float32[], Int32[]
    placed = zeros(Int32, length(adjacency))
    stack = [(Int32(root), _Arc(0, NaN32, NaN32), false)]
    while !isempty(stack)
        node, arc, done = pop!(stack)
        if done
            subtree_end[placed[node]] = length(parent)
            continue
        end
        push!(parent, from[node] == 0 ? Int32(0) : placed[from[node]])
        push!(subtree_end, 0)
        push!(dS, arc.dS)
        push!(dN, arc.dN)
        push!(leaf, forest_leaf[node])
        placed[node] = length(parent)
        push!(stack, (node, arc, true))
        below = filter(next -> next.to != from[node], adjacency[node])
        for next in sort!(below; by = next -> key[next.to], rev = true)
            push!(stack, (next.to, next, false))
        end
    end
    return (; parent, subtree_end, dS, dN, leaf)
end

#= Ages =#

"""
Given a preorder tree and per-branch `lengths`, date every node by the mean path
length from it down to the leaves below it, raised where needed to its oldest
child's age. Returns a `Vector{Float64}`: 0 at leaves, `NaN` when lengths are.
"""
function _node_ages(parent::Vector{Int32}, subtree_end::Vector{Int32}, lengths)
    n_nodes = length(parent)
    depth = zeros(Float64, n_nodes)
    for node = 2:n_nodes
        depth[node] = depth[parent[node]] + lengths[node]
    end
    depth_sums, n_leaves = zeros(Float64, n_nodes), zeros(Int, n_nodes)
    for node in reverse(1:n_nodes)
        subtree_end[node] == node && (depth_sums[node] = depth[node]; n_leaves[node] = 1)
        up = parent[node]
        up == 0 && continue
        depth_sums[up] += depth_sums[node]
        n_leaves[up] += n_leaves[node]
    end
    ages = depth_sums ./ n_leaves .- depth
    # Reverse preorder reaches every child before its parent, so each age is final when it
    # raises its parent's.
    for node in reverse(2:n_nodes)
        ages[parent[node]] = max(ages[parent[node]], ages[node])
    end
    return ages
end

#= Annotating =#

"""
Given child species masks, decide their parent's event pair by pair. Returns
`:duplication` when every pair shares a species, `:speciation` when none does,
else `:mixed`.
"""
function _pair_event(masks::AbstractVector{UInt64})
    n_pairs, n_shared = 0, 0
    for first_index in eachindex(masks), second_index = (first_index+1):lastindex(masks)
        n_pairs += 1
        masks[first_index] & masks[second_index] != 0 && (n_shared += 1)
    end
    return n_shared == 0 ? :speciation : n_shared == n_pairs ? :duplication : :mixed
end

"""
Given a forest tree and the forest leaves codeml kept, count them under any node
of that tree. Returns a function `node -> count`.
"""
function _kept_counter(forest::GeneForest, tree::Integer, kept_leaves)
    offset = first(_tree_nodes(forest, tree)) - 1
    present = zeros(Int32, length(_tree_nodes(forest, tree)) + 1)
    foreach(leaf -> present[leaf-offset+1] = 1, kept_leaves)
    cumsum!(present, present)
    return node -> present[forest.subtree_end[node]-offset+1] - present[node-offset]
end

"""
Given a rooted layout, map each node onto the forest: its LCA there, and whether
the forest node holds exactly its codeml leaves. Returns `(forest_node, exact)`;
throws when a codeml clade splits an OrthoFinder one ([`_check_compatible`](@ref)).
"""
function _forest_mapping(forest::GeneForest, tree::Integer, layout)
    n_nodes = length(layout.parent)
    low, high = copy(layout.leaf), copy(layout.leaf)
    n_leaves = Int32.(layout.leaf .!= 0)
    for node = n_nodes:-1:2
        parent = layout.parent[node]
        low[parent] = low[parent] == 0 ? low[node] : min(low[parent], low[node])
        high[parent] = max(high[parent], high[node])
        n_leaves[parent] += n_leaves[node]
    end
    kept_under = _kept_counter(forest, tree, filter(!iszero, layout.leaf))
    forest_node = [_forest_lca(forest, low[node], high[node]) for node = 1:n_nodes]
    exact = [kept_under(forest_node[node]) == n_leaves[node] for node = 1:n_nodes]
    _check_compatible(forest, layout.parent, forest_node, exact, n_leaves, kept_under)
    return forest_node, exact
end

"""
Given each node's forest node, check every codeml-only clade is a union of whole
children of its forest node (counting only kept leaves). Returns `nothing`;
throws an `ArgumentError` when one splits an OrthoFinder clade.
"""
function _check_compatible(forest, parent, forest_node, exact, n_leaves, kept_under)
    touched = Dict{Int32,Set{Int32}}()
    for node in reverse(eachindex(parent))
        up = parent[node]
        (up == 0 || exact[up]) && continue
        tops = get!(Set{Int32}, touched, up)
        if forest_node[node] == forest_node[up]
            union!(tops, touched[node])
        else
            push!(tops, _child_holding(forest, forest_node[up], forest_node[node]))
        end
    end
    for (node, tops) in touched
        sum(kept_under, tops) == n_leaves[node] ||
            throw(ArgumentError("the codeml tree splits a clade of OrthoFinder's tree"))
    end
    return nothing
end

"""
Given a rooted layout and its forest nodes, find each node's lineage mask: the
species of the child of its parent's forest node that holds it, so genes codeml
dropped still count; the OR of its children's where it shares that forest node
(a codeml-only split); the forest node's own at the root. Returns the masks.
"""
function _lineage_masks(
    forest::GeneForest,
    parent::Vector{Int32},
    forest_node::Vector{Int32},
)
    species, below = zeros(UInt64, length(parent)), zeros(UInt64, length(parent))
    for node in reverse(eachindex(parent))
        up = parent[node]
        species[node] = if up == 0
            forest.species[forest_node[node]]
        elseif forest_node[node] == forest_node[up]
            below[node]
        else
            forest.species[_child_holding(forest, forest_node[up], forest_node[node])]
        end
        up == 0 || (below[up] |= species[node])
    end
    return species
end

"""
Given a preorder tree's `parent`/`subtree_end` and the masks deciding events,
decide every node's event. Returns a `Symbol` per node: `:leaf`, or by
[`_pair_event`](@ref) over its children's masks; `:unresolved` when `unresolved`.
"""
function _node_events(parent, subtree_end, masks::Vector{UInt64}; unresolved::Bool = false)
    return map(eachindex(parent)) do node
        subtree_end[node] == node && return :leaf
        unresolved && return :unresolved
        _pair_event([masks[child] for child in _Children(subtree_end, node)])
    end
end

"""
    divergence_tree(clade, orthogroup; dS, dN = nothing)
    divergence_tree(clade, trees::CodemlTrees, orthogroup)

Given an orthogroup and its codeml Newick trees (`dS`, optionally `dN`), or a
[`CodemlTrees`](@ref) directory to read them from, root the dS tree where
OrthoFinder roots the orthogroup and map it onto the clade. Returns a
[`DivergenceTree`](@ref). Throws on a leaf the orthogroup lacks, or a tree that
conflicts with OrthoFinder's. Rules in `docs/divergence_tree.md`.
"""
function divergence_tree(
    clade::Clade,
    orthogroup::AbstractString;
    dS::AbstractString,
    dN::Union{Nothing,AbstractString} = nothing,
)
    forest = clade.gene_trees
    tree = gene_tree(clade, orthogroup).orthogroup
    codeml, dN_lengths = _codeml_layout(dS, dN)
    forest_leaf = _codeml_leaves(forest, tree, codeml)
    adjacency = _undirected(codeml, dN_lengths, forest_leaf)
    root = _place_root!(adjacency, forest_leaf, forest)
    layout = _rooted_layout(adjacency, root, forest_leaf)
    forest_node, exact = _forest_mapping(forest, tree, layout)
    species = _lineage_masks(forest, layout.parent, forest_node)
    unresolved = forest.source[tree] == TREE_PAIR
    event = _node_events(layout.parent, layout.subtree_end, species; unresolved)
    kept = Set(layout.leaf)
    dropped = Int32[
        node for
        node in _tree_nodes(forest, tree) if _is_leaf(forest, node) && !(node in kept)
    ]
    return DivergenceTree(
        clade,
        tree,
        layout.parent,
        layout.subtree_end,
        layout.dS,
        layout.dN,
        _node_ages(layout.parent, layout.subtree_end, layout.dS),
        _node_ages(layout.parent, layout.subtree_end, layout.dN),
        forest_node,
        exact,
        species,
        event,
        dropped,
    )
end

function divergence_tree(clade::Clade, trees::CodemlTrees, orthogroup::AbstractString)
    orthogroup in trees ||
        throw(ArgumentError("no codeml tree for `$orthogroup` in $(trees.directory)"))
    dN_path = joinpath(trees.directory, "$(orthogroup)_dN.nwk")
    return divergence_tree(
        clade,
        orthogroup;
        dS = read(joinpath(trees.directory, "$(orthogroup)_dS.nwk"), String),
        dN = isfile(dN_path) ? read(dN_path, String) : nothing,
    )
end

#= Pruning =#

"""Given a branch length, read an absent one (`NaN`) as 0."""
_or_zero(length::Float32) = isnan(length) ? 0.0f0 : length

"""
Given a divergence tree and a species mask, count each node's leaves in it and
decide which nodes survive pruning: kept leaves, and nodes with two or more
children holding some. Returns that `Vector{Bool}`.
"""
function _kept_nodes(tree::DivergenceTree, wanted::UInt64)
    forest = tree.clade.gene_trees
    n_kept = zeros(Int, length(tree.parent))
    n_branches = zeros(Int, length(tree.parent))
    for node in reverse(eachindex(tree.parent))
        _is_leaf(tree, node) &&
            (n_kept[node] = Int(forest.species[tree.forest_node[node]] & wanted != 0))
        parent = tree.parent[node]
        (parent == 0 || n_kept[node] == 0) && continue
        n_kept[parent] += n_kept[node]
        n_branches[parent] += 1
    end
    n_kept[1] >= 2 || throw(ArgumentError("pruning leaves fewer than two genes"))
    return [
        n_kept[node] > 0 && (_is_leaf(tree, node) || n_branches[node] >= 2) for
        node in eachindex(tree.parent)
    ]
end

"""
    prune(tree::DivergenceTree, species)

Given a divergence tree and species names (one or several), keep only their
leaves. Returns a new `DivergenceTree`: a chain left with one child becomes one
branch (lengths summed; the root's branch is its stem) whose `species` is the
full-tree lineage's at its top, and events are re-decided over the surviving
children. Throws when fewer than two genes remain.
"""
function prune(tree::DivergenceTree, species)
    wanted = _species_mask(tree.clade, species)
    kept = _kept_nodes(tree, wanted)
    n_nodes = length(tree.parent)
    anchor, top = zeros(Int32, n_nodes), collect(Int32, 1:n_nodes)
    above_dS, above_dN = zeros(Float32, n_nodes), zeros(Float32, n_nodes)
    selected, parent, placed = Int32[], Int32[], zeros(Int32, n_nodes)
    dS, dN, lineage = Float32[], Float32[], UInt64[]
    for node = 1:n_nodes
        up = tree.parent[node]
        if up != 0 && !kept[up]
            anchor[node], top[node] = anchor[up], top[up]
            above_dS[node] = above_dS[up] + _or_zero(tree.dS[up])
            above_dN[node] = above_dN[up] + _or_zero(tree.dN[up])
        elseif up != 0
            anchor[node] = placed[up]
        end
        kept[node] || continue
        push!(selected, node)
        placed[node] = length(selected)
        push!(parent, anchor[node])
        own_dS = node == 1 ? tree.dS[node] : _or_zero(tree.dS[node])
        own_dN = node == 1 ? tree.dN[node] : _or_zero(tree.dN[node])
        push!(dS, own_dS + above_dS[node])
        push!(dN, own_dN + above_dN[node])
        push!(lineage, tree.species[top[node]])
    end
    subtree_end = collect(Int32, eachindex(selected))
    for index in reverse(eachindex(selected))
        parent[index] == 0 || (
            subtree_end[parent[index]] =
                max(subtree_end[parent[index]], subtree_end[index])
        )
    end
    unresolved = tree.clade.gene_trees.source[tree.orthogroup] == TREE_PAIR
    forest = tree.clade.gene_trees
    return DivergenceTree(
        tree.clade,
        tree.orthogroup,
        parent,
        subtree_end,
        dS,
        dN,
        tree.dS_age[selected],
        tree.dN_age[selected],
        tree.forest_node[selected],
        tree.exact[selected],
        lineage,
        _node_events(parent, subtree_end, lineage; unresolved),
        filter(leaf -> forest.species[leaf] & wanted != 0, tree.dropped),
    )
end

#= Access =#

"""
Given a divergence tree and a node, check the node is the tree's. Returns it as
an `Int32`; throws an `ArgumentError` otherwise.
"""
function _checked(tree::DivergenceTree, node::Integer)
    1 <= node <= length(tree.parent) || throw(ArgumentError("node $node is not in $(tree)"))
    return Int32(node)
end

"""Given a divergence tree node, tell whether it is a leaf (unchecked)."""
_is_leaf(tree::DivergenceTree, node::Integer) = tree.subtree_end[node] == node

"""Given a divergence tree, its nodes, root (1) first, in preorder."""
tree_nodes(tree::DivergenceTree) = Int32(1):Int32(length(tree.parent))

"""Given a divergence tree, its root node: 1."""
root(tree::DivergenceTree) = Int32(1)

"""Given a divergence tree node, tell whether it is a leaf."""
is_leaf(tree::DivergenceTree, node::Integer) = _is_leaf(tree, _checked(tree, node))

"""Given a divergence tree node, its children, leaves in OrthoFinder's order."""
children(tree::DivergenceTree, node::Integer) =
    collect(_Children(tree.subtree_end, _checked(tree, node)))

"""Given a divergence tree node, its parent, or `nothing` at the root."""
function Base.parent(tree::DivergenceTree, node::Integer)
    parent = tree.parent[_checked(tree, node)]
    return parent == 0 ? nothing : parent
end

"""Given a divergence tree node (the root by default), the leaves under it, in order."""
function leaves(tree::DivergenceTree, node::Integer = root(tree))
    node = _checked(tree, node)
    return [leaf for leaf = node:tree.subtree_end[node] if _is_leaf(tree, leaf)]
end

"""Given a divergence tree and a rate (`:dS` or `:dN`), its per-branch lengths."""
function _rates(tree::DivergenceTree, rate::Symbol)
    rate === :dS && return tree.dS
    rate === :dN && return tree.dN
    throw(ArgumentError("`rate` must be :dS or :dN (got :$rate)"))
end

"""
Given a divergence tree node, the `rate` (`:dS` or `:dN`) of the branch above it.
Returns `NaN` at the root, unless pruning gave it a stem.
"""
branch_length(tree::DivergenceTree, node::Integer; rate::Symbol = :dS) =
    _rates(tree, rate)[_checked(tree, node)]

"""Given a divergence tree leaf, its forest leaf. Throws on an internal node."""
function _forest_leaf(tree::DivergenceTree, leaf::Integer)
    is_leaf(tree, leaf) || throw(ArgumentError("node $leaf is not a leaf"))
    return tree.forest_node[leaf]
end

"""Given a divergence tree leaf, its species' name."""
leaf_species(tree::DivergenceTree, leaf::Integer) =
    tree.clade.species[_leaf_position(tree.clade.gene_trees, _forest_leaf(tree, leaf))].name

"""Given a divergence tree leaf, its protein accession."""
leaf_protein(tree::DivergenceTree, leaf::Integer) =
    tree.clade.gene_trees.protein[_forest_leaf(tree, leaf)]

"""Given a divergence tree leaf, its gene ID, or `missing` when it is placed on none."""
leaf_gene(tree::DivergenceTree, leaf::Integer) =
    _gene_id(tree.clade, _forest_leaf(tree, leaf))

"""
Given a divergence tree node, its event: `:leaf`, `:speciation`, `:duplication`,
`:mixed` (a polytomy whose child pairs disagree) or `:unresolved` (a two-gene
orthogroup).
"""
node_event(tree::DivergenceTree, node::Integer) = tree.event[_checked(tree, node)]

"""
Given a divergence tree node, the species-tree node it maps to: the deepest one
holding its OrthoFinder clade's species (its children's lineages', at a
codeml-only node). Returns that node's label: a species name at a leaf.
"""
function species_node(tree::DivergenceTree, node::Integer)
    node = _checked(tree, node)
    mask =
        tree.exact[node] ? tree.clade.gene_trees.species[tree.forest_node[node]] :
        tree.species[node]
    return _species_label(tree.clade, mask)
end

"""Given a clade and a species mask, the label of the deepest species-tree node holding it."""
_species_label(clade::Clade, mask::UInt64) =
    clade.species_tree.labels[_species_node(clade.species_tree, mask)]

"""
Given a divergence tree, sum each node's `rate` (`:dS` or `:dN`) from the root.
Returns a `Vector{Float64}` over its nodes, the root at 0 (its stem excluded).
"""
function root_distances(tree::DivergenceTree; rate::Symbol = :dS)
    lengths = _rates(tree, rate)
    distances = zeros(Float64, length(tree.parent))
    for node = 2:length(tree.parent)
        distances[node] = distances[tree.parent[node]] + lengths[node]
    end
    return distances
end

"""
Given two divergence tree nodes, sum `rate` (`:dS` or `:dN`) along the path
between them. Returns it; for two leaves, codeml's pairwise `tree_dS`/`tree_dN`.
"""
function patristic(
    tree::DivergenceTree,
    first_node::Integer,
    second_node::Integer;
    rate::Symbol = :dS,
)
    lengths = _rates(tree, rate)
    low, high = minmax(_checked(tree, first_node), _checked(tree, second_node))
    meet = low
    while tree.subtree_end[meet] < high
        meet = tree.parent[meet]
    end
    total = 0.0
    for start in (low, high)
        node = start
        while node != meet
            total += lengths[node]
            node = tree.parent[node]
        end
    end
    return total
end

"""
Given a divergence tree, each node's age in `rate` (`:dS` or `:dN`): the mean
path length down to the leaves below it, at least its oldest child's. Returns a
`Vector{Float64}`, 0 at leaves; dated on the unpruned tree, so `prune` keeps them.
"""
function node_ages(tree::DivergenceTree; rate::Symbol = :dS)
    rate === :dS && return copy(tree.dS_age)
    rate === :dN && return copy(tree.dN_age)
    throw(ArgumentError("`rate` must be :dS or :dN (got :$rate)"))
end

"""Given a divergence tree, the name of its orthogroup."""
orthogroup(tree::DivergenceTree) = tree.clade.gene_trees.orthogroups[tree.orthogroup]

#= OrthoFinder gene trees =#

"""
Given a gene tree node, its event by the per-child-pair rule (`docs/clade.md`):
`:leaf`, `:speciation`, `:duplication`, `:mixed` (a polytomy whose pairs
disagree) or `:unresolved` (a two-gene orthogroup).
"""
function node_event(tree::GeneTree, node::Integer)
    forest = tree.clade.gene_trees
    node = _checked(tree, node)
    _is_leaf(forest, node) && return :leaf
    forest.source[tree.orthogroup] == TREE_PAIR && return :unresolved
    return _pair_event([
        forest.species[child] for child in _Children(forest.subtree_end, node)
    ])
end

"""Given a gene tree node, the label of the species-tree node it maps to."""
species_node(tree::GeneTree, node::Integer) =
    _species_label(tree.clade, tree.clade.gene_trees.species[_checked(tree, node)])

"""
    orthogroup(clade, id)

Given a gene ID or protein accession of any species, find its orthogroup.
Returns the name, or `nothing` when no species has it in an orthogroup.
"""
function orthogroup(clade::Clade, id::AbstractString)
    for position in eachindex(clade.species)
        leaf = _leaf(clade, position, id)
        leaf == 0 || return clade.gene_trees.orthogroups[_tree_of(clade.gene_trees, leaf)]
    end
    return nothing
end

Base.show(io::IO, tree::DivergenceTree) = print(
    io,
    "DivergenceTree(\"",
    orthogroup(tree),
    "\", ",
    count(node -> _is_leaf(tree, node), eachindex(tree.parent)),
    " leaves)",
)

Base.show(io::IO, trees::CodemlTrees) = print(
    io,
    "CodemlTrees(",
    length(trees.orthogroups),
    " orthogroups in ",
    trees.directory,
    ")",
)

#= Nearest orthologs =#

"""
Given two forest nodes of one tree, sum OrthoFinder's branch lengths along the
path between them. Returns it, `NaN` when a branch has none (a star).
"""
function _forest_patristic(forest::GeneForest, first_node::Integer, second_node::Integer)
    meet = _forest_lca(forest, first_node, second_node)
    total = 0.0
    for start in (first_node, second_node)
        node = start
        while node != meet
            total += forest.branch_length[node]
            node = forest.parent[node]
        end
    end
    return total
end

"""
Given an orthogroup's forest tree and codeml trees (or `nothing`), build its
divergence tree. Returns `(tree, local)`, `local` mapping forest leaves to its
leaves; `nothing` without a codeml tree or when it cannot be built.
"""
function _divergence_of(clade::Clade, codeml, tree::Integer)
    isnothing(codeml) && return nothing
    name = clade.gene_trees.orthogroups[tree]
    name in codeml || return nothing
    divergence = try
        divergence_tree(clade, codeml, name)
    catch error
        error isa ArgumentError || rethrow()
        return nothing
    end
    by_forest = Dict(divergence.forest_node[leaf] => leaf for leaf in leaves(divergence))
    return divergence, by_forest
end

"""
Given a leaf and its ortholog `candidates`, pick the nearest: by codeml dS when
the divergence tree holds them all, else by OrthoFinder's branch lengths, else
the first. Returns `(leaf, distance, basis)`, `basis` `:dS`, `:tree` or `:none`.
"""
function _nearest(forest::GeneForest, divergence, leaf::Int32, candidates::Vector{Int32})
    if !isnothing(divergence)
        tree, by_forest = divergence
        if haskey(by_forest, leaf) &&
           all(candidate -> haskey(by_forest, candidate), candidates)
            distances = [patristic(tree, by_forest[leaf], by_forest[c]) for c in candidates]
            best = argmin(distances)
            return candidates[best], distances[best], :dS
        end
    end
    distances = [_forest_patristic(forest, leaf, candidate) for candidate in candidates]
    all(isfinite, distances) || return first(candidates), NaN, :none
    best = argmin(distances)
    return candidates[best], distances[best], :tree
end

"""
    nearest_orthologs(clade, species; codeml = nothing, among = nothing)

Given a species of `clade`, find each other gene's nearest ortholog in it: by dS
on its orthogroup's codeml tree (a [`CodemlTrees`](@ref)) when that holds the
gene and all its orthologs there, else by OrthoFinder's branch lengths. `among`
(gene IDs or protein accessions of `species`) limits the orthologs considered to
those, so a gene whose orthologs are all outside it gets no row. Returns a
`DataFrame` (`species`, `gene`, `protein`, `ortholog_gene`, `ortholog_protein`,
`n_orthologs`, `distance`, `basis`). Rules in `docs/divergence_tree.md`.
"""
function nearest_orthologs(
    clade::Clade,
    species::AbstractString;
    codeml = nothing,
    among = nothing,
)
    forest = clade.gene_trees
    target = _bit(_species_position(clade, species))
    allowed = isnothing(among) ? nothing : Set{String}(among)
    eligible(leaf) =
        forest.protein[leaf] in allowed || coalesce(_gene_id(clade, leaf), "") in allowed
    from, to, counts, distances, bases = Int32[], Int32[], Int[], Float64[], Symbol[]
    for tree in eachindex(forest.orthogroups)
        top = first(_tree_nodes(forest, tree))
        mask = forest.species[top]
        (forest.source[tree] == TREE_PAIR || mask & target == 0 || mask == target) &&
            continue
        divergence = _divergence_of(clade, codeml, tree)
        for leaf in _tree_nodes(forest, tree)
            (_is_leaf(forest, leaf) && forest.species[leaf] != target) || continue
            candidates = _ortholog_leaves(forest, leaf, target)
            isnothing(allowed) || filter!(eligible, candidates)
            isempty(candidates) && continue
            nearest, distance, basis = _nearest(forest, divergence, leaf, candidates)
            push!(from, leaf)
            push!(to, nearest)
            push!(counts, length(candidates))
            push!(distances, distance)
            push!(bases, basis)
        end
    end
    return DataFrame(
        species = species_names(clade)[[_leaf_position(forest, leaf) for leaf in from]],
        gene = Union{Missing,String}[_gene_id(clade, leaf) for leaf in from],
        protein = forest.protein[from],
        ortholog_gene = Union{Missing,String}[_gene_id(clade, leaf) for leaf in to],
        ortholog_protein = forest.protein[to],
        n_orthologs = counts,
        distance = distances,
        basis = bases,
    )
end

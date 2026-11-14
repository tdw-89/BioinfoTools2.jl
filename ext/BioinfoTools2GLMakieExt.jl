"""
`gene_family_explorer`: an interactive GLMakie window over a clade's codeml
trees, drawn with `Plotting.gene_family_plot!`.

An extension so that GLMakie (and OpenGL) is loaded only by those who explore;
static figures need only CairoMakie.
"""
module BioinfoTools2GLMakieExt

using GLMakie

using BioinfoTools2: Orthologs
using BioinfoTools2.Orthologs:
    Clade, CodemlTrees, divergence_tree, leaf_gene, leaf_protein, leaves, orthogroup, prune
using BioinfoTools2.Plotting: gene_family_plot!, species_palette
import BioinfoTools2.Plotting: gene_family_explorer

"""
Given a clade, its codeml trees and a species mask, list the orthogroups to step
through: those with a codeml tree whose gene tree holds two or more genes of the
masked species. Returns their names, sorted.
"""
function _explorable(clade::Clade, trees::CodemlTrees, wanted::UInt64)
    forest = clade.gene_trees
    return filter(trees.orthogroups) do name
        position = searchsortedfirst(forest.orthogroups, name)
        position <= length(forest.orthogroups) && forest.orthogroups[position] == name ||
            return false
        n_genes = count(
            node ->
                forest.subtree_end[node] == node && forest.species[node] & wanted != 0,
            Orthologs._tree_nodes(forest, position),
        )
        return n_genes >= 2
    end
end

"""
Given an orthogroup name or any gene ID/protein, find its place among `choices`.
Returns `(index, message)`: the index, or `nothing` and why not.
"""
function _resolve(clade::Clade, trees::CodemlTrees, choices, query::AbstractString)
    name = query in trees ? String(query) : orthogroup(clade, query)
    isnothing(name) && return nothing, "nothing found for `$query`"
    position = searchsortedfirst(choices, name)
    position <= length(choices) && choices[position] == name && return position, ""
    return nothing,
    "$name has no codeml tree, or fewer than two genes of the chosen species"
end

"""
Given an orthogroup, build its divergence tree, pruned to `species` unless
`nothing`. Returns it, or the error message when it cannot be built.
"""
function _load(clade::Clade, trees::CodemlTrees, name::AbstractString, species)
    try
        tree = divergence_tree(clade, trees, name)
        return isnothing(species) ? tree : prune(tree, species)
    catch error
        error isa ArgumentError || rethrow()
        return sprint(showerror, error)
    end
end

"""
Given a layout, delete every block in it and detach every nested layout, which
would otherwise linger empty and skew the column widths. Returns `nothing`.
"""
function _clear!(grid::GridLayout)
    for item in copy(contents(grid))
        if item isa GridLayout
            _clear!(item)
            GLMakie.Makie.GridLayoutBase.remove_from_gridlayout!(
                GLMakie.Makie.GridLayoutBase.gridcontent(item),
            )
        else
            delete!(item)
        end
    end
    trim!(grid)
    return nothing
end

"""
Given the plotted tree's axis and its node `positions` (`(x, y)`), ring the
leaves whose gene ID or protein is `query`. Returns `nothing`.
"""
function _highlight!(axis, tree, query::AbstractString, positions)
    hits = filter(
        leaf ->
            isequal(leaf_gene(tree, leaf), query) || leaf_protein(tree, leaf) == query,
        leaves(tree),
    )
    isempty(hits) && return nothing
    x, y = positions
    scatter!(
        axis,
        x[hits],
        y[hits];
        markersize = 22,
        color = :transparent,
        strokecolor = :black,
        strokewidth = 2,
        inspectable = false,
    )
    return nothing
end

"""
What the explorer is stepping through: the orthogroups (`choices`) holding two
or more genes of the checked species, and the species each tree is pruned to
(`nothing` when every species is checked).
"""
mutable struct _Selection
    choices::Vector{String}
    species::Union{Nothing,Vector{String}}
end

"""
Given the checked species and all of the clade's, build the [`_Selection`](@ref)
over `trees`. Returns it; its `choices` are empty when nothing qualifies.
"""
function _Selection(clade::Clade, trees::CodemlTrees, checked, names)
    isempty(checked) && return _Selection(String[], String[])
    choices = _explorable(clade, trees, Orthologs._species_mask(clade, checked))
    return _Selection(
        choices,
        length(checked) == length(names) ? nothing : collect(checked),
    )
end

"""
Given the explorer's plot layout and status label, draw orthogroup `index` of the
selection, ringing `query`'s leaf. Returns `nothing`; a tree that cannot be built,
or an empty selection, leaves the layout empty and says why in `status`.
"""
function _draw!(
    plot_grid,
    status,
    clade,
    trees,
    selection::_Selection,
    index,
    query,
    options,
)
    _clear!(plot_grid)
    if isempty(selection.choices)
        status.text[] = "no orthogroup has two genes of the checked species"
        return nothing
    end
    name, place = selection.choices[index], "$index of $(length(selection.choices))"
    tree = _load(clade, trees, name, selection.species)
    if tree isa String
        status.text[] = "$name ($place): $tree"
        return nothing
    end
    drawn = gene_family_plot!(plot_grid, tree; options...)
    isempty(query) || _highlight!(drawn.tree_axis, tree, query, drawn.positions)
    status.text[] = "$name ($place)"
    return nothing
end

"""
Given a grid cell, add a checkbox per species, in its plot colour, checked when
in `checked`, wrapping onto a new row before `max_width` pixels (estimated from
the names' lengths). Returns a `name => Checkbox` `Dict`.
"""
function _species_boxes!(cell, names, checked, colors; max_width::Real)
    grid = GridLayout(cell; tellwidth = false, halign = :left)
    Label(grid[1, 1], "Species:"; halign = :left)
    boxes = Dict{String,Checkbox}()
    row, column, used = 1, 2, 90.0
    for name in names
        width = 50 + 9 * length(name)
        if used + width > max_width && column > 2
            row, column, used = row + 1, 2, 90.0
        end
        boxes[name] = Checkbox(
            grid[row, column];
            checked = name in checked,
            checkboxcolor_checked = colors[name],
            checkboxstrokecolor_checked = colors[name],
        )
        Label(grid[row, column+1], name; halign = :left)
        row, column, used = row, column + 2, used + width
    end
    return boxes
end

"""Given a `species` argument (a name, several, or `nothing` for all), list the names."""
_checked_names(names, ::Nothing) = copy(names)
_checked_names(names, species::AbstractString) = [String(species)]
_checked_names(names, species) = String.(collect(species))

function gene_family_explorer(
    clade::Clade,
    trees::CodemlTrees;
    species = nothing,
    orthogroup::AbstractString = "",
    layout::Symbol = :chronogram,
    size = (1300, 950),
    screen = GLMakie.Screen(),
    plot_options...,
)
    names = Orthologs.species_names(clade)
    checked = _checked_names(names, species)
    Orthologs._species_mask(clade, checked)   # throws on a species the clade lacks
    selection = _Selection(clade, trees, checked, names)
    isempty(selection.choices) && throw(
        ArgumentError(
            "no orthogroup has a codeml tree and two genes of the chosen species",
        ),
    )
    figure = Figure(; size)
    controls = GridLayout(figure[1, 1]; tellwidth = false, halign = :left)
    previous = Button(controls[1, 1]; label = "◀")
    textbox =
        Textbox(controls[1, 2]; placeholder = "orthogroup, gene ID or protein", width = 280)
    next = Button(controls[1, 3]; label = "▶")
    aligned = Toggle(controls[1, 4]; active = layout === :chronogram)
    Label(controls[1, 5], "align genes (ages)"; halign = :left)
    status = Label(controls[1, 6], ""; halign = :left, tellwidth = false)
    colors = get(plot_options, :species_colors, species_palette(clade))
    species_boxes =
        _species_boxes!(figure[2, 1], names, checked, colors; max_width = first(size) - 60)
    plot_grid = GridLayout(figure[3, 1])
    current, query = Observable(1), Observable("")
    redraw() = _draw!(
        plot_grid,
        status,
        clade,
        trees,
        selection,
        current[],
        query[],
        (; plot_options..., layout = aligned.active[] ? :chronogram : :phylogram),
    )
    on(_ -> redraw(), current)
    on(_ -> redraw(), aligned.active)

    go_to(text) = begin
        index, message = _resolve(clade, trees, selection.choices, strip(text))
        isnothing(index) && return status.text[] = message
        query[] = strip(text)
        current[] = index
    end
    step(offset) = begin
        isempty(selection.choices) && return nothing
        query[] = ""
        current[] = mod1(current[] + offset, length(selection.choices))
    end
    # Re-select on any checkbox, staying on the current orthogroup (or the next
    # one still qualifying).
    reselect() = begin
        shown = isempty(selection.choices) ? "" : selection.choices[current[]]
        updated = _Selection(
            clade,
            trees,
            filter(name -> species_boxes[name].checked[], names),
            names,
        )
        selection.choices, selection.species = updated.choices, updated.species
        current[] = clamp(
            searchsortedfirst(selection.choices, shown),
            1,
            max(length(selection.choices), 1),
        )
    end
    on(text -> isnothing(text) || go_to(text), textbox.stored_string)
    on(_ -> step(-1), previous.clicks)
    on(_ -> step(1), next.clicks)
    foreach(box -> on(_ -> reselect(), box.checked), values(species_boxes))
    isempty(orthogroup) ? notify(current) : go_to(orthogroup)
    DataInspector(figure)
    isnothing(screen) || display(screen, figure)
    return (;
        figure,
        screen,
        textbox,
        previous,
        next,
        aligned,
        status,
        current,
        species_boxes,
    )
end
end

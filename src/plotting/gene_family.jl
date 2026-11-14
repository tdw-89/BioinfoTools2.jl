#= Gene-family phylograms: a DivergenceTree drawn with dS as distance, plus
   per-gene heatmap tracks. Included into `Plotting`. =#

"""
Given a preorder tree's `subtree_end`, place its nodes vertically for a
rectangular phylogram. Returns a `Vector{Float64}`: leaves at 1, 2, … in
preorder, each internal node midway between its first and last child.
"""
function phylogram_layout(subtree_end::AbstractVector{<:Integer})
    y = zeros(Float64, length(subtree_end))
    n_leaves = 0
    for node in eachindex(subtree_end)
        subtree_end[node] == node && (y[node] = (n_leaves += 1))
    end
    for node in reverse(eachindex(subtree_end))
        subtree_end[node] == node && continue
        last_child = node + 1
        while subtree_end[last_child] < subtree_end[node]
            last_child = subtree_end[last_child] + 1
        end
        y[node] = (y[node+1] + y[last_child]) / 2
    end
    return y
end

"""
Given a preorder tree's `parent`, per-branch `lengths` and node positions `x`,
`y`, list a rectangular phylogram's segments. Returns `(solid, dashed)` point
vectors for `linesegments!`: each branch horizontal at its node's `y`, dashed
past `saturation`, and each parent's vertical connector, solid.
"""
function phylogram_segments(parent, lengths, x, y; saturation::Real = Inf)
    solid, dashed = Point2f[], Point2f[]
    low, high = copy(y), copy(y)
    for node in reverse(eachindex(parent))
        up = parent[node]
        up == 0 && continue
        target = lengths[node] > saturation ? dashed : solid
        push!(target, Point2f(x[up], y[node]), Point2f(x[node], y[node]))
        low[up], high[up] = min(low[up], y[node]), max(high[up], y[node])
    end
    for node in eachindex(parent)
        low[node] < high[node] &&
            push!(solid, Point2f(x[node], low[node]), Point2f(x[node], high[node]))
    end
    return solid, dashed
end

"""The ten categorical colours [`species_palette`](@ref) cycles through."""
const _TAB10 = CairoMakie.Makie.to_colormap(:tab10)

"""
Given a clade (or species names), give each species a fixed colour from the
categorical `:tab10` scheme, cycling past ten. Returns a `name => color` `Dict`,
the same for every tree of the clade.
"""
species_palette(names::AbstractVector{<:AbstractString}) =
    Dict(name => _TAB10[mod1(index, length(_TAB10))] for (index, name) in enumerate(names))

species_palette(clade::Orthologs.Clade) = species_palette(Orthologs.species_names(clade))

"""How each node event is drawn: `(label, marker, fill, stroke)`."""
const EVENT_STYLES = (
    duplication = ("duplication", :diamond, :firebrick, :firebrick),
    speciation = ("speciation", :circle, :white, :black),
    mixed = ("mixed polytomy", :rect, :darkorange, :darkorange),
    unresolved = ("unresolved pair", :xcross, :gray40, :gray40),
)

"""
Given a divergence tree leaf, its label: the gene ID less NCBI's `gene-` prefix,
else the protein accession.
"""
function default_leaf_label(tree::Orthologs.DivergenceTree, leaf::Integer)
    gene = Orthologs.leaf_gene(tree, leaf)
    return ismissing(gene) ? Orthologs.leaf_protein(tree, leaf) :
           String(chopprefix(gene, "gene-"))
end

"""
One heatmap column of [`gene_family_plot!`](@ref): per species name, a `gene ID
or protein accession => value` map (`values`), drawn with its own `colormap` and
`colorrange`. Tracks sharing `label`, colormap and range share a colorbar.
`notes` (the same shape, holding strings) are added to the leaves' tooltips, and
`opacity` (the same shape, in `[0, 1]`) fades each cell; a gene without one is opaque.
A cell without a value is [`MISSING_COLOR`](@ref), unless `why_missing` (the same
shape, holding strings) gives it a reason that `missing_colors` (`reason => color`
pairs) colours. A track with `missing_colors` lists them in the legend, after
[`MISSING_COLOR`](@ref) as `missing_label`.
"""
struct GeneTrack
    name::String
    values::Dict{String,Dict{String,Float64}}
    colormap::Any
    colorrange::Tuple{Float64,Float64}
    label::String
    notes::Dict{String,Dict{String,String}}
    opacity::Dict{String,Dict{String,Float64}}
    why_missing::Dict{String,Dict{String,String}}
    missing_colors::Vector{Pair{String,Any}}
    missing_label::String
end

"""
    GeneTrack(name, values; colormap = :viridis, colorrange = (0, 1), label = name,
              notes = Dict(), opacity = Dict(), why_missing = Dict(),
              missing_colors = [], missing_label = "no value")

Given a column name and its `species => (key => value)` map, build a
[`GeneTrack`](@ref). Returns it, the maps copied with `String` keys.
"""
GeneTrack(
    name::AbstractString,
    values::AbstractDict;
    colormap = :viridis,
    colorrange = (0.0, 1.0),
    label::AbstractString = name,
    notes::AbstractDict = Dict{String,Dict{String,String}}(),
    opacity::AbstractDict = Dict{String,Dict{String,Float64}}(),
    why_missing::AbstractDict = Dict{String,Dict{String,String}}(),
    missing_colors = Pair{String,Any}[],
    missing_label::AbstractString = "no value",
) = GeneTrack(
    String(name),
    _nested(Float64, values),
    colormap,
    Float64.(Tuple(colorrange)),
    String(label),
    _nested(String, notes),
    _nested(Float64, opacity),
    _nested(String, why_missing),
    Pair{String,Any}[String(reason) => color for (reason, color) in missing_colors],
    String(missing_label),
)

"""Given a `species => (key => x)` map, copy it with `String` keys and `T` values."""
_nested(::Type{T}, maps::AbstractDict) where {T} = Dict{String,Dict{String,T}}(
    String(species) => Dict{String,T}(String(key) => T(x) for (key, x) in by_key) for
    (species, by_key) in maps
)

"""
Given a track, or a `name => (species => (key => value))` pair taking the plot's
defaults, return it as a [`GeneTrack`](@ref).
"""
_track(track::GeneTrack; _...) = track
_track(track::Pair; colormap = :viridis, colorrange = (0.0, 1.0), label = first(track)) =
    GeneTrack(first(track), last(track); colormap, colorrange, label)

"""
Given a `species => (key => x)` map and a leaf, look the leaf up by species, then
gene ID, then protein accession. Returns its entry, or `default`.
"""
function _leaf_entry(maps::AbstractDict, tree::Orthologs.DivergenceTree, leaf, default)
    by_key = get(maps, Orthologs.leaf_species(tree, leaf), nothing)
    isnothing(by_key) && return default
    gene = Orthologs.leaf_gene(tree, leaf)
    !ismissing(gene) && haskey(by_key, gene) && return by_key[gene]
    return get(by_key, Orthologs.leaf_protein(tree, leaf), default)
end

"""
Given a tree's leaves and `tracks` ([`GeneTrack`](@ref)s or `name => (species
=> (key => value))` pairs), look each leaf's values up, by gene ID or else
protein accession. Returns a `n_tracks × n_leaves` matrix, `NaN` where missing.
"""
function track_values(tree::Orthologs.DivergenceTree, leaves, tracks)
    values = fill(NaN, length(tracks), length(leaves))
    for (row, track) in enumerate(_track.(tracks)), (column, leaf) in enumerate(leaves)
        values[row, column] = _leaf_entry(track.values, tree, leaf, NaN)
    end
    return values
end

"""Given a number, format it for a tooltip: 3 significant digits, `n/a` for `NaN`."""
_tip(value::Real) = isnan(value) ? "n/a" : string(round(value; sigdigits = 3))

"""
Given a divergence tree node, describe it for its tooltip: gene, protein,
species and `track_rows` (`(name, value, note)`) at a leaf; event, species-tree node,
genes below and dS age (see `node_ages`) elsewhere; then its branch dS/dN and its
`x` (dS) from the root.
"""
function node_tooltip(tree::Orthologs.DivergenceTree, node::Integer, x, track_rows = ())
    lines = String[]
    if Orthologs.is_leaf(tree, node)
        gene = Orthologs.leaf_gene(tree, node)
        push!(lines, ismissing(gene) ? "(no gene)" : gene)
        protein, species =
            Orthologs.leaf_protein(tree, node), Orthologs.leaf_species(tree, node)
        push!(lines, "$protein · $species")
        for (name, value, note) in track_rows
            push!(
                lines,
                isempty(note) ? "$name: $(_tip(value))" : "$name: $(_tip(value)) ($note)",
            )
        end
    else
        event = first(getproperty(EVENT_STYLES, Orthologs.node_event(tree, node)))
        origin = tree.exact[node] ? "OrthoFinder node" : "codeml-only split"
        push!(lines, "$event at $(Orthologs.species_node(tree, node)) ($origin)")
        n_genes = length(Orthologs.leaves(tree, node))
        push!(
            lines,
            "$n_genes genes below; age (mean dS to them) $(_tip(tree.dS_age[node]))",
        )
    end
    dS = Orthologs.branch_length(tree, node)
    dN = Orthologs.branch_length(tree, node; rate = :dN)
    push!(lines, "branch dS $(_tip(dS)), dN $(_tip(dN)); from root $(_tip(x[node]))")
    return join(lines, '\n')
end

"""
Given a tree axis, draw a divergence tree's branches (dashed past `saturation`)
and its stem when pruning left one. Returns `nothing`.
"""
function _draw_branches!(axis, tree, lengths, x, y; saturation)
    solid, dashed = phylogram_segments(tree.parent, lengths, x, y; saturation)
    stem = lengths[1]
    isfinite(stem) && stem > 0 && push!(solid, Point2f(-stem, y[1]), Point2f(0, y[1]))
    linesegments!(axis, solid; color = :gray25, linewidth = 1.2, inspectable = false)
    isempty(dashed) || linesegments!(
        axis,
        dashed;
        color = :gray25,
        linewidth = 1.2,
        linestyle = :dash,
        inspectable = false,
    )
    return nothing
end

"""The colormap for relative branch rates: slow blue, average grey, fast red."""
const RATE_COLORMAP = cgrad([:steelblue, :gray35, :firebrick])

"""The relative rates a chronogram colours, as log2: a quarter to four times the mean."""
const RATE_RANGE = (-2.0, 2.0)

"""
Given a branch's length and the age it spans, its rate relative to the mean the
ages assume, as log2 clamped to [`RATE_RANGE`](@ref). Returns it: 0 for a
zero-length branch over no time, the fastest for length over none, `NaN` if unknown.
"""
function relative_rate(length::Real, span::Real)
    (isfinite(length) && isfinite(span)) || return NaN
    if span <= eps()
        return length <= eps() ? 0.0 : last(RATE_RANGE)
    end
    length <= 0 && return first(RATE_RANGE)
    return clamp(log2(length / span), RATE_RANGE...)
end

"""
Given a tree axis, draw a divergence tree as a chronogram: each branch from its
parent's age to its own (`ages`), coloured by [`relative_rate`](@ref) and dashed
past `saturation`, grey connectors, and a short stub for a pruned tree's stem.
Returns `nothing`.
"""
function _draw_chronogram_branches!(axis, tree, lengths, ages, y; saturation)
    branches = (solid = (Point2f[], Float64[]), dashed = (Point2f[], Float64[]))
    connectors = Point2f[]
    low, high = copy(y), copy(y)
    for node in reverse(eachindex(tree.parent))
        up = tree.parent[node]
        up == 0 && continue
        points, rates = lengths[node] > saturation ? branches.dashed : branches.solid
        push!(points, Point2f(ages[up], y[node]), Point2f(ages[node], y[node]))
        append!(rates, fill(relative_rate(lengths[node], ages[up] - ages[node]), 2))
        low[up], high[up] = min(low[up], y[node]), max(high[up], y[node])
    end
    for node in eachindex(tree.parent)
        low[node] < high[node] && push!(
            connectors,
            Point2f(ages[node], low[node]),
            Point2f(ages[node], high[node]),
        )
    end
    stem = lengths[1]
    if isfinite(stem) && stem > 0
        stub = max(ages[1], eps()) * 0.06
        push!(connectors, Point2f(ages[1] + stub, y[1]), Point2f(ages[1], y[1]))
    end
    linesegments!(axis, connectors; color = :gray55, linewidth = 1, inspectable = false)
    for (style, (points, rates)) in pairs(branches)
        isempty(points) && continue
        linesegments!(
            axis,
            points;
            color = rates,
            colormap = RATE_COLORMAP,
            colorrange = RATE_RANGE,
            nan_color = :gray70,
            linewidth = 1.6,
            linestyle = style === :dashed ? :dash : :solid,
            inspectable = false,
        )
    end
    return nothing
end

"""
Given a tree axis, mark a divergence tree's internal nodes by event, faded where
codeml split what OrthoFinder left a polytomy. Returns the `Scatter` plots.
"""
function _draw_events!(axis, tree, x, y, tooltips)
    plots = Any[]
    for (event, (_, marker, fill, stroke)) in pairs(EVENT_STYLES)
        nodes = findall(
            node -> tree.subtree_end[node] != node && tree.event[node] == event,
            eachindex(tree.parent),
        )
        isempty(nodes) && continue
        fade = [tree.exact[node] ? 1.0 : 0.35 for node in nodes]
        plot = scatter!(
            axis,
            x[nodes],
            y[nodes];
            marker,
            markersize = 11,
            color = [(fill, alpha) for alpha in fade],
            strokecolor = [(stroke, alpha) for alpha in fade],
            strokewidth = 1.2,
            inspector_label = (_, index, _) -> tooltips[nodes[index]],
        )
        push!(plots, plot)
    end
    return plots
end

"""
Given a grid cell, draw [`GeneTrack`](@ref)s as heatmap columns beside the
leaves, each on its own scale, with the leaf `labels` as right-hand ticks.
Returns the `Axis` and the track values ([`track_values`](@ref)).
"""
function _draw_tracks!(cell, tree, leaves, tracks::Vector{GeneTrack}, labels)
    n_tracks, n_leaves = length(tracks), length(leaves)
    axis = Axis(
        cell;
        width = max(18 * n_tracks, 1),
        yreversed = true,
        yaxisposition = :right,
        yticks = (1:n_leaves, labels),
        yticklabelsize = 10,
        xticks = (1:max(n_tracks, 1), n_tracks == 0 ? [""] : [t.name for t in tracks]),
        xticklabelrotation = π / 3,
        xgridvisible = false,
        ygridvisible = false,
    )
    hidespines!(axis)
    values = track_values(tree, leaves, tracks)
    for (row, track) in enumerate(tracks)
        reasons = [
            isfinite(value) ? "" : _leaf_entry(track.why_missing, tree, leaf, "") for
            (leaf, value) in zip(leaves, values[row, :])
        ]
        heatmap!(
            axis,
            [row - 0.5, row + 0.5],
            collect(0.5:1:(n_leaves+0.5)),
            reshape(_track_colors(track, tree, leaves, values[row, :]), 1, n_leaves);
            inspector_label = (_, (_, column), _) -> string(
                "$(track.name): $(_tip(values[row, column]))",
                isempty(reasons[column]) ? "" : " ($(reasons[column]))",
                "\n$(labels[column])",
            ),
        )
    end
    n_tracks == 0 || xlims!(axis, 0.5, n_tracks + 0.5)
    return axis, values
end

"""The colour of a track cell with no value."""
const MISSING_COLOR = CairoMakie.Makie.RGBAf(0.88, 0.88, 0.88, 1)

"""
Given a track, its leaves and their `values`, colour each cell from the track's
colormap (clamped to its range) at the leaf's `opacity`. Returns `RGBAf`s; a
missing value takes its reason's `missing_colors` entry, else [`MISSING_COLOR`](@ref).
"""
function _track_colors(track::GeneTrack, tree, leaves, values)
    colormap = CairoMakie.Makie.to_colormap(track.colormap)
    low, high = track.colorrange
    reason_colors = Dict(
        reason => CairoMakie.Makie.to_color(color) for (reason, color) in track.missing_colors
    )
    return map(leaves, values) do leaf, value
        if !isfinite(value)
            reason = _leaf_entry(track.why_missing, tree, leaf, "")
            return CairoMakie.Makie.RGBAf(get(reason_colors, reason, MISSING_COLOR))
        end
        color = CairoMakie.Makie.interpolated_getindex(
            colormap,
            clamp(value, low, high),
            (low, high),
        )
        alpha = clamp(_leaf_entry(track.opacity, tree, leaf, 1.0), 0.0, 1.0)
        CairoMakie.Makie.RGBAf(color.r, color.g, color.b, alpha)
    end
end

"""
Given a grid cell and tracks, add one `Colorbar` per distinct `(label, colormap,
colorrange)`, side by side, then the branch-rate scale when `rates`. Returns them,
empty when there are none.
"""
function _colorbars!(cell, tracks::Vector{GeneTrack}; rates::Bool)
    scales = Any[
        (; label, colormap, limits, ticks = CairoMakie.Makie.automatic) for
        (label, colormap, limits) in
        unique((track.label, track.colormap, track.colorrange) for track in tracks)
    ]
    rates && push!(
        scales,
        (
            label = "branch rate (relative to the ages)",
            colormap = RATE_COLORMAP,
            limits = RATE_RANGE,
            ticks = (-2:2, ["¼×", "½×", "1×", "2×", "4×"]),
        ),
    )
    isempty(scales) && return Colorbar[]
    grid = GridLayout(cell)
    return [
        Colorbar(
            grid[1, column];
            scale.colormap,
            scale.limits,
            scale.label,
            scale.ticks,
            height = Relative(0.6),
            labelsize = 12,
        ) for (column, scale) in enumerate(scales)
    ]
end

"""
Given tracks, the `(label, color)` legend entries of those with `missing_colors`:
each one's `missing_label` for [`MISSING_COLOR`](@ref), then its reasons, without repeats.
"""
_missing_entries(tracks::Vector{GeneTrack}) = unique(
    entry for track in tracks if !isempty(track.missing_colors) for entry in [
        (track.missing_label, MISSING_COLOR)
        [(reason, CairoMakie.Makie.to_color(color)) for (reason, color) in track.missing_colors]
    ]
)

"""
Given a grid cell, the species and events present, whether the tree has
codeml-only splits (`faded`) or saturated branches (`dashed`), and the colours of
missing track cells (`missing_cells`, `(label, color)`s), add its legend.
Returns the `Legend`.
"""
function _family_legend!(
    cell,
    colors,
    species,
    events;
    faded,
    dashed,
    saturation,
    missing_cells = (),
)
    species_entries = [
        MarkerElement(marker = :circle, color = colors[name], markersize = 11) for
        name in species
    ]
    event_entries, event_labels = Any[], String[]
    for event in events
        label, marker, fill, stroke = getproperty(EVENT_STYLES, event)
        push!(
            event_entries,
            MarkerElement(; marker, color = fill, strokecolor = stroke, strokewidth = 1.2),
        )
        push!(event_labels, label)
    end
    if faded
        push!(event_entries, MarkerElement(marker = :diamond, color = (:firebrick, 0.35)))
        push!(event_labels, "codeml-only split (faded)")
    end
    if dashed
        push!(event_entries, LineElement(color = :gray25, linestyle = :dash))
        push!(event_labels, "branch > $saturation")
    end
    groups, labels = [species_entries, event_entries], [collect(species), event_labels]
    titles = ["Species", "Nodes"]
    if !isempty(missing_cells)
        groups = [groups..., [PolyElement(color = color) for (_, color) in missing_cells]]
        labels = [labels..., [label for (label, _) in missing_cells]]
        push!(titles, "Missing track values")
    end
    return Legend(
        cell,
        groups,
        labels,
        titles;
        framevisible = false,
        labelsize = 11,
        titlesize = 12,
        rowgap = 1,
    )
end

"""
    gene_family_plot!(position, tree; layout = :chronogram, tracks = [], rate = :dS,
                      saturation = 2.0, colormap = :viridis, colorrange = (0, 1),
                      track_label, species_colors = species_palette(tree.clade),
                      max_labels = 80, leaf_label = default_leaf_label, title)

Given a grid position (or a `GridLayout` to fill) and a [`DivergenceTree`](@ref),
draw it with heatmap `tracks` beside its leaves: [`GeneTrack`](@ref)s, or `name =>
(species => (key => value))` pairs drawn with `colormap`, `colorrange` and
`track_label`. `layout` `:chronogram` aligns the genes at age 0, nodes at their
`rate` age ([`node_ages`](@ref)) and branches coloured by [`relative_rate`](@ref);
`:phylogram` puts x at the `rate` from the root. Returns `(tree_axis, track_axis,
legend, colorbars, positions)`, `positions` the nodes' `(x, y)`. Labels are hidden
past `max_labels` leaves; tooltips need a `DataInspector`.
"""
function gene_family_plot!(
    position,
    tree::Orthologs.DivergenceTree;
    layout::Symbol = :chronogram,
    tracks = Pair{String,Dict{String,Dict{String,Float64}}}[],
    rate::Symbol = :dS,
    saturation::Real = 2.0,
    colormap = :viridis,
    colorrange = (0.0, 1.0),
    track_label = "promoter methylation",
    species_colors = species_palette(tree.clade),
    max_labels::Integer = 80,
    leaf_label = default_leaf_label,
    title = nothing,
)
    layout in (:chronogram, :phylogram) ||
        throw(ArgumentError("`layout` must be :chronogram or :phylogram (got :$layout)"))
    chronogram = layout === :chronogram
    grid = position isa GridLayout ? position : GridLayout(position)
    lengths = rate === :dS ? tree.dS : tree.dN
    from_root, y = Orthologs.root_distances(tree), phylogram_layout(tree.subtree_end)
    x = chronogram ? Orthologs.node_ages(tree; rate) : Orthologs.root_distances(tree; rate)
    leaves = Orthologs.leaves(tree)
    labels =
        length(leaves) <= max_labels ? [leaf_label(tree, leaf) for leaf in leaves] :
        fill("", length(leaves))
    tracks = GeneTrack[
        _track(track; colormap, colorrange, label = track_label) for track in tracks
    ]
    track_axis, values = _draw_tracks!(grid[1, 2], tree, leaves, tracks, labels)
    track_rows(column) = [
        (
            track.name,
            values[row, column],
            _leaf_entry(track.notes, tree, leaves[column], ""),
        ) for (row, track) in enumerate(tracks)
    ]
    tooltips = [node_tooltip(tree, node, from_root) for node in eachindex(tree.parent)]
    tooltips[leaves] .= [
        node_tooltip(tree, leaf, from_root, track_rows(column)) for
        (column, leaf) in enumerate(leaves)
    ]

    tree_axis = Axis(
        grid[1, 1];
        yreversed = true,
        xreversed = chronogram,
        xlabel = chronogram ? "age: mean $rate down to the genes  (older → now)" :
                 "$rate from the root  (older → younger)",
        title = something(title, _family_title(tree)),
        titlealign = :left,
        ygridvisible = false,
    )
    hideydecorations!(tree_axis)
    hidespines!(tree_axis, :l, :r, :t)
    if chronogram
        _draw_chronogram_branches!(tree_axis, tree, lengths, x, y; saturation)
    else
        _draw_branches!(tree_axis, tree, lengths, x, y; saturation)
    end
    _draw_events!(tree_axis, tree, x, y, tooltips)
    species = [Orthologs.leaf_species(tree, leaf) for leaf in leaves]
    scatter!(
        tree_axis,
        x[leaves],
        y[leaves];
        color = [species_colors[name] for name in species],
        markersize = 9,
        strokecolor = :white,
        strokewidth = 0.5,
        inspector_label = (_, index, _) -> tooltips[leaves[index]],
    )
    linkyaxes!(tree_axis, track_axis)
    ylims!(tree_axis, length(leaves) + 0.5, 0.5)
    chronogram || _clip_stem!(tree_axis, lengths[1], maximum(x))
    colorbars = _colorbars!(grid[1, 3], tracks; rates = chronogram)
    present = filter(name -> name in species, Orthologs.species_names(tree.clade))
    events = filter(event -> event in tree.event, collect(keys(EVENT_STYLES)))
    faded, dashed = !all(tree.exact), any(>(saturation), filter(isfinite, lengths))
    legend_cell = grid[1, isempty(colorbars) ? 3 : 4]
    legend = _family_legend!(
        legend_cell,
        species_colors,
        present,
        events;
        faded,
        dashed,
        saturation,
        missing_cells = _missing_entries(tracks),
    )
    colgap!(grid, 1, 4)
    return (; tree_axis, track_axis, legend, colorbars, positions = (x, y))
end

"""
Given the tree axis, a pruned tree's stem and its depth, start the x axis at most
a tenth of the depth left of the root, so a long stem runs off the edge rather
than squeezing the tree. Returns `nothing`; its full length is in the root's tooltip.
"""
function _clip_stem!(axis, stem::Real, depth::Real)
    (isfinite(stem) && stem > 0 && depth > 0) || return nothing
    stem > depth / 10 && xlims!(axis, -depth / 10, depth * 1.03)
    return nothing
end

"""Given a divergence tree, its default title: orthogroup, gene and species counts."""
function _family_title(tree::Orthologs.DivergenceTree)
    leaves = Orthologs.leaves(tree)
    n_species = length(unique(Orthologs.leaf_species(tree, leaf) for leaf in leaves))
    dropped = isempty(tree.dropped) ? "" : ", $(length(tree.dropped)) not in codeml's tree"
    return "$(Orthologs.orthogroup(tree)): $(length(leaves)) genes, $n_species species$dropped"
end

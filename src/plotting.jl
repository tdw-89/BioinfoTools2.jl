"""
Figures for the package's data types: metagene line profiles and heatmaps, and
per-quantile distributions.

Every function draws into a caller-supplied grid position, so figure size and
layout stay the caller's.
"""
module Plotting

using CairoMakie

"""
Given a flank width and the number of body columns, place a metagene's four
landmark ticks. Returns `(positions, labels)` — region start, TSS, TES and
region end — ready for an axis's `xticks`.
"""
metagene_xticks(flank::Integer, body_columns::Integer) = (
    [1.0, flank + 0.5, flank + body_columns + 0.5, Float64(2 * flank + body_columns)],
    ["-$(flank) bp", "TSS", "TES", "+$(flank) bp"],
)

"""
Given a flank width and the number of body columns, locate the TSS and TES.
Returns the two positions, each between a flank column and a body column.
"""
metagene_boundaries(flank::Integer, body_columns::Integer) =
    [flank + 0.5, flank + body_columns + 0.5]

"""
Given a metagene matrix (profiles along its rows), widen its body from
`body_bins` to `display_bins` columns, leaving the `flank`-wide flanks untouched.
Returns the stretched matrix — or `matrix` itself when the widths already agree.

**NOTE:** columns are repeated, not interpolated: a bin is drawn wider, never
invented, so the true resolution stays `body_bins`.
"""
function stretch_body_columns(
    matrix::AbstractMatrix,
    flank::Integer,
    body_bins::Integer,
    display_bins::Integer,
)
    display_bins == body_bins && return matrix
    n_columns = size(matrix, 2)
    n_columns == 2 * flank + body_bins || throw(
        DimensionMismatch("expected $(2 * flank + body_bins) columns, got $n_columns"),
    )

    stretched = similar(matrix, size(matrix, 1), 2 * flank + display_bins)
    stretched[:, 1:flank] .= view(matrix, :, 1:flank)
    stretched[:, (end-flank+1):end] .= view(matrix, :, (n_columns-flank+1):n_columns)
    for slot = 1:display_bins
        source = clamp(ceil(Int, slot * body_bins / display_bins), 1, body_bins)
        stretched[:, flank+slot] .= view(matrix, :, flank + source)
    end
    return stretched
end

"""
Given a metagene profile and its `flank`, draw it as a line with the TSS and TES
dashed in (unless `mark_boundaries = false`). Returns the new `Axis`, placed at
`position` (e.g. `figure[1, 1]`); the body width is read off the profile's
length.
"""
function metagene_lines!(
    position,
    profile::AbstractVector;
    flank::Integer,
    title::AbstractString = "",
    ylabel::AbstractString = "",
    mark_boundaries::Bool = true,
)
    body_columns = length(profile) - 2 * flank
    axis = Axis(
        position;
        title = title,
        ylabel = ylabel,
        xticks = metagene_xticks(flank, body_columns),
    )
    lines!(axis, eachindex(profile), profile)
    mark_boundaries && vlines!(
        axis,
        metagene_boundaries(flank, body_columns);
        color = :gray,
        linestyle = :dash,
    )
    return axis
end

"""
Given a `groups × positions` metagene matrix (as `Exploration.quantile_profiles`
returns, perhaps through [`stretch_body_columns`](@ref)) and its `flank`, draw it
as a heatmap, groups up the y axis. Returns `(axis, plot)` placed at `position`;
pass `plot` to a `Colorbar` in a neighbouring cell.

`NaN` cells — unmeasured positions — are left blank. Set `mark_boundaries =
false` to leave off the dashed TSS/TES lines.
"""
function metagene_heatmap!(
    position,
    matrix::AbstractMatrix;
    flank::Integer,
    title::AbstractString = "",
    ylabel::AbstractString = "",
    colormap = :viridis,
    colorrange = CairoMakie.Makie.automatic,
    mark_boundaries::Bool = true,
)
    n_groups, n_columns = size(matrix)
    body_columns = n_columns - 2 * flank
    axis = Axis(
        position;
        title = title,
        ylabel = ylabel,
        xticks = metagene_xticks(flank, body_columns),
        yticks = 1:n_groups,
    )
    plot = heatmap!(
        axis,
        1:n_columns,
        1:n_groups,
        permutedims(matrix);
        colormap = colormap,
        colorrange = colorrange,
    )
    mark_boundaries && vlines!(
        axis,
        metagene_boundaries(flank, body_columns);
        color = :white,
        linestyle = :dash,
    )
    return axis, plot
end

"""
Given `values` grouped by integer `bins` (as `Exploration.values_by_quantile`
returns), draw each group as a violin with a narrow boxplot over it, keeping the
median and quartiles legible. Returns the new `Axis`, placed at `position`.
"""
function violin_box!(
    position,
    bins::AbstractVector{<:Integer},
    values::AbstractVector{<:Real};
    title::AbstractString = "",
    xlabel::AbstractString = "quantile",
    ylabel::AbstractString = "",
)
    axis = Axis(
        position;
        title = title,
        xlabel = xlabel,
        ylabel = ylabel,
        xticks = 1:maximum(bins; init = 0),
    )
    isempty(bins) && return axis
    violin!(
        axis,
        bins,
        values;
        width = 0.8,
        color = (:steelblue, 0.4),
        strokecolor = :steelblue,
        strokewidth = 1,
    )
    boxplot!(
        axis,
        bins,
        values;
        width = 0.2,
        color = :white,
        strokecolor = :black,
        strokewidth = 1,
        show_outliers = false,
    )
    return axis
end

export metagene_xticks,
    metagene_boundaries,
    stretch_body_columns,
    metagene_lines!,
    metagene_heatmap!,
    violin_box!

end

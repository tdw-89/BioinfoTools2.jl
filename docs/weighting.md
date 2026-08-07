# Weighting in methylation profiles

Where read depth becomes a statistical weight, and why each site does what it
does. All code is in `src/exploration.jl`.

There are **three** places a methylation level gets averaged. Only the last two
are user-controllable; the first is deliberately fixed.

---

## 1. Pooling reads at one cytosine — always linear in depth

`feature_frequency(genome, feature, ::MethylationData)`

```julia
fraction = meth_fraction(call)
isnan(fraction) && continue
# Always linear in depth: this pools reads at one cytosine
# back into that cytosine's own methylation fraction, which
# only comes out right weighted by read count.
weighted_fraction += fraction * Int(depth)
depth_total += Int(depth)
```

Two calls can land on one base — the same position reported on opposite strands,
or two datasets combined by `merge_calls`. Pooling them is not a modelling
choice: `fraction × depth` recovers the methylated read count, so the result is
that cytosine's true fraction over its combined reads. A transform here would
produce a number that is not the methylation fraction of anything, so
`weight_transform` is not applied at this step.

The stored `FeatureLevels.weights` are therefore raw summed depths, which keeps
a `MethylationFrequency` independent of any weighting scheme — you can re-profile
the same object under `identity`, `sqrt` and `log` without rebuilding it.

---

## 2. Pooling bases into a body bin — `weight_transform`

`gene_profile(::FeatureLevels)`

```julia
depth = nonzero_weights[entry]
depth == 0 && continue
weight = weight_by_depth ? weight_of(depth, weight_transform) : 1.0
weight == 0 && continue
...
weighted_sums[slot] += Float64(feature_levels.levels[base]) * weight
profile_weights[slot] += weight
```

The body is binned rather than interpolated, so several cytosines usually fall
in one bin. Each contributes `weight_transform(depth)`. `identity` (the default)
weights linearly by depth; `sqrt` or `log` compress the advantage of deeply
covered sites, which matters when coverage varies by an order of magnitude
across a library.

`weight_by_depth = false` gives every measured base a weight of 1 and ignores
`weight_transform`.

The returned `weights` vector doubles as the coverage mask: **a nonzero weight,
not a nonzero level, is what marks a position as measured**, since a base
measured at 0% methylation is a structural zero in `levels`.

---

## 3. Averaging genes into a metagene — same transform, applied once

`mean_gene_profile(::MethylationFrequency)`

```julia
# `gene_profile` already applied `weight_transform` to each base, so
# this sum is used as-is: transforming it again would give
# `f(Σ f(depth))`, which is not a weight on any observation.
profile.weights[slot] == 0 && continue
weight = weight_by_depth ? profile.weights[slot] : 1.0
weighted_sums[slot] += profile.levels[slot] * weight
weight_totals[slot] += weight
```

A gene contributes to a position only if it was measured there, so a position is
the mean level among genes with data rather than being dragged toward zero by
genes with no nearby cytosine. Positions no gene measured return `NaN`.

**Depth does drive this average.** A gene's weight at a position is the slot
weight step 2 computed, which for a flank position — one base per slot — is
exactly `weight_transform(depth)` of that gene's cytosine there. Two genes
measured at the same offset at 40× and 5× contribute in a 40:5 ratio under
`identity`, and in a compressed ratio under `sqrt` or `log`.

What this step does *not* do is apply the transform a second time. Every raw
depth becomes a weight exactly once, in step 2; everything above that aggregates
linearly. Re-transforming here would give `f(Σ f(depth))`, a transform of
something that is no longer a depth.

One consequence to be aware of: for a *body bin*, which can hold several
cytosines, the slot weight is their **sum**, so at equal depth a CpG-dense gene
outweighs a sparse one in proportion to its cytosine count. That is a
deliberate "more measurements, more information" reading, but it does let
cytosine density back into the weighting after being excluded from the value.
Use the mean rather than the sum if that is not wanted.

---

## Validation

`weight_of` guards the transform:

```julia
@inline function weight_of(depth::Real, weight_transform)
    weight = Float64(weight_transform(Float64(depth)))
    (isfinite(weight) && weight >= 0) || throw(
        ArgumentError(
            "`weight_transform` returned $weight for depth $depth; weights must be finite and non-negative",
        ),
    )
    return weight
end
```

Depth is converted to `Float64` first because depths are stored unsigned, where
`-depth` wraps to a huge positive number rather than going negative — a bad
transform would pass unnoticed.

Two gotchas worth knowing:

- `log(1) == 0`, so a depth-1 site under `log` gets weight 0 and drops out
  entirely. Prefer `log1p` unless that is what you want. With the default
  `min_depth = 5` it cannot arise.
- The transform must be monotonic for the result to be interpretable as a
  weighting; that is not checkable, so it is on the caller.

---

## Summary

| Step | Function | Weight | Configurable |
|---|---|---|---|
| Reads → one cytosine | `feature_frequency` | read depth, linear | no, by design |
| Cytosines → flank base or body bin | `gene_profile` | `weight_transform(depth)` per base | yes |
| Genes → metagene | `mean_gene_profile` | that slot's weight — `weight_transform(depth)` for a flank base, the sum over the gene's cytosines for a body bin | yes, via the same transform |

Defaults reproduce plain linear depth weighting at every step.

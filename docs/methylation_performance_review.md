# Methylation performance review

A review of `Data.Methylation` (`src/data/methylation.jl`), its consumers in
`Exploration` (`feature_frequency`/`gene_profile`/`_group_profiles`), and how the
analysis notebooks use them (`zebrafish_and_stickleback`:
`methylation_vs_ds.qmd`, `expression_vs_methylation.qmd`, `tissue_comparison.qmd`;
`primates`: `methylation_vs_ds.qmd`). Scope is strictly *optimization* — nothing
below changes what any analysis computes unless explicitly flagged.

## Status (implemented)

Measured on real zebrafish data (3 × `.cov.gz`, 8 threads), outputs hashed
against the pre-change code:

| Step | Before | After |
|---|---|---|
| `feature_frequency` (methylation) | 1.02 s | 0.16 s |
| `quantile_profiles` (methylation) | 3.55 s | 0.33 s |
| `tss_window` (methylation) | 0.37 s | 0.02 s |
| `merge_calls`, 3 datasets | 0.60 s | 0.14 s |
| `load_bismark_cov`, 3 files, cached | 3.36 s | 0.95 s |

All 52 813 probes are bit-identical except `mean_gene_profile`, which moves by
≤ 1.1e-16: its per-gene sum now runs in a different `Dict` order.

- **Done:** §1 shared-index `FeatureLevels`; §2 as `load_bismark_cov(...; cache)`;
  §4 per-scaffold merge, plus the k-way merge as opt-in `single_rounding`
  (sort-based); §5 threaded `feature_frequency` (both methods), presized
  `_region_levels`; §6 presized `aggregate_keys!`, column-indexed
  `merge_scaffold`; §7 `materialize` (17 % faster `feature_frequency` on
  materialized columns); §9 notebook items.
- **Found beyond the review:** `gene_profile` did not specialize on
  `weight_transform`, boxing every weight — most of the `quantile_profiles` gain.
- **Not done:** `SortingAlgorithms.RadixSort` — Base already radix-sorts `UInt64`
  (20 M keys in 0.17 s); porting `load_bismark` to the chunk machinery — deferred
  by this review itself until extractor files are used at scale; §3 and §8 carry
  no code changes.

## Benchmark setup

Apple M3 (8 cores, 4P+4E), **16 GB** unified memory, Julia 1.13, real inputs:

- `primates/.../chimp/merged.cov` — 1.3 GB plain text, 32.5 M sites, 26 scaffolds.
- one zebrafish `.cov.gz` — 232 MB gzipped (~1.1 GB inflated), 40.3 M sites.
- The zebrafish set is **21 such files (4.8 GB gzipped)**, and each notebook
  render parses the whole set **twice** (once per tissue group, once pooled).

Measured (hot, after compile; `@time`):

| Operation | 8 threads | 1 thread |
|---|---|---|
| `load_bismark_cov`, plain 1.3 GB `.cov` | **0.48 s** | 1.96 s |
| `load_bismark_cov`, one 232 MB `.cov.gz` | **2.43 s** | 4.29 s |
| `gzcat` (Apple zlib) of that `.gz`, baseline | 0.98 s | — |
| `merge_calls`, one pairwise whole-genome merge | 0.23 s (serial) | — |
| `write_methylation` (`:zstd`) of 40 M sites | 0.27 s | — |
| `read_methylation`, zstd (170 MB on disk) | **0.27 s** | — |
| `read_methylation`, uncompressed (322 MB, mmap) | **0.002 s** | — |
| 12k × `find_calls_in_range` sweeps, in-memory `Vector` | 0.8 ms | — |
| same sweeps on mmap'd Arrow columns | 1.6 ms | — |

Two structural facts fall out immediately:

1. The plain-text parse is excellent — ~2.7 GB/s at 8 threads, 4.1× scaling —
   and is **not** the bottleneck.
2. The `.gz` path scales only 1.8× because CodecZlib's serial inflate (~2 s per
   file, vs ~1 s for Apple's accelerated zlib) is the floor. Everything else
   waits on it.

---

## 1. Data structures: SoA vs AoS — verdict: already right

`AggregatedCall` as 8 packed bytes in a `StructArray` is the correct layout and
should not change:

- `calls.pos` is a contiguous `Vector{UInt32}`, so `find_calls_in_range`'s
  binary search touches only the position column — an AoS layout would drag the
  payloads through cache for every probe.
- 8 B/site means a whole-genome, all-sample zebrafish dataset is ~320 MB — it
  fits in 16 GB with room for the genome and the sparse frequency structures.
- The level-plus-depth payload (vs two counts) is a space decision already made
  correctly; nothing to revisit.

The one structure worth shrinking is **`Exploration.FeatureLevels`**: two
parallel `SparseVector{_, Int}`s over the *same* sparsity pattern cost
`(8+4) + (8+4) = 24 B` per stored base, and the two `Int64` index vectors are
duplicates of each other.

- Cheap fix: `Int32` indices (`sparsevec` with `Vector{Int32}` indices — region
  lengths are far below `typemax(Int32)`) → 16 B/entry, −33 %.
- Full fix: a small custom struct holding **one** shared `Vector{Int32}` of
  indices plus `Vector{Float32}` levels and `Vector{UInt32}` weights →
  12 B/entry, −50 %. `gene_profile(::FeatureLevels)` already walks the two
  stored-entry lists in step, so it would actually get *simpler* (one list, no
  dual-cursor sync), and the "weight ≠ 0 marks measured" convention is easier to
  enforce with one pattern. For ~25 k genes × hundreds of covered bases each,
  this is a few hundred MB saved per `MethylationFrequency` — and the notebooks
  hold one per tissue loop iteration plus a pooled one.

---

## 2. The single biggest win: stop re-parsing `.cov.gz` — use your own Arrow layer

`write_methylation`/`read_methylation` exist, are well designed (mmap-able SoA
columns, per-scaffold files, name manifest) — and **no notebook uses them**.
Meanwhile a render of `methylation_vs_ds.qmd` parses 4.8 GB of gzip twice
(~21 files × 2.4 s × 2 ≈ 100 s of pure I/O per render, re-paid on every
iteration of a figure).

Recommended workflow (bit-identical results):

1. One-off conversion script: for each sample,
   `write_methylation(dir, load_bismark_cov(path))`. Cost ≈ one current render's
   I/O. ~170 MB/sample zstd or ~320 MB uncompressed (the original `.gz` is
   232 MB, so uncompressed Arrow is only ~1.4× the gzip on disk).
2. Notebooks then build per-tissue and pooled datasets as
   `merge_calls(read_methylation(dir) for dir in sample_dirs)` **in the same
   `paths` order** as today. `load_bismark_cov(paths)` itself merges per file in
   `paths` order, so folding the same per-file datasets in the same order is
   exactly the same arithmetic — bit-identical output, no analysis change.
3. Per-file reload cost drops from 2.4 s to 0.27 s (zstd) or ~0 s + lazy page-in
   (uncompressed). The double parse disappears as a concern entirely.

Compression choice: keep `:zstd` for archival copies; prefer **uncompressed**
for the working set you re-render against — it memory-maps (0.002 s "load",
pages evicted under pressure instead of living on the Julia heap, a real benefit
on 16 GB), and the ~2× slower binary-search sweep on mapped columns (1.6 ms vs
0.8 ms for 12 k queries) is noise.

Do **not** replace the pooled `load_bismark_cov(all_files)` with
`merge_calls(per_tissue_datasets)`: that changes the fold order, and
`merge_calls` is not perfectly associative once quantization bites (documented
in the module). The Arrow route above keeps the fold order and is exact.

## 3. gzip is the serial bottleneck of the current path

If you keep loading from `.cov.gz` at all: CodecZlib inflates at ~0.5 GB/s on
one core while the parser can eat 2.7 GB/s. Options, in order of preference:

- Convert to Arrow (above) and the problem vanishes.
- Store the working `.cov`s uncompressed — the plain path scales 4× and runs at
  0.5 s/GB (disk is not the limit on Apple SSDs).
- If gzip must stay, `max_concurrent_files` (default 4) is what hides the serial
  inflate; on an 8-core M3 with 21 files the current default is about right —
  raising it would help wall time slightly but each in-flight file holds ~2 GB
  of transient allocations (measured), which 16 GB cannot afford at 6+.

## 4. Parallelize the merge across scaffolds (exact, easy)

`merge_calls` folds datasets serially *and* merges scaffolds serially within
each fold. Scaffolds are independent — spawn `merge_scaffold` per scaffold
(exactly like `load_cov_stream` already spawns `collapse_cov!`):

```julia
# inside merge_calls, per dataset: collect the scaffold pairs first, then
@sync for (name, calls) in pairs_to_merge
    Threads.@spawn merged[name] = merge_scaffold(existing[name], calls)
end
```

Same pairwise arithmetic, same fold order, bit-identical output; a 21-file merge
goes from ~20 × 0.23 s serial to ~20 × 0.23/min(8, n_scaffolds) s. It also
overlaps the merge with the gated file loads instead of making fetched results
queue behind a single-threaded fold — which matters on 16 GB, because completed
`Task`s each pin a ~320 MB dataset until the fold consumes them.

A k-way (heap) merge that sums all N counts per site before packing once would
cut the 20 × 0.5 GB of intermediate-array churn as well — but it rounds once
instead of N−1 times, i.e. output differs at the quantization level. Worth
having only as an opt-in flag; the per-scaffold spawn is the safe default.

## 5. Thread `Exploration.feature_frequency` (exact, high leverage)

The methylation `feature_frequency` is a serial double loop over scaffolds ×
features, and the notebooks call it ~12× per render (once per tissue + pooled).
Each feature's work (binary search + run-pool + `sparsevec`) is independent.
Parallelize per scaffold into per-task `Dict`s and merge (feature IDs are unique
across scaffolds, so the merge is a plain `merge!`):

```julia
results = [Threads.@spawn _scaffold_levels(genome, tree, calls, ...) for ...]
features = reduce(merge!, fetch.(results); init = Dict{String,FeatureLevels}())
```

Same per-feature arithmetic, order-independent by construction, bit-identical.
`calculate_frequency` already sets the precedent with `Threads.@threads` over
scaffolds.

Minor, same function: `_region_levels` could `sizehint!` its three vectors to
`length(region_calls)`, and `feature_frequency` could pass the indices/values
straight into the shared-index `FeatureLevels` of §1 instead of building two
`SparseVector`s (each `sparsevec` call re-sorts and copies).

## 6. Julia-specific points in the module itself

The `.cov` hot path is already close to optimal — byte-level fields, no
`SubString`/`tryparse`, `@inbounds` where it counts, `findnext(==(0x0a), ...)`
hits the memchr fast path, pooled 4 MiB chunks, `CovChunk` as byte ranges.
Remaining items, none huge:

- **`load_bismark` (extractor path) is the laggard**: `eachline` allocates a
  `String` per line and `parse_bismark_line` works through `SubString`s. Fine
  for the micro fixtures; if multi-GB extractor files ever return, port it to
  the `.cov` chunk machinery. Also `aggregate_keys!` should
  `sizehint!(positions/payloads, n_unique_estimate)` and its `sort!(keys)` is
  the textbook case for `SortingAlgorithms.RadixSort` (UInt64 keys, ~5× on
  large arrays). Not worth touching until the path is actually used at scale.
- `parse_cov_parallel`'s ~1.4 k lock conflicts per gz file (measured) come from
  the bounded `Channel` — harmless at current chunk sizes; if you ever shrink
  `COV_CHUNK_BYTES`, batch chunks instead.
- `merge_scaffold` re-reads `left[left_index]`/`right[right_index]` as whole
  structs; indexing the `pos`/`payload` columns directly avoids constructing an
  `AggregatedCall` per comparison. Micro — the merge is already 0.23 s per
  genome pass.
- Keep returning `Union{Nothing, T}` from the parsers — small-union dispatch is
  free here; no need for sentinel values.

## 7. Arrow usage — effective; two refinements

The design (per-scaffold files + manifest, columns wrapped without copying,
zero-copy mmap when uncompressed) is exactly how Arrow.jl is meant to be used.
Refinements:

- Consider `Arrow.write(...; ntasks)` default is fine, but for the conversion
  script write samples **concurrently** (one task per sample) — writes are
  0.27 s each, the script will be dominated by the initial `.gz` parses anyway.
- `read_methylation_arrow` wraps `Arrow.Primitive` columns directly in the
  `StructArray`; per-element access is marginally slower than `Vector`. For a
  dataset you're about to hammer with millions of point queries, `copy(table.pos)`
  /`copy(table.payload)` materialises to plain `Vector`s at ~memcpy cost — worth
  a `materialize::Bool = false` keyword rather than a behavioural change.

## 8. GPU (Metal.jl) — honest verdict: not yet

Measured hot kernels are already sub-second per whole genome on the CPU
(parse 0.48 s, merge 0.23 s, collapse ~no-op on sorted input), and the true
bottlenecks are serial gzip inflate and redundant re-parsing — neither is GPU
work. Specific to this codebase and an M3:

- The byte-level parser is branchy, variable-length work — a poor GPU fit, and
  at 2.7 GB/s it already approaches the SoC's shared memory bandwidth, which is
  the same resource the GPU would compete for.
- The natural GPU shapes here are segmented reductions (`collapse_cov!`,
  `merge_scaffold`) and scatter-adds (`_group_profiles` / `quantile_heatmap`
  accumulation). Unified memory + Metal.jl's shared buffers would make the
  transfer cost near-zero — but the CPU versions are 0.2–0.5 s, so the ceiling
  on total savings is a few seconds per render *after* the I/O fixes above.
- If profile accumulation ever grows 10× (many more groups × positions ×
  samples), that scatter-add is the first candidate — note Metal has no
  `Float64`, so the `Float64` accumulators in `gene_profile`/`_group_profiles`
  would need `Float32` accumulation, which *is* an analysis-visible change and
  should be gated the same way as the k-way merge.

Priority order on this machine: threads everywhere (free), Arrow caching
(§2), threaded merge + `feature_frequency` (§4–5), and only then think about
Metal.

## 9. Notebook-level items (no analysis changes)

- **`primates/notebooks/exploration/` has no `_metadata.yml`** — unlike the
  fish repo, whose `_metadata.yml` passes `--threads=auto`. Unless
  `JULIA_NUM_THREADS` is set globally, primates renders run **single-threaded**:
  the measured difference is 2–4.5× on every load. Add the same `_metadata.yml`
  (minus the sysimage line, or build one). Likewise set the VS Code Julia
  extension's `julia.NumThreads` to `"auto"` so interactive REPL work gets the
  threaded paths.
- `expression_vs_methylation.qmd` keeps `methylation_cov_z` (the pooled ~320 MB
  dataset) alive for the entire notebook although nothing after
  `feature_frequency` uses it — `methylation_cov_z = nothing; GC.gc()` right
  after, as the per-tissue loops already do. On 16 GB this is the difference
  between staying resident and paging during the per-tissue loop.
- Both `methylation_vs_ds.qmd` and `expression_vs_methylation.qmd` traverse all
  21 files twice (tissue groups ∪ pooled). With the §2 Arrow cache this becomes
  cheap; until then, it is the single largest time cost of a render.
- The per-tissue `tissue_data = nothing; tissue_freq = nothing; GC.gc()`
  pattern is correct and worth keeping — it is what bounds peak memory today.

## 10. What is already right (don't churn)

- 8-byte packed SoA calls; view-returning binary-searched range queries.
- Chunked producer/consumer with a recycling `ChunkPool`; byte-range chunks;
  allocation-free per-record parsing with the two-layer scaffold-name cache.
- Chunk-order concatenation and `paths`-order merging (scheduling-independent
  output) — preserve these invariants in any parallelization from §4–5.
- Per-scaffold spawned `collapse_cov!`; semaphore-gated concurrent files.
- Arrow persistence layer design (it just needs to actually be used).

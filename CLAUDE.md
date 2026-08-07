# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

BioinfoTools2.jl is a pure-Julia bioinformatics package for genome annotations (GFF3), interval data (BED), and per-sample tabular measurements, plus exploration/statistics on top of them.

## Setup

Depends on **forked** packages that are not on the public registry. Add them before building or testing:

```julia
using Pkg
Pkg.add([
    PackageSpec(url="https://github.com/tdw-89/Indexes.jl"),
    PackageSpec(url="https://github.com/tdw-89/GenomicFeatures.jl"),
    PackageSpec(url="https://github.com/tdw-89/GFF3.jl.git"),
    PackageSpec(url="https://github.com/tdw-89/BED.jl.git"),
])
```

Exact CI setup: `.github/workflows/main.yml` (Julia 1.12, ubuntu-latest).

## Common commands

Run from the repo root against the package's own environment:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'   # full test suite
./format.sh                                    # jlfmt, config baked into the script
```

`test/TESTS_*.jl` are `include`d by `test/runtests.jl` and assume `using BioinfoTools2` and `using Test` are already active — they are not standalone scripts. To run one in isolation: `julia --project=test`, then `using BioinfoTools2, Test; include("test/TESTS_reference.jl")`.

Fixtures live in `test/data/`. `micro_CpG_OT.txt` is a 406-record excerpt of a Bismark CpG_OT file (chromosome 15, a 5 kb window over 76 sites); `micro.cov` is the coverage counterpart over the same window (151 cytosines, i.e. both strands of those CpGs). The overlap is deliberate: the 76 forward-strand sites in `micro.cov` carry counts identical to the aggregate of `micro_CpG_OT.txt`, which is what `TESTS_methylation.jl` cross-checks the two loaders against. The multi-GB originals are gitignored (`*.deduplicated.txt`, `*.bismark.cov`); tests must use only the micro fixtures.

`jlfmt` via `format.sh` is the only linter.

### Testing conventions

Never put `@test` inside a `for` loop — it inflates the reported test count without adding coverage and hides which iteration failed. Collect each iteration's boolean into a `Vector{Bool}`, then assert once after the loop with `@test all(...)` (or `@test any(...)` when one match suffices). Aim for maximum meaningful coverage per assertion, not a high test count.

### Style

- **Docstrings, not comments, on definitions.** Every function, struct and meaningful constant gets a `"""..."""` docstring, exported or not. Reserve `#` for explaining a step *inside* a body.
- **Concise and minimal-but-informative.** State what is required, what happens to the input, and any pitfall. Name code elements rather than paraphrasing them. Do not justify or advertise a design decision — a directive ("don't replace X with a `Dict`") belongs; a paragraph defending X does not.
- **Descriptive snake_case names.** Prefer `meth_count`, `line_start`, `n_records` over `m`, `lo`, `n`. This covers loop and index variables too (`for offset in ...`, not `for k in ...`); one- and two-letter abbreviations are not acceptable just because a scope is short.

## Architecture

`src/BioinfoTools2.jl` is a thin top-level module that `include`s and re-exports the submodules. Read them in this order; each builds on the last.

0. **`BitCodes`** (`src/bit_codes.jl`) — the bit-packing vocabulary every other module encodes against. Holds the 2-bit strand codes (`STRAND_FWD` = 0, `STRAND_REV` = 1, `STRAND_BOTH` = 2, `STRAND_NA` = 3), conversions to/from `GenomicFeatures.Strand` and strand characters (`strand_code`/`decode_strand`/`strand_char`/`get_strand`), the generic field helpers (`get_field`/`set_field`/`field_mask`/`clamp_field`), and the pack/parse helpers for `Reference`'s 64-bit feature metadata code (`pack_metadata`/`parse_index`/`parse_strand`/`parse_so_term`, re-exported as aliases by `Reference`). Every layout uses the 2-bit strand encoding — the narrowest slot, in the methylation payload — and wider slots leave their upper bits zero. **`0x00` means forward, not unknown**: a zeroed code decodes to `+`, so always write a strand explicitly, `STRAND_NA` when there isn't one. Add new packed fields here rather than hand-rolling shifts at the call site.

1. **`SOTerms`** (`src/so_terms.jl`) — internal, not exported at the top level. Loads `assets/SOFA.json` at build time into a `SOTermLookup`, a tuple-based binary-search table mapping SO term short codes (`UInt16`), full ontology IDs, and Symbol labels (`:gene`). Refresh the snapshot with `./update_assets.sh`.

2. **`Reference`** (`src/reference.jl`) — the in-memory genome: `Species` → `Genome` → `Scaffold` → interval tree of features. Each feature interval carries a packed 64-bit `code` (layout diagram in the `Genome` docstring) holding an SO term, a strand, and an index into a string-interning metadata store (`vocab`/`meta_offsets`/`meta_blob`). Don't replace that design with a `Dict`-of-structs; it is what makes whole-genome GFF3 loads fit in memory. `add_features!` parses GFF3 (optionally gzipped) producer/consumer style: the main task batches records onto a `Channel`, a spawned task commits them. `Genome`/`Scaffold` `getindex` takes a single ID, a `Vector` of IDs, or a metadata index, and is an O(features) linear scan by design — no ID index is maintained (rationale in the docstrings). Lookups return a `FeatureRecord`, carrying the raw `code` alongside ID/type/coordinates/metadata so callers needing the packed fields don't re-resolve them. `IntervalSimple` is a bare `(start_pos, end_pos, code)` struct for where no scaffold/vocab is needed (e.g. `Paralogs.GeneFamily`'s per-scaffold `StructArray`s).

3. **`Data`** (`src/data.jl`) — sample-level data loaded *against* a `Genome`. `BedData` wraps per-scaffold interval trees from a BED file. `TabularData{T}` wraps a numeric matrix whose rows match genome features by ID. `Experiment{T}` is a trie of `Variable`s built from a sample sheet (a `DataFrame` whose last column is file paths and preceding columns categorical variables); samples load concurrently via `Threads.@spawn` and fold into an `ImmutableDict` tree in `_build_variables`. Also provides `intersect` (against a `Genome`/`Scaffold`/`BedData`), `leftjoin` (two interval trees, on metadata/start/end/interval), and `merge_segments`.

3a. **`Data.Methylation`** (`src/data/methylation.jl`) — single-base methylation calls: a sample-level signal, hence its place under `Data`, but a purely *positional* one. It touches no `Genome`, `Scaffold` or interval type and depends on nothing but `BitCodes`; keep it that way, so a whole-genome call set stays one flat, memory-mappable block. `Data` re-exports the submodule, so `BioinfoTools2.Methylation` reaches it alongside `BioinfoTools2.Data.Methylation`.

   - `AggregatedCall` is exactly 8 bytes: `pos::UInt32` plus a packed `payload::UInt32` holding a 16-bit depth (saturating at `MAX_COUNT` = 65535), an 8-bit methylation level, a 2-bit context (`CTX_CPG`/`CTX_CHG`/`CTX_CHH`/`CTX_UNKNOWN`), a 2-bit strand, and 4 reserved bits.
   - **The payload stores a level and a depth, not two counts** — two 12-bit counts would cap each at 4095, losing exactly the deepest sites. The level maps `[0, 100]` percent linearly onto `[0, 255]` (`encode_percent`/`decode_percent`, `MAX_PERCENT_CODE` = 255). `get_depth`/`meth_percent`/`meth_fraction` read stored fields; `get_meth`/`get_unmeth` **reconstruct** counts as `round(depth * percent)` — exact for depth ≤ 255 (verified exhaustively in the tests), off by a read or two beyond, with `get_meth + get_unmeth == get_depth` always holding. A saturating depth never corrupts the level, which is encoded before the clamp. Build payloads with `pack_payload(meth, unmeth, …)` or `pack_percent_payload(percent, depth, …)`.
   - `MethylationData` holds calls per scaffold as `StructArray{AggregatedCall}`, sorted by position, so `find_calls_in_range` binary searches the contiguous `pos` column and returns a **view**, not a copy.
   - `load_bismark` reads methylation-extractor output (optionally gzipped), discarding read IDs and summing per (position, context, strand) via a sort-and-scan over packed `UInt64` keys. The **strand comes from the file name** (`OT`/`CTOT` forward, `OB`/`CTOB` reverse — `infer_strand`); Bismark's second column is the methylation state, not a strand.
   - `load_bismark_cov` reads the *already aggregated* `.cov` output of `bismark2bedGraph`/`coverage2cytosine` (one 1-based line per cytosine with methylated/unmethylated counts), so it only buckets by scaffold and collapses. `collapse_cov!` sorts and sums repeated positions; both are no-ops for a normal `.cov`, which is already sorted and unique. A `.cov` records **neither context nor strand**: `context` defaults to `CTX_CPG` (Bismark's default coverage output is CpG-only, and a mixed `--CX` file cannot be split by context afterwards) and `strand` to `infer_strand(path)`, i.e. `STRAND_NA` for the usual strand-pooled whole-sample file.
   - `.cov` parsing is parallel at two levels, `.cov` being the multi-GB input in practice.
     - *Within a file* (`parse_cov_parallel`): one task reads while `n_workers` tasks parse off a bounded channel into their own `CovBuffer`s. `chunk_lines!` emits a `CovChunk` — a byte **range** into a `ChunkPool` buffer, not a trimmed copy — so a chunk costs no allocation, and the worker recycles the buffer when done. A pooled buffer holds stale bytes past `stop_index`; `parse_cov_chunk` must stop there, never at the next newline it happens to find. Only a line straddling two reads is copied, into a short chunk of its own that the pool declines to recycle.
     - *Across files* (`load_bismark_cov(paths)`): up to `max_concurrent_files` load at once behind a `Base.Semaphore`, splitting `nthreads()` between them via `load_cov_stream`'s `n_workers`. This exists for gzip, whose decompression is serial per stream and so caps what one file can reach.
     - Results concatenate **in chunk order**, and files merge in `paths` order — don't drop either; that is what keeps a load scheduling-independent and leaves an already-sorted `.cov` sorted.
     - Per-line work is byte-level throughout: `parse_cov_record`/`parse_uint32` read fields straight out of the buffer, never via `SubString` or `tryparse`, and a two-layer scaffold-name cache (last name seen, then a bounded scan, both byte-comparing via `name_matches`) avoids a `String` per record. `parse_uint32` takes digits only — no sign, no surrounding whitespace — so malformed-field tests need non-digit garbage rather than `"+7"`.
   - Both loaders produce the same `MethylationData` and combine via `merge_calls`.
   - `write_methylation_arrow`/`read_methylation_arrow` persist one scaffold as Arrow (Zstd by default; compression and zero-copy `mmap` are a tradeoff — uncompressed files map, compressed ones decompress into memory). `write_methylation`/`read_methylation` do a whole dataset as one file per scaffold plus a `scaffolds.tsv` manifest, since scaffold names are not always valid file names.

4. **`Paralogs`** (`src/paralogs.jl`) — depends on `Reference` (`Genome`/`IntervalSimple`) plus `Graphs`, `SimpleWeightedGraphs`, `SparseArrays`, `StructArrays`. Centres on `GeneFamily`, built from a `Genome` and a paralog-pair `DataFrame` (2–6 columns: mandatory query/subject ID, then optional `dN`/`dS`/`id_subject_query`/`id_query_subject`, matched *by column name*, case-insensitively). Genes are resolved against the genome, grouped by scaffold, and given a contiguous linear index (`scaffold_ranges`/`intervals`/`id_to_index`). Every relation is a `SparseMatrixCSC` over that index, oriented so the axis you look up by is column-major: `topology` is a symmetric `Bool` adjacency; `dN`/`dS` are symmetric and stored **once** in the upper triangle (`[min(a,b), max(a,b)]`); `id_subject_query`/`id_query_subject` are directed (subjects and queries on the columns respectively). Rows with a `NaN` in any matched value column are dropped up front with a warning, since `NaN != NaN` breaks the `*_graph` symmetry checks. There is exactly one *public* constructor (an unexported `_RawFields` marker gates the raw-fields one used by sub-family `getindex`), plus `getindex` (`gf[i]`/`gf["id"]` → a column dict; `gf[indices]`/`gf[ids]` → a sub-family) and `Base.show`. The `*_graph` functions convert relations to graph views; dN/dS are symmetrised as `m + permutedims(m)` because `SimpleWeightedGraph` requires full symmetry, while the directed pair pass straight to `SimpleWeightedDiGraph`. `rbh` has two methods: `rbh(df::DataFrame; scoring)` scores a flat pair table by percent identity (max/mean/double_max, plus `rbh_ds` for dS) with no genome context; `rbh(gf::GeneFamily; scoring)` finds `gf.topology`'s connected components and keeps each one's single lowest-distance edge (ascending `dS`, then `dN`, then `%ID`, skipping levels `gf` lacks), warning about and deterministically resolving ties.

5. **`Plotting`** (`src/plotting.jl`) — stub (`using CairoMakie`, no functions).

6. **`Exploration`** (`src/exploration.jl`) — statistics over `Reference`/`Data`: per-feature `coverage` fractions from `BedData` or mean methylation levels from `MethylationData`; `kde` density estimates (on either, or on `TabularData`); three `quantiles` overloads — value-based on `TabularData` (collapses each row via `merge`, default `mean`, then bins by the *value* quantiles of those scalars), rank-based on `TabularData` (`quantiles(data, ranking::Vector{String})`, cascading tie-breaks over named variables then original sample order), and rank-based on a `DataFrame` (`quantiles(pairs, ranking)`, same scheme over a paralog-pair table shaped like `GeneFamily`'s input or `rbh`'s output, appending a `"quantile"` column); and a genome-wide coverage pipeline: `calculate_frequency` (per-base overlap counts across many `BedData`, `Threads.@threads` over scaffolds using a difference-array sweep) → `feature_frequency` (projects it onto each feature ± a flank, strand-oriented so index 1 is the 5' end) → `gene_profile`/`mean_gene_profile` (interpolates each body to a fixed bin count and averages into a metagene profile).

   The same three functions have a **methylation** path that deliberately does *not* share the `BedData` one's semantics, and the difference is the whole point — don't "unify" them:
   - `feature_frequency(genome, feature, ::MethylationData)` skips `calculate_frequency` (a `MethylationData` is already per-base) and returns a `MethylationFrequency` of `FeatureLevels`, not a `FeatureFrequency`. Each feature stores two sparse vectors: `levels` (depth-weighted fraction in `[0, 1]`) and `weights` (total depth). **A nonzero weight, not a nonzero level, is what marks a base as measured** — a base measured at 0% and a base with no cytosine are both structural zeros in `levels`. Calls are filtered by `min_depth` (default `DEFAULT_MIN_DEPTH` = 5) and `context` (default `CTX_CPG`, `nothing` for all); several calls on one base combine in proportion to their depths. The region is *not* clipped at the scaffold start, so index `flank + 1` is always the feature's first base.
   - `gene_profile(::FeatureLevels)` returns a `(levels, weights)` NamedTuple, and **bins** the body rather than interpolating it — interpolating would invent levels across the long uncovered stretches between cytosines.
   - `mean_gene_profile(::MethylationFrequency)` averages each position only over the features actually measured there, so the curve is mean level among measured cytosines rather than level × cytosine density; positions nothing measured come back as `NaN`, not `0`.
   - Weighting: `weight_by_depth` (default `true`) and `weight_transform` (default `identity`; `log`/`sqrt` compress deep coverage) are read once per base by `weight_of`, and everything above that aggregates linearly — don't re-apply the transform to an already-summed weight. Pooling several calls on *one* cytosine stays linear in depth regardless, since that step reconstructs the site's own methylation fraction. `docs/weighting.md` documents all three steps; keep it in step with the code.
   - Everything is on the `[0, 1]` fraction scale, `coverage` included (sum of per-base levels over the feature *length*, so uncovered bases drag it down).

7. **`Modeling`** (`src/modeling.jl`) — stub (imports GLM/HypothesisTests/MultipleTesting/StatsBase, no functions).

### Cross-cutting conventions

- Genomic intervals are `IntervalTree{UInt32, IntervalValue{UInt32,UInt64}}` (aliased `IntervalTreeM64`), 1-based closed `[start, end]`. BED's 0-based half-open coordinates are converted on load (`load_bed`).
- An interval's 64-bit `value` is always a packed metadata code. Unpack with `BitCodes.parse_index`/`parse_strand`/`parse_so_term` (re-exported by `Reference`; gene features) or `Data.parse_bed_strand` (BED features), not inline bit-twiddling. For methylation payloads use `Methylation.get_depth`/`meth_percent`/`get_meth`/`get_unmeth`/`get_context`/`is_forward`, preferring `get_depth`/`meth_percent` — those are stored, the counts are reconstructions.
- Strand is encoded identically everywhere via `BitCodes`: gene and BED codes keep an 8-bit slot, methylation payloads a 2-bit one, all storing the same values. Encode with `strand_code`, decode with `decode_strand`; don't reintroduce a local strand mapping.
- Feature/sample IDs resolve lazily against the interned vocab (`Reference.get_metadata_id`) instead of a separate index. This memory/lookup-speed tradeoff repeats across `Reference` and `Data`; preserve it rather than adding ID `Dict`s.
- Parallel hot paths: `add_features!`'s channel producer/consumer, `Experiment`'s per-sample `Threads.@spawn`, `calculate_frequency`'s `Threads.@threads` over scaffolds, and `load_bismark_cov`'s per-file semaphore, chunked producer/consumer, and per-scaffold `collapse_cov!`. Run `julia --project=. -t auto` to exercise them. Each is order-independent by construction and the tests assert it; don't let output depend on the thread count.
- Most fixtures are smaller than one `.cov` chunk, so anything touching `chunk_lines!` or `ChunkPool` needs a test exceeding `COV_CHUNK_BYTES` (`TESTS_methylation.jl`'s "multi-chunk input") and one that reuses a pooled buffer ("chunk buffer pooling").
- `quantiles` lives entirely in `Exploration`, including the `DataFrame` overload used for paralog-pair tables. Two sibling modules each declaring and exporting `quantiles` produces an ambiguous binding as soon as both are loaded in one session (as `runtests.jl` does), and it breaks at call time, not compile time. Don't reintroduce a same-named function in another module without either extending the existing generic or picking a different name.

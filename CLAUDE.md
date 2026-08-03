# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

BioinfoTools2.jl is a pure-Julia bioinformatics package for working with genome annotations (GFF3), interval data (BED), and per-sample tabular measurements, plus downstream exploration/statistics on top of them.

## Setup

This package depends on **forked** versions of several packages that are not on the public registry. They must be added from the maintainer's GitHub forks before the package will build or run tests:

```julia
using Pkg
Pkg.add([
    PackageSpec(url="https://github.com/tdw-89/Indexes.jl"),
    PackageSpec(url="https://github.com/tdw-89/GenomicFeatures.jl"),
    PackageSpec(url="https://github.com/tdw-89/GFF3.jl.git"),
    PackageSpec(url="https://github.com/tdw-89/BED.jl.git"),
])
```

See `.github/workflows/main.yml` for the exact CI setup (Julia 1.12, ubuntu-latest).

## Common commands

Run from the repo root, using the package's own environment (`--project=.`):

```sh
# Run the full test suite
julia --project=. -e 'using Pkg; Pkg.test()'

# Or from within the Julia REPL
julia --project=.
julia> using Pkg; Pkg.test()

# Run a single test file directly (faster iteration)
julia --project=. test/TESTS_exploration.jl   # note: these are includes, not standalone — see below

# Format code (uses jlfmt, config baked into the invocation)
./format.sh
```

Each `test/TESTS_*.jl` file is `include`d from `test/runtests.jl` and assumes `using BioinfoTools2` and `using Test` are already active — they are not independently runnable scripts. To run just one test file in isolation, start a REPL with `julia --project=test`, `using BioinfoTools2, Test`, then `include("test/TESTS_reference.jl")`.

Test fixture files (small GFF3/BED/narrowPeak/Bismark samples) live under `test/data/`. `micro_CpG_OT.txt` is a real 406-record excerpt of a Bismark CpG_OT file (chromosome 15, a 5 kb window covering 76 sites). The full ~4.5 GB dataset it was cut from (`CpG_OT_SRR12173567_1_val_1_bismark_bt2_pe.deduplicated.txt`) is far too large to commit and is gitignored via `*.deduplicated.txt` — tests must only use the micro fixture.

`micro.cov` is the Bismark **coverage** counterpart, cut from the same sample's `SRR12173567_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov` (~1.1 GB, gitignored via `*.bismark.cov`) over exactly the same chromosome-15 window: 151 cytosines, i.e. both strands of the same CpGs. The overlap is deliberate — the 76 forward-strand sites in it carry counts identical to the aggregate of `micro_CpG_OT.txt`, which is what `TESTS_methylation.jl` cross-checks the two loaders against.

There is no linter configured beyond `jlfmt` via `format.sh`.

## Architecture

The package is a thin top-level module (`src/BioinfoTools2.jl`) that just `include`s and re-exports its submodules. Read them in this order — each one builds on the last:

0. **`BitCodes`** (`src/bit_codes.jl`) — the shared bit-packing vocabulary every other module encodes against. Holds the canonical 2-bit strand codes (`STRAND_FWD` = 0, `STRAND_REV` = 1, `STRAND_BOTH` = 2, `STRAND_NA` = 3), conversions to/from `GenomicFeatures.Strand` and strand characters (`strand_code`/`decode_strand`/`strand_char`/`get_strand`), and the generic field helpers (`get_field`/`set_field`/`field_mask`/`clamp_field`) used to read and write packed codes. Two bits is the narrowest strand field any layout has (the methylation payload), so every layout uses that encoding; wider slots just leave their upper bits zero. **Note `0x00` means forward, not unknown** — a zeroed code decodes to `+`, so always write a strand explicitly (`STRAND_NA` when there isn't one). Add new packed fields here rather than hand-rolling shifts at the call site.

1. **`SOTerms`** (`src/so_terms.jl`) — internal only, not exported at the top level. Loads `assets/SOFA.json` (Sequence Ontology Feature Annotation) at package build time and builds a `SOTermLookup`, a tuple-based binary-search table mapping between SO term short codes (`UInt16`), full ontology IDs, and Symbol labels (e.g. `:gene`). Refresh the ontology snapshot with `./update_assets.sh`.

2. **`Reference`** (`src/reference.jl`) — the core in-memory genome representation. Key types: `Species` → `Genome` → `Scaffold` → interval tree of features. Each feature interval carries a packed 64-bit `code` (see the bit-layout diagram in the `Genome` docstring) encoding an SO term + strand + an index into a separate string-interning metadata store (`vocab`/`meta_offsets`/`meta_blob`). This packed-code + interned-string design is deliberate for memory efficiency when loading whole-genome GFF3 annotations — don't casually replace it with a `Dict`-of-structs without understanding why. `add_features!` parses a GFF3 file (optionally gzipped) using a producer/consumer pattern: the main task parses records into batches on a `Channel`, a spawned task commits them into the `Genome`. `Genome`/`Scaffold` support ID-based lookup via `Base.getindex`, which is an O(features) linear scan by design (no ID index is maintained) — see the docstrings for the rationale.

3. **`Data`** (`src/data.jl`) — sample-level data loaded *against* a `Genome`. `BedData` wraps per-scaffold interval trees from a BED file. `TabularData{T}` wraps a numeric matrix whose rows are matched to genome features by ID. `Experiment{T}` is a trie of `Variable`s built from a sample sheet (`DataFrame` whose last column is file paths, preceding columns are categorical variables); samples are loaded concurrently via `Threads.@spawn` and folded into an `ImmutableDict`-based tree by `_build_variables`. Also provides `intersect` (interval-tree intersection against a `Genome`/`Scaffold`/`BedData`), `leftjoin` (join two interval trees on metadata/start/end/interval), and `merge_segments`.

3a. **`Data.Methylation`** (`src/data/methylation.jl`) — single-base methylation calls: another sample-level signal, hence its home under `Data`, but a purely *positional* one. It touches no `Genome`, `Scaffold` or interval type, and depends on nothing but `BitCodes` — that independence is what lets a whole-genome call set stay one flat, memory-mappable block, so keep it. `Data` re-exports the submodule, so `BioinfoTools2.Methylation` and `using BioinfoTools2.Methylation` still reach it alongside the fully-qualified `BioinfoTools2.Data.Methylation`. `AggregatedCall` is exactly 8 bytes — `pos::UInt32` plus a bit-packed `payload::UInt32` holding a 16-bit total read depth (saturating at `MAX_COUNT` = 65535), an 8-bit methylation level, a 2-bit context (`CTX_CPG`/`CTX_CHG`/`CTX_CHH`/`CTX_UNKNOWN`), a 2-bit strand, and 4 reserved bits. **The payload stores a level and a depth, not two counts** — two 12-bit counts would cap each at 4095 and destroy the depth of exactly the deep sites whose statistics are most trustworthy, so the depth is kept whole and the spare byte holds the level, linearly mapped from `[0, 100]` percent onto `[0, 255]` (`encode_percent`/`decode_percent`, `MAX_PERCENT_CODE` = 255). `get_depth`/`meth_percent`/`meth_fraction` read what is stored; `get_meth`/`get_unmeth` **reconstruct** counts as `round(depth * percent)` and are exact for any depth ≤ 255 (verified exhaustively in the tests), drifting by at most a read or two beyond that while `get_meth + get_unmeth == get_depth` always holds. Saturating the depth never corrupts the level, since the level is encoded from the counts before the depth is clamped. Build payloads with `pack_payload(meth, unmeth, …)` from counts or `pack_percent_payload(percent, depth, …)` from a level. Calls are held per scaffold in `MethylationData` as `StructArray{AggregatedCall}` (struct-of-arrays), sorted by position so `find_calls_in_range` is a binary search over the contiguous `pos` column and returns a **view**, not a copy. `load_bismark` parses Bismark methylation-extractor output (optionally gzipped), discarding read IDs and summing calls per (position, context, strand) via a sort-and-scan over packed `UInt64` keys; the **strand comes from the file name** (`OT`/`CTOT` = forward, `OB`/`CTOB` = reverse, see `infer_strand`) because Bismark's second column is the methylation state, not a strand. `load_bismark_cov` reads the other, *already aggregated* Bismark output — `.cov` coverage files from `bismark2bedGraph`/`coverage2cytosine`, one 1-based line per cytosine with methylated/unmethylated counts — so it only has to bucket by scaffold and collapse (`collapse_cov!` sorts and sums repeated positions; a normal `.cov` is already sorted and unique, so both are no-ops). A `.cov` records **neither context nor strand**, so both are caller-supplied kwargs: `context` defaults to `CTX_CPG` (Bismark's default coverage output is CpG-only; a mixed `--CX` file cannot be split by context after the fact) and `strand` to `infer_strand(path)`, which is `STRAND_NA` for the usual strand-pooled whole-sample file. Both loaders produce the same `MethylationData`, so their results can be combined with `merge_calls`. `write_methylation_arrow`/`read_methylation_arrow` persist one scaffold as Arrow (Zstd by default; note compression and zero-copy `mmap` are a tradeoff — uncompressed files map, compressed ones decompress into memory), and `write_methylation`/`read_methylation` do a whole dataset as one file per scaffold plus a `scaffolds.tsv` manifest (scaffold names are not always valid file names).

4. **`Paralogs`** (`src/paralogs.jl`) — standalone; only depends on `DataFrames`, not on `Reference`/`Data`. Reciprocal-best-hit (RBH) detection over a paralog-pair `DataFrame`, via `rbh` (percent-identity scoring: max/mean/double_max) and `rbh_ds` (dS-based scoring).

5. **`Plotting`** (`src/plotting.jl`) — currently a near-empty stub (`using CairoMakie` only, no functions yet).

6. **`Exploration`** (`src/exploration.jl`) — statistics built on `Reference`/`Data`: per-feature `coverage` fractions from `BedData`, `kde` density estimates (on coverage fractions or on `TabularData`), `quantiles` binning of samples, and a genome-wide-coverage pipeline: `calculate_frequency` (per-base overlap counts across many `BedData` measurements, parallelized over scaffolds with `Threads.@threads` using a difference-array sweep) → `feature_frequency` (projects that per-base frequency onto each feature ± a flank, strand-oriented so index 1 is the 5' end) → `gene_profile`/`mean_gene_profile` (interpolates each feature's body to a fixed bin count and averages into a metagene profile).

7. **`Modeling`** (`src/modeling.jl`) — currently a stub (imports GLM/HypothesisTests/MultipleTesting/StatsBase, no functions yet).

### Cross-cutting conventions

- Genomic intervals are consistently represented as `IntervalTree{UInt32, IntervalValue{UInt32,UInt64}}` (aliased `IntervalTreeM64`), 1-based closed `[start, end]` coordinates. BED's native 0-based half-open coordinates are converted to this convention on load (`load_bed`).
- The 64-bit `value` on an interval is always a packed metadata code; use `Reference.parse_index`/`parse_strand`/`parse_so_term` (gene features) or `Data.parse_bed_strand` (BED features) to unpack it rather than bit-twiddling inline. Methylation payloads are unpacked the same way, via `Methylation.get_depth`/`meth_percent`/`get_meth`/`get_unmeth`/`get_context`/`is_forward` — prefer `get_depth`/`meth_percent` there, as those are the stored fields and the counts are reconstructions.
- Strand is encoded identically everywhere via `BitCodes` (see above): gene and BED codes keep an 8-bit slot and methylation payloads a 2-bit one, but all three store the same values. Encode with `strand_code`, decode with `decode_strand`; don't reintroduce a local strand mapping.
- Feature/sample IDs are resolved lazily against the interned string vocab (`Reference.get_metadata_id`) rather than kept in a separate index — this is a deliberate memory/lookup-speed tradeoff repeated across `Reference` and `Data`; preserve it rather than adding ID `Dict`s, as the tradeoff is deliberate.
- Parallelism is used deliberately in a few hot paths (`add_features!`'s channel-based producer/consumer, `Experiment`'s per-sample `Threads.@spawn`, `calculate_frequency`'s `Threads.@threads` over scaffolds) — run Julia with multiple threads (`julia --project=. -t auto`) to exercise these paths meaningfully.

# Mappability from UniMap.jl

How a UniMap `mappability.zarr` becomes something BioinfoTools2 can query lazily,
and why it is converted rather than read in place. Code: `src/data/mappability.jl`
(`Data.Mappability`), `ext/BioinfoTools2ZarrExt.jl` (`convert_unimap`), and the
mappability methods in `src/exploration.jl` and `src/data.jl`.

## The workflow

```julia
using BioinfoTools2, Zarr          # Zarr only for converting
species = Species("zebrafish")
add_features!("GRCz12tu.gff", species.genome)

convert_unimap("scores/mappability.zarr", "mappability_cache", species.genome;
               channels = [(150, 300, 2), (150, 300, 10)],
               loci = "simulation/loci.tsv")      # optional: exact genome check

data = read_mappability("mappability_cache"; L = 150, F = 300, mapq_cutoff = 2)
frequency = feature_frequency(species.genome, :gene, data)   # MappabilityFrequency
mean_gene_profile(frequency); body_mean(frequency); tss_window(frequency)
low = mappability_intervals(data, species.genome; below = 0.5)  # BedData
clean = mask_calls(methylation, data; below = 0.9)              # MethylationData
```

A **channel** is one grid point `(L, F, mapq_cutoff)`: read length, fragment
length, and the MAPQ a fragment needed to count. Each is converted and read on
its own.

## What UniMap writes, and why it isn't queried directly

`mappability.zarr` holds exact per-base counts for every channel:
`n_simulated[pos, F, L]`, `n_correct[pos, F, cutoff, L]` and
`n_ectopic[pos, F, cutoff, L]`. It stays the archive: exact, and readable from
Python/xarray. Its layout suits writing it, not BioinfoTools2's queries:

- **Loci are laid end to end** along one position axis, located through a JSON
  `loci` attribute. There is no index by scaffold and coordinate.
- **Overlapping loci are stored twice**, with identical values on shared bases.
- **A chunk holds every channel** (16,384 positions × all F × all cutoffs). One
  channel across many genes would decompress ~78× what it uses.
- **Blosc is not thread-safe:** Blosc.jl decompresses through global state, and
  BioinfoTools2 reads one scaffold per task.
- **Gene IDs differ between the packages.** UniMap strips NCBI's `gene-`;
  `sanitize_id` strips Ensembl's `gene:`. Matching is therefore by coordinates.

## The cache

```
mappability_cache/
  manifest.toml       format, source store, genome_checksum (FASTA sha256), grid,
                      n_loci, log_ratio_pseudocount, compress, genome_match,
                      channels converted and their run counts
  scaffolds.tsv       scaffold name -> file stem (names aren't always file names)
  L150_F300_q2/<stem>.arrow   one file per scaffold per channel
```

Each file is a `StructArray{MappabilityRun}`: one 16-byte row per run of
consecutive bases with identical `(n_simulated, n_correct, n_ectopic)`,
1-based closed `start:stop`, sorted and non-overlapping. That is **lossless**:
expanding the runs reproduces the store's counts exactly, base for base
(`TESTS_mappability.jl`). A base in no run had all three counts zero (outside
every locus, or an N), and reads as unmeasured (`NaN`). Metrics come from the
counts on reading, with UniMap's formulas: `map_back = n_correct / n_simulated`,
`ectopic_rate = n_ectopic / n_simulated`, and `log2_ratio`, using the store's
pseudocount.

- **Compressed or mapped.** `compress = :zstd` (default) gives small files that
  are decompressed into memory per scaffold on reading. `compress = nothing` gives
  larger files that are memory-mapped, so nothing is resident until touched.
- **Adding channels.** Converting into an existing cache of the same store adds
  only the channels it lacks, with that cache's compression. `overwrite = true`
  replaces a cache. A directory that isn't a cache is never touched.
- **Partial conversions.** The manifest is written last, so an interrupted
  conversion leaves a cache that `read_mappability` refuses (a new one) or that
  doesn't list the half-written channels (an added one).

## Sizes: why channels are chosen, not all converted

Runs are short wherever the counts change often. For most channels that is
0.2–1.2% of bases. But fragment starts are `step` bp apart (10 in the pipeline's
default grid), so when `F` is not a multiple of `step` the number of fragments
covering a base alternates with period `step`. With `F = 75`, a base is covered
by 7 or 8 fragments, switching every 5 bp. Those channels have a run every ~5
bases, ~20%, and are ~10× larger.

Measured on 50 zebrafish genes (1.52 Mb), and scaled to all 28,407 genes
(~1,060 Mb):

| channels | runs (% of bases) | Zstd size per channel, full genome |
|---|---|---|
| F a multiple of `step` (e.g. 100, 200, 300) | 0.2–1.2% | ~70–140 MB |
| F not a multiple (75, 125, 175 at step 10) | ~20% | ~1.3 GB |

All 156 channels of that grid would be ~55–60 GB, so `convert_unimap` makes you
name them (`channels = :all` is still accepted). Every call reads the whole store
(~30 min for zebrafish, one reader task), so convert the channels you need
together.

## Genome matching

BioinfoTools2's `Genome` holds annotation only (no lengths or checksum), so a
"matching genome" is checked through the loci. A conversion:

- **throws** if any locus's scaffold is missing from `genome`;
- with `loci = ` UniMap's `loci.tsv`, matches each locus's exact gene span and
  strand to a `:gene` feature; without it, looks for a same-strand gene inside
  each locus span;
- **throws** if no locus matches (wrong genome or annotation), **warns** if only
  some do;
- records the result, and the store's `genome_checksum`, in `manifest.toml`.

## Analysis on top

- `feature_frequency(genome, feature, data; flank, metric)` returns a
  `MappabilityFrequency` of `FeatureRuns`: each feature's padded, strand-oriented
  region kept as runs. Mappability is dense, so per-base `FeatureLevels` would
  cost ~12 bytes × every base; `FeatureLevels(runs)` expands one when needed.
  `levels` is the metric and `weights` is `n_simulated`.
- `MappabilityFrequency` and `MethylationFrequency` are both `LevelFrequency`s.
  `mean_gene_profile`, `quantile_profiles`, `tss_window` and `body_mean` take
  either. `gene_profile` and the window sums visit each base as the per-base path
  does, so results equal those on the expanded `FeatureLevels` exactly.
- `mappability_intervals` merges low runs into a `BedData`; `mask_calls` filters
  methylation calls in one pass per scaffold; `binned_mappability` gives
  `n_simulated`-weighted means per fixed bin.

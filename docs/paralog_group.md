# `ParalogGroup` reference

Storage conventions behind `Homologs.Paralogs.ParalogGroup`, split out of the
constructor's docstring.

## Input columns

`pairs` has 2-8 columns. The first two are required and taken **by position**;
the rest are optional and matched to the relation matrices **by name**
(case-insensitive):

| # | Column             | Type     | Becomes                                    |
|---|--------------------|----------|--------------------------------------------|
| 1 | query ID           | `String` | a vertex                                   |
| 2 | subject ID         | `String` | a vertex                                   |
| 3 | `dN`               | `Real`   | `pg.dN`                                    |
| 4 | `dS`               | `Real`   | `pg.dS`                                    |
| 5 | `id_subject_query` | `Real`   | `pg.id_subject_query`                      |
| 6 | `id_query_subject` | `Real`   | `pg.id_query_subject`                      |
| 7 | `lca`              | `String` | `pg.lca` (a code into `pg.lca_labels`)     |
| 8 | `lca_depth`        | `Real`   | `pg.lca_depth`                             |

A column matching no name is ignored.

## Matrix layout

Every matrix is `n_genes × n_genes` over the group's linear gene indexing
(`pg.id_to_index`), oriented so that the axis you look up by is the column-major
fast axis:

- `topology` — symmetric `Bool` adjacency. A pair sets both `[a, b]` and
  `[b, a]`, so the relation is undirected. **`topology` is the authority on
  which pairs exist**; the other matrices only carry values.
- `dN`, `dS`, `lca_depth`, `lca` — symmetric by definition, so stored **once**
  in the upper triangle: genes `a`, `b` live at `[min(a, b), max(a, b)]`.
- `id_subject_query` — %ID subject → query with **subjects on the columns**
  (`[query, subject]`), so `M[:, s]` gathers subject `s`'s scores in O(nnz).
- `id_query_subject` — %ID query → subject with **queries on the columns**
  (`[subject, query]`), so `M[:, q]` gathers query `q`'s scores in O(nnz).

`lca` holds a 1-based code into `pg.lca_labels` rather than a value; read it back
with `lca_label`. Code 0 is unreachable, so an unstored cell means "unknown".

An explicit zero is *kept*: a dN or dS of 0 between recent duplicates survives
construction, sub-group slicing and `rbh`. Note, though, that the `dN_graph` /
`dS_graph` views symmetrise with `m + permutedims(m)`, which prunes zeros, and a
weighted graph cannot hold a zero-weight edge — so such a pair has no edge there.

## `lca_depth`

The depth of the speciation event preceding the pair's duplication, higher being
more recent. A paralog group lives in **one** species, whose ancestors form a
single lineage, so every duplication maps onto that one chain and depths *are*
comparable across pairs — equal depth means the same ancestral lineage.

Number the root 1. A 0 is reserved for "not recorded", and `rbh` refuses to rank
by it.

## Which genes are indexed

Only genes that end up **weakly connected**. A row is kept when both IDs resolve
against the genome and the pair shows some similarity, and a gene is indexed only
if it appears in a kept row. Two rules drop rows up front, each warning with a
count:

- a `NaN` in any matched numeric column — `NaN != NaN` would break the symmetry
  the `*_graph` views check;
- a 0 in *every* matched %ID column, i.e. no alignment in either direction. A
  pair with no %ID column at all is taken on trust.

Repeated or reciprocal entries for one cell collapse to a single value.

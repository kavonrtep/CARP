# Feature request (dante_ltr): deterministic library clustering

**Status:** ✅ RESOLVED in **dante_ltr 0.5.2.0** (petrnovak channel); CARP's
`envs/tidecluster.yaml` pin bumped 0.5.1.1 → 0.5.2.0.
**Affects:** `dante_ltr_to_library` → `utils/mmseq_clustering.R` (and, second-order, `utils/extract_fasta.R`).
**CARP tracking:** part of `fix/deterministic-library-clustering`.

## Summary

`dante_ltr_to_library` clusters the extracted transposable-element sequences
(`TE_all.fasta`) with `mmseqs easy-cluster`, and takes the cluster
representatives verbatim as the LTR-RT library. **`mmseqs easy-cluster` is
order-sensitive**: the same set of input sequences presented in a different
order yields a *different set of representative consensi and a different cluster
count*. Because `TE_all.fasta` is emitted in the row order of `DANTE_LTR.gff3`
(rtracklayer `import` preserves file order), and that order is not canonical
(see below), the LTR-RT library — and therefore the downstream RepeatMasker
annotation — is not reproducible run-to-run / across machines.

This is the same class of bug, with the same fix, that TideCluster already
addressed in its comparative analysis (TideCluster `changelog.md`, 1.13.1,
issue #4: "Deterministic comparative analysis").

## Evidence (measured in CARP)

- `mmseqs easy-cluster` on a fixed LINE input, same file, at 1 and 16 threads →
  **byte-identical** representatives. So it is **not** thread-related.
- The *same* sequences **shuffled** → **~19 % of representatives differ** and the
  **cluster count changes** (e.g. 2250 → 2260). Re-sorting the input into one
  canonical order before clustering restores byte-identical representatives.
- On one genome (GCA_964200825.2) the Ty1_copia/Angela library consensi differ
  ~49 % between pipeline runs; the elements themselves are stable (same count),
  so the churn is representative-election over a stable set — i.e. exactly this
  order-sensitivity.
- `DANTE_LTR.gff3` row order is **not canonical** on the multi-chunk path: chunk
  results are concatenated in chunk-index order, and the chunk count depends on
  genome size / open-file limit / machine, so `TE_all.fasta` order (hence the
  clustered library) varies across environments.

## Requested fix (in dante_ltr)

Make the library a deterministic function of the input **set**, independent of
record order. Preferred: sort `TE_all.fasta` into a canonical order **by
sequence content** immediately before clustering, inside
`utils/mmseq_clustering.R` (so every caller of `dante_ltr_to_library` benefits).
An out-of-core sort keeps it cheap on large inputs, e.g.

```sh
seqkit fx2tab TE_all.fasta \
  | LC_ALL=C sort -t$'\t' -k2,2 -k1,1 -S <buf> --parallel=<threads> -T <tmp> \
  | seqkit tab2fx > TE_all.sorted.fasta
```

Sorting by **sequence** (not coordinate/ID) makes the order invariant to
upstream coordinate jitter and to chunk grouping. Optionally also expose a
`--deterministic` switch that additionally forces `mmseqs --threads 1` for a
byte-identical `*_rep_seq.fasta` (matching TideCluster's flag), though sorting
alone already fixes the representative set.

Performance note: the sort is O(n log n) over bytes the clustering reads anyway
and is disk-backed, so it is a small fraction of clustering wall-time even on
30–90 Gbp genomes with large TE sets.

## What CARP does in the meantime

CARP cannot fix this without modifying the dependency, so it does **not** work
around it inside `make_library_of_ltrs`. CARP *does* independently sort the
inputs to the clustering steps it drives itself (`dante_line`, `reduce_library`
CAP3/mmseqs, `make_tir_combined_library`), so the LTR consensi are re-clustered
deterministically at the `reduce_library` stage. Full reproducibility of the
Angela/LTR layer requires this dante_ltr fix to land as well.

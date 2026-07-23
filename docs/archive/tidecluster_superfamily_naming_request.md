# TideCluster TRC-superfamily output: naming inconsistency + request

> **RESOLVED in TideCluster 1.16.2 (2026-07-01).** All three requests below were
> implemented: `make_empty_outputs` now writes the canonical
> `<prefix>_trc_superfamilies.csv` with the `Superfamily,TRC,fallback` header
> even when empty (0 data rows), and a `write_superfamily_manifest()` was added.
> CARP pins `tidecluster=1.16.2` (`envs/tidecluster_run.yaml`); the both-name
> lookup in the `tidecluster_reannotate` rule is retained for backward
> compatibility with 1.16.1. This document is kept as the record of the finding.


CARP's superfamily-aware RM-on-TideCluster filter
(`scripts/tc_reannotate_sf_filter.py`, config
`tidecluster_reannotate_superfamily_merge`) consumes TideCluster's per-TRC
superfamily map. While integrating it we found a **naming inconsistency in
TideCluster** that makes the file's basename depend on whether any superfamily
was found.

## Finding (TideCluster 1.16.1, `tarean/compare_trc_by_blast.R`)

| case | code | file written |
|------|------|--------------|
| superfamilies **found** (normal) | line 670: `write.csv(tbl_out, paste(opt$prefix, "_trc_superfamilies.csv", ...))` | `<prefix>_trc_superfamilies.csv` (header `"Superfamily","TRC","fallback"` + rows) |
| **no** superfamilies (empty BLAST) | line 325 (`make_empty_outputs`): `cat("", file = paste(prefix, "_superfamilies.csv", ...))` | `<prefix>_superfamilies.csv` — **drops `trc_`**, 0 bytes, no header |

So the *same* logical output has two different basenames, and the empty case is
an unheadered 0-byte file. Observed directly:

- A genome with superfamilies (Merodon `GCA_937616625.2`, A. thaliana T2T):
  `TideCluster_trc_superfamilies.csv`, populated.
- A genome without (tiny_pea, 9 TRCs, no cross-TRC BLAST edges):
  `TideCluster_superfamilies.csv`, empty.

The HTML sibling is *always* `<prefix>_trc_superfamilies.html` (line 321 and
611), so only the CSV basename diverges. There is a related note referenced in
the source: `docs/tidecluster_bug_superfamily_fallback_empty_blast.md`.

## Impact on downstream consumers

A consumer that looks for the canonical `<prefix>_trc_superfamilies.csv` silently
finds nothing on no-superfamily genomes and cannot distinguish "TideCluster ran
and found none" from "file missing / renamed". CARP now probes both names and
takes the first non-empty one, but that is a workaround for a producer-side
inconsistency.

## Request to the TideCluster developer

1. **Use one basename.** `make_empty_outputs` should write
   `<prefix>_trc_superfamilies.csv` (with `trc_`), not `<prefix>_superfamilies.csv`.
2. **Always write a header**, even when empty:
   `"Superfamily","TRC","fallback"` with zero rows — so the file is always
   present with a stable schema and "empty" is unambiguous (0 data rows, not 0
   bytes / absent / differently named).
3. **(Optional) emit a small run manifest** (e.g. `<prefix>_manifest.json` or
   `.tsv`) listing the produced artefacts and their paths, so downstream tools
   key on a declared contract rather than on guessed filenames.

Keeping the `Superfamily,TRC,fallback` schema stable across releases is the other
thing CARP relies on.

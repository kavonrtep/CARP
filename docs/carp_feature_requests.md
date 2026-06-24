# CARP feature requests (from the carp_jbrowse_server portal)

Context for the CARP developer: we publish CARP repeat-annotation output to a public
**JBrowse 2** browser. The genome-wide density BigWigs CARP emits are served as JBrowse
wiggle tracks. Requests below — (1) shrink the BigWigs (they're mostly zeros), and
(2) emit two derived density BigWigs we currently have to compute ourselves so the
pipeline stays the single source of truth — **both DELIVERED in rc6** — plus a new
(3) emit a machine-readable **output manifest** (schema version + file map) so the
server stops guessing CARP's layout and can serve multiple CARP-version annotations of
one genome side by side. **All three are now delivered — FR-1/FR-2 in rc6, FR-3 in rc7.**

Observed against **CARP 0.9.0rc4** (`tiny_pea`, ~50 Mbp). Numbers scale with genome size.

> **Status — DELIVERED in CARP 0.9.0rc6.** Both requests below shipped, and the
> portal now consumes the native CARP output directly (no publish-side synthesis or
> per-family recomputation). The publish pipeline targets the rc6 density layout
> **only** — a pre-rc6 run is refused at publish and must be re-annotated (the
> `Repeat_density_by_class_bigwig/` marker dir). Code: `domain/carp_outputs.py`,
> `publish/{extract,jbrowse,preflight}.py`. Per-request notes inline below.

---

## FR-1 — Emit sparse (run-length-merged) density BigWigs

> **DELIVERED (0.9.0rc6).** All density BigWigs are now sparse (run-length-merged),
> Unified-sourced, and live under the renamed dirs `Repeat_density/`
> (whole-genome total) + `Repeat_density_by_class_bigwig/` (per class × {10k,100k}).
> The portal serves them verbatim.

### Current behaviour
The per-class / per-cluster density BigWigs (e.g. `*_100k.bw`, `*_10k.bw`, the
`Repeat_Annotation_NoSat_split_by_class_bigwig/**` set, the TideCluster
`TRC_<n>_*.bw`) are written as **fixed-width windows with one entry per window across
the whole genome, including every all-zero window**.

### Problem
Most windows are zero, so the files are dominated by zero entries. Measured on the
spike genome:

| Track | % windows = 0 | file now | with zeros merged/dropped |
|---|---|---|---|
| `TRC_22_100k.bw` (a localized satellite) | 100% | 68 KB | ~1 KB |
| `TRC_1_100k.bw` | 97% | 69 KB | ~3 KB |
| `Tandem_repeats_100k.bw` (satellite aggregate) | 46% | 81 KB | ~14 KB |
| `All_Ty3_Gypsy_100k.bw` (dense, genome-wide) | 1% | 108 KB | ~68 KB |

Satellite/per-family tracks are ~**20–60× larger than they need to be**; across the
served set we measured ~**82% reduction**. On Gbp-scale plant genomes (our target) this
is the difference between tens of MB and ~1 MB per genome, ×thousands of genomes.

### Requested change
Write the density BigWigs **run-length-merged**: collapse consecutive windows that have
the **same value into a single interval** `[start_of_first_window, end_of_last_window,
value]`. BigWig natively supports variable-width intervals, so this is lossless and
needs no consumer change. In particular this merges long zero-runs into one interval (or
omit all-zero runs entirely — JBrowse renders absent data as zero on a 0-baseline).

Key point: **only merge windows whose values are exactly equal** — keep full per-window
resolution wherever the density actually changes. The displayed track is byte-identical
where it's non-zero; only the redundant repeated/zero windows are coalesced.

### Acceptance criteria
- Non-zero density values are unchanged at every position.
- Adjacent equal-value windows (incl. zeros) are emitted as one interval.
- Applies to all density BigWig outputs (per-class, per-superfamily, TideCluster
  per-cluster, totals — both window sizes).
- (Nice to have) the BigWig header summary (`mean`/`nBasesCovered`) is computed
  consistently; we don't rely on it, but other consumers might.

---

## FR-2 — Emit two derived density BigWigs CARP doesn't currently produce

> **DELIVERED (0.9.0rc6).** Both gaps are filled, so the publish step no longer
> computes either one — it serves CARP's output directly (the `merge_bigwigs` TIR
> synthesis and the per-family density recomputation were removed):
> - **FR-2a** — CARP now emits a native `Class_II.Subclass_1.TIR_100k.bw` rollup;
>   the portal serves it as the TIR row in the by-class overview.
> - **FR-2b** — superseded by a *unified* per-family split. CARP emits per-family
>   sparse BigWigs for two flavours under
>   `Tandem_repeats_unified_split_by_family_bigwig/100k/` (RM ∪ structural — shown
>   by default) and `Tandem_repeats_TideCluster_split_by_family_bigwig/100k/`
>   (TideCluster structural — published, off by default), named `TRC_<n>_100k.bw`.
>   The portal copies the top-N of each verbatim rather than splitting the
>   RepeatMasker tandem GFF3 itself. (The standalone RepeatMasker per-family
>   density track this asked for is dropped in favour of the unified one.)

We presently compute these two at publish time because they're missing from CARP output.
We'd prefer CARP to be the single source of truth (reproducibility + so they get FR-1's
compression for free). Both are "split an annotation by a class/family and write a
per-group density BigWig", which CARP already does for other classes — these are just
two gaps.

### FR-2a — A TIR-level density rollup
- **Now:** CARP emits per-superfamily TIR BigWigs
  (`Class_II.Subclass_1.TIR.EnSpm_CACTA_100k.bw`, `…MuDR_Mutator…`, `…PIF_Harbinger…`,
  `…Tc1_Mariner…`, `…hAT…`) but **no single TIR-level rollup**, whereas it *does* emit
  class rollups for others (e.g. `All_Ty1_Copia_100k.bw`, `All_Ty3_Gypsy_100k.bw`,
  `Class_I.LTR.Ty1_copia_100k.bw`).
- **We currently:** sum the five per-superfamily TIR BigWigs window-by-window into one
  `TIR_100k.bw`.
- **Request:** emit a `Class_II.Subclass_1.TIR_100k.bw` (and the matching `_10k`)
  density rollup directly, analogous to the existing superfamily/class rollups, so TIR
  reads as one row alongside Copia/Gypsy/LINE/Helitron.

### FR-2b — Per-family density for the RepeatMasker tandem pass
- **Now:** CARP emits per-cluster density BigWigs for the **TideCluster** pass
  (`TideCluster/default/TideCluster_clustering_split_files_bigwig/100k/TRC_<n>_100k.bw`)
  but **not for the RepeatMasker remasking pass** (`Tandem_repeats_RepeatMasker.gff3`,
  whose features also carry `Name=TRC_<n>`). Only the combined `Tandem_repeats_100k.bw`
  exists for that pass.
- **We currently:** split `Tandem_repeats_RepeatMasker.gff3` by `Name=TRC_<n>` and
  compute a per-family fraction-covered density BigWig for the top families, so the
  browser can show the RepeatMasker satellite distribution per family next to the
  TideCluster one.
- **Request:** emit per-family density BigWigs for the RepeatMasker tandem pass — the
  same split-by-family + per-window density CARP already does for the TideCluster pass,
  applied to `Tandem_repeats_RepeatMasker.gff3` (suggested naming mirroring TideCluster,
  e.g. `Tandem_repeats_RepeatMasker_split_files_bigwig/100k/TRC_<n>_100k.bw`). This lets
  us compare the TideCluster calls vs. the genome-wide RepeatMasker remasking per family.

### Acceptance criteria
- FR-2a: a TIR-level density BigWig at both window sizes, values = the per-window union
  density of the TIR superfamilies (equivalently the sum where superfamilies don't
  overlap), consistent with the other class rollups.
- FR-2b: one density BigWig per RepeatMasker tandem family (`TRC_<n>`), same window
  binning and fraction-covered value semantics as the TideCluster per-cluster BigWigs.
- Both benefit from FR-1 (sparse output) automatically.

---

## FR-3 — Emit a machine-readable output manifest (schema version + file map)

> **DELIVERED (0.9.0rc7).** CARP now writes `carp_manifest.json` at the output root on
> every run, via the wrapper (`run_pipeline.py`): `exit_status:"running"` before
> snakemake, flipped to `completed`/`failed` after — so it is present on success **and**
> failure. It carries `manifest_version`, `schema_version` (currently **"2"** — the rc6
> density layout; adding the manifest is additive and did **not** bump the schema),
> `carp_version`, `carp_git_ref`, `produced_at`, `exit_status`, and an `outputs` map of
> stable logical name → path-relative-to-output-root. It is a **separate artifact** from
> `run_provenance.json` and is **not** a completion gate. Code: `scripts/manifest.py`;
> `schema_version` changelog: [`docs/output_schema.md`](output_schema.md). The portal can
> resolve files by logical name and key its extractor on `schema_version` instead of
> sniffing marker directories.

### Current behaviour
CARP's output layout is **implicit**. The server discovers what a run produced by
hard-coding paths and sniffing marker directories — e.g. it treats the presence of
`Repeat_density_by_class_bigwig/` as "this is the rc6 layout" (`domain/carp_outputs.py`).
There is no declared, machine-readable statement of *which output contract* a run
conforms to, or *where* each consumed file lives.

### Problem
- A **breaking layout change** (rc4 → rc6 already happened; rc6 → rc7 is coming) silently
  moves or renames paths. The server can only guess which schema it got, so a mismatch
  surfaces as a deep extraction failure rather than a clean "unsupported output version"
  at intake.
- We are adding **multiple annotations of one genome from different CARP versions, served
  simultaneously** (a newer "primary" plus retained "legacy" versions). The publish step
  must pick the **right extractor per annotation**, which needs each run's output schema
  stated reliably — not sniffed.
- Reproducibility: an annotation should carry, *in its own output*, the contract it
  satisfies.

### Requested change
Write a small JSON manifest at the output root — `carp_manifest.json` — on every run:

```json
{
  "manifest_version": 1,
  "schema_version": "2",          // OUTPUT-CONTRACT version; bump ONLY on a breaking
                                  //   layout/content change (NOT on every CARP release)
  "carp_version": "0.9.0rc7",     // the CARP release/tag that ran
  "carp_git_ref": "…",            // optional, if cheap to include
  "produced_at": "2026-06-23T12:00:00Z",
  "exit_status": "completed",
  "outputs": {                    // stable logical name -> path relative to the output root
    "unified_gff3":            "Repeat_Annotation_Unified.gff3",
    "cleaned_fasta":           "genome_cleaned.fasta",
    "summary_statistics":      "summary_statistics.csv",
    "density_total_dir":       "Repeat_density/",
    "density_by_class_dir":    "Repeat_density_by_class_bigwig/",
    "tandem_unified_split_dir":"Tandem_repeats_unified_split_by_family_bigwig/100k/",
    "report_dir":              "report/"
  }
}
```

Two parts, both load-bearing:

1. **`schema_version`** — one value the server keys its extractor on. It changes **only**
   when the output layout/contents change in a way a consumer must adapt to — *not* on
   every CARP release. Keep a one-line changelog in the CARP repo for each bump.
2. **`outputs`** — a map from **stable logical names** (the things the server consumes) to
   their **relative paths** in this run. Logical names stay constant across schema
   versions where the *meaning* is unchanged; only the paths move. The server resolves
   files by logical name and never hard-codes a path.

### Acceptance criteria
- `carp_manifest.json` is written at the output root on every run (success or fail; on
  fail, `exit_status != "completed"` and `outputs` may be partial).
- `schema_version` is present and bumped on any breaking output change, with a changelog
  entry per bump.
- `outputs` lists every file/dir the portal consumes (current set: `domain/carp_outputs.py`)
  by a stable logical name → relative path.
- Logical names are stable across versions for unchanged semantics; a renamed/removed
  output shows up as a logical-name change **plus** a `schema_version` bump.

### Why this helps both sides
- **Server:** picks the right versioned extractor, refuses an unknown `schema_version`
  cleanly at intake, and serves old + new annotations of one genome together. No layout
  sniffing.
- **CARP:** the output contract becomes explicit and versioned — changes are intentional
  and visible, not silent breakages downstream.

---

### Priority
All three delivered. FR-1 (rc6) was the high-value one (file size at genome scale); FR-2
(rc6) was convenience/correctness; **FR-3 (rc7)** unblocks serving multiple CARP-version
annotations per genome and removes the layout-sniffing coupling — the server can now pick
its extractor from `carp_manifest.json`'s `schema_version` and resolve files by logical
name instead of sniffing marker directories.

# Density BigWig redesign — plan & data flow

**Status:** IMPLEMENTED in 0.9.0rc6 (all phases). This document is retained as
the design rationale and old→new path reference.
**Scope:** density BigWig outputs and their source annotations only. GFF3
annotations and the masking pipeline are unchanged except where noted.

This document is the agreed plan for reworking the density BigWig tracks so
that (a) every annotation-derived density comes from the **Unified**
annotation, (b) file/dir names stop overloading the misleading `NoSat`
prefix, and (c) each tandem-repeat family gets **two** complementary
per-family tracks (structural vs. unified). It records the target data
flow and an explicit diff vs. the previous version for downstream
consumers (JBrowse, etc.).

---

## 1. Why (problems in the current version)

1. **Partial Unified switch.** `summary_statistics.csv` and the per-class
   split GFFs/BigWigs are derived from `Repeat_Annotation_Unified.gff3`, but
   the genome-wide *total* density BigWig
   (`RepeatMasker/Repeat_Annotation_NoSat_{10k,100k}.bw`) is still built from
   `Repeat_Annotation_NoSat.gff3` — **RepeatMasker-only, tandems excluded** —
   and is in fact *orphaned* (its only consumer, `make_summary_plots`, takes
   it purely as a DAG trigger; nothing reads its values).
2. **`NoSat` name overload.** Two unrelated things share the prefix:
   - `Repeat_Annotation_NoSat.gff3` — genuinely RM-only, no tandems;
   - `Repeat_Annotation_NoSat_split_by_class_*/` — actually **Unified**-derived.
   Same prefix, opposite meaning. This is what made the README wrong.
3. **FR-2b used the wrong source.** The new per-family RepeatMasker tandem
   tracks were split from the raw, untrimmed `RM_on_TideCluster_Library.gff3`
   (neither Unified nor NoSat).
4. **Empty-input crash.** `make_rm_tandem_per_family_bigwig` fails whenever
   `RM_on_TideCluster_Library.gff3` is 0 bytes — a legitimate, expected state
   (`tidecluster_reannotate` does `: > {output}` when there is no dimer
   library). `rtracklayer::import()` dies on an empty file in
   `.sniffGFFVersion` (`argument is of length zero`). This breaks CI on any
   genome with no RM-remasked tandems (the medium fixture).

---

## 2. Decisions (locked)

- **Unified is the authoritative source** for every annotation-derived
  density track (total, per-class, per-superfamily, tandem aggregate, tandem
  per-family).
- **Keep the structural TideCluster tracks** (total + per-family). They are a
  deliberate *subset* view: structural-TC ⊂ Unified tandem union (RM ∪ TC).
- **Per-family tandem tracks come in two flavours**, so each `TRC_<n>` family
  has two BigWigs:
  - **structural** — split from the TideCluster clustering (TC only);
  - **unified** — split from `Repeat_Annotation_Unified.gff3` filtered to
    `Name=TRC_<n>` (structural ∪ similarity, conflict-resolved by tier).
- **Rename freely.** Breaking backward compatibility is acceptable **for
  BigWig files/directories only**. GFF3 outputs keep their current paths
  (including the historical `Repeat_Annotation_NoSat_split_by_class_gff3/`
  name) to avoid breaking non-BigWig consumers — see §7.
- **Fix the empty-input crash** as part of this work.
- All tracks keep FR-1 run-length-merged (sparse) encoding and the existing
  value semantics (fraction-covered, 10-bin smoothing, `_10k`=1 kb bins /
  `_100k`=10 kb bins).

---

## 3. Target data flow

```
                         Repeat_Annotation_Unified.gff3   (authoritative, non-overlapping by tier)
                                      │
        ┌─────────────────────────────┼───────────────────────────────────────────┐
        │                             │                                            │
        ▼                             ▼                                            ▼
 (all features)            calculate_statistics_and_make_groups.R         (filter Name=TRC_*)
 Repeat_density/           → per-class / superfamily / rollup GFF3s        Tandem_repeats_unified
 Repeat_density_total        + Tandem_repeats aggregate GFF3                  split_by_family GFF3s
 _{10k,100k}.bw              + Class_II.Subclass_1.TIR rollup                       │
   [Unified TOTAL]                    │                                            ▼
                                      ▼                                  calculate_density_batch.R
                            calculate_density_batch.R                     → Tandem_repeats_unified_
                            → Repeat_density_by_class_bigwig/                split_by_family_bigwig/
                                {10k,100k}/<Class|Superfam|rollup>_*.bw       {10k,100k}/TRC_<n>_*.bw
                                incl. Tandem_repeats_*.bw  (UNION aggregate)   [UNIFIED per-family]


 TideCluster/default/TideCluster_clustering.gff3   (structural TC only — a SUBSET of the union)
        │
        ├── calculate_density.R  → Tandem_repeats_TideCluster_{10k,100k}.bw          [STRUCTURAL total]
        └── split_files + calculate_density_batch.R
              → Tandem_repeats_TideCluster_split_by_family_bigwig/{10k,100k}/TRC_<n>_*.bw
                                                                                     [STRUCTURAL per-family]
```

All density values: fraction of window covered (0–1), 10-bin moving-average
smoothed, run-length-merged. `coverage()` of non-overlapping Unified features
⇒ value ∈ [0,1].

### Per-family pairing (the key new capability)
For each tandem family `TRC_<n>`:

| View | Source | Path |
|---|---|---|
| structural | TideCluster clustering (TC only) | `Tandem_repeats_TideCluster_split_by_family_bigwig/{10k,100k}/TRC_<n>_*.bw` |
| unified | Unified `Name=TRC_<n>` (TC ∪ RM-on-TC, tier-resolved) | `Tandem_repeats_unified_split_by_family_bigwig/{10k,100k}/TRC_<n>_*.bw` |

Note: RM satellites named `Satellite/Unknown` and TideHunter residuals
contribute to the **aggregate** union but have no `TRC_<n>` name, so they do
not appear in the per-family splits. This is expected and documented.

---

## 4. Target output layout (BigWig)

```
<output_dir>/
├── Repeat_density/                                   # NEW dir — genome-wide totals
│   ├── Repeat_density_total_10k.bw                   # NEW — Unified, all repeats
│   └── Repeat_density_total_100k.bw                  # NEW
├── Repeat_density_by_class_bigwig/                    # RENAMED (was Repeat_Annotation_NoSat_split_by_class_bigwig/)
│   ├── 10k/  <Class>_10k.bw, All_Ty1_Copia_10k.bw, Class_II.Subclass_1.TIR_10k.bw, Tandem_repeats_10k.bw, …
│   └── 100k/ …
├── Tandem_repeats_TideCluster_10k.bw                 # RENAMED (was TideCluster/default/TideCluster_clustering_10k.bw)
├── Tandem_repeats_TideCluster_100k.bw                # RENAMED
├── Tandem_repeats_TideCluster_split_by_family_bigwig/ # RENAMED (was TideCluster/default/TideCluster_clustering_split_files_bigwig/)
│   ├── 10k/  TRC_<n>_10k.bw
│   └── 100k/ TRC_<n>_100k.bw
└── Tandem_repeats_unified_split_by_family_bigwig/     # NEW (replaces Tandem_repeats_RepeatMasker_split_files_bigwig/)
    ├── 10k/  TRC_<n>_10k.bw
    └── 100k/ TRC_<n>_100k.bw
```

The aggregate union tandem density remains
`Repeat_density_by_class_bigwig/{10k,100k}/Tandem_repeats_*.bw` (it is the
`Tandem_repeats` class produced by the by-class splitter from Unified).

---

## 5. Old → New mapping (what downstream must repoint)

| Previous (≤ current) | New | Change | Source change |
|---|---|---|---|
| `RepeatMasker/Repeat_Annotation_NoSat_{10k,100k}.bw` | `Repeat_density/Repeat_density_total_{10k,100k}.bw` | renamed + relocated | **NoSat (RM-only, no tandems) → Unified (all repeats)** |
| `Repeat_Annotation_NoSat_split_by_class_bigwig/` | `Repeat_density_by_class_bigwig/` | renamed | none (already Unified) |
| `TideCluster/default/TideCluster_clustering_{10k,100k}.bw` | `Tandem_repeats_TideCluster_{10k,100k}.bw` | renamed + relocated | none (structural TC) |
| `TideCluster/default/TideCluster_clustering_split_files_bigwig/` | `Tandem_repeats_TideCluster_split_by_family_bigwig/` | renamed + relocated | none (structural TC) |
| `Tandem_repeats_RepeatMasker_split_files_bigwig/` (FR-2b, rc-only) | `Tandem_repeats_unified_split_by_family_bigwig/` | **removed/replaced** | **raw RM-on-TC → Unified `TRC_<n>`** |

Unchanged BigWigs: none in name (all density BigWig dirs are renamed). The
GFF3 split dir `Repeat_Annotation_NoSat_split_by_class_gff3/` is **unchanged**
(GFF3 backcompat retained).

**JBrowse / external impact:** every wiggle track that points at a density
BigWig must be repointed to the new path. The values are unchanged for all
tracks except the total track, whose meaning changes from "RM dispersed
repeats" to "all Unified repeats incl. tandems".

---

## 6. Implementation plan (rule/script level)

### 6.1 Robustness fix (do first — unblocks CI)
- `scripts/split_gff_by_name.R`: before `import()`, detect empty/0-byte/
  header-only input; if no feature lines, create no output and exit 0.
  (Guard on `file.size()==0` and on `length(import(...))==0`.)
- `scripts/calculate_density.R` and `calculate_density_batch.R` already handle
  zero-feature input; no change needed there.
- Rules that consume possibly-empty tandem GFF3s must tolerate "no families"
  (the `if ls *.gff3` guard already does; keep it).

### 6.2 New / changed Snakemake rules
1. **`calculate_bigwig_density` (modify).** Stop emitting
   `Repeat_Annotation_NoSat_{10k,100k}.bw`. Emit
   `Repeat_density/Repeat_density_total_*.bw` from
   `Repeat_Annotation_Unified.gff3` instead. Keep building the structural
   TideCluster total but write it as `Tandem_repeats_TideCluster_*.bw`.
2. **`make_bigwig_density` (modify).** Output dir →
   `Repeat_density_by_class_bigwig/` (input GFF3 dir unchanged). No logic
   change otherwise — the aggregate `Tandem_repeats_*.bw` and the
   `Class_II.Subclass_1.TIR_*.bw` rollup already come through here.
3. **TideCluster per-family structural — new dedicated rule** (decision Q1).
   Add `make_tidecluster_tandem_per_family_bigwig` that reads the existing
   `TideCluster/default/TideCluster_clustering_split_files/` and writes
   `Tandem_repeats_TideCluster_split_by_family_bigwig/{10k,100k}/TRC_<n>_*.bw`.
   `tidecluster_long` is left untouched (it still produces the split GFF3s).
   **Default run only** (decision Q2) — no short-monomer per-family set.
4. **`make_rm_tandem_per_family_bigwig` → `make_unified_tandem_per_family_bigwig`
   (replace).** New rule:
   - input: `Repeat_Annotation_Unified.gff3`, `genome_seqlengths.rds`;
   - extract features with `Name` matching `^TRC_` (selection by Name prefix,
     not a classification rewrite — allowed), split by Name into
     `Tandem_repeats_unified_split_by_family_gff3/TRC_<n>.gff3`;
   - run `calculate_density_batch.R` →
     `Tandem_repeats_unified_split_by_family_bigwig/{10k,100k}/TRC_<n>_*.bw`;
   - no-op cleanly when there are no TRC features.
   - Delete the old FR-2b rule and its `split_gff_by_name.R`-on-raw-RM usage.
5. **`rule all` / `localrules`:** update target paths to the new names; drop
   the removed ones.
6. **`make_summary_plots` (modify).** Its `bw1` trigger (old NoSat total) →
   point at `Repeat_density_total_10k.bw` (or the by-class `.done`). Update
   `make_summary_plots.R` track paths (`dir_100k_RA`, `dir_100k_TC`, the
   `TideCluster_clustering_100k.bw` reference) to the new names.
7. **`make_repeat_report.R` (modify).** Update the `split_by_class_bigwig`
   directory reference and the TideCluster per-cluster dir
   (`TideCluster_clustering_split_files_bigwig`) to the new paths.

### 6.3 Scripts touched
- `split_gff_by_name.R` — empty-input guard (and reuse for the Unified TRC
  split; optionally add a `--name-prefix TRC_` filter so the rule can pass the
  whole Unified file and split only TRC families).
- `calculate_statistics_and_make_groups.R` — unchanged (already Unified, already
  emits the `Tandem_repeats` aggregate and TIR rollup).
- `make_summary_plots.R`, `make_repeat_report.R` — path renames only.
- No conda env changes (all reuse `envs/tidecluster.yaml`); Singularity recipe
  needs no change (scripts ship via the `scripts/` directory copy).

### 6.4 Dependencies / ordering
- The Unified per-family rule depends on `Repeat_Annotation_Unified.gff3`
  (already a late artifact). The total and by-class rules already depend on
  Unified (via summary stats). The structural TideCluster per-family depends on
  `tidecluster_long`. No new cross-stage barriers.

---

## 7. Backward-compatibility & migration

- **BigWig paths change** (intended). Publish a short migration note in the
  README and CHANGELOG; provide the §5 mapping table verbatim.
- **GFF3 paths unchanged**, including `Repeat_Annotation_NoSat_split_by_class_gff3/`.
  Rationale: the user constrained breakage to BigWigs; several consumers
  hard-code the GFF3 dir. **Known wart:** that GFF3 dir keeps the misleading
  `NoSat` name while being Unified-derived. Out of scope here; flagged for a
  future major version.
- `Repeat_Annotation_NoSat.gff3` itself **stays** (still consumed by
  `make_track_for_masking` for the masking BED — a masking union, not an
  annotation summary, so leaving it is correct).

---

## 8. Validation plan
- Medium fixture (the empty-RM-tandem case) must pass end-to-end — this is the
  current CI failure; the empty-input guard is the regression test.
- Re-run the FR-1 losslessness check on a renamed track (per-base scores
  unchanged; zero-runs collapsed).
- Per-family: confirm structural and unified TRC_<n> tracks both exist and that
  unified ⊇ structural per base on a genome where RM-on-TC is non-empty.
- `snakemake -n` DAG builds; new targets scheduled once; removed targets gone.
- Reports/plots render with the new paths (no missing-file fallbacks).

---

## 9. Resolved decisions
1. **Structural TideCluster per-family** → **new dedicated rule** reading
   `TideCluster_clustering_split_files/`; `tidecluster_long` left untouched.
2. **Short-monomer per-family** → **no** — default TideCluster run only.
3. **Total track location** → under a new **`Repeat_density/`** dir
   (`Repeat_density/Repeat_density_total_{10k,100k}.bw`).
4. **FR-2b raw-RM per-family comparison** → **acceptable to drop**. Per-family
   is Unified-only; the raw TideCluster-vs-RepeatMasker per-family comparison
   is intentionally not retained.
```

# Implementation plan for `docs/todo.md`

Three independent workstreams. Each is separately shippable and testable.
File:line references below were verified against the current tree (2026-07).

Decisions taken with the user:
- **Cleanup** exposed as a config key `cleanup_intermediates: minimal | maximal |
  none` (default **minimal**, i.e. on by default), plus a `run_pipeline.py
  --keep-all` CLI override that forces `none`. Minimal-default with a heavier
  `maximal` mode.
- **LTR_RT_TR member tag**: Level-2 children get **both**
  `in_structure=LTR_RT_TR` and `member_of=<container ID>`.
- Sequencing: A (bugfix) first, then C, then B.

---

## A. Bugfix — satellite density label shows `TRC_1(?bp)`

### Root cause
`scripts/make_repeat_report.R` `discover_trc_bw_files()`:
- `:576` defaults `mon_sizes <- setNames(rep("?", N), tnames)` (the literal `?`).
- `:564` reads `TideCluster/default/TideCluster_kite/monomer_size_best_estimate_stat.csv`
  and `:580,:584` looks up column `position`. **Both are the old TideCluster
  schema.** Current TideCluster (1.16.x) writes
  `monomer_size_top3_estimats.csv` with column `monomer_size` (and a final
  `consensus_period_bp`). The old file never exists → `?` sticks.
- Label built at `:588-589`: `paste0(tnames, " (", mon_sizes, "bp)")`.

Same stale-schema bug in `scripts/make_summary_plots.R` (PDF): reads
`monomer_size_best_estimate_stat.csv` (`:195`), checks `position`/`TRC_ID`
(`:213`), reads `$position` (`:218`); fallback `NA` at `:221` drops the `(…bp)`
entirely. Note it already loads the correct file into `monomer_estimates_kite`
(`:191-192`) but never uses it for the label.

### Changes
1. `make_repeat_report.R`: point the lookup at
   `TideCluster/<run>/TideCluster_kite/monomer_size_top3_estimats.csv`, column
   `monomer_size`, aggregated per `TRC_ID` by **mode** (same aggregation the old
   code did over `position`). Do it for the default **and** short-monomer TRC
   sets (check whether short TRCs have their own
   `TideCluster/short_monomer/TideCluster_kite/`).
2. `make_summary_plots.R`: same fix; reuse the already-loaded
   `monomer_estimates_kite`.
3. If the file/row is genuinely absent, **drop the `(…bp)` suffix** rather than
   printing `?`.

### Test
Extract the monomer-size lookup into a small pure helper and unit-test it
against a `monomer_size_top3_estimats.csv` fixture (both scripts can share it),
OR validate by running the report rule on a fixture and asserting the label. The
report scripts need the bioconductor stack, so a extracted-helper unit test is
preferred (fits a light env; wire into `unit.yml` + the coverage guard).

### Notes
- No config/manifest/format change → no docs gate impact. Add a CHANGELOG bullet.

---

## C. Post-run output cleanup

### Config surface
- New config key `cleanup_intermediates: minimal | maximal | none` (default
  `minimal`). Read by `run_pipeline.py` (not a Snakemake rule).
- CLI `run_pipeline.py --keep-all` → forces `none` (overrides config).
- Recorded in `run_provenance.json` (config is already captured there).
- **Docs gate**: `run_pipeline.py` "reads" this key, so it must be added to
  `docs/configuration.md` + the README config table + `config*.yaml` examples,
  or `tests/test_config_docs.py` fails. Use the `config-docs` agent.

### Implementation
- New `scripts/cleanup_outputs.py` (standalone, unit-testable). `chmod +x`,
  `#!/usr/bin/env python3` shebang. `scripts/` is already in Singularity
  `%files`, so container sync is automatic; still confirm the pre-commit hook.
- **Keep-list derived by importing `manifest.py:OUTPUTS`** (so it can't drift),
  PLUS: every symlink target, the CI count-checked files
  (`DANTE_TIR/DANTE_TIR_final.fasta`, `DANTE_TIR/DANTE_TIR_combined.gff3`,
  `DANTE_LINE/DANTE_LINE.gff3`, `DANTE_LTR/LTR_RTs_library.fasta` + its target
  `DANTE_LTR/library/mmseqs2/mmseqs_representative_seq_clean.fasta` and
  `…_rm_compatible.fasta`), `run_provenance.json`, `carp_manifest.json`, and the
  `.classifications_validated` marker.
- **Never delete a symlink target.** Resolve every top-level symlink first and
  add its target to the keep-set.

### Hook point
`run_pipeline.py`: after the manifest/provenance finalize block (`:251-258`),
before `if rc != 0: sys.exit(rc)` (`:260`), guarded
`if rc == 0 and not is_dry_run and cleanup_mode != "none":`. `rc == 0` is the
established success gate (already used by `finalise_provenance` /
`finalise_manifest`). `--keep-incomplete` preserves partials on failure, which
`rc==0` gating respects.

### Minimal mode (default) — delete these (all verified unconsumed by any rule)
- `DANTE/DANTE_filtered.gff3.tmp.gff3` (within-rule scratch, `Snakefile:1594-1595`)
- `DANTE_TIR/DANTE_TIR.RData` (~35–53 MB debug workspace)
- `DANTE_TIR_FALLBACK/TPase_5prime_alignment.tsv`, `TPase_3prime_alignment.tsv`
- `RepeatMasker/genome_cleaned.fasta` — staged genome copy (~150–210 MB).
  **Exact path only** — never the top-level `genome_cleaned.fasta` (manifest
  `cleaned_fasta`, consumed by ~7 rules).
- `TideCluster/genome_cleaned.fasta` — staged copy from `tidecluster_reannotate`
  (`Snakefile:1102`).
- `DANTE_LTR/library/TE*.fasta`, `DANTE_LTR/library/mmseqs2/mmseqs_all_seqs.fasta`,
  `DANTE_LTR/library/mmseqs2/partitioned_s900_w1000.fasta` — **with symlink-target
  guard** (preserve `mmseqs_representative_seq_clean.fasta` / `…_rm_compatible.fasta`).
- `DANTE_LTR/LTR_RTs_library.fasta.reformatted`, `.blast.csv`, `.filtered_ids`
  (`filter_ltr_rt_library` within-rule scratch).

### Maximal mode (additionally) — after a report-deep-link check
- `TideCluster/*/TideCluster_tarean/` (~206 MB), `TideCluster_kite/`,
  `TideCluster_consensus/`. **Precondition:** verify `TideCluster_report.html`
  (and its generator) do NOT hyperlink into these before enabling by default in
  maximal. If they do, keep the linked assets.
- Workdirs: `Libraries/workdir/`, `Libraries/containment_workdir/`,
  `RepeatMasker/workdir/`, `DANTE_TIR/mmseqs_combined/` (+ `tmp/`),
  `DANTE_TIR/fallback_library_workdir/`, `DANTE_TIR/mmseqs2/`, `DANTE_LINE/mmseqs/`.
- `DANTE_TIR/log/`, `DANTE_TIR/img/`, `DANTE_TIR/report.html` (internal).

### Safety / behaviour
- Cleanup only runs on `rc == 0`, non-dry-run. On any error → keep everything
  (matches `--keep-incomplete`).
- Log every deletion (path + bytes freed) to the run log; print a summary total.
- Metacentrum: `run_pipeline.py` returns before the wrapper's scratch→NFS
  copy-back, so cleanup runs on the scratch copy (good — less to copy). Confirm
  `annotate_repeats_metacentrum*.sh` copy-back doesn't expect deleted files.

### Test
`tests/test_cleanup_outputs.py`: build a fake output tree (with symlinks +
manifest), run each mode, assert keep-list survives (incl. symlink targets) and
the right scratch is gone. Wire into `unit.yml` + coverage guard.

---

## B. LTR_RT_TR reporting & density

Current state (verified):
- Containers created in `scripts/resolve_ltr_tandems.py:139-153`
  (`structure=LTR_RT_TR`, `copy_number`, `members=`), ingested by
  `make_unified_annotation.R:167-228` `load_tier1_ltr` (drops `members`, sets
  Name/classification = lineage path, tier 1). Final L1 attrs: `structure`,
  `copy_number`, no `members`/`Parent`.
- Members = Level-2, `Parent`→container only; `structure`/`copy_number` are `NA`
  (`make_unified_annotation.R:205-213`). **No membership tag today.**
- No `LTR_RT_TR` density track; per-class split keys on `~Name` (lineage),
  L1-only (`calculate_statistics_and_make_groups.R:36-37,91-96`).
- "Complete TEs" count in the report reads `DANTE_LTR.gff3` Rank
  (`make_repeat_report.R:115-134`), per lineage; no tandem awareness.
- Container/member double-count is **already avoided** (L1-filter + `reduce()`
  union in `calculate_density*.R`); the report composition also uses L1-only.

### B3 — child membership attribute (do this first; others depend on it)
- In `make_unified_annotation.R` `meta_el`/member handling: set
  `in_structure="LTR_RT_TR"` and `member_of=<container ID>` on members of an
  LTR_RT_TR array. The container ID is known at ingest (before `finalise_output`
  re-IDs to `UA_L1_…`) — must map the pre-ID container to its final `UA_L1_…`
  and stamp `member_of` with the **final** ID (resolve during/after parent
  assignment in `finalise_output:969-999`, where the Parent is already computed;
  `member_of` == the Parent for these, so set it there).
- Add `in_structure`, `member_of` to `.META_COLS`
  (`make_unified_annotation.R:122-124`).
- **Contract updates** (all must move together):
  - `scripts/validate_unified_gff3.py`: allow `in_structure` / `member_of` on
    Level-2 DANTE_LTR members; keep `structure`/`copy_number` container-only.
  - `docs/unified_annotation_gff3_spec.md:132-141`: document the two new attrs.
  - `tests/test_unified_gff3_spec.py:84-94`: extend the member example + add
    drift checks.

### B1 — LTR_RT_TR density track (roll-up)
- In `calculate_statistics_and_make_groups.R` (near the TIR roll-up at
  `:154-157`): write `LTR_RT_TR.gff3` from L1 features with
  `structure=="LTR_RT_TR"`. L1-only → containers only → no internal double count.
- `make_bigwig_density` (`Snakefile:1774-1803`) auto-emits its 10k/100k bigwigs.
- **It is a roll-up (overlapping) view, NOT part of the disjoint partition** —
  document alongside `All_Ty1_Copia` etc. in `classification_vocabulary.yaml`
  aggregation_buckets (`:207-229`) and in the report legend. Must NOT be summed
  with the partition tracks.

### B2 — "Complete TEs" parenthetical
- In `make_repeat_report.R`: per lineage, compute M = number of complete copies
  that are LTR_RT_TR members. Source options: parse `members=` in
  `DANTE_LTR/DANTE_LTR_tandems.gff3` (count member IDs per container lineage),
  or count `UA_L2_` features whose parent `structure==LTR_RT_TR` in the unified
  GFF3. Prefer the tandems file (already lineage-tagged, no join needed).
- Render `N (M in LTR_RT_TR)` in the Complete-TEs column (`build_comp_tree`
  `:315-433`, `html_comp_table` `:1107-1159`, column `:1149`).
- **Legend/caption** (`:1157`, intro `:1393-1395`): explain what LTR_RT_TR means
  (a tandem array of complete LTR-RTs) and that M counts complete copies sitting
  inside such arrays.

### B4 — double-count verification (test only)
- Add a regression test asserting a container + its members do not double-count
  in (a) the per-class density (L1-filter) and (b) the total density (`reduce()`
  clamps ≤1). Guards live at `calculate_statistics_and_make_groups.R:36-37` and
  `calculate_density*.R` `reduce()`. Also spot-check the report composition bars.

### Test data
`run-000116` is on ceph (not reachable from the sandbox). Need either a copied
slice with LTR_RT_TR arrays, or validate on fixtures + a `tmp/` run that has
them. Confirm availability before B validation.

---

## Cross-cutting release hygiene
- New runtime scripts (`cleanup_outputs.py`) → Singularity `%files` (auto via
  `scripts/`), `chmod +x`, shebang; pre-commit hook check.
- New config key → `docs/configuration.md` + README table + `config*.yaml`
  (config-docs gate) — use the `config-docs` agent.
- New tests → wire into `.github/workflows/unit.yml` AND keep
  `tests/test_ci_test_coverage.py` green (run or exempt with reason).
- CHANGELOG bullets per change (via the `changelog` agent).

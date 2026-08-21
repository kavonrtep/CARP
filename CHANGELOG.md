# Changelog

## Unreleased

- **TideCluster 1.21.1** (from 1.20.1), closing out the cleanup thread above.
  **1.20.2** makes `tc_comparative_analysis.R` honour the run prefix for all six
  of its inputs — it hardcoded `tc_` for the two consensus FASTAs, so any run made
  with a different `-pr` (CARP uses `TideCluster`) aborted on contact.
  **1.21.0** adds `--cleanup`, a file-level purge inside
  `<prefix>_{kite,tarean,consensus}/`; CARP does not pass it, but keeps its own
  set in step with `tc_utils.CLEANUP_PATTERNS`. **1.21.1** writes
  `<prefix>_consensus/consensus_sequences_all.fasta` during **clustering** rather
  than annotation: it previously appeared only when TideCluster was given a
  reference library, so a CARP run without `tandem_repeat_library` produced a
  directory that looked complete and then aborted TideCluster's cross-genome
  comparative analysis. Every CARP run is now comparative-ready, archived ones
  included (the analysis rebuilds the pool from `TRC_*_dimers.fasta` when it is
  missing).
- **`cleanup_intermediates: maximal` no longer deletes TideCluster's three
  working trees whole.** It removed the scratch, but it also took the inputs for
  TideCluster's per-array consensus (`tc_per_tra_consensus.py`) and report
  re-render (`tc_rerender_report.py`), plus the kite monomer-size CSV that
  `make_repeat_report.R` and `make_summary_plots.R` read — so a report rebuilt
  from an archived run silently lost its `(172 bp)` monomer labels. `maximal` now
  deletes the disposable **files** inside those trees instead: the per-array
  k-mer counts, `ggmin.RData` / `monomers.RData` (written by TAREAN and never
  read back by anything), the kite periodogram, the duplicated per-array FASTAs,
  and the RepeatMasker leftovers in `TideCluster_consensus`. Verified on a 94 Gbp
  run: **40.2 GB of the trees' 44.6 GB still freed (90 %)**, 4.4 GB kept, all
  three capabilities preserved. Raised upstream as issue #3; the glob set is
  kept in step with TideCluster's own `tc_utils.CLEANUP_PATTERNS` (its
  `--cleanup`, added in 1.21.0), so one definition governs both and a
  divergence shows up as a plain diff. CARP keeps applying it itself rather
  than passing `--cleanup`, so cleanup stays one post-run step behind one
  config key.

## 1.6.0

> **Annotation outputs are unchanged.** The memory work caps *concurrency* only:
> both fixtures were run in pairs differing in `--max-memory` (~31×),
> thread count and `PYTHONHASHSEED`, and the manifest outputs came back
> byte-identical, as did `make_unified` single-batch vs multi-batch. The two
> dependency bumps do change behaviour in one visible way: **TideCluster now
> fails loudly where it used to exit 0 with missing results**, so a run that
> previously "succeeded" while a step had died will now stop.

- **`max_memory_gb` / `run_pipeline.py -m <GB>`: one knob for every memory gate,
  and a budget that reflects the job rather than the host.** CARP ships as a
  `.sif` and is usually run under PBS/Slurm, where `/proc/meminfo MemAvailable`
  reports the **node** — the kernel does not namespace it — and the cgroup limit
  that will actually kill the job normally sits on an *ancestor* job scope that
  may be unreachable from inside the container (controllers not delegated,
  `/sys/fs/cgroup` not mounted, or a scheduler enforcing by polling). Detection
  therefore could not distinguish "no limit" from "a limit I cannot see", and
  `run_pipeline.py` had no memory flag at all, so a default
  `singularity run … -c config.yaml -t 96` carried **no** memory information into
  the pipeline: every gated pool sized itself against the whole node and the run
  was OOM-killed hours in with no earlier symptom (the failure behind
  [TideCluster #6](https://github.com/kavonrtep/TideCluster/issues/6), where a
  128 GB job believed it had 1.6 TB). The new value is the job's **allocation**;
  it reaches the tools through `resources.mem_mb`, so an HPC profile or
  `--set-resources` still overrides it per rule, and it now sizes the
  `make_unified_annotation` worker pool, the three BigWig density pools and
  DANTE_TIR's concurrent CAP3 assemblies from a single place. Default `0` keeps
  every rule byte-identical to before.
- **Memory-budget resolution is now scheduler-aware, and says where the number
  came from.** `scripts/mem_utils.R` gained the rungs that survive a container
  boundary and a Python mirror (`scripts/mem_utils.py`), sharing one chain and
  one set of source labels with TideCluster 1.20.0's
  `tc_utils.memory_budget_mb()`: explicit budget → `AGENT_MEMORY` →
  `PBS_RESC_MEM` / `SLURM_MEM_PER_NODE` / `LSB_MAX_MEM_RUSAGE` /
  `SLURM_MEM_PER_CPU` × `SLURM_CPUS_ON_NODE` → the tighter of the cgroup
  headroom (walked leaf → root, v1 + v2) and `MemAvailable` → nothing. Allocation
  sources carry 80 % headroom; the measured readings are used as-is, so runs that
  detect a budget today are unaffected. Every run now prints
  `[mem] budget … from <source>` as its second line and records
  `resources.memory_budget_mb` / `memory_budget_source` in
  `run_provenance.json`; when the budget came from host memory while
  `$PBS_JOBID`/`$SLURM_JOB_ID` is set, it warns at startup instead of failing in
  hour five.
- **DANTE_LTR 0.5.4.0** (from 0.5.3.0), closing the last blind spot in the stack
  (upstream issue #13, a request from this project —
  `docs/archive/dante_ltr_max_memory_request.md`). 0.5.3.0 gates its chunk pool at
  `0.8 × MemAvailable / per-chunk peak RSS`; the per-chunk term is measured exactly
  (largest chunk probed first, `ru_maxrss`) but `MemAvailable` is not namespaced, so
  inside a `.sif` it is the node's free memory. `pool_size` is linear in that budget,
  so a 128 GB job on a node reporting ~1.6 TB oversubscribes 12.5× and the one
  safeguard between a large `-c` and an OOM stops binding — on the 94 Gbp genome a
  single chunk in flight already peaked at 108.9 GB. 0.5.4.0 adds `--max_memory <GB>`
  and the same resolution chain as TideCluster 1.20.0 and `scripts/mem_utils.py`,
  names the winning source in the pool-sizing line, and warns when the budget is
  host-wide while a scheduler job id is set. The `dante_ltr` rule now passes the
  job's allocation explicitly via `resources.mem_mb`. Output is unchanged: per-chunk
  results stay index-keyed and are concatenated in index order, so pool size cannot
  reach the GFF3 — verified on tiny_pea, where 0.5.3.0 and 0.5.4.0 on identical
  inputs give a byte-identical `DANTE_LTR.gff3` (26,150 features) and
  `DANTE_LTR_statistics.csv` apart from the `##DANTE_LTR version` header.
- **TideCluster 1.20.1** (from 1.19.0). **1.20.1** fixes two silent failures,
  one of which hits a CARP feature directly (upstream issue #7, reported from
  this project): `compare_trc_by_blast.R` named each dotplot after every member
  TRC id joined together, so a superfamily of ~61 or more three-digit members
  exceeded `NAME_MAX` and killed the script on its **first** iteration —
  superfamilies are ordered by decreasing member count, so the largest one is
  always first and the loss was total. The superfamily CSV went with it, and
  that CSV is what `tidecluster_reannotate_superfamily_merge` reads: on a genome
  with one big satellite superfamily the merge silently fell back to the plain
  per-TRC filter — the exact failure it exists to prevent — while `run_all`
  exited 0. Dotplots are now named by rank and a failed image can no longer
  discard the CSV. 1.20.1 also makes failed pipeline steps abort non-zero
  instead of leading to a report that denies its own missing results; a run that
  previously "succeeded" with missing TideCluster outputs will now fail the rule
  loudly. **1.20.0** fixes the same failure class as the work above
  (its [issue #6](https://github.com/kavonrtep/TideCluster/issues/6), reported
  from this project). TideHunter's pool size and TAREAN's thread cap were
  resolved from `MemAvailable`, so a 128 GB PBS job read a ~1.6 TB budget, ran 32
  concurrent TideHunter parts wanting ~310 GB and was OOM-killed. 1.20.0 adds
  `--max_memory <GB>` and one shared resolution chain (the same one CARP now
  uses), logs the winning source, and warns when it falls through to
  `MemAvailable` under a scheduler. The `tidecluster_long` / `tidecluster_short`
  rules pass the job's allocation explicitly via `resources.mem_mb`, so
  TideCluster never has to guess. `tc_reannotate` has no such flag upstream (its
  RepeatMasker chunk pool is CPU-sized), so `tidecluster_reannotate` passes
  nothing. On a plain host with no scheduler or container the budget is unchanged
  (`MemAvailable` × 0.8); where a source does hit, the only difference is a
  smaller pool — TideHunter and TAREAN output are thread-count-independent.
- **`dante_tir_cap3_max_memory_gb` now outranks `resources.mem_mb`** (it was the
  other way round), so the global `max_memory_gb` cannot silently override a CAP3
  budget pinned for that tool. No effect at the default `0`.

## 1.5.1

> **All BigWig tracks are unchanged** — every `.bw` verified byte-identical on
> tiny_pea and on a 2.82 Gbp genome, across all four density rules. What changes
> is the memory the density layer needs: it no longer scales with genome size.

- **BigWig density: build only the bins that are kept.** `calculate_density.R` and
  `calculate_density_batch.R` used to materialise the **whole-genome** `tileGenome`
  grid for *every* track and then discard the bins the track did not occupy. On the
  94.26 Gbp run-000156 that is 94.26 M bins per task at `step=1000` (plus a 94.26 M
  character vector for the `%in%` filter), and 3,206 such tasks ran concurrently:
  `make_tidecluster_tandem_per_family_bigwig` peaked at **568.5 GB — 74 % of a
  768 GB host** — to produce 1,603 small per-family tracks. Across those families
  the occupied territory sums to 7.25 Tbp against 151.1 Tbp of grid actually built,
  i.e. **20.8× of the tiling was thrown away**; the median family occupies one
  2.14 Gb chromosome and the commonest is a single feature on a single scaffold.
  A previous attempt at "tile the occupied sequences only" was correctly rejected,
  because `tileGenome` derives its effective tile width from the whole-genome
  total, so subset tiling shifts the bins. The new `density_track()` instead
  reproduces that phase in closed form and streams **one sequence at a time**, so
  peak grid memory is `O(longest sequence / step)` — ~2.15 M bins here, independent
  of genome size *and* of how much of the genome a track covers. Every downstream
  step was already per-seqname (`binnedAverage` is per-bin independent, the
  smoothing splits by seqname, the run-length merge never spans seqnames), so the
  output is unchanged: `tests/test_density_tiling.R` pins the new grid arithmetic
  against the real `tileGenome` over 3,068 randomised and boundary cases
  (including a genome larger than `.Machine$integer.max`), `tests/test_density_utils.R`
  pins `density_track` against verbatim copies of both implementations it replaces,
  and all four density rules were re-run on two real genomes with byte-identical
  `.bw` output. This also collapses two near-duplicate density implementations
  (`density_per_family` / `get_density2`) into one.
- **New config parameter `bigwig_max_workers` (default `0` = no ceiling) + a
  measured, cgroup-aware memory gate on the BigWig worker pool.** The density rules
  ran `mc.cores = workflow.cores` (96 on the HPC profile) with no memory bound. Per-
  task cost varies by orders of magnitude between a 1-feature family and a 52 M-feature
  class, so rather than guess it the heaviest task now runs **first, on its own**, and
  the pool is `budget × 0.8 / its peak RSS` — the budget being `resources.mem_mb`
  when a profile sets it, else the tightest of the cgroup limit (walking **up** the
  hierarchy, so a PBS/Slurm job-scope limit is honoured rather than the node's free
  memory) and `/proc/meminfo MemAvailable`. Remaining tasks are dispatched
  longest-processing-time-first so a dominant class file cannot start late and
  strand the pool. Concurrency only — each task writes its own `.bw`, so no track
  can change. `[mem]` lines now report the probe peak, the gate decision and the
  max/median worker peak.
- **`calculate_density_batch.R`: detect dead workers and missing tracks.** A
  signal-killed `mclapply` child returns `NULL`, not a `try-error`, so scanning the
  returned messages for `^ERROR` could not see it — and because these rules are
  checkpointed by a `.done` marker rather than by the `.bw` files, a silently
  missing track would have gone unnoticed. Every task is now validated positively
  (a result was returned, and the `.bw` it claims to have written exists).
- Internal: the memory-instrumentation and cgroup-budget helpers move from
  `scripts/make_unified_annotation.R` to a shared `scripts/mem_utils.R` (pure move,
  no behaviour change).

## 1.5.0

> **Annotations are unchanged** — verified byte-identical on tiny_pea end-to-end,
> and at 1 vs 3 workers and 1 vs 5 workers with batching held fixed. What changes
> is that `make_unified_annotation` now caps its own concurrency (default 8
> workers, previously up to `workflow.cores`), and that a partial worker failure
> aborts instead of silently writing a truncated annotation. Large genomes are the
> reason: on a 94.26 Gbp assembly this rule was OOM-killed after 5 h 50 m at 96
> threads, and completes in under 2 h once the pool is capped.

- **New config parameter `make_unified_max_workers` (default `8`) + a cgroup-aware
  memory gate on `make_unified_annotation`'s worker pool.** Both cap **concurrency
  only** — batch composition still comes from `threads`/`batch_size`, so neither can
  change the annotation (verified byte-identical at 1 vs 5 workers).
  **Why:** `mclapply` forks, and R's GC dirties the inherited parent heap, so each
  worker's peak RSS converges on the *parent's* regardless of its own batch size —
  measured on a 94 Gbp genome at **48.3 GB per worker for every one of 55 batches**
  (min = median = max) against a 48.4 GB parent, for batches spanning 143 Mb to
  2.15 Gb of sequence. At `threads: workflow.cores` = 96 that demanded ~2.66 TB on a
  768 GB host: 388 M minor page faults (≈1.5 TB of copy-on-write copying) *at only 4
  workers*, 15 workers OOM-killed, rule dead after 5 h 50 m. The identical work at 4
  workers took **885 s**, and produced the full annotation (69.7 M features, 80.44 %
  genome coverage, zero overlap violations). Because per-worker cost tracks the
  parent heap rather than the batch, splitting oversized batches would not have
  helped — the cap is the fix.
  The memory gate uses `budget × 0.8 / parent RSS`, where the budget is
  `resources.mem_mb` when a scheduler/profile sets it, else the tightest of the
  **cgroup limit — walking up the hierarchy, so a PBS/Slurm job-scope limit is
  honoured rather than the node's free memory** (v2 `memory.max`/`memory.current`,
  v1 `limit_in_bytes`/`usage_in_bytes`, plus the namespaced-container case) — and
  `/proc/meminfo MemAvailable`. `--mem_budget_gb` overrides. Note this rule gains
  little from wide parallelism anyway: tier resolution was 885 s of a 6,949 s run,
  the rest (65 min loading, 20 min GFF3 export) being serial.
- **`make_unified_annotation`: detect dead `mclapply` workers instead of writing a
  truncated annotation.** On a 94 Gbp genome (run-000156) 15 of 55 forked workers
  were killed; the guard only tested `inherits(r, "try-error")`, which a
  signal-killed child never returns, so the run continued, concatenated 40 list
  elements instead of millions of features, and died much later in
  `finalise_output` with a misleading `seqinfo … signature '"list"'`. Had the
  survivors concatenated cleanly it would have written an annotation missing
  **27 % of the genome and exited 0**. The check is now positive — expected result
  count, and every element carrying `GRanges` `level1`/`level2` — and names the
  failed batches.
- **`make_unified_annotation`: per-phase and per-worker RSS logging.** `/proc`-based
  `[mem]` lines at load / fork / resolve / combine / final, with each worker
  reporting its own `VmHWM` (which covers inherited copy-on-write pages, so it is
  the per-worker footprint to size a pool against) and the parent summarising
  min/median/max. The rule's benchmark file is only written on success, so a failed
  run previously left no memory data at all — which is exactly what happened on
  run-000156. Degrades silently off Linux.
- **`tests/test_resolve_tier1_overlaps.R` now tests the shipped code.** It compared
  two *local copies* of the resolver carrying a "keep in sync" comment, so it passed
  regardless of what `make_unified_annotation.R` actually did. It now extracts the
  real `trim_to_nonoverlap` / `resolve_tier1_overlaps` from the script.

## 1.4.1

> **No output change.** This release is dependency pins only — both upstream
> releases are output-neutral, verified on tiny_pea against the 1.4.0 baseline
> under a different `PYTHONHASHSEED` and thread count: all **196 manifest outputs
> byte-identical**. It is a performance release for large genomes; on small ones
> nothing changes but the wall clock of three rules.

- **DANTE_LTR 0.5.2.0 → 0.5.3.0 and TideCluster 1.18.0 → 1.19.0** (pin-only; both
  dependency sets are unchanged, both envs re-solve clean). Two upstream releases
  driven by a parallelism audit of a 94 Gbp run on 96 cores, where the pipeline
  spent **199 h of wall clock at an average of well under 10 of 96 cores**. Both
  are reported output-neutral, and TideCluster verified byte-identical to 1.18.0 on
  the deterministic outputs.
  - **`dante_ltr` was strictly serial**: `for i in range(number_of_temp_files):
    subprocess.check_call(...)` ran one `detect_putative_ltr.R` at a time, so
    `dante_ltr -c 96` used **1.0 core for 24.2 h** over 1,889 chunks. 0.5.3.0 runs
    the chunks in a pool — largest chunk first as a memory probe, the rest
    longest-first via `imap_unordered`, children at `-c 1` — sized by
    `min(cpu, 0.8 × MemAvailable / per-chunk peak RSS)`. Per-chunk outputs are still
    index-keyed and concatenated in index order, so completion order cannot reach
    the output, and the fd-budget guard from the too-many-open-files fix is retained.
  - **TideCluster `tc_reannotate` was bracketed by serial phases**: its 96-worker
    pool was correct, but the parent then re-read every chunk `.out`, built a
    `Gff3Feature` per hit, accumulated all hits and sorted — holding the stage to
    **2.3 of 96 cores over 43.3 h**. 1.19.0 parses and remaps in the workers and
    k-way merges sorted fragments in the parent; parent memory now scales with
    chunk count rather than hit count.
  - **TideCluster TAREAN was gated by one straggler**: jobs were always launched at
    `-n 1`, so on the 94 Gbp genome 441 of 445 finished within ~2 h and `TRC_1` then
    ran **~55 h alone on a single core** — 74 % of that stage — while 95 cores idled.
    1.19.0 dispatches longest-first and gives large TRCs multiple threads
    (memory-gated); `tarean.R` output is thread-count-independent.
  - Write-ups: `docs/archive/dante_ltr_chunk_parallelism_request.md`,
    `docs/archive/tidecluster_parallelism_request.md`.

## 1.4.0

> **Upgrade note — annotations change.** Every release from 1.0.0 through 1.3.0
> produced a RepeatMasker library containing **zero `Class_II` sequences**, so DNA
> transposons were annotated only where DANTE_TIR found a complete structural
> element and never by similarity. Re-running is advisable if DNA-transposon
> content matters for your assembly. Measured end-to-end on an 80 Mbp genome:
> unified `Class_II/Subclass_1/TIR/EnSpm_CACTA` **481,962 → 880,445 bp (+83%)**,
> `Class_II` total **+25%**, whole annotation 66.16% → 66.65% of the genome; every
> other `Class_II` subfamily moved exactly 0 bp and the largest non-target shift was
> `Tandem_repeats` −6.8 kb. On a 2.82 Gbp genome the library gained 23 TIR consensi
> (114 kb) and the newly-functional `Class_II` screen removed 26 LTR-RT consensi.

- **Fix: the DANTE_TIR library was empty on every run.**
  `filter_dante_tir_by_multiplicity.py` joined the GFF3 `ID`
  (`Class_II_Subclass_1_TIR_hAT_2`) against the FASTA name (`hAT_2`) verbatim.
  DANTE_TIR derives the FASTA name by stripping that classification prefix, so the
  comparison was unsatisfiable and the filter kept **0 records on every input**
  (0/7233 on a 90 Gb genome; 0/N on all 18 runs checked locally). Because
  `dante_tir_min_multiplicity` defaults to `3`, the filter runs by default, so
  `all_representative_elements_combined.fasta` was **0 bytes in every run** since the
  script was introduced (1.0.0 and all 0.9.0rc*). The rule then `touch`ed an empty
  library and exited 0, making it silent. Downstream this also emptied
  `class_ii_library.fasta`, left **zero `Class_II` sequences** in the RepeatMasker
  library, and degraded `filter_ltr_rt_library` to a plain `cp` — so the Class_II
  screen that removes DNA-transposon-like sequence from the LTR library had never
  actually run. The join now keys on `Name=` (DANTE_TIR ≥ 0.3.1, which publishes the
  FASTA name explicitly) and falls back to stripping the prefix off `ID=` for older
  outputs.
- **Fail loudly on a total join failure.** A join matching nothing is a bug, not an
  empty result; the script now reports a join census (rows parsed / joined /
  unjoined / missing `Multiplicity`) and exits non-zero when no record joins, so the
  calling rule cannot silently emit an empty library.
- **`Multiplicity` absent is no longer silently treated as 1.** DANTE_TIR < 0.3.1
  clustered before `round4()`, so round-4 elements carried no `Multiplicity` at all
  (25% of elements on a real run) and were dropped at any threshold ≥ 2 for a reason
  unrelated to their copy number. Missing is now a distinct state: still dropped by
  default, but counted, warned about, and recoverable with
  `--keep-missing-multiplicity`.
- **DANTE_TIR 0.3.0 → 0.3.1** (pin-only; same dependency set). Upstream release
  driven by this bug: `DANTE_TIR_final.gff3` is now byte-deterministic at fixed
  input+version (the fragment order handed to CAP3 no longer inherits mmseqs'
  thread-dependent output order); every element carries `Multiplicity`
  (`cluster_tir_sequences()` re-runs on the final set after `round4()`); and each
  `sequence_feature` row carries `Name=`, an explicit GFF3↔FASTA join key. Two
  deliberate result changes: `Multiplicity` values rise because they are now computed
  over the complete element set (maize: elements passing a `>= 3` floor 806 → 978),
  and the element set shifts slightly with the fixed fragment order (maize: −4 curated
  loci of 1927, −0.2 points) — the cost of reproducibility.
- **New test** `tests/test_filter_dante_tir_multiplicity.py` (wired into `unit.yml`):
  covers the join on both the raw and canonicalised header forms, the `Name=` path,
  numeric threshold comparison, missing-`Multiplicity` handling and the fail-loud
  behaviour. The `output_medium` fixture cannot catch this regression — all its
  records are `Multiplicity=1`, so the broken and correct results are both 0 kept.
- **The `filter_ltr_rt_library` Class_II screen runs again.** Because
  `class_ii_library.fasta` was empty, the rule took its `[ ! -s ]` branch and degraded to
  a plain `cp` — the screen that removes DNA-transposon-like sequence from the LTR
  library (a documented "Critical step") had never actually executed. Validated on a
  2.82 Gbp genome: 58 blastn hits, **26 LTR-RT consensi removed** across nine lineages
  (Tekay, Tork, Ale, Reina, SIRE, CRM, Retand, Ikeros, Alesia), e.g. a `Ty1_copia/SIRE`
  consensus matching `PIF_Harbinger` at 92.8% identity over 332 bp (e-140). Small
  genomes may see no hits (the 80 Mbp test genome returns zero).

## 1.3.0

- **TE-derived TRC detection: domain-rhythm gate.** A Tier-3 satellite is tagged
  `TE_origin` (TE-derived, satellite wins the region) only if — in addition to the
  existing coverage test — the TE's DANTE domains recur **through the tandem**:
  `identify_te_derived_trcs` now tiles each array at the TRC's true monomer period
  and requires `domain_occupancy ≥ 0.5` across `frac_arrays_in_rhythm ≥ 0.5`
  (`te_domain_rhythm`). This distinguishes a genuine TE-derived tandem (the TE
  repeats as the monomer) from a plain satellite merely *interrupted* by a few TE
  insertions — coverage alone cannot, since a degraded genuine tandem is also only
  partially covered, whereas domain rhythm is degradation-robust. A candidate that
  fails is not tagged and is split around the interrupting TEs by normal tier
  resolution. Calibrated on 3 genomes (100 TE_origin TRCs).
- **Correct tandem period source: `trc_table.tsv`.** The per-TRC monomer now comes
  from TideCluster's report table (`monomer_tarean` → `monomer_kite` founder
  fallback) via `--tc_trc_table_default`/`--tc_trc_table_short`, replacing the kite
  `monomer_size` CSV. The kite peak collapses to SSR sub-periods (e.g. 79 bp for a
  genuine 13134 bp TIR-derived monomer), which mis-tiled the occupancy test *and*
  showed the wrong monomer in the CSV; `trc_table.tsv` is also `cleanup_intermediates:
  maximal`-safe (the kite tree is not).
- **New output: `Repeat_Annotation_Unified.te_derived_trc.csv`** — a per-TRC
  summary of the tandem-repeat clusters the unified annotation tags `TE_origin`.
  Written by `make_unified_annotation.R` beside the unified GFF3, columns `trc_id,
  run, n_arrays, total_array_bp, monomer_length_bp, te_classification,
  te_origin_structure, protein_domains, n_complete_elements, n_expected_monomers,
  complete_bp_fraction, domain_occupancy, frac_arrays_in_rhythm` (header-only when
  no such TRCs exist), registered in `scripts/manifest.py` OUTPUTS (determinism
  gate). See [`docs/output_reference.md`](docs/output_reference.md).
- **HTML report: "Tandem repeats derived from transposable elements" table** in
  the Tandem repeats (TideCluster) section of `make_repeat_report.R`, populated
  from the new CSV (now a rule input) and shown only when TE-derived TRCs are
  present.
- **Fix (determinism + quality): tier-1/2 overlap resolution corrupted by the
  te_sat pre-pass.** The TE-derived-satellite pre-pass trimmed `t1`/`t2` against the
  satellite regions with `trim_to_nonoverlap` **before** those tiers' own overlap
  resolution ran. `trim_to_nonoverlap` does `disjoin(c(lower, higher))`, which
  silently decomposes `lower`'s *internal* overlaps into disjoint pieces — so the
  subsequent greedy `resolve_tier1_overlaps` saw an already-non-overlapping `t1` and
  did nothing. Because `te_sat` is non-empty only in a batch that contains a
  TE-derived-TRC array (and batching is thread-count-dependent), overlapping
  structural TEs were resolved greedy-longest-first in some runs and naively
  disjoined (over-fragmented, arbitrary metadata) in others — 816/2.82 M tier-1
  features differed run116 threads 1 vs N, and the result was also lower quality
  when it triggered. Fixes: (a) run `resolve_tier1_overlaps` **before** the te_sat
  trim so `t1` is non-overlapping when carved; (b) fold the te_sat carve of `t2`
  into Step 2 (trim vs `reduce(level1)` = te_sat ∪ tier-1, disjoined once by
  `resolve_within_tier`) instead of a separate pre-pass trim — identical to today
  when there are no TE-derived TRCs; (c) extend the greedy's deterministic tie-break
  to `(seqname,start,end,strand,classification,source_tool)`; (d) canonicalise the
  GFF3 mcols column order before export, so column-9 attribute order (e.g. `ID`
  position) is byte-identical across thread counts, not just data-identical. Full
  A/B on run116: 816 → 0 differing features, byte-identical unified GFF3. (Supersedes
  the earlier `(seqname,start,end,strand)`-only tie-break; see
  `docs/archive/tier1_resolution_determinism_audit.md`.)
- **Fix (determinism): `trim_to_nonoverlap` always decomposes internal overlaps.**
  It had two `return(lower)` short-circuits that left `lower`'s internal overlaps
  intact, while the main path disjoined them — and which ran was batch-dependent
  (`higher = reduce(level1)`). So same-tier features overlapping only each other
  (e.g. two DANTE domains far from any tier-1) came out with min-index metadata /
  kept strand in one batching but LCA / strand `*` in another (a tier-2/4/5 twin of
  the te_sat pre-disjoin bug; run116 never hit it — tier-1 is dense — but the
  `unified_multibatch` fixture does). It now always decomposes `lower` (disjoining
  `lower` alone when nothing overlaps `higher`, so the fast path stays cheap),
  matching what a single-batch `threads=1` run already produced — run116 output is
  unchanged, the multi-batch fixture is now byte-identical.
- **Fix (determinism): TE-derived TRC decision is now global, not per-batch.**
  `identify_te_derived_trcs` ran once per processing batch, and batch composition
  is thread-count-dependent (threads=1 → one batch), so whether a multi-sequence
  TRC was tagged `TE_origin` could differ run-to-run across thread counts. It now
  runs once over every array of each TRC genome-wide. This also required a **robust
  LCA**: the shared-lineage test now ignores rare stray overlaps (a class must
  cover ≥ `TE_ORIGIN_LCA_MIN_SHARE` = 10 % of the covered bp), so a single inserted
  element of an unrelated family (e.g. 1 Ty1_copia/SIRE among 73 Ty3_gypsy/CRM)
  no longer collapses the LCA and drops a genuine TE-derived TRC.
- **CI (determinism): multi-batch gate for `make_unified_annotation`.** The
  full-pipeline determinism gate runs only single-batch fixtures, so it cannot
  catch per-batch resolution bugs (like the te_sat pre-disjoin above). New test
  `tests/test_unified_multibatch_determinism.sh` drives make_unified directly on
  `tests/fixtures/unified_multibatch/` — a fixture whose `.fai` declares 3× 2 Mb
  sequences (make_unified reads only lengths) so `--batch_size` forces a 3-batch
  split at `--threads 3` while `--threads 1` stays single-batch — and asserts the
  unified GFF3 is byte-identical across the two. Wired into the `determinism` job;
  FAILS against the pre-fix code, PASSES on the fix.
- **CI: runner-deps guard.** New `tests/test_ci_runner_deps.py` (wired into the
  `unit` job and the pre-commit hook) asserts every test is wired into a CI job
  whose environment actually provides the packages it imports — deriving each
  job's package set from its `create-args`. Catches "a test runs on a runner
  without its dependency" (e.g. a `library(GenomicRanges)` R test in the
  lightweight `carp-unit` env) at commit / PR time instead of as a red CI job.
  Also moved the two GenomicRanges R unit tests into a dedicated bioconductor
  `unit_r_bioc` job so they actually run.

## 1.2.1

- **Fix: non-deterministic RepeatMasker library order (Class_I/LINE base counts
  drifted run-to-run).** `reduce_library_size.py` streams mmseqs cluster
  representatives in mmseqs' native, thread-scheduling-dependent order (to stay
  byte-identical to its R reference), so `combined_library_reduced.fasta` — and
  the containment-reduced library RepeatMasker actually indexes — had stable
  **content** but an unstable **order**. Because `makeblastdb` assigns OIDs by
  input order and RepeatMasker tie-breaks equal-scoring HSPs by OID, ~600 bp of
  `Class_I/LINE` masking flipped between two full-genome runs of identical input.
  The `reduce_library_containment` rule now canonically sorts its output
  (`canonical_sort_fasta.py`) — the single choke point immediately before
  RepeatMasker — so masking is reproducible. Audited every other library→BLAST
  path: the TideCluster dimer library is already canonically ordered
  (`reduce_dimer_library.py`), all source libraries are sorted at
  `concatenate_libraries`, and the blastn filters make order-independent
  set decisions, so this was the only affected path.

- **Regression guards for the library-order fix.** The `determinism` CI gate now
  runs over both the `small` and `medium` fixtures (medium is the only fixture
  that builds a multi-representative mmseqs LINE library, i.e. the path that
  drifted). A new deterministic check, `scripts/assert_library_canonical_order.py`
  (unit-tested), asserts the RepeatMasker library equals its own canonical sort
  and runs in the `fixture` job — it fails 100% if the sort is ever removed,
  independent of whether a two-run diff happens to trigger the drift.

- **Determinism is now enforced in CI, not just intended.** A `determinism` job
  in `.github/workflows/pipeline.yml` runs the small fixture twice under a
  different `PYTHONHASHSEED` and thread count and fails the build if any
  manifest data output differs between the two runs
  (`scripts/assert_run_determinism.py`, comparator unit-tested). A deeper
  full-size two-run check is available on demand
  (`.claude/skills/check-determinism/`), and a new "Determinism" section in
  CLAUDE.md documents the rules that keep the pipeline reproducible.

## 1.2.0

- **Dependency bump: DANTE_TIR 0.2.8 → 0.3.0** (`envs/dante_tir.yaml`), with a new
  CAP3 memory guard. 0.3.0 pulls three new transitive dependencies —
  `bioconductor-rsamtools`, `bioconductor-genomeinfodbdata` and `numpy` — all
  resolvable from the channels the env already lists, so only the version pin
  changed. 0.3.0 also sizes its CAP3 pool to a memory budget (it auto-detects the
  cgroup limit, falling back to node RAM, at 60 %); to keep that safe on
  schedulers that enforce no cgroup memory limit, the `dante_tir` rule now
  declares `resources: mem_mb` and passes `--cap3_max_memory` = 60 % of the
  allocation, with a new `dante_tir_cap3_max_memory_gb` config knob (default `0` =
  auto) as the fallback.

- **Bounded memory in the shared all-vs-all alignment engine (very high-copy
  families).** The engine used by both `dante_line` and `dante_tir_fallback` no
  longer materialises the whole run in RAM: `run_all_vs_all_alignment` streams
  each clustering group's records to disk as the group finishes (instead of
  accumulating every record across all groups before one write),
  `_compare_sequences` bounds the number of in-flight worker futures, and
  `analyze_alignment_lengths` keeps only a size-N min-heap per group rather than
  every alignment length. Peak memory now scales with the largest group / a fixed
  window, not the whole family, so very high-copy LINE / TIR families no longer
  OOM. Output is unchanged (verified byte-identical on single- and multi-group
  paths; the length-threshold rewrite is equivalence-tested against the old
  algorithm over randomised inputs).

- **`dante_line` failures no longer silently drop the LINE layer.** The rule
  wrapper previously turned *any* `dante_line.py` non-zero exit into empty LINE
  outputs, so a genuine crash looked like a successful run with no LINEs.
  `dante_line.py` now exits with a distinct code (3) only for the benign
  "too few LINE features / no valid patterns" case; the Snakefile creates empty
  outputs solely for that code and **fails loudly** on any other error. Also
  hardens the shared all-vs-all alignment engine used by both `dante_line` and
  `dante_tir_fallback`: intermediate tables are written atomically (a killed run
  can no longer leave a truncated file that a checkpoint trusts), the mmseqs
  clustering checkpoint is gated on a completion marker, the alignment TSV rows
  are emitted in a deterministic order (checksum-verifiable), and the deprecated
  `tempfile.mktemp` calls were replaced with `mkstemp`. No annotation output
  changes.

- **Reproducible RepeatMasker under culling: pin the culling `rmblastn` to a
  single thread.** NCBI BLAST's `-culling_limit` is non-deterministic across
  threads (per-thread HSPs merge in completion order, then equal-scoring
  enveloped HSPs are culled in that order), and RepeatMasker launches `rmblastn`
  with `-num_threads 4` even at `-pa 1` — so culling + multithreading jittered
  ~0.03–0.07% of hits run-to-run (alignment boundaries, tied sub-consensus
  choice, divergence %); family-level classification totals were already stable.
  The RMBLAST_DIR shim (`scripts/rmblast_culling_shim.py`) now strips
  RepeatMasker's `-num_threads N` and forces `-num_threads 1` on real search
  calls (in addition to appending `-culling_limit N`; `-version` still passes
  through), making culling bit-reproducible. This fixes **both** RepeatMasker
  layers through the one shared shim — the main `repeatmasker` rule and the
  `tidecluster_reannotate` rule (whose internal RepeatMasker inherits
  `RMBLAST_DIR`; TideCluster itself is not patched). Parallelism is recovered at
  the chunk level (`scripts/repeatmasker_wrapper.py` now dispatches chunks
  largest-first via LPT scheduling, `Pool.starmap(chunksize=1)` — dispatch order
  only; the merged `.out` is byte-identical), and on large genomes this also
  removes the previous 16-workers × 4-threads CPU oversubscription. Only affects
  runs with culling enabled (`repeatmasker_culling_limit` > 0 or
  `tidecluster_reannotate_culling_limit` > 0; culling-off runs were already
  deterministic) — no config parameters changed. **NOTE:** outputs shift slightly
  (~0.3% of hits) versus prior releases and become reproducible. Verified on a
  190 Mb genome: two runs byte-identical (were not before), RepeatMasker
  wall-time unchanged-to-faster (322 s vs 334–379 s baseline), and TideCluster
  reannotate byte-identical across runs (was not before).

## 1.1.5

- **Dependency bump: TideCluster 1.17.0 → 1.18.0** (`envs/tidecluster_run.yaml`).
  A large-genome performance / scalability update, **byte-identical to 1.17.0 on
  the deterministic outputs** (no behaviour, CLI, or output change): TAREAN/rDNA
  now use a random-access `IndexedFasta` reader instead of loading the whole
  genome into RAM (which was OOMing large assemblies), O(1) TideHunter
  coordinate remap, linear duplicate-array filtering, `tidehunter --long`
  streams rounds to disk, and a memory-gated parallel per-part TideHunter —
  directly relevant to the 30–90 Gbp genomes. The self-contained report / maximal
  purge safety from 1.17.0 was reconfirmed end-to-end on the small fixture (0 of 9
  report PNG references missing after a maximal purge).

## 1.1.4

- **Reproducible repeat libraries: canonically sort every clustering input.**
  The greedy clustering CARP drives — `mmseqs easy-cluster` and CAP3 — is
  order-sensitive: the *same* sequences in a different order yield different
  representative consensi and a different cluster count (measured: ~19% of mmseqs
  reps and the cluster count itself change under a pure input shuffle; ~3.7% for
  CAP3; mmseqs is otherwise thread-count-invariant). Upstream tools hand CARP the
  same sequences in a run-varying order (DANTE's ~0.1% boundary jitter;
  environment-sensitive chunk grouping), which this amplified into large
  run-to-run library churn — on one genome the Ty1_copia/Angela library differed
  ~49% between runs and the LINE quantity moved, while total repeat content stayed
  constant (a pure reclassification at the LINE↔Ty1_copia RT-ambiguity boundary).
  Fix: a new out-of-core helper `scripts/canonical_sort_fasta.py` (GNU `sort`,
  disk-backed, so it stays cheap on the multi-GB intermediates of 30–90 Gbp
  genomes) sorts each clustering input by **sequence content** immediately before
  clustering — invariant to upstream record order — at every site CARP invokes
  the tool itself: `dante_line.py` (LINE), `concatenate_libraries` (one sort of
  `combined_library.fasta`, which canonicalises every `reduce_library` per-class
  CAP3/mmseqs input), `make_tir_combined_library` (TIR), and
  `build_fallback_tir_library.py`. The reduced library — and the RepeatMasker
  annotation — is now a reproducible function of the input *set*, not its order.
  Regression test `tests/test_canonical_sort_fasta.py`.
  - The Ty1_copia/Angela clustering and the DANTE boundary-jitter seed live in the
    upstream dependencies; both fixes have now **landed** and are pinned here:
    **dante_ltr 0.5.2.0** (canonically sorts `TE_all.fasta` before its
    `mmseqs easy-cluster`, `envs/tidecluster.yaml`) and **dante 0.2.12**
    (byte-deterministic domain output, `envs/dante.yaml`). With these plus the
    CARP-side sorts above, the LINE/Angela libraries — and the RepeatMasker
    annotation — are reproducible run-to-run. See
    `docs/dante_ltr_deterministic_clustering_request.md` and
    `docs/dante_deterministic_output_request.md`.
- **Dependency bump: DANTE_TIR 0.2.7 → 0.2.8** (`envs/dante_tir.yaml`). Verified
  annotation-neutral on the primary-TIR medium fixture: the only output change is
  the `##DANTE_TIR version` stamp in `DANTE_TIR_final.gff3`; the TIR feature set,
  the TIR library, the reduced RepeatMasker library and the final unified
  annotation are all byte-identical to 0.2.7, and repeated runs are deterministic.
- **Dependency bump: TideCluster 1.16.6 → 1.17.0** (`envs/tidecluster_run.yaml`).
  1.17.0 makes the HTML report self-contained — it vendors its per-TRC PNGs into
  `<prefix>_report/img/` instead of referencing them via `../` inside the
  `TideCluster_kite` / `TideCluster_tarean` working trees. This makes
  `cleanup_intermediates: maximal` (which purges those multi-GB trees) safe:
  verified end-to-end on the small fixture that after the maximal purge, 0 of 9
  report PNG references are missing (the report survives intact). Resolves
  `docs/tidecluster_self_contained_report_request.md`.

## 1.1.3

- **Bugfix: satellite density labels in the HTML report showed `TRC_n (?bp)`.**
  Both `make_repeat_report.R` and `make_summary_plots.R` read the monomer size
  from the old TideCluster kite file `monomer_size_best_estimate_stat.csv` /
  column `position`; current TideCluster (>=1.15) writes
  `monomer_size_top3_estimats.csv` / column `monomer_size`, so the lookup always
  missed and the `?` default stuck. Now reads the current file/column (per-TRC
  mode) and drops the `(bp)` suffix when a size is genuinely unavailable instead
  of printing `?`. Regression test `tests/test_trc_monomer_label.R`.
- **LTR_RT_TR (tandem arrays of complete LTR-RTs) surfaced in the report + GFF3.**
  - Report density panel: a new **`LTR_RT_TR`** roll-up track (L1 containers only;
    overlaps the Ty1/copia + Ty3/gypsy lineage tracks, so it is not part of the
    disjoint partition), rendered directly below the Mobile-elements aggregate.
    Backed by `Repeat_density_by_class_bigwig/LTR_RT_TR_{10k,100k}.bw`, written by
    `calculate_statistics_and_make_groups.R`.
  - Classification table: the **Complete TEs** column now shows `N (M in LTR_RT_TR)`
    per lineage — of `N` complete LTR-RTs, `M` sit inside tandem arrays (both
    counts are complete-only: DANTE_LTR ranks `DL/DLT/DLP/DLTP`, never `D`
    fragments). A darkened, asterisked legend explains it. No double count: the bp
    column counts the container span once (L1-only), members are Level-2.
  - Unified GFF3: Level-2 tandem-member copies now carry
    `in_structure=LTR_RT_TR;member_of=<container UA_L1 id>` (== their `Parent`), so
    a consumer can identify a tandem member directly. New attributes added to the
    contract (`validate_unified_gff3.py`, `docs/unified_annotation_gff3_spec.md`,
    `tests/test_unified_gff3_spec.py`). Regression test
    `tests/test_ltr_rt_tr_no_double_count.R` locks in the L1-filter + `reduce()`
    no-double-count guards.
- **Post-run cleanup of intermediate scratch (config `cleanup_intermediates`,
  default `minimal`).** After a successful run (rc == 0, non-dry-run),
  `run_pipeline.py` deletes per-tool scratch no downstream rule consumes:
  `minimal` (default) removes staged genome copies (`RepeatMasker`/`TideCluster`
  `genome_cleaned.fasta`), `DANTE_TIR.RData`, the `DANTE_LTR/library` mmseqs /
  `TE*.fasta` scratch, `filter_ltr_rt_library` side-files and the DANTE tmp GFF3
  (multiple GB on a large genome); `maximal` additionally purges the big
  `TideCluster_tarean` / `TideCluster_kite` / `TideCluster_consensus` trees and
  tool workdirs; `none` (or the CLI `--keep-all`) keeps everything. New
  `scripts/cleanup_outputs.py` derives its keep-set from `manifest.py:OUTPUTS` +
  every top-level symlink target + the CI/count-checked files + run metadata, so
  a manifest output (or the real file a deliverable symlink points at) is never
  deleted — validated by a real maximal cleanup of a fixture leaving every
  manifest output intact. The applied mode is recorded in `run_provenance.json`
  (`cleanup_intermediates`). **Tested end-to-end in CI**: the fixture jobs run
  through `run_pipeline.py`, so cleanup fires — the small fixture uses `maximal`
  and the medium `minimal`, and a new `scripts/assert_cleanup.py` step (in both
  `pipeline.yml` and the in-container `release.yml`) reads the applied mode from
  provenance and asserts the mode's scratch is gone and every manifest output
  survived. Unit tests: `tests/test_cleanup_outputs.py`.

## 1.1.2

- **DANTE_LTR 0.5.1.0 → 0.5.1.1** (`envs/tidecluster.yaml`): a newer DANTE_LTR
  build with improved efficiency on very large genomes. Dependency-identical to
  0.5.1.0 (still `r-dplyr 1.0.7.*` → r-base 4.1), so the env shape and the whole
  R/report stack are unchanged — a drop-in bump. Validated end-to-end on the
  small fixture (env rebuilt from scratch, full DAG green).
- **TideCluster 1.16.5 → 1.16.6** (`envs/tidecluster_run.yaml`): small upstream
  update. Dependency-identical to 1.16.5 (same tidehunter 1.4.3 / mmseqs2 /
  repeatmasker 4.1.2.p1 / r-igraph 2.0.3), so the env shape and the four
  `tidecluster_*` behaviours the pipeline consumes are unchanged — a drop-in bump.
- **CI: every test file is now enforced to actually run** (dev-infra). Added
  `tests/test_ci_test_coverage.py` (wired into `unit.yml`): CI goes red if any
  `tests/test_*.{py,R,sh}` is neither referenced by a workflow nor in an explicit
  exemption list. Also wired the two stdlib-only regression tests
  (`test_flank_index.py`, `test_repeatmasker_wrapper_streaming.py`) that were
  present but never executed. Closes the "test on disk but CI never runs it" gap.

## 1.1.1

- **Hotfix: `make_unified_annotation` crashed at the end of a real pipeline run
  (both CI fixtures red on the 1.1.0 tag).** The 1.1.0 streaming provenance-header
  prepend (in `finalise_output`) guarded its connection cleanup with
  `on.exit(if (isOpen(con)) close(con))`. In R, `isOpen()` **raises**
  `"invalid connection"` on an already-closed connection — it does not return
  `FALSE` — so after the explicit `close()` on the success path, the `on.exit`
  handler threw at function exit, *after* the unified GFF3 had been written
  correctly. R exited non-zero → snakemake failed the rule. The branch only runs
  when `<output_dir>/run_provenance.json` exists, i.e. every `run_pipeline.py` /
  container run (both CI fixtures, and every real user run), but **not** a bare
  `snakemake` invocation — which is why 1.1.0's local validation missed it. Fixed
  by using `try(close(con), silent = TRUE)` in `on.exit` instead of an `isOpen()`
  guard. Regression test `tests/test_provenance_prepend.R` (wired into `unit.yml`)
  exercises the success / missing / malformed-provenance paths and fails against
  the old `isOpen()` cleanup.

## 1.1.0

- **Large-genome (~90 Gbp) scaling pass across CARP's own scripts.** A sweep of
  the hangs / OOMs / quadratics that surfaced running the pipeline on a ~90 Gbp,
  >10k-contig assembly. Every change is output-identical to the prior behaviour
  and ships with a regression test:
  - `dante_line.py` + `dante_tir_fallback.py`: flank clipping rescanned every
    DANTE domain **and** every raw-TideHunter mask feature (millions on a large
    genome) per pattern/anchor — O(patterns × features). Replaced with a shared
    bisect `FeatureIndex` → O(patterns × log N) (`tests/test_flank_index.py`).
  - `repeatmasker_wrapper.py`: `split_fasta_to_chunks` / `split_fasta_to_files`
    loaded the whole genome into a Python dict (~90 GB each, twice) → stream one
    record at a time via `iter_fasta_records()` (also closes two leaked handles;
    `tests/test_repeatmasker_wrapper_streaming.py`).
  - `clean_rm_output.R` + `merge_repeat_annotations.R`: `gff_cleanup` did
    `as.list(revmap)` inside an `mclapply` that COW-forked the whole GRanges
    (296 GB max_rss on a 3.9 Gb genome) → single vectorized `extractList` +
    `unstrsplit`, no fork; `clean_rm_output.R` also reads only the 5 used
    RepeatMasker `.out` columns via `colClasses` (`tests/test_gff_cleanup.R`).
  - `make_unified_annotation.R`: dropped a per-batch full-Rle→character coercion
    (`seqnames(gr) %in% seqs`), stopped `readLines()`-ing the whole unified GFF3
    to prepend a 3-line provenance header (streaming prepend), and made the
    tier-1 overlap resolver O(k²) instead of dragging the whole batch's tier-1
    set through a greedy loop for a single overlapping pair
    (`tests/test_resolve_tier1_overlaps.R`).
  - `density_utils.R`: `rle_merge_granges` scanned all tiles once per seqlevel
    (O(tiles × seqnames), effectively a hang at 90M tiles × >10k contigs) → a
    single pass over the sorted seqnames-Rle runs (`tests/test_density_utils.R`).
  - `reduce_dimer_library.py`: parallelised the per-TRC reduction (one group per
    worker, each mmseqs single-threaded; the reduction is thread-invariant),
    with per-group temp cleanup (`tests/test_reduce_dimer_parallel.sh`).
  - `calculate_statistics_and_make_groups.R`: `readDNAStringSet(genome)`
    (~90 GB just to sum widths) → `fasta.seqlengths()`.
  - `Snakefile`: added the `ulimit -n` open-file stopgap to `tidecluster_long`,
    `tidecluster_short` and `repeatmasker` (mirrors the existing `dante_ltr` /
    `tidecluster_reannotate` stopgaps). Full audit in
    `docs/carp_large_genome_audit.md`.
- **Dependency modernization to the released upstream tools that now ship CARP's
  large-genome fixes: `dante` 0.2.11, `dante_ltr` 0.5.1.0, `dante_tir` 0.2.7,
  `TideCluster` 1.16.5.** These releases contain the fd-budget / merge-tail /
  streaming-flanks / bounded-handle fixes contributed upstream from this project.
  `dante` 0.2.11 moved to `r-base 4.2.3` while `dante_ltr` 0.5.1.0 still pins
  `r-dplyr 1.0.7` (`r-base <4.2`), so DANTE can no longer share
  `envs/tidecluster.yaml`: it splits into its own **`envs/dante.yaml`**
  (r-base 4.2.3), and the two rules that run DANTE-package tools (`dante`,
  `filter_dante` → `dante_gff_output_filtering.py`) point at it. Validated
  end-to-end on tiny_pea (29/29 jobs, 0 errors; dante_ltr's 0.4→0.5 GFF3 schema
  intact; `validate_classifications` green).
- **Fixed misplaced plot / y-axis labels in the benchmark report's 2nd–4th
  charts.** The four inline SVGs (from cairo's `svg()`) carried document-local
  ids (`clip-N`, `glyph-N-M`) that collide when inlined into one page, so later
  charts resolved their glyph/clip references to the first chart's definitions.
  `make_benchmark_report.R` now namespaces each SVG's ids + references with a
  unique per-chart prefix.

- **Pin CARP's own conda env dependencies to the validated 1.0.7 SIF versions
  (reproducibility).** The env files listed several deps unpinned
  (`envs/dante_line.yaml`: `seqkit`, `parasail-python`, `mmseqs2`, `blast`,
  `pyyaml`; `envs/tidecluster.yaml` / `envs/tidecluster_run.yaml`: the
  `r-jsonlite`/`r-yaml`/`r-optparse`/`pyyaml`/`jq` helper deps and
  `bioconductor-rtracklayer`/`mmseqs2`), so each build re-solved them to whatever
  the channels served that day. They are now pinned to the exact versions that
  the validated 1.0.7 container shipped (a known co-solvable set), so the CARP
  envs build reproducibly. The r-* helper pins deliberately differ between the
  two TideCluster envs because they track each env's `r-base` (4.1.3 vs 4.3.3).
  The upstream tools' transitive deps are still unconstrained — that is the
  subject of `docs/dante_tidecluster_dependency_pinning_request.md`.

## 1.0.7
- **Pin the Singularity base image and harden `%post` against base-image drift
  (fixes the SIF build).** The `continuumio/miniconda3` base was unpinned, so
  every build pulled `:latest`; between the working 1.0.4 build and the next one
  it drifted twice and broke the release build (the `%post` and `envs/` were
  byte-identical the whole time — `git diff 1.0.4 1.0.5 -- Singularity envs/` is
  empty). First a newer conda began enforcing Anaconda's channel Terms of
  Service on the default `pkgs/main`/`pkgs/r` channels
  (`CondaToSNonInteractiveError`) — and that conda build enforced it but didn't
  ship the `conda tos` CLI to accept it. Then the base python jumped to 3.14 with
  an Anaconda helper package (`anaconda-channel-guide`) that hard-pinned it,
  blocking the downgrade to 3.11 (`LibMambaUnsatisfiableError`). Fixes:
  (1) **pin** `From: continuumio/miniconda3:24.9.2-0` (a pre-drift conda 24.9.2
  build — the same conda the project's sandbox uses); (2) build only from
  **conda-forge + bioconda**, never Anaconda's default channels (drops `defaults`
  from the config + `--override-channels` on the bootstrap installs) — avoids the
  ToS gate on any conda and keeps CARP off Anaconda's commercially-licensed
  channels; (3) stop downgrading the base env's python — snakemake 8.12.0 is a
  noarch package, so it installs onto the base python (3.12) as-is; the old
  `conda install python=3.11` rewrote the base env and broke the base `conda`
  package, whose next operation then failed plugin discovery
  (`PluginError: Conflicting post_solves plugins: signature-verification`);
  (4) remove the Anaconda telemetry helper packages (`anaconda-anon-usage`,
  `anaconda-channel-guide`) so `conda info --json` stays clean for Snakemake.
  (The earlier 1.0.4 `jq` system-install fix for the Bioconductor
  `genomeinfodbdata` post-link is retained.)
- **HTML report: fix large-genome crash and make report generation non-fatal.**
  `make_repeat_report.R` aborted on a large genome while building the density
  panels: a bin midpoint on a chromosome > ~1.07 Gbp overflowed 32-bit integer
  arithmetic (`starts + ends` → `NA`), and those `NA`s drove the recursive
  `rdp_simplify` line-simplifier into O(n) recursion → `Error: C stack usage is
  too close to the limit`. Fixes: (1) do the midpoint sum in double precision so
  it never overflows; (2) rewrite `rdp_simplify` with an explicit stack (no
  recursion → cannot overflow the C stack), verified identical to the old output
  on 300 random cases. Robustness: (3) each chart section is now built via
  `safe_build()` so a failure in one panel leaves a visible "not generated"
  placeholder and the rest of the report still renders; (4) the `make_repeat_report`
  rule writes a minimal placeholder HTML if the script fails entirely, so a
  report failure can never abort the pipeline (the annotation outputs are
  unaffected).

## 1.0.4
- **Fix Singularity/SIF build failure (bioconda `genomeinfodbdata` post-link),
  properly this time.** `bioconductor-genomeinfodbdata`'s post-link script runs
  `installBiocDataPackage.sh` → `yq` → `jq`; with `jq` absent the link errored
  and the whole SIF build failed (`yq: Error starting jq`). 1.0.3 tried pinning
  `jq` inside the two Bioconductor env yamls, but that did **not** fix it: conda
  gives no ordering guarantee that the env-level `jq` links before
  `genomeinfodbdata`'s post-link runs (no dependency between them), so `jq` can
  still be missing at post-link time. The authoritative fix is to install `jq`
  at the **system level** in the `Singularity` `%post` (`apt-get install jq`,
  before the conda env creation), so it is on `PATH` for every post-link
  regardless of link order. The env-yaml `jq` is kept as redundancy. (1.0.2 and
  1.0.3 published no artifacts — their builds never completed — so this is the
  first release to ship the flank-alignment grouping change below.)
- **Bounded, deterministic grouping for the DANTE_TIR_FALLBACK / DANTE_LINE flank
  alignment (OOM fix on large genomes).** The all-vs-all flank aligner
  (`scripts/global_local_aln.py`) is O(N²) in memory and compute; on a big
  genome a single high-copy TIR subtype (e.g. CACTA/Mutator) has thousands of
  near-identical TPase anchors, so even after the mmseqs prefilter the candidate
  pairs stay ~O(N²) and the run OOMs. It now mirrors DANTE_TIR's `--max_class_size`
  strategy: when a subtype (fallback) or the LINE pattern set exceeds
  `dante_tir_fallback_max_group_size` / `dante_line_max_group_size` (default
  `1000`), the sequences are split into deterministic, mmseqs-clustering-based
  groups of at most that size and aligned per group, dropping cost to O(N·group).
  Clustering keeps similar flanks together, so the split is near-lossless
  (verified byte-identical `Selected_Length` vs the ungrouped path on distinct
  families); at or below the threshold the code path is unchanged. New unit test
  `tests/test_aln_grouping.py`.

## 1.0.1
- **TideCluster 1.16.1 → 1.16.3.** 1.16.2 fixed the TRC-superfamily
  empty-fallback naming (canonical `<prefix>_trc_superfamilies.csv` with header,
  even when empty), giving the map consumed by
  `tidecluster_reannotate_superfamily_merge` a stable name/schema. 1.16.3 fixes
  the `tc_reannotate` chunk-pool RepeatMasker library-build race that could
  silently truncate the reannotation on a fresh env's first run (RM libraries
  are now built once before the pool, and a failed chunk fails the run instead
  of being silently dropped).
- **New `tidecluster_reannotate_superfamily_merge` option (default `True`).**
  Superfamily-aware array recovery for the RM-on-TideCluster reannotation.
  `tc_reannotate`'s array-length filter is strictly per-TRC, so a real tandem
  array tiled by ≥2 TRCs of the **same superfamily** is fragmented below
  threshold and lost — the underlying TE then prevails, and enabling
  `tidecluster_reannotate_culling_limit` only masks this by collapsing each locus
  to one TRC (making the annotation culling-dependent). The pipeline now runs
  `tc_reannotate --debug` and re-filters its raw hits
  (`scripts/tc_reannotate_sf_filter.py`) grouping sibling TRCs by superfamily, so
  those arrays are recovered **deterministically, independent of culling**, while
  every feature keeps its bare `TRC_<n>` `Name`.
- **New `rm_tc_tandem_gate` option (default `True`).** A Tier-4 RM-on-TideCluster
  satellite may override a Tier-5 TE call only where it has independent tandem
  evidence (raw TideHunter); an unsupported RM_TC array — a short AT-rich
  consensus tiling a genuinely non-tandem TE (e.g. a Tekay LTR-RT) — is demoted
  below the TE, preventing spurious TE→satellite over-masking. Genuine satellites
  and RM_TC over non-TE sequence are unaffected. Measured on a Merodon genome:
  recovers ~7 Mb of TE that culling had re-labelled satellite, while keeping the
  ~22 Mb of real TE-derived tandem arrays.

## 1.0.0
- **New `repeatmasker_culling_limit` option (default `0` = off).** Passes rmblastn
  `-culling_limit N` to the RepeatMasker step to cap the redundant per-locus HSP
  explosion a de-novo library causes (each genomic TE copy aligns to hundreds of
  near-identical consensi). `2` → ~3× faster RepeatMasker at ~−0.7 % masked bp,
  classification preserved. Injected via an `RMBLAST_DIR` shim
  (`scripts/rmblast_culling_shim.py`) — no RepeatMasker/container patch, works in
  the existing image.
- **New `tidecluster_reannotate_culling_limit` option (default `0` = off).** The
  same culling for TideCluster's internal RepeatMasker in `tidecluster_reannotate`
  — the strongest culling target (a consensus dimer matches a satellite array at
  every rotational phase). `2` → ~3.7× faster reannotation at ~−0.8 % masked bp.
  Independent of `repeatmasker_culling_limit`; validate masked bp on
  large-satellite genomes before enabling (may under-mask megasatellite arrays).
- **`dante_tir_fallback` performance fix.** Reverted a subtype-parallelisation
  that regressed the rule ~3.6× on genomes with size-imbalanced TIR subtypes (the
  dominant subtype was starved to `threads/N` while small subtypes left cores
  idle). Subtypes now run serially with the full thread budget; output unchanged.
- **Tool bumps:** TideCluster 1.16.0 → 1.16.1, dante_tir 0.2.5 → 0.2.6.
- **rDNA classification restructured (breaking; output schema v3).** The flat
  top-level classes `rDNA_45S` / `rDNA_5S` are now nested under a single `rDNA`
  parent — `rDNA/45S_rDNA`, `rDNA/45S_rDNA/{18S,25S,5.8S,IGS,ITS1,ITS2}`,
  `rDNA/5S_rDNA`, `rDNA/5S_rDNA/5S` (`classification_vocabulary.yaml` +
  `data/rdna_library.fasta` headers + the TideCluster rDNA map in
  `make_unified_annotation.R`). This **renames the per-class rDNA GFF3/BigWig
  files** (`rDNA_45S.* → rDNA.45S_rDNA.*`, etc.); the `rDNA` **rollups**
  (`rDNA_RepeatMasker.gff3`, `rDNA_{10k,100k}.bw`) and the REQUIRED outputs keep
  their names. No backward-compat symlinks; the manifest `schema_version` bump to
  `3` is the migration signal. Consumers that resolve by manifest logical name
  are unaffected; those globbing `rDNA_45S*` filenames must update. See
  [`docs/output_schema.md`](docs/output_schema.md).
- **HTML report:** classification table and pie now use a **fixed category
  order** (Class_I, Class_II, rDNA, Tandem_repeats, Simple_repeat,
  Low_complexity, Unknown) instead of ordering by content; the pie is
  **genome-relative** (inner ring = Repeats vs grey Non-repetitive, so every
  percentage matches the table); overview cards show **LTR-RT / DNA transposons /
  LINE as % of genome** and the "Sequences in charts" card was removed. The
  table+pie show **rDNA at the array level** (45S / 5S; subunit detail kept in
  the GFF3/density outputs but not broken out). The "Unspecified" row tag was
  renamed **"Unclassified"**. The per-sequence repeat-content panel now draws
  **two reference lines** (shown-sequence average + whole-genome average) with a
  note on how many contigs/bp are omitted from the bars. Fixed an off-by-one that
  dropped the last chromosome separator line in `summary_plots.pdf`.
- **Per-sequence repeat-content panel reworked + Level-1-only statistics.** The
  per-sequence composition bars are now computed at base precision from the
  unified GFF3 Level-1 union (not the smoothed density BigWigs, which inflated
  the totals), so each bar sums to the sequence's true repeat fraction. All
  sub-threshold contigs are aggregated into one prominent "Other (M contigs,
  X Mb)" bar, so the bars cover the whole genome; the redundant "shown avg" line
  is dropped (one "genome avg" line remains). `summary_statistics.csv` and the
  per-class split now count **Level-1 features only** — Level-2 nested children
  (LTR_RT_TR member copies, simple repeats nested in satellites) are no longer
  double-counted, so per-class numbers and the total equal the true union and
  the per-class partition is exactly disjoint. The Class_II bar/legend label is
  now just "Class II" (it includes Helitron, not only TIR).
- **Fixed tier-4 rDNA mislabelling in the unified annotation:** rDNA arrays
  annotated only via the RM-on-TideCluster pass were labelled
  `Satellite/TideCluster/TRC_<n>` instead of rDNA; the authoritative TideCluster
  `<prefix>_rdna.tsv` now drives rDNA labelling across all TideCluster tiers.

## 1.0.0rc3
- Density BigWig tracks no longer exceed 1.0: overlapping features are merged into a **strand-agnostic union** before coverage (`calculate_density.R` / `calculate_density_batch.R`). The Unified annotation tolerates overlap (L1 `Simple_repeat`/`Low_complexity` over a TE; nested L2 children), which previously stacked the total track to ~3.5× and per-class tracks to ~2.2×; every density track is now a true union fraction in `[0, 1]`. Validated across three full assemblies (Boechera, Dunaliella, *A. thaliana* Col-CC T2T).
- **Tandem LTR-RT (`LTR_RT_TR`)**: head-to-tail, same-lineage LTR-RT arrays that share boundary LTRs are collapsed to one Level-1 container with the member copies as Level-2 children (`scripts/resolve_ltr_tandems.py`), so a shared-LTR array is annotated once instead of double-counting the overlapping LTRs. Containers are also written to a new top-level `DANTE_LTR_tandems.gff3` (header-only when none; logical name `dante_ltr_tandems_gff3`). See [`docs/dante_ltr_tandem_feature_request.md`](dante_ltr_tandem_feature_request.md).
- **TE-derived satellite conflict** resolved: where a TideCluster satellite is a tandem of complete LTR-RTs, the satellite wins the region and is tagged `TE_origin` (plus `TE_origin_structure=tandem_LTR_RT` for the full-LTR-RT case), and the underlying structural elements/members are trimmed out of the unified file — eliminating the double annotation. The unified-GFF3 drift guard and spec were extended to cover the new attributes.

## 1.0.0rc2
- Output contract for `Repeat_Annotation_Unified.gff3` written down ([`docs/unified_annotation_gff3_spec.md`](unified_annotation_gff3_spec.md)) and enforced by an executable drift guard (`scripts/validate_unified_gff3.py` + `tests/test_unified_gff3_spec.py`, run in CI and the release gate).
- Fixed a satellite `Name` regression that silently emptied the per-family BigWig outputs — every TideCluster satellite `Name` stays the bare `TRC_<n>` (downstream apps and `split_gff_by_name.R --name-prefix TRC_` key on it); rDNA is routed by `classification` instead.
- New `/release` skill: one-command version bump + cheap-CI gate + tag.

## 1.0.0rc1
- **TideCluster upgraded to 1.16.0**, adding three default-on behaviours the unified annotation consumes: array-level **rDNA identification** (`rDNA_45S`/`rDNA_5S`, labelled clearly as rDNA while still counted as a tandem family), cross-TRC overlap resolution, and a chunked/pooled `tc_reannotate` RepeatMasker (fixes `-pa`-with-custom-`-lib` under-parallelism). New config knobs `tidecluster_detect_rdna` / `tidecluster_rdna_library` / `tidecluster_keep_trc_overlaps` / `tidecluster_chunk_size`.
- Unified annotation now **resolves TR-from-TE overlaps**: a tandem array built from multiple same-family structural TEs is reported once as the satellite with a `TE_origin` tag rather than double-annotated; a non-fatal `Repeat_Annotation_Unified.overlaps.tsv` reports any residual Level-1 overlaps.
- Major performance/memory work: fixed a `reduce_library` merge OOM (~296 GB RSS → bounded), parallelised the density batch / DANTE_TIR fallback / `reduce_library` BLAST+CAP3, and packed small scaffolds into shared RepeatMasker chunks.

## 0.9.0rc10
- New **second-round containment reduction** of the RepeatMasker library (`reduce_library_containment`, default True, independent of `reduce_library`). After the per-class CAP3/mmseqs reduction, a `reduce_library_containment` rule drops short repeat fragments fully contained in a longer element of the **same class** (greedy blastn self-comparison, `scripts/containment_reduce_library.py`); RepeatMasker masks their genomic copies via the container, so masking **and** classification are preserved. Validated masked-bp-lossless with RepeatMasker on the Pisum pangenome over two genome regions (−0.09…−0.12%) while cutting ~22% of library bp and ~30% of RepeatMasker wall-time. Thresholds via `containment_min_identity` (80) / `containment_min_coverage` (0.90); set `reduce_library_containment: False` to feed RepeatMasker the full per-class-reduced library.

## 0.9.0rc9
- Fixed a hard pipeline failure in `tidecluster_reannotate`: the TideCluster dimer-library reducer (`reduce_dimer_library.py`) ran `mmseqs easy-search` with the default nucleotide k-mer (15), which **segfaults on short tandem monomers** (the query is the first-half monomer, down to ~21 bp). Pinned `-k 7`; groups with monomers below the prefilter floor are skipped and any residual `mmseqs` error retries single-threaded then keeps the group unreduced — so the reducer can no longer abort a run.
- Dimer reduction made lossless on more families: replaced single-linkage clustering with greedy/star clustering (each dropped dimer aligns directly to a kept representative) plus a length guard (nested/harmonic periods no longer collapse onto the longest rep). Validated masked-bp-lossless on tiny_pea (−0.02%) and GCA_041296365.1. Known limitation: large satellite families (dimer ≳ a few kb) can still under-mask because RepeatMasker needs the full consensus diversity — set `reduce_tidecluster_library: False` to mask with the full dimer library.

## 0.9.0rc8
- Per-family / by-class density BigWigs build ~140× faster on assemblies with many scaffolds: `calculate_density_batch.R` now parallelises across families (`-t`, wired to `workflow.cores`) and computes density only on the scaffolds each family occupies instead of re-binning the whole genome per family. Track values are byte-identical to rc7; empty scaffolds are no longer written to the BigWig header (occupied-only). On a 2 Gbp / 1888-scaffold assembly with 477 tandem families the step drops from ~6.5 h to ~3 min.
- CI now asserts `carp_manifest.json` matches the produced output tree — every declared logical-name→path must exist after a run (`scripts/assert_manifest_outputs.py`, run on the release-gate fixture) — so the manifest can no longer drift out of sync with the pipeline's outputs.

## 0.9.0rc7
- Emit `carp_manifest.json` at the output root (machine-readable output contract: `schema_version` + logical-name→path map; written by `run_pipeline.py` on success and failure). See [`docs/output_schema.md`](output_schema.md).
- Faster `tidecluster_reannotate`: rotation-invariant reduction of the TideCluster consensus dimer library (aligns each monomer against the already-doubled dimers, collapsing redundant phase variants per TRC). Decoupled from `reduce_library` via the new `reduce_tidecluster_library` option (default True). Lossless for masking (validated; masked bp unchanged within ±0.15%), library ~5–20× smaller, remasking ~1.5–5× faster.
- `summary_plots.pdf` hardened: empty repeat categories render a placeholder panel instead of crashing; any render failure (e.g. data too large) falls back to a one-page placeholder PDF so the rule never fails.
- Snakefile-derived workflow schematic ([`figs/workflow_overview.svg`](../figs/workflow_overview.svg), regenerated by `scripts/make_workflow_diagram.py`).

## 0.9.0rc6
- Density BigWigs reworked: all annotation-derived tracks now sourced from `Repeat_Annotation_Unified.gff3`; BigWig outputs renamed/relocated (breaking — see "Migration" in [`docs/output_reference.md`](output_reference.md)); genome-wide total now includes tandems (`Repeat_density/Repeat_density_total_*.bw`); per-family tandem tracks come in two flavours per `TRC_<n>` (structural TideCluster + Unified union); fixed crash on empty `RM_on_TideCluster_Library.gff3`.

## 0.9.0rc5
- Density BigWigs written run-length-merged (sparse, lossless, much smaller); new `Class_II.Subclass_1.TIR` density rollup; new per-family BigWigs for the RepeatMasker tandem pass (superseded by rc6).

## 0.8.0
- Bug fix in DANTE_LINE, filtering of tandem repeats from DANTE_LINE added.

## 0.7.4
- TideCluster updated to v1.8.0 with `--long` option added to detect tandem repeats with monomer up to 25 kb. Bugfix in subtracting tandem repeats from dispersed repeats.

## 0.7.2
- DANTE_LINE added.

## 0.7.1
- DANTE_TIR added.

## 0.6.7
- More efficient calculation of bigwig files.

## 0.6.6
- DANTE_LTR updated to 0.6.0.4, TideCluster updated to 1.6.

## 0.6.5
- Bugfix in gff3 merging.

## 0.6.4
- dante_ltr runs on smaller chunks (50000000) → better memory usage.

## 0.6.3
- DANTE update to 0.2.5 — bugfix.

## 0.6.2
- Bugfix in bigwig calculation.

## 0.6.1
- DANTE_LTR update to 0.4.0.3 (bugfix).

## 0.6.0
- REXdb Viridiplantae v4.0, library size reduction added, RepeatMasker parallelization added, missing full LTR-RT handling added.

## 0.5.2
- RepeatMasker sensitivity can be set.

## 0.5.1
- Graphical output to PDF added.

# Large-genome (~90 Gbp) scaling analysis — plan + findings

**Status:** analysis complete (Phases 0–4 done). Upstream feature requests
written; pipeline-side stopgaps applied and DAG-validated. See
**Findings & outcomes** at the bottom.
**Trigger:** running CARP on a ~90 Gbp assembly (N50 ≈ 2.1 Gb, >10k contigs,
`-c 32`, `-S 50000000`) surfaces two hard failures and raises the question of
whether sibling modules share the same weaknesses.

## 0. Problem statement

Two confirmed problems on the 90 Gbp run, plus an audit request:

1. **DANTE — single-threaded post-search stage.** The `lastal` similarity
   search parallelises across `-c {threads}`, but the stage *after* the search
   (hit dedup → top-hit / LCA classification → domain filtering → GFF3
   assembly) runs in **one process** and takes roughly as long as the search
   itself. On a 90 Gbp genome that serial tail dominates wall-time and leaves
   the machine ~idle. **Goal:** parallelise (or stream) the parse/assemble
   stage so CPU utilisation stays high after the search.

2. **DANTE_LTR — "Too many open files".** On the big genome DANTE_LTR aborts
   with a file-descriptor-exhaustion error. **Goal:** find the fd source
   (leak vs. genuine concurrency need vs. too-low `ulimit -n`) and fix it so
   the tool completes on large / fragmented assemblies.

3. **Audit the siblings for the same two failure classes** (serial-tail CPU
   waste + fd exhaustion), and any large-genome-specific memory/temp-file
   blowups:
   - upstream: **DANTE_TIR**, **TideCluster**
   - in-repo (ours, directly fixable): **`scripts/dante_line.py`**,
     **`scripts/dante_tir_fallback.py`**

## Deliverables

- One upstream **feature-request markdown per tool that needs upstream
  changes**, following the existing convention
  (`docs/tidecluster_reannotate_chunk_race_request.md`,
  `docs/tidecluster_superfamily_naming_request.md`,
  `docs/dante_ltr_tandem_feature_request.md`): a `docs/<tool>_<topic>_request.md`
  with Summary → Evidence (numbers/traceback) → Root cause → Proposed fix →
  Pipeline-side stopgap. Candidates: `dante`, `dante_ltr`, `dante_tir`,
  `TideCluster`.
- **Direct in-repo fixes** (PR branch) for `dante_line.py` /
  `dante_tir_fallback.py` where the problem is ours to fix.
- A **pipeline-side stopgap section** (things CARP can do *now* in the
  Snakefile — raise `ulimit -n` in a rule, tune `-S`/`-wd`/`--chunk_size`,
  size `TMPDIR`) so the 90 Gbp job can keep running before upstream releases
  land.

## Versions under analysis (must analyse the exact shipped versions)

| tool | version in pipeline | source | how invoked |
|------|--------------------|--------|-------------|
| DANTE | 0.2.5 | `github.com/kavonrtep/dante` | `dante -q IN -o OUT -c {threads}` (`Snakefile:421`) |
| DANTE_LTR | 0.4.0.4 | `github.com/kavonrtep/dante_ltr` | `dante_ltr -o P -s FA -g GFF -c {threads} -M 1 -S 50000000` (`Snakefile:793`) |
| DANTE_TIR | 0.2.6 | `github.com/kavonrtep/dante_tir` | `dante_tir.py -g GFF -f FA -o DIR -c {threads}` (`Snakefile:450`) + `dante_tir_summary.R` |
| TideCluster | 1.16.3 | `github.com/kavonrtep/TideCluster` | `run_all` / `tc_reannotate.py` (rules `tidecluster_long/short/reannotate`) |
| dante_line (ours) | in-repo | `scripts/dante_line.py` | `dante_line.py … -t {threads}` (`Snakefile:758`) |
| dante_tir_fallback (ours) | in-repo | `scripts/dante_tir_fallback.py` | `dante_tir_fallback.py … -t {threads}` (`Snakefile:524`) |

## Method — the same 5-step protocol per tool

For each of the six targets, run the identical protocol so findings are
comparable and each yields a self-contained request/fix.

**S1 — Pin & clone.** `git clone` the repo and `git checkout` the tag that
matches the conda version above (fall back to the release commit if the tag is
missing). Clone into scratch, not the project tree. Analyse *that* revision —
not `master` — because the fix must apply to what the SIF actually ships.

**S2 — Map stages, mark serial vs. parallel, enumerate fd sites.** Read the
entry point and trace the pipeline. For every target, produce two lists:
  - *Compute stages*, each tagged parallel (worker pool / mclapply / mmseqs
    `--threads` / lastal `-c`) or **serial** — and for each serial stage, whether
    it is embarrassingly parallel by chunk / contig / lineage / superfamily.
  - *File-descriptor sites*: `open()` / `subprocess`/`Popen` / `tempfile.mktemp`
    / R `file()`/`pipe()`/`gzfile()`/`textConnection` / blast+`makeblastdb`
    temp DBs / per-chunk temp files — and whether the count scales with genome
    size, chunk count, or contig count, and whether handles are held
    concurrently (the fd-exhaustion condition) vs. opened-closed serially.

**S3 — Characterise scaling for 90 Gbp** (order-of-magnitude, to confirm):
  - DANTE windows at `-wd` 10 Mb default → ~9,000 windows.
  - DANTE_LTR chunks at `-S` 50 Mb → ~1,800 chunks.
  - TideCluster at `--chunk_size` 50 Mb → ~1,800 chunks.
  - Contig count is unknown and matters for fd (a fragmented 90 Gbp assembly
    can have >10⁵ scaffolds). **Ask the user** for `grep -c '>' genome.fa`
    and the assembly N50 — this changes which hypothesis is dominant.

**S4 — Reproduce at reduced scale + profile (confirm, don't guess).** We can't
run 90 Gbp here, so build the *smallest input that exercises the same code
path*:
  - fd exhaustion → a genome of *many small contigs* (e.g. 50k × 2 kb) run
    under a deliberately low `ulimit -n 256`, watching peak concurrent fds via
    `ls /proc/<pid>/fd | wc -l` sampling or `strace -f -e trace=openat,close`.
  - serial-tail CPU waste → a repeat-dense chunk timed with `/usr/bin/time -v`;
    locate the hot serial loop with `py-spy record`/`py-spy dump` (Python) or
    `Rprof()` (R). Measure the serial stage's share of wall-time to size the
    win.
  - **Authoritative evidence comes from the user's real run** — request the
    DANTE_LTR traceback (the exact failing call), the `benchmarks/dante.tsv` /
    `benchmarks/dante_ltr.tsv` rows, and the log tails. Reduced-scale repro
    *validates the mechanism*; the real logs *prove it's the one that bit*.

**S5 — Design the fix + stopgap.** Write the concrete upstream change (which
function, which stage, what strategy) and a throwaway pipeline-side mitigation.

## Preliminary hypotheses (from reconnaissance already done)

Starting points to confirm/refute in S2–S4 — not conclusions.

- **DANTE 0.2.5.** `dante` parallelises `lastal` per window, then a single
  process runs hit dedup + top-hit(80%)/LCA classification (Fig-1 B/C of Novák
  et al. 2024) and `dante_gff_output_filtering.py`. Hypothesis: that tail is
  embarrassingly parallel **per window/contig** (classification is local to a
  hit-region) and the fix is to fan it across the existing pool then concat, or
  to pipeline it so parsing of window *i* overlaps the search of window *i+1*.
  Confirm where the serial boundary is (whole-genome sort/merge? a global
  structure?) — that determines whether per-chunk fan-out is valid.

- **DANTE_LTR 0.4.0.4.** R-based, chunked by `-S`. fd-exhaustion candidates:
  (a) per-chunk BLAST DBs / temp files opened inside an `mclapply`/`foreach`
  pool and not released before the next chunk; (b) an R connection
  (`file()`/`pipe()`/`gzfile()`) opened in a loop without `close()`; (c) simply
  the default soft `ulimit -n` (often 1024) being below `chunks × per-chunk
  handles`. Grep the R source for connection-opens inside loops and the
  parallel backend; count worst-case concurrent handles. The right fix differs
  per cause (close-in-loop vs. bounded pool vs. `ulimit` bump + doc).

- **DANTE_TIR 0.2.6.** Python + `dante_tir_summary.R` (mmseqs-backed). Check
  the per-superfamily worker layout for the same **thread-starvation** pattern
  the in-repo fallback already documents (see below), and fd from per-element
  `seqkit`/`blast` subprocesses across many TPase domains.

- **TideCluster 1.16.3.** Already chunked/pooled (we drove that via the two
  existing requests). Audit `tc_reannotate` + clustering fd behaviour across
  ~1,800 chunks and a high contig count; check TideHunter/mmseqs thread use and
  temp-file cleanup at this scale.

- **`scripts/dante_line.py` + `scripts/dante_tir_fallback.py` (ours).** Mostly
  **serial** Python glue around `seqkit`/`mmseqs` subprocesses, using
  `tempfile.mktemp` heavily (`dante_line.py:402,453,563,648-669`;
  `dante_tir_fallback.py:311,342-368,530`). Two things to nail down:
  1. `run_prime_alignments(..., threads=)` (`dante_line.py:1145`) — does the
     parasail alignment step actually parallelise, or is it a serial
     per-pattern loop that should use a pool? (Same failure class as DANTE's
     serial tail, but ours to fix.)
  2. `dante_tir_fallback.py:1088-1091` already carries a comment about
     `threads_per_subtype = threads // n_workers` **starving the one big
     subtype** — verify and fix the thread-allocation so the largest subtype
     isn't throttled.
  3. `mktemp` usage: files land in `TMPDIR` and (a) are a known race/predictable-
     name smell and (b) at 90 Gbp scale × thousands of groups can accumulate /
     exhaust temp space or fds if not unlinked. Audit cleanup; consider
     `NamedTemporaryFile`/`mkstemp`+context managers.

## Cross-cutting: immediate pipeline-side stopgaps (so the 90 Gbp job runs now)

Independent of upstream releases, evaluate and (where safe) apply:
- `ulimit -n <high>` prepended in the `dante_ltr` rule shell (and any other rule
  that trips it) — unblocks fd exhaustion immediately if the cause is (c).
- Chunk-size tuning: larger `-S` / `--chunk_size` = fewer chunks = fewer
  concurrent handles (trades memory); smaller `-wd` changes DANTE's window
  count. Quantify the trade.
- `TMPDIR` on a large, fast filesystem sized for the temp-file peak.

## Sequencing

- **Phase 0** — this document. *(get sign-off)*
- **Phase 1** — S1–S2 for all six targets (clone + static map). Cheap, parallelisable.
- **Phase 2** — S3–S4 on the two **blocking** problems first (DANTE serial
  tail, DANTE_LTR fd), using reduced-scale repro + the user's real logs.
- **Phase 3** — write the upstream feature requests + the in-repo fixes.
- **Phase 4** — land the pipeline-side stopgaps in the Snakefile.

## Inputs needed from the user

1. Confirm I may `git clone` the four `kavonrtep/*` repos here (network), or
   drop local clones into scratch.
2. The **DANTE_LTR traceback** ("too many open files" — the exact failing
   line/call) and, if kept, `benchmarks/dante_ltr.tsv`.
3. The **DANTE benchmark** row (`benchmarks/dante.tsv`) or search-vs-parse
   wall-time split, to size the serial-tail win.
4. Assembly shape: `grep -c '>' genome.fa` (contig count) and N50.
5. Preference: upstream findings as `docs/*_request.md` (default, matches
   convention) or filed as GitHub issues.

---

## Findings & outcomes (2026-07-14)

Sources analysed at the exact shipped tags (cloned under `scratch/`, gitignored):
dante 0.2.5, dante_ltr 0.4.0.4, dante_tir 0.2.6, TideCluster 1.16.3, plus the
in-repo `scripts/dante_line.py` / `scripts/dante_tir_fallback.py`.

### Upstream feature requests written

| tool | doc | headline finding | class |
|------|-----|------------------|-------|
| DANTE_LTR | `dante_ltr_too_many_open_files_request.md` | `[open(f,'w') for f in temp_files]` at `dante_ltr:565/585` opens ~1800 handles at once (90 Gbp / 50 Mb) — the reported crash. `dante` already fixed this via open-append. | **fd — confirmed crash** |
| DANTE | `dante_parse_parallelization_request.md` | single-threaded post-`Pool` merge tail (recalc → global `sort` → dedup, wrapper 576-607); amplified by `chunk_size=500000` → ~180k chunks; `lastal` run twice/chunk | **serial CPU** |
| TideCluster | `tidecluster_large_genome_request.md` | **same fd bug** at `tc_utils.py:1040` on the `tc_reannotate` path (the pipeline uses it) **+** clustering loads whole genome (~90 GB) that is then unused (`TideCluster.py:623`) | **fd — confirmed crash + OOM** |
| DANTE_TIR | `dante_tir_large_genome_request.md` | whole genome into a Python dict ~90–180 GB (`dt_utils.py:168`, OOM) + single-threaded Round-1 `beast()` loop (`dt_utils.R:900/1170`) | **memory + serial CPU** |

### In-repo scripts — audited CLEAN (no code change)

`dante_line.py` + `dante_tir_fallback.py` have **none** of the three failure
classes: sequence extraction streams via `seqkit subseq`/`faidx` (no
whole-genome dict), no concurrent-handle fan-out, the parasail all-vs-all is
already parallel (`global_local_aln.py:_compare_sequences`,
`ThreadPoolExecutor` + `max_group_size` O(N²) guard), and temp files are cleaned
in `finally`. The subtype-serial strategy in `dante_tir_fallback.py:1083-1098`
is a deliberate, benchmarked choice (concurrent variant was 3.6× slower). Only
nit: `tempfile.mktemp` is deprecated (hygiene, not a scaling bug) — left as-is.

### Pipeline-side stopgaps applied (Snakefile)

`ulimit -n "$(ulimit -Hn)"` added to the **`dante_ltr`** and
**`tidecluster_reannotate`** rule shells — both open ~1800 chunk handles at 90
Gbp and both hit the default 1024 soft limit. DAG re-validated (41 jobs build
clean). This unblocks the run only where the hard limit exceeds the chunk count;
the durable fix is upstream (open-append / bounded-handle). **No stopgap exists**
for the DANTE_TIR / TideCluster-clustering whole-genome memory loads — those need
the upstream fix.

### Not yet done (needs the user)

- File the four requests upstream (kavonrtep/*) — or hand them to the maintainer.
- The two OOM findings (DANTE_TIR dict, TideCluster clustering) will still fail a
  90 Gbp run even with the fd stopgaps; they are the next hard blockers.

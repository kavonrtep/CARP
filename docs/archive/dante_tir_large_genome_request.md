# DANTE_TIR: large-genome memory and single-threaded Round-1 bottleneck

> Feature request for `github.com/kavonrtep/dante_tir` — audited on **0.2.6**.
> Invoked as `dante_tir.py -g DANTE.gff3 -f genome.fasta -o OUTDIR -c 32`.

## Summary

DANTE_TIR does **not** have the `dante_ltr`-style "too many open files" bug —
every `open()` is context-managed and no handle count scales with genome or
contig count. But on a ~90 Gbp / >10k-contig assembly it has two real hazards:

1. **Memory (most serious): the whole genome is loaded into a Python `dict`**,
   transiently as a list of line-strings then joined — a ~2× peak of roughly
   **90–180 GB RAM** on a 90 Gbp genome, incurred *before* any TIR work. Likely
   OOM.
2. **CPU: Round 1 of the R detector is single-threaded** and contains the
   heaviest per-item operation in the tool (an `Rbeast::beast()` MCMC change-point
   fit per contig region), while the later, lighter rounds already use all cores.

## Finding 1 — whole-genome in-memory dict (`dt_utils.py`, `dante_tir.py`)

```python
# dante_tir.py:140
genome = dt.fasta_to_dict(args.fasta)          # entire assembly into RAM

# dt_utils.py:168-185
def fasta_to_dict(fasta_file):
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                fasta_dict[seq_name] = [""]        # list of per-line strings ...
            else:
                fasta_dict[seq_name].append(line.strip())
    fasta_dict = {k: "".join(v) for k, v in fasta_dict.items()}   # ... then joined -> ~2x peak
    return fasta_dict
```

Only ±6 kb flanks around each TPase domain are actually needed
(`extract_flanking_regions`, `dt_utils.py:724-791`), yet the whole genome is
resident. On the R side the genome is additionally `readDNAStringSet`-loaded up
to three times (`detect_tirs.R:144`, `dt_utils.R:980`, and in
`dante_tir_summary.R:93`) — 2-bit-packed (~22 GB for 90 Gbp) but re-read rather
than reused.

**Proposed fix.** Replace the whole-genome dict + random-access pattern with an
**indexed reader** (`pysam.FastaFile` / faidx) that fetches only the flanks per
domain (`dt_utils.py:753-773`). This removes the 90–180 GB peak outright. If the
dict is kept, at least avoid the transient list-of-lines (append to a running
buffer) to halve the peak. On the R side, read the genome once and pass the
object into `cluster_tir_sequences` instead of re-reading at `dt_utils.R:980`.

## Finding 2 — Round 1 is single-threaded (`dt_utils.R`)

`round1` (`dt_utils.R:1103`) and its callee `process_region_files`
(`dt_utils.R:887`) take **no `threads`/`mc.cores` argument** (call site
`detect_tirs.R:56` passes none), unlike Round 3 and Round 4 which parallelise
detection with `mclapply(mc.cores = threads)` (`dt_utils.R:1755`, `1950`). Two
serial hot loops:

```r
# dt_utils.R:900  — one beast() MCMC per contig region file, on ONE core
for (f in file_list) {
    ... find_switch_point(...) -> Rbeast::beast(...)   # heaviest op in the tool
}
# dt_utils.R:1170 — per common-TIR-ID pairwise eval, on ONE core
for (id in both_side_id) {
    ... eval_aln_length_alt(...)   # inner sapply(1:L, ... binom.test ...) is O(L) per pair
}
```

So the **earliest and heaviest** detection round runs single-threaded while the
later rounds use all `-c` cores — the worst asymmetry for wall-time on a big
genome with many TPase domains. A secondary serial spot is the Round-2 post-BLAST
pairwise loop (`dt_utils.R:1486-1502`; blastn itself is threaded at `1333`).

**Proposed fix.** Thread `threads` through `round1(...)` /
`process_region_files(...)` (`dt_utils.R:1103,887`; call site `detect_tirs.R:56`)
and convert the `for (f in file_list)` (`900`) and `for (id in both_side_id)`
(`1170`) loops to `mclapply(..., mc.cores = threads)`, mirroring Round 3. Same
treatment for the Round-2 loop (`1486-1502`).

> **⚠️ Blocked: `Rbeast::beast()` is not fork-safe.** We implemented exactly
> this (`mclapply` over `process_region_files`, `mc.cores` threaded through) and
> it **silently changed detection results** at `threads > 1`: on the `short`
> test the run dropped from **14 → 10** TIR records at `-c 2`, with **no error
> or warning**. `-c 1` (serial `lapply`) matched stock exactly (14). The loss is
> **not** RNG-related (`mc.set.seed = FALSE` → still 10) and **not** OpenMP
> (`OMP_NUM_THREADS=1` → still 10): `find_switch_point` passes a fixed
> `mcmc.seed`, so each `beast()` call should be deterministic, yet concurrent
> forked execution makes some return `NA` (→ the domain is silently skipped),
> pointing to shared global/temp state inside Rbeast that is unsafe under `fork`.
> **A naive `mclapply` here is a silent-data-loss bug and must not be shipped.**
> A safe parallelisation needs either an upstream Rbeast fork-safety fix, or a
> **process-isolated** backend (a PSOCK `parLapply` cluster — separate R
> processes, no shared fork state) with a result-count guard that fails loudly
> if any worker returns fewer detections than serial. Until then Round 1 stays
> serial (the memory fix above is independent and already landed).

## Finding 3 — secondary (bounded, lower priority)

- **CAP3-parse pool ignores `-c`.** `ProcessPoolExecutor(...)` at
  `dt_utils.py:512` is created with no `max_workers`, so it spawns
  `os.cpu_count()` workers (each inheriting a copy of the part's FASTA), defeating
  the `-c 32` throttle on a many-core node. The `ncpus` parameter of
  `parse_cap3_aln` (`dt_utils.py:444`) is accepted but unused. Fix: pass
  `max_workers = ncpus`.
- **Per-contig temp-FASTA inode growth.** One FASTA is written per assembled CAP3
  contig (`dante_tir.py:292-315`), accumulating until the final `shutil.rmtree`
  (and never cleaned in `--debug`). Not a correctness or fd bug — inode/space
  pressure only; relieved if written into per-class subdirs or streamed.

## Priority for the 90 Gbp use case

Finding 1 (memory) is the likely hard failure and should land first; Finding 2
(Round-1 parallelisation) is the wall-time win once it runs. Finding 3 is
polish. No fd-exhaustion fix is needed here.

#!/usr/bin/env python3
"""Reduce a CARP combined repeat library by per-class CAP3 / mmseqs2.

Drop-in CLI-compatible replacement for ``scripts/reduce_library_size.R``
producing **byte-identical** output. The R script has been retained
for parity testing; the Snakefile rule ``reduce_library`` continues to
call the R script until a separate cutover commit flips it to this one
after parity has been demonstrated.

Why this exists
---------------
The R + Biostrings stack used by the original script consumes 5–8 GB
of resident memory on import alone. Combined with ``mclapply``'s
fork-then-COW model, the four parallel workers each carry a
copy-on-write snapshot of the parent's heap. Snakemake's benchmark
collector aggregates RSS across the process tree, summing those
shared pages independently per child — producing 30 GB+ peak
``max_rss`` numbers on inputs of a few hundred KB. The actual
physical memory was lower (≈ parent RSS × 1–2), but high enough that
GHA hosted runners (16 GB) reaped the job before the heavy phase
completed.

This Python rewrite:

* Streams the input FASTA — full library never resides in memory.
* Uses ``multiprocessing`` with ``spawn`` start method — each worker
  is a fresh ~30 MB Python process with no inherited heap. No COW
  over-counting in benchmark output.
* Adds two-pass scheduling: large LTR classes go through CAP3
  sequentially in Phase 1 (caps peak at 1× CAP3 working set); small
  classes (CAP3 or mmseqs2) run in parallel in Phase 2.
* Keeps the CAP3 / mmseqs2 / blastn invocations character-for-character
  identical to the R script — same args, same input bytes, same
  output bytes. The only difference between the two implementations is
  the Python orchestrator's memory profile.

The byte-identity contract is enforced by
``tests/test_reduce_library_parity.sh`` (Commit B of the rollout).
Until that test is green on the reference fixture, the Snakefile
should not switch to this script.
"""
from __future__ import annotations

import argparse
import errno
import hashlib
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


# ──────────────────────────────────────────────────────────────────────────
# Constants — chosen to match the R script exactly
# ──────────────────────────────────────────────────────────────────────────
FASTA_LINE_WIDTH = 80                     # Biostrings writeXStringSet default
CAP3_ARGS = ("-p", "80", "-o", "50")      # identical to R: cap3 input -p 80 -o 50
MMSEQS_THREADS = "1"                       # R passes --threads 1 to easy-cluster
BLASTN_ARGS = (
    "-task", "blastn",
    "-outfmt",
    "6 qseqid sseqid pident length qstart qend evalue bitscore qlen slen qcovs",
    "-evalue", "1e-20",
    "-perc_identity", "95",
    "-word_size", "9",
    "-max_target_seqs", "20",
    "-gapextend", "1", "-gapopen", "2",
    "-reward", "1", "-penalty", "-1",
)
BLAST_LENGTH_MIN = 50                      # alignment rows shorter than this
                                            # don't count toward coverage; R's
                                            # calculate_total_coverage skip
QCOV_DROP_THRESHOLD = 0.98                 # singlets with qcov > 0.98 dropped;
                                            # strict > matches R's qcov[[i]] > 0.98


# ──────────────────────────────────────────────────────────────────────────
# FASTA I/O — streaming, no Biopython dependency
# ──────────────────────────────────────────────────────────────────────────
def _header_to_name(header_line: str) -> str:
    """Convert ``>foo bar baz\\n`` → ``foo bar baz``.

    Strips only the leading ``>`` and the trailing newline (and CR if
    present). **Preserves all internal and trailing whitespace** —
    Biostrings ``readDNAStringSet`` keeps the full header verbatim,
    and ``writeXStringSet`` writes it back the same way, so we have to
    too for byte-identity. mmseqs2's ``cluster_rep_seq.fasta`` adds a
    trailing space to some headers; if we stripped it, output would
    diverge from the R reference by exactly one byte per affected
    record.
    """
    return header_line[1:].rstrip("\n").rstrip("\r")


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield ``(name, sequence)`` per FASTA record.

    The ``name`` is the entire header line minus the leading ``>`` and
    trailing newline — matches Biostrings' read/write round-trip so
    output bytes match. Sequence is the concatenation of body lines
    with line breaks stripped, and **uppercased**: Biostrings'
    ``readDNAStringSet`` stores sequences in its ACGT alphabet with
    lowercase letters coerced to uppercase, so ``writeXStringSet``
    emits uppercase. CAP3's contigs/singlets carry mixed case (lower
    = soft-masked / low-confidence consensus); without ``.upper()``
    here, output diverges from the R reference on the lowercase bases
    that CAP3 introduces.
    """
    name: str | None = None
    seq_parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_parts).upper()
                name = _header_to_name(line)
                seq_parts = []
            else:
                seq_parts.append(line.rstrip("\n").rstrip("\r"))
        if name is not None:
            yield name, "".join(seq_parts).upper()


def write_fasta_record(name: str, seq: str, fh,
                        line_width: int = FASTA_LINE_WIDTH) -> None:
    """Write one FASTA record matching ``writeXStringSet(width=80)``.

    Header is ``>{name}\\n`` (no trailing whitespace, no description
    field). Sequence is wrapped at ``line_width`` with a trailing
    newline after every chunk including the last. For an empty
    sequence the output is just the header line.
    """
    fh.write(f">{name}\n")
    for i in range(0, len(seq), line_width):
        fh.write(seq[i:i + line_width])
        fh.write("\n")


def parse_classification(name: str) -> str:
    """Extract the classification component of a CARP-style FASTA name.

    Mirrors the R script's ``gsub('.+#', '', name)``: returns
    everything after the **last** ``#``. Names without ``#`` are
    returned unchanged.
    """
    return name.rsplit("#", 1)[-1]


# ──────────────────────────────────────────────────────────────────────────
# Stream-split by classification
# ──────────────────────────────────────────────────────────────────────────
@dataclass
class ClassJob:
    """Per-class job record. Pickle-safe so it crosses spawn workers."""
    index: int                  # 1-based, alphabetical class order (matches R)
    class_name: str             # the classification path, e.g. Class_I/LTR/Ty1_copia
    fasta_path: str             # input.fasta inside the class workdir
    workdir: str                # the class workdir itself
    size_bp: int                # total bp in fasta_path; used for big/small split
    # Filled in by workers:
    mode: str = ""              # "cap3" or "mmseqs"
    contig_file: str = ""       # CAP3 only
    singlet_file: str = ""      # CAP3 only
    aln_file: str = ""          # CAP3 only
    mmseqs_rep_file: str = ""   # mmseqs only
    singlets_final_file: str = ""   # set by Phase 3 BLAST filter (or = singlet_file)


def stream_split_by_class(input_path: Path, workdir: Path) -> list[ClassJob]:
    """Stream the input FASTA, partition into per-class files on disk.

    Open one file handle per unique class encountered; close all at end.
    Then sort classes alphabetically (mirrors R's ``split()`` factor-
    levels-default ordering) and rename the per-class files to indexed
    directories ``workdir/1/``, ``workdir/2/``, ... so downstream code
    matches the R script's per-index workdir layout.

    Returns one ``ClassJob`` per class, in alphabetical order.
    """
    workdir.mkdir(parents=True, exist_ok=True)

    # Stage 1: stream-write to per-class temp files (named by hashed
    # class so we never have to escape special chars for filenames).
    handles: dict[str, tuple[Path, object]] = {}
    sizes: dict[str, int] = {}
    try:
        for name, seq in iter_fasta(input_path):
            cls = parse_classification(name)
            if cls not in handles:
                tag = hashlib.sha1(cls.encode("utf-8")).hexdigest()[:16]
                tmppath = workdir / f"_class_{tag}.fasta"
                fh = open(tmppath, "w")
                handles[cls] = (tmppath, fh)
                sizes[cls] = 0
            tmppath, fh = handles[cls]
            write_fasta_record(name, seq, fh)
            sizes[cls] += len(seq)
    finally:
        for _, fh in handles.values():
            fh.close()

    # Stage 2: sort alphabetically, move temp files to indexed dirs.
    sorted_classes = sorted(handles)
    jobs: list[ClassJob] = []
    for idx, cls in enumerate(sorted_classes, start=1):
        class_dir = workdir / str(idx)
        class_dir.mkdir(parents=True, exist_ok=True)
        final_path = class_dir / "input.fasta"
        tmppath, _ = handles[cls]
        os.replace(tmppath, final_path)
        jobs.append(ClassJob(
            index=idx,
            class_name=cls,
            fasta_path=str(final_path),
            workdir=str(class_dir),
            size_bp=sizes[cls],
        ))
    return jobs


# ──────────────────────────────────────────────────────────────────────────
# External-tool wrappers — module-level so they're picklable for spawn
# ──────────────────────────────────────────────────────────────────────────
def _run_cap3(job: ClassJob) -> ClassJob:
    """Execute ``cap3 input.fasta -p 80 -o 50`` for one class.

    Reproduces the R script's invocation exactly: stdout and stderr
    are merged into the .aln file as a logging artefact, while CAP3's
    own implicit outputs (``input.fasta.cap.contigs``,
    ``input.fasta.cap.singlets``) are produced next to the input.
    """
    job.contig_file = job.fasta_path + ".cap.contigs"
    job.singlet_file = job.fasta_path + ".cap.singlets"
    job.aln_file = job.fasta_path + ".cap.aln"
    # Clean stale outputs from a previous invocation (R does the same)
    for p in (job.contig_file, job.singlet_file, job.aln_file):
        try:
            os.remove(p)
        except FileNotFoundError:
            pass
    cmd = ["cap3", job.fasta_path, *CAP3_ARGS]
    with open(job.aln_file, "w") as aln:
        rc = subprocess.run(cmd, stdout=aln, stderr=aln).returncode
    if rc != 0:
        raise RuntimeError(
            f"[{job.class_name}] CAP3 failed with status {rc} "
            f"(see {job.aln_file})"
        )
    job.mode = "cap3"
    job.singlets_final_file = job.singlet_file  # may be replaced by BLAST filter
    return job


def _run_mmseqs(job: ClassJob) -> ClassJob:
    """Execute ``mmseqs easy-cluster --threads 1`` for one class."""
    mmseqs_prefix = os.path.join(job.workdir, "mmseqs_cluster")
    mmseqs_tmp = os.path.join(job.workdir, "mmseqs_tmp")
    mmseqs_log = os.path.join(job.workdir, "mmseqs.log")
    # Clean stale outputs (R does the same via list.files+unlink)
    for stale in os.listdir(job.workdir):
        if stale.startswith("mmseqs_cluster") or stale == "mmseqs_tmp":
            full = os.path.join(job.workdir, stale)
            if os.path.isdir(full):
                shutil.rmtree(full, ignore_errors=True)
            else:
                try:
                    os.remove(full)
                except FileNotFoundError:
                    pass
    os.makedirs(mmseqs_tmp, exist_ok=True)
    cmd = [
        "mmseqs", "easy-cluster",
        job.fasta_path, mmseqs_prefix, mmseqs_tmp,
        "--threads", MMSEQS_THREADS,
    ]
    with open(mmseqs_log, "w") as log:
        rc = subprocess.run(cmd, stdout=log, stderr=log).returncode
    rep_path = mmseqs_prefix + "_rep_seq.fasta"
    if rc != 0 or not os.path.exists(rep_path):
        raise RuntimeError(
            f"[{job.class_name}] mmseqs2 failed with status {rc} "
            f"(see {mmseqs_log})"
        )
    job.mode = "mmseqs"
    job.mmseqs_rep_file = rep_path
    return job


def _should_use_cap3(class_name: str) -> bool:
    """Return True iff the class name triggers the CAP3 branch.

    Matches R's ``startsWith(classification_name, 'Class_I/LTR')``.
    """
    return class_name.startswith("Class_I/LTR")


def process_class(job: ClassJob) -> ClassJob:
    """Worker entry point. Dispatches to CAP3 or mmseqs2."""
    if _should_use_cap3(job.class_name):
        return _run_cap3(job)
    return _run_mmseqs(job)


# ──────────────────────────────────────────────────────────────────────────
# BLAST post-pass — runs in the main process after the CAP3 phases.
# Drops singletons that are ≥98 % covered by a contig from the same class.
# ──────────────────────────────────────────────────────────────────────────
def analyze_blast(blast_path: Path) -> dict[str, float]:
    """Compute per-query coverage from a blastn outfmt-6 result.

    Faithful port of the R ``analyze_blast`` + ``calculate_total_coverage``:

    * Empty file → empty dict (``return(0)`` in R, which silently means
      "no hits to drop").
    * Group rows by qseqid; for each query, allocate a length-qlen
      bitmap and mark covered positions.
    * Skip alignment rows whose ``length`` < 50 — same threshold the
      R script applies. Without this filter the coverage metric over-
      counts very short spurious hits.
    * ``qstart``/``qend`` are 1-based inclusive in BLAST. Python
      bitmap is 0-based; ``[qstart-1:qend]`` gives the same set of
      positions (count = qend - qstart + 1).
    * Final value is ``sum(covered) / qlen`` — fraction of query bases
      covered by at least one ≥50 bp hit.
    """
    if not os.path.exists(blast_path) or os.path.getsize(blast_path) == 0:
        return {}
    by_query: dict[str, list[tuple[int, int, int, int]]] = {}
    with open(blast_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11:
                continue
            qseqid = cols[0]
            length = int(cols[3])
            qstart = int(cols[4])
            qend = int(cols[5])
            qlen = int(cols[8])
            by_query.setdefault(qseqid, []).append((length, qstart, qend, qlen))
    qcov: dict[str, float] = {}
    for qseqid, rows in by_query.items():
        qlen = rows[0][3]  # all rows for the same query share qlen
        covered = bytearray(qlen)
        for length, qstart, qend, _ in rows:
            if length < BLAST_LENGTH_MIN:
                continue
            # qstart/qend 1-based inclusive → Python slice [qstart-1:qend]
            for i in range(qstart - 1, qend):
                covered[i] = 1
        qcov[qseqid] = sum(covered) / qlen if qlen > 0 else 0.0
    return qcov


def _filter_singlets(
    singlets_path: str,
    out_path: str,
    drop_ids: set[str],
) -> None:
    """Stream-copy a FASTA, omitting records whose full header name is
    in ``drop_ids``.

    Comparison uses the FULL header (matching ``iter_fasta`` /
    Biostrings ``names(...)``). ``drop_ids`` is built from BLAST
    ``qseqid`` which is the first whitespace-token, so a record
    survives the filter unless its full header equals a first-token —
    in practice this means trailing whitespace in the singlet header
    keeps the record. R's ``%in%`` exhibits the same behaviour, so
    matching it preserves byte-identity.
    """
    with open(singlets_path) as fin, open(out_path, "w") as fout:
        keep = False
        for line in fin:
            if line.startswith(">"):
                name = _header_to_name(line)
                keep = name not in drop_ids
            if keep:
                fout.write(line)


def run_blast_filter(jobs: list[ClassJob], threads: int) -> None:
    """For each CAP3 job that produced both contigs and singlets,
    run ``makeblastdb`` + ``blastn`` and drop singlets covered ≥ 98 %.

    This is the post-pass equivalent of the R script's loop after
    ``mclapply``. We run sequentially (matches R's ``lapply``); BLAST
    itself is multi-threaded via ``-num_threads``.
    """
    candidates = [
        j for j in jobs
        if j.mode == "cap3"
        and os.path.exists(j.contig_file) and os.path.getsize(j.contig_file) > 0
        and os.path.exists(j.singlet_file) and os.path.getsize(j.singlet_file) > 0
    ]
    if not candidates:
        return
    for job in candidates:
        # makeblastdb on contigs
        rc = subprocess.run(
            ["makeblastdb", "-in", job.contig_file, "-dbtype", "nucl"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        ).returncode
        if rc != 0:
            raise RuntimeError(
                f"[{job.class_name}] makeblastdb failed (status {rc})"
            )
        blast_out = job.singlet_file + ".blastn"
        cmd = [
            "blastn",
            "-query", job.singlet_file,
            "-db", job.contig_file,
            *BLASTN_ARGS,
            "-num_threads", str(threads),
            "-out", blast_out,
        ]
        rc = subprocess.run(cmd).returncode
        if rc != 0:
            raise RuntimeError(
                f"[{job.class_name}] blastn failed (status {rc})"
            )
        qcov = analyze_blast(Path(blast_out))
        drop_ids = {q for q, c in qcov.items() if c > QCOV_DROP_THRESHOLD}
        if drop_ids:
            filtered_path = job.singlet_file + "_filtered"
            _filter_singlets(job.singlet_file, filtered_path, drop_ids)
            job.singlets_final_file = filtered_path


# ──────────────────────────────────────────────────────────────────────────
# Output assembly — stream-concatenate in canonical order
# ──────────────────────────────────────────────────────────────────────────
def assemble_output(jobs: list[ClassJob], output_path: Path) -> int:
    """Write the final reduced library by streaming through the per-class
    artefacts in the **exact** order the R script emits them.

    Order (matches the ``c(do.call(c, contigs), …, do.call(c, mmseqs_reps))``
    construction at the bottom of the R script):

    1. CAP3 contigs across all classes — class-alphabetical order;
       names rewritten as ``<orig>_<idx>#<class>``.
    2. CAP3 (filtered) singlets across all classes — class-alphabetical
       order; names left untouched.
    3. mmseqs2 cluster representatives — class-alphabetical order;
       names left untouched.

    Returns total bp written (for the optional summary log line).
    """
    total_bp = 0
    # Atomic write: stage to ``output_path.tmp`` then rename. Avoids
    # leaving a half-written file if the script is killed.
    output_tmp = Path(str(output_path) + ".tmp")
    with open(output_tmp, "w") as fh:
        # 1. Contigs
        for job in jobs:
            if job.mode != "cap3":
                continue
            if not (os.path.exists(job.contig_file)
                    and os.path.getsize(job.contig_file) > 0):
                continue
            for name, seq in iter_fasta(Path(job.contig_file)):
                renamed = f"{name}_{job.index}#{job.class_name}"
                write_fasta_record(renamed, seq, fh)
                total_bp += len(seq)
        # 2. Singlets (filtered if BLAST trimmed any)
        for job in jobs:
            if job.mode != "cap3":
                continue
            singlets = job.singlets_final_file or job.singlet_file
            if not (os.path.exists(singlets) and os.path.getsize(singlets) > 0):
                continue
            for name, seq in iter_fasta(Path(singlets)):
                write_fasta_record(name, seq, fh)
                total_bp += len(seq)
        # 3. mmseqs reps
        for job in jobs:
            if job.mode != "mmseqs":
                continue
            rep = job.mmseqs_rep_file
            if not (os.path.exists(rep) and os.path.getsize(rep) > 0):
                continue
            for name, seq in iter_fasta(Path(rep)):
                write_fasta_record(name, seq, fh)
                total_bp += len(seq)
    os.replace(output_tmp, output_path)
    return total_bp


# ──────────────────────────────────────────────────────────────────────────
# main()
# ──────────────────────────────────────────────────────────────────────────
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("-i", "--input", required=True, type=Path,
                    help="Input library FASTA (combined library to reduce).")
    ap.add_argument("-o", "--output", required=True, type=Path,
                    help="Output reduced library FASTA.")
    ap.add_argument("-t", "--threads", type=int, default=4,
                    help="Worker count for the small-class phase (default 4). "
                         "Phase 1 (big CAP3 classes) always uses 1 worker.")
    ap.add_argument("-d", "--directory", type=Path, default=None,
                    help="Working directory for per-class staging "
                         "(default: tempfile.mkdtemp()).")
    ap.add_argument("-m", "--max-parallel-bp", type=int, default=50_000_000,
                    help="LTR classes whose input FASTA total is at least this "
                         "many bp run sequentially in Phase 1 (default 50 MB). "
                         "Lower this on memory-constrained machines.")
    ap.add_argument("--start-method", default="spawn",
                    choices=("spawn", "fork", "forkserver"),
                    help="multiprocessing start method (default: spawn). "
                         "spawn gives the cleanest memory profile; fork is "
                         "faster to start workers but inherits the parent's "
                         "heap (the very thing we're avoiding here).")
    return ap.parse_args(argv)


def _phase_partition(jobs: list[ClassJob], max_parallel_bp: int
                     ) -> tuple[list[ClassJob], list[ClassJob]]:
    """Split jobs into (big CAP3, everything else).

    A job is "big" iff it triggers the CAP3 branch AND its input
    library exceeds the threshold. Such jobs go through Phase 1 one
    at a time. Everything else (small CAP3 + all mmseqs2) runs in
    parallel in Phase 2.
    """
    big = [j for j in jobs
           if _should_use_cap3(j.class_name) and j.size_bp >= max_parallel_bp]
    small = [j for j in jobs if j not in big]
    return big, small


def _run_phase(jobs: list[ClassJob], threads: int,
               start_method: str, label: str) -> list[ClassJob]:
    """Run a phase via ``ProcessPoolExecutor`` with the chosen start
    method. Returns the worker-updated job records, in input order.
    """
    if not jobs:
        return []
    sys.stderr.write(
        f"[reduce_library] {label}: {len(jobs)} class(es), workers={threads}, "
        f"start_method={start_method}\n"
    )
    ctx = multiprocessing.get_context(start_method)
    with ProcessPoolExecutor(max_workers=threads, mp_context=ctx) as pool:
        return list(pool.map(process_class, jobs))


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if not args.input.exists():
        sys.stderr.write(f"reduce_library_size.py: input not found: {args.input}\n")
        return 2
    if args.input.stat().st_size == 0:
        # R reads an empty FASTA into DNAStringSet() and produces an empty
        # output FASTA. Mirror that.
        args.output.write_text("")
        sys.stderr.write("reduce_library_size.py: input empty; output empty\n")
        return 0

    workdir = (Path(args.directory) if args.directory
               else Path(tempfile.mkdtemp(prefix="reduce_library_")))
    workdir.mkdir(parents=True, exist_ok=True)

    sys.stderr.write(
        f"[reduce_library] input: {args.input}  ({args.input.stat().st_size} B)\n"
        f"[reduce_library] workdir: {workdir}\n"
        f"[reduce_library] threads: {args.threads}  "
        f"max_parallel_bp: {args.max_parallel_bp}  "
        f"start_method: {args.start_method}\n"
    )

    # Phase 0: stream-split input by class
    jobs = stream_split_by_class(args.input, workdir)
    sys.stderr.write(
        f"[reduce_library] split: {len(jobs)} class(es)\n"
    )
    for j in jobs:
        sys.stderr.write(
            f"  [{j.index:3d}] {j.class_name}  "
            f"({j.size_bp} bp, {'CAP3' if _should_use_cap3(j.class_name) else 'mmseqs2'})\n"
        )

    # Phase split: big CAP3 vs everything else
    big_jobs, small_jobs = _phase_partition(jobs, args.max_parallel_bp)

    # Phase 1: big CAP3 classes, sequential
    big_done = _run_phase(big_jobs, threads=1,
                          start_method=args.start_method,
                          label="phase 1 (big CAP3, sequential)")
    # Phase 2: small remainder, parallel
    small_done = _run_phase(small_jobs, threads=args.threads,
                            start_method=args.start_method,
                            label="phase 2 (small classes, parallel)")

    # Reassemble in original (alphabetical, 1-based) class order so
    # downstream BLAST + assembly see the same ordering as R's
    # ``group_results``.
    completed: dict[int, ClassJob] = {}
    for j in big_done + small_done:
        completed[j.index] = j
    jobs = [completed[i + 1] for i in range(len(jobs))]

    # Phase 3: BLAST filter (singletons vs contigs, qcov > 0.98 → drop).
    # Matches the R loop after mclapply.
    sys.stderr.write("[reduce_library] phase 3: BLAST filter on CAP3 classes\n")
    run_blast_filter(jobs, threads=args.threads)

    # Phase 4: assemble output
    sys.stderr.write("[reduce_library] phase 4: assemble output\n")
    output_bp = assemble_output(jobs, args.output)

    # Final summary line — equivalent to R's three message() calls at end
    input_bp = sum(j.size_bp for j in jobs)
    pct = (output_bp / input_bp * 100) if input_bp > 0 else 0.0
    sys.stderr.write(
        f"-----------------------------------------------\n"
        f"Input library size: {input_bp} bp\n"
        f"Output library size: {output_bp} bp\n"
        f"Reduction: {pct:.2f}%\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

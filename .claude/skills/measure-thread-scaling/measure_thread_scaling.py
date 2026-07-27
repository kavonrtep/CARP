#!/usr/bin/env python3
"""PR-C (item 5): measure thread scaling + peak RSS of the shared all-vs-all
alignment engine (``run_all_vs_all_alignment`` in scripts/global_local_aln.py).

Why: threads only add throughput if the work releases the GIL. parasail's Python
binding is ctypes (it DOES release the GIL during the C alignment), so we expect
real scaling — but the per-pair Python overhead (building result dicts, the
score-threshold filter) is GIL-bound and can dominate for very short flanks. The
sister-pipeline finding was that blastn barely scaled with -num_threads; this
harness answers the same question for THIS loop, on real data, before anyone
changes the threading.

Each (thread-count, rep) runs in its OWN subprocess so peak RSS is isolated
(RUSAGE_SELF ru_maxrss). Reports, per thread count:
  wall_s, cpu_self_s, eff_cores (=cpu_self/wall), speedup + parallel efficiency
  vs the smallest thread count, peak RSS (MB), record count,
plus a determinism check (output must be byte-identical across ALL thread counts
— guards item 3). The peak-RSS column also serves as the PR-B (items 1/2/4)
memory confirmation: point it at a genuinely high-copy family's prime/flank FASTA
and add --max-group-size to exercise the grouped/streaming path.

Run under the dante_line conda env (needs parasail; mmseqs only with --prefilter
or --max-group-size):

  conda run -p <dante_line_env> python measure_thread_scaling.py \
      -f <output>/DANTE_LINE/ENDO_RT_5prime.fasta -t 1,4,8,16,32 --end 5 --reps 2

Good inputs = the all-vs-all inputs a real run produces, e.g. dante_line's
ENDO_RT_5prime.fasta / ENDO_RT_3prime.fasta / ENDO_RT_RH_*prime.fasta, or the
fallback's TPase_5prime.fasta / TPase_3prime.fasta. For a high-copy family one of
these carries thousands of sequences — that is the O(N^2) step to profile.
"""
import argparse
import hashlib
import importlib.util
import json
import os
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def load_engine():
    """Import scripts/global_local_aln.py from the enclosing repo checkout."""
    here = Path(__file__).resolve()
    for parent in here.parents:
        cand = parent / "scripts" / "global_local_aln.py"
        if cand.exists():
            spec = importlib.util.spec_from_file_location("global_local_aln", cand)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module
    sys.exit("measure_thread_scaling: could not find scripts/global_local_aln.py "
             "above this script")


def run_one(fasta, out, end, threads, max_group_size, prefilter):
    """One alignment run; returns timing/memory as a dict (child-process body)."""
    engine = load_engine()
    t0 = time.perf_counter()
    records = engine.run_all_vs_all_alignment(
        fasta_file=fasta, output_file=out, end=end, threads=threads,
        verbose=False, use_prefilter=prefilter, max_group_size=max_group_size,
    )
    wall = time.perf_counter() - t0
    ru = resource.getrusage(resource.RUSAGE_SELF)
    ruc = resource.getrusage(resource.RUSAGE_CHILDREN)
    md5 = (hashlib.md5(open(out, "rb").read()).hexdigest()
           if os.path.exists(out) else "NOFILE")
    n = records if isinstance(records, int) else len(records or [])
    return {
        "threads": threads,
        "wall": wall,
        "cpu_self": ru.ru_utime + ru.ru_stime,
        "cpu_children": ruc.ru_utime + ruc.ru_stime,   # mmseqs subprocesses, if any
        "maxrss_mb": ru.ru_maxrss / 1024.0,            # Linux: ru_maxrss is KB
        "records": n,
        "md5": md5,
    }


def main():
    ap = argparse.ArgumentParser(
        description="Measure all-vs-all alignment thread scaling + peak RSS (PR-C / item 5)")
    ap.add_argument("-f", "--fasta", required=True,
                    help="all-vs-all input FASTA (a prime/flank FASTA from a real run)")
    ap.add_argument("-t", "--threads", default="1,4,8,16",
                    help="comma-separated thread counts (default 1,4,8,16)")
    ap.add_argument("--end", default="5", choices=["5", "3"])
    ap.add_argument("--reps", type=int, default=2,
                    help="repetitions per thread count; the median is reported (default 2)")
    ap.add_argument("--max-group-size", type=int, default=None,
                    help="bound the all-vs-all into clustering groups (exercises the "
                         "grouped/streaming path; default: single pass)")
    ap.add_argument("--prefilter", action="store_true",
                    help="enable the mmseqs prefilter (default off, to isolate the "
                         "parasail alignment scaling)")
    ap.add_argument("--outdir", default=None,
                    help="dir for the per-run TSVs (default: a temp dir)")
    # hidden single-run child mode
    ap.add_argument("--_child", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--out", help=argparse.SUPPRESS)
    ap.add_argument("--child-threads", type=int, help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args._child:
        res = run_one(args.fasta, args.out, args.end, args.child_threads,
                      args.max_group_size, args.prefilter)
        print("CHILD_JSON " + json.dumps(res))
        return

    if not Path(args.fasta).exists():
        sys.exit(f"FASTA not found: {args.fasta}")
    thread_list = [int(x) for x in args.threads.split(",") if x.strip()]
    outdir = args.outdir or tempfile.mkdtemp(prefix="thread_scaling_")
    os.makedirs(outdir, exist_ok=True)
    nseq = sum(1 for line in open(args.fasta) if line.startswith(">"))
    print(f"# input: {args.fasta}  ({nseq} sequences)  end={args.end}  reps={args.reps}"
          f"  max_group_size={args.max_group_size}  prefilter={args.prefilter}")

    agg = {}
    all_md5 = set()
    any_children = False
    for tc in thread_list:
        for rep in range(args.reps):
            out = os.path.join(outdir, f"t{tc}_r{rep}.tsv")
            cmd = [sys.executable, str(Path(__file__).resolve()), "--_child",
                   "-f", args.fasta, "--out", out, "--end", args.end,
                   "--child-threads", str(tc)]
            if args.max_group_size is not None:
                cmd += ["--max-group-size", str(args.max_group_size)]
            if args.prefilter:
                cmd += ["--prefilter"]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            line = next((ln for ln in proc.stdout.splitlines()
                         if ln.startswith("CHILD_JSON ")), None)
            if line is None:
                sys.exit(f"child run failed (threads={tc} rep={rep}):\n"
                         f"{proc.stdout}\n{proc.stderr}")
            d = json.loads(line[len("CHILD_JSON "):])
            agg.setdefault(tc, []).append(d)
            all_md5.add(d["md5"])
            any_children = any_children or d["cpu_children"] > 0.1

    base_tc = min(thread_list)
    base_wall = statistics.median(d["wall"] for d in agg[base_tc])

    print()
    print(f"{'threads':>7} {'wall_s':>9} {'cpu_self':>9} {'eff_cores':>9} "
          f"{'speedup':>8} {'effic%':>7} {'peakRSS_MB':>11} {'records':>9}")
    for tc in thread_list:
        ds = agg[tc]
        wall = statistics.median(d["wall"] for d in ds)
        cpu_self = statistics.median(d["cpu_self"] for d in ds)
        rss = max(d["maxrss_mb"] for d in ds)
        rec = ds[0]["records"]
        eff_cores = cpu_self / wall if wall > 0 else 0.0
        speedup = base_wall / wall if wall > 0 else 0.0
        effic = 100.0 * speedup / (tc / base_tc) if tc > 0 else 0.0
        print(f"{tc:>7} {wall:>9.2f} {cpu_self:>9.2f} {eff_cores:>9.2f} "
              f"{speedup:>8.2f} {effic:>7.0f} {rss:>11.1f} {rec:>9}")

    print()
    print("eff_cores = cpu_self/wall (in-process parallelism of the parasail loop; "
          "~threads = good scaling, ~1 = GIL/overhead-bound).")
    if any_children:
        print("NOTE: mmseqs subprocess CPU (grouping/prefilter) is excluded from "
              "cpu_self/eff_cores — those measure the alignment loop only.")
    if len(all_md5) == 1:
        print(f"determinism: OK — output byte-identical across all thread counts/reps "
              f"({next(iter(all_md5))[:12]})")
    else:
        print(f"determinism: WARNING — {len(all_md5)} distinct output checksums: "
              f"{sorted(m[:12] for m in all_md5)}")
    print(f"(per-run TSVs in {outdir})")


if __name__ == "__main__":
    main()

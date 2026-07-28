# Full-size determinism check (manual)

CARP's determinism is enforced automatically on every push/PR by the
`determinism` job in `.github/workflows/pipeline.yml` — it runs the **small and
medium** fixtures twice each (different `PYTHONHASHSEED` + thread count) and
diffs the manifest data outputs via `scripts/assert_run_determinism.py` (plus a
canonical-order check on the RepeatMasker library in the `fixture` job). That
catches most regressions cheaply.

But several real non-determinism bugs (clustering order, culling × threads,
`as_completed` row order, `set`/hash-seed leaks) only became visible on a
**full-size genome**. This skill runs the same two-run diff on a real genome, on
demand — do it before trusting a release or a large change to any rule that
clusters, aligns all-vs-all, or masks.

## Run (on the host / HPC, where `run_pipeline.py` works)
```bash
.claude/skills/check-determinism/check_determinism.sh -c <your_genome.yaml>
# optional: -a <threads run 1>  -b <threads run 2>  (default: all cores, then half)
```
It runs the whole pipeline twice on the same config — pass 1 `PYTHONHASHSEED=1 -t
<all>`, pass 2 `PYTHONHASHSEED=2 -t <half>` — into `<output_dir>_det_a` and
`_det_b`, then compares. A full-size genome is ~30–60 min per pass (~1–2 h total).
This is a **manual** gate; never wire it into per-PR CI.

## What it compares
`scripts/assert_run_determinism.py` diffs the manifest `OUTPUTS` (the
downstream-consumed set — unified + per-class GFF3s, `summary_statistics.csv`,
BigWig tracks, libraries, BEDs). GFF3 `#` header lines (e.g. `##date`) are
normalized out; the HTML reports and `run_provenance.json` are excluded (they
embed timestamps). Everything else must be **byte-identical** between the two
runs. Exit 0 = deterministic; exit 1 lists the outputs that drifted.

## When it fails
A listed output changed between two runs of the same input — a determinism
regression. Find the source using the checklist in the **Determinism** section of
`CLAUDE.md` (canonical-sort clustering inputs; `-num_threads 1` for culling; sort
intermediate row order; no `set`-driven output order; no wall-clock in data
outputs; atomic writes). The two `_det_a` / `_det_b` trees are left in place so
you can `diff` the offending file directly.

## Notes
- Lives under `.claude/` → not in the release SIF, not a pipeline rule.
- Uses `run_pipeline.py` (the real entrypoint), so it also exercises the
  provenance/manifest path — though those files themselves are excluded from the
  comparison.

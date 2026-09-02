# k-mer prefilter test — run this on the exec host

## Commands

```bash
cd /nfsroot/projects/darwin/runs/tmp3/kmer_test
./RUN.sh --check      # verify environment + inputs, runs nothing   <-- DO THIS FIRST
./RUN.sh --setup      # only if --check says the environment was not found
./RUN.sh              # run the calibration arm and score it
```

`--check` costs a second and verifies three things before you commit hours to
the run: an environment with the right tools, every input file readable, and
that the scripts actually **load** in that environment and carry the
`--prefilter-kmer` flag. That last one matters — an earlier version passed all
its environment checks and then every job died instantly on a missing module.

If no environment is found, `--setup` builds one with conda/mamba into
`./conda_env` from the pinned `dante_line_env.yaml` (seqkit, parasail-python,
mmseqs2, pyyaml — the same versions the validated container was built with).

That runs the calibration arm (3 genomes x 2 settings, 6 jobs in parallel at 32
threads each) and prints the scored result against the pre-registered rule.
Expect **roughly 1-2 hours**, set by one genome: GCA_977018505 needed ~102,500
core-seconds in its production run, which is ~55 min at 32 threads, and the k=7
arm does more alignment work than that.

Other modes:

```bash
./RUN.sh --list       # what would run, what is already done
./RUN.sh --score      # re-print the scores without re-running
./RUN.sh --validate   # the held-out genomes -- ONLY if calibration passes
```

Knobs: `THREADS=32` per job, `SIF=/path/to/image.sif`, `ENV_PREFIX=/path/for/env`.

If a job fails, its exit code and the last lines of its log are printed and its
`DONE_` marker is removed, so simply re-running retries only what failed.

Safe to re-run: finished jobs are skipped via `results/DONE_*`. To redo one,
delete its `results/DONE_<tag>_k<k>` file.

**boechera is pre-seeded** — both its arms were already run and are in
`results/`, so `./RUN.sh` will skip them and only run the two large genomes.
Delete `results/DONE_boechera_k*` if you would rather it re-ran there.

## What is being tested, in one paragraph

`dante_line` works out how far a LINE extends past its protein-domain core by
aligning many copies' flanking sequence against each other. Two flanks are only
ever aligned if mmseqs first finds a shared exact k-mer in the 30 nt next to the
core. At mmseqs' default k-mer length that window is too short for most related
pairs to qualify, so they are discarded before identity is ever checked, and
**71% of LINE elements end up with no extension at all**. On the small test
genome, passing `-k 7` raised admitted pairs from 171 to 907 and extended
elements from 28.6% to 66.5%.

The question is not whether it extends more. It is whether the extra extension
is REAL element. More extension is the same lever that produced the 16-22 kb
chimeric "LINE" consensi this whole investigation started from.

## How it is judged

Ground truth is LINE loci where BOTH ends were actually observed -- a poly-A
tail found and a target-site duplication clearing a measured chance floor. For
those loci we know the true extension, so:

- **PURITY**  = emitted extension bp lying inside the true element / emitted bp.
  The metric that matters. An over-extended consensus mislabels whatever its
  tail matches genome-wide; a truncated one only loses a little masking.
- **COVERAGE** = emitted bp inside the true element / true extension bp.
  How much real element was recovered.

`PREREGISTERED.md` was committed BEFORE any of this was run and fixes the genome
split, the metrics and the accept/reject bar. Judge the **worst** genome, never
the average -- an average lets one good genome hide a bad one.

**Adopt k=7 only if, on every calibration genome, coverage rises by >= 10
percentage points, purity stays >= 95%, and purity drops by no more than 1
point. Then it must reproduce on the held-out genomes with no retuning. Any
failure means keep the mmseqs default.**

## Result so far (the small genome, already run)

```
boechera (25 confirmed loci)
  k-default    purity 100.0%   coverage 13.9%
  k-7          purity 100.0%   coverage 16.9%
```

Purity is untouched, but coverage rose only 3.0 points against a >= 10 bar, so
boechera **fails** the rule. It is also the least representative genome here --
343 LINE loci against 12-13k in the other two -- which is why the other
calibration genomes are worth running rather than rejecting on this one alone.
If they also land near +3, the rejection is solid and the default stays.

## Why this directory exists at all

The shipped container does **not** have the `--prefilter-kmer` flag. The
`scripts/` here are patched copies that do. `RUN.sh` uses the container only for
its tools (python3 + parasail, mmseqs, seqkit) and runs these scripts instead of
the ones baked into the image.

## Contents

| path | what |
|------|------|
| `RUN.sh` | the runner |
| `PREREGISTERED.md` | the rule, fixed in advance |
| `scripts/` | the full patched CARP scripts directory (dante_line.py pulls in several modules via importlib, so shipping a hand-picked subset does not work) |
| `classification_vocabulary.yaml` | needed by the patched scripts |
| `dante_line_env.yaml` | pinned env spec used by `--setup` |
| `gt/gt_<genome>.tsv` | confirmed-boundary ground truth |
| `score_extensions.py` | purity/coverage scorer |
| `results/` | per-run output, logs and scores |

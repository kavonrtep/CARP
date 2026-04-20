# CI test fixtures

Three genome slices of increasing size, carved out of the Pisum JI2822 run
(`tmp/output/genome_cleaned.fasta`) so that every pipeline stage sees
non-trivial input while the full end-to-end run fits inside a CI job's time
budget.

| dataset | FASTA | total size | sequences | coverage highlights |
|---|---|---:|---:|---|
| **micro** | `genome_micro.fasta` | 200 kb | 2 | densest hAT TPase block (CM055026.1) + Ale-rich LTR block — smoke test, not all pipeline stages produce output |
| **tiny** | `genome_tiny.fasta` | 500 kb | 5 | adds EnSpm-dense tile, Ivana LTR tile and TRC_2 satellite tile |
| **small** | `genome_small.fasta` | 1 Mb | 10 | full REXdb coverage: 4 LTR lineages with ≥3 copies each + 4 TIR families + satellites + LINE; exercises library reduction and multi-scaffold batching. **~5 min wall time with 4 cores** |
| **medium** | `genome_medium.fasta` | 2.4 Mb | 11 | 200-kb windows around 10 real primary-DANTE_TIR hits (5 hAT + 5 EnSpm/CACTA). Triggers substantive DANTE_TIR **fallback** detection (≥18 combined TIR features). Primary DANTE_TIR detection still returns 0 — the algorithm needs hundreds of TPase copies for statistical TIR-consensus convergence, which is unreachable at CI scale. |

Matching `config_*.yaml` files are sibling to each FASTA. Invoke from the
repo root:

```bash
python run_pipeline.py -c tests/fixtures/config_micro.yaml -t 4
python run_pipeline.py -c tests/fixtures/config_tiny.yaml  -t 4
python run_pipeline.py -c tests/fixtures/config_small.yaml -t 8
```

## Reproducing the slices

Windows were chosen by tiling the reference genome into 100-kb bins and
ranking each bin by a diversity score (# of distinct LTR families × 3 +
# of distinct TIR families × 3 + # of satellite clusters × 2 + LINE
presence × 2). The exact window coordinates are encoded in the sequence
names (`<chrom>_<start1>-<end1>`), so a slice can be re-extracted with:

```bash
seqkit subseq --bed <(awk -F'[_-]' -v OFS='\t' '/^>/ {print $1, $2-1, $3}' \
  tests/fixtures/genome_micro.fasta) \
  path/to/source/genome.fasta
```

## Adding more fixtures

Pick a new bin from the reference unified annotation:

```bash
python3 scripts/pick_fixture_windows.py --out-bed /tmp/extra.bed \
  --target-bp 250000 --include-satellite --include-line
seqkit subseq --bed /tmp/extra.bed tmp/output/genome_cleaned.fasta \
  > tests/fixtures/genome_extra.fasta
```

(The helper script is not checked in yet; the logic used to build the
current fixtures lives in the git history of this directory.)

## What NOT to commit into this directory

- Output directories (`output_micro/`, `output_tiny/`, `output_small/`)
  — they're generated on every run and would balloon the repo.
- Large intermediate FASTAs (cluster libraries, mmseqs dbs, etc.).
- The reference genome slice source (`tmp/output/genome_cleaned.fasta`).

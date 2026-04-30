# tests/

Lightweight test scripts for ad-hoc validation. Not (yet) integrated
into a formal test framework — each script is self-contained, prints
a clear OK / FAIL signal, and is intended to be runnable from the
repo root in either a developer checkout or a CI step.

## Classification round-trip

| File | Purpose |
|---|---|
| `test_classification.py` | Verifies `scripts/classification.py` parses every case in `classification_cases.tsv` and round-trips through canonicalise → is_canonical. |
| `test_classification.R` | Same, for the R mirror in `scripts/classification.R`. |
| `classification_cases.tsv` | Shared test vector used by both Python and R parsers. |

Run after any change to `classification_vocabulary.yaml` or the
`classification` modules:

```bash
python3 tests/test_classification.py
Rscript tests/test_classification.R
```

Both should report `72/72 cases passed` (the row count grows when
new entries are added).

## reduce_library Python ↔ R parity

| File | Purpose |
|---|---|
| `test_reduce_library_parity.sh` | Runs `scripts/reduce_library_size.R` and `scripts/reduce_library_size.py` on the same input and asserts the output FASTAs are byte-identical. |
| `test_reduce_library_memprofile.sh` | Profiles peak RSS / wall time / CPU time of both implementations on the same input via `/usr/bin/time -v`. |

Both expect `cap3`, `mmseqs`, `makeblastdb`, and `blastn` on PATH —
activate the `tidecluster` conda env (or whichever env hosts those
tools) before running. The Python implementation is a drop-in
replacement for the R version; production behaviour will not flip
until `test_reduce_library_parity.sh` runs green on the medium
fixture and a separate Snakefile cutover commit lands.

Recommended development sweep — three fixtures of increasing
complexity, exercising different code paths:

```bash
# mmseqs-only path (4 classes, no LTR, no CAP3)
./tests/test_reduce_library_parity.sh tests/fixtures/output_micro/Libraries/combined_library.fasta

# CAP3 + BLAST + mmseqs (9 classes, 11 LTR records)
./tests/test_reduce_library_parity.sh tests/fixtures/output_small/Libraries/combined_library.fasta

# realistic mix (the fixture release.yml's test-in-container uses)
./tests/test_reduce_library_parity.sh tests/fixtures/output_medium/Libraries/combined_library.fasta
```

The medium fixture is the one to watch: passing on it is the
green-light signal that the Snakefile rule can be flipped from R to
Python in a follow-up commit.

For a memory-profile snapshot (run on real hardware, not in a
container or sandbox without `/usr/bin/time`):

```bash
./tests/test_reduce_library_memprofile.sh tests/fixtures/output_medium/Libraries/combined_library.fasta
```

Output is a 2-row table comparing peak RSS / wall / user / sys time
plus the R-peak / Py-peak ratio. The expected difference is a 5–10x
peak RSS reduction in favour of Python — this is the headline result
backing the rewrite.

## Fixtures

| Path | Notes |
|---|---|
| `fixtures/output_tiny/`   | Smallest fixture, spike-out for quick sanity checks. |
| `fixtures/output_micro/`  | Single-genome rDNA-only library; mmseqs path. |
| `fixtures/output_small/`  | Adds CAP3 / BLAST coverage (LTR retro records). |
| `fixtures/output_medium/` | The realistic CI fixture — used by `pipeline.yml` and `release.yml`. |

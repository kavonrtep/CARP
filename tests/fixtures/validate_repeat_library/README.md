# Fixtures for `scripts/validate_repeat_library.py`

Golden inputs and expected outputs for the library check / conversion tool.
`tests/test_validate_repeat_library.py` reruns the tool on each input and
requires byte-identical FASTA and TSV output (the HTML report is not compared:
it carries the CARP version).

| file | what it is |
|---|---|
| `tandem_input.fasta` | Real-world tandem library with non-CARP headers (22 *Arabidopsis* sequences from GenBank: CEN180, 5S/45S rDNA, telomere, minisatellites, ...; class written without the `Satellite/` prefix). |
| `tandem_expected.fasta` / `.tsv` | `-t tandem` result: every sequence converted (`Satellite/<name>`, `rDNA/5S_rDNA`, `rDNA/45S_rDNA`), exit 0. |
| `custom_input.fasta` | Synthetic custom library covering each conversion path: RepeatMasker/Dfam names, REXdb/DANTE/DANTE_TIR spellings, bare lineage name, Wicker codes, FASTA hygiene (description, gaps, `U`, duplicate ID), plus three unconvertible records (`tRNA`, no class, protein sequence). |
| `custom_expected.tsv` + `custom_expected_empty.fasta` | `-t custom --out-fasta ...` result: 3 errors, exit 2, FASTA created empty (the Galaxy contract). |
| `custom_dropped_expected.fasta` | `-t custom --drop-invalid` result: the 15 convertible records, exit 0. |

Regenerate after an intended behaviour change (then review the diff):

```bash
F=tests/fixtures/validate_repeat_library; V=scripts/validate_repeat_library.py
$V $F/tandem_input.fasta -t tandem --out-fasta $F/tandem_expected.fasta --out-tsv $F/tandem_expected.tsv --out-html /tmp/r.html
$V $F/custom_input.fasta -t custom --out-fasta $F/custom_expected_empty.fasta --out-tsv $F/custom_expected.tsv --out-html /tmp/r.html
$V $F/custom_input.fasta -t custom --drop-invalid --out-fasta $F/custom_dropped_expected.fasta --out-tsv /tmp/r.tsv --out-html /tmp/r.html
```

# Checking and converting a repeat library

CARP accepts two optional user libraries (see
[configuration.md](configuration.md)):

- `tandem_repeat_library` — satellites and other tandem repeats. TideCluster
  uses it to name the tandem repeat clusters it finds.
- `custom_library` — any repeats, usually transposable elements. They are added
  to the RepeatMasker library.

Both are FASTA files where each header carries a classification:
`>ID#CLASS`. The class must use CARP's vocabulary
([`classification_vocabulary.yaml`](../classification_vocabulary.yaml)).
A library with other class names (`#CEN180`, `#LTR/Gypsy`, `#RLC`) stops the
pipeline at the `validate_classifications` step. That step runs after DANTE,
DANTE_LTR, DANTE_TIR and DANTE_LINE, so a bad library is only reported hours
into a run.

`validate_repeat_library.py` checks a library before the run. It rewrites the
headers it can convert safely, and it explains the ones it cannot.

## Usage

The tool is in the CARP container:

```bash
singularity exec -B $PWD carp.sif \
    validate_repeat_library.py my_satellites.fasta -t tandem -o TR_lib

singularity exec -B $PWD carp.sif \
    validate_repeat_library.py my_TEs.fasta -t custom -o custom_lib
```

Then use `TR_lib.fasta` as `tandem_repeat_library`, or `custom_lib.fasta` as
`custom_library`. From a repo checkout, run `scripts/validate_repeat_library.py`
the same way. It needs only Python 3 and PyYAML.

| option | meaning |
|---|---|
| `-t tandem` / `-t custom` | which library the file is meant to be (required) |
| `-o PREFIX` | writes `PREFIX.fasta`, `PREFIX.tsv`, `PREFIX.html` (default prefix: `<input name>.carp_<type>`) |
| `--out-fasta`, `--out-tsv`, `--out-html` | set each output path on its own (overrides `-o`; this is what a Galaxy wrapper uses) |
| `--drop-invalid` | write the library without the sequences that cannot be converted |

## Outputs

- **Library FASTA.** The converted library. Sequences are written 60 bp per
  line, in input order, and the letter case of the input is kept.
- **TSV.** One row per input sequence. Columns: `index`, `input_id`,
  `input_class`, `output_id`, `output_class`, `status` (`unchanged`,
  `converted`, `error`, `dropped`), `length_bp`, `rule` (which rule set the
  class), `changes`, `errors`.
- **HTML report.** A single self-contained page with the verdict, counts, the
  table of every sequence (filter: all / changed / errors), and the format
  guideline for that library type. For a custom library, the guideline lists
  every valid class.

### When some sequences cannot be converted

By default no library is produced. Dropping sequences changes what the genome
is annotated with, so it is never done silently:

| situation | library FASTA | exit status |
|---|---|---|
| all sequences valid, or converted | written | `0` |
| some cannot be converted, `--drop-invalid` given | written without them (marked `dropped` in the report) | `0` |
| some cannot be converted, run with `-o PREFIX` | not written (an older file at that path is removed) | `2` |
| some cannot be converted, run with `--out-fasta PATH` | created **empty** | `2` |
| input missing, unreadable, or has no FASTA records | none | `3` |

A non-zero exit status always means "do not use the FASTA". The empty file in
the `--out-fasta` case exists only because a Galaxy tool must create every
output it declares. A Galaxy wrapper should treat exit status `2` as a failed
job, so the empty file is never passed on to a CARP run. The HTML report
explains why.

## What is converted

Before the class is looked at, both library types get the same clean-up:

| problem | what the tool does |
|---|---|
| text after the first space in a header | dropped (RepeatMasker reads only the first word) |
| characters other than letters, digits and `. _ : - + \| ( ) [ ]` in the ID | replaced with `_` |
| duplicate IDs | the second copy becomes `ID_dup2`, the third `ID_dup3`, ... |
| alignment gaps `-` `.` in the sequence | removed |
| RNA base `U` | changed to `T` |
| trailing `?` on the class (RepeatMasker "uncertain") | removed |
| Windows line endings | changed to Unix line endings |

These cannot be fixed, so the sequence is reported as an error: an empty
sequence, non-DNA letters (for example a protein FASTA), an empty ID, or more
than one `#` in a header.

### Tandem library (`-t tandem`)

Only `Satellite`, `Satellite/<name>` and `rDNA/...` classes are allowed. The
tool trusts the user that the sequences are tandem repeats. It does not check
for tandem structure, because a library entry is often a single monomer, which
has no internal repeat to find.

| input header | output header | rule |
|---|---|---|
| `>X04322.1#CEN180` | `>X04322.1#Satellite/CEN180` | a bare name becomes the satellite name |
| `>t#telomere` | `>t#Satellite/telomere` | telomeric repeats count as satellites |
| `>s#satellite`, `>s#Satellite_x` | `>s#Satellite`, `>s#Satellite/x` | spelling normalised |
| `>M65137.1#5S_rDNA` | `>M65137.1#rDNA/5S_rDNA` | rDNA name |
| `>X52322.1#45S_rDNA` | `>X52322.1#rDNA/45S_rDNA` | rDNA name (`45S`, `35S`, `47S`, `48S`) |
| `>r#18S`, `#ITS1`, `#IGS`, `#25S` | `>r#rDNA/45S_rDNA/18S`, ... | rDNA subunit name |
| `>CEN180_a` (no `#`) | `>CEN180_a#Satellite/CEN180_a` | the ID is used as the name |
| `>a#LTR/Gypsy`, `>a#RLC`, `>a#Class_I/...` | — error | a TE class: move the sequence to `custom_library` |

rDNA names are matched as whole tokens, and `45S` is tested before `5S`. So
`45S_rDNA` is never read as 5S, and `AthSat5Sx` is not read as rDNA at all.

### Custom library (`-t custom`)

Every class must be a CARP classification. These rules are tried in order:

1. **The pipeline's own normaliser** (`classification.py`). It converts
   DANTE-style `Class_I|LTR|Ty3/gypsy|chromovirus|Tekay` and DANTE_TIR-style
   `Class_II_Subclass_1_TIR_hAT` to slash form.
2. **`Satellite` spellings**, for example `satellite/centromeric` becomes
   `Satellite/centromeric`.
3. **A CARP class name or path ending.** `Athila` becomes
   `Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Athila`, and `Ty1/copia/Ale`
   becomes `Class_I/LTR/Ty1_copia/Ale`. This is used only when exactly one CARP
   class matches.
4. **The dictionary** of RepeatMasker/Dfam names and Wicker codes (below).
5. **rDNA names**, the same rules as for a tandem library.

Anything else is an error. The tool never guesses a class.

#### Dictionary

The dictionary is the `library_class_aliases` section of
[`classification_vocabulary.yaml`](../classification_vocabulary.yaml). Each
rule is a regular expression matched against the whole class name (ignoring
case), and the first match wins. Every target is checked against the
vocabulary when the tool starts. The pipeline itself never reads this section.

When the foreign name is more specific than anything in CARP, the rule maps to
the most specific CARP class that is still correct. For example ERV and BEL/Pao
become `Class_I/LTR`. Names with no CARP equivalent (`tRNA`, `snRNA`, `scRNA`,
...) are left out on purpose, so they are reported instead of being mapped to a
wrong class.

**RepeatMasker / Dfam names (examples)**

| foreign class | CARP class |
|---|---|
| `LTR/Copia`, `LTR/Gypsy` | `Class_I/LTR/Ty1_copia`, `Class_I/LTR/Ty3_gypsy` |
| `LTR/ERV*`, `LTR/Pao`, `LTR/Unknown`, `LTR` | `Class_I/LTR` |
| `LTR/Caulimovirus` | `Class_I/pararetrovirus` |
| `LINE/*` (L1, RTE, ...), `SINE/*` | `Class_I/LINE`, `Class_I/SINE` |
| `LINE/Penelope`, `PLE` | `Class_I/Penelope` |
| `RC/Helitron` | `Class_II/Subclass_2/Helitron` |
| `DNA/hAT*`, `DNA/CMC-EnSpm`, `DNA/MULE-MuDR`, `DNA/PIF-Harbinger`, `DNA/TcMar*` | `Class_II/Subclass_1/TIR/` + hAT, EnSpm_CACTA, MuDR_Mutator, PIF_Harbinger, Tc1_Mariner |
| `DNA/Maverick` / `DNA/Crypton` / `DNA` | `Class_II/Subclass_2` / `Class_II/Subclass_1` / `Class_II` |
| `rRNA` | `rDNA` |

**Wicker et al. (2007) three-letter codes**

Wicker T. et al. (2007) *A unified classification system for eukaryotic
transposable elements.* Nat Rev Genet 8:973–982. A code is accepted alone
(`RLC`) or followed by a family name after `_`, `-` or `/` (`RLC_Angela`,
`DTC-CACTA_1`). The family part is not used: Wicker family names are specific
to one species and do not match REXdb lineages, so the class stops at the
superfamily. A code glued to other text (`DTA1`) is not read as a code.

| code | CARP class |
|---|---|
| `RLC` / `RLG` | `Class_I/LTR/Ty1_copia` / `Class_I/LTR/Ty3_gypsy` |
| `RLB`, `RLR`, `RLE`, `RLX` | `Class_I/LTR` |
| `RYD`, `RYN`, `RYV`, `RYX` | `Class_I/DIRS` |
| `RPP`, `RPX` | `Class_I/Penelope` |
| `RIR`, `RIT`, `RIJ`, `RIL`, `RII`, `RIX` | `Class_I/LINE` |
| `RST`, `RSL`, `RSS`, `RSX` | `Class_I/SINE` |
| `RXX` | `Class_I` |
| `DTT`, `DTA`, `DTM`, `DTE`, `DTP`, `DTB`, `DTH`, `DTC` | `Class_II/Subclass_1/TIR/` + Tc1_Mariner, hAT, MuDR_Mutator, Merlin, P, PiggyBac, PIF_Harbinger, EnSpm_CACTA |
| `DTR`, `DTX` | `Class_II/Subclass_1/TIR` |
| `DYC`, `DYX` | `Class_II/Subclass_1` |
| `DHH`, `DHX` | `Class_II/Subclass_2/Helitron` |
| `DMM`, `DMX` | `Class_II/Subclass_2` |
| `DXX` | `Class_II` |

To support a new foreign name, add a rule to `library_class_aliases`, put it
before any more general rule it overlaps, and add a case to
`tests/test_validate_repeat_library.py`.

## Running from Galaxy

The tool runs from the CARP image as-is. Inside the container it is found by
name on `PATH` (`/opt/pipeline/scripts`), and it finds the vocabulary file in
`/opt/pipeline` the same way the pipeline does. A wrapper should:

- pass `--out-fasta`, `--out-tsv` and `--out-html` with the paths of Galaxy's
  output datasets (the `-o` prefix does not fit Galaxy's dataset naming);
- treat exit status `0` as success, and `2` and `3` as failure;
- put the HTML output on the HTML allowlist, so the report's style and its
  filter buttons are kept. Without the allowlist the table is still complete,
  but unstyled.

## Tests

`tests/test_validate_repeat_library.py` (run in the CI `unit` job) covers each
conversion rule, the error cases, the exit codes, the `--out-*` contract, and
that the output is the same under different `PYTHONHASHSEED` values. It also
reruns the tool on the golden fixtures in `tests/fixtures/validate_repeat_library/`
and requires byte-identical FASTA and TSV output.

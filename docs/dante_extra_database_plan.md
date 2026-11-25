# Implementation plan — optional DANTE extra protein database

Status: implemented. Request: issue #8 (CARP 1.9.9 user, Asteraceae).

Adds one optional config parameter, `dante_extra_database`: a FASTA file of
protein domains that DANTE appends to REXdb (`dante -e`). It is for genomes
whose transposable elements REXdb covers poorly.

**Invariants**

* (a) A run without the parameter gives the same output as today. The rule
  graph (DAG) of such a run does not change.
* (b) No tool version bump and no new conda environment. DANTE 0.2.12, already
  pinned in `envs/dante.yaml`, has `-e`.
* (c) No change to `classification_vocabulary.yaml`, the manifest or
  `OUTPUT_SCHEMA_VERSION`.
* (d) A bad file stops the run before DANTE starts, with a message naming the
  record.

## What was tested

DANTE 0.2.12 on `tests/fixtures/genome_small.fasta` (about 1 Mb; 294 domain
records without `-e`). Each test file held up to 60 proteins that DANTE had found
in that same genome. These proteins match it almost perfectly, so the size of
each effect is an upper-end example, not an estimate for a real genome.

| Question | Result |
|---|---|
| Does `-e` itself change the search? | No. An empty file gives byte-identical feature lines. |
| Is `-e` reproducible? | Yes. Four runs (three with different `PYTHONHASHSEED` and 3, 6 or 9 threads) gave byte-identical feature lines. Only the position of repeated `##` header lines differs, which also happens without `-e`. |
| Does record order matter? | No. A shuffled file gave identical feature lines. |
| Does a failing `lastdb` pass silently? | No. `dante` does not check `lastdb`, but `lastal` then fails and `dante` exits 1. |
| Scratch and time cost | 82 MB for a copy of REXdb plus its index, in `$TMPDIR`; `lastdb` takes about 4 s. |
| Does the vocabulary cover REXdb? | Yes. All 56 classes of Viridiplantae_v4.0 are canonical, and none has a child class in the vocabulary. |

**Bad input is the main risk.** DANTE exits 0 for every case below except the
last one:

| Content of the extra file | Effect |
|---|---|
| New lineage (`...\|Ale\|Cicho`) | Adds an unknown class and loses LTR-RTs (next table) |
| Internal class (`Class_I\|LTR\|Ty1/copia`) | Lowers matching REXdb hits to that level |
| `Ty1_copia` instead of `Ty1/copia` | Lowers matching hits to `Class_I\|LTR` |
| Unknown domain name (`FOO`) | `Ambiguous_domain` records, which the DANTE filter drops |
| Same ID on several records | Classes are mixed up: `NONE(Annotations from different classes)` |
| ID already used by REXdb (`REXdb_ID1`) | The ID's class is overwritten, so records get a wrong class |
| `;` in the ID | Breaks the GFF3 attribute column |
| DNA instead of protein | Accepted without a warning |
| Two spaces between header fields | `ValueError`, exit 1 |

DANTE_LTR on the same data (lineage mode, `-M 1`; 0.6.4.0 and the pinned 0.6.5.0
gave the same counts). Complete means `Rank` other than `D`:

| Extra database (60 proteins) | Complete LTR-RTs | Partial (`Rank=D`) |
|---|---|---|
| none | 24 | 10 |
| existing REXdb classes | 24 | 10 |
| 9 of the 60 proteins with an internal class | 15 | 20 |
| new lineage below an existing one | 21 | 9 |
| all LTR proteins with an internal class | 13 | 9 |
| `Ty1_copia` / `Ty3_gypsy` spelling | 13 | 9 |

Reason: lineage mode accepts only lineages listed in its own
`lineage_domain_order.csv`, and all domains of an element must agree.
DANTE assigns the common prefix of all good hits, so one less specific protein
makes the whole call less specific.

## Design decisions

1. **Accept the 56 classes that REXdb itself uses, and the internal nodes of its
   tree, spelled exactly as in REXdb.** New lineages and the underscore spelling
   are errors. Internal nodes give a warning, because they can lower a lineage
   call (15 instead of 24 complete LTR-RTs in the test when 9 of 60 proteins had
   an internal class). Decision of the project
   owner, confirmed with these counts; the alternative, rejecting them, would
   also reject the example in
   `dante --help` (`...|OTA|Tat`).
2. **Reject, do not repair, a wrong class.** The message prints the closest
   valid class. Repair would need `classification.py`, and PyYAML is not in
   `envs/dante.yaml`.
3. **Read the allowed classes, domains and IDs from the installed DANTE**
   (`tool-data/protein_domains/Viridiplantae_v4.0_{class,pdb}`), not from a list
   in CARP. They then always match the pinned version.
4. **Validate in a separate rule.** It does not depend on the genome, so it fails
   in seconds. `validate_classifications` runs hours later and sees only classes
   that reached the output.

## Changes

### C1 — config (`Snakefile`)

Default `""` (off), and must be a string, as for `tidecluster_rdna_library`. An
empty value and a missing key mean the same.

### C2 — validator (`scripts/prepare_dante_extra_database.py`, new, `chmod +x`)

Standard library only. Arguments: input, `--reference-class`, `--reference-pdb`,
`-o`. Checks every record and reports all errors at once (exit 2):

* header has exactly three fields: `>ID DOMAIN CLASS`;
* ID matches `[A-Za-z0-9_.:-]+`, and is unique and absent from REXdb after
  DANTE's own `-` → `_` replacement;
* DOMAIN is one of the 13 REXdb domain types;
* CLASS is one of the 56 REXdb classes or an internal node (warning);
* sequence is not empty, uses letters plus `*` and `.` (`lastdb` refuses `-`, `/`
  and digits), and is not DNA
  (over 90 % `ACGTUN` in a sequence of at least 50 letters);
* file has at least one record.

A DOMAIN never seen with that CLASS in REXdb (for example CHD in Ty1/copia) gives
a warning. Output: the same records in the same order with Windows line
endings, blank lines and repeated spaces cleaned; written atomically. The log
reports warnings. A summary TSV (source file name, SHA-256, REXdb version, record
counts per class and domain) feeds the HTML report.

### C3 — rules (`Snakefile`)

* `prepare_dante_extra_database` (only defined when the parameter is set) writes
  `DANTE/extra_database.fasta` and `DANTE/extra_database_summary.tsv`, in
  `envs/dante.yaml`.
* `dante` takes the checked copy as an input and gets `-e <copy>`; both are empty
  when the parameter is not set. It also gets `-D Viridiplantae_v4.0`.
* `make_repeat_report` gets `--dante_extra_db_summary` when the parameter is set.

`-D Viridiplantae_v4.0` is already the default. Stating it stops a future DANTE
default change from making the validator check the wrong REXdb version.
The file is a rule `input`, so editing it reruns DANTE. The cleaned copy stays in
`DANTE/` as a record of what was used; `cleanup_outputs.py` must not list it.
Use the standard PATH-prepend shell pattern.

### C4 — launcher and provenance

* `run_pipeline.py`: check that the file exists, and add its directory to the
  printed `-B` options.
* Same function: `show_singularity_settings` reads `custom_library` and
  `tandem_repeat_library` without checking that they are set. A config without
  them raises `KeyError` instead of printing the hint. Fix it here.
* `scripts/record_provenance.py`: add the key to `_filter_config`.
* `scripts/make_repeat_report.R`: a note in the Overview section with the file
  name, REXdb version, number of proteins and a table of classes. When the report
  is rendered by hand, provenance decides whether the note is shown, because a
  later run without the parameter leaves the summary file behind.

### C5 — documentation

* `docs/configuration.md`: parameter row, header format, the class rule, and a
  note that `dante_ltr_mode: core` is the alternative when no proteins exist.
* `config_full.yaml`, README table, `docs/output_reference.md`
  (`DANTE/extra_database.fasta`), and `CHANGELOG.md` under `## Unreleased`.

### C6 — tests (each must be added to `.github/workflows/unit.yml`)

* `tests/test_prepare_dante_extra_database.py`: one good file, one case per row
  of the bad-input table, CRLF and blank lines, and identical output under two
  `PYTHONHASHSEED` values. Uses a small made-up reference, so no DANTE is needed.
* `tests/test_dante_extra_database_wiring.py`: default is off; the flag is empty
  when off and `-e <file>` when on; the `dante` shell uses it. Same style as
  `test_dante_ltr_mode_wiring.py`.
* `tests/test_report_extra_db.R`: content of the note, and when it is shown.

## Acceptance

1. Parameter unset: dry-run job list unchanged, and the small and medium fixtures
   give manifest outputs identical to 1.9.10.
2. Parameter set, small fixture: full run passes, twice with different hash seed
   and thread count, and `assert_run_determinism.py` reports no difference.
3. Each bad file stops the run before `dante` starts.
4. Container run with the file in its own `-B` directory.

Items 1 to 3 pass. Item 4 is open: it needs a newly built image.

Not tested: a full pipeline run on a genome that REXdb covers poorly. Only the
validator and the DANTE step were run on such a genome (a 30 Mb Asteraceae
window, no internal classes in the file). The effect of internal classes on real
data is unknown.

## Limits

* No new lineages. They would need a vocabulary entry, a row in DANTE_LTR's
  `lineage_domain_order.csv` (upstream), and report changes.
* A wrong class on an otherwise valid protein cannot be detected. The result
  depends on the quality of the user's classification.
* Viridiplantae only. Choosing another REXdb (`dante -D`) is a separate feature.
* The DANTE GFF3 header shows `pdb.fasta` as the database version when `-e` is
  used. Nothing in CARP reads that line; provenance records the real file.

# Robust repeat-classification handling — implementation plan

## Goal

Replace the ad-hoc sed/regex translation of repeat classification strings
scattered across Snakefile rules and analysis scripts with a single,
vocabulary-driven normalisation module. After the refactor every
classification that flows through the pipeline is validated against the
authoritative list in `classification_vocabulary.yaml`.

The plan follows the **B+A** option agreed earlier: a canonical vocabulary
(data) plus a centralised normalisation module (code), used from both
Python and R. No refactor of the underlying pipeline logic, no structured
(tuple) representation — just string in, canonical string out, with
validation.

## Success criteria

- Adding a new classification is a one-line edit of `classification_vocabulary.yaml`.
- No `sed` / `awk` on classification strings remains in `Snakefile`.
- No bespoke regex normalisation remains in `scripts/*.R` or `scripts/*.py`.
- A single function is the only code path that turns a tool-native string
  into canonical form; callers never touch separators.
- Any classification that cannot be resolved against the vocabulary produces
  a hard error (with the offending string, source tool, and the closest
  canonical match logged), not a silent pass-through.
- Running the pipeline on the existing `JI2822_2026-02-02` output should
  reproduce the same `summary_statistics.csv` classifications (verified by
  diff), with zero normalisation warnings.

## Deliverables

### D1. `classification_vocabulary.yaml`  (already created)
Authoritative data: canonical list, leaf aliases, special classes,
aggregation buckets, tool dialects.

### D2. `scripts/classification.py`
Python module. Public API:

```python
load_vocabulary(path=None) -> Vocabulary        # cached
canonicalise(s, source=None) -> str             # raise on unresolvable
is_canonical(s) -> bool
strip_similarity_suffix(s) -> str               # TideCluster "(NN.N%25)"
iter_canonical() -> Iterator[str]
```

Internals:
- Loads the YAML once, builds three lookups: canonical set, leaf-alias
  map, tool-prefix table.
- `canonicalise()` dispatches on `source` if given (DANTE / DANTE_TIR /
  TideCluster / RepeatMasker / custom_library); otherwise auto-detects by
  looking for `|`, underscore-prefix, or similarity suffix in that order.
- DANTE path: strip optional suffix, split on `|`, apply leaf_aliases
  per-component, join with `/`.
- DANTE_TIR path: match longest hierarchy_prefix from the vocabulary,
  substitute `/` for the matched prefix underscores, leave the remainder
  untouched (this is how `Class_II_Subclass_1_TIR_Tc1_Mariner` becomes
  `Class_II/Subclass_1/TIR/Tc1_Mariner` without touching the `Tc1_Mariner`
  leaf).
- RepeatMasker / TideCluster path: strip similarity suffix, validate.
- After normalisation the result must be in the canonical set OR match
  an allowed Satellite sub-path OR match an allowed special class. If
  not, raise `UnknownClassification`.

CLI entry point (`python -m classification normalise ...`) so the
Snakefile can shell out when Python isn't embedded.

### D3. `scripts/classification.R`
Thin R implementation of the same API (`canonicalise`, `is_canonical`,
`strip_similarity_suffix`). Reads the same YAML via the `yaml` R package.
R and Python must pass a shared test vector (see D7).

Not a reticulate bridge — a direct reimplementation. Two reasons: (1) the
Snakemake-managed conda envs would each need `reticulate` + a Python env
with the module, doubling deployment complexity; (2) the algorithm is
~50 lines and the vocabulary is the source of truth, so drift is checked
by the shared test vector.

Add `r-yaml` to `envs/tidecluster.yaml` (already carries `r-jsonlite`, so
adding one more trivially-sized package is fine).

### D4. Retire Snakefile sed calls

| Snakefile line | Current | Replacement |
|---|---|---|
| 338 (`make_tir_combined_library`) | `sed -E '/^>/ s@(#Class_II/Subclass_[12]/TIR)_@\\1/@'` | `python -m classification canonicalise-fasta-headers --source DANTE_TIR` |
| 749 (`filter_ltr_rt_library`) | `sed 's/\|/_/g' \| sed 's/|/\//g'` | `python -m classification canonicalise-fasta-headers --source DANTE_LTR` |

Both become one-liner shell calls into the module. The 692 `TRC_` → `TRC_S_`
rename is unrelated to classifications and stays as-is.

### D5. Retire per-script regex normalisation

Replace in-place:

| File | Current fragment | Replacement |
|---|---|---|
| `scripts/clean_DANTE_names.R` L6-11 | `gsub("|", "/", gsub("/", "_", ...))` | `vapply(x, canonicalise, character(1), source="DANTE")` |
| `scripts/clean_rm_output.R` L69-90 | `convert_names()` helper + `Unspecified` fallback | `canonicalise(..., source="RepeatMasker")` ; keep the `Unspecified` fallback logic since it extracts info from seq name |
| `scripts/make_unified_annotation.R` L78-94 (`fix_sep`, `convert_tir_cls`) | Per-tool regexes | `canonicalise(x, source=<tool>)` per loader |
| `scripts/make_repeat_report.R` L70-98 | Tool-specific splits | Same |
| `scripts/dante_tir_fallback.py` L647 | Hard-coded `Class_II/Subclass_1/TIR/{subtype}` | `canonicalise(f"Class_II_Subclass_1_TIR_{subtype}", source="DANTE_TIR")` — this also validates that `subtype` produces a known leaf |

Each replacement is a ~3-line diff. Keep behaviour-preserving; no logic
changes in the surrounding pipeline.

### D6. Vocabulary-level assertions at pipeline boundaries

Add one validation step per upstream tool: after that tool's output lands
in `output/`, run `python -m classification validate <file> --source <tool>`.
The validator scans the relevant GFF3 attribute, canonicalises every value,
and fails the rule if any value is unresolvable.

This is cheap (seconds per GFF3) and catches tool version upgrades that
silently introduce new classifications, which was the failure mode behind
the two recent bugfix commits.

### D7. Shared test vector

Create `tests/classification_cases.tsv` with columns
`source | raw_input | expected_canonical | notes`, populated from the
real data we inspected plus synthetic edge cases:

- DANTE_TIR with every leaf (hAT, Tc1_Mariner, EnSpm_CACTA, MuDR_Mutator,
  PIF_Harbinger, MITE, MITE/Stowaway)
- DANTE_LTR with Ty1/copia and Ty3/gypsy leaves
- TideCluster annotation with similarity suffix
- RepeatMasker canonical passthrough
- Intermediate classifications (`Class_I/LTR`, `Class_II/Subclass_1/TIR`)
- Satellite sub-paths (`Satellite/FabTR-1`)
- `Simple_repeat`, `Low_complexity`, `Unknown`
- Negative cases that should raise (`Class_III`, `Ty5_foo`)

Both the Python and R modules run against this TSV in CI; divergence
fails the build.

### D8. Documentation

Update `README.md`:
- Expand the classification list to reflect the full vocabulary
  (add `Class_I/pararetrovirus`, `rDNA_45S/IGS`, `rDNA_45S/ITS1`,
  `rDNA_45S/ITS2`, all intermediate nodes).
- Point maintainers at `classification_vocabulary.yaml` as the place to
  add new classifications.

Update `CLAUDE.md`:
- Section "Classification handling" describing the module and the
  "never sed a classification string" rule.

## Execution order

1. **D1** — YAML vocabulary (done).
2. **D2** — Python module + CLI.
3. **D7** — test vector (drives D2 and D3).
4. **D3** — R module against the same vector.
5. **D5** — replace in-script regexes (small, reviewable diffs).
6. **D4** — replace Snakefile sed calls.
7. **D6** — add validator rules.
8. **D8** — docs.

Each step is independently shippable and reversible; if D5 breaks a
downstream script the diff to revert is tiny. The validator (D6) goes in
last so earlier refactor steps don't trip their own checks mid-conversion.

## Risk + mitigations

| Risk | Mitigation |
|---|---|
| New DANTE / DANTE_TIR release introduces a leaf not in vocabulary | D6 validator fails loudly at the rule boundary, pointing at the offending value. Add the leaf to YAML; one-line fix. |
| R and Python normalisers drift | D7 shared test vector, run in CI. |
| Existing output directories contain legacy non-canonical strings | Normaliser is forgiving on read (accepts underscores-as-separator for known DANTE_TIR prefixes); strict only on write. |
| User-supplied custom library has unexpected classification | D6 validator runs on `custom_library` too — early hard failure is better than downstream mis-annotation. |
| `Satellite` vs `Satellites` confusion | Explicitly separated in vocabulary (`special_classes` vs `aggregation_buckets`); `distinct_from` field makes the relationship searchable. |

## Out of scope

- Representing classification as a structured list/tuple internally
  (Option C from the earlier review). Could be revisited if B+A proves
  insufficient, but the current string-based interface is compatible with
  all existing code.
- Replacing `jsonlite` with a different serialisation for the R side.
- Changing RepeatMasker library header conventions (still
  `>name#class/subclass`).
- TideCluster's similarity-suffix format (upstream behaviour; we strip it
  on read).

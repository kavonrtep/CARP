# TideCluster rDNA library: hard-coded class names make `--rdna_library` silently fail

TideCluster 1.16.0 added array-level rDNA identification, which CARP consumes
(`load_rdna_map()` in `scripts/make_unified_annotation.R` reads
`<prefix>_rdna.tsv`). It works well on the bundled library. While tracing an
rDNA misannotation we found that the feature's **user-supplied-library option
(`--rdna_library`) only works if the library uses TideCluster's exact
classification spelling**, and fails *silently* — producing zero rDNA calls —
for any other convention.

Observed in TideCluster 1.21.1 (`envs/tidecluster_run.yaml`).

## Finding

`identify_rdna` derives an rDNA "top type" from each library entry's
RepeatMasker-style `name#class` header, then discards anything that is not one
of two hard-coded literals.

`tc_utils.py:563`:

```python
RDNA_TOP_TYPES = ("rDNA_45S", "rDNA_5S")
```

`tc_utils.py` in `blastn_rdna_reference_coverage()`:

```python
cls = h.split("#", 1)[1] if "#" in h else ""
ref_top[h] = cls.split("/")[0]          # top type = first '/'-separated token
...
top = ref_top.get(qseqid)
if top not in RDNA_TOP_TYPES:
    continue                            # silently skipped — no warning
```

So the top type is the **first `/`-separated component** of the class, and it
must equal `rDNA_45S` or `rDNA_5S` exactly:

| library class | derived top type | matches? |
|---|---|---|
| `rDNA_45S/18S` (bundled) | `rDNA_45S` | yes |
| `rDNA_5S/5S` (bundled) | `rDNA_5S` | yes |
| `rDNA/45S_rDNA/18S` | `rDNA` | **no — skipped** |
| `45S_rDNA/18S` | `45S_rDNA` | **no — skipped** |

When nothing matches, `blastn_rdna_reference_coverage()` returns an empty dict,
`assign_rdna_to_trcs()` gets nothing, and `<prefix>_rdna.tsv` is written with
its header and **no data rows**. The run exits 0.

## Impact on downstream consumers

`--rdna_library` is offered as a user-facing option (`TideCluster.py:1517` and
`:1648`, "classes `rDNA_45S/*` and `rDNA_5S/*`"), but the constraint is enforced
only by that silent `continue`. A user pointing it at a library that uses any
other naming gets an empty `_rdna.tsv`, which is **indistinguishable from
"this genome has no rDNA arrays"**. There is no warning, no non-zero exit, and
nothing in the log that names the cause.

This is not hypothetical for CARP. CARP ships `data/rdna_library.fasta`
containing the **same 117 sequences as TideCluster's bundled library** —
109 unique, every one md5-identical, same per-class counts — relabelled to
CARP's canonical vocabulary (`classification_vocabulary.yaml`):

| CARP class | TideCluster class | n |
|---|---|---:|
| `rDNA/45S_rDNA/18S` | `rDNA_45S/18S` | 23 |
| `rDNA/45S_rDNA/25S` | `rDNA_45S/25S` | 17 |
| `rDNA/45S_rDNA/5.8S` | `rDNA_45S/5.8S` | 18 |
| `rDNA/5S_rDNA/5S` | `rDNA_5S/5S` | 59 |

The two files are the same database that has diverged in labelling only. CARP
exposes `tidecluster_rdna_library` as a config knob, so pointing it at CARP's
own copy — the natural thing to do, since the sequences are identical — would
**silently disable rDNA detection** for the whole run. CARP currently avoids
this only by leaving the knob empty so the bundled library is used.

The same trap catches real user libraries. A curated *Pisum* library in
production use here carries `rDNA_45S/{18S,25S,5.8S,IGS,ITS1,ITS2}`, which
happens to match; a library using the fuller hierarchical form would not.

## Requests

1. **Fail loudly instead of silently.** If no library reference yields a top
   type in `RDNA_TOP_TYPES`, emit a warning naming the classes actually found
   and the ones expected, or abort. A user-supplied library that contributes
   zero usable references is a configuration error, not a genome with no rDNA.
   This mirrors the fix already made in 1.20.1, where failed steps stopped
   producing output that denied its own missing results.

2. **Accept the subunit hierarchy more permissively.** Rather than matching the
   first token exactly, match on whether any token of the class path identifies
   the subunit family — e.g. treat a class as 45S if it contains `45S` and as 5S
   if it contains `5S`. That would accept `rDNA_45S/18S`, `rDNA/45S_rDNA/18S`
   and `45S_rDNA/18S` alike, and needs no change to the bundled library. If a
   strict contract is preferred instead, a `--rdna_class_map` option or a
   documented normalisation step would serve.

3. **Document the requirement precisely** in the `--rdna_library` help text: not
   just "classes `rDNA_45S/*` and `rDNA_5S/*`" but that the class token **before
   the first `/`** must be exactly `rDNA_45S` or `rDNA_5S`, and that other
   spellings are ignored.

Item 1 alone would have turned this from a silent wrong answer into an obvious
one, and is the smallest useful change.

## Reproduction

```bash
# Take the bundled library and relabel it to an equivalent hierarchical form.
sed 's/#rDNA_45S\//#rDNA\/45S_rDNA\//; s/#rDNA_5S\//#rDNA\/5S_rDNA\//' \
    <tidecluster>/data/rdna_library.fasta > rdna_renamed.fasta

# Same sequences, same genome, only the class labels differ.
TideCluster.py run_all -f genome.fasta -pr TC_bundled
TideCluster.py run_all -f genome.fasta -pr TC_renamed --rdna_library rdna_renamed.fasta

# TC_bundled_rdna.tsv  -> rDNA TRCs listed
# TC_renamed_rdna.tsv  -> header only, no rows, exit 0, no warning
```

## Evidence from the CARP side

Re-running scaffold_9 of the *Pisum sativum* CDC_Dakota assembly (2,061,187 bp,
~90 % 45S rDNA array) with the bundled library and TideCluster 1.21.1:

```
TideCluster_rdna.tsv
TRC     rDNA_type   coverage
TRC_1   45S         1.0
TRC_2   45S         1.0
TRC_3   45S         1.0
```

Coverage 1.0 on all three arrays, and CARP's unified annotation then labels
90.33 % of the scaffold `rDNA/45S_rDNA` — agreeing with RepeatMasker's
independent 89.6 % to within 0.7 %. The detection is decisive when the library
naming matches, which is what makes the silent-empty failure mode on a mismatch
so easy to miss.

# DANTE_LTR feature request — detect & represent tandem LTR-RT (`LTR_RT_TR`) arrays

**Summary.** When same-lineage LTR-RTs are arranged head-to-tail sharing a
boundary LTR (`LTR-INT-LTR-INT-LTR…`), DANTE_LTR currently reports each
domain-set-flanked-by-LTRs as a **separate intact `transposable_element`**, so
consecutive copies **overlap by one shared LTR**. Downstream, that shared LTR is
annotated/counted twice. We propose DANTE_LTR detect these arrays and represent
them once (a tandem container with member copies), emitting a structural tag
(`LTR_RT_TR`). This is the established "tandem LTR-RT" structure (Macko-Podgórni
et al., *Mobile DNA* 2025; tool IDENTAM) — common and Gypsy/Tekay/Athila-skewed
in plant genomes — so native support would benefit every DANTE_LTR user.

## Evidence (real DANTE_LTR output)

*Boechera stricta* `GCA_018361405.1`, `CM031550.1:7,521,409–7,605,371`. DANTE_LTR
(`-M 1`) on the region emits **four overlapping complete elements**, all
`Class_I|LTR|Ty3/gypsy|chromovirus|Tekay`, all `−` strand:

| element | span (abs) | overlaps next by |
|---|---|---|
| TE1 | 7,545,752–7,555,702 | 3,212 bp |
| TE2 | 7,552,491–7,562,948 | 3,188 bp |
| TE3 | 7,559,761–7,569,930 | 3,087 bp |
| TE4 | 7,566,844–7,576,995 | — |

Each overlap equals a shared LTR: e.g. TE1's 3′LTR `7,553,136–7,555,702`
coincides with TE2's 5′LTR `7,552,491–7,555,702` (`LTR_Identity` 90–98 % across
the array). Net: a 4-copy tandem (`LTR-INT-LTR-INT-LTR-INT-LTR-INT-LTR`, 5 LTRs /
4 internals), **~9.5 kb of shared-LTR sequence double-annotated**. (A partial 5th
copy abuts it; a separate Tekay element follows after a gap — correctly distinct.)

## Why it isn't recombination of two insertions / a nested insertion

Per the 2025 paper (worth mirroring in DANTE_LTR's logic): high flanking-sequence
identity and a **single set of TSDs** rule out recombination of two independent
insertions; absence of extra TSDs rules out nested insertion. The favoured origin
is a **single retrotransposition event** (template switching on a read-through
tandem transcript, or illegitimate LTR recombination) — *not* gene conversion. So
the array is one biological unit and is best annotated as such.

## Proposed detection

A maximal chain of **consecutive same-lineage, same-strand complete elements**
where each adjacent pair's terminal LTRs coincide — element A's rightmost
`long_terminal_repeat` reciprocally overlaps element B's leftmost LTR (≥ ~0.5 of
the shorter LTR). The shared-LTR test distinguishes a tandem (boundary LTR
shared) from a nested insertion (deep overlap, LTRs don't align) and self-limits
array length, so **no arbitrary kb locus cap is needed** (Gypsy elements can be
>25 kb; the example array is 31 kb). DANTE_LTR already has every input it needs
(intact elements, per-element LTR coordinates, domains, lineage).

## Proposed representation

One **container** `transposable_element` spanning the array
(`structure=LTR_RT_TR`, `copy_number=N`, `Final_Classification` = the shared
lineage), with the member copies kept as **children** (`Parent=<container>`),
each retaining its own LTR/domain children. Non-tandem elements unchanged. This
makes the annotation non-overlapping while preserving every copy's detail, and
exposes the tandem structure explicitly for downstream tools.

## Interim workaround (this pipeline)

Until DANTE_LTR supports this natively, the `resolve_ltr_tandems` rule (between
`dante_ltr` and the unified annotation, `scripts/resolve_ltr_tandems.py`) applies
the detection above and writes a **small companion file**,
`DANTE_LTR/DANTE_LTR_tandems.gff3`, holding only the derived containers (each
listing its member element IDs). **`DANTE_LTR.gff3` is left untouched**, so the
LTR library, masking track, repeat report and `dante_line` keep seeing every
individual element exactly as DANTE_LTR emitted them. `make_unified_annotation.R`
reads both files, counts each array once as a container, and nests the member
copies as Level-2 children. On the Boechera region this turns the 4 overlapping
elements into 1 container + 4 nested members (overlap eliminated); a genome with
no tandems yields an empty companion file. Validation on a complete assembly
(copy-number range, cross-lineage edge cases, degenerate `LTR-RT-related`
members) is planned.

## Validation suggested for any DANTE_LTR change

Before/after on a TR/Gypsy-rich reference: confirm (a) the per-array union bp is
unchanged (no base lost, only de-duplicated), (b) member coordinates match the
former per-element calls, (c) genomes without tandems are byte-identical.

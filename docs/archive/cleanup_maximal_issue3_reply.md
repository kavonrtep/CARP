# Reply to issue #3 — "Align cleanup_intermediates: maximal with TideCluster's purge contract"

> Draft reply to `kavonrtep/assembly_repeat_annotation_pipeline#3`, written
> 2026-08-21 against TideCluster **1.21.1** and CARP `e9e8494`. **Not yet
> posted** — the sandbox has no `gh` and the remote is SSH-only, so paste the
> body below as a comment from the host.
>
> Rewritten after TideCluster shipped 1.20.2, 1.21.0 and 1.21.1 the same day: the
> two corrections this draft originally carried (the hardcoded `tc_` prefix, and
> comparative analysis needing `-l`) are **already fixed upstream**, so the reply
> now confirms rather than corrects. Everything below was verified against the
> released 1.21.1 package and a live 94 Gbp run (`run-000170`).

---

Adopted, and thanks — this was a good catch, and 1.20.2/1.21.0/1.21.1 landed
faster than our reply did.

Shipped as option **(a)**, the CARP-side glob set, in `e9e8494`. `maximal` now
deletes the disposable files inside `TideCluster_{tarean,kite,consensus}/`
instead of the three trees, and the patterns are kept character-for-character in
step with `tc_utils.CLEANUP_PATTERNS`.

## Verification on our side

Measured by dry-running the new rule over a live 94 Gbp run (*drVisAlbu1.1*,
5,096 sequences) rather than a fixture:

| | |
|---|---|
| three trees | **44.59 GB** (tarean 39 GB, kite 3.5 GB, consensus 134 MB) |
| freed by the file-level set | **40.23 GB across 3,566 items = 90.2 %** |
| kept | 4.36 GB, and all three capabilities |

Each claim checked independently before adopting it:

| Item | Size | Check |
|---|---|---|
| `monomers.RData` | 25.02 GB | `save()`d at `tarean/methods.R:654-655`; **no `load()` anywhere** |
| `ggmin.RData` | 2.76 GB | same |
| `*.kmers` | 8.43 GB | written by `kmer_counting.py`, never read back |
| `kitehor.periodogram` | 3.46 GB | input to `kite_heatmaps.R`, already rendered |
| per-TRC `TRC_*.fasta` | 0.55 GB | byte-identical to `tarean/fasta/` (confirmed with `cmp`) |

Your 800 Mbp figure was 82 %; at genome scale it is **90 %**, because the
write-only `.RData` pair alone is 62 % of the trees. The bigger the assembly, the
better the trade — 4.4 GB is a rounding error beside a 95 GB genome. Good to see
that number already in the 1.21.0 notes.

We hit the `_tarean/*/TRC_*.fasta` trap you documented in
`CLEANUP_PROTECTED_DIRS` while writing our set, independently and for the same
reason. Our test asserts it: it fails with
`maximal deleted a TideCluster capability file: .../fasta/TRC_1.fasta` if the
glob is loosened.

## One possible omission in `CLEANUP_PATTERNS`

The long-period re-search leaves three files. `--cleanup` removes
`_kite_input_longext.fasta` and `_rescored_longext.*` but keeps
`_kite_peaks_longext.tsv`, which is written and consumed inside the same
function (`tc_utils.py:3458`) and looks equally disposable. It is 146 KB, so this
is tidiness rather than disk — flagging it in case it was an oversight. We match
your set exactly, including keeping that file, rather than using a broader
`_*longext*`: being wider than the owner's contract on paths we do not own is how
this breaks silently later.

## Your open question 2 — yes, there was one

`make_repeat_report.R:646` and `make_summary_plots.R:192` both read

```
TideCluster/<run>/TideCluster_kite/monomer_size_top3_estimats.csv
```

to label a tandem family with its monomer length (`TRC_7 (172 bp)`). Under the
old tree-level purge this degraded silently — cleanup runs after the report, so a
single run was fine, but re-rendering from an archived run quietly dropped every
`(bp)` label. Your `CLEANUP_PATTERNS` keeps it, so nothing is needed on your
side; recording it so it stays in the keep set.

Nothing else: no CARP manifest output lives inside the three trees, and neither
`TideCluster_clustering.gff3_1.gff3` nor the kite `_*longext*` files are
referenced anywhere in CARP.

## On the comparative-analysis argument — now moot, and thank you

Our draft was going to push back on this one. When the issue was filed it did not
apply to CARP, for two reasons that had nothing to do with cleanup:

- three of the six inputs are written only by the annotation step, which
  `run_all` skips without `-l`, and CARP passes `-l` only when the optional
  `tandem_repeat_library` key is set — so on our runs they never existed;
- `get_seq_files()` hardcoded a `tc_` prefix (with a `# TODO`), so it could not
  read a `TideCluster_`-prefixed directory even when they did.

**1.20.2 fixed the second and 1.21.1 fixed the first** — moving
`consensus_sequences_all.fasta` into the clustering step, plus the rebuild
fallback for archived runs, is a better fix than the one we were going to ask
for, because it makes every existing CARP run comparative-ready rather than
requiring a config change and a re-run. Most of our runs do not set
`tandem_repeat_library`, so that is the difference between "no CARP run can be
compared" and "all of them can". Worth calling out in the release notes for
downstream users, if it is not already.

## Your open question 3 — no

Dropping `tarean/fasta/` as well would add **0.61 GB of 44.59 GB** on our run —
not worth losing per-array consensus for 1.4 %.

## On (a) vs (b)

We went with (a) because `--cleanup` did not exist when we implemented it, but we
expect to **keep doing this CARP-side now that it does**. CARP's cleanup is one
post-run step behind one config key (`cleanup_intermediates: minimal|maximal|none`,
with `--keep-all` overriding it) that also covers DANTE / DANTE_LTR / DANTE_TIR /
RepeatMasker / mmseqs scratch in the same pass; splitting the TideCluster part
onto a flag on one tool would make the contract harder to explain, not easier.

Your point about our set being a list of paths we do not own is the real risk, and
we have taken it seriously: the globs are documented one reason each, they name
`tc_utils.CLEANUP_PATTERNS` as the upstream definition so a divergence shows up as
a plain diff, and a test fails loudly if the layout moves under us.

That also settles your "please avoid" note — we will not be passing `--cleanup`,
so both purges can never be applied to the same directory.

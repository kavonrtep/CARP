# Reply to issue #3 — "Align cleanup_intermediates: maximal with TideCluster's purge contract"

> Draft reply to `kavonrtep/assembly_repeat_annotation_pipeline#3`, written
> 2026-08-21 against TideCluster 1.20.1 and CARP `e9e8494`. **Not yet posted** —
> the sandbox has no `gh` and the remote is SSH-only, so paste the body below as
> a comment from the host. Everything in it was verified against the code and
> against a live 94 Gbp run (`run-000170`); the numbers are measurements, not
> estimates.

---

Adopted — thanks, this was a good catch and the analysis held up everywhere we
could check it.

Shipped as option **(a)**, the CARP-side glob set, in `e9e8494`. `maximal` now
deletes the disposable files inside `TideCluster_{tarean,kite,consensus}/`
instead of the three trees.

## Verification on our side

Measured on a live 94 Gbp run (*drVisAlbu1.1*, 5,096 sequences), dry-running the
new rule over the real tree rather than a fixture:

| | |
|---|---|
| three trees | **44.59 GB** (tarean 39 GB, kite 3.5 GB, consensus 134 MB) |
| freed by the file-level set | **40.23 GB across 3,567 items = 90.2 %** |
| kept | 4.36 GB, and all three capabilities |

Breakdown, with each claim checked:

| Item | Size | Check |
|---|---|---|
| `monomers.RData` | 25.02 GB | `save()`d at `tarean/methods.R:654-655`; **no `load()` of it anywhere** in TideCluster |
| `ggmin.RData` | 2.76 GB | same |
| `*.kmers` | 8.43 GB | written by `kmer_counting.py`, never read back |
| `kitehor.periodogram` | 3.46 GB | input to `kite_heatmaps.R`, already rendered |
| per-array `TRC_*.fasta` | 0.55 GB | byte-identical to `tarean/fasta/TRC_*.fasta` (confirmed with `cmp`) |

Your tomato run gave 82.3 %; at genome scale it is **90 %**, because the
write-only `.RData` pair alone is 62 % of the trees. The larger the assembly, the
better the trade — 4.4 GB is a rounding error next to a 95 GB genome.

## One warning if you implement `--cleanup` with globs

TideCluster stores each array's FASTA **twice**:

```
TideCluster_tarean/TRC_1.fasta_tarean/TRC_1.fasta   <- disposable copy
TideCluster_tarean/fasta/TRC_1.fasta                <- needed by tc_per_tra_consensus.py
```

The natural pattern `TideCluster_tarean/*/TRC_*.fasta` matches **both** — it
would delete the file the whole change exists to preserve. We had to write
`*_tarean/TRC_*.fasta`, and our test asserts it (it fails with
`maximal deleted a TideCluster capability file: .../fasta/TRC_1.fasta` if the
glob is loosened). Worth a test on your side too.

## Your open question 2 — yes, there is one

`make_repeat_report.R:646` and `make_summary_plots.R:192` both read

```
TideCluster/<run>/TideCluster_kite/monomer_size_top3_estimats.csv
```

to label a tandem family with its monomer length (`TRC_7 (172 bp)`). Under the
old tree-level purge this degraded silently: cleanup runs after the report, so a
single run was fine, but re-rendering a report from an archived run quietly
dropped every `(bp)` label. Your proposed set already keeps this file, so (a)
fixed it as a side effect — please keep it in the `--cleanup` keep set too.

Nothing else: no CARP manifest output lives inside the three trees, and neither
`TideCluster_clustering.gff3_1.gff3` nor the kite `_*longext*` files are
referenced anywhere in CARP.

## Correction: argument 1 (comparative analysis) does not apply to CARP

The other two capability losses are real and are why we took the change. The
comparative-analysis one is not — for two reasons, neither of them cleanup.

**a) Three of the six inputs are never written on a CARP run.**
`TideCluster_annotation.gff3`, `TideCluster_annotation.tsv` and
`TideCluster_consensus/consensus_sequences_all.fasta` all come from the
annotation step, which `run_all` skips unless a library is supplied:

```python
if cmd_args.library:
    annotation(cmd_args.prefix, cmd_args.library, cpu=cmd_args.cpu)
```

CARP passes `-l` only when the optional `tandem_repeat_library` config key is
set. On the 94 Gbp run above all three are absent right now, mid-run, long
before any cleanup could touch them.

**b) Even with a library, `tc_comparative_analysis.R` cannot read a CARP tree.**
`get_seq_files()` hardcodes a `tc_` prefix:

```r
get_seq_files <- function(input_dirs, prefix){
  # TODO - handle tc_ prefixes
  tarean_consensus_path <- "tc_consensus_dimer_library.fasta"
  consensus_groups_path <- "tc_consensus/consensus_sequences_all.fasta"
```

`prefix` is passed in and ignored, while lines 1361-1363 of the same script *do*
build their paths from it. CARP names everything `TideCluster_*`, so those two
`readDNAStringSet` calls fail on any CARP directory regardless of cleanup.

So: had we changed nothing, comparative analysis still would not run on CARP
output. If you want that capability to work for CARP runs, the fix is your
`get_seq_files()` TODO, not our purge set — and it would be worth stating in the
docs that comparative analysis needs a run made **with** a reference library.

## Your open question 3 — no

Dropping `tarean/fasta/` as well would add **0.61 GB of 44.59 GB** on our run.
Not worth giving up per-array consensus for 1.4 %.

## On (a) vs (b)

We went with (a) because `--cleanup` does not exist in 1.20.1 yet. To set
expectations: we expect to **keep doing this CARP-side even once the flag ships**.
CARP's cleanup is one post-run step the user controls through a single config key
(`cleanup_intermediates: minimal|maximal|none`, with `--keep-all` overriding it),
and it covers DANTE / DANTE_LTR / DANTE_TIR / RepeatMasker / mmseqs scratch in the
same pass; splitting the TideCluster part off into a flag on one tool would make
the contract harder to explain, not easier. Your point about our set being a list
of paths we do not own is fair — the mitigation is that the globs are narrow,
documented with a reason each, and covered by a test that fails loudly if the
layout moves under us.

That also settles your "please avoid" note: we will not be passing `--cleanup`,
so there is no risk of both purges being applied.

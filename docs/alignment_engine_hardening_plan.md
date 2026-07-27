# Implementation plan — hardening the shared all-vs-all alignment engine

Addresses the 9 findings in [`carp_fallback_feature_requests.md`](carp_fallback_feature_requests.md).

**Scope note (read first).** These are not fallback-only. `dante_line` and
`dante_tir_fallback` share one engine — `scripts/global_local_aln.py`
(`run_all_vs_all_alignment`) plus two helpers *defined in* `scripts/dante_line.py`
(`analyze_alignment_lengths`, `run_mmseqs_clustering`) that the fallback imports.
So **every fix below changes both rules and must be validated on both** (the
medium fixture exercises each). See memory `project_alignment_engine_shared`.

**Invariants to preserve.** (a) Annotation output unchanged — these are
robustness / memory / determinism fixes, not algorithm changes. (b) Don't regress
the 1.1.6 determinism work: outputs stay content-stable. (c) Keep the
importlib-fallback dual-context import pattern intact.

Phases are ordered by the agreed priority. Each is independently shippable.

---

## Phase 1 — Stop silent-wrong-output (items 6, 7) — HIGHEST

Correctness beats everything: today a killed/《failed run can look successful with
missing data.

- **Item 6 — atomic intermediates + honest checkpoints.** Write every derived
  file to `<name>.tmp` and `os.replace()` onto the final name only on success, so
  a partial file can never sit under the final name. Then the `.exists()`-skip
  checkpoints (`dante_line.py:1247`, `dante_tir_fallback.py:383`, and
  `run_mmseqs_clustering`'s checkpoint at `dante_line.py:1125`) only ever accept
  complete outputs. Treat zero-length as absent.
  - Files: `global_local_aln.py` (`write_results_table`), `dante_line.py`
    (`analyze_alignment_lengths` writer, `run_mmseqs_clustering`).
- **Item 7 — surface failures instead of degrading quietly.** Collect per-fasta /
  per-subtype failures (currently `except Exception: print("Error…"); continue`
  at `dante_line.py:1234/1256`, `dante_tir_fallback.py:404/885`) and print a
  summary block at the end (`WARNING: N of M inputs produced no alignments: …`);
  write a marker file (or non-zero exit) when a rule yields nothing.
  - **Rule-level (dante_line only, important):** `Snakefile:770` turns *any*
    dante_line non-zero exit into "empty LINE outputs". Split the legitimate
    "too few LINE features" path (expected, → empty outputs OK) from an
    unexpected crash (should fail loudly). Simplest: have `dante_line.py` exit
    with a distinct code / write a sentinel for the benign case, and let the
    wrapper create empty outputs only for that code.
- **Tests:** truncate an intermediate then rerun → it is redone, not trusted; a
  forced subtype failure → summary line + marker; the dante_line "too few
  features" path still produces empty outputs without failing the pipeline.
- **Effort:** small. **Risk:** low. **Rules:** both (+ Snakefile for dante_line).

## Phase 2 — Bounded memory for high-copy families (items 1, 4) — SCALING

The two that decide whether a high-copy family (LINE or TIR) fits in RAM on a
30–90 Gbp genome.

- **Item 1 — stream the alignment table.** `run_all_vs_all_alignment` currently
  holds `results` for *all* groups then writes once (`global_local_aln.py:1009/
  1031/1044`). Open the output once, write each group's records as that group
  finishes, keep only one group in memory. Single-group path writes directly.
  (Row order handled in Phase 3.)
- **Item 4 — bounded length tracking.** `analyze_alignment_lengths`
  (`dante_line.py:1021`) stores every length per group then sorts. Replace with a
  per-group **min-heap of size N** (`heapq.heappushpop`) for `Selected_Length`,
  plus running counters for `Num_Shorter` (it needs `total_count` and
  `count_strictly_greater_than_selected` — both derivable from counters, no stored
  list). Input is already read line-by-line.
- **Tests:** old-vs-new byte-identical output on a fixture; capture the rule's
  `benchmark:` peak RSS before/after on a high-copy input.
- **Effort:** small–moderate, one function each. **Risk:** low (fixture
  equivalence test). **Rules:** both.

## Phase 3 — Determinism / verifiability (items 3, 8)

Cheap, and they restore the checksum-the-intermediate check we leaned on in 1.1.6.

- **Item 3 — stable alignment-TSV row order.** Rows are appended in
  `as_completed` order (`global_local_aln.py:907`) → non-deterministic file.
  Collect results in submission order (index-and-sort, or `executor.map(...,
  chunksize=…)` which preserves input order) and, combined with Phase-1 streaming,
  emit each group in a fixed order. Result: byte-stable intermediate. (Downstream
  was already stable via the sort at `dante_line.py:1064`; this makes the
  intermediate itself checksummable.)
- **Item 8 — guarantee mmseqs input order.** mmseqs clustering is order-sensitive
  (measured: shuffling the same seqs moved 258/5,983 clusters). The FASTA order
  through the BED/seqkit path looks deterministic but is **not asserted**, and a
  `set` in the chain (e.g. `dante_line.py:207`) would make it vary per process
  via `PYTHONHASHSEED`. Audit `group_sequences_for_alignment` / the
  `LINE_regions_extended.fasta` build for `set`→ordered fixes, and add a test
  that builds the extracted FASTA twice under different `PYTHONHASHSEED` and
  asserts byte-equality (mirror DANTE_TIR's `tests/test_aa_fasta_order.py`; needs
  no data).
- **Tests:** two-run checksum of the alignment TSV; the hash-seed FASTA test.
- **Effort:** small. **Risk:** none to results. **Rules:** both.

## Phase 4 — Thread-scaling measurement (item 5) — MEASURE FIRST

- Run one large family through `run_all_vs_all_alignment` at `-t 4` vs `-t 32`;
  compare wall time and CPU-minutes. **Nuance:** unlike the blastn case that
  motivated this finding, parasail's Python binding is **ctypes** (releases the
  GIL during the C alignment), so threads should help the compute; the cap is the
  per-pair Python overhead (result-dict build + `score_threshold` filter at
  `global_local_aln.py:885/902`), which for very short LINE-end flanks can
  dominate. So do not assume the blastn numbers transfer.
- **If wall time barely moves:** switch the inner loop to process-level
  concurrency (several groups at once, or a `ProcessPoolExecutor`) rather than
  more threads per group.
- **Effort:** measurement, then possibly moderate. **Risk:** n/a until measured.

## Phase 5 — Housekeeping (items 2, 9)

- **Item 2 — bound in-flight futures.** `_compare_sequences` submits all pairs
  up front (`global_local_aln.py:893-905`). Submit in a bounded window / use
  `executor.map(chunksize=…)` so memory scales with the window. Largely subsumed
  by Phase 3's ordering rework — do them together.
- **Item 9 — drop `tempfile.mktemp`.** Replace at `global_local_aln.py:552` and
  `:1026` with `mkstemp()` / `NamedTemporaryFile(delete=False)`.
- **Effort:** trivial. **Risk:** none.

---

## Cross-cutting

- **Validate on both rules every time.** Run the medium fixture (exercises
  dante_line and the fallback) and assert unchanged annotation outputs after each
  phase; `scripts/assert_fixture_outputs.sh medium` + a diff of the split GFF3s.
- **Confirmation experiment (validates Phases 1–2 & 4).** A CARP run on a genome
  with a genuinely high-copy LINE family (and/or a high-copy TIR family) —
  capture the rule `benchmark:` peak RSS and `len(results)` per group. This is
  the missing measurement the feature-request doc flags; nothing here has been
  reproduced on CARP at scale yet.
- **Suggested PR grouping:** (A) Phase 1 + Phase 3 + item 9 — low-risk
  correctness/determinism, ship first; (B) Phase 2 + item 2 — memory scaling,
  gated on the fixture-equivalence test + an RSS measurement; (C) Phase 4 —
  measurement, separate.
- **Docs/changelog:** user-visible only where behaviour changes (Phase 1's
  loud-failure / marker + the dante_line exit-code split) — add a CHANGELOG
  bullet for those; the rest are internal.

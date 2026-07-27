# Feature requests for CARP's `dante_tir_fallback`

Findings from reading `kavonrtep/CARP` (cloned at CHANGELOG 1.1.6) after fixing
the equivalent problems in DANTE_TIR's Round 3. The fallback rule solves the
same shape of problem — an all-vs-all comparison over the flanks of every
copy of a high-copy family — so several of the failure modes we measured on
run-000129 apply here too.

Scope: `scripts/dante_tir_fallback.py`, plus the shared machinery it calls in
`scripts/global_local_aln.py` and `scripts/dante_line.py`. Nothing here has
been reproduced *on CARP* — these are code-reading findings backed by measured
experience in the sister pipeline, so each item says what would confirm it.

## What CARP already gets right

Worth stating, because it means the remaining items are narrow:

- **Bounded, cluster-aware grouping.** `max_group_size` (default 1000) with
  mmseqs-based grouping, and a deterministic id-sorted fallback when mmseqs is
  unavailable, is exactly the design DANTE_TIR needed and lacked. The comment
  at `global_local_aln.py:1013` shows the O(N²) blow-up was understood.
- **Streaming sequence extraction.** seqkit-based extraction never loads the
  genome. DANTE_TIR had two `readDNAStringSet(genome)` calls that asked for
  ~94 GB on a 94.3 Gbp assembly; CARP does not have that class of bug here.
- **Correct subprocess usage.** Every call in the fallback path uses
  `subprocess.run(...)`, which drains its pipes. DANTE_TIR had
  `check_call(..., stdout=PIPE, stderr=PIPE)`, which deadlocks once the child
  writes past the 64 KB pipe buffer — that cost a real run. CARP is clean;
  `canonical_sort_fasta.py:103` uses `check_call` but redirects to a file.
- **Measured, not assumed, parallelism.** The comment at
  `dante_tir_fallback.py:1070` documents a parallelisation that was benchmarked
  (2:11 → 7:48 wall) and reverted. That is the right instinct, and item 5 below
  suggests extending the same scepticism one level down.

---

## 1. The alignment table is materialised in RAM before anything is written

**Where:** `global_local_aln.py:1009–1044` (`run_all_vs_all_alignment`).

`results = []` accumulates every surviving alignment record across *all*
groups, and `write_results_table(results, output_file)` writes once at the end.
Peak memory therefore scales with the total number of records in the whole
run, not with the largest group.

**Why it matters:** this is precisely the pattern that made DANTE_TIR's Round 3
unrunnable. There, the equivalent table reached 2.54 × 10⁹ rows (41.5 GB *after*
filtering) on a 168,012-copy EnSpm/CACTA family, and the run died. The grouping
in CARP bounds the *per-group* work but not the accumulation across groups: a
family of 168k copies at `max_group_size=1000` is 168 groups, all of whose
records are held simultaneously as Python dicts (roughly 300–500 bytes each).

**Suggested fix:** open the output once and append each group's records as that
group finishes, keeping only one group in memory. The single-group path can
write directly too. Output is unchanged except for row ordering, which item 3
addresses anyway.

**Effort:** small, contained in one function. **Risk:** low.

**Confirm it first:** log `len(results)` per group on the largest family in a
real run; multiply by ~400 bytes.

## 2. Every pair of a group is held as a live `Future`

**Where:** `global_local_aln.py:892–914`.

All pairs are submitted up front (`futures.append(future)` for the entire
`pair_generator`) before any result is consumed. At `max_group_size=1000` that
is up to ~500k `Future` objects alive at once, each retaining its result until
`as_completed` reaches it.

**Suggested fix:** bound the in-flight window — submit in chunks and drain, or
use `executor.map(..., chunksize=...)`. Memory then scales with the window, not
the group.

**Effort:** small. **Risk:** low.

## 3. The alignment TSV row order is non-deterministic

**Where:** `global_local_aln.py:907` (`for future in as_completed(futures)`).

Rows are appended in completion order, so the file differs between runs of
identical input. The derived `*_aln_length.tsv` sorts per group, so the
*downstream* result appears order-independent — but the intermediate cannot be
checksum-compared, which removes the cheapest way to verify a rerun.

**Why we care:** in DANTE_TIR we leaned on byte-identical intermediates
repeatedly to prove that a refactor changed nothing. Losing that makes every
future change harder to validate.

**Suggested fix:** preserve submission order (collect into a preallocated list
indexed by submission, or sort by `(query_id, ref_id)` before writing).

**Effort:** small. **Risk:** none to results.

## 4. `analyze_alignment_lengths` keeps every alignment length, but needs only the Nth largest

**Where:** `dante_line.py:1021–1092`.

`group_lengths[query_id].append(...)` and the same for `ref_id`, for every row
of the alignment TSV — then the function sorts each list and takes element
`min_num_alignments - 1` (default **3**).

So the peak is two Python ints per alignment record for the whole file, when
the answer needs at most N per group.

**Suggested fix:** a bounded min-heap of size N per group
(`heapq.heappushpop`). Identical output, memory proportional to
groups × N instead of records. `Num_Shorter` needs a counter alongside, which
is a running count rather than a stored list.

**Effort:** small, ~15 lines. **Risk:** low; easily covered by a test comparing
old and new on a fixture.

## 5. Thread scaling of the inner all-vs-all is worth measuring

**Where:** `global_local_aln.py:892` (`ThreadPoolExecutor(max_workers=threads)`).

Python threads only help if the alignment work releases the GIL. If parasail is
called through a binding that holds it, `-t 96` buys nothing and the rule is
effectively serial inside each group.

**Why we raise it:** we assumed blastn scaled with `-num_threads` and it does
not — measured 24.7 effective cores of 96 on run-000129, 6.1 of 14 locally, and
efficiency *rises* as threads fall (61 % at 7). Splitting the work into several
smaller concurrent searches gave **3.0× on fewer total threads**. The same
question — "does adding threads actually add throughput here?" — has not been
answered for this loop.

**Suggested check:** run one large subtype at `-t 4` and `-t 32` and compare
CPU-minutes and wall time. If wall time barely moves, the fix is process-level
concurrency (several groups at once) rather than more threads per group.

**Effort:** measurement first, then possibly moderate. **Risk:** n/a.

## 6. Cached intermediates are trusted without an integrity check

**Where:** `dante_tir_fallback.py:383` (`if output_path.exists(): continue`),
and the same idiom for `*_aln_length.tsv`.

A run killed midway leaves a truncated TSV that the next run accepts as
finished.

**Why we care:** DANTE_TIR had the identical idiom around CAP3's output. When
CAP3 segfaulted it left a zero-byte `.cap.aln`, the pipeline treated it as a
completed assembly, and an entire superfamily silently lost its primary
detection round — on a real run, unnoticed until we went looking months later.

**Suggested fix:** write to `<name>.tmp` and rename on success, so a partial
file cannot exist under the final name. Optionally treat zero-length outputs as
absent.

**Effort:** small. **Risk:** none.

## 7. Failures degrade the result with only a log line

**Where:** `dante_tir_fallback.py:404` (`Error analyzing ...`), `:885`
(`Warning: mmseqs clustering failed for {subtype}`).

Both continue with reduced output. In a long log on a multi-day run, a single
warning is easy to miss, and the result looks like a successful run with fewer
elements.

**Suggested fix:** collect failures and print a summary block at the end
(`WARNING: N of M subtypes produced no alignments: ...`), and consider a
non-zero exit or a marker file when a subtype yields nothing. DANTE_TIR now
does this for CAP3.

**Effort:** small. **Risk:** none.

## 8. mmseqs2 clustering is order-sensitive — is the input order guaranteed?

**Where:** `group_sequences_for_alignment` and `run_mmseqs_prefilter` in
`global_local_aln.py`; `run_mmseqs_clustering` in `dante_line.py`.

**Measured** on run-000129's MuDR domains (11,260 sequences, mmseqs
`easy-cluster`, `--min-seq-id 0.8 -c 0.8 -s 7.0`):

| test | result |
|---|---|
| same input, same threads, rerun | identical |
| same input, 1 thread vs 4 | identical |
| **same sequences, shuffled order** | **differs — 258 of 5,983 clusters changed** |

So grouping is reproducible only as long as the FASTA handed to mmseqs is
byte-stable. In CARP that order comes from the BED/GFF3 path through seqkit,
which looks deterministic, but it is not asserted anywhere. The hazard is a
`set` entering that chain later: Python randomises string hashing per process,
so set iteration order varies between runs of the same code, and the symptom —
clusters differing between identical runs — is not something anyone thinks to
check.

**Suggested fix:** a unit test that builds the extracted FASTA twice in
separate interpreters under different `PYTHONHASHSEED` values and asserts the
bytes match. DANTE_TIR added `tests/test_aa_fasta_order.py` for exactly this;
it needs no data and runs instantly.

**Effort:** small. **Risk:** none.

## 9. `tempfile.mktemp()` is deprecated and racy

**Where:** `global_local_aln.py:552` and `:1026`.

`mktemp` returns a name, not an open handle, so it is a documented race. Use
`tempfile.mkstemp()` or `NamedTemporaryFile(delete=False)`.

**Effort:** trivial. **Risk:** none.

---

## Suggested order

1. Item 6 (atomic writes) and item 7 (failure summary) — smallest, and they
   protect against silently wrong output, which is the worst failure mode.
2. Item 1 and item 4 (streaming writes, bounded length tracking) — the two that
   decide whether a very high-copy family fits in memory.
3. Item 8 (determinism test) and item 3 (stable row order) — cheap, and they
   make everything afterwards verifiable.
4. Item 5 (thread scaling measurement) — potentially the largest speed win, but
   measure before changing anything.
5. Items 2 and 9 — housekeeping.

## What we would need to confirm any of this on CARP itself

A run of `dante_tir_fallback` on a genome with a genuinely high-copy TIR family
— run-000129's `drVisAlbu1.1` (168,012 EnSpm/CACTA copies) is the obvious
candidate, and its `DANTE_TIR_FALLBACK/` directory is empty, so the rule has
not been exercised at that scale there. Peak RSS from the rule's `benchmark:`
entry plus `len(results)` per group would settle items 1, 2 and 4 immediately.

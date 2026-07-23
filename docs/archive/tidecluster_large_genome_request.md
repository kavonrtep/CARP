# TideCluster: large-genome fd exhaustion + whole-genome memory load

> Feature request for `github.com/kavonrtep/TideCluster` — audited on **1.16.3**.
> Two issues each break a ~90 Gbp run outright; both fixes are low-risk.

## Summary

1. **fd exhaustion (same class as the DANTE_LTR bug):**
   `split_fasta_to_chunk_files` opens **one handle per chunk file, all at once**
   (`tc_utils.py:1040`). At `--chunk_size 50000000` a 90 Gbp genome makes ~1800
   chunks ⇒ ~1800 concurrent handles vs. the default `ulimit -n` of 1024 ⇒
   `OSError: [Errno 24] Too many open files`. This is on the
   `tc_reannotate.py --ref_seq` → `run_repeatmasker_genome_chunked` path
   (`tc_reannotate.py:195` → `tc_utils.py:1163`), so it is deterministic for any
   genome larger than ~50 Gbp (or smaller with a reduced `--chunk_size`).
2. **memory (whole genome loaded, then discarded):** the `clustering` stage
   loads the entire reference FASTA into a Python dict (~90 GB for 90 Gbp) via
   `gff3_to_fasta` → `fasta_to_dict` (`TideCluster.py:623` → `tc_utils.py:1538`),
   but the extracted genomic sequence **is never used** — the loop consumes only
   the `consensus_sequence` GFF attribute. Likely OOM, for no benefit.

## Finding A — fd exhaustion (`tc_utils.py:1040`)

```python
# tc_utils.py:1040, inside split_fasta_to_chunk_files (chunk_size default 50_000_000)
open_handles = {p: open(p, "w") for p in file_paths}     # <-- all handles open at once
try:
    with open(fasta_file) as fh:
        for header, sequence in read_single_fasta_as_generator(fh):
            for _orig, _i, start, end, tok in pieces_by_header.get(header, []):
                out_fh = open_handles[token_to_file[tok]]
                ...
finally:
    for fh in open_handles.values():
        fh.close()
```

`len(file_paths) ≈ genome_bp / chunk_size` (bin-packing at
`tc_utils.py:1028-1034` places ~1 piece per file, since each large-sequence piece
is already `chunk_size + overlap`). 90 Gbp / 50 Mb ⇒ ~1800 handles, opened
before any RepeatMasker runs. Identical canonical pattern to the DANTE_LTR bug
(`[open(f,'w') for f in temp_files]`).

**Call path (the pipeline uses exactly this):**
`tc_reannotate.py:195 run_repeatmasker_genome_chunked(...)` →
`tc_utils.py:1163 split_fasta_to_chunk_files(...)`. CARP's `tidecluster_reannotate`
rule runs `tc_reannotate.py --ref_seq`, so this fires on any large-genome run.

**Proposed fix.** Bound the number of concurrently open handles. Simplest:
an LRU of ≤256 append-mode handles — on a cache miss `open(path, "a")`, and when
over the cap `close()` the least-recently-used (correctness preserved because
every write is append-only). Alternatively raise `RLIMIT_NOFILE` toward the hard
limit at start-up (`resource.setrlimit`) as a defensive complement, but the
handle cap is the real fix (it holds even when the chunk count exceeds the hard
limit). Add a regression test to `tests/test_chunked_repeatmasker.py` with a
synthetic input that forces > 1024 chunk files.

## Finding C1 — clustering loads the whole genome, then ignores it (`TideCluster.py:623`)

```python
# TideCluster.py:623  — `seq` is bound but never read in the loop body
for seq_id, seq, cons in tc.gff3_to_fasta(gff3, fasta, "consensus_sequence"):
    mult = round(1 + 10000 / len(cons))
    consensus = cons * mult          # uses only `cons` (the consensus_sequence attribute)
    consensus_dimers = cons * 4
    ...
    f.write(F">{seq_id}\n{consensus}\n")     # and only `seq_id`
```

`gff3_to_fasta` (`tc_utils.py:1538`) begins with `fasta_dict = fasta_to_dict(fasta_file)`,
slurping the entire reference FASTA into `{name: str}` — ~90 GB for a 90 Gbp
assembly (more with object overhead) — and slices a genomic sub-sequence per
feature. Here that genomic sequence (`seq`) is discarded; only the
`consensus_sequence` attribute is needed.

**Proposed fix (near-zero risk).** Don't load the genome for this call. Either
add a lightweight generator that yields `(ID, consensus_sequence_attr)` straight
from the GFF3 (no FASTA read at all), or give `gff3_to_fasta` a flag to skip
`fasta_to_dict`/`get_seq_from_fasta` when only the attribute is requested. This
removes the ~90 GB load that OOMs every large clustering run.

## Secondary findings (scaling, lower priority)

- **C2 — rDNA genomic fallback reloads the whole genome.** `identify_rdna` →
  `_write_trc_region_subject` (`tc_utils.py:691`) calls
  `gff3_to_fasta(tmp_gff, fasta, "Name")` → full `fasta_to_dict(genome)`, firing
  whenever a TRC lacks a usable consensus. Replace with a streaming extractor
  that indexes only the seqids in `tmp_gff` and reads the FASTA once via
  `read_single_fasta_as_generator`.
- **B1 — `resolve_trc_overlaps` is O(F²) per seqid and single-threaded**
  (`tc_utils.py:2337`): the per-segment `covering = [... for ... in fl if ...]`
  re-scans all features on the seqid. Replace with an active-set sweep
  (O(F log F)); parallelise the per-seqid loop. Dominant serial cost after the
  parallel mmseqs/blast search on satellite-rich chromosomes.
- **B2 — chunked-RM serial reduce tail** (`tc_utils.py:1210-1241` +
  `tc_reannotate.py:203-242`): all `.out` hits parsed and accumulated into one
  `records` list on a single core, then sorted/merged. Move `.out` parsing into
  `_repeatmasker_chunk_worker` (return parsed records) and merge per-seqid. Also
  caps the C3 in-RAM `records` growth.
- **B4 — `split_gff3_by_cluster_name` open/close per feature**
  (`tc_utils.py:351-362`): millions of `open`+`close` syscalls on a large genome.
  Bucket lines by `Name` in memory (bounded by TRC count), write each file once.

## Priority

**Finding A (fd) and Finding C1 (memory)** each break a 90 Gbp run and should
land first — A crashes `tc_reannotate`, C1 OOMs `clustering`. B1/B2/C2/C3 are the
wall-time/RAM scaling improvements for the satellite-dense case.

## Pipeline-side stopgap (CARP)

Because Finding A is on the `tc_reannotate` path, the CARP `tidecluster_reannotate`
rule needs the **same `ulimit -n` bump** as the `dante_ltr` rule until a fixed
TideCluster release ships:

```bash
ulimit -n "$(ulimit -Hn)"     # tc_reannotate opens ~1800 chunk handles at 90 Gbp / 50 Mb
```

The memory issue (C1) has no pipeline-side mitigation — it needs the upstream fix
(or a much smaller genome). Raising `--chunk_size` reduces the fd handle count
(fewer chunks) as a partial stopgap for Finding A.

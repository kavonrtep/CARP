# DANTE_LTR: `Too many open files` on large / fragmented genomes

> Feature request for `github.com/kavonrtep/dante_ltr` — observed on **0.4.0.4**.

## Summary

On a large genome (~90 Gbp), `dante_ltr` aborts at start-up with
`OSError: [Errno 24] Too many open files`. The chunk-splitting step opens **one
file handle per chunk and holds them all open simultaneously**; the number of
chunks is `total_size / max_chunk_size`, so it grows without bound with genome
size and exceeds the default `ulimit -n` (1024) on any genome larger than
~`1024 × max_chunk_size` (≈ 51 Gbp at the pipeline's `-S 50000000`, less once
other fds are counted). It is independent of assembly quality — a contiguous or
a fragmented 90 Gbp assembly both trip it.

## Evidence (DANTE_LTR 0.4.0.4, ~90 Gbp genome, `-c 32 -M 1 -S 50000000`)

```
running analysis on chunks of ~ 50000000 bp
Traceback (most recent call last):
  File ".../bin/dante_ltr", line 700, in <module>
    main()
  File ".../bin/dante_ltr", line 565, in main
    file_handles = [open(temp_file, 'w') for temp_file in temp_files_fasta]
  File ".../bin/dante_ltr", line 565, in <listcomp>
    file_handles = [open(temp_file, 'w') for temp_file in temp_files_fasta]
OSError: [Errno 24] Too many open files: '.../scratch/tmp/tmpulsinpc1'
```

At `-S 50000000`, 90 Gbp ⇒ `number_of_temp_files = int(90e9 / 50e6) + 1 = 1801`
handles opened at once — well over the default soft limit of 1024.

## Root cause (`dante_ltr`, `main()`)

Two spots open a handle per chunk in one list comprehension and keep them open
to distribute sequences round-robin in a single streaming pass:

```python
# line 560
number_of_temp_files = int(total_size / args.max_chunk_size) + 1
if number_of_temp_files > number_of_sequences:
    number_of_temp_files = number_of_sequences
...
# line 564-567  (FASTA split)
temp_files_fasta = make_temp_files(number_of_temp_files)
file_handles = [open(temp_file, 'w') for temp_file in temp_files_fasta]   # <-- fd blow-up
seq_id_file_handle_dict = dict(zip(seq_id_size_sorted, cycle(file_handles)))
...
# line 584-587  (GFF split — identical pattern, same bug)
temp_files_gff = make_temp_files(number_of_temp_files)
file_handles = [open(temp_file, 'w') for temp_file in temp_files_gff]
seq_id_file_handle_dict = dict(zip(seq_id_size_sorted, cycle(file_handles)))
```

`number_of_temp_files` is `total_size / max_chunk_size` (capped only by the
sequence count), so **concurrent open handles scale linearly with genome size**.
Any genome over ~1024 chunks fails.

The chunks are consumed **serially** afterwards
(`for i in range(number_of_temp_files): subprocess.check_call([... detect_putative_ltr.R ...])`,
line 603), and `max_chunk_size` is the per-chunk **memory** control for the R
detector — so the chunk count cannot simply be capped without changing the
memory profile the user tuned `-S` for.

## The fix already exists in the sibling tool `dante`

`kavonrtep/dante` had this exact bug and fixed it — it now distributes with
**open-append per sequence** instead of holding all handles open. From
`dante` (the `dante` wrapper, ~line 544):

```python
temp_files_fasta = make_temp_files(number_of_temp_files)
# do not use file handles as it can cause problem with maximum number of open files
seq_id_file_path_dict = dict(zip(seq_id_size_sorted, cycle(temp_files_fasta)))
with open(fasta_file_split, 'r') as f:
    for line in f:
        if line[0] == '>':
            header = line.strip().split(' ')[0][1:]
        with open(seq_id_file_path_dict[header], 'a') as fout:   # open-append, 1 fd at a time
            fout.write(line)
```

The same treatment applied to `dante_ltr`'s two split loops removes the bug
with no change to chunking behaviour or memory profile.

## Proposed fix (layered — recommend both)

**1. Raise the soft fd limit at start-up (cheap, zero behavioural change).**
Before the split, lift `RLIMIT_NOFILE` toward the hard limit:

```python
import resource
soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
need = number_of_temp_files + 64          # + stdio / subprocess headroom
if soft < need:
    resource.setrlimit(resource.RLIMIT_NOFILE, (min(max(need, soft), hard), hard))
```

This alone fixes every case where the chunk count fits under the hard limit
(typically 4096–1048576), which covers the reported 90 Gbp run (1801 < hard).

**2. Don't hold all handles open (robust for any size / low hard limit).**
Adopt `dante`'s open-append pattern in both split loops (FASTA at line 565, GFF
at line 585). This is O(sequences) `open()`/`close()` calls but bounded to one
handle at a time, so it can never exhaust fds regardless of genome size or a
capped hard limit. (If per-sequence open/close proves too slow for
many-contig inputs, a bounded-batch variant — write at most `B ≤ soft_limit-64`
chunks per streaming pass — keeps it to `ceil(N/B)` passes with the same fd
guarantee.)

Layer 1 makes the common case a no-op; layer 2 guarantees correctness on the
extreme tail (300 Gbp+, or hosts pinning a low hard limit).

## Pipeline-side stopgap (until a fixed release ships)

Raise the soft limit to the hard limit in the `dante_ltr` rule shell — this
unblocks the reported case immediately because 1801 < the hard limit:

```bash
ulimit -n "$(ulimit -Hn)"     # or an explicit high value, e.g. 65536
dante_ltr -o ... -s ... -g ... -c {threads} -M 1 -S 50000000
```

This is a stopgap only; the durable fix is upstream (layer 2 guarantees it for
genomes whose chunk count exceeds even the hard limit).

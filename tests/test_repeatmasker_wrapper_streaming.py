#!/usr/bin/env python3
"""Equivalence test for the streaming rewrite of repeatmasker_wrapper.py.

split_fasta_to_chunks and split_fasta_to_files used to slurp the whole genome
into a Python dict (read_single_fasta_to_dictionary) — ~90 GB on a 90 Gbp
assembly, materialised twice. They now stream one record at a time via
iter_fasta_records(). This test proves the streamed output is byte-identical to
the original dict-based output across: a long sequence split into overlapping
chunks, several short sequences packed together, line-wrapped input, and a
sequence long enough to span multiple RM bins.

Run: python3 tests/test_repeatmasker_wrapper_streaming.py
"""
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import repeatmasker_wrapper as rw  # noqa: E402


def make_seq(n, salt=0):
    return "".join("ACGT"[(i * 7 + salt) % 4] for i in range(n))


def write_wrapped_fasta(path, seqs, width=60):
    with open(path, "w") as fh:
        for name, seq in seqs:
            fh.write(">%s some description\n" % name)
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


# --- reference (original dict-based) write passes -----------------------------
def old_write_chunks(fasta_file, matching_table, out_path):
    fasta_dict = rw.read_single_fasta_to_dictionary(open(fasta_file, "r"))
    rows_by_header = {}
    for row in matching_table:
        rows_by_header.setdefault(row[0], []).append(row)
    with open(out_path, "w") as fh_out:
        for header in fasta_dict:
            for header2, i, start, end, new_header in rows_by_header.get(header, []):
                fh_out.write(">" + new_header + "\n")
                fh_out.write(fasta_dict[header][start:end] + "\n")


def old_split_to_files(fasta_file, workdir, chunk_size):
    fasta_dict = rw.read_single_fasta_to_dictionary(open(fasta_file, "r"))
    files = []
    cur_fh = None
    cur_size = 0
    for header, sequence in fasta_dict.items():
        seq_len = len(sequence)
        if cur_fh is not None and cur_size > 0 and cur_size + seq_len > chunk_size:
            cur_fh.close()
            cur_fh = None
        if cur_fh is None:
            path = tempfile.NamedTemporaryFile(delete=False, dir=workdir).name
            files.append(path)
            cur_fh = open(path, "w")
            cur_size = 0
        cur_fh.write(">" + header + "\n")
        cur_fh.write(sequence + "\n")
        cur_size += seq_len
    if cur_fh is not None:
        cur_fh.close()
    return files


def main():
    tmp = tempfile.mkdtemp(prefix="rmw_stream_test_")
    genome = os.path.join(tmp, "genome.fasta")
    # long seq (splits into overlapping chunks), several short, one medium.
    seqs = [
        ("chrA", make_seq(2500, 1)),   # > 2*chunk_size(500) -> multiple chunks
        ("scaf1", make_seq(120, 2)),
        ("scaf2", make_seq(300, 3)),
        ("scaf3", make_seq(900, 4)),   # > 2*chunk_size -> chunks too
        ("tiny", make_seq(30, 5)),
    ]
    write_wrapped_fasta(genome, seqs)

    chunk_size, overlap = 500, 50

    # 1) split_fasta_to_chunks: new (streaming) vs old (dict) write, same matching_table
    new_split, matching_table = rw.split_fasta_to_chunks(
        genome, chunk_size=chunk_size, overlap=overlap, temp_dir=tmp)
    ref_split = os.path.join(tmp, "ref_split.fasta")
    old_write_chunks(genome, matching_table, ref_split)
    assert open(new_split).read() == open(ref_split).read(), \
        "split_fasta_to_chunks: streamed output != dict output"
    print("  split_fasta_to_chunks: streamed == dict (%d chunk records)" % len(matching_table))

    # 2) split_fasta_to_files: new (streaming) vs old (dict) binning
    rm_chunk = 800
    d_new = os.path.join(tmp, "new_bins"); os.makedirs(d_new)
    d_old = os.path.join(tmp, "old_bins"); os.makedirs(d_old)
    new_files = rw.split_fasta_to_files(new_split, d_new, chunk_size=rm_chunk)
    old_files = old_split_to_files(new_split, d_old, chunk_size=rm_chunk)
    assert len(new_files) == len(old_files), \
        "different number of RM bins: %d vs %d" % (len(new_files), len(old_files))
    for nf, of in zip(new_files, old_files):
        assert open(nf).read() == open(of).read(), "bin content differs"
    print("  split_fasta_to_files: streamed == dict (%d bins)" % len(new_files))

    # 3) round-trip sanity: concatenated chunk slices reconstruct each original
    #    sequence's covered span (chunks tile [0, size) with overlap).
    seqmap = dict(seqs)
    chunk_seqs = dict(rw.iter_fasta_records(new_split))
    for header, i, start, end, new_header in matching_table:
        assert chunk_seqs[new_header] == seqmap[header][start:end], \
            "chunk %s slice mismatch" % new_header
    print("  chunk slices match the source genome exactly")
    print("test_repeatmasker_wrapper_streaming: PASSED")


if __name__ == "__main__":
    main()

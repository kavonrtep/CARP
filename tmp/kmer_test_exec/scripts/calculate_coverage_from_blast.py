#!/usr/bin/env python3
"""
Input: BLAST tabular output and subject FASTA file
Output: Coverage (as number of hits) of subject sequences based on BLAST hits, number of hits for each nucleotide position
Output format - tab delimited:
subject_id1    position1    number_of_hits
subject_id1    position2    number_of_hits
...
subject_idN    position1    number_of_hits
..
subject_idN    positionL    number_of_hits


tabular output from blast - outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path


def parse_fasta_lengths(fasta_file):
    """
    Parse FASTA file to get sequence lengths.
    Returns dict: {sequence_id: length}
    """
    seq_lengths = {}
    current_seq = None
    current_length = 0

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq is not None:
                    seq_lengths[current_seq] = current_length
                current_seq = line[1:].split()[0]  # Take only the first part of header
                current_length = 0
            elif current_seq is not None:
                current_length += len(line)

    if current_seq is not None:
        seq_lengths[current_seq] = current_length

    return seq_lengths


def parse_blast_output(blast_file):
    """
    Parse BLAST tabular output (outfmt 6).
    Returns list of hits with subject coordinates.
    """
    hits = []

    with open(blast_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue

            fields = line.split('\t')
            if len(fields) < 12:
                continue

            qseqid = fields[0]
            sseqid = fields[1]
            pident = float(fields[2])
            length = int(fields[3])
            mismatch = int(fields[4])
            gapopen = int(fields[5])
            qstart = int(fields[6])
            qend = int(fields[7])
            sstart = int(fields[8])
            send = int(fields[9])
            evalue = float(fields[10])
            bitscore = float(fields[11])

            # Ensure sstart <= send for consistent processing
            if sstart > send:
                sstart, send = send, sstart

            hits.append({
                'qseqid': qseqid,
                'sseqid': sseqid,
                'sstart': sstart,
                'send': send,
                'length': length,
                'pident': pident,
                'evalue': evalue,
                'bitscore': bitscore
            })

    return hits


def calculate_position_coverage(hits, seq_lengths):
    """
    Calculate coverage for each position of each subject sequence.
    Returns dict: {subject_id: {position: hit_count}}
    """
    coverage = defaultdict(lambda: defaultdict(int))

    for hit in hits:
        sseqid = hit['sseqid']
        sstart = hit['sstart']
        send = hit['send']

        # Only process if sequence is in our FASTA file
        if sseqid not in seq_lengths:
            continue

        # Count hits for each position in the alignment
        for pos in range(sstart, send + 1):
            if 1 <= pos <= seq_lengths[sseqid]:  # Ensure position is within sequence bounds
                coverage[sseqid][pos] += 1

    return coverage


def write_coverage_output(coverage, seq_lengths, output_file):
    """
    Write position-wise coverage to output file.
    """
    with open(output_file, 'w') as f:
        for sseqid in sorted(coverage.keys()):
            seq_len = seq_lengths[sseqid]

            # Write coverage for each position (1-based coordinates)
            for pos in range(1, seq_len + 1):
                hit_count = coverage[sseqid].get(pos, 0)
                f.write(f"{sseqid}\t{pos}\t{hit_count}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate position-wise coverage from BLAST output",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python calculate_coverage_from_blast.py -b blast_output.txt -f subjects.fasta -o coverage.txt
  python calculate_coverage_from_blast.py -b hits.blast -f database.fa -o pos_coverage.tsv
        """
    )

    parser.add_argument('-b', '--blast', required=True,
                       help='BLAST tabular output file (outfmt 6)')
    parser.add_argument('-f', '--fasta', required=True,
                       help='Subject FASTA file')
    parser.add_argument('-o', '--output', required=True,
                       help='Output coverage file (tab-delimited)')
    parser.add_argument('--min-identity', type=float, default=0.0,
                       help='Minimum percent identity to include hit (default: 0.0)')
    parser.add_argument('--max-evalue', type=float, default=float('inf'),
                       help='Maximum e-value to include hit (default: no limit)')

    args = parser.parse_args()

    # Check input files
    if not Path(args.blast).exists():
        print(f"Error: BLAST file {args.blast} not found")
        sys.exit(1)

    if not Path(args.fasta).exists():
        print(f"Error: FASTA file {args.fasta} not found")
        sys.exit(1)

    print(f"Parsing FASTA file: {args.fasta}")
    seq_lengths = parse_fasta_lengths(args.fasta)

    if not seq_lengths:
        print("Error: No sequences found in FASTA file")
        sys.exit(1)

    print(f"Found {len(seq_lengths)} sequences in FASTA file")

    print(f"Parsing BLAST output: {args.blast}")
    hits = parse_blast_output(args.blast)

    if not hits:
        print("Error: No BLAST hits found")
        sys.exit(1)

    print(f"Found {len(hits)} BLAST hits")

    # Apply filters
    filtered_hits = []
    for hit in hits:
        if (hit['pident'] >= args.min_identity and
            hit['evalue'] <= args.max_evalue):
            filtered_hits.append(hit)

    print(f"After filtering: {len(filtered_hits)} hits")

    if not filtered_hits:
        print("Warning: No hits passed filtering criteria")
        # Still create output file with zero coverage
        hits = []
    else:
        hits = filtered_hits

    print("Calculating position-wise coverage...")
    coverage = calculate_position_coverage(hits, seq_lengths)

    print(f"Writing coverage to: {args.output}")
    write_coverage_output(coverage, seq_lengths, args.output)

    # Summary statistics
    total_positions = sum(seq_lengths.values())
    covered_positions = sum(len([pos for pos, count in seq_cov.items() if count > 0])
                           for seq_cov in coverage.values())

    print(f"Summary:")
    print(f"  Total positions: {total_positions}")
    print(f"  Covered positions: {covered_positions}")
    print(f"  Coverage percentage: {covered_positions/total_positions*100:.2f}%")
    print("Coverage calculation completed!")


if __name__ == "__main__":
    main()
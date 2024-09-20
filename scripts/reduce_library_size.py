#!/usr/bin/env python

"""
This script takes library of repeats formated for repeatmasker and reduce the size of
library using cap3 assembly program, requires cap3 to be installed, resulting contigs
are osed only if all input sequences have same classification. If not, contig is split
to individual sequences.

input: repeat library in fasta format
output: reduced repeat library in fasta format
"""

import argparse
import subprocess


def parse_ace(file_path):
    contig_reads = {}
    current_contig = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()

            # Look for the start of a new contig (CO section)
            if line.startswith('CO '):
                parts = line.split()
                current_contig = parts[1]
                contig_reads[current_contig] = []

            # Look for AF (read alignment) lines after a contig has been identified
            elif line.startswith('AF ') and current_contig:
                parts = line.split()
                read_name = parts[1]
                contig_reads[current_contig].append(read_name)
            elif line.startswith('BS ') and current_contig:
                parts = line.split()
                read_name = parts[3]
                contig_reads[current_contig].append(read_name)

    return contig_reads

def evaluate_classification(contigs):
    # classification if the string after `#` check is
    # all reads in contig have the same classification
    # if not contig is split to individual reads
    contigs_to_keep = {}
    reads_to_keep = []
    for contig, reads in contigs.items():
        # reads is list of reads in contig
        # make set of classifications
        classifications = set()
        for read in reads:
            classification = read.split("#")[-1]
            classifications.add(classification)
        if len(classifications) == 1:
            cls = classifications.pop()
            contigs_to_keep[contig] = cls
        else:
            reads_to_keep.extend(reads)

    return contigs_to_keep, reads_to_keep


def filter_fasta_contigs(contig_file, contigs_to_keep, output_file, append=False):
    contig_ok = False
    with open(contig_file, 'r') as file, open(output_file, 'a' if append else 'w') as out:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                contig_name = line[1:]
                if contig_name in contigs_to_keep:
                    # print this contig with classification and subsequence sequences
                    out.write(f">{contig_name}#{contigs_to_keep[contig_name]}\n")
                    contig_ok = True
                else:
                    contig_ok = False
            else:
                if contig_ok:
                    out.write(line + "\n")

def filte_fasta_reads(read_file, reads_to_keep, output_file, append=False):
    read_ok = False
    with open(read_file, 'r') as file, open(output_file, 'a' if append else 'w') as out:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                read_name = line[1:]
                if read_name in reads_to_keep:
                    out.write(line + "\n")
                    read_ok = True
                else:
                    read_ok = False
            else:
                if read_ok:
                    out.write(line + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", help="input fasta file")
    parser.add_argument("--output", help="output fasta file")
    parser.add_argument("--cap3", default="cap3", help="path to cap3 executable")
    args = parser.parse_args()

    ace_output = args.input + ".cap.ace"
    aln_output = args.input + ".cap.aln"
    contigs_output = args.input + ".cap.contigs"
    singlets_output = args.input + ".cap.singlets"

    # run cap3
    cmd = f"{args.cap3} {args.input} -p 80 -o 50"
    subprocess.check_call(cmd, shell=True, stdout=open(aln_output, 'w'))

    # parse ace file
    contigs = parse_ace(ace_output)
    contigs_to_keep, reads_to_keep = evaluate_classification(contigs)

    # write output
    filter_fasta_contigs(contigs_output, contigs_to_keep, args.output)
    filte_fasta_reads(args.input, reads_to_keep, args.output, append=True)
    # append singlets without filtering
    with open(singlets_output, 'r') as file, open(args.output, 'a') as out:
        for line in file:
            out.write(line)



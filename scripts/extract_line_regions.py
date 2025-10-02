#!/usr/bin/env python3
"""
Sophisticated script for extracting DNA sequences based on DANTE annotation patterns.
Extracts regions containing specific domain patterns (ENDO-RT or ENDO-RT-RH) from LINE elements
with strand-aware ordering and distance constraints.

Features:
- Pattern matching: ENDO-RT (2 features) or ENDO-RT-RH (3 features)
- Strand-aware ordering: + strand (ENDO->RT->RH), - strand (RH->RT->ENDO)
- Distance constraints between features (default 2000bp)
- Flanking region extraction with overlap detection
- Reverse complement for minus strand features
- Grouped output with matching sequence and GFF3 features
"""

import argparse
import sys
import subprocess
import tempfile
import os
import re
from pathlib import Path
from collections import defaultdict, namedtuple
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Set


@dataclass
class GFF3Feature:
    seqname: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: str

    def get_attribute(self, key: str) -> Optional[str]:
        """Extract specific attribute value from attributes string."""
        pattern = f'{key}=([^;]+)'
        match = re.search(pattern, self.attributes)
        return match.group(1) if match else None

    def get_name(self) -> Optional[str]:
        """Get the Name attribute."""
        return self.get_attribute('Name')

    def has_line_classification(self) -> bool:
        """Check if feature has Final_Classification=Class_I|LINE."""
        return 'Final_Classification=Class_I|LINE' in self.attributes


@dataclass
class FeatureGroup:
    group_id: str
    seqname: str
    strand: str
    features: List[GFF3Feature]
    pattern_type: str  # "ENDO-RT" or "ENDO-RT-RH"

    def get_region_bounds(self) -> Tuple[int, int]:
        """Get start and end coordinates of the entire group."""
        starts = [f.start for f in self.features]
        ends = [f.end for f in self.features]
        return min(starts), max(ends)


def parse_gff3_features(gff3_file: str) -> List[GFF3Feature]:
    """Parse GFF3 file and return sorted features."""
    features = []

    with open(gff3_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if line.startswith('#') or not line:
                continue

            fields = line.split('\t')
            if len(fields) != 9:
                print(f"Warning: Invalid GFF3 line {line_num}: {line}", file=sys.stderr)
                continue

            try:
                feature = GFF3Feature(
                    seqname=fields[0],
                    source=fields[1],
                    feature=fields[2],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    score=fields[5],
                    strand=fields[6],
                    phase=fields[7],
                    attributes=fields[8]
                )
                features.append(feature)
            except ValueError as e:
                print(f"Warning: Error parsing line {line_num}: {e}", file=sys.stderr)
                continue

    # Sort by seqname, then by start position
    features.sort(key=lambda x: (x.seqname, x.start))
    return features


def get_line_features(features: List[GFF3Feature]) -> List[GFF3Feature]:
    """Filter features that have Final_Classification=Class_I|LINE."""
    return [f for f in features if f.has_line_classification()]


def group_features_by_sequence_and_strand(features: List[GFF3Feature]) -> Dict[Tuple[str, str], List[GFF3Feature]]:
    """Group features by sequence name and strand."""
    groups = defaultdict(list)
    for feature in features:
        key = (feature.seqname, feature.strand)
        groups[key].append(feature)
    return dict(groups)


def find_domain_patterns(features: List[GFF3Feature], max_distance: int) -> List[FeatureGroup]:
    """Find ENDO-RT and ENDO-RT-RH patterns within distance constraints."""
    patterns = []
    group_counter = 0

    # Group by sequence and strand
    seq_strand_groups = group_features_by_sequence_and_strand(features)

    for (seqname, strand), seq_features in seq_strand_groups.items():
        # Sort by position
        seq_features.sort(key=lambda x: x.start)

        # Get features with Name attributes
        named_features = [(i, f) for i, f in enumerate(seq_features) if f.get_name()]

        # Find patterns
        used_indices = set()

        for i in range(len(named_features)):
            if named_features[i][0] in used_indices:
                continue

            # Try to find 3-feature pattern first (ENDO-RT-RH)
            pattern_3 = find_three_feature_pattern(named_features, i, strand, max_distance)
            if pattern_3:
                group_counter += 1
                group_id = f"LINE_group_{group_counter:04d}"
                features_list = [named_features[j][1] for j in pattern_3]
                patterns.append(FeatureGroup(
                    group_id=group_id,
                    seqname=seqname,
                    strand=strand,
                    features=features_list,
                    pattern_type="ENDO-RT-RH"
                ))
                used_indices.update(pattern_3)
                continue

            # Try to find 2-feature pattern (ENDO-RT)
            pattern_2 = find_two_feature_pattern(named_features, i, strand, max_distance, used_indices)
            if pattern_2:
                group_counter += 1
                group_id = f"LINE_group_{group_counter:04d}"
                features_list = [named_features[j][1] for j in pattern_2]
                patterns.append(FeatureGroup(
                    group_id=group_id,
                    seqname=seqname,
                    strand=strand,
                    features=features_list,
                    pattern_type="ENDO-RT"
                ))
                used_indices.update(pattern_2)

    return patterns


def find_three_feature_pattern(named_features: List[Tuple[int, GFF3Feature]],
                              start_idx: int, strand: str, max_distance: int) -> Optional[List[int]]:
    """Find ENDO-RT-RH pattern (3 features)."""
    if start_idx + 2 >= len(named_features):
        return None

    # Get the three consecutive features
    f1_idx, f1 = named_features[start_idx]
    f2_idx, f2 = named_features[start_idx + 1]
    f3_idx, f3 = named_features[start_idx + 2]

    # Check if they are consecutive in the original list
    if f2_idx != f1_idx + 1 or f3_idx != f2_idx + 1:
        return None

    # Get names
    name1, name2, name3 = f1.get_name(), f2.get_name(), f3.get_name()

    # Check pattern based on strand
    if strand == '+':
        expected_pattern = ['ENDO', 'RT', 'RH']
    else:
        expected_pattern = ['RH', 'RT', 'ENDO']

    if [name1, name2, name3] != expected_pattern:
        return None

    # Check distances
    dist1 = f2.start - f1.end
    dist2 = f3.start - f2.end

    if dist1 > max_distance or dist2 > max_distance:
        return None

    return [f1_idx, f2_idx, f3_idx]


def find_two_feature_pattern(named_features: List[Tuple[int, GFF3Feature]],
                            start_idx: int, strand: str, max_distance: int,
                            used_indices: Set[int]) -> Optional[List[int]]:
    """Find ENDO-RT pattern (2 features)."""
    if start_idx + 1 >= len(named_features):
        return None

    # Get the two consecutive features
    f1_idx, f1 = named_features[start_idx]
    f2_idx, f2 = named_features[start_idx + 1]

    # Check if already used
    if f1_idx in used_indices or f2_idx in used_indices:
        return None

    # Check if they are consecutive in the original list
    if f2_idx != f1_idx + 1:
        return None

    # Get names
    name1, name2 = f1.get_name(), f2.get_name()

    # Check pattern based on strand
    if strand == '+':
        expected_pattern = ['ENDO', 'RT']
    else:
        expected_pattern = ['RT', 'ENDO']

    if [name1, name2] != expected_pattern:
        return None

    # Check distance
    distance = f2.start - f1.end
    if distance > max_distance:
        return None

    return [f1_idx, f2_idx]


def get_flanking_regions(pattern: FeatureGroup, all_features: List[GFF3Feature],
                        flank_size: int) -> Tuple[int, int]:
    """Calculate flanking regions with overlap detection."""
    region_start, region_end = pattern.get_region_bounds()

    # Initial flanking coordinates
    flank_start = region_start - flank_size
    flank_end = region_end + flank_size

    # Find overlapping features on the same sequence
    seq_features = [f for f in all_features if f.seqname == pattern.seqname]

    # Adjust left flank
    for feature in seq_features:
        if (feature.end < region_start and feature.end > flank_start and
            feature not in pattern.features):
            flank_start = feature.end + 1

    # Adjust right flank
    for feature in seq_features:
        if (feature.start > region_end and feature.start < flank_end and
            feature not in pattern.features):
            flank_end = feature.start - 1

    # Ensure we don't go below 1
    flank_start = max(1, flank_start)

    return flank_start, flank_end


def create_extraction_bed(patterns: List[FeatureGroup], all_features: List[GFF3Feature],
                         flank_size: int, bed_file: str, include_flanking: bool = False) -> None:
    """Create BED file for sequence extraction.

    Args:
        patterns: List of feature groups to extract
        all_features: All GFF3 features for overlap detection
        flank_size: Size of flanking regions (only used if include_flanking=True)
        bed_file: Output BED file path
        include_flanking: If True, include flanking regions; if False, extract only core region
    """
    with open(bed_file, 'w') as f:
        for pattern in patterns:
            if include_flanking:
                start, end = get_flanking_regions(pattern, all_features, flank_size)
            else:
                # Extract only the core region defined by the features
                start, end = pattern.get_region_bounds()

            # BED format: chr, start, end, name
            # seqkit will use the name field as the sequence ID
            f.write(f"{pattern.seqname}\t{start}\t{end}\t{pattern.group_id}\n")


def extract_all_sequences(genome_fasta: str, bed_file: str, patterns: List[FeatureGroup],
                         all_features: List[GFF3Feature], output_files: Dict, flank_size: int) -> bool:
    """Extract all required sequences: full regions, grouped by pattern, and 5'/3' prime sequences."""
    try:
        # Extract main sequences with flanking regions
        temp_fasta = tempfile.mktemp(suffix='.fasta')
        cmd = ['seqkit', 'subseq', '--bed', bed_file, genome_fasta]
        with open(temp_fasta, 'w') as outf:
            result = subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE,
                                  text=True, check=True)

        # Process main sequences for strand orientation
        process_strand_orientation_with_seqkit(temp_fasta, str(output_files['main_fasta']), patterns)

        # Separate patterns by type
        endo_rt_patterns = [p for p in patterns if p.pattern_type == "ENDO-RT"]
        endo_rt_rh_patterns = [p for p in patterns if p.pattern_type == "ENDO-RT-RH"]

        # Extract grouped sequences from the processed main FASTA
        if endo_rt_patterns:
            extract_grouped_sequences(str(output_files['main_fasta']), str(output_files['endo_rt_fasta']),
                                    endo_rt_patterns, patterns)

        if endo_rt_rh_patterns:
            extract_grouped_sequences(str(output_files['main_fasta']), str(output_files['endo_rt_rh_fasta']),
                                    endo_rt_rh_patterns, patterns)

        # Extract 5' and 3' prime sequences
        if endo_rt_patterns:
            extract_prime_sequences(genome_fasta, endo_rt_patterns, all_features, flank_size,
                                   str(output_files['endo_rt_5prime']), str(output_files['endo_rt_3prime']))

        if endo_rt_rh_patterns:
            extract_prime_sequences(genome_fasta, endo_rt_rh_patterns, all_features, flank_size,
                                   str(output_files['endo_rt_rh_5prime']), str(output_files['endo_rt_rh_3prime']))

        # Clean up
        os.remove(temp_fasta)
        return True

    except subprocess.CalledProcessError as e:
        print(f"Error running seqkit: {e}")
        print(f"stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print("Error: seqkit not found. Please ensure seqkit is installed and in PATH.")
        return False


def extract_sequences_with_seqkit(genome_fasta: str, bed_file: str, output_fasta: str,
                                 patterns: List[FeatureGroup]) -> bool:
    """Extract sequences and handle reverse complement for minus strand using seqkit."""
    try:
        # First extract all sequences
        temp_fasta = tempfile.mktemp(suffix='.fasta')

        cmd = ['seqkit', 'subseq', '--bed', bed_file, genome_fasta]
        with open(temp_fasta, 'w') as outf:
            result = subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE,
                                  text=True, check=True)

        # Process sequences for strand orientation using seqkit
        process_strand_orientation_with_seqkit(temp_fasta, output_fasta, patterns)

        # Clean up
        os.remove(temp_fasta)
        return True

    except subprocess.CalledProcessError as e:
        print(f"Error running seqkit: {e}")
        print(f"stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print("Error: seqkit not found. Please ensure seqkit is installed and in PATH.")
        return False


def extract_group_id_from_header(header: str) -> str:
    """Extract group_id from seqkit header format.

    seqkit adds coordinates like: >chr1:100-200 LINE_group_0001
    We need to extract just LINE_group_0001
    """
    # Remove leading '>'
    header = header.lstrip('>')
    # Split by whitespace and get the last part which should be our group_id
    parts = header.split()
    if parts:
        # The group_id might have suffixes like _5prime, _3prime
        # Extract the base group_id or the full ID with suffix
        return parts[-1]
    return header


def process_strand_orientation_with_seqkit(input_fasta: str, output_fasta: str,
                                          patterns: List[FeatureGroup]) -> None:
    """Process sequences to ensure forward orientation using seqkit for reverse complement."""
    # Create lookup for strand info using group_id only
    strand_lookup = {p.group_id: p.strand for p in patterns}

    # Separate plus and minus strand sequences
    plus_sequences = []
    minus_sequences = []

    with open(input_fasta, 'r') as f:
        current_header = None
        current_seq = []

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Process previous sequence if any
                if current_header and current_seq:
                    seq_id = extract_group_id_from_header(current_header)
                    # Remove suffix for strand lookup
                    base_seq_id = seq_id.replace('_5prime', '').replace('_3prime', '')

                    # Use base_seq_id (without suffix) for the header
                    sequence_data = (f">{base_seq_id}", ''.join(current_seq))

                    if base_seq_id in strand_lookup and strand_lookup[base_seq_id] == '-':
                        minus_sequences.append(sequence_data)
                    else:
                        plus_sequences.append(sequence_data)

                current_header = line
                current_seq = []
            else:
                current_seq.append(line)

        # Process last sequence
        if current_header and current_seq:
            seq_id = extract_group_id_from_header(current_header)
            # Remove suffix for strand lookup
            base_seq_id = seq_id.replace('_5prime', '').replace('_3prime', '')

            # Use base_seq_id (without suffix) for the header
            sequence_data = (f">{base_seq_id}", ''.join(current_seq))

            if base_seq_id in strand_lookup and strand_lookup[base_seq_id] == '-':
                minus_sequences.append(sequence_data)
            else:
                plus_sequences.append(sequence_data)

    # Write output, processing minus strand sequences with seqkit
    with open(output_fasta, 'w') as outf:
        # Write plus strand sequences directly
        for header, sequence in plus_sequences:
            outf.write(f"{header}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(sequence), 80):
                outf.write(f"{sequence[i:i+80]}\n")

        # Process minus strand sequences with seqkit revcom
        if minus_sequences:
            process_minus_strand_sequences(minus_sequences, outf)


def process_minus_strand_sequences(minus_sequences: List[Tuple[str, str]], output_file) -> None:
    """Process minus strand sequences using seqkit revcom."""
    if not minus_sequences:
        return

    # Create temporary file for minus strand sequences
    temp_minus = tempfile.mktemp(suffix='.fasta')
    temp_revcomp = tempfile.mktemp(suffix='.fasta')

    try:
        # Write minus strand sequences to temporary file
        with open(temp_minus, 'w') as f:
            for header, sequence in minus_sequences:
                f.write(f"{header}\n")
                # Write sequence in 80-character lines
                for i in range(0, len(sequence), 80):
                    f.write(f"{sequence[i:i+80]}\n")

        # Use seqkit to reverse complement
        cmd = ['seqkit', 'seq', '--reverse', '--complement', temp_minus]
        with open(temp_revcomp, 'w') as outf:
            result = subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE,
                                  text=True, check=True)

        # Read reverse complemented sequences and add _revcomp suffix
        with open(temp_revcomp, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Add revcomp suffix to header
                    output_file.write(f"{line}_revcomp\n")
                else:
                    output_file.write(f"{line}\n")

    except subprocess.CalledProcessError as e:
        print(f"Error running seqkit revcom: {e}")
        print(f"stderr: {e.stderr}")
        raise

    finally:
        # Clean up temporary files
        for temp_file in [temp_minus, temp_revcomp]:
            if os.path.exists(temp_file):
                os.remove(temp_file)


def extract_grouped_sequences(main_fasta: str, output_fasta: str,
                            target_patterns: List[FeatureGroup], all_patterns: List[FeatureGroup]) -> None:
    """Extract sequences for specific pattern groups from main FASTA."""
    # Use only group_id for matching
    target_ids = {p.group_id for p in target_patterns}

    with open(main_fasta, 'r') as inf, open(output_fasta, 'w') as outf:
        current_header = None
        current_seq = []
        write_sequence = False

        for line in inf:
            line = line.strip()
            if line.startswith('>'):
                # Process previous sequence if needed
                if write_sequence and current_header and current_seq:
                    outf.write(f"{current_header}\n")
                    for seq_line in current_seq:
                        outf.write(f"{seq_line}\n")

                # Check if this sequence should be included
                seq_id = line[1:].split()[0]
                # Handle both original IDs and those with _revcomp suffix
                base_seq_id = seq_id.replace('_revcomp', '')
                write_sequence = base_seq_id in target_ids
                current_header = line
                current_seq = []
            else:
                if write_sequence:
                    current_seq.append(line)

        # Process last sequence
        if write_sequence and current_header and current_seq:
            outf.write(f"{current_header}\n")
            for seq_line in current_seq:
                outf.write(f"{seq_line}\n")


def extract_prime_sequences(genome_fasta: str, patterns: List[FeatureGroup],
                          all_features: List[GFF3Feature], flank_size: int,
                          output_5prime: str, output_3prime: str) -> None:
    """Extract 5' and 3' prime sequences around ENDO and RT/RH domains."""
    # Create BED files for 5' and 3' sequences
    bed_5prime = tempfile.mktemp(suffix='_5prime.bed')
    bed_3prime = tempfile.mktemp(suffix='_3prime.bed')

    try:
        create_prime_bed_files(patterns, flank_size, bed_5prime, bed_3prime)

        # Extract 5' prime sequences
        if os.path.getsize(bed_5prime) > 0:
            temp_5prime = tempfile.mktemp(suffix='_5prime.fasta')
            cmd = ['seqkit', 'subseq', '--bed', bed_5prime, genome_fasta]
            with open(temp_5prime, 'w') as outf:
                subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, text=True, check=True)

            process_strand_orientation_with_seqkit(temp_5prime, output_5prime, patterns)
            os.remove(temp_5prime)

        # Extract 3' prime sequences
        if os.path.getsize(bed_3prime) > 0:
            temp_3prime = tempfile.mktemp(suffix='_3prime.fasta')
            cmd = ['seqkit', 'subseq', '--bed', bed_3prime, genome_fasta]
            with open(temp_3prime, 'w') as outf:
                subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, text=True, check=True)

            process_strand_orientation_with_seqkit(temp_3prime, output_3prime, patterns)
            os.remove(temp_3prime)

    finally:
        # Clean up BED files
        for bed_file in [bed_5prime, bed_3prime]:
            if os.path.exists(bed_file):
                os.remove(bed_file)


def create_prime_bed_files(patterns: List[FeatureGroup], flank_size: int,
                          bed_5prime: str, bed_3prime: str) -> None:
    """Create BED files for 5' and 3' prime sequence extraction.

    For plus strand: 5' is upstream, 3' is downstream (relative to feature direction)
    For minus strand: After reverse complement, what was downstream becomes 5', upstream becomes 3'
    So we need to swap the labels for minus strand to maintain biological orientation.
    """
    with open(bed_5prime, 'w') as f5, open(bed_3prime, 'w') as f3:
        for pattern in patterns:
            # Find ENDO domain and last domain (RT or RH)
            if pattern.strand == '+':
                # Plus strand: ENDO is first, RT/RH is last
                endo_feature = pattern.features[0]  # Should be ENDO
                last_feature = pattern.features[-1]  # Should be RT or RH

                # 5' upstream of ENDO (will remain 5' after extraction)
                prime5_start = max(1, endo_feature.start - flank_size)
                prime5_end = endo_feature.start - 1

                # 3' downstream of RT/RH (will remain 3' after extraction)
                prime3_start = last_feature.end + 1
                prime3_end = last_feature.end + flank_size

            else:
                # Minus strand: RH/RT is first, ENDO is last
                # After reverse complement: what was downstream (after ENDO) becomes 5'
                #                          what was upstream (before RH/RT) becomes 3'
                first_feature = pattern.features[0]  # Should be RH or RT
                endo_feature = pattern.features[-1]  # Should be ENDO

                # Downstream of ENDO → will become 5' after revcomp
                prime5_start = endo_feature.end + 1
                prime5_end = endo_feature.end + flank_size

                # Upstream of first feature → will become 3' after revcomp
                prime3_start = max(1, first_feature.start - flank_size)
                prime3_end = first_feature.start - 1

            # Write BED entries if regions are valid
            if prime5_end >= prime5_start:
                f5.write(f"{pattern.seqname}\t{prime5_start}\t{prime5_end}\t{pattern.group_id}_5prime\n")

            if prime3_end >= prime3_start:
                f3.write(f"{pattern.seqname}\t{prime3_start}\t{prime3_end}\t{pattern.group_id}_3prime\n")



def write_grouped_gff(patterns: List[FeatureGroup], output_gff: str) -> None:
    """Write GFF3 file with grouped features."""
    with open(output_gff, 'w') as f:
        f.write("##gff-version 3\n")

        for pattern in patterns:
            # Write all features in the group
            for feature in pattern.features:
                # Add group ID to attributes
                new_attributes = f"{feature.attributes};Group_ID={pattern.group_id};Pattern_Type={pattern.pattern_type}"

                f.write(f"{feature.seqname}\t{feature.source}\t{feature.feature}\t"
                       f"{feature.start}\t{feature.end}\t{feature.score}\t"
                       f"{feature.strand}\t{feature.phase}\t{new_attributes}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Sophisticated extraction of LINE regions with domain patterns",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Pattern Requirements:
  - Two-feature pattern: ENDO-RT (strand +) or RT-ENDO (strand -)
  - Three-feature pattern: ENDO-RT-RH (strand +) or RH-RT-ENDO (strand -)
  - Three-feature patterns take precedence over two-feature patterns
  - Features must be consecutive without interruption
  - Distance between features must be <= max_distance

Examples:
  python extract_line_regions.py -g genome.fa -a dante.gff3 -o sequences.fa --gff-out features.gff3
  python extract_line_regions.py -g genome.fa -a dante.gff3 -o seq.fa -d 1500 -f 5000
        """
    )

    parser.add_argument('-g', '--genome', required=True,
                       help='Input genome FASTA file')
    parser.add_argument('-a', '--annotations', required=True,
                       help='Input GFF3 file with DANTE annotations')
    parser.add_argument('-o', '--output-dir', required=True,
                       help='Output directory for all result files')
    parser.add_argument('-d', '--max-distance', type=int, default=2000,
                       help='Maximum distance between features in pattern (default: 2000)')
    parser.add_argument('-f', '--flank', type=int, default=10000,
                       help='Flanking region size in bp (default: 10000)')
    parser.add_argument('-b', '--bed-output',
                       help='Optional: save BED file with coordinates')
    parser.add_argument('--keep-bed', action='store_true',
                       help='Keep the intermediate BED file')

    args = parser.parse_args()

    # Validate input files
    for file_path, name in [(args.genome, "Genome"), (args.annotations, "Annotations")]:
        if not Path(file_path).exists():
            print(f"Error: {name} file {file_path} not found")
            sys.exit(1)

    print(f"Parsing GFF3 file: {args.annotations}")
    all_features = parse_gff3_features(args.annotations)
    print(f"Found {len(all_features)} total features")

    print("Filtering LINE features...")
    line_features = get_line_features(all_features)
    print(f"Found {len(line_features)} LINE features")

    if not line_features:
        print("No LINE features found with Final_Classification=Class_I|LINE")
        sys.exit(1)

    print(f"Finding domain patterns (max distance: {args.max_distance}bp)...")
    patterns = find_domain_patterns(line_features, args.max_distance)

    if not patterns:
        print("No valid domain patterns found")
        sys.exit(1)

    print(f"Found {len(patterns)} valid patterns:")
    pattern_counts = defaultdict(int)
    for p in patterns:
        pattern_counts[p.pattern_type] += 1

    for pattern_type, count in pattern_counts.items():
        print(f"  {pattern_type}: {count}")

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Define output files
    output_files = {
        'main_fasta': output_dir / 'LINE_regions.fasta',
        'endo_rt_fasta': output_dir / 'ENDO_RT_regions.fasta',
        'endo_rt_rh_fasta': output_dir / 'ENDO_RT_RH_regions.fasta',
        'endo_rt_5prime': output_dir / 'ENDO_RT_5prime.fasta',
        'endo_rt_3prime': output_dir / 'ENDO_RT_3prime.fasta',
        'endo_rt_rh_5prime': output_dir / 'ENDO_RT_RH_5prime.fasta',
        'endo_rt_rh_3prime': output_dir / 'ENDO_RT_RH_3prime.fasta',
        'gff_out': output_dir / 'grouped_features.gff3'
    }

    # Set up BED file
    if args.bed_output:
        bed_file = args.bed_output
        keep_bed = True
    else:
        bed_file = tempfile.mktemp(suffix='.bed')
        keep_bed = args.keep_bed

    print(f"Creating BED file for core regions (no flanking)...")
    create_extraction_bed(patterns, all_features, args.flank, bed_file, include_flanking=False)

    print(f"Extracting sequences from {args.genome}...")
    success = extract_all_sequences(args.genome, bed_file, patterns, all_features,
                                   output_files, args.flank)

    print(f"Writing grouped GFF3 to {output_files['gff_out']}...")
    write_grouped_gff(patterns, str(output_files['gff_out']))

    # Clean up BED file
    if not keep_bed and os.path.exists(bed_file):
        os.remove(bed_file)
        print("Temporary BED file removed")
    elif keep_bed:
        print(f"BED file saved: {bed_file}")

    if success:
        print("Extraction completed successfully!")
        print(f"Output files in {output_dir}:")
        for name, path in output_files.items():
            if path.exists():
                print(f"  {name}: {path.name}")
        sys.exit(0)
    else:
        print("Extraction failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
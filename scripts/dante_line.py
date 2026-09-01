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
- All-vs-all alignment analysis of flanking sequences
"""

import argparse
import sys
import subprocess
import tempfile
import os
import re
import heapq
import math
from pathlib import Path
import bisect
from collections import defaultdict, namedtuple
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Set

sys.path.insert(0, str(Path(__file__).resolve().parent))
import classification  # noqa: E402

# Import alignment function from global_local_aln module
try:
    from global_local_aln import run_all_vs_all_alignment
except ImportError:
    # Try importing from same directory
    import importlib.util
    spec = importlib.util.spec_from_file_location("global_local_aln",
                                                   Path(__file__).parent / "global_local_aln.py")
    global_local_aln = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(global_local_aln)
    run_all_vs_all_alignment = global_local_aln.run_all_vs_all_alignment

# Import the canonical FASTA sorter (deterministic clustering input ordering)
try:
    from canonical_sort_fasta import sort_fasta_by_sequence
except ImportError:
    import importlib.util
    _spec = importlib.util.spec_from_file_location(
        "canonical_sort_fasta", Path(__file__).parent / "canonical_sort_fasta.py")
    canonical_sort_fasta = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(canonical_sort_fasta)
    sort_fasta_by_sequence = canonical_sort_fasta.sort_fasta_by_sequence


def load_sequence_lengths(genome_fasta: str) -> Dict[str, int]:
    """Load sequence lengths from FASTA index file (.fai).

    If .fai file doesn't exist, it will be created using seqkit faidx.

    Parameters
    ----------
    genome_fasta : str
        Path to genome FASTA file

    Returns
    -------
    Dict[str, int]
        Dictionary mapping sequence names to their lengths
    """
    fai_file = f"{genome_fasta}.fai"

    # Check if .fai file exists, if not create it
    if not Path(fai_file).exists():
        print(f"FASTA index not found, creating {fai_file}...")
        try:
            subprocess.run(['seqkit', 'faidx', genome_fasta],
                         check=True, capture_output=True, text=True)
            print(f"  → Created {Path(fai_file).name}")
        except subprocess.CalledProcessError as e:
            print(f"Error creating FASTA index: {e}", file=sys.stderr)
            print(f"stderr: {e.stderr}", file=sys.stderr)
            raise
        except FileNotFoundError:
            print("Error: seqkit not found. Please ensure seqkit is installed and in PATH.", file=sys.stderr)
            raise

    # Read .fai file
    seq_lengths = {}
    with open(fai_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                seq_name = fields[0]
                seq_length = int(fields[1])
                seq_lengths[seq_name] = seq_length

    return seq_lengths


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


def parse_gff3_intervals(gff3_file: str) -> List[Tuple[str, int, int]]:
    """Parse a GFF3 into bare ``(seqname, start, end)`` triples.

    Used for ``--mask-gff3`` inputs, which are only ever asked where the nearest
    feature boundary is. Keeping the attribute column would dominate memory: a
    DANTE_LTR or DANTE GFF3 stores the protein sequence of every domain in
    column 9.
    """
    intervals = []
    with open(gff3_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                if line.strip():
                    print(f"Warning: Invalid GFF3 line {line_num} in {gff3_file}",
                          file=sys.stderr)
                continue
            try:
                intervals.append((fields[0], int(fields[3]), int(fields[4])))
            except ValueError as e:
                print(f"Warning: Error parsing line {line_num} in {gff3_file}: {e}",
                      file=sys.stderr)
    return intervals


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


class FeatureIndex:
    """Per-seqname index for O(log n) nearest-boundary lookups.

    Flank clipping used to rescan the whole feature list per pattern
    (``[f for f in all_features if f.seqname == ...]`` inside the per-pattern
    loop), i.e. O(patterns x features). On a large genome the mask set (raw
    TideHunter arrays plus the DANTE_LTR / DANTE_TIR elements) holds millions of
    records, so that quadratic scan dominates wall-time. This indexes features once and answers
    the two flank-clipping queries by binary search.

    Both queries reproduce the exact result of the original linear loops:
      * ``nearest_end_below(seqname, ref, floor)`` -> the largest feature end
        that is ``< ref`` (and ``> floor``), plus 1; otherwise ``floor``.
        (Original: ``for f: if f.end < ref and f.end > floor: floor = f.end+1``.)
      * ``nearest_start_above(seqname, ref, ceil)`` -> the smallest feature start
        that is ``> ref`` (and ``< ceil``), minus 1; otherwise ``ceil``.
        (Original: ``for f: if f.start > ref and f.start < ceil: ceil = f.start-1``.)
    """

    def __init__(self, features: List[GFF3Feature]):
        self._build((f.seqname, f.start, f.end) for f in features)

    @classmethod
    def from_intervals(cls, intervals) -> "FeatureIndex":
        """Build from bare ``(seqname, start, end)`` triples.

        Mask files are only ever queried for boundaries, so they are parsed with
        ``parse_gff3_intervals`` rather than into ``GFF3Feature`` objects -- a
        DANTE_LTR GFF3 carries protein sequences in column 9 and is multiple GB
        on a large genome, none of which this index needs.
        """
        index = cls.__new__(cls)
        index._build(intervals)
        return index

    def _build(self, intervals) -> None:
        ends = defaultdict(list)
        starts = defaultdict(list)
        for seqname, start, end in intervals:
            ends[seqname].append(end)
            starts[seqname].append(start)
        self._ends = {sn: sorted(v) for sn, v in ends.items()}
        self._starts = {sn: sorted(v) for sn, v in starts.items()}

    def nearest_end_below(self, seqname: str, ref: int, floor: int) -> int:
        ends = self._ends.get(seqname)
        if not ends:
            return floor
        i = bisect.bisect_left(ends, ref)   # first end >= ref
        if i == 0:
            return floor                     # no end < ref
        best = ends[i - 1]                   # largest end < ref
        return best + 1 if best > floor else floor

    def nearest_start_above(self, seqname: str, ref: int, ceil: int) -> int:
        starts = self._starts.get(seqname)
        if not starts:
            return ceil
        i = bisect.bisect_right(starts, ref)  # first start > ref
        if i == len(starts):
            return ceil                        # no start > ref
        best = starts[i]                       # smallest start > ref
        return best - 1 if best < ceil else ceil


def get_flanking_regions(pattern: FeatureGroup, all_index: "FeatureIndex",
                        flank_size: int, mask_index: "FeatureIndex" = None) -> Tuple[int, int]:
    """Calculate flanking regions with overlap detection.

    Parameters
    ----------
    pattern : FeatureGroup
        The feature group to calculate flanking regions for
    all_index : FeatureIndex
        Index over all DANTE annotation features
    flank_size : int
        Maximum size of flanking regions
    mask_index : FeatureIndex, optional
        Index over additional masking features that limit flanking regions

    Returns
    -------
    Tuple[int, int]
        Start and end coordinates of flanking region
    """
    region_start, region_end = pattern.get_region_bounds()

    # Initial flanking coordinates
    flank_start = region_start - flank_size
    flank_end = region_end + flank_size

    # Clip against neighbouring DANTE features, then mask features (mask uses the
    # DANTE-updated flanks as its floor/ceiling, matching the original order).
    # The pattern's own features never satisfy end<region_start / start>region_end
    # (they lie within [region_start, region_end]), so the old `not in
    # pattern.features` guard was a no-op and is dropped.
    sn = pattern.seqname
    flank_start = all_index.nearest_end_below(sn, region_start, flank_start)
    flank_end = all_index.nearest_start_above(sn, region_end, flank_end)
    if mask_index is not None:
        flank_start = mask_index.nearest_end_below(sn, region_start, flank_start)
        flank_end = mask_index.nearest_start_above(sn, region_end, flank_end)

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
    # Build the feature index once (not per pattern) when flanking is needed.
    all_index = FeatureIndex(all_features) if include_flanking else None
    with open(bed_file, 'w') as f:
        for pattern in patterns:
            if include_flanking:
                start, end = get_flanking_regions(pattern, all_index, flank_size)
            else:
                # Extract only the core region defined by the features
                start, end = pattern.get_region_bounds()

            # BED format: chr, start, end, name
            # seqkit will use the name field as the sequence ID
            f.write(f"{pattern.seqname}\t{start}\t{end}\t{pattern.group_id}\n")


def extract_all_sequences(genome_fasta: str, bed_file: str, patterns: List[FeatureGroup],
                         all_features: List[GFF3Feature], output_files: Dict, flank_size: int,
                         mask_index: "FeatureIndex" = None, seq_lengths: Dict[str, int] = None) -> bool:
    """Extract all required sequences: full regions, grouped by pattern, and 5'/3' prime sequences."""
    # Index the DANTE features once: create_prime_bed_files is called for both
    # pattern types, and on a large genome rebuilding it per call is expensive.
    all_index = FeatureIndex(all_features)
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

        # Extract 5' and 3' prime sequences with mask support and sequence boundaries
        if endo_rt_patterns:
            extract_prime_sequences(genome_fasta, endo_rt_patterns, all_index, flank_size,
                                   str(output_files['endo_rt_5prime']), str(output_files['endo_rt_3prime']),
                                   mask_index, seq_lengths)

        if endo_rt_rh_patterns:
            extract_prime_sequences(genome_fasta, endo_rt_rh_patterns, all_index, flank_size,
                                   str(output_files['endo_rt_rh_5prime']), str(output_files['endo_rt_rh_3prime']),
                                   mask_index, seq_lengths)

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
                          all_index: "FeatureIndex", flank_size: int,
                          output_5prime: str, output_3prime: str,
                          mask_index: "FeatureIndex" = None,
                          seq_lengths: Dict[str, int] = None) -> None:
    """Extract 5' and 3' prime sequences around ENDO and RT/RH domains."""
    # Create BED files for 5' and 3' sequences
    bed_5prime = tempfile.mktemp(suffix='_5prime.bed')
    bed_3prime = tempfile.mktemp(suffix='_3prime.bed')

    try:
        create_prime_bed_files(patterns, all_index, flank_size, bed_5prime, bed_3prime, mask_index, seq_lengths)

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


def create_prime_bed_files(patterns: List[FeatureGroup], all_index: "FeatureIndex",
                          flank_size: int, bed_5prime: str, bed_3prime: str,
                          mask_index: "FeatureIndex" = None,
                          seq_lengths: Dict[str, int] = None) -> None:
    """Create BED files for 5' and 3' prime sequence extraction with mask support.

    For plus strand: 5' is upstream, 3' is downstream (relative to feature direction)
    For minus strand: After reverse complement, what was downstream becomes 5', upstream becomes 3'
    So we need to swap the labels for minus strand to maintain biological orientation.

    Flanking regions are limited by neighboring DANTE features, optional mask features, and sequence boundaries.
    """
    with open(bed_5prime, 'w') as f5, open(bed_3prime, 'w') as f3:
        for pattern in patterns:
            # Get sequence length if available
            seq_length = seq_lengths.get(pattern.seqname) if seq_lengths else None
            sn = pattern.seqname
            # Find ENDO domain and last domain (RT or RH)
            if pattern.strand == '+':
                # Plus strand: ENDO is first, RT/RH is last
                endo_feature = pattern.features[0]  # Should be ENDO
                last_feature = pattern.features[-1]  # Should be RT or RH

                # Initial 5' upstream of ENDO
                prime5_start = max(1, endo_feature.start - flank_size)
                prime5_end = endo_feature.start - 1

                # Initial 3' downstream of RT/RH
                prime3_start = last_feature.end + 1
                prime3_end = last_feature.end + flank_size

                # Clip 5' (left of ENDO) and 3' (right of last) by DANTE, then mask.
                prime5_start = all_index.nearest_end_below(sn, endo_feature.start, prime5_start)
                prime3_end = all_index.nearest_start_above(sn, last_feature.end, prime3_end)
                if mask_index is not None:
                    prime5_start = mask_index.nearest_end_below(sn, endo_feature.start, prime5_start)
                    prime3_end = mask_index.nearest_start_above(sn, last_feature.end, prime3_end)


            else:
                # Minus strand: RH/RT is first, ENDO is last
                first_feature = pattern.features[0]  # Should be RH or RT
                endo_feature = pattern.features[-1]  # Should be ENDO

                # Initial downstream of ENDO → will become 5' after revcomp
                prime5_start = endo_feature.end + 1
                prime5_end = endo_feature.end + flank_size

                # Initial upstream of first feature → will become 3' after revcomp
                prime3_start = max(1, first_feature.start - flank_size)
                prime3_end = first_feature.start - 1

                # Clip 5' (right of ENDO) and 3' (left of first) by DANTE, then mask.
                prime5_end = all_index.nearest_start_above(sn, endo_feature.end, prime5_end)
                prime3_start = all_index.nearest_end_below(sn, first_feature.start, prime3_start)
                if mask_index is not None:
                    prime5_end = mask_index.nearest_start_above(sn, endo_feature.end, prime5_end)
                    prime3_start = mask_index.nearest_end_below(sn, first_feature.start, prime3_start)

            # Constrain to sequence boundaries if seq_length is available
            if seq_length:
                prime5_end = min(prime5_end, seq_length)
                prime3_end = min(prime3_end, seq_length)

            # Write BED entries if regions are valid (start < end)
            if prime5_end > prime5_start:
                f5.write(f"{pattern.seqname}\t{prime5_start}\t{prime5_end}\t{pattern.group_id}_5prime\n")

            if prime3_end > prime3_start:
                f3.write(f"{pattern.seqname}\t{prime3_start}\t{prime3_end}\t{pattern.group_id}_3prime\n")



@dataclass
class LineElement:
    """Represents a complete LINE element with boundaries."""
    group_id: str
    seqname: str
    strand: str
    start: int
    end: int
    pattern_type: str
    extension_5prime: int = 0
    extension_3prime: int = 0
    # What the flank alignment proposed before the length caps were applied.
    # Equal to the applied values when nothing was capped.
    extension_5prime_inferred: int = 0
    extension_3prime_inferred: int = 0


def load_alignment_lengths(output_dir: Path) -> Dict[str, Dict[str, int]]:
    """Load alignment length data from TSV files.

    Returns dict: {group_id: {'5prime': length, '3prime': length}}
    """
    alignment_lengths = defaultdict(dict)

    length_files = [
        'ENDO_RT_5prime_aln_length.tsv',
        'ENDO_RT_3prime_aln_length.tsv',
        'ENDO_RT_RH_5prime_aln_length.tsv',
        'ENDO_RT_RH_3prime_aln_length.tsv',
    ]

    for filename in length_files:
        filepath = output_dir / filename
        if not filepath.exists():
            continue

        # Determine if this is 5prime or 3prime
        prime_type = '5prime' if '5prime' in filename else '3prime'

        try:
            with open(filepath, 'r') as f:
                # Read header
                header = f.readline().strip().split('\t')
                group_id_idx = header.index('Group_ID')
                selected_length_idx = header.index('Selected_Length')

                # Read data rows
                for line in f:
                    fields = line.strip().split('\t')
                    group_id = fields[group_id_idx]
                    selected_length = int(fields[selected_length_idx])
                    alignment_lengths[group_id][prime_type] = selected_length
        except Exception as e:
            print(f"Warning: Could not load {filename}: {e}", file=sys.stderr)

    return dict(alignment_lengths)


def line_max_extension(side: str, explicit: int = 0) -> int:
    """Per-side flank bound in bp for Class_I/LINE; 0 means unbounded.

    ``explicit`` (the CLI --max-extension) wins when set, so the old symmetric
    behaviour stays reachable. Otherwise the bound comes from
    ``max_extension_per_side`` in classification_vocabulary.yaml, which is
    ASYMMETRIC for LINE: the domain core is ORF2 (ENDO+RT), so the 5' side still
    has inter-ORF + ORF1 + 5'UTR to cover (~2 kb) while the 3' side has only the
    ORF2 C-terminus + 3'UTR + polyA (~0.8 kb). The previous shared 1500 was
    therefore too tight 5' and too loose 3' at the same time.
    """
    if explicit:
        return explicit
    try:
        return classification.max_extension_for_class("Class_I/LINE", side)
    except Exception:
        return 0


def line_core_to_3prime_end(explicit: int = 0) -> int:
    """Bound on (core length + 3' extension) for Class_I/LINE; 0 = unbounded.

    A non-zero ``explicit`` (the CLI --max-extension) disables this bound as
    well as the per-side ones, so --max-extension is a COMPLETE manual override:
    the number given is the cap, with no vocabulary bound quietly tightening it.
    """
    if explicit:
        return 0
    try:
        return classification.max_core_to_3prime_end_for_class("Class_I/LINE")
    except Exception:
        return 0


def cap_extensions(core_len: int, ext_5prime: int, ext_3prime: int,
                   max_extension: int = 0, max_element_length: int = 0,
                   max_ext_5prime: int = 0, max_ext_3prime: int = 0,
                   max_core_to_3prime_end: int = 0) -> Tuple[int, int]:
    """Clamp a pair of inferred extensions to the biological length bounds.

    The flank alignment can propose an extension all the way to ``--flank``
    (10 kb by default) on each side, which produces "LINE" consensi of 16-22 kb.
    No plant LINE is that long -- full-length elements run 4-7 kb around an
    ENDO+RT core of roughly 2.1 kb -- so the extensions are capped rather than
    the elements discarded: the domain core is always kept, only the appended
    flank is trimmed.

    ``max_ext_5prime`` / ``max_ext_3prime`` bound the two sides SEPARATELY and
    take precedence over ``max_extension`` on the side they are set for; they
    come from ``max_extension_per_side`` in classification_vocabulary.yaml, which
    is per superfamily for TIR (medians span 6x: PIF_Harbinger 951 bp vs
    EnSpm_CACTA 3544 bp) and asymmetric for LINE (the core is ORF2, so the 5'
    side has ORF1 + 5'UTR to cover but the 3' side only 3'UTR + polyA). A single
    shared cap cut into the MEDIAN real element for four of five TIR
    superfamilies.

    ``max_core_to_3prime_end`` bounds ``core_len + ext_3prime`` -- the distance
    from the START of the domain core to the element's 3' end. For LINEs this is
    the bound that works: measured over 148 confirmed loci, the 3' extension
    varies 13.8x (p90/p10) while this span varies 1.16x, because the extension is
    only the remainder after the domain annotation stops (r = -0.936 between core
    length and 3' extension). It is applied after the per-side caps.

    ``max_extension`` bounds each side independently. ``max_element_length``
    then bounds the whole element; when the two sides together exceed the
    remaining budget it is split evenly, except that a side asking for less than
    its half keeps what it asked for and yields the rest to the other side. An
    odd base pair goes to the 3' side. A core already at or over
    ``max_element_length`` yields no extension at all -- the core itself is never
    trimmed.

    Either bound is disabled by passing 0.
    """
    cap5 = max_ext_5prime if max_ext_5prime > 0 else max_extension
    cap3 = max_ext_3prime if max_ext_3prime > 0 else max_extension
    if cap5 > 0:
        ext_5prime = min(ext_5prime, cap5)
    if cap3 > 0:
        ext_3prime = min(ext_3prime, cap3)

    if max_core_to_3prime_end > 0:
        ext_3prime = min(ext_3prime, max(0, max_core_to_3prime_end - core_len))

    if max_element_length > 0:
        budget = max(0, max_element_length - core_len)
        if ext_5prime + ext_3prime > budget:
            half = budget // 2
            if ext_5prime <= half:
                ext_3prime = budget - ext_5prime
            elif ext_3prime <= half:
                ext_5prime = budget - ext_3prime
            else:
                ext_5prime = half
                ext_3prime = budget - half

    return ext_5prime, ext_3prime


def create_line_elements(patterns: List[FeatureGroup], alignment_lengths: Dict[str, Dict[str, int]],
                         max_extension: int = 0, max_element_length: int = 0) -> List[LineElement]:
    """Create LINE_element features with boundaries based on alignment lengths.

    Parameters
    ----------
    patterns : List[FeatureGroup]
        List of feature groups (ENDO-RT or ENDO-RT-RH patterns)
    alignment_lengths : Dict[str, Dict[str, int]]
        Dictionary mapping group_id to 5prime/3prime alignment lengths
    max_extension : int, optional
        Cap on each side's extension in bp; 0 disables
    max_element_length : int, optional
        Cap on the whole element in bp; 0 disables

    Returns
    -------
    List[LineElement]
        List of LINE elements with calculated boundaries
    """
    line_elements = []

    for pattern in patterns:
        # Get base boundaries from features (ENDO to RT/RH)
        base_start, base_end = pattern.get_region_bounds()

        # Get alignment extensions if available
        # For minus strand, the group_id in TSV has _revcomp suffix
        lookup_id = f"{pattern.group_id}_revcomp" if pattern.strand == '-' else pattern.group_id
        extensions = alignment_lengths.get(lookup_id, {})
        raw_5prime = extensions.get('5prime', 0)
        raw_3prime = extensions.get('3prime', 0)

        ext_5prime, ext_3prime = cap_extensions(
            base_end - base_start + 1, raw_5prime, raw_3prime,
            max_extension=max_extension, max_element_length=max_element_length,
            max_ext_5prime=line_max_extension("5prime", max_extension),
            max_ext_3prime=line_max_extension("3prime", max_extension),
            max_core_to_3prime_end=line_core_to_3prime_end(max_extension))

        # Calculate extended boundaries based on strand
        if pattern.strand == '+':
            # Plus strand: 5' is upstream (subtract from start), 3' is downstream (add to end)
            element_start = base_start - ext_5prime if ext_5prime > 0 else base_start
            element_end = base_end + ext_3prime if ext_3prime > 0 else base_end
        else:
            # Minus strand: 5' is downstream (add to end), 3' is upstream (subtract from start)
            # After reverse complement: 5' extension goes to what was downstream (higher coords)
            #                          3' extension goes to what was upstream (lower coords)
            element_start = base_start - ext_3prime if ext_3prime > 0 else base_start
            element_end = base_end + ext_5prime if ext_5prime > 0 else base_end

        # Ensure start <= end and >= 1
        element_start = max(1, element_start)
        element_end = max(element_start, element_end)

        line_elements.append(LineElement(
            group_id=pattern.group_id,
            seqname=pattern.seqname,
            strand=pattern.strand,
            start=element_start,
            end=element_end,
            pattern_type=pattern.pattern_type,
            extension_5prime=ext_5prime,
            extension_3prime=ext_3prime,
            extension_5prime_inferred=raw_5prime,
            extension_3prime_inferred=raw_3prime
        ))

    return line_elements


def extract_extended_line_regions(genome_fasta: str, line_elements: List[LineElement],
                                 output_fasta: str) -> bool:
    """Extract extended LINE regions based on LINE_element boundaries.

    Parameters
    ----------
    genome_fasta : str
        Path to genome FASTA file
    line_elements : List[LineElement]
        List of LINE elements with extended boundaries
    output_fasta : str
        Output FASTA file for extended sequences

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    if not line_elements:
        print("No LINE elements to extract")
        return False

    try:
        # Create BED file for extended regions
        bed_file = tempfile.mktemp(suffix='_extended.bed')

        with open(bed_file, 'w') as f:
            for line_element in line_elements:
                f.write(f"{line_element.seqname}\t{line_element.start}\t{line_element.end}\t{line_element.group_id}\n")

        # Extract sequences using seqkit
        temp_fasta = tempfile.mktemp(suffix='_extended.fasta')
        cmd = ['seqkit', 'subseq', '--bed', bed_file, genome_fasta]
        with open(temp_fasta, 'w') as outf:
            subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, text=True, check=True)

        # Create a list of FeatureGroups from LineElements for strand processing
        pseudo_patterns = [
            FeatureGroup(
                group_id=le.group_id,
                seqname=le.seqname,
                strand=le.strand,
                features=[],  # Not needed for this operation
                pattern_type=le.pattern_type
            )
            for le in line_elements
        ]

        # Process strand orientation
        process_strand_orientation_with_seqkit(temp_fasta, output_fasta, pseudo_patterns)

        # Clean up
        os.remove(bed_file)
        os.remove(temp_fasta)
        return True

    except subprocess.CalledProcessError as e:
        print(f"Error extracting extended regions: {e}")
        print(f"stderr: {e.stderr}")
        return False
    except Exception as e:
        print(f"Error extracting extended regions: {e}")
        return False


def write_grouped_gff(patterns: List[FeatureGroup], output_gff: str, line_elements: List[LineElement] = None) -> None:
    """Write GFF3 file with grouped features and LINE_element features.

    Parameters
    ----------
    patterns : List[FeatureGroup]
        List of feature groups
    output_gff : str
        Output GFF3 file path
    line_elements : List[LineElement], optional
        List of LINE elements with extended boundaries
    """
    with open(output_gff, 'w') as f:
        f.write("##gff-version 3\n")

        for pattern in patterns:
            # Find corresponding LINE_element if available
            line_element = None
            if line_elements:
                line_element = next((le for le in line_elements if le.group_id == pattern.group_id), None)

            # Write LINE_element feature first (parent)
            if line_element:
                attributes = f"ID={pattern.group_id};Pattern_Type={pattern.pattern_type}"
                if line_element.extension_5prime > 0 or line_element.extension_3prime > 0:
                    attributes += f";Extension_5prime={line_element.extension_5prime};Extension_3prime={line_element.extension_3prime}"
                # Only emitted when a cap actually bound, so an uncapped element's
                # attributes are unchanged from previous releases.
                capped = [
                    side for side, applied, inferred in (
                        ("5prime", line_element.extension_5prime, line_element.extension_5prime_inferred),
                        ("3prime", line_element.extension_3prime, line_element.extension_3prime_inferred),
                    ) if applied != inferred
                ]
                if capped:
                    attributes += (f";Extension_capped={'+'.join(capped)}"
                                   f";Extension_5prime_inferred={line_element.extension_5prime_inferred}"
                                   f";Extension_3prime_inferred={line_element.extension_3prime_inferred}")

                f.write(f"{line_element.seqname}\tDANTE\tLINE_element\t"
                       f"{line_element.start}\t{line_element.end}\t.\t"
                       f"{line_element.strand}\t.\t{attributes}\n")

            # Write protein domain features with Parent attribute
            for feature in pattern.features:
                # Add Parent and Group_ID to attributes
                new_attributes = f"{feature.attributes};Parent={pattern.group_id};Group_ID={pattern.group_id};Pattern_Type={pattern.pattern_type}"

                f.write(f"{feature.seqname}\t{feature.source}\t{feature.feature}\t"
                       f"{feature.start}\t{feature.end}\t{feature.score}\t"
                       f"{feature.strand}\t{feature.phase}\t{new_attributes}\n")


def analyze_alignment_lengths(alignment_tsv: Path, output_tsv: Path,
                              min_num_alignments: int,
                              support_fraction: float = 0.0,
                              min_group_alignments: int = 0) -> None:
    """Choose, per group, how far its flank is supported by the other copies.

    For each ``Group_ID`` the alignment table holds one length per pairwise flank
    alignment the group takes part in. The extension kept for the group is the
    ``k``-th largest of those lengths, i.e. the distance out to which at least
    ``k`` partners still align.

    Why ``k`` is not a constant
    ---------------------------
    ``k`` used to be ``min_num_alignments`` (default 3) regardless of how many
    partners a group had, which makes the threshold satisfiable by construction:
    a group with exactly three alignments has its third-largest *equal to its
    minimum*, so nothing is filtered at all. Measured on a wheat run, most
    runaway extensions came from groups of 3-4, and the ones that did not came
    from a handful of outlier pairs in a large group -- ``LINE_group_2147`` took
    a 5,968 bp extension from 3 of its 51 partners while the 25th percentile of
    that same group was 751 bp.

    So ``k`` scales with the evidence available::

        k = max(min_num_alignments, ceil(support_fraction * n_alignments))

    and a group with fewer than ``min_group_alignments`` alignments is not
    reported at all -- too few partners to say anything, so the element keeps its
    domain core and no extension.

    ``support_fraction=0.0`` reproduces the old fixed-``k`` behaviour exactly, and
    takes a **single-pass** code path: ``k`` is then known up front, so the counting
    pass is unnecessary. That matters because ``dante_tir_fallback`` shares this
    function and its alignment tables are far larger than DANTE_LINE's -- on a
    14.5 Gb wheat run the EnSpm_CACTA 3' table alone is **43 GB** (every row carries
    its aligned sequences), so a gratuitous second read costs tens of GB of I/O.

    Parameters
    ----------
    alignment_tsv : Path
        Input alignment TSV file from all-vs-all alignment
    output_tsv : Path
        Output TSV file with length thresholds per group
    min_num_alignments : int
        Floor on ``k``; also the value of ``k`` when ``support_fraction`` is 0
    support_fraction : float, optional
        Fraction of a group's alignments that must reach the kept length
    min_group_alignments : int, optional
        Groups with fewer alignments than this get no extension
    """
    if not 0.0 <= support_fraction <= 1.0:
        raise ValueError("support_fraction must be between 0 and 1")

    def read_lengths():
        """Yield (group_id, length) for both partners of every alignment row."""
        with open(alignment_tsv, 'r') as f:
            header = f.readline().strip().split('\t')
            query_id_idx = header.index('query_id')
            ref_id_idx = header.index('ref_id')
            degapped_query_len_idx = header.index('degapped_query_len')
            degapped_ref_len_idx = header.index('degapped_ref_len')
            for line in f:
                fields = line.strip().split('\t')
                yield fields[query_id_idx], int(fields[degapped_query_len_idx])
                yield fields[ref_id_idx], int(fields[degapped_ref_len_idx])

    group_count = defaultdict(int)
    group_heap = defaultdict(list)

    if support_fraction > 0:
        # k scales with the group, so it cannot be known while streaming the
        # lengths: count first, then keep the k largest per group.
        for group_id, _ in read_lengths():
            group_count[group_id] += 1

        group_k = {}
        for group_id, total in group_count.items():
            k = max(min_num_alignments, math.ceil(total * support_fraction))
            if total < k or total < min_group_alignments:
                continue
            group_k[group_id] = k

        for group_id, length in read_lengths():
            k = group_k.get(group_id)
            if k is None:
                continue
            heap = group_heap[group_id]
            if len(heap) < k:
                heapq.heappush(heap, length)
            else:
                heapq.heappushpop(heap, length)  # keep the k largest
    else:
        # k is a constant, so one pass suffices: a size-k min-heap per group
        # holds the k largest lengths and the count comes along for free.
        k = min_num_alignments
        for group_id, length in read_lengths():
            group_count[group_id] += 1
            heap = group_heap[group_id]
            if len(heap) < k:
                heapq.heappush(heap, length)
            else:
                heapq.heappushpop(heap, length)
        group_k = {g: k for g, total in group_count.items()
                   if total >= k and total >= min_group_alignments}

    # The selected length is the k-th largest (= the heap's min). Num_Shorter is
    # (# lengths <= selected) - 1; every length strictly greater than the k-th
    # largest is itself among the k largest (so it is in the heap), hence
    # (# > selected) = k - (# heap entries == selected).
    results = []
    for group_id in sorted(group_k.keys()):
        total = group_count[group_id]
        k = group_k[group_id]
        heap = group_heap[group_id]
        selected_length = heap[0]
        num_greater = k - sum(1 for x in heap if x == selected_length)
        num_shorter_or_equal = total - num_greater - 1

        results.append({
            'Group_ID': group_id,
            'Selected_Length': selected_length,
            'Num_Shorter': num_shorter_or_equal
        })

    # Write results to TSV atomically (.tmp + rename) so a run killed mid-write
    # cannot leave a truncated file that the caller's `.exists()` checkpoint would
    # accept as a completed length analysis.
    tmp_out = output_tsv.with_name(output_tsv.name + '.tmp')
    with open(tmp_out, 'w') as f:
        # Write header
        f.write('Group_ID\tSelected_Length\tNum_Shorter\n')

        # Write data rows
        for result in results:
            f.write(f"{result['Group_ID']}\t{result['Selected_Length']}\t{result['Num_Shorter']}\n")
    os.replace(tmp_out, output_tsv)


def run_mmseqs_clustering(input_fasta: str, output_dir: Path, threads: int = 1) -> bool:
    """Run MMseqs2 easy-cluster on extended LINE sequences.

    Parameters
    ----------
    input_fasta : str
        Path to input FASTA file (LINE_regions_extended.fasta)
    output_dir : Path
        Output directory for MMseqs2 results
    threads : int, optional
        Number of threads for clustering, default 1

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    if not Path(input_fasta).exists():
        print(f"Error: Input FASTA file {input_fasta} not found", file=sys.stderr)
        return False

    # Create mmseqs subdirectory
    mmseqs_dir = output_dir / 'mmseqs'
    mmseqs_dir.mkdir(parents=True, exist_ok=True)

    # Define output files
    output_prefix = mmseqs_dir / 'cluster'
    rep_seq_fasta = mmseqs_dir / 'cluster_rep_seq.fasta'
    cluster_tsv = mmseqs_dir / 'cluster.tsv'
    tmp_dir = mmseqs_dir / 'tmp'
    done_marker = mmseqs_dir / 'clustering.done'

    # Checkpoint: only trust a run explicitly marked complete. A killed mmseqs
    # can leave a partial cluster_rep_seq.fasta / cluster.tsv, which must NOT be
    # reused as a finished clustering.
    if done_marker.exists() and rep_seq_fasta.exists() and cluster_tsv.exists():
        print(f"  MMseqs2 clustering already exists in {mmseqs_dir}")
        return True

    print(f"\nRunning MMseqs2 clustering on extended LINE sequences...")
    print(f"  Output directory: {mmseqs_dir}")

    try:
        # Deterministic clustering: mmseqs easy-cluster is order-sensitive
        # (same sequences in a different order -> different representatives and
        # a different cluster count). Canonically sort the input by sequence so
        # the LINE library is a reproducible function of the input set, not its
        # (run-varying) order. See scripts/canonical_sort_fasta.py.
        sorted_input = mmseqs_dir / 'input_sorted.fasta'
        sort_fasta_by_sequence(input_fasta, str(sorted_input),
                               threads=threads, tmpdir=str(mmseqs_dir))

        # Run MMseqs2 easy-cluster
        cmd = [
            'mmseqs', 'easy-cluster',
            str(sorted_input),
            str(output_prefix),
            str(tmp_dir),
            '--threads', str(threads)
        ]

        print(f"  Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # Check if output files were created
        if rep_seq_fasta.exists():
            print(f"  → cluster_rep_seq.fasta")
        if cluster_tsv.exists():
            print(f"  → cluster.tsv")

        # Mark complete only once both expected outputs exist, so the checkpoint
        # above distinguishes a finished run from a killed/partial one.
        if rep_seq_fasta.exists() and cluster_tsv.exists():
            done_marker.write_text("ok\n")
        return True

    except subprocess.CalledProcessError as e:
        print(f"Error running MMseqs2: {e}", file=sys.stderr)
        print(f"stdout: {e.stdout}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        return False
    except FileNotFoundError:
        print("Error: mmseqs not found. Please ensure MMseqs2 is installed and in PATH.", file=sys.stderr)
        return False


def run_prime_alignments(output_dir: Path, threads: int = 1, min_num_alignments: int = 3,
                         verbose: bool = False, max_group_size: Optional[int] = None,
                         support_fraction: float = 0.0, min_group_alignments: int = 0) -> None:
    """Run all-vs-all alignment analysis on prime sequences.

    For 5' sequences: use --end 3 (3' end fixed, analyze 5' similarity)
    For 3' sequences: use --end 5 (5' end fixed, analyze 3' similarity)

    Parameters
    ----------
    output_dir : Path
        Output directory containing prime sequence FASTA files
    threads : int, optional
        Number of threads for parallel processing, default 1
    min_num_alignments : int, optional
        Minimum number of alignments for length threshold calculation, default 3
    verbose : bool, optional
        Print verbose progress messages for alignment analysis, default False
    support_fraction : float, optional
        Fraction of a group's alignments that must reach the kept extension
    min_group_alignments : int, optional
        Groups with fewer alignments than this get no extension
    """
    alignment_files = [
        ('ENDO_RT_5prime.fasta', 'ENDO_RT_5prime_alignment.tsv', '3'),  # 5' seqs -> fix 3' end
        ('ENDO_RT_3prime.fasta', 'ENDO_RT_3prime_alignment.tsv', '5'),  # 3' seqs -> fix 5' end
        ('ENDO_RT_RH_5prime.fasta', 'ENDO_RT_RH_5prime_alignment.tsv', '3'),
        ('ENDO_RT_RH_3prime.fasta', 'ENDO_RT_RH_3prime_alignment.tsv', '5'),
    ]

    print("\nRunning all-vs-all alignment analysis on prime sequences...")

    for fasta_name, output_name, end_param in alignment_files:
        fasta_path = output_dir / fasta_name
        output_path = output_dir / output_name

        # Skip if output TSV already exists (checkpoint)
        if output_path.exists():
            print(f"  Skipping {fasta_name} (alignment file {output_name} already exists)")
            continue

        # Skip if FASTA file doesn't exist or is empty
        if not fasta_path.exists():
            print(f"  Skipping {fasta_name} (file not found)")
            continue

        if fasta_path.stat().st_size == 0:
            print(f"  Skipping {fasta_name} (empty file)")
            continue

        print(f"  Analyzing {fasta_name} (--end {end_param})...")

        try:
            run_all_vs_all_alignment(
                fasta_file=str(fasta_path),
                output_file=str(output_path),
                end=end_param,
                gap_open=12,
                gap_extend=3,
                match=2,
                mismatch=-2,
                score_threshold=20,
                threads=threads,
                verbose=verbose,
                max_group_size=max_group_size
            )
            print(f"    → {output_name}")
        except Exception as e:
            print(f"    Error analyzing {fasta_name}: {e}", file=sys.stderr)

    # Analyze alignment lengths for each TSV file
    print("\nAnalyzing alignment lengths and calculating thresholds...")
    for fasta_name, alignment_name, end_param in alignment_files:
        alignment_path = output_dir / alignment_name
        length_output = output_dir / alignment_name.replace('_alignment.tsv', '_aln_length.tsv')

        # Skip if alignment file doesn't exist
        if not alignment_path.exists():
            continue

        # Skip if length file already exists (checkpoint)
        if length_output.exists():
            print(f"  Skipping {alignment_name} (length file already exists)")
            continue

        print(f"  Processing {alignment_name}...")
        try:
            analyze_alignment_lengths(alignment_path, length_output, min_num_alignments,
                                      support_fraction=support_fraction,
                                      min_group_alignments=min_group_alignments)
            print(f"    → {length_output.name}")
        except Exception as e:
            print(f"    Error processing {alignment_name}: {e}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Extraction of LINE regions with domain patterns",
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
    parser.add_argument('-t', '--threads', type=int, default=1,
                       help='Number of threads for alignment analysis (default: 1)')
    parser.add_argument('--min-num-alignments', type=int, default=3,
                       help='Floor on the number of partner alignments that must reach the '
                            'kept extension (default: 3)')
    parser.add_argument('--support-fraction', type=float, default=0.5,
                       help='Fraction of a group\'s flank alignments that must reach the kept '
                            'extension. Scales the support requirement with the evidence '
                            'available, so a few near-identical outlier pairs cannot set the '
                            'boundary for the whole group. 0 restores the old fixed '
                            '--min-num-alignments behaviour (default: 0.5)')
    parser.add_argument('--min-group-alignments', type=int, default=5,
                       help='Groups with fewer flank alignments than this get no extension at '
                            'all: too few partners to place a boundary, so the element keeps '
                            'its domain core. 0 disables (default: 5)')
    parser.add_argument('--max-extension', type=int, default=0,
                       help='SYMMETRIC ESCAPE HATCH: a value > 0 caps BOTH sides at that many '
                            'bp, overriding the per-side bounds. At the default 0 the bounds '
                            'come from max_extension_per_side in classification_vocabulary.yaml, '
                            'which is asymmetric for LINE (2000 bp 5\', 800 bp 3\'): the core is '
                            'ORF2, so ORF1 + 5\'UTR lie outside it 5\' but only 3\'UTR + polyA 3\'. '
                            'Deliberately tighter than a full-length LINE needs either way: a '
                            'truncated consensus still masks its family, an over-extended one '
                            'mislabels its neighbour genome-wide (default: 0)')
    parser.add_argument('--library-source', choices=['core', 'element'], default='core',
                       help="Which span of each element seeds the RepeatMasker library. "
                            "'core' (default) uses only the span DANTE actually annotated "
                            "-- the ENDO..RT/RH domain region -- and leaves the inferred "
                            "flanks out of the similarity search entirely, while "
                            "DANTE_LINE.gff3 still reports the full element. 'element' "
                            "restores the previous behaviour of clustering the extended "
                            "regions. The flanks are an inference, and a wrong one is "
                            "amplified by every RepeatMasker search that follows.")
    parser.add_argument('--max-element-length', type=int, default=8000,
                       help='Cap on the total LINE element length in bp; the domain core is '
                            'never trimmed, only the appended flanks. 0 disables (default: 8000)')
    parser.add_argument('--mask-gff3', action='append', metavar='GFF3',
                       help='GFF3 whose features terminate a flank. Repeatable: pass '
                            'the TideHunter tandem arrays plus DANTE_LTR / DANTE_TIR '
                            'elements so a flank cannot run through a neighbouring '
                            'element that another tool has already annotated.')
    parser.add_argument('--max-group-size', type=int, default=None,
                       help='When more than this many LINE patterns are found, split them '
                            'into deterministic clustering-based groups of at most this size '
                            'and run the all-vs-all flank alignment per group (bounds O(N^2) '
                            'memory). Below the threshold the run is byte-identical to no '
                            'grouping. Default: unset.')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Print verbose progress messages for alignment analysis (default: False)')

    args = parser.parse_args()
    if args.max_group_size is not None and args.max_group_size < 2:
        parser.error("--max-group-size must be >= 2")
    if not 0.0 <= args.support_fraction <= 1.0:
        parser.error("--support-fraction must be between 0 and 1")
    for name, value in (("--min-group-alignments", args.min_group_alignments),
                        ("--max-extension", args.max_extension),
                        ("--max-element-length", args.max_element_length)):
        if value < 0:
            parser.error(f"{name} must be >= 0")

    # Validate input files
    for file_path, name in [(args.genome, "Genome"), (args.annotations, "Annotations")]:
        if not Path(file_path).exists():
            print(f"Error: {name} file {file_path} not found")
            sys.exit(1)

    print(f"Parsing GFF3 file: {args.annotations}")
    all_features = parse_gff3_features(args.annotations)
    print(f"Found {len(all_features)} total features")

    # Parse mask GFF3s if provided. Several may be given; they are merged into a
    # single boundary index, so the flank stops at whichever comes first.
    mask_index = None
    if args.mask_gff3:
        mask_intervals = []
        for mask_path in args.mask_gff3:
            if not Path(mask_path).exists():
                print(f"Error: Mask GFF3 file {mask_path} not found")
                sys.exit(1)
            print(f"Parsing mask GFF3 file: {mask_path}")
            found = parse_gff3_intervals(mask_path)
            print(f"  {len(found)} mask features")
            mask_intervals.extend(found)
        print(f"Found {len(mask_intervals)} mask features in total")
        mask_index = FeatureIndex.from_intervals(mask_intervals)
        del mask_intervals

    print("Filtering LINE features...")
    line_features = get_line_features(all_features)
    print(f"Found {len(line_features)} LINE features")

    if not line_features:
        # Exit 3 = "no LINE content" (benign): the genome legitimately has too
        # few LINE features. The dante_line Snakefile rule treats ONLY exit 3 as
        # the empty-output case; every other non-zero exit is a real failure that
        # must stop the pipeline (rather than being silently turned into empty
        # LINE outputs).
        print("No LINE features found with Final_Classification=Class_I|LINE")
        sys.exit(3)

    print(f"Finding domain patterns (max distance: {args.max_distance}bp)...")
    patterns = find_domain_patterns(line_features, args.max_distance)

    if not patterns:
        print("No valid domain patterns found")
        sys.exit(3)  # benign "no LINE content" — see above

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
        'gff_out': output_dir / 'DANTE_LINE.gff3',
        'line_rep_lib': output_dir / 'LINE_rep_lib.fasta'
    }

    # Set up BED file
    if args.bed_output:
        bed_file = args.bed_output
        keep_bed = True
    else:
        bed_file = tempfile.mktemp(suffix='.bed')
        keep_bed = args.keep_bed

    # Load sequence lengths from genome FASTA index
    print(f"Loading sequence lengths from genome index...")
    seq_lengths = load_sequence_lengths(args.genome)
    print(f"  Loaded lengths for {len(seq_lengths)} sequences")

    print(f"Creating BED file for core regions (no flanking)...")
    create_extraction_bed(patterns, all_features, args.flank, bed_file, include_flanking=False)

    print(f"Extracting sequences from {args.genome}...")
    success = extract_all_sequences(args.genome, bed_file, patterns, all_features,
                                   output_files, args.flank, mask_index, seq_lengths)

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

        # Run alignment analysis on prime sequences
        run_prime_alignments(output_dir, threads=args.threads,
                             min_num_alignments=args.min_num_alignments,
                             verbose=args.verbose, max_group_size=args.max_group_size,
                             support_fraction=args.support_fraction,
                             min_group_alignments=args.min_group_alignments)

        # Load alignment lengths and create LINE elements
        print("\nCreating LINE_element features from alignment data...")
        alignment_lengths = load_alignment_lengths(output_dir)
        line_elements = create_line_elements(patterns, alignment_lengths,
                                             max_extension=args.max_extension,
                                             max_element_length=args.max_element_length)
        n_capped = sum(1 for le in line_elements
                       if le.extension_5prime != le.extension_5prime_inferred
                       or le.extension_3prime != le.extension_3prime_inferred)
        if n_capped:
            print(f"  {n_capped} of {len(line_elements)} elements had a boundary capped "
                  f"(--max-extension {args.max_extension}, "
                  f"--max-element-length {args.max_element_length})")

        # Write GFF3 with LINE_element features
        print(f"Writing GFF3 with LINE_element features to {output_files['gff_out']}...")
        write_grouped_gff(patterns, str(output_files['gff_out']), line_elements)

        # Extract extended LINE regions based on LINE_element boundaries. This
        # file is always written -- it is the record of the element extent the
        # boundary search inferred -- but see --library-source for whether it is
        # what seeds RepeatMasker.
        extended_fasta_path = output_dir / 'LINE_regions_extended.fasta'
        print(f"\nExtracting extended LINE regions to {extended_fasta_path}...")
        extended_success = extract_extended_line_regions(args.genome, line_elements, str(extended_fasta_path))

        if extended_success:
            print(f"  → LINE_regions_extended.fasta")

            # The library is built from the DANTE-annotated core by default: the
            # domain span is observed, the flanks are inferred, and a wrong flank
            # in the library is amplified across the whole RepeatMasker search.
            # LINE_regions.fasta is that core span, already strand-normalised.
            if args.library_source == 'core':
                library_source_path = output_files['main_fasta']
                print(f"  Library source: LINE_regions.fasta (DANTE domain core; "
                      f"inferred flanks excluded from the RepeatMasker library)")
            else:
                library_source_path = extended_fasta_path
                print(f"  Library source: LINE_regions_extended.fasta (core + inferred flanks)")

            # Run MMseqs2 clustering on the chosen library source
            clustering_success = run_mmseqs_clustering(str(library_source_path), output_dir, threads=args.threads)

            if clustering_success:
                # Process cluster representative sequences to create repeat library
                print(f"\nCreating LINE repeat library from cluster representatives...")
                mmseqs_rep_fasta = output_dir / 'mmseqs' / 'cluster_rep_seq.fasta'
                line_rep_lib = output_files['line_rep_lib']
                if mmseqs_rep_fasta.exists():
                    try:
                        with open(mmseqs_rep_fasta, 'r') as inf, open(line_rep_lib, 'w') as outf:
                            for line in inf:
                                if line.startswith('>'):
                                    # Add #Class_I/LINE suffix to sequence name
                                    seq_name = line.strip()
                                    outf.write(f"{seq_name}#Class_I/LINE\n")
                                else:
                                    outf.write(line)
                        print(f"  → LINE_rep_lib.fasta ({line_rep_lib})")
                    except Exception as e:
                        print(f"  Warning: Failed to create LINE repeat library: {e}", file=sys.stderr)
                else:
                    print(f"  Warning: Cluster representative file not found: {mmseqs_rep_fasta}", file=sys.stderr)
            else:
                print("  Warning: MMseqs2 clustering failed", file=sys.stderr)
        else:
            print("  Warning: Failed to extract extended LINE regions", file=sys.stderr)

        print("\nAll processing completed successfully!")
        sys.exit(0)
    else:
        print("Extraction failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
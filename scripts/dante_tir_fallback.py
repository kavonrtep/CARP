#!/usr/bin/env python3
"""
Fallback extraction pipeline for DANTE TIR annotations.

This script mirrors the LINE workflow in dante_line.py, but it operates on
TPase-only TIR annotations:

- keep only features with Final_Classification=Class_II|Subclass_1|TIR|*
- process each TIR subtype independently
- use the TPase feature as the anchor
- extract strand-aware 5' and 3' flanks
- run all-vs-all flank alignments to infer extension lengths
- extract extended regions and cluster them with mmseqs easy-cluster

The output is a fallback annotation and library, not a validated TIR caller.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from dante_line import (
        GFF3Feature,
        analyze_alignment_lengths,
        load_sequence_lengths,
        parse_gff3_features,
        process_strand_orientation_with_seqkit,
        run_mmseqs_clustering,
    )
except ImportError:
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "dante_line", Path(__file__).parent / "dante_line.py"
    )
    dante_line = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dante_line)
    GFF3Feature = dante_line.GFF3Feature
    analyze_alignment_lengths = dante_line.analyze_alignment_lengths
    load_sequence_lengths = dante_line.load_sequence_lengths
    parse_gff3_features = dante_line.parse_gff3_features
    process_strand_orientation_with_seqkit = dante_line.process_strand_orientation_with_seqkit
    run_mmseqs_clustering = dante_line.run_mmseqs_clustering

try:
    from global_local_aln import run_all_vs_all_alignment
except ImportError:
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "global_local_aln", Path(__file__).parent / "global_local_aln.py"
    )
    global_local_aln = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(global_local_aln)
    run_all_vs_all_alignment = global_local_aln.run_all_vs_all_alignment


SOURCE_NAME = "DANTE_TIR_FALLBACK"
GFF_PARENT_TYPE = "sequence_feature"


@dataclass
class TIRAnchor:
    group_id: str
    seqname: str
    strand: str
    subtype: str
    feature: GFF3Feature


@dataclass
class TIRFallbackElement:
    group_id: str
    seqname: str
    strand: str
    subtype: str
    start: int
    end: int
    extension_5prime: int = 0
    extension_3prime: int = 0


def sanitize_label(label: str) -> str:
    """Convert a classification label into a filesystem-safe token."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", label).strip("_")


def get_tir_subtype(feature: GFF3Feature) -> Optional[str]:
    """Return the TIR subtype from Final_Classification if present."""
    name = feature.get_name()
    if name != "TPase":
        return None

    classification = feature.get_attribute("Final_Classification")
    if not classification or not classification.startswith("Class_II|Subclass_1|TIR|"):
        return None

    parts = classification.split("|")
    if len(parts) < 4:
        return None
    return parts[-1]


def build_tir_anchors(features: List[GFF3Feature]) -> Dict[str, List[TIRAnchor]]:
    """Group TPase features by TIR subtype and assign stable group IDs."""
    subtype_groups: Dict[str, List[GFF3Feature]] = defaultdict(list)
    for feature in features:
        subtype = get_tir_subtype(feature)
        if subtype:
            subtype_groups[subtype].append(feature)

    anchors: Dict[str, List[TIRAnchor]] = {}
    for subtype, subtype_features in subtype_groups.items():
        subtype_features.sort(key=lambda x: (x.seqname, x.start, x.end, x.strand))
        subtype_slug = sanitize_label(subtype)
        classification_id = f"Class_II_Subclass_1_TIR_{subtype_slug}"
        subtype_anchors: List[TIRAnchor] = []
        for idx, feature in enumerate(subtype_features, 1):
            # _partial suffix is embedded in the ID from the start so it flows
            # through BED → seqkit → FASTA → GFF3 uniformly; downstream code
            # (merge_tir_fallback.py) no longer needs to append it.
            group_id = f"{classification_id}_{idx}_partial"
            subtype_anchors.append(
                TIRAnchor(
                    group_id=group_id,
                    seqname=feature.seqname,
                    strand=feature.strand,
                    subtype=subtype,
                    feature=feature,
                )
            )
        anchors[subtype] = subtype_anchors

    return anchors


def make_core_bed(anchors: List[TIRAnchor], bed_file: str) -> None:
    """Write BED coordinates for TPase anchor regions."""
    with open(bed_file, "w") as handle:
        for anchor in anchors:
            handle.write(
                f"{anchor.seqname}\t{anchor.feature.start}\t{anchor.feature.end}\t{anchor.group_id}\n"
            )


def create_prime_bed_files(
    anchors: List[TIRAnchor],
    all_features: List[GFF3Feature],
    flank_size: int,
    bed_5prime: str,
    bed_3prime: str,
    mask_features: List[GFF3Feature] = None,
    seq_lengths: Dict[str, int] = None,
) -> None:
    """Create strand-aware BED files for 5' and 3' flanking regions."""
    with open(bed_5prime, "w") as f5, open(bed_3prime, "w") as f3:
        for anchor in anchors:
            feature = anchor.feature
            seq_length = seq_lengths.get(anchor.seqname) if seq_lengths else None

            if anchor.strand == "+":
                prime5_start = max(1, feature.start - flank_size)
                prime5_end = feature.start - 1
                prime3_start = feature.end + 1
                prime3_end = feature.end + flank_size

                seq_features = [f for f in all_features if f.seqname == anchor.seqname]
                for other in seq_features:
                    if other is feature:
                        continue
                    if other.end < feature.start and other.end > prime5_start:
                        prime5_start = other.end + 1
                for other in seq_features:
                    if other is feature:
                        continue
                    if other.start > feature.end and other.start < prime3_end:
                        prime3_end = other.start - 1

                if mask_features:
                    mask_seq_features = [f for f in mask_features if f.seqname == anchor.seqname]
                    for other in mask_seq_features:
                        if other.end < feature.start and other.end > prime5_start:
                            prime5_start = other.end + 1
                    for other in mask_seq_features:
                        if other.start > feature.end and other.start < prime3_end:
                            prime3_end = other.start - 1
            else:
                prime5_start = feature.end + 1
                prime5_end = feature.end + flank_size
                prime3_start = max(1, feature.start - flank_size)
                prime3_end = feature.start - 1

                seq_features = [f for f in all_features if f.seqname == anchor.seqname]
                for other in seq_features:
                    if other is feature:
                        continue
                    if other.start > feature.end and other.start < prime5_end:
                        prime5_end = other.start - 1
                for other in seq_features:
                    if other is feature:
                        continue
                    if other.end < feature.start and other.end > prime3_start:
                        prime3_start = other.end + 1

                if mask_features:
                    mask_seq_features = [f for f in mask_features if f.seqname == anchor.seqname]
                    for other in mask_seq_features:
                        if other.start > feature.end and other.start < prime5_end:
                            prime5_end = other.start - 1
                    for other in mask_seq_features:
                        if other.end < feature.start and other.end > prime3_start:
                            prime3_start = other.end + 1

            if seq_length:
                prime5_end = min(prime5_end, seq_length)
                prime3_end = min(prime3_end, seq_length)

            if prime5_end > prime5_start:
                f5.write(
                    f"{anchor.seqname}\t{prime5_start}\t{prime5_end}\t{anchor.group_id}_5prime\n"
                )
            if prime3_end > prime3_start:
                f3.write(
                    f"{anchor.seqname}\t{prime3_start}\t{prime3_end}\t{anchor.group_id}_3prime\n"
                )


def extract_sequences_with_seqkit(
    genome_fasta: str,
    bed_file: str,
    output_fasta: str,
    anchors: List[TIRAnchor],
) -> bool:
    """Extract BED regions and normalize minus-strand sequences."""
    temp_fasta = tempfile.mktemp(suffix=".fasta")
    try:
        cmd = ["seqkit", "subseq", "--bed", bed_file, genome_fasta]
        with open(temp_fasta, "w") as outf:
            subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, text=True, check=True)

        process_strand_orientation_with_seqkit(temp_fasta, output_fasta, anchors)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running seqkit: {e}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        return False
    except FileNotFoundError:
        print("Error: seqkit not found. Please ensure seqkit is installed and in PATH.", file=sys.stderr)
        return False
    finally:
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)


def extract_prime_sequences(
    genome_fasta: str,
    anchors: List[TIRAnchor],
    all_features: List[GFF3Feature],
    flank_size: int,
    output_5prime: str,
    output_3prime: str,
    mask_features: List[GFF3Feature] = None,
    seq_lengths: Dict[str, int] = None,
) -> bool:
    """Extract 5' and 3' flanking sequences around TPase anchors."""
    bed_5prime = tempfile.mktemp(suffix="_5prime.bed")
    bed_3prime = tempfile.mktemp(suffix="_3prime.bed")
    temp_5prime = None
    temp_3prime = None

    try:
        create_prime_bed_files(
            anchors,
            all_features,
            flank_size,
            bed_5prime,
            bed_3prime,
            mask_features,
            seq_lengths,
        )

        if os.path.getsize(bed_5prime) > 0:
            temp_5prime = tempfile.mktemp(suffix="_5prime.fasta")
            cmd = ["seqkit", "subseq", "--bed", bed_5prime, genome_fasta]
            with open(temp_5prime, "w") as outf:
                subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, text=True, check=True)
            process_strand_orientation_with_seqkit(temp_5prime, output_5prime, anchors)

        if os.path.getsize(bed_3prime) > 0:
            temp_3prime = tempfile.mktemp(suffix="_3prime.fasta")
            cmd = ["seqkit", "subseq", "--bed", bed_3prime, genome_fasta]
            with open(temp_3prime, "w") as outf:
                subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, text=True, check=True)
            process_strand_orientation_with_seqkit(temp_3prime, output_3prime, anchors)

        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running seqkit: {e}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        return False
    except FileNotFoundError:
        print("Error: seqkit not found. Please ensure seqkit is installed and in PATH.", file=sys.stderr)
        return False
    finally:
        for bed_file in [bed_5prime, bed_3prime, temp_5prime, temp_3prime]:
            if bed_file and os.path.exists(bed_file):
                os.remove(bed_file)


def run_prime_alignments(
    subtype_dir: Path,
    threads: int = 1,
    min_num_alignments: int = 3,
    verbose: bool = False,
) -> None:
    """Run all-vs-all alignment analysis for TIR flanks."""
    alignment_files = [
        ("TPase_5prime.fasta", "TPase_5prime_alignment.tsv", "3"),
        ("TPase_3prime.fasta", "TPase_3prime_alignment.tsv", "5"),
    ]

    for fasta_name, output_name, end_param in alignment_files:
        fasta_path = subtype_dir / fasta_name
        output_path = subtype_dir / output_name

        if output_path.exists():
            continue
        if not fasta_path.exists() or fasta_path.stat().st_size == 0:
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
            )
            print(f"    -> {output_name}")
        except Exception as e:
            print(f"    Error analyzing {fasta_name}: {e}", file=sys.stderr)

    for alignment_name in ["TPase_5prime_alignment.tsv", "TPase_3prime_alignment.tsv"]:
        alignment_path = subtype_dir / alignment_name
        length_output = subtype_dir / alignment_name.replace("_alignment.tsv", "_aln_length.tsv")
        if not alignment_path.exists() or length_output.exists():
            continue

        try:
            analyze_alignment_lengths(alignment_path, length_output, min_num_alignments)
            print(f"    -> {length_output.name}")
        except Exception as e:
            print(f"    Error processing {alignment_name}: {e}", file=sys.stderr)


def load_prime_alignment_lengths(subtype_dir: Path) -> Dict[str, Dict[str, int]]:
    """Load inferred 5'/3' lengths for TPase anchors."""
    alignment_lengths = defaultdict(dict)

    for prime_type in ["5prime", "3prime"]:
        filepath = subtype_dir / f"TPase_{prime_type}_aln_length.tsv"
        if not filepath.exists():
            continue

        try:
            with open(filepath, "r") as handle:
                header = handle.readline().strip().split("\t")
                group_id_idx = header.index("Group_ID")
                selected_length_idx = header.index("Selected_Length")

                for line in handle:
                    fields = line.strip().split("\t")
                    if len(fields) <= max(group_id_idx, selected_length_idx):
                        continue
                    group_id = fields[group_id_idx]
                    selected_length = int(fields[selected_length_idx])
                    alignment_lengths[group_id][prime_type] = selected_length
        except Exception as e:
            print(f"Warning: Could not load {filepath.name}: {e}", file=sys.stderr)

    return dict(alignment_lengths)


def create_tir_elements(
    anchors: List[TIRAnchor],
    alignment_lengths: Dict[str, Dict[str, int]],
) -> List[TIRFallbackElement]:
    """Create extended fallback elements from inferred flank lengths."""
    tir_elements: List[TIRFallbackElement] = []

    for anchor in anchors:
        base_start = anchor.feature.start
        base_end = anchor.feature.end
        lookup_id = f"{anchor.group_id}_revcomp" if anchor.strand == "-" else anchor.group_id
        extensions = alignment_lengths.get(lookup_id, {})
        ext_5prime = extensions.get("5prime", 0)
        ext_3prime = extensions.get("3prime", 0)

        if anchor.strand == "+":
            element_start = base_start - ext_5prime if ext_5prime > 0 else base_start
            element_end = base_end + ext_3prime if ext_3prime > 0 else base_end
        else:
            element_start = base_start - ext_3prime if ext_3prime > 0 else base_start
            element_end = base_end + ext_5prime if ext_5prime > 0 else base_end

        element_start = max(1, element_start)
        element_end = max(element_start, element_end)

        tir_elements.append(
            TIRFallbackElement(
                group_id=anchor.group_id,
                seqname=anchor.seqname,
                strand=anchor.strand,
                subtype=anchor.subtype,
                start=element_start,
                end=element_end,
                extension_5prime=ext_5prime,
                extension_3prime=ext_3prime,
            )
        )

    return tir_elements


def extract_extended_tir_regions(
    genome_fasta: str,
    tir_elements: List[TIRFallbackElement],
    output_fasta: str,
) -> bool:
    """Extract extended fallback regions and normalize strand orientation."""
    if not tir_elements:
        return False

    bed_file = tempfile.mktemp(suffix="_extended.bed")
    temp_fasta = tempfile.mktemp(suffix="_extended.fasta")

    try:
        with open(bed_file, "w") as handle:
            for element in tir_elements:
                handle.write(
                    f"{element.seqname}\t{element.start}\t{element.end}\t{element.group_id}\n"
                )

        cmd = ["seqkit", "subseq", "--bed", bed_file, genome_fasta]
        with open(temp_fasta, "w") as outf:
            subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, text=True, check=True)

        process_strand_orientation_with_seqkit(temp_fasta, output_fasta, tir_elements)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error extracting extended regions: {e}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error extracting extended regions: {e}", file=sys.stderr)
        return False
    finally:
        for path in [bed_file, temp_fasta]:
            if os.path.exists(path):
                os.remove(path)


def write_tir_gff(
    anchors: List[TIRAnchor],
    tir_elements: List[TIRFallbackElement],
    output_gff: str,
) -> None:
    """Write a fallback GFF3 with inferred parent elements and TPase anchors."""
    element_map = {element.group_id: element for element in tir_elements}

    with open(output_gff, "w") as handle:
        handle.write("##gff-version 3\n")
        handle.write(f"##{SOURCE_NAME.lower()} version 0.1\n")

        for anchor in anchors:
            element = element_map.get(anchor.group_id)
            subtype_slug = sanitize_label(anchor.subtype)
            classification_id = f"Class_II_Subclass_1_TIR_{subtype_slug}"

            if element:
                attributes = (
                    f"ID={anchor.group_id};Classification={classification_id};"
                    f"Subtype={anchor.subtype};Method={SOURCE_NAME};Anchor=TPase"
                )
                if element.extension_5prime > 0 or element.extension_3prime > 0:
                    attributes += (
                        f";Extension_5prime={element.extension_5prime}"
                        f";Extension_3prime={element.extension_3prime}"
                    )
                handle.write(
                    f"{element.seqname}\t{SOURCE_NAME}\t{GFF_PARENT_TYPE}\t"
                    f"{element.start}\t{element.end}\t.\t{element.strand}\t.\t{attributes}\n"
                )

            child_attributes = (
                f"{anchor.feature.attributes};Parent={anchor.group_id};"
                f"Group_ID={anchor.group_id};Subtype={anchor.subtype};Method={SOURCE_NAME}"
            )
            handle.write(
                f"{anchor.feature.seqname}\t{anchor.feature.source}\t{anchor.feature.feature}\t"
                f"{anchor.feature.start}\t{anchor.feature.end}\t{anchor.feature.score}\t"
                f"{anchor.feature.strand}\t{anchor.feature.phase}\t{child_attributes}\n"
            )


def append_gff_body(source_gff: Path, destination_gff_handle) -> None:
    """Append a GFF file body without its header lines."""
    with open(source_gff, "r") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            destination_gff_handle.write(line)


def append_rep_library(rep_fasta: Path, output_handle, classification_suffix: str) -> None:
    """Append cluster representatives to the combined repeat library."""
    with open(rep_fasta, "r") as handle:
        for line in handle:
            if line.startswith(">"):
                output_handle.write(f"{line.strip()}#{classification_suffix}\n")
            else:
                output_handle.write(line)


def annotate_extended_fasta_headers(fasta_path: str, classification_with_slashes: str) -> None:
    """Rewrite headers in a fallback extended FASTA to append ``#classification``.

    The seqkit → strand-orientation pipeline writes headers like
    ``>Class_II_Subclass_1_TIR_EnSpm_CACTA_123_partial`` without any
    classification suffix. reduce_library_size.R relies on the
    ``>unique_id#classification`` naming convention to group sequences
    for per-class reduction, so we append the slash-form classification
    here. Idempotent — strips any existing ``#...`` tail before adding.
    """
    if not os.path.exists(fasta_path) or os.path.getsize(fasta_path) == 0:
        return
    tmp = tempfile.mktemp(suffix=".fasta")
    with open(fasta_path, "r") as fin, open(tmp, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                header = line.rstrip("\n").rstrip("\r")
                base = header.split("#", 1)[0]
                fout.write(f"{base}#{classification_with_slashes}\n")
            else:
                fout.write(line)
    os.replace(tmp, fasta_path)


def process_subtype(
    subtype: str,
    anchors: List[TIRAnchor],
    genome_fasta: str,
    all_features: List[GFF3Feature],
    output_dir: Path,
    flank_size: int,
    threads: int,
    min_num_alignments: int,
    mask_features: List[GFF3Feature] = None,
    seq_lengths: Dict[str, int] = None,
    verbose: bool = False,
) -> Tuple[Optional[Path], Optional[Path], int]:
    """Run the fallback workflow for one TIR subtype."""
    subtype_slug = sanitize_label(subtype)
    subtype_dir = output_dir / subtype_slug
    subtype_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nProcessing subtype {subtype} ({len(anchors)} TPase anchors)")
    print(f"  Output directory: {subtype_dir}")

    output_files = {
        "core_fasta": subtype_dir / "TPase_regions.fasta",
        "5prime_fasta": subtype_dir / "TPase_5prime.fasta",
        "3prime_fasta": subtype_dir / "TPase_3prime.fasta",
        "gff_out": subtype_dir / "DANTE_TIR_FALLBACK.gff3",
        "extended_fasta": subtype_dir / "TPase_regions_extended.fasta",
        "rep_lib": subtype_dir / "TIR_fallback_rep_lib.fasta",
    }

    core_bed = subtype_dir / "TPase_regions.bed"
    make_core_bed(anchors, str(core_bed))

    try:
        success = extract_sequences_with_seqkit(
            genome_fasta,
            str(core_bed),
            str(output_files["core_fasta"]),
            anchors,
        )
        if not success:
            print(f"  Warning: failed to extract TPase regions for {subtype}", file=sys.stderr)
            return None, None, 0

        prime_success = extract_prime_sequences(
            genome_fasta,
            anchors,
            all_features,
            flank_size,
            str(output_files["5prime_fasta"]),
            str(output_files["3prime_fasta"]),
            mask_features,
            seq_lengths,
        )
        if not prime_success:
            print(f"  Warning: failed to extract prime regions for {subtype}", file=sys.stderr)

        run_prime_alignments(
            subtype_dir,
            threads=threads,
            min_num_alignments=min_num_alignments,
            verbose=verbose,
        )

        alignment_lengths = load_prime_alignment_lengths(subtype_dir)
        tir_elements = create_tir_elements(anchors, alignment_lengths)
        write_tir_gff(anchors, tir_elements, str(output_files["gff_out"]))

        extended_success = extract_extended_tir_regions(
            genome_fasta,
            tir_elements,
            str(output_files["extended_fasta"]),
        )
        if not extended_success:
            print(f"  Warning: failed to extract extended regions for {subtype}", file=sys.stderr)
            return output_files["gff_out"], None, len(tir_elements)

        # Append #classification suffix to extended-FASTA headers so that
        # reduce_library_size.R can group these sequences correctly. Use the
        # canonical slash form "Class_II/Subclass_1/TIR/<subtype>" so that
        # similarity-based hits line up with structure-based DANTE annotations
        # (clean_DANTE_names.R converts Class_II|Subclass_1|TIR|hAT to
        # Class_II/Subclass_1/TIR/hAT). Primary DANTE_TIR_final.fasta headers
        # are normalized to the same slash form by make_tir_combined_library.
        classification_with_slashes = f"Class_II/Subclass_1/TIR/{subtype_slug}"
        annotate_extended_fasta_headers(
            str(output_files["extended_fasta"]),
            classification_with_slashes,
        )

        clustering_success = run_mmseqs_clustering(
            str(output_files["extended_fasta"]),
            subtype_dir,
            threads=threads,
        )
        if not clustering_success:
            print(f"  Warning: mmseqs clustering failed for {subtype}", file=sys.stderr)
            return output_files["gff_out"], None, len(tir_elements)

        mmseqs_rep_fasta = subtype_dir / "mmseqs" / "cluster_rep_seq.fasta"
        if not mmseqs_rep_fasta.exists():
            print(f"  Warning: cluster representatives missing for {subtype}", file=sys.stderr)
            return output_files["gff_out"], None, len(tir_elements)

        with open(output_files["rep_lib"], "w") as rep_out:
            append_rep_library(mmseqs_rep_fasta, rep_out, classification_with_slashes)

        return output_files["gff_out"], output_files["rep_lib"], len(tir_elements)
    finally:
        if os.path.exists(core_bed):
            os.remove(core_bed)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fallback extraction of TIR regions anchored on TPase hits",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
The workflow processes each subtype independently. Only features with
Final_Classification=Class_II|Subclass_1|TIR|* and Name=TPase are used.
        """,
    )

    parser.add_argument("-g", "--genome", required=True, help="Input genome FASTA file")
    parser.add_argument("-a", "--annotations", required=True, help="Input GFF3 file with DANTE annotations")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory for result files")
    parser.add_argument(
        "-f",
        "--flank",
        type=int,
        default=10000,
        help="Flanking region size in bp (default: 10000)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of threads for alignment analysis (default: 1)",
    )
    parser.add_argument(
        "--min-num-alignments",
        type=int,
        default=3,
        help="Minimum number of alignments for length threshold calculation (default: 3)",
    )
    parser.add_argument("--mask-gff3", help="Optional GFF3 file with features that can limit flanking regions")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print verbose progress messages for alignment analysis",
    )

    args = parser.parse_args()

    for file_path, name in [(args.genome, "Genome"), (args.annotations, "Annotations")]:
        if not Path(file_path).exists():
            print(f"Error: {name} file {file_path} not found", file=sys.stderr)
            sys.exit(1)

    print(f"Parsing GFF3 file: {args.annotations}")
    all_features = parse_gff3_features(args.annotations)
    print(f"Found {len(all_features)} total features")

    mask_features = None
    if args.mask_gff3:
        if not Path(args.mask_gff3).exists():
            print(f"Error: Mask GFF3 file {args.mask_gff3} not found", file=sys.stderr)
            sys.exit(1)
        print(f"Parsing mask GFF3 file: {args.mask_gff3}")
        mask_features = parse_gff3_features(args.mask_gff3)
        print(f"Found {len(mask_features)} mask features")

    print("Filtering TIR TPase features...")
    tir_anchors_by_subtype = build_tir_anchors(all_features)
    if not tir_anchors_by_subtype:
        print("No features found with Name=TPase and Final_Classification=Class_II|Subclass_1|TIR|*")
        print("Creating empty output files and exiting.")
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        for fname in ["DANTE_TIR_FALLBACK.gff3", "TIR_fallback_rep_lib.fasta", "TIR_fallback_extended.fasta"]:
            Path(output_dir / fname).touch()
        sys.exit(0)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")

    combined_gff = output_dir / "DANTE_TIR_FALLBACK.gff3"
    combined_rep_lib = output_dir / "TIR_fallback_rep_lib.fasta"
    combined_extended = output_dir / "TIR_fallback_extended.fasta"

    with open(combined_gff, "w") as gff_handle:
        gff_handle.write("##gff-version 3\n")
        gff_handle.write(f"##{SOURCE_NAME.lower()} version 0.1\n")

    with open(combined_rep_lib, "w"):
        pass
    with open(combined_extended, "w"):
        pass

    print("Loading sequence lengths from genome index...")
    seq_lengths = load_sequence_lengths(args.genome)
    print(f"  Loaded lengths for {len(seq_lengths)} sequences")

    total_elements = 0
    processed_subtypes = []

    for subtype in sorted(tir_anchors_by_subtype.keys()):
        anchors = tir_anchors_by_subtype[subtype]
        gff_path, rep_path, n_elements = process_subtype(
            subtype=subtype,
            anchors=anchors,
            genome_fasta=args.genome,
            all_features=all_features,
            output_dir=output_dir,
            flank_size=args.flank,
            threads=args.threads,
            min_num_alignments=args.min_num_alignments,
            mask_features=mask_features,
            seq_lengths=seq_lengths,
            verbose=args.verbose,
        )

        if gff_path and Path(gff_path).exists():
            with open(combined_gff, "a") as gff_handle:
                append_gff_body(Path(gff_path), gff_handle)

        if rep_path and Path(rep_path).exists():
            with open(combined_rep_lib, "a") as rep_handle:
                with open(rep_path, "r") as rep_in:
                    rep_handle.write(rep_in.read())

        subtype_slug = sanitize_label(subtype)
        extended_fasta = output_dir / subtype_slug / "TPase_regions_extended.fasta"
        if extended_fasta.exists() and extended_fasta.stat().st_size > 0:
            with open(combined_extended, "a") as ext_handle:
                with open(extended_fasta, "r") as ext_in:
                    ext_handle.write(ext_in.read())

        if gff_path:
            processed_subtypes.append(subtype)
            total_elements += n_elements

    print("\nAll processing completed successfully!")
    print(f"  Subtypes processed: {len(processed_subtypes)}")
    print(f"  Total fallback elements: {total_elements}")
    print(f"  Combined GFF: {combined_gff}")
    print(f"  Combined repeat library: {combined_rep_lib}")
    print(f"  Combined extended FASTA: {combined_extended}")


if __name__ == "__main__":
    main()

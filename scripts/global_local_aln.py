#!/usr/bin/env python3
import argparse, sys
import parasail
from concurrent.futures import ThreadPoolExecutor, as_completed

def read_first_fasta_seq(path):
    seq, seen = [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if seen: break
                seen = True
                continue
            seq.append(line.strip())
    s = ''.join(seq).upper()
    if not s:
        raise SystemExit(f"No sequence found in {path}")
    return s


def read_fasta_sequences(path):
    """
    Read all sequences from a FASTA file.

    Parameters
    ----------
    path : str
        Path to FASTA file

    Returns
    -------
    list[tuple]
        List of (sequence_id, sequence) tuples
    """
    sequences = []
    current_id = None
    current_seq = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_id is not None and current_seq:
                    sequences.append((current_id, ''.join(current_seq).upper()))

                # Start new sequence
                current_id = line[1:].split()[0]  # Take only the first part of header
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_id is not None and current_seq:
            sequences.append((current_id, ''.join(current_seq).upper()))

    if not sequences:
        raise SystemExit(f"No sequences found in {path}")

    return sequences

def pick_func(end_anchor: str):
    # Fix 5′ (left) on both -> free 3′ ends on both
    # Fix 3′ (right) on both -> free 5′ ends on both
    if end_anchor == "5":
        return parasail.sg_qe_de_trace_striped_16
    else:
        return parasail.sg_qb_db_trace_striped_16
    # See Parasail naming table. :contentReference[oaicite:3]{index=3}

def make_matrix(match: int, mismatch: int):
    # Build a simple DNA matrix; strong negative mismatch stops extension.
    # (Alternative is parasail.nuc44, but simple is clearer/tunable.)
    return parasail.matrix_create("ACGT", match, mismatch)

def softclip_cigar(cigar_str, q_beg, q_end, q_len):
    left = q_beg
    right = q_len - (q_end + 1)
    parts = []
    if left > 0:
        parts.append(f"{left}S")
    parts.append(cigar_str)
    if right > 0:
        parts.append(f"{right}S")
    return ''.join(parts)

def per_column_scores(
    result,
    match=2,
    mismatch=-6,
    gap_open=14,
    gap_extend=3,
    free_side="right",   # "right" for sg_qe_de (fixed 5′); "left" for sg_qb_db (fixed 3′)
):
    """
    Compute per-column score contributions for a Parasail traceback (semi-global with one free side).

    Parameters
    ----------
    result : parasail.Result
        Returned by sg_qe_de_trace_striped_16 or sg_qb_db_trace_striped_16 (with _trace_*).
    match, mismatch : int
        Simple DNA match/mismatch scores for aligned letters.
    gap_open, gap_extend : int
        Affine gap penalties. Opening a gap of length L costs (gap_open + L*gap_extend).
    free_side : {"left","right"}
        Which end is *free* (unpenalized) in this semi-global setup:
        - "right": both 3′ ends are free (use with sg_qe_de… → fixed 5′).
        - "left" : both 5′ ends are free (use with sg_qb_db… → fixed 3′).

    Returns
    -------
    scores : list[int]
        One integer per printed alignment column (including gap columns).
        Terminal gap columns on the *free* side contribute 0 by definition of free-end gaps.
    """
    tb = result.traceback
    q = tb.query  # alignment string with '-' for gaps
    r = tb.ref

    n = len(q)
    assert n == len(r), "Traceback strings must have equal length"

    # Identify terminal gap columns on the *free* side → zero them out
    left_trim = 0
    right_trim_from = n

    if free_side == "left":
        # Count leading columns where either sequence has a gap
        i = 0
        while i < n and (q[i] == '-' or r[i] == '-'):
            i += 1
        left_trim = i
    elif free_side == "right":
        # Count trailing columns where either sequence has a gap
        j = n - 1
        while j >= 0 and (q[j] == '-' or r[j] == '-'):
            j -= 1
        right_trim_from = j + 1
    else:
        raise ValueError("free_side must be 'left' or 'right'")

    scores = [0] * n

    # Walk the alignment once and compute contributions with affine gaps.
    # We track whether we're inside a gap run in query (I) or ref (D).
    in_gap_q = False  # gap in query (deletion relative to ref → q=='-')
    in_gap_r = False  # gap in ref (insertion relative to ref → r=='-')

    for i in range(n):
        qc = q[i]
        rc = r[i]

        # Terminal gap columns on the FREE side are unpenalized by definition
        if free_side == "left" and i < left_trim:
            scores[i] = 0
            continue
        if free_side == "right" and i >= right_trim_from:
            scores[i] = 0
            continue

        if qc != '-' and rc != '-':
            # Closing any ongoing gaps
            in_gap_q = False
            in_gap_r = False
            # Match/mismatch
            scores[i] = match if qc == rc else mismatch
        elif qc == '-' and rc != '-':
            # Gap in query (D). Apply open at first, extend otherwise.
            if not in_gap_q:
                scores[i] = -(gap_open + gap_extend)
                in_gap_q = True
            else:
                scores[i] = -gap_extend
        elif qc != '-' and rc == '-':
            # Gap in ref (I)
            if not in_gap_r:
                scores[i] = -(gap_open + gap_extend)
                in_gap_r = True
            else:
                scores[i] = -gap_extend
        else:
            # Both '-' should not occur in a valid alignment
            raise ValueError("Invalid alignment column: both are gaps")

    return scores


def per_column_scores_alt(
    result,
    match=2,
    mismatch=-6,
    gap_open=14,
    gap_extend=3,
    cumulative=False
):
    """
    Alternative function to compute per-column score contributions using the comp string format.

    Parameters
    ----------
    result : parasail.Result
        Returned by alignment function with traceback information.
    match : int
        Score for matching positions (default 2).
    mismatch : int
        Score for mismatching positions (default -6).
    gap_open : int
        Gap opening penalty (default 14).
    gap_extend : int
        Gap extension penalty (default 3).
    cumulative : bool, optional
        If True, return cumulative scores instead of per-column scores (default False).

    Returns
    -------
    scores : list[int]
        One integer per alignment column based on comp string:
        - '|' = match (positive score)
        - '.' = mismatch (negative score)
        - ' ' = gap (gap_open for first gap, gap_extend for extensions)
    """
    tb = result.traceback
    comp = tb.comp  # comparison string: "|" = match, "." = mismatch, " " = gap
    query = tb.query  # alignment string with '-' for gaps
    ref = tb.ref      # reference alignment string with '-' for gaps

    n = len(comp)
    assert n == len(query) == len(ref), "Traceback strings must have equal length"

    scores = [0] * n

    # Track gap states to handle gap opening vs extension
    in_gap_query = False    # True when query has gap (query[i] == '-')
    in_gap_ref = False      # True when ref has gap (ref[i] == '-')

    for i in range(n):
        comp_char = comp[i]
        query_char = query[i]
        ref_char = ref[i]

        if comp_char == '|':
            # Match - reset gap states
            in_gap_query = False
            in_gap_ref = False
            scores[i] = match

        elif comp_char == '.':
            # Mismatch - reset gap states
            in_gap_query = False
            in_gap_ref = False
            scores[i] = mismatch

        elif comp_char == ' ':
            # Gap - determine if in query or reference
            if query_char == '-' and ref_char != '-':
                # Gap in query (deletion from query perspective)
                if not in_gap_query:
                    # Gap opening
                    scores[i] = -(gap_open + gap_extend)
                    in_gap_query = True
                else:
                    # Gap extension
                    scores[i] = -gap_extend
                # Reset ref gap state
                in_gap_ref = False

            elif query_char != '-' and ref_char == '-':
                # Gap in reference (insertion from query perspective)
                if not in_gap_ref:
                    # Gap opening
                    scores[i] = -(gap_open + gap_extend)
                    in_gap_ref = True
                else:
                    # Gap extension
                    scores[i] = -gap_extend
                # Reset query gap state
                in_gap_query = False

            else:
                # Both gaps or neither - shouldn't happen in valid alignment
                raise ValueError(f"Invalid alignment at position {i}: query='{query_char}', ref='{ref_char}', comp='{comp_char}'")
        else:
            # Unknown comparison character
            raise ValueError(f"Unknown comparison character '{comp_char}' at position {i}")

    if cumulative:
        for i in range(1, n):
            scores[i] += scores[i-1]
    return scores


def extract_optimal_alignment(result, scores, end="5"):
    """
    Find the position with maximum cumulative score and extract trimmed sequences.

    Parameters
    ----------
    result : parasail.Result
        Returned by alignment function with traceback information.
    scores : list[int]
        Per-column scores (not cumulative) from per_column_scores_alt or similar.
    end : str, optional
        Which end is fixed ("5" or "3"). Affects direction of cumulative scoring and trimming.
        - "5": Cumulative scoring from left (5' end), trim right side to max position
        - "3": Cumulative scoring from right (3' end), trim left side to max position

    Returns
    -------
    dict
        Dictionary containing:
        - 'max_pos': Position of maximum cumulative score (0-based)
        - 'max_score': Maximum cumulative score value
        - 'trimmed_query': Query alignment string trimmed to max position
        - 'trimmed_ref': Reference alignment string trimmed to max position
        - 'degapped_query': Trimmed query with gaps removed
        - 'degapped_ref': Trimmed reference with gaps removed
        - 'degapped_query_len': Length of degapped query
        - 'degapped_ref_len': Length of degapped reference
        - 'end_mode': Which end was used for scoring direction
    """
    # Extract traceback sequences
    tb = result.traceback
    query = tb.query
    ref = tb.ref
    n = len(scores)

    if end == "5":
        # 5' end fixed: cumulative scoring from left (5' -> 3')
        cumulative_scores = []
        cumsum = 0
        for score in scores:
            cumsum += score
            cumulative_scores.append(cumsum)

        # Find position of maximum cumulative score
        max_score = max(cumulative_scores)
        max_pos = cumulative_scores.index(max_score)

        # Trim sequences from start to maximum score position (inclusive)
        trimmed_query = query[:max_pos + 1]
        trimmed_ref = ref[:max_pos + 1]

    elif end == "3":
        # 3' end fixed: cumulative scoring from right (3' -> 5')
        # Reverse the scores and calculate cumulative from the end
        reversed_scores = scores[::-1]
        cumulative_scores = []
        cumsum = 0
        for score in reversed_scores:
            cumsum += score
            cumulative_scores.append(cumsum)

        # Find position of maximum cumulative score in reversed array
        max_score = max(cumulative_scores)
        max_pos_reversed = cumulative_scores.index(max_score)

        # Convert back to original array position
        max_pos = n - 1 - max_pos_reversed

        # Trim sequences from maximum score position to end (inclusive)
        trimmed_query = query[max_pos:]
        trimmed_ref = ref[max_pos:]

    else:
        raise ValueError("end parameter must be '5' or '3'")

    # Create degapped versions by removing gaps
    degapped_query = trimmed_query.replace('-', '')
    degapped_ref = trimmed_ref.replace('-', '')

    # Calculate lengths
    degapped_query_len = len(degapped_query)
    degapped_ref_len = len(degapped_ref)

    return {
        'max_pos': max_pos,
        'max_score': max_score,
        'trimmed_query': trimmed_query,
        'trimmed_ref': trimmed_ref,
        'degapped_query': degapped_query,
        'degapped_ref': degapped_ref,
        'degapped_query_len': degapped_query_len,
        'degapped_ref_len': degapped_ref_len,
        'end_mode': end
    }


def run_alignment_comparison(seq1_id, seq1, seq2_id, seq2, args):
    """
    Run alignment between two sequences and return results.

    Parameters
    ----------
    seq1_id : str
        ID of first sequence
    seq1 : str
        First sequence
    seq2_id : str
        ID of second sequence
    seq2 : str
        Second sequence
    args : argparse.Namespace
        Command line arguments

    Returns
    -------
    dict
        Alignment results including optimal alignment data
    """
    func = pick_func(args.end)
    mat = make_matrix(args.match, args.mismatch)

    res = func(seq1, seq2, args.open, args.extend, mat)

    scores = per_column_scores_alt(
        res,
        match=args.match,
        mismatch=args.mismatch,
        gap_open=args.open,
        gap_extend=args.extend,
        cumulative=False
    )

    result_data = extract_optimal_alignment(res, scores, args.end)

    # Add sequence IDs and original alignment score
    result_data['query_id'] = seq1_id
    result_data['ref_id'] = seq2_id
    result_data['parasail_score'] = res.score
    result_data['query_len'] = len(seq1)
    result_data['ref_len'] = len(seq2)

    return result_data


def generate_sequence_pairs(sequences):
    """
    Generate all unique sequence pairs for comparison.

    Parameters
    ----------
    sequences : list[tuple]
        List of (sequence_id, sequence) tuples

    Yields
    ------
    tuple
        (i, j, seq1_id, seq1, seq2_id, seq2) for each unique pair
    """
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            seq1_id, seq1 = sequences[i]
            seq2_id, seq2 = sequences[j]
            yield (i, j, seq1_id, seq1, seq2_id, seq2)


def run_single_comparison(pair_data, args):
    """
    Worker function for running a single pairwise comparison.

    Parameters
    ----------
    pair_data : tuple
        (i, j, seq1_id, seq1, seq2_id, seq2) from generate_sequence_pairs
    args : argparse.Namespace
        Command line arguments

    Returns
    -------
    dict or None
        Alignment result if above threshold, None otherwise
    """
    i, j, seq1_id, seq1, seq2_id, seq2 = pair_data

    # Run alignment
    result = run_alignment_comparison(seq1_id, seq1, seq2_id, seq2, args)

    # Apply score threshold filter
    if result['max_score'] >= args.score_threshold:
        return result
    else:
        return None


def write_results_table(results, output_file):
    """
    Write alignment results to tab-delimited file.

    Parameters
    ----------
    results : list[dict]
        List of alignment result dictionaries
    output_file : str
        Output file path
    """
    if not results:
        return

    # Define column order
    columns = [
        'query_id', 'ref_id', 'query_len', 'ref_len',
        'parasail_score', 'max_score', 'max_pos', 'end_mode',
        'degapped_query_len', 'degapped_ref_len',
        'trimmed_query', 'trimmed_ref',
        'degapped_query', 'degapped_ref'
    ]

    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join(columns) + '\n')

        # Write data rows
        for result in results:
            row = [str(result.get(col, '')) for col in columns]
            f.write('\t'.join(row) + '\n')


def extract_end_region(fasta_file, end, region_length=30, output_file=None):
    """
    Extract a fixed-length region from the 5' or 3' end of sequences.

    Parameters
    ----------
    fasta_file : str
        Path to input FASTA file
    end : str
        Which end to extract ("5" or "3")
    region_length : int, optional
        Length of region to extract, default 30
    output_file : str, optional
        Output FASTA file. If None, creates a temporary file

    Returns
    -------
    str
        Path to output FASTA file with extracted regions
    """
    import tempfile

    if output_file is None:
        output_file = tempfile.mktemp(suffix='_end_region.fasta')

    sequences = read_fasta_sequences(fasta_file)

    with open(output_file, 'w') as f:
        for seq_id, seq in sequences:
            if end == "5":
                # Extract first region_length nucleotides
                region = seq[:region_length] if len(seq) >= region_length else seq
            else:  # end == "3"
                # Extract last region_length nucleotides
                region = seq[-region_length:] if len(seq) >= region_length else seq

            f.write(f">{seq_id}\n{region}\n")

    return output_file


def run_mmseqs_prefilter(fasta_file, end, min_identity=0.8, threads=1, verbose=True):
    """
    Run MMseqs2 easy-search to identify candidate pairs for alignment.

    Uses a 30nt region from the correct end (matching the fixed end in alignment).

    Parameters
    ----------
    fasta_file : str
        Path to input FASTA file
    end : str
        Which end is fixed ("5" or "3")
    min_identity : float, optional
        Minimum sequence identity (0-1), default 0.8 (80%)
    threads : int, optional
        Number of threads, default 1
    verbose : bool, optional
        Print progress messages, default True

    Returns
    -------
    set of tuple
        Set of (seq1_id, seq2_id) pairs that passed prefiltering
    """
    import subprocess
    import tempfile
    import os

    # Extract 30nt from the correct end
    if verbose:
        print(f"  Extracting 30nt from {end}' end for MMseqs2 prefiltering...")

    end_region_file = extract_end_region(fasta_file, end, region_length=30)

    # Create temporary directory for MMseqs2
    tmp_dir = tempfile.mkdtemp(prefix='mmseqs_prefilter_')
    output_file = os.path.join(tmp_dir, 'search_results.tsv')

    try:
        # Run MMseqs2 easy-search
        # Format: query target identity alnlen mismatch gapopen qstart qend tstart tend evalue bits
        if verbose:
            print(f"  Running MMseqs2 easy-search with {min_identity*100:.0f}% identity threshold...")

        cmd = [
            'mmseqs', 'easy-search',
            end_region_file,
            end_region_file,  # Search against itself
            output_file,
            tmp_dir,
            '--min-seq-id', str(min_identity),
            '--threads', str(threads),
            '--format-output', 'query,target',
            '--search-type', '3'   #Nucleotide search
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # Parse results to get candidate pairs
        candidate_pairs = set()
        with open(output_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    query_id = fields[0]
                    target_id = fields[1]

                    # Skip self-hits and ensure consistent ordering
                    if query_id != target_id:
                        pair = tuple(sorted([query_id, target_id]))
                        candidate_pairs.add(pair)

        if verbose:
            print(f"  MMseqs2 prefiltering identified {len(candidate_pairs)} candidate pairs")

        return candidate_pairs

    except subprocess.CalledProcessError as e:
        print(f"Warning: MMseqs2 prefiltering failed: {e}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        print("Falling back to all-vs-all comparison without prefiltering", file=sys.stderr)
        return None
    except FileNotFoundError:
        print("Warning: mmseqs not found. Falling back to all-vs-all comparison without prefiltering",
              file=sys.stderr)
        return None
    finally:
        # Clean up temporary files
        try:
            os.remove(end_region_file)
        except:
            pass
        try:
            import shutil
            shutil.rmtree(tmp_dir)
        except:
            pass


def generate_sequence_pairs_filtered(sequences, candidate_pairs):
    """
    Generate sequence pairs from a filtered set of candidates.

    Parameters
    ----------
    sequences : list of tuple
        List of (seq_id, sequence) tuples
    candidate_pairs : set of tuple
        Set of (seq1_id, seq2_id) pairs to process

    Yields
    ------
    tuple
        (seq1_id, seq1, seq2_id, seq2) for each candidate pair
    """
    # Create lookup dict for fast access
    seq_dict = {seq_id: seq for seq_id, seq in sequences}

    for seq1_id, seq2_id in candidate_pairs:
        if seq1_id in seq_dict and seq2_id in seq_dict:
            yield (seq1_id, seq_dict[seq1_id], seq2_id, seq_dict[seq2_id])


def run_all_vs_all_alignment(
    fasta_file,
    output_file,
    end="5",
    gap_open=12,
    gap_extend=3,
    match=2,
    mismatch=-2,
    score_threshold=20,
    threads=1,
    verbose=True,
    use_prefilter=True,
    prefilter_identity=0.8
):
    """
    Run all-vs-all pairwise alignment analysis on sequences in a FASTA file.

    Parameters
    ----------
    fasta_file : str
        Path to input FASTA file with multiple sequences
    output_file : str
        Path to output tab-delimited file with alignment results
    end : str, optional
        Which end is FIXED for BOTH sequences ("5" or "3"), default "5"
    gap_open : int, optional
        Gap opening penalty, default 12
    gap_extend : int, optional
        Gap extension penalty, default 3
    match : int, optional
        Match score, default 2
    mismatch : int, optional
        Mismatch score, default -2
    score_threshold : int, optional
        Minimum max_score for alignment to be reported, default 20
    threads : int, optional
        Number of threads for parallel processing, default 1
    verbose : bool, optional
        Print progress messages, default True
    use_prefilter : bool, optional
        Use MMseqs2 prefiltering to speed up alignment, default True
    prefilter_identity : float, optional
        Minimum identity threshold for MMseqs2 prefiltering (0-1), default 0.8 (80%)

    Returns
    -------
    list[dict]
        List of alignment result dictionaries
    """
    # Create a namespace object to mimic argparse args
    class Args:
        pass

    args = Args()
    args.fasta = fasta_file
    args.output = output_file
    args.end = end
    args.open = gap_open
    args.extend = gap_extend
    args.match = match
    args.mismatch = mismatch
    args.score_threshold = score_threshold
    args.threads = threads

    # Read all sequences from FASTA file
    if verbose:
        print(f"Reading sequences from {fasta_file}...")
    sequences = read_fasta_sequences(fasta_file)
    if verbose:
        print(f"Found {len(sequences)} sequences")

    # Calculate total possible pairs
    total_possible_pairs = len(sequences) * (len(sequences) - 1) // 2
    if verbose:
        print(f"Total possible pairwise comparisons: {total_possible_pairs}")

    # Run MMseqs2 prefiltering if enabled and sequence count is sufficient
    candidate_pairs = None
    if use_prefilter:
        if len(sequences) < 100:
            # Skip prefiltering for small datasets
            if verbose:
                print(f"\nSkipping prefiltering (fewer than 100 sequences)")
            use_prefilter = False
        else:
            if verbose:
                print(f"\nRunning MMseqs2 prefiltering...")
            candidate_pairs = run_mmseqs_prefilter(
                fasta_file,
                end,
                min_identity=prefilter_identity,
                threads=threads,
                verbose=verbose
            )

            # If prefiltering failed, fall back to all-vs-all
            if candidate_pairs is None:
                use_prefilter = False
            elif verbose:
                # Report filtering statistics
                num_filtered = len(candidate_pairs)
                reduction_pct = (1 - num_filtered / total_possible_pairs) * 100 if total_possible_pairs > 0 else 0
                print(f"  Prefiltering reduced comparisons by {reduction_pct:.1f}% ({total_possible_pairs} → {num_filtered})")

    # Determine comparison strategy
    if use_prefilter and candidate_pairs is not None:
        # Use prefiltered pairs
        total_comparisons = len(candidate_pairs)
        if verbose:
            print(f"\nPerforming {total_comparisons} prefiltered pairwise alignments using {threads} threads...")
        pair_generator = generate_sequence_pairs_filtered(sequences, candidate_pairs)
    else:
        # Use all-vs-all pairs
        total_comparisons = len(sequences) * (len(sequences) - 1) // 2
        if verbose:
            print(f"\nPerforming {total_comparisons} pairwise alignments using {threads} threads...")
        pair_generator = generate_sequence_pairs(sequences)

    # Perform comparisons
    results = []

    if threads == 1:
        # Serial processing for single thread
        comparison_count = 0
        for pair_data in pair_generator:
            comparison_count += 1
            if verbose and (comparison_count % 100 == 0 or comparison_count == total_comparisons):
                print(f"  Progress: {comparison_count}/{total_comparisons} comparisons")

            # Handle different tuple formats from generators
            if len(pair_data) == 6:
                # From generate_sequence_pairs: (i, j, seq1_id, seq1, seq2_id, seq2)
                result = run_single_comparison(pair_data, args)
            else:
                # From generate_sequence_pairs_filtered: (seq1_id, seq1, seq2_id, seq2)
                seq1_id, seq1, seq2_id, seq2 = pair_data
                result = run_alignment_comparison(seq1_id, seq1, seq2_id, seq2, args)
                # Apply score threshold
                if result['max_score'] < args.score_threshold:
                    result = None

            if result is not None:
                results.append(result)
    else:
        # Parallel processing
        comparison_count = 0
        with ThreadPoolExecutor(max_workers=threads) as executor:
            # Submit all tasks - need to handle different tuple formats
            futures = []
            for pair_data in pair_generator:
                if len(pair_data) == 6:
                    # From generate_sequence_pairs
                    future = executor.submit(run_single_comparison, pair_data, args)
                else:
                    # From generate_sequence_pairs_filtered
                    seq1_id, seq1, seq2_id, seq2 = pair_data
                    # Create a wrapper function that applies threshold
                    def align_and_filter(s1_id, s1, s2_id, s2, a):
                        res = run_alignment_comparison(s1_id, s1, s2_id, s2, a)
                        return res if res['max_score'] >= a.score_threshold else None
                    future = executor.submit(align_and_filter, seq1_id, seq1, seq2_id, seq2, args)
                futures.append(future)

            # Collect results as they complete
            for future in as_completed(futures):
                comparison_count += 1
                if verbose and (comparison_count % 100 == 0 or comparison_count == total_comparisons):
                    print(f"  Progress: {comparison_count}/{total_comparisons} comparisons")

                result = future.result()
                if result is not None:
                    results.append(result)

    # Write results to output file
    if verbose:
        print(f"Writing results to {output_file}...")
    write_results_table(results, output_file)

    if verbose:
        print(f"Analysis complete! Generated {len(results)} alignment records from {total_comparisons} comparisons.")

    return results


def main():
    ap = argparse.ArgumentParser(
        description="All-to-all anchored semi-global alignment with Parasail (DNA). "
                    "Fix BOTH sequences at 5′ or 3′; the opposite ends are free."
    )
    ap.add_argument("fasta", help="Input FASTA file with multiple sequences")
    ap.add_argument("-o", "--output", required=True,
                    help="Output tab-delimited file with alignment results")
    ap.add_argument("--end", choices=["5","3"], required=True,
                    help="Which end is FIXED for BOTH sequences (5 or 3).")
    ap.add_argument("--open", type=int, default=12, help="Gap open (default 12)")
    ap.add_argument("--extend", type=int, default=3, help="Gap extend (default 3)")
    ap.add_argument("--match", type=int, default=2, help="Match score (default 2)")
    ap.add_argument("--mismatch", type=int, default=-2, help="Mismatch score (default -2)")
    ap.add_argument("--score-threshold", type=int, default=20,
                    help="Minimum max_score for alignment to be reported (default 20)")
    ap.add_argument("-t", "--threads", type=int, default=1,
                    help="Number of threads for parallel processing (default 1)")
    ap.add_argument("--no-prefilter", action="store_true",
                    help="Disable MMseqs2 prefiltering (use full all-vs-all comparison)")
    ap.add_argument("--prefilter-identity", type=float, default=0.8,
                    help="Minimum identity for MMseqs2 prefiltering (0-1, default 0.8)")
    ap.add_argument("-v", "--verbose", action="store_true",
                    help="Print progress messages (default: False)")
    args = ap.parse_args()

    if args.verbose:
        print(args)

    # Call the main function with parsed arguments
    run_all_vs_all_alignment(
        fasta_file=args.fasta,
        output_file=args.output,
        end=args.end,
        gap_open=args.open,
        gap_extend=args.extend,
        match=args.match,
        mismatch=args.mismatch,
        score_threshold=args.score_threshold,
        threads=args.threads,
        verbose=args.verbose,
        use_prefilter=not args.no_prefilter,
        prefilter_identity=args.prefilter_identity
    )


if __name__ == "__main__":
    main()

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
    args = ap.parse_args()
    print(args)

    # Read all sequences from FASTA file
    print(f"Reading sequences from {args.fasta}...")
    sequences = read_fasta_sequences(args.fasta)
    print(f"Found {len(sequences)} sequences")

    # Perform all-to-all comparisons
    results = []
    total_comparisons = len(sequences) * (len(sequences) - 1) // 2

    print(f"Performing {total_comparisons} pairwise alignments using {args.threads} threads...")

    if args.threads == 1:
        # Serial processing for single thread
        comparison_count = 0
        for pair_data in generate_sequence_pairs(sequences):
            comparison_count += 1
            if comparison_count % 100 == 0 or comparison_count == total_comparisons:
                print(f"  Progress: {comparison_count}/{total_comparisons} comparisons")

            result = run_single_comparison(pair_data, args)
            if result is not None:
                results.append(result)
    else:
        # Parallel processing
        comparison_count = 0
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            # Submit all tasks
            future_to_pair = {executor.submit(run_single_comparison, pair_data, args): pair_data
                            for pair_data in generate_sequence_pairs(sequences)}

            # Collect results as they complete
            for future in as_completed(future_to_pair):
                comparison_count += 1
                if comparison_count % 100 == 0 or comparison_count == total_comparisons:
                    print(f"  Progress: {comparison_count}/{total_comparisons} comparisons")

                result = future.result()
                if result is not None:
                    results.append(result)

    # Write results to output file
    print(f"Writing results to {args.output}...")
    write_results_table(results, args.output)

    print(f"Analysis complete! Generated {len(results)} alignment records from {total_comparisons} comparisons.")


if __name__ == "__main__":
    main()

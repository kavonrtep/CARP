#!/usr/bin/env python3
import argparse, sys
import parasail

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



def main():
    ap = argparse.ArgumentParser(
        description="Anchored semi-global alignment with Parasail (DNA). "
                    "Fix BOTH sequences at 5′ or 3′; the opposite ends are free."
    )
    ap.add_argument("fasta1")
    ap.add_argument("fasta2")
    ap.add_argument("--end", choices=["5","3"], required=True,
                    help="Which end is FIXED for BOTH sequences (5 or 3).")
    ap.add_argument("--open", type=int, default=12, help="Gap open (default 12)")
    ap.add_argument("--extend", type=int, default=3, help="Gap extend (default 2)")
    ap.add_argument("--match", type=int, default=2, help="Match score (default 2)")
    ap.add_argument("--mismatch", type=int, default=-2, help="Mismatch score (default -4)")
    ap.add_argument("--show-aln", action="store_true", help="Print only the aligned core strings")
    args = ap.parse_args()

    s1 = read_first_fasta_seq(args.fasta1)
    s2 = read_first_fasta_seq(args.fasta2)

    func = pick_func(args.end)
    mat = make_matrix(args.match, args.mismatch)

    res = func(s1, s2, args.open, args.extend, mat)

    scores = per_column_scores_alt(
        res,
        match=args.match,
        mismatch=args.mismatch,
        gap_open=args.open,
        gap_extend=args.extend,
        cumulative=False
    )
    print("Cumulative scores (alt method):")
    print(scores)
    result_data = extract_optimal_alignment(res, scores, args.end)
    print(f"Max cumulative score {result_data['max_score']} at position {result_data['max_pos']}")
    print("Trimmed alignment to max score position:")
    print(result_data['trimmed_query'])
    print(result_data['trimmed_ref'])
    print("Degapped trimmed sequences:")
    print(f">query_len={result_data['degapped_query_len']}")
    print(result_data['degapped_query'])
    print(f">ref_len={result_data['degapped_ref_len']}")
    print(result_data['degapped_ref'])


    exit()
    cumulative_scores = [sum(scores[:i+1]) for i in range(len(scores))]
    print("Cumulative scores:")
    print(cumulative_scores)

    print("Per-column scores:")
    print(args)
    print(scores)
    print(len(scores), "columns, total =", sum(scores))

    if res is None or res.cigar is None:
        sys.exit("Parasail returned no traceback/CIGAR; check inputs and parameters.")

    cig = res.cigar
    cigar_str = cig.decode.decode()  # SAM-style (M/=/X/I/D etc.)



    q_beg, r_beg = cig.beg_query, cig.beg_ref          # 0-based inclusive
    q_end, r_end = res.end_query, res.end_ref          # 0-based inclusive

    # Build soft-clipped CIGAR for the query so you can see hidden flanks
    cigar_q_soft = softclip_cigar(cigar_str, q_beg, q_end, len(s1))

    print("== Anchored semi-global alignment ==")
    print(f"Mode: end {args.end}′ fixed on BOTH "
          f"({'sg_qe_de' if args.end=='5' else 'sg_qb_db'})")  # :contentReference[oaicite:4]{index=4}
    print(f"Score: {res.score}")
    print(f"Query coords (1-based): {q_beg+1}-{q_end+1} / len={len(s1)}")
    print(f"Target coords (1-based): {r_beg+1}-{r_end+1} / len={len(s2)}")
    print(f"Query CIGAR (with soft-clips): {cigar_q_soft}")

    if args.show_aln:
        core_q = s1[q_beg:q_end+1]
        core_r = s2[r_beg:r_end+1]
        tb = res.traceback  # only the aligned core; no unaligned flanks
        print("\nAligned core (no flanks):")
        print(tb.query)  # aligned core with gaps
        print(tb.comp)
        print(tb.ref)
        print("\nCore slices (raw, ungapped):")
        print(">query_core")
        print(core_q)
        print(">target_core")
        print(core_r)

if __name__ == "__main__":
    main()

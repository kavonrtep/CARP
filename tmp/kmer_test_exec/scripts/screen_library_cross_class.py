#!/usr/bin/env python3
"""
Cross-class contamination screen for the RepeatMasker library.

Why
---
Every library reduction CARP runs is *within* a class: ``reduce_library_size.py``
clusters per classification and ``containment_reduce_library.py`` only drops a
fragment into a container of the SAME class. So a consensus that is part one
class and part another -- a LINE whose inferred boundary ran out into the
neighbouring Retand array, say -- is invisible to both: it looks like a
perfectly ordinary, slightly long LINE. RepeatMasker then masks the foreign part
genome-wide under the wrong label.

``filter_ltr_rt_library`` already does a special case of this (LTR library
screened against Class_II, whole sequence dropped on any hit). This generalises
it to every class and, instead of dropping, **trims**: the contaminated span is
removed and the rest of the consensus is kept, so a real element that merely
picked up a foreign tail keeps masking its own copies.

What counts as contamination
----------------------------
A blastn hit between two library consensi whose classifications diverge near the
root of the classification tree. "Near the root" is the important part: sibling
LTR lineages (``Ty1_copia/Ale`` vs ``Ty1_copia/Angela``) share real sequence
homology and must NOT be trimmed against each other, so two classes count as
incompatible only when their lowest common ancestor is at depth
``--max-shared-depth`` or shallower (default 1). With that default:

    Class_I/LINE            vs Class_I/LTR/.../Retand   LCA Class_I        (1)  CONFLICT
    Class_I/...             vs Class_II/...             LCA -              (0)  CONFLICT
    Class_I/LINE            vs rDNA/45S_rDNA            LCA -              (0)  CONFLICT
    Class_I/LTR/Ty1_copia   vs Class_I/LTR/Ty3_gypsy    LCA Class_I/LTR    (2)  ok
    Ty1_copia/Ale           vs Ty1_copia/Angela         LCA .../Ty1_copia  (3)  ok

An ancestor/descendant pair (``Class_I/LTR`` vs ``Class_I/LTR/Ty1_copia``) is
always compatible.

Which side is the chimera
-------------------------
A blastn hit is symmetric, so "these two consensi share sequence they should not"
does not by itself say which one is wrong. Trimming both would be worse than
trimming neither: on a wheat library the naive symmetric rule trimmed 1,114
Ty1_copia/Angela consensi against 86 LINE consensi -- i.e. it damaged the correct
library because of the broken one.

The region is assigned to whichever consensus it makes up a **larger fraction
of**. A chimeric LINE of 8.5 kb carrying 2 kb of Angela at its 3' end gives
``q_frac = 0.23`` against the Angela consensus's ``s_frac = 0.40``: the sequence
is more of what Angela *is* than of what the LINE is, so it is trimmed out of the
LINE and left in Angela. The same hit seen from the other direction reaches the
same verdict, so the decision is consistent whichever way round BLAST reports it.
Ties, and margins below ``--ownership-margin``, are left alone.

Truncate to the clean part
--------------------------
The consensus is cut down to its **longest span containing no conflicting
material at all**, and dropped if too little of it survives.

An earlier version only peeled conflict blocks inward from the two ends, on the
grounds that cutting an internal block fragments a sequence whose middle might be
contested for a legitimate reason (a nested insertion). On real data that leaves
most of the contamination in place: the foreign material in a chimeric consensus
is not one clean block glued to an end but a mosaic, because the genomic region
the boundary ran out into is itself a patchwork of decayed insertions with
unclassifiable rubble between them. Peeling stops at the first gap. Measured on a
wheat library, end-peeling removed 46% of the foreign base pairs it had already
identified and left 54% sitting inside consensi that RepeatMasker would then use.

Taking the longest clean span instead leaves **zero** identified foreign base
pairs. It costs more library sequence -- the clean fragments stranded between two
conflict blocks go too -- but a truncated consensus still masks its family, while
a chimeric one mislabels whatever the foreign part matches, genome-wide. On the
same library nothing fell below the drop thresholds (shortest survivor 829 bp).

Also reports (never trims on) consensi longer than their class plausibly allows,
per ``max_consensus_length`` in ``classification_vocabulary.yaml``.

Lossless fallback: if blastn is unavailable or fails, the library is copied
through unchanged with a warning. The screen is a correctness improvement, but a
missing aligner must not take the run down.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

AUDIT_HEADER = [
    "name", "classification", "original_length", "retained_length",
    "action", "conflict_bp", "conflict_classes", "detail",
]


def blast_id_to_name(seqid: str) -> str:
    """Strip the ``#<classification>`` suffix BLAST carries through as part of
    the sequence id.

    ``makeblastdb`` does not split a FASTA header on ``#``, so ``qseqid`` for
    ``>LINE_group_2147#Class_I/LINE`` comes back as the whole string. Looking the
    class up by the bare name against un-normalised BLAST ids silently matches
    nothing, and the screen then reports a clean library on a contaminated one.
    """
    return seqid.partition("#")[0]


def read_fasta(path):
    """[(name, classification, sequence)] preserving file order."""
    records = []
    name = cls = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, cls, "".join(chunks)))
                header = line[1:].strip()
                name, _, cls = header.partition("#")
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        records.append((name, cls, "".join(chunks)))
    return records


def lca_depth(a: str, b: str) -> int:
    """Number of leading path components two classifications share."""
    pa, pb = a.split("/"), b.split("/")
    n = 0
    for x, y in zip(pa, pb):
        if x != y:
            break
        n += 1
    return n


def is_conflict(a: str, b: str, max_shared_depth: int) -> bool:
    """True when a hit between classes `a` and `b` is cross-class contamination."""
    if not a or not b or a == b:
        return False
    pa, pb = a.split("/"), b.split("/")
    depth = lca_depth(a, b)
    # Ancestor/descendant is never a conflict: one is just a less specific
    # label for the same lineage.
    if depth == len(pa) or depth == len(pb):
        return False
    return depth <= max_shared_depth


def run_self_blast(library: str, workdir: Path, threads: int, evalue: str = "1e-10"):
    """makeblastdb + self blastn; return the outfmt-6 hits path, or None."""
    workdir.mkdir(parents=True, exist_ok=True)
    db = str(workdir / "libdb")
    hits = workdir / "cross_class_hits.tsv"
    try:
        subprocess.run(["makeblastdb", "-in", library, "-dbtype", "nucl",
                        "-out", db],
                       check=True, capture_output=True, text=True)
        with open(hits, "w") as out:
            subprocess.run(
                ["blastn", "-query", library, "-db", db,
                 "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend",
                 "-evalue", evalue, "-num_threads", str(threads)],
                check=True, stdout=out, stderr=subprocess.PIPE, text=True)
    except (subprocess.CalledProcessError, FileNotFoundError, OSError) as exc:
        sys.stderr.write(f"WARNING: blastn cross-class screen failed ({exc}); "
                         f"copying library through unscreened.\n")
        return None
    return hits


def merge_intervals(intervals):
    out = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return out


def longest_clean_span(length: int, conflicts):
    """Largest 1-based span of the consensus containing no conflict block.

    Returns ``(start, end, trimmed_5, trimmed_3)``; an empty span is reported as
    ``end < start``. `conflicts` must be merged and sorted.
    """
    best_start, best_end = 1, 0
    pos = 1
    for c_start, c_end in conflicts:
        if c_start - 1 >= pos and (c_start - 1 - pos) > (best_end - best_start):
            best_start, best_end = pos, c_start - 1
        pos = max(pos, c_end + 1)
    if length - pos > (best_end - best_start):
        best_start, best_end = pos, length
    return best_start, best_end, max(0, best_start - 1), max(0, length - best_end)


def load_max_consensus_length(vocabulary_path=None):
    """{classification prefix: max bp} from the vocabulary, or {} if absent."""
    try:
        from classification import load_vocabulary
        vocab = load_vocabulary(vocabulary_path)
        return dict(getattr(vocab, "max_consensus_length", {}) or {})
    except Exception as exc:  # vocabulary is advisory here, never fatal
        sys.stderr.write(f"NOTE: max_consensus_length unavailable ({exc}); "
                         f"length reporting disabled.\n")
        return {}


def max_length_for(cls: str, table) -> int:
    """Longest-prefix lookup of the class length bound; 0 = unbounded."""
    best = 0
    best_depth = -1
    for prefix, value in table.items():
        if cls == prefix or cls.startswith(prefix + "/"):
            depth = len(prefix.split("/"))
            if depth > best_depth:
                best, best_depth = value, depth
    return best


def screen(records, hits_path, opts):
    """Return (kept_records, audit_rows)."""
    by_class = {name: cls for name, cls, _ in records}
    lengths = {name: len(seq) for name, _, seq in records}
    conflicts = {}
    conflict_classes = {}

    if hits_path is not None:
        with open(hits_path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                qid, sid, pident, alen, qlen, slen, qstart, qend = parts[:8]
                qid, sid = blast_id_to_name(qid), blast_id_to_name(sid)
                if qid == sid:
                    continue
                if float(pident) < opts.min_identity or int(alen) < opts.min_length:
                    continue
                qcls, scls = by_class.get(qid, ""), by_class.get(sid, "")
                if not is_conflict(qcls, scls, opts.max_shared_depth):
                    continue
                # Ownership: the shared region belongs to whichever consensus it
                # makes up more of. Only the other one is contaminated by it.
                q_len_i, s_len_i = int(qlen), int(slen)
                if not q_len_i or not s_len_i:
                    continue
                q_frac = int(alen) / q_len_i
                s_frac = int(alen) / s_len_i
                if s_frac < q_frac + opts.ownership_margin:
                    continue
                a, b = int(qstart), int(qend)
                if a > b:
                    a, b = b, a
                conflicts.setdefault(qid, []).append((a, b))
                conflict_classes.setdefault(qid, set()).add(scls)

    max_len_table = load_max_consensus_length(opts.vocabulary)

    kept = []
    audit = []
    for name, cls, seq in records:
        length = len(seq)
        blocks = merge_intervals(conflicts.get(name, []))
        conflict_bp = sum(e - s + 1 for s, e in blocks)
        classes = ",".join(sorted(conflict_classes.get(name, ())))
        over_long = max_length_for(cls, max_len_table)
        notes = []
        if over_long and length > over_long:
            notes.append(f"exceeds max_consensus_length for {cls} ({length} > {over_long})")

        if not blocks:
            kept.append((name, cls, seq))
            if notes:
                audit.append([name, cls, length, length, "kept", 0, "", "; ".join(notes)])
            continue

        start, end, trim5, trim3 = longest_clean_span(length, blocks)
        retained = max(0, end - start + 1)

        if retained < opts.min_retained_length or (
                length and retained < opts.min_retained_fraction * length):
            audit.append([name, cls, length, 0, "dropped", conflict_bp, classes,
                          "; ".join(notes + ["retained sequence too short"])])
            continue

        kept.append((name, cls, seq[start - 1:end]))
        notes.append(f"trimmed {trim5} bp from 5' and {trim3} bp from 3'")
        audit.append([name, cls, length, retained, "trimmed", conflict_bp, classes,
                      "; ".join(notes)])

    return kept, audit


def write_fasta(records, path, width=60):
    tmp = str(path) + ".tmp"
    with open(tmp, "w") as fh:
        for name, cls, seq in records:
            fh.write(f">{name}#{cls}\n" if cls else f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")
    os.replace(tmp, path)


def write_audit(rows, path):
    tmp = str(path) + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(AUDIT_HEADER) + "\n")
        for row in sorted(rows, key=lambda r: str(r[0])):
            fh.write("\t".join(str(x) for x in row) + "\n")
    os.replace(tmp, path)


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Trim cross-class contamination out of the RepeatMasker library.")
    p.add_argument("-i", "--input", required=True, help="Input library FASTA")
    p.add_argument("-o", "--output", required=True, help="Output screened FASTA")
    p.add_argument("-a", "--audit", help="TSV of every decision (default: <output>.audit.tsv)")
    p.add_argument("-d", "--workdir", default="cross_class_workdir",
                   help="Directory for the blast database and hit table")
    p.add_argument("-t", "--threads", type=int, default=1)
    p.add_argument("--min-identity", type=float, default=80.0,
                   help="Minimum percent identity of a conflicting hit (default 80)")
    p.add_argument("--min-length", type=int, default=200,
                   help="Minimum alignment length of a conflicting hit (default 200)")
    p.add_argument("--max-shared-depth", type=int, default=1,
                   help="Two classes conflict only when their lowest common ancestor is "
                        "at this depth or shallower. 1 keeps sibling LTR lineages (which "
                        "share real homology) safe while catching LINE-vs-LTR and "
                        "Class_I-vs-Class_II (default 1)")
    p.add_argument("--ownership-margin", type=float, default=0.10,
                   help="A shared region is treated as belonging to the subject (and so "
                        "trimmed out of the query) only when it covers this much more of "
                        "the subject than of the query. Stops a chimera dragging its "
                        "correctly-built neighbour down with it (default 0.10)")
    p.add_argument("--min-retained-length", type=int, default=300,
                   help="Drop the consensus if fewer bp survive the trim (default 300)")
    p.add_argument("--min-retained-fraction", type=float, default=0.2,
                   help="Drop the consensus if less than this fraction survives (default 0.2)")
    p.add_argument("--vocabulary", help="classification_vocabulary.yaml (default: autodetect)")
    p.add_argument("--disabled", action="store_true",
                   help="Copy the library through unchanged (screen turned off in config)")
    args = p.parse_args(argv)

    audit_path = args.audit or (args.output + ".audit.tsv")

    if args.disabled:
        shutil.copyfile(args.input, args.output)
        write_audit([], audit_path)
        print("cross-class screen disabled; library copied through unchanged")
        return 0

    records = read_fasta(args.input)
    if not records:
        shutil.copyfile(args.input, args.output)
        write_audit([], audit_path)
        print("empty library; nothing to screen")
        return 0

    hits = run_self_blast(args.input, Path(args.workdir), args.threads)
    kept, audit = screen(records, hits, args)

    write_fasta(kept, args.output)
    write_audit(audit, audit_path)

    trimmed = sum(1 for r in audit if r[4] == "trimmed")
    dropped = sum(1 for r in audit if r[4] == "dropped")
    bp_in = sum(len(s) for _, _, s in records)
    bp_out = sum(len(s) for _, _, s in kept)
    print(f"cross-class screen: {len(records)} consensi in, {len(kept)} out "
          f"({trimmed} trimmed, {dropped} dropped); "
          f"{bp_in:,} -> {bp_out:,} bp ({100.0 * (bp_in - bp_out) / bp_in:.2f}% removed)")
    print(f"  audit: {audit_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

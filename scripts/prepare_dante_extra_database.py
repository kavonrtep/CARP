#!/usr/bin/env python3
"""Check a user protein database for ``dante -e`` (config dante_extra_database).

DANTE appends the file to REXdb and accepts almost anything: duplicate IDs, IDs
already used by REXdb, unknown domains, misspelled classes and DNA sequences
all finish with exit status 0 and change the annotation silently. This script
checks every record before DANTE runs and reports all errors at once.

Header format, as DANTE defines it (``DANTE_filtered.fasta`` already has it)::

    >ID DOMAIN CLASS
    >chr1:1000-1400 RT Class_I|LTR|Ty1/copia|Ale

* ID: unique, not used by REXdb, letters, digits and ``_ . : -`` only.
* DOMAIN: a REXdb domain type (RT, RH, INT, ...).
* CLASS: a REXdb class in REXdb spelling (``Ty1/copia``, not ``Ty1_copia``), or
  an internal node of the REXdb tree (``Class_I|LTR|Ty1/copia``). Internal nodes
  are accepted with a warning: DANTE classifies a domain by the common prefix of
  its good hits, so such a protein lowers calls that would otherwise reach a
  lineage.

The allowed classes, domains and reserved IDs are read from the installed REXdb,
so they always match the DANTE version in use. Standard library only (the DANTE
environment has no PyYAML).

Outputs: the cleaned FASTA (same records and order; CRLF, blank lines and extra
spaces removed) and a TSV summary that the HTML report reads.

Exit status: 0 = usable, 2 = records need a fix, 3 = a file cannot be read.
"""
import argparse
import difflib
import hashlib
import os
import re
import sys

ID_RE = re.compile(r"[A-Za-z0-9_.:-]+")
# lastdb aborts with "bad symbol in sequence" on '-', '/' and digits.
SEQ_RE = re.compile(r"[A-Za-z*.]+")
DNA_RE = re.compile(r"[ACGTUN]", re.I)
DNA_MIN_LENGTH = 50
DNA_MIN_FRACTION = 0.9
MAX_REPORTED = 50


class Reference:
    """Classes, internal nodes, domains and IDs of the installed REXdb."""

    def __init__(self, class_path, pdb_path):
        self.class_of = {}
        with open(class_path) as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) > 1:
                    self.class_of[fields[0]] = "|".join(fields[1:])
        self.classes = set(self.class_of.values())
        self.internal = {"|".join(c.split("|")[:i])
                         for c in self.classes
                         for i in range(1, c.count("|") + 1)} - self.classes
        self.domains = set()
        self.domains_of = {}          # class or internal node -> domain types
        with open(pdb_path) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                # >Lineage-DOMAIN__ID, split the way dante.py splits it
                prefix, _, seq_id = line[1:].split()[0].partition("__")
                domain = prefix.split("-")[1]
                self.domains.add(domain)
                parts = self.class_of.get(seq_id, "").split("|")
                for i in range(1, len(parts) + 1):
                    self.domains_of.setdefault(
                        "|".join(parts[:i]), set()).add(domain)
        if not self.classes or not self.domains:
            raise ValueError("no REXdb classes or domains found in "
                             f"{class_path} / {pdb_path}")


def read_fasta(path):
    """[(line number, header without '>', sequence)]; ValueError if not FASTA."""
    records = []
    with open(path, encoding="utf-8-sig") as fh:     # tolerate a BOM; CRLF -> LF
        for number, line in enumerate(fh, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                records.append((number, line[1:].strip(), []))
            elif records:
                records[-1][2].append("".join(line.split()))
            else:
                raise ValueError(f"line {number}: sequence before the first '>' "
                                 "header; this is not a FASTA file")
    return [(n, h, "".join(s)) for n, h, s in records]


def class_hint(value, ref):
    """What to tell the user about a CLASS that REXdb does not have."""
    def key(s):
        return re.sub(r"[/|_\s-]", "", s).lower()
    allowed = sorted(ref.classes | ref.internal)
    same = [c for c in allowed if key(c) == key(value)]
    if same:
        return f"REXdb spells it '{same[0]}'"
    parts = value.split("|")
    for i in range(len(parts) - 1, 0, -1):
        parent = "|".join(parts[:i])
        if parent in ref.classes:
            return (f"'{parent}' is the most detailed class here; new lineages "
                    "are not supported")
    close = difflib.get_close_matches(value, allowed, n=1, cutoff=0.6)
    if close:
        return f"did you mean '{close[0]}'?"
    return "it must be a class of the installed REXdb"


def is_dna(seq):
    return (len(seq) >= DNA_MIN_LENGTH
            and len(DNA_RE.findall(seq)) / len(seq) > DNA_MIN_FRACTION)


def check(records, ref):
    """Return (rows, errors, warnings); rows = (id, domain, class, sequence)."""
    rows, errors, warnings = [], [], []
    reserved = set(ref.class_of)
    seen = {}
    for number, header, seq in records:
        where = f"line {number}: >{header}"
        fields = header.split()
        if len(fields) != 3:
            errors.append(f"{where}: expected 3 fields, '>ID DOMAIN CLASS', "
                          f"found {len(fields)}")
            continue
        seq_id, domain, cls = fields

        if not ID_RE.fullmatch(seq_id):
            errors.append(f"{where}: ID may contain only letters, digits and _ . : -")
        stored = seq_id.replace("-", "_")      # DANTE stores the ID like this
        if stored in reserved:
            errors.append(f"{where}: ID is already used in REXdb")
        if stored in seen:
            errors.append(f"{where}: same ID as line {seen[stored]} "
                          "(DANTE treats '-' and '_' as equal)")
        else:
            seen[stored] = number

        if domain not in ref.domains:
            errors.append(f"{where}: unknown domain '{domain}'; "
                          f"REXdb has {', '.join(sorted(ref.domains))}")
        if cls in ref.internal:
            warnings.append(f"{where}: internal class '{cls}'; domains that match "
                            "this protein may lose their lineage")
        elif cls not in ref.classes:
            errors.append(f"{where}: unknown class '{cls}'; {class_hint(cls, ref)}")
        if (domain in ref.domains and cls in ref.domains_of
                and domain not in ref.domains_of[cls]):
            warnings.append(f"{where}: REXdb has no {domain} domain in '{cls}'")

        if not seq:
            errors.append(f"{where}: empty sequence")
        elif not SEQ_RE.fullmatch(seq):
            bad = "".join(sorted(set(re.sub(SEQ_RE, "", seq))))
            errors.append(f"{where}: sequence contains '{bad}'; use letters and "
                          "'*' only (remove gaps and frameshift marks)")
        elif is_dna(seq):
            errors.append(f"{where}: this looks like DNA; a protein is needed")
        rows.append((seq_id, domain, cls, seq))
    if not records:
        errors.append("the file has no sequences")
    return rows, errors, warnings


def write_atomic(path, text):
    tmp = f"{path}.tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(text)
    os.replace(tmp, path)


def fasta_text(rows):
    out = []
    for seq_id, domain, cls, seq in rows:
        out.append(f">{seq_id} {domain} {cls}\n")
        out.extend(f"{seq[i:i + 60]}\n" for i in range(0, len(seq), 60))
    return "".join(out)


def summary_text(rows, ref, source, rexdb):
    """'#key<TAB>value' lines, then one row per (class, domain), sorted."""
    with open(source, "rb") as fh:
        digest = hashlib.sha256(fh.read()).hexdigest()
    counts = {}
    for _, domain, cls, _ in rows:
        counts[(cls, domain)] = counts.get((cls, domain), 0) + 1
    out = [f"#source\t{os.path.basename(source)}\n", f"#sha256\t{digest}\n",
           f"#rexdb\t{rexdb}\n", f"#proteins\t{len(rows)}\n",
           "classification\tdomain\tproteins\tlevel\n"]
    for (cls, domain), n in sorted(counts.items()):
        level = "internal" if cls in ref.internal else "terminal"
        out.append(f"{cls}\t{domain}\t{n}\t{level}\n")
    return "".join(out)


def report(kind, messages):
    for message in messages[:MAX_REPORTED]:
        print(f"{kind}: {message}", file=sys.stderr)
    if len(messages) > MAX_REPORTED:
        print(f"{kind}: ... and {len(messages) - MAX_REPORTED} more",
              file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Check a protein FASTA for dante -e and write a clean copy.",
        epilog="Exit status: 0 usable, 2 records need a fix, 3 unreadable file.")
    ap.add_argument("input", help="protein FASTA, headers '>ID DOMAIN CLASS'")
    ap.add_argument("--reference-class", required=True, help="REXdb *_class file")
    ap.add_argument("--reference-pdb", required=True, help="REXdb *_pdb FASTA")
    ap.add_argument("-o", "--output", required=True, help="cleaned FASTA")
    ap.add_argument("--summary", required=True, help="summary TSV")
    a = ap.parse_args(argv)
    try:
        ref = Reference(a.reference_class, a.reference_pdb)
        records = read_fasta(a.input)
    except (OSError, ValueError) as exc:        # UnicodeDecodeError is a ValueError
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    rows, errors, warnings = check(records, ref)
    report("WARNING", warnings)
    report("ERROR", errors)
    if errors:
        print(f"FAILED: {len(errors)} error(s) in {a.input} "
              f"({len(records)} records)", file=sys.stderr)
        return 2
    rexdb = re.sub(r"_class$", "", os.path.basename(a.reference_class))
    write_atomic(a.summary, summary_text(rows, ref, a.input, rexdb))
    write_atomic(a.output, fasta_text(rows))
    n_internal = sum(r[2] in ref.internal for r in rows)
    print(f"OK: {len(rows)} proteins will be added to {rexdb} "
          f"({n_internal} of them with an internal class)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

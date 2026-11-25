#!/usr/bin/env python3
"""prepare_dante_extra_database.py: what it accepts, rejects and writes.

DANTE accepts nearly every malformed record with exit status 0 and changes the
annotation silently (measured with 0.2.12; see docs/dante_extra_database_plan.md),
so each of those cases must be an error here. Uses a small made-up REXdb; no
DANTE needed.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "prepare_dante_extra_database.py"

ALE = "Class_I|LTR|Ty1/copia|Ale"
CRM = "Class_I|LTR|Ty3/gypsy|chromovirus|CRM"
CLASS_TABLE = (f"REXdb_ID1\t{ALE}\nREXdb_ID2\t{CRM}\nREXdb_ID3\tClass_I|LINE\n"
               "REXdb_ID4\tClass_II|Subclass_1|TIR|hAT\n").replace("|", "\t")
PDB = (">Ale-RT__REXdb_ID1\nMKV\n>Ale-GAG__REXdb_ID1\nMKV\n"
       ">CRM-CHDCR__REXdb_ID2\nMKV\n>CRM-RT__REXdb_ID2\nMKV\n"
       ">LINE-ENDO__REXdb_ID3\nMKV\n>hAT-TPase__REXdb_ID4\nMKV\n")
SEQ = "MKVLDWNTYRHIKCEFGQSW" * 4
failures = []


def check(label, ok, detail=""):
    if not ok:
        failures.append(f"{label}: {detail}")


def run(data, seed="0"):
    """(return code, stderr, fasta or None, summary or None) for one input."""
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        (d / "X_v1_class").write_text(CLASS_TABLE)
        (d / "X_v1_pdb").write_text(PDB)
        (d / "extra.fa").write_bytes(data if isinstance(data, bytes) else data.encode())
        out, summary = d / "out.fa", d / "out.tsv"
        p = subprocess.run(
            [sys.executable, str(SCRIPT), str(d / "extra.fa"),
             "--reference-class", str(d / "X_v1_class"),
             "--reference-pdb", str(d / "X_v1_pdb"),
             "-o", str(out), "--summary", str(summary)],
            capture_output=True, text=True, env={**os.environ, "PYTHONHASHSEED": seed})
        check("no temporary file left", not list(d.glob("*.tmp")))
        return (p.returncode, p.stderr,
                out.read_bytes().decode() if out.exists() else None,
                summary.read_text() if summary.exists() else None)


def rec(header, seq=SEQ):
    return f">{header}\n{seq}\n"


# ── accepted ──────────────────────────────────────────────────────────────
good = rec(f"chr1:10-400 RT {ALE}") + rec(f"b RT {CRM}") + rec("c.2 TPase Class_II|Subclass_1|TIR|hAT")
rc, err, fasta, summary = run(good)
check("good file", rc == 0 and "WARNING" not in err, err)
check("good headers kept in order",
      fasta and [l for l in fasta.split("\n") if l.startswith(">")] ==
      [l for l in good.split("\n") if l.startswith(">")], fasta)
check("summary rows", summary and f"{ALE}\tRT\t1\tterminal\n" in summary
      and "#rexdb\tX_v1\n" in summary and "#proteins\t3\n" in summary
      and "#source\textra.fa\n" in summary, summary)

# CRLF, BOM, blank lines, repeated spaces, no final newline -> same clean FASTA
messy = b"\xef\xbb\xbf" + (good.replace("\n>", "\n\n>").replace(" RT ", "  RT\t")
                           .rstrip("\n").replace("\n", "\r\n").encode())
rc, err, messy_fasta, _ = run(messy)
check("messy file cleaned", rc == 0 and messy_fasta == fasta, err)
check("60 letters per line", all(len(l) <= 60 for l in fasta.split("\n")
                                 if not l.startswith(">")))
# the cleaned copy is itself valid and unchanged
check("idempotent", run(fasta)[2] == fasta)
# lower case and stop codons are normal in DANTE protein output
rc, err, _, _ = run(rec(f"a RT {ALE}", "mkvLDW*NTYR"))
check("lower case and stop codon", rc == 0, err)

# internal classes: accepted, warned about, marked in the summary
for internal in ("Class_I|LTR|Ty1/copia", "Class_I|LTR|Ty3/gypsy|chromovirus", "Class_I|LTR"):
    rc, err, fa, tsv = run(rec(f"a RT {internal}"))
    check(f"internal {internal}", rc == 0 and fa and "WARNING" in err
          and "internal class" in err and f"{internal}\tRT\t1\tinternal" in tsv, err)

# a domain REXdb does not have in that class is suspicious but not fatal
rc, err, _, _ = run(rec(f"a CHDCR {ALE}"))
check("unusual domain warns", rc == 0 and "no CHDCR domain" in err, err)

# ── rejected (exit 2, nothing written, the reason is named) ────────────────
BAD = {
    "new lineage": (rec(f"a RT {ALE}|Cicho"), "new lineages are not supported"),
    "underscore spelling": (rec("a RT Class_I|LTR|Ty1_copia|Ale"), f"spells it '{ALE}'"),
    "slash separators": (rec("a RT Class_I/LTR/Ty1_copia/Ale"), f"spells it '{ALE}'"),
    "typo in class": (rec("a RT Class_I|LTR|Ty1/copia|Alle"), "unknown class"),
    "unknown domain": (rec(f"a FOO {ALE}"), "unknown domain 'FOO'"),
    "domain suffix": (rec(f"a RT-x {ALE}"), "unknown domain 'RT-x'"),
    "two fields": (rec(f"a {ALE}"), "expected 3 fields"),
    "four fields": (rec(f"a RT {ALE} note"), "expected 3 fields"),
    "duplicate ID": (rec(f"a RT {ALE}") + rec(f"a RT {CRM}"), "same ID as line 1"),
    "duplicate after DANTE's - to _": (rec(f"a-1 RT {ALE}") + rec(f"a_1 RT {ALE}"), "same ID"),
    "ID used by REXdb": (rec(f"REXdb_ID1 RT {CRM}"), "already used in REXdb"),
    "REXdb ID after - to _": (rec(f"REXdb-ID1 RT {CRM}"), "already used in REXdb"),
    "semicolon in ID": (rec(f"a;x=1 RT {ALE}"), "ID may contain only"),
    "pipe in ID": (rec(f"a|b RT {ALE}"), "ID may contain only"),
    "DNA": (rec(f"a RT {ALE}", "ACGTN" * 20), "looks like DNA"),
    "gap": (rec(f"a RT {ALE}", "MKV-LDW"), "contains '-'"),
    "frameshift": (rec(f"a RT {ALE}", "MKV/LDW\\K"), "sequence contains"),
    "digit": (rec(f"a RT {ALE}", "MKV1LDW"), "contains '1'"),
    "empty sequence": (f">a RT {ALE}\n>b RT {ALE}\nMKV\n", "empty sequence"),
    "empty file": ("", "no sequences"),
}
for label, (data, reason) in BAD.items():
    rc, err, fa, tsv = run(data)
    check(label, rc == 2 and reason in err, f"rc={rc} {err}")
    check(f"{label}: nothing written", fa is None and tsv is None)

# every bad record is reported in one go, not only the first
rc, err, _, _ = run(rec(f"a FOO {ALE}") + rec(f"b;c RT {ALE}") + rec("d RT Nope"))
check("all errors reported", err.count("ERROR:") == 3, err)

# short sequences made of A/C/G/T letters are legitimate peptides
check("short ACGT peptide", run(rec(f"a RT {ALE}", "ACGTACGT"))[0] == 0)

# ── unreadable (exit 3) ────────────────────────────────────────────────────
check("not FASTA", run("hello\n>a RT x\nMKV\n")[0] == 3)
check("binary", run(b"\x1f\x8b\xff\xfe\x00")[0] == 3)

# ── determinism ────────────────────────────────────────────────────────────
many = "".join(rec(f"id{i} {'RT' if i % 2 else 'GAG'} {ALE}") for i in range(40))
check("hash seed invariance", run(many, "1") == run(many, "2"))

if failures:
    print("FAIL: prepare_dante_extra_database")
    for f in failures:
        print(f"  - {f}")
    sys.exit(1)
print(f"OK: prepare_dante_extra_database ({len(BAD)} bad inputs rejected)")

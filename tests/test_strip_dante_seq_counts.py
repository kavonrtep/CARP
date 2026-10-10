#!/usr/bin/env python3
"""strip_dante_seq_counts.py (issue #5): drops exactly DANTE's per-sequence
SEQ/DOMAIN/COUNTS header table and nothing else.

That table is ~1 row per sequence x domain (1.1 M rows on a 305k-sequence
assembly) and rtracklayer's import.gff3 is quadratic in comment lines, so the
filter_dante rule cuts it out. The cut must never touch a feature line or the
other header lines (thresholds, classification counts)."""
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))
from strip_dante_seq_counts import strip_lines  # noqa: E402

SEP = "##-----------------------------------------------\n"
HEAD = [
    "##gff-version 3\n",
    "##DANTE version 0.2.12\n",
    "##Filtering thresholdss: min identity: 0.35 ...\n",
    "##CLASSIFICATION\tORIGINAL_COUNTS\tFILTERED_COUNTS\n",
    "##Class_I|LTR|Ty3/gypsy\t12\t10\n",
    SEP,
]
TABLE = ["##SEQ\tDOMAIN\tCOUNTS\n",
         "##CM116228.1\tCHD\t2075\n",
         "##JBDLLP010305167.1\taRH\t1\n",
         SEP]
FEATS = ["chr1\tdante\tprotein_domain\t5850\t6120\t180\t+\t.\tName=GAG\n",
         "chr2\tdante\tprotein_domain\t100\t400\t90\t-\t.\tName=RT\n"]

failures = []


def check(label, got, want):
    if got != want:
        failures.append(label)
        print(f"FAIL {label}\n  got:  {got!r}\n  want: {want!r}")


# 1. The DANTE layout: table removed, everything else byte-identical in order.
check("dante layout", list(strip_lines(HEAD + TABLE + FEATS)), HEAD + FEATS)

# 2. Empty-table form DANTE writes when there are no domains.
check("no-domains form",
      list(strip_lines(HEAD + ["##SEQ\tDOMAIN\tCOUNTS\n", "##NO DOMAINS\n", SEP])),
      HEAD)

# 3. File without the table (e.g. RepeatMasker GFF3): untouched.
plain = ["##gff-version 3\n", "##sequence-region chr1 1 1000\n", SEP] + FEATS
check("no table -> no-op", list(strip_lines(plain)), plain)

# 4. Unterminated table: features after it are never dropped.
check("unterminated table keeps features",
      list(strip_lines(HEAD + TABLE[:-1] + FEATS)), HEAD + FEATS)

# 5. CRLF line endings are recognised.
crlf = [l.replace("\n", "\r\n") for l in HEAD + TABLE + FEATS]
check("CRLF", list(strip_lines(crlf)),
      [l.replace("\n", "\r\n") for l in HEAD + FEATS])

# 6. CLI: several inputs concatenated in order.
with tempfile.TemporaryDirectory() as d:
    a, b, out = Path(d, "a.gff3"), Path(d, "b.gff3"), Path(d, "out.gff3")
    a.write_text("".join(plain))
    b.write_text("".join(HEAD + TABLE + FEATS))
    subprocess.run([sys.executable, str(SCRIPTS / "strip_dante_seq_counts.py"),
                    str(a), str(b), "-o", str(out)], check=True)
    check("CLI concat", out.read_text(), "".join(plain + HEAD + FEATS))
    check("CLI leaves no .tmp", sorted(p.name for p in Path(d).iterdir()),
          ["a.gff3", "b.gff3", "out.gff3"])
    r = subprocess.run([sys.executable, str(SCRIPTS / "strip_dante_seq_counts.py"),
                        str(Path(d, "missing.gff3")), "-o", str(out)],
                       capture_output=True)
    check("CLI missing input fails", r.returncode != 0, True)

if failures:
    sys.exit(f"test_strip_dante_seq_counts: {len(failures)} FAILED")
print("test_strip_dante_seq_counts: PASSED")

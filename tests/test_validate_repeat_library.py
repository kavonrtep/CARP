#!/usr/bin/env python3
"""Tests for scripts/validate_repeat_library.py (user library check + conversion).

Covers the two library types the pipeline accepts:
  * tandem_repeat_library: bare names -> Satellite/<name>, rDNA names -> rDNA
    classes (45S before 5S), TE classes refused (they belong in custom_library);
  * custom_library: shared canonicaliser (DANTE pipes, DANTE_TIR underscores),
    bare CARP lineage names, the RepeatMasker/Dfam dictionary
    (library_class_aliases in classification_vocabulary.yaml), unknown classes
    refused, never guessed.
Plus FASTA hygiene (descriptions, gaps, U->T, protein input, duplicate IDs),
exit codes, --drop-invalid, and that every written library passes the
pipeline's own validate_classifications gate. Data-free; runs in carp-unit.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import classification  # noqa: E402
import validate_repeat_library as vrl  # noqa: E402

FAILS = 0


def eq(got, want, label):
    global FAILS
    if got != want:
        FAILS += 1
        print(f"FAIL {label}: got {got!r}, want {want!r}")


def run(fasta_text, kind, *extra):
    d = Path(tempfile.mkdtemp())
    src = d / "in.fa"
    src.write_text(fasta_text)
    rc = vrl.main([str(src), "-t", kind, "-o", str(d / "out"), *extra])
    rows = (d / "out.tsv").read_text().splitlines()
    hdr = rows[0].split("\t")
    recs = [dict(zip(hdr, r.split("\t"))) for r in rows[1:]]
    fa = d / "out.fasta"
    return rc, recs, (fa.read_text() if fa.exists() else None), d


def classes(recs):
    return [r["output_class"] for r in recs]


# ---- 1. tandem: the example-library shapes --------------------------------
rc, recs, fa, d = run(
    ">X04322.1#CEN180\nACGT\n>M65137.1#5S_rDNA\nACGT\n>X52322.1#45S_rDNA\nACGT\n"
    ">t#telomere\nTTTAGGG\n>s#Satellite/PisTR-B\nACGT\n>r#18S\nACGT\n"
    ">noclass\nACGT\n>sat#satellite\nACGT\n>q#AthSat500?\nACGT\n", "tandem")
eq(rc, 0, "tandem: converted -> library written, exit 0")
eq(classes(recs), ["Satellite/CEN180", "rDNA/5S_rDNA", "rDNA/45S_rDNA",
                   "Satellite/telomere", "Satellite/PisTR-B", "rDNA/45S_rDNA/18S",
                   "Satellite/noclass", "Satellite", "Satellite/AthSat500"],
   "tandem: class conversion")
eq(recs[4]["status"], "unchanged", "tandem: canonical Satellite left alone")
eq(classification.validate_values(classification.collect_fasta_classes(d / "out.fasta"),
                                  "RepeatMasker"), [], "tandem: output passes the pipeline gate")

# ---- 2. rDNA token boundaries ---------------------------------------------
eq(vrl.rdna_class("45S_rDNA"), "rDNA/45S_rDNA", "rdna: 45S not read as 5S")
eq(vrl.rdna_class("5S"), "rDNA/5S_rDNA", "rdna: 5S")
eq(vrl.rdna_class("rDNA_5.8S"), "rDNA/45S_rDNA/5.8S", "rdna: 5.8S not read as 5S")
eq(vrl.rdna_class("CEN180"), None, "rdna: satellite name untouched")
eq(vrl.rdna_class("AthSat5Sx"), None, "rdna: 5S inside a word ignored")

# ---- 3. tandem: TE classes refused ----------------------------------------
rc, recs, fa, _ = run(">a#Class_I/LTR/Ty1_copia\nACGT\n>b#LTR/Gypsy\nACGT\n"
                      ">c#Satellite/A\nACGT\n", "tandem")
eq(rc, 2, "tandem: TE class -> needs manual fix")
eq(fa, None, "tandem: no library written on error")
eq([r["status"] for r in recs], ["error", "error", "unchanged"], "tandem: TE rows are errors")
rc, recs, fa, _ = run(">a#LTR/Gypsy\nACGT\n>c#Satellite/A\nACGT\n", "tandem", "--drop-invalid")
eq(rc, 0, "tandem --drop-invalid: filtered library written, exit 0")
eq(fa, ">c#Satellite/A\nACGT\n", "tandem --drop-invalid keeps only valid")

# ---- 4. custom: translation paths ------------------------------------------
rc, recs, fa, _ = run(
    ">a#LTR/Gypsy\nACGT\n>b#DNA/hAT-Ac?\nACGT\n>c#LINE/L1\nACGT\n>d#Ty1/copia/Ale\nACGT\n"
    ">e#Athila\nACGT\n>f#RC/Helitron\nACGT\n>g#Class_I|LTR|Ty3/gypsy|chromovirus|Tekay\nACGT\n"
    ">h#Class_II_Subclass_1_TIR_hAT\nACGT\n>j#satellite/cen\nACGT\n>k#rRNA\nACGT\n"
    ">l#DNA/CMC-EnSpm\nACGT\n>m#Class_I/LTR/Ty1_copia/SIRE\nACGT\n", "custom")
eq(rc, 0, "custom: converted exit code")
eq(classes(recs), [
    "Class_I/LTR/Ty3_gypsy", "Class_II/Subclass_1/TIR/hAT", "Class_I/LINE",
    "Class_I/LTR/Ty1_copia/Ale", "Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Athila",
    "Class_II/Subclass_2/Helitron", "Class_I/LTR/Ty3_gypsy/chromovirus/Tekay",
    "Class_II/Subclass_1/TIR/hAT", "Satellite/cen", "rDNA",
    "Class_II/Subclass_1/TIR/EnSpm_CACTA", "Class_I/LTR/Ty1_copia/SIRE"],
   "custom: class conversion")

# ---- 4b. Wicker et al. 2007 three-letter codes (superfamily level only) ----
R0 = vrl.Resolver(None)
for code, want in [("RLC", "Class_I/LTR/Ty1_copia"), ("RLG_Athila", "Class_I/LTR/Ty3_gypsy"),
                   ("RLX", "Class_I/LTR"), ("RIL", "Class_I/LINE"), ("RST", "Class_I/SINE"),
                   ("RXX", "Class_I"), ("DTC-CACTA_1", "Class_II/Subclass_1/TIR/EnSpm_CACTA"),
                   ("DTA", "Class_II/Subclass_1/TIR/hAT"), ("DTX", "Class_II/Subclass_1/TIR"),
                   ("DHH", "Class_II/Subclass_2/Helitron"), ("DXX", "Class_II")]:
    eq(R0.by_dictionary(code)[0], want, f"wicker: {code}")
eq(R0.by_dictionary("DTA1"), None, "wicker: code glued to a name is not a code")
rc, recs, fa, _ = run(">a#RLC\nACGT\n", "tandem")
eq(recs[0]["status"], "error", "wicker: TE code refused in a tandem library")

# ---- 5. custom: refusals ----------------------------------------------------
rc, recs, fa, _ = run(">i#tRNA\nACGT\n>m\nACGT\n>n#LTR/Copia\nMKLVPQ\n>e#CEN180\nACGT\n", "custom")
eq(rc, 2, "custom: unknown/missing/protein -> manual fix")
eq([r["status"] for r in recs], ["error"] * 4, "custom: all four refused")

# ---- 6. FASTA hygiene + duplicate IDs ---------------------------------------
rc, recs, fa, _ = run(">a#Satellite/x desc words\nAC-G.Tu\n>a#Satellite/y\nacgt\n", "tandem")
eq(rc, 0, "hygiene: converted")
eq(fa, ">a#Satellite/x\nACGTt\n>a_dup2#Satellite/y\nacgt\n", "hygiene: gaps, U, case, dedupe")

# ---- 7. already valid -> exit 0, input reproduced ---------------------------
rc, recs, fa, _ = run(">a#Satellite/x\nACGT\n>b#rDNA/5S_rDNA\nACGT\n", "tandem")
eq(rc, 0, "valid: exit 0")
eq(fa, ">a#Satellite/x\nACGT\n>b#rDNA/5S_rDNA\nACGT\n", "valid: output identical")

# ---- 7b. explicit output paths (Galaxy): all three outputs always exist ----
d = Path(tempfile.mkdtemp())
(d / "in.fa").write_text(">a#LTR/Gypsy\nACGT\n>b#CEN\nACGT\n")
args = [str(d / "in.fa"), "-t", "tandem", "--out-fasta", str(d / "lib.dat"),
        "--out-tsv", str(d / "rep.dat"), "--out-html", str(d / "rep_html.dat")]
eq(vrl.main(args), 2, "galaxy: errors -> exit 2")
eq([(d / f).exists() for f in ("lib.dat", "rep.dat", "rep_html.dat")], [True] * 3,
   "galaxy: every declared output created")
eq((d / "lib.dat").read_text(), "", "galaxy: no library -> empty FASTA, not a partial one")
eq(vrl.main(args + ["--drop-invalid"]), 0, "galaxy: --drop-invalid -> exit 0")
eq((d / "lib.dat").read_text(), ">b#Satellite/CEN\nACGT\n", "galaxy: filtered library")

# ---- 7c. golden fixtures: byte-identical FASTA + TSV -------------------------
FX = ROOT / "tests" / "fixtures" / "validate_repeat_library"
for kind, src, extra, want_rc, want_fa, want_tsv in [
        ("tandem", "tandem_input.fasta", [], 0, "tandem_expected.fasta", "tandem_expected.tsv"),
        ("custom", "custom_input.fasta", [], 2, "custom_expected_empty.fasta", "custom_expected.tsv"),
        ("custom", "custom_input.fasta", ["--drop-invalid"], 0, "custom_dropped_expected.fasta", None)]:
    d = Path(tempfile.mkdtemp())
    got_rc = vrl.main([str(FX / src), "-t", kind, *extra, "--out-fasta", str(d / "o.fasta"),
                       "--out-tsv", str(d / "o.tsv"), "--out-html", str(d / "o.html")])
    label = f"golden {src} {' '.join(extra)}".strip()
    eq(got_rc, want_rc, f"{label}: exit code")
    eq((d / "o.fasta").read_text(), (FX / want_fa).read_text(), f"{label}: FASTA")
    if want_tsv:
        eq((d / "o.tsv").read_text(), (FX / want_tsv).read_text(), f"{label}: TSV")

# ---- 8. dictionary targets are canonical (also enforced at load) -----------
R = vrl.Resolver(None)
eq(all(R.vocab.is_canonical(t) for _, _, t in R.aliases), True, "dictionary targets canonical")

# ---- 9. HTML written and self-contained -------------------------------------
_, _, _, d = run(">a#CEN180\nACGT\n", "tandem")
page = (d / "out.html").read_text()
eq("<script src" in page or "<link" in page, False, "html: no external resources")
eq("Satellite/CEN180" in page, True, "html: shows converted header")

# ---- 10. same output under different PYTHONHASHSEED ------------------------
outs = []
for seed in ("0", "1", "12345"):
    d = Path(tempfile.mkdtemp())
    (d / "in.fa").write_text(">b#LTR/Gypsy\nACGT\n>a#CEN\nACGT\n>a#DNA/hAT\nACGT\n")
    subprocess.run([sys.executable, str(ROOT / "scripts" / "validate_repeat_library.py"),
                    str(d / "in.fa"), "-t", "custom", "-o", str(d / "o"), "--drop-invalid"],
                   env={"PYTHONHASHSEED": seed, "PATH": "/usr/bin:/bin"},
                   capture_output=True, check=False)
    outs.append((d / "o.fasta").read_text() + (d / "o.tsv").read_text())
eq(len(set(outs)), 1, "deterministic across PYTHONHASHSEED")

if FAILS:
    print(f"test_validate_repeat_library: {FAILS} FAILED")
    sys.exit(1)
print("test_validate_repeat_library: PASSED")

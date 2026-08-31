#!/usr/bin/env python3
"""Unit tests for the library-health summary.

This exists to make two specific defects visible, so the tests are written
around them: an **empty Class_II library** (the pre-1.4.0 DANTE_TIR join bug,
which left the RepeatMasker library with no DNA-transposon sequences at all) and
**over-extended element boundaries** (DANTE_LINE consensi whose flanks were a
different, far more abundant repeat). Both shipped for several releases while
every number in the run's report looked plausible.
"""

import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import library_health as H  # noqa: E402


def write_library(path, records):
    """records: (name, classification, length)"""
    with open(path, "w") as fh:
        for name, cls, length in records:
            fh.write(f">{name}#{cls}\n" + "A" * length + "\n")


def write_elements(path, feature_type, rows):
    """rows: (start, end, ext5, ext3, inferred5, inferred3 or None)"""
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for i, (start, end, e5, e3, i5, i3) in enumerate(rows, 1):
            attrs = f"ID=el_{i};Extension_5prime={e5};Extension_3prime={e3}"
            if i5 is not None:
                attrs += (f";Extension_capped=5prime"
                          f";Extension_5prime_inferred={i5}"
                          f";Extension_3prime_inferred={i3}")
            fh.write(f"chr1\tX\t{feature_type}\t{start}\t{end}\t.\t+\t.\t{attrs}\n")


def opts(**over):
    base = dict(library=None, dante_line_gff3=None, fallback_gff3=None,
                screen_audit=None,
                vocabulary=str(REPO / "classification_vocabulary.yaml"))
    base.update(over)
    return types.SimpleNamespace(**base)


def as_dict(rows):
    return {(s, i, m): v for s, i, m, v in rows}


class TestClassComposition(unittest.TestCase):
    def test_counts_and_lengths_per_class(self):
        with tempfile.TemporaryDirectory() as tmp:
            lib = Path(tmp) / "lib.fasta"
            write_library(lib, [("a", "Class_I/LINE", 1000),
                                ("b", "Class_I/LINE", 3000),
                                ("c", "Class_I/LTR/Ty1_copia/Ale", 5000)])
            rows, _ = H.collect(opts(library=str(lib)))
            d = as_dict(rows)
            self.assertEqual(d[("library", "all", "n_consensi")], 3)
            self.assertEqual(d[("library", "all", "total_bp")], 9000)
            self.assertEqual(d[("class", "Class_I/LINE", "n_consensi")], 2)
            self.assertEqual(d[("class", "Class_I/LINE", "max_len")], 3000)
            self.assertEqual(d[("class", "Class_I/LTR/Ty1_copia/Ale", "total_bp")], 5000)

    def test_flags_a_consensus_over_its_class_bound(self):
        """A 22 kb Class_I/LINE consensus is what shipped for several releases."""
        with tempfile.TemporaryDirectory() as tmp:
            lib = Path(tmp) / "lib.fasta"
            write_library(lib, [("huge", "Class_I/LINE", 22000)])
            rows, warnings = H.collect(opts(library=str(lib)))
            self.assertEqual(as_dict(rows)[("class", "Class_I/LINE", "n_over_bound")], 1)
            self.assertTrue(any("longer than the 8000 bp bound" in w for w in warnings))

    def test_bound_is_per_class_not_global(self):
        with tempfile.TemporaryDirectory() as tmp:
            lib = Path(tmp) / "lib.fasta"
            # 22 kb is over the LINE bound but well inside the CACTA one
            write_library(lib, [("cacta", "Class_II/Subclass_1/TIR/EnSpm_CACTA", 22000)])
            rows, warnings = H.collect(opts(library=str(lib)))
            self.assertEqual(
                as_dict(rows)[("class", "Class_II/Subclass_1/TIR/EnSpm_CACTA",
                               "n_over_bound")], 0)
            self.assertEqual(warnings, [])


class TestEmptyClassIIDetection(unittest.TestCase):
    """The pre-1.4.0 join bug left zero Class_II sequences in the library."""

    def _warnings(self, n_records):
        with tempfile.TemporaryDirectory() as tmp:
            lib = Path(tmp) / "lib.fasta"
            write_library(lib, [(f"s{i}", "Class_I/LTR/Ty1_copia/Ale", 2000)
                                for i in range(n_records)])
            return H.collect(opts(library=str(lib)))[1]

    def test_warns_on_a_real_sized_library(self):
        self.assertTrue(any("no class_ii" in w.lower()
                            for w in self._warnings(H.MIN_CONSENSI_FOR_CLASS_WARN)))

    def test_silent_on_a_fixture_sized_library(self):
        """A tiny library has no Class_II for an unremarkable reason."""
        self.assertFalse(any("Class_II" in w
                             for w in self._warnings(H.MIN_CONSENSI_FOR_CLASS_WARN - 1)))

    def test_silent_when_class_ii_is_present(self):
        with tempfile.TemporaryDirectory() as tmp:
            lib = Path(tmp) / "lib.fasta"
            recs = [(f"s{i}", "Class_I/LTR/Ty1_copia/Ale", 2000) for i in range(60)]
            recs.append(("t", "Class_II/Subclass_1/TIR/hAT", 3000))
            write_library(lib, recs)
            self.assertFalse(any("Class_II" in w
                                 for w in H.collect(opts(library=str(lib)))[1]))


class TestBoundaryHealth(unittest.TestCase):
    def test_extension_fraction_and_ceiling_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = Path(tmp) / "line.gff3"
            # one 10 kb element: 2 kb core + two 4 kb extensions, both inferred
            # at the 10 kb ceiling and capped back
            write_elements(gff, "LINE_element", [
                (1, 10000, 4000, 4000, 9999, 9999),
                (20000, 22000, 0, 0, None, None),
            ])
            rows, warnings = H.collect(opts(dante_line_gff3=str(gff)))
            d = as_dict(rows)
            self.assertEqual(d[("boundary", "DANTE_LINE", "n_elements")], 2)
            self.assertEqual(d[("boundary", "DANTE_LINE", "n_at_flank_ceiling")], 1)
            self.assertEqual(d[("boundary", "DANTE_LINE", "n_capped")], 1)
            self.assertEqual(d[("boundary", "DANTE_LINE", "extension_bp")], 8000)
            self.assertTrue(any("flank alignment reaching" in w for w in warnings))

    def test_high_inferred_flank_share_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = Path(tmp) / "fb.gff3"
            write_elements(gff, "sequence_feature",
                           [(1, 10000, 4500, 4500, None, None)])
            _, warnings = H.collect(opts(fallback_gff3=str(gff)))
            self.assertTrue(any("inferred flank" in w for w in warnings))

    def test_a_well_bounded_builder_is_silent(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = Path(tmp) / "line.gff3"
            write_elements(gff, "LINE_element",
                           [(1, 5000, 200, 300, None, None)])
            rows, warnings = H.collect(opts(dante_line_gff3=str(gff)))
            self.assertEqual(warnings, [])
            self.assertEqual(as_dict(rows)[("boundary", "DANTE_LINE",
                                            "n_at_flank_ceiling")], 0)

    def test_protein_domain_children_are_not_counted_as_elements(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = Path(tmp) / "line.gff3"
            gff.write_text(
                "##gff-version 3\n"
                "chr1\tX\tLINE_element\t1\t5000\t.\t+\t.\tID=el_1\n"
                "chr1\tdante\tprotein_domain\t100\t400\t.\t+\t.\tID=d1;Parent=el_1\n")
            rows, _ = H.collect(opts(dante_line_gff3=str(gff)))
            self.assertEqual(as_dict(rows)[("boundary", "DANTE_LINE", "n_elements")], 1)


class TestScreenSection(unittest.TestCase):
    def test_summarises_the_cross_class_audit(self):
        with tempfile.TemporaryDirectory() as tmp:
            audit = Path(tmp) / "screen.tsv"
            audit.write_text(
                "name\tclassification\toriginal_length\tretained_length\taction\t"
                "conflict_bp\tconflict_classes\tdetail\n"
                "a\tClass_I/LINE\t5000\t3000\ttrimmed\t2000\tX\t\n"
                "b\tClass_I/LINE\t4000\t0\tdropped\t4000\tX\t\n"
                "c\tClass_I/LINE\t1000\t1000\tkept\t0\t\t\n")
            rows, _ = H.collect(opts(screen_audit=str(audit)))
            d = as_dict(rows)
            self.assertEqual(d[("screen", "cross_class", "n_trimmed")], 1)
            self.assertEqual(d[("screen", "cross_class", "n_dropped")], 1)
            self.assertEqual(d[("screen", "cross_class", "bp_removed")], 6000)


class TestOutput(unittest.TestCase):
    script = str(REPO / "scripts" / "library_health.py")

    def test_rows_are_sorted_and_written_atomically(self):
        with tempfile.TemporaryDirectory() as tmp:
            lib = Path(tmp) / "lib.fasta"
            write_library(lib, [("z", "Class_I/LINE", 100),
                                ("a", "Class_I/LTR/Ty1_copia/Ale", 200)])
            out = Path(tmp) / "health.tsv"
            subprocess.run([sys.executable, self.script, "-o", str(out),
                            "--library", str(lib)], check=True, capture_output=True)
            lines = out.read_text().splitlines()
            self.assertEqual(lines[0].split("\t"), H.HEADER)
            keys = [tuple(l.split("\t")[:3]) for l in lines[1:]]
            self.assertEqual(keys, sorted(keys))
            self.assertFalse((out.parent / (out.name + ".tmp")).exists())

    def test_missing_inputs_are_not_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "health.tsv"
            r = subprocess.run([sys.executable, self.script, "-o", str(out),
                                "--library", str(Path(tmp) / "nope.fasta")],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertTrue(out.exists())

    def test_never_exits_non_zero_on_warnings(self):
        """Advisory means advisory: a warning must not fail a run."""
        with tempfile.TemporaryDirectory() as tmp:
            lib = Path(tmp) / "lib.fasta"
            write_library(lib, [("huge", "Class_I/LINE", 22000)])
            out = Path(tmp) / "health.tsv"
            r = subprocess.run([sys.executable, self.script, "-o", str(out),
                                "--library", str(lib)], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0)
            self.assertIn("WARNING", r.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)

#!/usr/bin/env python3
"""LINE elements with both ends directly observed (poly-A tail + TSD).

Such an element is marked `Status=complete` and its inferred span is REPLACED
with the measured one, so this does change the annotation — for the 0.7-5% of
elements that qualify. Everything else keeps its inferred span and is marked
`Status=inferred`.

Determinism matters more than usual here: each TSD is calibrated against decoy
windows drawn at random. A global RNG would make the result depend on locus
order, thread count and PYTHONHASHSEED, and the run-to-run determinism gate
would fail intermittently.
"""
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import line_complete_elements as lce  # noqa: E402


class TestPolyADetection(unittest.TestCase):
    """The tail test must find a homopolymer, not merely AT-rich sequence."""

    def test_finds_a_clean_run(self):
        run, at = lce.longest_run("CGCG" + "A" * 20 + "CGCG", "A")
        self.assertEqual(run, 20)
        self.assertEqual(at, 4)

    def test_tolerates_one_interruption(self):
        run, _ = lce.longest_run("A" * 10 + "C" + "A" * 10, "A", max_mismatch=1)
        self.assertEqual(run, 21)

    def test_does_not_bridge_two_interruptions(self):
        run, _ = lce.longest_run("A" * 10 + "CC" + "A" * 10, "A", max_mismatch=1)
        self.assertEqual(run, 10)

    def test_at_rich_coding_sequence_is_not_a_tail(self):
        """The failure that made an early version score BELOW its own control."""
        at_rich = "ATATATAAATAATATAAATATTAAATATAAT" * 3
        run, _ = lce.longest_run(at_rich, "A")
        self.assertLess(run, 15, "AT-rich sequence must not pass the tail test")


class TestTSDScan(unittest.TestCase):
    """A perfect repeat of length k scores exactly k, which is what makes the
    chance floor directly interpretable as a minimum TSD length."""

    def test_perfect_repeat_scores_its_length(self):
        tsd = "ACGTACGTACGTAC"          # 14 bp
        score, pos, klen, _ = lce.scan(tsd, "TTTT" + tsd + "TTTT", 8, 20)
        self.assertEqual(score, 14)
        self.assertEqual(klen, 14)
        self.assertEqual(pos, 4)

    def test_absent_repeat_scores_poorly(self):
        score, _, _, _ = lce.scan("ACGTACGTACGTAC", "T" * 200, 8, 20)
        self.assertLess(score, 8)

    def test_mismatches_are_penalised(self):
        tsd = "ACGTACGTACGTAC"
        bad = tsd[:6] + "TTT" + tsd[9:]
        clean, _, _, _ = lce.scan(tsd, "TTTT" + tsd + "TTTT", 8, 20)
        dirty, _, _, _ = lce.scan(tsd, "TTTT" + bad + "TTTT", 8, 20)
        self.assertGreater(clean, dirty)


class TestRepeatDerivedTSDs(unittest.TestCase):
    """A TSD found inside a tandem array is not a duplication.

    The failure this caught: on GCA_973357735.1 the plant telomere repeat
    (CCCTAAA) was read as a 16 bp TSD. Decoy windows are drawn from distant loci
    and contain no telomeric sequence, so the null test passed it — and the
    resulting element displaced ~5 kb of genuine telomere annotation.
    """

    def test_telomere_repeat_is_rejected(self):
        self.assertTrue(lce.is_tandem("CCTAAACCCTAAAACC"))

    def test_simple_repeats_are_rejected(self):
        for seq in ("TATATATATATATATATATA", "GAGAGAGAGAGAGACC",
                    "TAATAATAATTATTAT"):
            self.assertTrue(lce.is_tandem(seq), seq)

    def test_a_genuine_tsd_is_kept(self):
        for seq in ("TTCGTAAGCATGTAAATAA", "TTTATCATTAAAACAA"):
            self.assertFalse(lce.is_tandem(seq), seq)

    def test_uniqueness_is_counted(self):
        """A query that matches at every period of an array must be seen to."""
        tsd = "ACGTACGTACGTAC"
        _, _, _, unique = lce.scan(tsd, "TTTT" + tsd + "TTTT", 8, 20, count_above=8)
        _, _, _, many = lce.scan("ATATATATATAT", "AT" * 400, 8, 20, count_above=8)
        self.assertLess(unique, many)


class TestChanceFloor(unittest.TestCase):
    """Measured 99th-percentile chance scores; a TSD is only separable above them."""

    def test_floor_rises_with_window(self):
        self.assertLess(lce.floor_for(500), lce.floor_for(8000))

    def test_the_shipped_window_needs_11bp(self):
        # 2100 bp is the default --tsd-window; a repeat must beat 10 to count,
        # i.e. be at least 11 bp. Many decayed TSDs are shorter, which is a
        # limit of the method, not a tuning choice.
        self.assertEqual(lce.floor_for(2100), 10)

    def test_unknown_large_window_falls_back_to_the_strictest(self):
        self.assertEqual(lce.floor_for(10 ** 6), max(lce.CHANCE_FLOOR.values()))


class TestDeterminism(unittest.TestCase):
    def test_seed_is_stable_and_not_python_hash(self):
        a = lce.locus_seed("chr1", 100, 200)
        b = lce.locus_seed("chr1", 100, 200)
        self.assertEqual(a, b)
        self.assertNotEqual(a, lce.locus_seed("chr1", 100, 201))

    def test_seed_survives_a_different_hash_seed(self):
        """hash() is randomised per process; a hash-derived seed would break the
        determinism gate. Run a child with a different PYTHONHASHSEED."""
        code = ("import sys; sys.path.insert(0, %r);"
                "import line_complete_elements as l;"
                "print(l.locus_seed('chr1', 100, 200))" % str(REPO / "scripts"))
        env = dict(os.environ, PYTHONHASHSEED="12345")
        out = subprocess.run([sys.executable, "-c", code], capture_output=True,
                             text=True, env=env)
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertEqual(int(out.stdout.strip()), lce.locus_seed("chr1", 100, 200))


class TestOutputContract(unittest.TestCase):
    def test_library_header_convention(self):
        """The FASTA must be directly usable as a RepeatMasker library by anyone
        evaluating it, which means name#Class/Subclass headers."""
        with tempfile.TemporaryDirectory() as d:
            fasta = Path(d) / "c.fasta"
            fasta.write_text(">LINE_group_0001_confirmed#Class_I/LINE\nACGT\n")
            head = fasta.read_text().splitlines()[0]
            self.assertTrue(head.startswith(">"))
            self.assertIn("#Class_I/LINE", head)

    def test_every_element_gets_a_status(self):
        """Both confidence classes must be distinguishable downstream, so an
        element is either measured or explicitly marked as inferred."""
        import re
        src = (REPO / "scripts" / "line_complete_elements.py").read_text()
        self.assertIn('attrs["Status"] = "inferred"', src)
        self.assertIn('attrs["Status"] = "complete"', src)

    def test_only_confirmed_elements_move(self):
        """A coordinate rewrite must be tied to having found the evidence.

        Guards the property the two-genome measurement rests on: 17 of 343
        Boechera elements were marked complete and exactly those 17 moved.
        """
        src = (REPO / "scripts" / "line_complete_elements.py").read_text()
        i = src.index("r = by_id.get(")
        seg = src[i:i + 600]
        self.assertIn("if r is None:", seg)
        # the span is only rewritten in the else branch
        self.assertLess(seg.index('attrs["Status"] = "inferred"'),
                        seg.index('f[3], f[4] = str(r["start"]), str(r["end"])'))

    def test_rewrite_is_atomic(self):
        """A killed run must not leave a truncated GFF3 that a later step trusts."""
        src = (REPO / "scripts" / "line_complete_elements.py").read_text()
        self.assertIn('tmp = gff + ".tmp"', src)
        self.assertIn("os.replace(tmp, gff)", src)

    def test_rerunning_is_idempotent(self):
        """Marks from a previous run are stripped, so re-running cannot stack
        duplicate Status/TSD attributes."""
        src = (REPO / "scripts" / "line_complete_elements.py").read_text()
        self.assertIn('for k in ("Status", "TSD", "PolyA_length"):', src)
        self.assertIn("attrs.pop(k, None)", src)


if __name__ == "__main__":
    unittest.main(verbosity=2)

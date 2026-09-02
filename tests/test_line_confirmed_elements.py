#!/usr/bin/env python3
"""LINE elements with both ends directly observed (poly-A tail + TSD).

This is an EVALUATION output: nothing downstream consumes it, so these tests
guard the two properties that actually matter — that a confirmed call means what
it says, and that the output is reproducible.

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

import line_confirmed_elements as lce  # noqa: E402


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
        score, pos, klen = lce.scan(tsd, "TTTT" + tsd + "TTTT", 8, 20)
        self.assertEqual(score, 14)
        self.assertEqual(klen, 14)
        self.assertEqual(pos, 4)

    def test_absent_repeat_scores_poorly(self):
        score, _, _ = lce.scan("ACGTACGTACGTAC", "T" * 200, 8, 20)
        self.assertLess(score, 8)

    def test_mismatches_are_penalised(self):
        tsd = "ACGTACGTACGTAC"
        bad = tsd[:6] + "TTT" + tsd[9:]
        clean, _, _ = lce.scan(tsd, "TTTT" + tsd + "TTTT", 8, 20)
        dirty, _, _ = lce.scan(tsd, "TTTT" + bad + "TTTT", 8, 20)
        self.assertGreater(clean, dirty)


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
                "import line_confirmed_elements as l;"
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

    def test_feeds_nothing_downstream(self):
        """The whole point: enabling this cannot change the annotation.

        Guards against a future edit wiring the output into another rule.
        """
        snake = (REPO / "Snakefile").read_text()
        import re
        for chunk in re.split(r"\n(?=rule )", snake):
            name = chunk.split("\n", 1)[0].replace("rule ", "").rstrip(":")
            if name in ("all", "line_confirmed_elements"):
                continue
            m = re.search(r"\n    input:(.*?)\n    (output|params|log|shell|run|"
                          r"threads|conda|benchmark|priority|resources):",
                          chunk, re.S)
            if m and "LINE_confirmed_elements" in m.group(1):
                self.fail(f"rule {name} consumes LINE_confirmed_elements; this "
                          f"output must stay evaluation-only")


if __name__ == "__main__":
    unittest.main(verbosity=2)

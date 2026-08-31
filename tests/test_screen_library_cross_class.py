#!/usr/bin/env python3
"""Unit tests for the cross-class library screen.

The screen exists because every other reduction CARP runs is within a class, so
a consensus that is part LINE and part LTR is invisible to all of them. The two
things most likely to go wrong here, and so most heavily tested:

  * **direction** -- a BLAST hit is symmetric, and a screen that cannot tell the
    chimera from its victim damages the correct library. On a real wheat library
    the naive symmetric rule trimmed 1,114 Ty1_copia/Angela consensi against 86
    LINE consensi. `TestOwnership` pins the rule that fixes it.
  * **sibling safety** -- Ty1_copia/Ale and Ty1_copia/Angela share genuine
    homology and must never be trimmed against each other.
"""

import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import screen_library_cross_class as S  # noqa: E402


def opts(**over):
    base = dict(min_identity=80.0, min_length=200, max_shared_depth=1,
                min_retained_length=300,
                min_retained_fraction=0.2, ownership_margin=0.10,
                vocabulary=str(REPO / "classification_vocabulary.yaml"))
    base.update(over)
    return types.SimpleNamespace(**base)


def hits_file(tmp, rows):
    """rows: (qid, sid, pident, alen, qlen, slen, qstart, qend)"""
    path = Path(tmp) / "hits.tsv"
    path.write_text("".join("\t".join(str(c) for c in r) + "\n" for r in rows))
    return path


class TestConflictRule(unittest.TestCase):
    def test_line_versus_ltr_is_a_conflict(self):
        self.assertTrue(S.is_conflict(
            "Class_I/LINE", "Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tat/Retand", 1))

    def test_class_i_versus_class_ii_is_a_conflict(self):
        self.assertTrue(S.is_conflict(
            "Class_I/LTR/Ty1_copia/Angela", "Class_II/Subclass_1/TIR/EnSpm_CACTA", 1))

    def test_te_versus_rdna_is_a_conflict(self):
        self.assertTrue(S.is_conflict("Class_I/LINE", "rDNA/45S_rDNA/18S", 1))

    def test_sibling_ltr_lineages_are_never_conflicts(self):
        """The single most destructive false positive available to this screen."""
        self.assertFalse(S.is_conflict(
            "Class_I/LTR/Ty1_copia/Ale", "Class_I/LTR/Ty1_copia/Angela", 1))
        self.assertFalse(S.is_conflict(
            "Class_I/LTR/Ty3_gypsy/chromovirus/Tekay",
            "Class_I/LTR/Ty3_gypsy/chromovirus/CRM", 1))

    def test_copia_versus_gypsy_is_not_a_conflict_at_the_default(self):
        self.assertFalse(S.is_conflict(
            "Class_I/LTR/Ty1_copia/Ale", "Class_I/LTR/Ty3_gypsy/chromovirus/CRM", 1))

    def test_ancestor_and_descendant_are_never_conflicts(self):
        self.assertFalse(S.is_conflict("Class_I/LTR", "Class_I/LTR/Ty1_copia", 1))
        self.assertFalse(S.is_conflict("Class_I/LTR/Ty1_copia", "Class_I/LTR", 1))
        self.assertFalse(S.is_conflict("Class_I/LINE", "Class_I/LINE", 1))

    def test_identical_and_empty_classes_are_safe(self):
        self.assertFalse(S.is_conflict("", "Class_I/LINE", 1))
        self.assertFalse(S.is_conflict("Class_I/LINE", "", 1))

    def test_depth_knob_widens_the_rule(self):
        self.assertFalse(S.is_conflict(
            "Class_I/LTR/Ty1_copia/Ale", "Class_I/LTR/Ty3_gypsy/chromovirus/CRM", 1))
        self.assertTrue(S.is_conflict(
            "Class_I/LTR/Ty1_copia/Ale", "Class_I/LTR/Ty3_gypsy/chromovirus/CRM", 2))


class TestBlastIdNormalisation(unittest.TestCase):
    def test_strips_the_class_suffix_blast_carries_through(self):
        self.assertEqual(
            S.blast_id_to_name("LINE_group_2147#Class_I/LINE"), "LINE_group_2147")
        self.assertEqual(S.blast_id_to_name("plain_name"), "plain_name")

    def test_unnormalised_ids_would_screen_nothing(self):
        """Regression: the first version looked classes up by un-normalised ids,
        found none, and reported a clean library on a contaminated one."""
        with tempfile.TemporaryDirectory() as tmp:
            records = [("chimera", "Class_I/LINE", "A" * 3000),
                       ("angela", "Class_I/LTR/Ty1_copia/Angela", "C" * 1000)]
            hits = hits_file(tmp, [
                ("chimera#Class_I/LINE", "angela#Class_I/LTR/Ty1_copia/Angela",
                 "95.0", 900, 3000, 1000, 2101, 3000),
            ])
            kept, audit = S.screen(records, hits, opts())
            self.assertEqual([r[4] for r in audit if r[0] == "chimera"], ["trimmed"])


class TestOwnership(unittest.TestCase):
    """A shared region belongs to whichever consensus it makes up more of."""

    def _pair(self, tmp, chimera_len, angela_len, hit_len):
        records = [("chimera", "Class_I/LINE", "A" * chimera_len),
                   ("angela", "Class_I/LTR/Ty1_copia/Angela", "C" * angela_len)]
        # the shared block sits at the 3' end of the chimera
        hits = hits_file(tmp, [
            ("chimera", "angela", "95.0", hit_len, chimera_len, angela_len,
             chimera_len - hit_len + 1, chimera_len),
            ("angela", "chimera", "95.0", hit_len, angela_len, chimera_len,
             1, hit_len),
        ])
        return records, hits

    def test_chimera_is_trimmed_and_its_victim_is_not(self):
        with tempfile.TemporaryDirectory() as tmp:
            records, hits = self._pair(tmp, chimera_len=8000, angela_len=2000,
                                       hit_len=1800)
            kept, audit = S.screen(records, hits, opts())
            actions = {r[0]: r[4] for r in audit}
            self.assertEqual(actions.get("chimera"), "trimmed")
            self.assertNotIn("angela", [r[0] for r in audit if r[4] != "kept"])
            lengths = {n: len(s) for n, _, s in kept}
            self.assertEqual(lengths["angela"], 2000)   # untouched
            self.assertLess(lengths["chimera"], 8000)

    def test_verdict_is_the_same_whichever_way_blast_reports_the_hit(self):
        with tempfile.TemporaryDirectory() as tmp:
            records, hits = self._pair(tmp, 8000, 2000, 1800)
            _, audit_a = S.screen(records, hits, opts())
            _, audit_b = S.screen(list(reversed(records)), hits, opts())
            self.assertEqual({r[0]: r[4] for r in audit_a},
                             {r[0]: r[4] for r in audit_b})

    def test_a_close_call_is_left_alone(self):
        """Comparable coverage on both sides -> no evidence which is wrong."""
        with tempfile.TemporaryDirectory() as tmp:
            records, hits = self._pair(tmp, chimera_len=2000, angela_len=2000,
                                       hit_len=900)
            kept, audit = S.screen(records, hits, opts())
            self.assertFalse([r for r in audit if r[4] in ("trimmed", "dropped")])
            self.assertEqual(len(kept), 2)


class TestTrimming(unittest.TestCase):
    def _one(self, tmp, length, blocks, **over):
        records = [("q", "Class_I/LINE", "A" * length),
                   ("s", "Class_II/Subclass_1/TIR/EnSpm_CACTA", "C" * 400)]
        rows = [("q", "s", "95.0", b - a + 1, length, 400, a, b) for a, b in blocks]
        return S.screen(records, hits_file(tmp, rows), opts(**over))

    def test_three_prime_block_is_trimmed(self):
        with tempfile.TemporaryDirectory() as tmp:
            kept, audit = self._one(tmp, 3000, [(2601, 3000)])
            self.assertEqual([r[4] for r in audit if r[0] == "q"], ["trimmed"])
            self.assertEqual(len(dict((n, s) for n, _, s in kept)["q"]), 2600)

    def test_five_prime_block_is_trimmed(self):
        with tempfile.TemporaryDirectory() as tmp:
            kept, audit = self._one(tmp, 3000, [(1, 400)])
            self.assertEqual(len(dict((n, s) for n, _, s in kept)["q"]), 2600)

    def test_internal_block_truncates_to_the_longer_clean_side(self):
        """Policy: better a truncated consensus than one carrying foreign sequence.

        An earlier version left internal conflicts in place, which on real data
        stranded 54% of the identified foreign base pairs inside the library.
        """
        with tempfile.TemporaryDirectory() as tmp:
            kept, audit = self._one(tmp, 3000, [(1400, 1800)])
            row = [r for r in audit if r[0] == "q"][0]
            self.assertEqual(row[4], "trimmed")
            # clean spans are 1-1399 (1399 bp) and 1801-3000 (1200 bp)
            self.assertEqual(len(dict((n, s) for n, _, s in kept)["q"]), 1399)

    def test_mosaic_of_blocks_keeps_only_the_largest_clean_run(self):
        """The real shape of a chimera: foreign fragments interleaved with rubble."""
        with tempfile.TemporaryDirectory() as tmp:
            kept, audit = self._one(
                tmp, 8000, [(2900, 3300), (4100, 4900), (5300, 6100), (6900, 8000)])
            # clean runs: 1-2899 (2899), 3301-4099 (799), 4901-5299, 6101-6899
            self.assertEqual(len(dict((n, s) for n, _, s in kept)["q"]), 2899)
            self.assertEqual([r[4] for r in audit if r[0] == "q"], ["trimmed"])

    def test_no_foreign_material_survives_a_trim(self):
        with tempfile.TemporaryDirectory() as tmp:
            blocks = [(500, 900), (2000, 2400)]
            kept, _ = self._one(tmp, 5000, blocks)
            seq_len = len(dict((n, s) for n, _, s in kept)["q"])
            start, end, _, _ = S.longest_clean_span(5000, blocks)
            self.assertEqual(seq_len, end - start + 1)
            for a, b in blocks:
                self.assertFalse(min(b, end) - max(a, start) + 1 > 0,
                                 "a conflict block survived inside the kept span")

    def test_blocks_at_both_ends_both_come_off(self):
        with tempfile.TemporaryDirectory() as tmp:
            kept, audit = self._one(tmp, 4000, [(1, 500), (3501, 4000)])
            self.assertEqual(len(dict((n, s) for n, _, s in kept)["q"]), 3000)

    def test_dropped_when_the_clean_remainder_is_tiny(self):
        """Discarding beats keeping a fragment too short to mask with."""
        with tempfile.TemporaryDirectory() as tmp:
            kept, audit = self._one(tmp, 3000, [(200, 2900)])
            self.assertEqual([r[4] for r in audit if r[0] == "q"], ["dropped"])

    def test_dropped_when_almost_nothing_survives(self):
        with tempfile.TemporaryDirectory() as tmp:
            kept, audit = self._one(tmp, 3000, [(1, 2900)])
            row = [r for r in audit if r[0] == "q"][0]
            self.assertEqual(row[4], "dropped")
            self.assertNotIn("q", [n for n, _, _ in kept])

    def test_weak_hits_are_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            kept, audit = self._one(tmp, 3000, [(2900, 3000)])          # 101 bp
            self.assertFalse([r for r in audit if r[4] != "kept"])
        with tempfile.TemporaryDirectory() as tmp:
            records = [("q", "Class_I/LINE", "A" * 3000),
                       ("s", "Class_II/Subclass_1/TIR/EnSpm_CACTA", "C" * 400)]
            rows = [("q", "s", "70.0", 400, 3000, 400, 2601, 3000)]     # 70% id
            kept, audit = S.screen(records, hits_file(tmp, rows), opts())
            self.assertFalse([r for r in audit if r[4] != "kept"])


class TestLongestCleanSpan(unittest.TestCase):
    def test_never_extends_past_the_sequence(self):
        for blocks in ([(1, 10)], [(90, 100)], [(1, 10), (90, 100)], [(40, 60)]):
            start, end, t5, t3 = S.longest_clean_span(100, blocks)
            self.assertGreaterEqual(start, 1)
            self.assertLessEqual(end, 100)
            self.assertGreaterEqual(t5, 0)
            self.assertGreaterEqual(t3, 0)

    def test_full_coverage_leaves_nothing(self):
        start, end, _, _ = S.longest_clean_span(100, [(1, 100)])
        self.assertLess(end, start)

    def test_picks_the_longest_of_several_clean_runs(self):
        start, end, _, _ = S.longest_clean_span(100, [(20, 25), (40, 45)])
        self.assertEqual((start, end), (46, 100))


class TestMaxConsensusLength(unittest.TestCase):
    def test_longest_prefix_wins(self):
        table = {"Class_I": 5000, "Class_I/LINE": 8000}
        self.assertEqual(S.max_length_for("Class_I/LINE", table), 8000)
        self.assertEqual(S.max_length_for("Class_I/LTR/Ty1_copia", table), 5000)
        self.assertEqual(S.max_length_for("Class_II/Subclass_1", table), 0)

    def test_vocabulary_ships_a_line_bound(self):
        table = S.load_max_consensus_length(str(REPO / "classification_vocabulary.yaml"))
        self.assertEqual(S.max_length_for("Class_I/LINE", table), 8000)

    def test_over_long_is_reported_never_trimmed(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [("big", "Class_I/LINE", "A" * 20000)]
            kept, audit = S.screen(records, None, opts())
            self.assertEqual(len(dict((n, s) for n, _, s in kept)["big"]), 20000)
            self.assertIn("exceeds max_consensus_length", audit[0][7])
            self.assertEqual(audit[0][4], "kept")


class TestCli(unittest.TestCase):
    script = str(REPO / "scripts" / "screen_library_cross_class.py")

    def test_disabled_copies_through(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "in.fasta"
            src.write_text(">a#Class_I/LINE\nACGT\n")
            dst = Path(tmp) / "out.fasta"
            subprocess.run([sys.executable, self.script, "-i", str(src),
                            "-o", str(dst), "-a", str(Path(tmp) / "audit.tsv"),
                            "--disabled"], check=True, capture_output=True)
            self.assertEqual(dst.read_text(), src.read_text())

    def test_empty_library_is_not_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "in.fasta"
            src.write_text("")
            dst = Path(tmp) / "out.fasta"
            audit = Path(tmp) / "audit.tsv"
            subprocess.run([sys.executable, self.script, "-i", str(src),
                            "-o", str(dst), "-a", str(audit)],
                           check=True, capture_output=True)
            self.assertEqual(dst.read_text(), "")
            self.assertTrue(audit.read_text().startswith("name\t"))

    def test_missing_blastn_copies_through_rather_than_failing(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "in.fasta"
            src.write_text(">a#Class_I/LINE\n" + "ACGT" * 100 + "\n")
            dst = Path(tmp) / "out.fasta"
            env = {"PATH": "/nonexistent", "HOME": tmp}
            r = subprocess.run([sys.executable, self.script, "-i", str(src),
                                "-o", str(dst), "-a", str(Path(tmp) / "audit.tsv"),
                                "-d", str(Path(tmp) / "wd")],
                               capture_output=True, text=True, env=env)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("WARNING", r.stderr)
            self.assertIn(">a#Class_I/LINE", dst.read_text())


class TestDeterminism(unittest.TestCase):
    def test_audit_is_sorted_and_output_order_follows_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [("zeta", "Class_I/LINE", "A" * 3000),
                       ("alpha", "Class_I/LINE", "A" * 3000),
                       ("s", "Class_II/Subclass_1/TIR/hAT", "C" * 400)]
            rows = [("zeta", "s", "95.0", 400, 3000, 400, 2601, 3000),
                    ("alpha", "s", "95.0", 400, 3000, 400, 1, 400)]
            kept, audit = S.screen(records, hits_file(tmp, rows), opts())
            self.assertEqual([n for n, _, _ in kept], ["zeta", "alpha", "s"])
            out = Path(tmp) / "audit.tsv"
            S.write_audit(audit, out)
            names = [l.split("\t")[0] for l in out.read_text().splitlines()[1:]]
            self.assertEqual(names, sorted(names))


if __name__ == "__main__":
    unittest.main(verbosity=2)

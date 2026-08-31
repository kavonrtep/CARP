#!/usr/bin/env python3
"""Unit tests for the DANTE_LINE consensus-boundary guards.

Covers the three L1 changes that stop `dante_line.py` building chimeric LINE
consensi out of the high-copy repeats that surround a LINE in a repeat-dense
genome:

  * repeatable ``--mask-gff3`` and the lightweight interval parser behind it,
  * ``cap_extensions`` -- the biological length bounds,
  * ``analyze_alignment_lengths`` -- the group-scaled support rule, including a
    parity check that the legacy settings reproduce the previous fixed-``k``
    implementation exactly.
"""

import math
import os
import sys
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import dante_line  # noqa: E402
from dante_line import (  # noqa: E402
    FeatureIndex,
    LineElement,
    analyze_alignment_lengths,
    cap_extensions,
    create_line_elements,
    parse_gff3_intervals,
)


ALIGNMENT_HEADER = (
    "query_id\tref_id\tquery_len\tref_len\tparasail_score\tmax_score\tmax_pos\t"
    "end_mode\tdegapped_query_len\tdegapped_ref_len\ttrimmed_query\ttrimmed_ref\t"
    "degapped_query\tdegapped_ref"
)


def write_alignment_tsv(path, pairs):
    """pairs: iterable of (query_id, ref_id, query_len, ref_len)."""
    with open(path, "w") as fh:
        fh.write(ALIGNMENT_HEADER + "\n")
        for q, r, ql, rl in pairs:
            fh.write(f"{q}\t{r}\t0\t0\t0\t0\t0\t5\t{ql}\t{rl}\tA\tA\tA\tA\n")


def read_result(path):
    rows = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        assert header == ["Group_ID", "Selected_Length", "Num_Shorter"]
        for line in fh:
            gid, sel, shorter = line.rstrip("\n").split("\t")
            rows[gid] = (int(sel), int(shorter))
    return rows


def legacy_analyze(alignment_tsv, output_tsv, min_num_alignments):
    """The pre-L1 implementation, kept here as the parity reference."""
    group_lengths = defaultdict(list)
    with open(alignment_tsv) as f:
        header = f.readline().strip().split("\t")
        qi = header.index("query_id")
        ri = header.index("ref_id")
        qli = header.index("degapped_query_len")
        rli = header.index("degapped_ref_len")
        for line in f:
            fields = line.strip().split("\t")
            group_lengths[fields[qi]].append(int(fields[qli]))
            group_lengths[fields[ri]].append(int(fields[rli]))

    results = []
    for group_id in sorted(group_lengths.keys()):
        lengths = sorted(group_lengths[group_id], reverse=True)
        if len(lengths) >= min_num_alignments:
            selected = lengths[min_num_alignments - 1]
            num_shorter = sum(1 for x in lengths if x <= selected) - 1
            results.append((group_id, selected, num_shorter))

    with open(output_tsv, "w") as f:
        f.write("Group_ID\tSelected_Length\tNum_Shorter\n")
        for gid, sel, shorter in results:
            f.write(f"{gid}\t{sel}\t{shorter}\n")


class TestParseGff3Intervals(unittest.TestCase):
    def test_reads_coordinates_and_skips_comments(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = Path(tmp) / "mask.gff3"
            gff.write_text(
                "##gff-version 3\n"
                "chr1\tDANTE_LTR\ttransposable_element\t100\t500\t.\t+\t.\tID=te1\n"
                "\n"
                "chr2\tDANTE_TIR\tsequence_feature\t7\t9\t.\t-\t.\tID=tir1\n"
            )
            self.assertEqual(
                parse_gff3_intervals(str(gff)),
                [("chr1", 100, 500), ("chr2", 7, 9)],
            )

    def test_bad_rows_are_skipped_not_fatal(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = Path(tmp) / "mask.gff3"
            gff.write_text(
                "chr1\tx\ty\tNOT_A_NUMBER\t500\t.\t+\t.\tID=a\n"
                "chr1\ttoo\tfew\tcols\n"
                "chr1\tx\ty\t10\t20\t.\t+\t.\tID=b\n"
            )
            self.assertEqual(parse_gff3_intervals(str(gff)), [("chr1", 10, 20)])


class TestMaskIndex(unittest.TestCase):
    def test_from_intervals_matches_feature_objects(self):
        rows = [("chr1", 10, 20), ("chr1", 100, 200), ("chr2", 5, 8)]
        features = [
            dante_line.GFF3Feature(sn, "s", "f", start, end, ".", "+", ".", "")
            for sn, start, end in rows
        ]
        from_features = FeatureIndex(features)
        from_intervals = FeatureIndex.from_intervals(rows)
        for seqname, ref, floor, ceil in (
            ("chr1", 150, 0, 10_000),
            ("chr1", 5, 0, 10_000),
            ("chr2", 7, 0, 10_000),
            ("missing", 50, 3, 999),
        ):
            self.assertEqual(
                from_features.nearest_end_below(seqname, ref, floor),
                from_intervals.nearest_end_below(seqname, ref, floor),
            )
            self.assertEqual(
                from_features.nearest_start_above(seqname, ref, ceil),
                from_intervals.nearest_start_above(seqname, ref, ceil),
            )

    def test_merged_masks_clip_to_whichever_comes_first(self):
        """A second mask file must be able to tighten the flank, never loosen it."""
        tidehunter = [("chr1", 1_000, 1_100)]
        dante_ltr = [("chr1", 1_500, 4_000)]
        one = FeatureIndex.from_intervals(tidehunter)
        both = FeatureIndex.from_intervals(tidehunter + dante_ltr)
        # element starts at 5000; flank would reach back to 1 without masks
        self.assertEqual(one.nearest_end_below("chr1", 5_000, 1), 1_101)
        self.assertEqual(both.nearest_end_below("chr1", 5_000, 1), 4_001)


class TestCapExtensions(unittest.TestCase):
    def test_disabled_by_zero(self):
        self.assertEqual(cap_extensions(2000, 9999, 9999, 0, 0), (9999, 9999))

    def test_per_side_cap(self):
        self.assertEqual(
            cap_extensions(2000, 9999, 400, max_extension=2500, max_element_length=0),
            (2500, 400),
        )

    def test_element_length_cap_splits_evenly(self):
        # core 2000, budget 8000-2000 = 6000, both sides want 2500 -> fits
        self.assertEqual(
            cap_extensions(2000, 9999, 9999, 2500, 8000), (2500, 2500)
        )
        # core 4000, budget 4000, both sides want 2500 -> split 2000/2000
        self.assertEqual(
            cap_extensions(4000, 9999, 9999, 2500, 8000), (2000, 2000)
        )

    def test_per_side_cap_is_applied_before_the_element_budget(self):
        # core 5000 -> budget 3000, but --max-extension 2500 binds the 3' side
        # first and 200+2500 already fits, so the budget never engages.
        self.assertEqual(cap_extensions(5000, 200, 9999, 2500, 8000), (200, 2500))

    def test_modest_side_keeps_what_it_asked_and_yields_the_rest(self):
        # core 6000 -> budget 2000; after the per-side cap the sides want
        # 200 + 2500 = 2700 > 2000. 5' is below half (1000) so it keeps 200
        # and the 3' side takes the remaining 1800.
        self.assertEqual(cap_extensions(6000, 200, 9999, 2500, 8000), (200, 1800))
        # mirrored
        self.assertEqual(cap_extensions(6000, 9999, 200, 2500, 8000), (1800, 200))

    def test_odd_budget_goes_to_three_prime(self):
        self.assertEqual(cap_extensions(4999, 9999, 9999, 2500, 8000), (1500, 1501))

    def test_core_at_or_over_the_cap_gets_no_extension_and_is_never_trimmed(self):
        self.assertEqual(cap_extensions(9000, 4000, 4000, 2500, 8000), (0, 0))

    def test_never_lengthens(self):
        for core in (500, 2000, 5000):
            for e5 in (0, 100, 3000, 9999):
                for e3 in (0, 100, 3000, 9999):
                    c5, c3 = cap_extensions(core, e5, e3, 2500, 8000)
                    self.assertLessEqual(c5, e5)
                    self.assertLessEqual(c3, e3)
                    self.assertGreaterEqual(c5, 0)
                    self.assertGreaterEqual(c3, 0)


class TestCreateLineElements(unittest.TestCase):
    def _pattern(self, strand, start, end):
        feat = dante_line.GFF3Feature("chr1", "dante", "protein_domain", start, end,
                                      ".", strand, ".", "Name=RT")
        return dante_line.FeatureGroup("LINE_group_0001", "chr1", strand, [feat], "ENDO-RT")

    def test_records_inferred_and_applied_extensions(self):
        patterns = [self._pattern("+", 10_000, 12_000)]
        lengths = {"LINE_group_0001": {"5prime": 9999, "3prime": 9999}}
        (element,) = create_line_elements(patterns, lengths,
                                          max_extension=2500, max_element_length=8000)
        self.assertEqual(element.extension_5prime_inferred, 9999)
        self.assertEqual(element.extension_3prime_inferred, 9999)
        # core is 2001 bp, budget 5999, each side capped at 2500 first
        self.assertEqual((element.extension_5prime, element.extension_3prime), (2500, 2500))
        self.assertEqual(element.end - element.start + 1, 2001 + 5000)

    def test_minus_strand_extensions_stay_on_the_right_side(self):
        patterns = [self._pattern("-", 10_000, 12_000)]
        lengths = {"LINE_group_0001_revcomp": {"5prime": 1000, "3prime": 300}}
        (element,) = create_line_elements(patterns, lengths, 2500, 8000)
        # 5' of a minus-strand element is at the higher coordinate
        self.assertEqual(element.end, 12_000 + 1000)
        self.assertEqual(element.start, 10_000 - 300)

    def test_uncapped_run_leaves_extensions_untouched(self):
        patterns = [self._pattern("+", 10_000, 12_000)]
        lengths = {"LINE_group_0001": {"5prime": 400, "3prime": 700}}
        (element,) = create_line_elements(patterns, lengths, 2500, 8000)
        self.assertEqual(element.extension_5prime, element.extension_5prime_inferred)
        self.assertEqual(element.extension_3prime, element.extension_3prime_inferred)


class TestSupportRule(unittest.TestCase):
    def run_analyze(self, pairs, **kwargs):
        with tempfile.TemporaryDirectory() as tmp:
            aln = Path(tmp) / "aln.tsv"
            out = Path(tmp) / "out.tsv"
            write_alignment_tsv(aln, pairs)
            analyze_alignment_lengths(aln, out, kwargs.pop("min_num_alignments", 3), **kwargs)
            return read_result(out)

    def test_legacy_settings_reproduce_the_previous_implementation(self):
        """support_fraction=0, min_group_alignments=0 must be byte-identical."""
        pairs = [
            ("g1", "g2", 900, 850), ("g1", "g3", 800, 700), ("g1", "g4", 100, 90),
            ("g2", "g3", 500, 500), ("g2", "g4", 400, 300), ("g3", "g4", 60, 50),
            ("g5", "g6", 10, 10),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            aln = Path(tmp) / "aln.tsv"
            new_out = Path(tmp) / "new.tsv"
            old_out = Path(tmp) / "old.tsv"
            write_alignment_tsv(aln, pairs)
            analyze_alignment_lengths(aln, new_out, 3,
                                      support_fraction=0.0, min_group_alignments=0)
            legacy_analyze(aln, old_out, 3)
            self.assertEqual(new_out.read_bytes(), old_out.read_bytes())

    def test_a_tiny_group_no_longer_sets_its_own_boundary(self):
        """The bug: with k fixed at 3, the 3rd largest of 3 IS the minimum."""
        pairs = [("g1", "p1", 9998, 9998), ("g1", "p2", 9995, 9995), ("g1", "p3", 9935, 9935)]
        legacy = self.run_analyze(pairs, support_fraction=0.0, min_group_alignments=0)
        self.assertEqual(legacy["g1"][0], 9935)  # the minimum -- no filtering at all
        guarded = self.run_analyze(pairs, support_fraction=0.5, min_group_alignments=5)
        self.assertNotIn("g1", guarded)  # too few partners -> core only

    def test_outlier_pairs_cannot_set_the_boundary_of_a_large_group(self):
        """Three long partners out of many must not carry the whole group."""
        pairs = [("g1", f"p{i}", 6000, 6000) for i in range(3)]
        pairs += [("g1", f"p{i}", 750, 750) for i in range(3, 51)]
        legacy = self.run_analyze(pairs, support_fraction=0.0, min_group_alignments=0)
        self.assertEqual(legacy["g1"][0], 6000)
        guarded = self.run_analyze(pairs, support_fraction=0.5, min_group_alignments=5)
        self.assertEqual(guarded["g1"][0], 750)

    def test_a_genuinely_supported_extension_survives(self):
        """If the whole group agrees, the long extension is kept."""
        pairs = [("g1", f"p{i}", 6000, 6000) for i in range(20)]
        guarded = self.run_analyze(pairs, support_fraction=0.5, min_group_alignments=5)
        self.assertEqual(guarded["g1"][0], 6000)

    def test_k_is_the_ceiling_of_the_fraction(self):
        # 5 alignments, fraction 0.5 -> k = 3 -> 3rd largest
        pairs = [("g1", f"p{i}", length, length)
                 for i, length in enumerate([500, 400, 300, 200, 100])]
        got = self.run_analyze(pairs, support_fraction=0.5, min_group_alignments=5)
        self.assertEqual(got["g1"][0], 300)
        self.assertEqual(math.ceil(5 * 0.5), 3)

    def test_min_num_alignments_is_a_floor_on_k(self):
        # 6 alignments at fraction 0.1 -> ceil(0.6)=1, but the floor keeps k=3
        pairs = [("g1", f"p{i}", length, length)
                 for i, length in enumerate([600, 500, 400, 300, 200, 100])]
        got = self.run_analyze(pairs, min_num_alignments=3, support_fraction=0.1,
                               min_group_alignments=0)
        self.assertEqual(got["g1"][0], 400)

    def test_rejects_out_of_range_fraction(self):
        with self.assertRaises(ValueError):
            self.run_analyze([("g1", "g2", 10, 10)], support_fraction=1.5)

    def test_output_is_sorted_and_deterministic(self):
        pairs = [("gz", "ga", 100, 200), ("gm", "ga", 300, 400), ("ga", "gz", 50, 60)]
        first = self.run_analyze(pairs, support_fraction=0.5, min_group_alignments=0)
        second = self.run_analyze(pairs, support_fraction=0.5, min_group_alignments=0)
        self.assertEqual(first, second)
        with tempfile.TemporaryDirectory() as tmp:
            aln = Path(tmp) / "a.tsv"
            out = Path(tmp) / "o.tsv"
            write_alignment_tsv(aln, pairs)
            analyze_alignment_lengths(aln, out, 3, support_fraction=0.5)
            names = [l.split("\t")[0] for l in out.read_text().splitlines()[1:]]
            self.assertEqual(names, sorted(names))

    def test_writes_atomically(self):
        pairs = [("g1", "g2", 10, 10), ("g1", "g3", 20, 20), ("g2", "g3", 30, 30)]
        with tempfile.TemporaryDirectory() as tmp:
            aln = Path(tmp) / "a.tsv"
            out = Path(tmp) / "o.tsv"
            write_alignment_tsv(aln, pairs)
            analyze_alignment_lengths(aln, out, 3)
            self.assertTrue(out.exists())
            self.assertFalse((out.parent / (out.name + ".tmp")).exists())


class TestSinglePassWhenFractionRuleIsOff(unittest.TestCase):
    """The counting pass must not run when k is a constant.

    `dante_tir_fallback` shares `analyze_alignment_lengths`, and its alignment
    tables are far larger than DANTE_LINE's -- on a 14.5 Gb wheat run the
    EnSpm_CACTA 3' table alone is 43 GB, because every row carries its aligned
    sequences. An unconditional second read costs tens of GB of I/O for nothing.
    """

    def _count_opens(self, path, **kwargs):
        import builtins
        real_open = builtins.open
        seen = {"n": 0}

        def counting(f, *args, **kw):
            if str(f) == str(path):
                seen["n"] += 1
            return real_open(f, *args, **kw)

        builtins.open = counting
        try:
            with tempfile.TemporaryDirectory() as tmp:
                dante_line.analyze_alignment_lengths(
                    path, Path(tmp) / "out.tsv", 3, **kwargs)
        finally:
            builtins.open = real_open
        return seen["n"]

    def test_read_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            aln = Path(tmp) / "aln.tsv"
            write_alignment_tsv(aln, [("g1", f"p{i}", 900 - i * 10, 900 - i * 10)
                                      for i in range(12)])
            self.assertEqual(self._count_opens(aln), 1)
            self.assertEqual(
                self._count_opens(aln, support_fraction=0.5, min_group_alignments=5), 2)

    def test_both_paths_agree_where_they_should(self):
        """min_group_alignments must behave identically on both code paths."""
        with tempfile.TemporaryDirectory() as tmp:
            aln = Path(tmp) / "aln.tsv"
            # g1 has 6 alignments, g2 has 2
            rows = [("g1", f"p{i}", 500, 500) for i in range(6)]
            rows += [("g2", f"q{i}", 400, 400) for i in range(2)]
            write_alignment_tsv(aln, rows)
            out = Path(tmp) / "o.tsv"
            analyze_alignment_lengths(aln, out, 3, min_group_alignments=5)
            names = {l.split("\t")[0] for l in out.read_text().splitlines()[1:]}
            self.assertIn("g1", names)
            self.assertNotIn("g2", names)   # under-supported, single-pass path


class TestFallbackSharesTheGuards(unittest.TestCase):
    """dante_tir_fallback runs the same engine and must carry the same guards."""

    def setUp(self):
        import dante_tir_fallback
        self.fallback = dante_tir_fallback

    def test_it_imports_the_shared_cap(self):
        self.assertIs(self.fallback.cap_extensions, dante_line.cap_extensions)

    def test_length_bound_is_per_tir_superfamily(self):
        """CACTA genuinely reaches ~20 kb; Tc1_Mariner cannot."""
        cacta = self.fallback.max_element_length_for_subtype("EnSpm/CACTA")
        mariner = self.fallback.max_element_length_for_subtype("Tc1/Mariner")
        self.assertGreater(cacta, mariner)
        self.assertEqual(cacta, 30000)
        self.assertEqual(mariner, 10000)

    def test_unknown_subtype_falls_back_to_the_generic_tir_bound(self):
        self.assertEqual(
            self.fallback.max_element_length_for_subtype("Something_New"), 25000)

    def test_explicit_value_overrides_the_vocabulary(self):
        self.assertEqual(
            self.fallback.max_element_length_for_subtype("EnSpm/CACTA", 4321), 4321)


if __name__ == "__main__":
    unittest.main(verbosity=2)

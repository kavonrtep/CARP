#!/usr/bin/env python3
"""Per-superfamily / per-side flank bounds (max_extension_per_side).

The regression these guard: DANTE_LINE and DANTE_TIR_FALLBACK shared ONE
symmetric extension cap (1500 bp). Measured against 13,760 element-sides from
DANTE_TIR primary successes -- elements with tir_seq5 + tir_seq3 + tsd, so the
boundary is exact -- the median true flank spans 6x between superfamilies
(PIF_Harbinger 951 bp, EnSpm_CACTA 3544 bp). A shared 1500 therefore cut into
the MEDIAN real element for four of the five superfamilies, while being four
times too loose for Harbinger. LINE is asymmetric for a structural reason: the
core is ORF2, so ORF1 + 5'UTR sit outside it 5' but only 3'UTR + polyA 3'.
"""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import classification
from dante_line import cap_extensions, line_max_extension, line_core_to_3prime_end
from dante_tir_fallback import max_extension_for_subtype

# medians measured from DANTE_TIR successes (wheat, Lycopus, Boechera)
MEASURED_MEDIAN = {
    "EnSpm/CACTA": 3544, "MuDR/Mutator": 2558, "Tc1/Mariner": 2099,
    "hAT": 1515, "PIF/Harbinger": 951,
}


class TestVocabularyLookup(unittest.TestCase):
    def test_line_bound_is_asymmetric(self):
        # MEASURED 2026-09-01 from 218 confirmed loci across 7 genomes; the
        # earlier structure-derived 2000/800 truncated 57% and 72% of them.
        self.assertEqual(classification.max_extension_for_class("Class_I/LINE", "5prime"), 3400)
        # 3' is a backstop only -- max_core_to_3prime_end is the real bound.
        self.assertEqual(classification.max_extension_for_class("Class_I/LINE", "3prime"), 2500)

    def test_longest_prefix_wins_over_generic(self):
        generic = classification.max_extension_for_class("Class_II/Subclass_1/TIR", "5prime")
        cacta = classification.max_extension_for_class(
            "Class_II/Subclass_1/TIR/EnSpm_CACTA", "5prime")
        self.assertEqual(generic, 6000)
        self.assertEqual(cacta, 5800)
        self.assertNotEqual(cacta, generic)

    def test_unlisted_superfamily_falls_back_to_generic(self):
        self.assertEqual(
            classification.max_extension_for_class(
                "Class_II/Subclass_1/TIR/Novosib", "5prime"), 6000)

    def test_layers_we_do_not_extend_are_unbounded(self):
        # LTR-RTs are found structurally by DANTE_LTR, never flank-extended.
        self.assertEqual(classification.max_extension_for_class("Class_I/LTR/Ty1_copia", "5prime"), 0)

    def test_rejects_bad_side(self):
        with self.assertRaises(ValueError):
            classification.max_extension_for_class("Class_I/LINE", "middle")

    def test_per_side_bound_never_exceeds_whole_element_bound(self):
        """The two tables must not contradict: two flanks + a core must fit."""
        vocab = classification.load_vocabulary()
        for cls, val in vocab.max_extension_per_side.items():
            sides = val.values() if isinstance(val, dict) else [val]
            whole = 0
            for prefix, w in vocab.max_consensus_length.items():
                if cls == prefix or cls.startswith(prefix + "/"):
                    whole = max(whole, w)
            if whole:
                self.assertLessEqual(
                    2 * max(int(x) for x in sides), whole,
                    f"{cls}: two flanks alone exceed max_consensus_length")


class TestCapExtensions(unittest.TestCase):
    def test_symmetric_cap_still_works(self):
        self.assertEqual(cap_extensions(2000, 5000, 5000, max_extension=1500), (1500, 1500))

    def test_per_side_caps_are_independent(self):
        self.assertEqual(
            cap_extensions(2000, 5000, 5000, max_ext_5prime=2000, max_ext_3prime=800),
            (2000, 800))

    def test_per_side_overrides_symmetric_on_the_side_it_is_set(self):
        self.assertEqual(
            cap_extensions(2000, 5000, 5000, max_extension=1200, max_ext_5prime=2000),
            (2000, 1200))

    def test_zero_disables_a_side(self):
        self.assertEqual(cap_extensions(2000, 5000, 4000, max_ext_5prime=0, max_ext_3prime=0),
                         (5000, 4000))

    def test_element_budget_still_applies_after_per_side_caps(self):
        ext5, ext3 = cap_extensions(2000, 5000, 5000, max_ext_5prime=2000,
                                    max_ext_3prime=800, max_element_length=3500)
        self.assertEqual(ext5 + ext3, 1500)

    def test_core_is_never_trimmed(self):
        self.assertEqual(cap_extensions(9000, 500, 500, max_element_length=8000), (0, 0))


class TestTheRegression(unittest.TestCase):
    def test_shared_1500_truncated_the_median_real_element(self):
        """This is why the shared cap had to go -- four of five superfamilies."""
        truncated = [sf for sf, med in MEASURED_MEDIAN.items() if med > 1500]
        self.assertEqual(len(truncated), 4)
        for sf in truncated:
            got5, _ = cap_extensions(1300, MEASURED_MEDIAN[sf], 0, max_extension=1500)
            self.assertLess(got5, MEASURED_MEDIAN[sf])

    def test_superfamily_bounds_admit_the_median_real_element(self):
        for sf, med in MEASURED_MEDIAN.items():
            cap = max_extension_for_subtype(sf, "5prime")
            self.assertGreaterEqual(
                cap, med, f"{sf}: bound {cap} truncates the median real flank {med}")

    def test_harbinger_is_not_given_cacta_room(self):
        self.assertLess(max_extension_for_subtype("PIF/Harbinger", "5prime"),
                        max_extension_for_subtype("EnSpm/CACTA", "5prime"))


class TestResolverPrecedence(unittest.TestCase):
    def test_explicit_cli_value_wins_for_line(self):
        self.assertEqual(line_max_extension("5prime", 1500), 1500)
        self.assertEqual(line_max_extension("3prime", 1500), 1500)

    def test_explicit_cli_value_wins_for_tir(self):
        self.assertEqual(max_extension_for_subtype("EnSpm/CACTA", "5prime", 999), 999)

    def test_default_reads_the_vocabulary(self):
        self.assertEqual(line_max_extension("5prime"), 3400)
        self.assertEqual(line_max_extension("3prime"), 2500)

    def test_subtype_label_is_sanitised(self):
        """DANTE writes 'EnSpm/CACTA'; the vocabulary key is 'EnSpm_CACTA'."""
        self.assertEqual(max_extension_for_subtype("EnSpm/CACTA", "5prime"), 5800)


class TestCoreToThreePrimeSpan(unittest.TestCase):
    """The 3' bound is on (core + extension), not on the extension alone.

    Measured over 148 confirmed LINE loci: Pearson r(core length, 3' extension)
    = -0.936, i.e. they are nearly a constant sum. The extension alone varies
    13.8x (p90/p10); the span varies 1.16x. A flat per-side 3' cap therefore
    limits a quantity that is not stable.
    """

    # two real Boechera elements that end ~3750 bp after their ENDO start
    SHORT_CORE = (2164, 1568)   # short core, long leftover extension
    LONG_CORE = (3893, 155)     # long core, almost nothing left over

    def _cap(self, core, want):
        return cap_extensions(core, 0, want, max_ext_5prime=3400,
                              max_ext_3prime=2500, max_core_to_3prime_end=4500)[1]

    def test_both_real_shapes_survive_the_span_rule(self):
        for core, ext in (self.SHORT_CORE, self.LONG_CORE):
            self.assertEqual(self._cap(core, ext), ext,
                             f"core {core} + ext {ext} should not be truncated")

    def test_flat_800_cap_would_have_truncated_the_short_core_element(self):
        """The regression: the shipped flat cap clipped 72% of confirmed loci."""
        core, ext = self.SHORT_CORE
        flat = cap_extensions(core, 0, ext, max_ext_3prime=800)[1]
        self.assertEqual(flat, 800)
        self.assertLess(flat, ext)

    def test_span_binds_when_the_inference_overruns(self):
        # long core asking for far more than the span allows
        self.assertEqual(self._cap(3893, 1500), 4500 - 3893)

    def test_core_longer_than_the_span_yields_no_extension(self):
        self.assertEqual(self._cap(4600, 500), 0)

    def test_core_is_never_trimmed_by_the_span_rule(self):
        e5, e3 = cap_extensions(4600, 300, 300, max_core_to_3prime_end=4500)
        self.assertEqual(e3, 0)
        self.assertEqual(e5, 300)      # the 5' side is untouched by this bound

    def test_explicit_max_extension_is_a_complete_override(self):
        """--max-extension must not be quietly tightened by a vocabulary bound."""
        self.assertEqual(line_core_to_3prime_end(2500), 0)
        self.assertEqual(line_core_to_3prime_end(), 4500)

    def test_vocabulary_is_line_only(self):
        self.assertEqual(
            classification.max_core_to_3prime_end_for_class(
                "Class_II/Subclass_1/TIR/EnSpm_CACTA"), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)

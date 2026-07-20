#!/usr/bin/env python3
"""Test scripts/cleanup_outputs.py: the right scratch is deleted, and nothing in
the keep-set (manifest outputs, their symlink targets, CI/count-checked files,
run metadata) is ever touched — for minimal, maximal, dry-run, and none."""
import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from cleanup_outputs import cleanup  # noqa: E402
from assert_cleanup import check_cleanup  # noqa: E402
from manifest import OUTPUTS  # noqa: E402


def touch(p: Path, content: str = "x"):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(content)


def build_tree(root: Path):
    # ── keep: manifest outputs / metadata / CI-checked ───────────────────────
    touch(root / "Repeat_Annotation_Unified.gff3")
    touch(root / "summary_statistics.csv")
    touch(root / "carp_manifest.json")
    touch(root / "run_provenance.json")
    touch(root / ".classifications_validated")
    touch(root / "DANTE_TIR/DANTE_TIR_final.fasta")     # CI count-checked
    touch(root / "DANTE_TIR/DANTE_TIR_combined.gff3")   # CI count-checked
    touch(root / "DANTE_LINE/DANTE_LINE.gff3")          # CI count-checked
    # a manifest directory output with a file inside
    touch(root / "Repeat_density_by_class_bigwig/100k/LTR_RT_TR_100k.bw")
    # top-level deliverable as a SYMLINK into a subdir (container layout)
    touch(root / "DANTE_LTR/DANTE_LTR.gff3")
    os.symlink("DANTE_LTR/DANTE_LTR.gff3", root / "DANTE_LTR_top.link")
    # LTR_RTs_library.fasta symlinked to its library target (both must survive)
    touch(root / "DANTE_LTR/library/mmseqs2/mmseqs_representative_seq_clean.fasta")
    os.symlink("library/mmseqs2/mmseqs_representative_seq_clean.fasta",
               root / "DANTE_LTR/LTR_RTs_library.fasta")
    # a scratch file that ALSO happens to be a top-level symlink target — the
    # keep-set must protect it even though the TE*.fasta glob matches it.
    touch(root / "DANTE_LTR/library/TE_special.fasta")
    os.symlink("DANTE_LTR/library/TE_special.fasta", root / "special.link")

    # ── delete (minimal) ─────────────────────────────────────────────────────
    touch(root / "DANTE/DANTE_filtered.gff3.tmp.gff3")
    touch(root / "DANTE_TIR/DANTE_TIR.RData")
    touch(root / "DANTE_TIR_FALLBACK/TPase_5prime_alignment.tsv")
    touch(root / "DANTE_TIR_FALLBACK/TPase_3prime_alignment.tsv")
    touch(root / "RepeatMasker/genome_cleaned.fasta")
    touch(root / "TideCluster/genome_cleaned.fasta")
    touch(root / "DANTE_LTR/library/mmseqs2/mmseqs_all_seqs.fasta")
    touch(root / "DANTE_LTR/library/mmseqs2/partitioned_s900_w1000.fasta")
    touch(root / "DANTE_LTR/library/TE_all.fasta")      # matched by TE*.fasta glob
    touch(root / "DANTE_LTR/LTR_RTs_library.fasta.reformatted")

    # ── delete (maximal only) ────────────────────────────────────────────────
    touch(root / "TideCluster/default/TideCluster_kite/monomer.csv")
    touch(root / "TideCluster/short_monomer/TideCluster_tarean/x.kmers")
    touch(root / "RepeatMasker/workdir/chunk1.fa")
    touch(root / "DANTE_TIR/mmseqs_combined/tmp/db")
    touch(root / "DANTE_LINE/mmseqs/clu.tsv")


KEEP = [
    "Repeat_Annotation_Unified.gff3", "summary_statistics.csv",
    "carp_manifest.json", "run_provenance.json", ".classifications_validated",
    "DANTE_TIR/DANTE_TIR_final.fasta", "DANTE_TIR/DANTE_TIR_combined.gff3",
    "DANTE_LINE/DANTE_LINE.gff3",
    "Repeat_density_by_class_bigwig/100k/LTR_RT_TR_100k.bw",
    "DANTE_LTR/DANTE_LTR.gff3",
    "DANTE_LTR/library/mmseqs2/mmseqs_representative_seq_clean.fasta",
    "DANTE_LTR/LTR_RTs_library.fasta",
    "DANTE_LTR/library/TE_special.fasta",               # protected via symlink
]
MINIMAL_GONE = [
    "DANTE/DANTE_filtered.gff3.tmp.gff3", "DANTE_TIR/DANTE_TIR.RData",
    "DANTE_TIR_FALLBACK/TPase_5prime_alignment.tsv",
    "DANTE_TIR_FALLBACK/TPase_3prime_alignment.tsv",
    "RepeatMasker/genome_cleaned.fasta", "TideCluster/genome_cleaned.fasta",
    "DANTE_LTR/library/mmseqs2/mmseqs_all_seqs.fasta",
    "DANTE_LTR/library/mmseqs2/partitioned_s900_w1000.fasta",
    "DANTE_LTR/library/TE_all.fasta",
    "DANTE_LTR/LTR_RTs_library.fasta.reformatted",
]
MAXIMAL_GONE = [
    "TideCluster/default/TideCluster_kite",
    "TideCluster/short_monomer/TideCluster_tarean",
    "RepeatMasker/workdir", "DANTE_TIR/mmseqs_combined", "DANTE_LINE/mmseqs",
]


class CleanupTest(unittest.TestCase):
    def setUp(self):
        self.d = Path(tempfile.mkdtemp(prefix="cleanup_test_"))
        build_tree(self.d)

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def _exists(self, rel):
        return (self.d / rel).exists() or (self.d / rel).is_symlink()

    def assert_keep_all(self):
        for rel in KEEP:
            self.assertTrue(self._exists(rel), f"keep-list file deleted: {rel}")

    def test_minimal(self):
        cleanup(self.d, "minimal", log=lambda m: None)
        self.assert_keep_all()
        for rel in MINIMAL_GONE:
            self.assertFalse(self._exists(rel), f"minimal did not delete: {rel}")
        for rel in MAXIMAL_GONE:               # maximal-only survives minimal
            self.assertTrue(self._exists(rel), f"minimal wrongly deleted: {rel}")

    def test_maximal(self):
        cleanup(self.d, "maximal", log=lambda m: None)
        self.assert_keep_all()
        for rel in MINIMAL_GONE + MAXIMAL_GONE:
            self.assertFalse(self._exists(rel), f"maximal did not delete: {rel}")

    def test_dry_run_deletes_nothing(self):
        cleanup(self.d, "maximal", dry_run=True, log=lambda m: None)
        self.assert_keep_all()
        for rel in MINIMAL_GONE + MAXIMAL_GONE:
            self.assertTrue(self._exists(rel), f"dry-run deleted: {rel}")

    def test_none_is_noop(self):
        cleanup(self.d, "none", log=lambda m: None)
        self.assert_keep_all()
        for rel in MINIMAL_GONE + MAXIMAL_GONE:
            self.assertTrue(self._exists(rel), f"none deleted: {rel}")

    def test_unknown_mode_is_noop(self):
        cleanup(self.d, "bogus", log=lambda m: None)
        for rel in MINIMAL_GONE:
            self.assertTrue(self._exists(rel), f"unknown mode deleted: {rel}")


class AssertCleanupTest(unittest.TestCase):
    """scripts/assert_cleanup.py — the CI check that cleanup ran (scratch gone,
    every manifest output kept), reading the applied mode from provenance."""

    def setUp(self):
        self.d = Path(tempfile.mkdtemp(prefix="assert_cleanup_test_"))
        build_tree(self.d)
        # create every manifest output so the presence check has a full tree
        for rel in OUTPUTS.values():
            p = self.d / rel.rstrip("/")
            if rel.endswith("/"):
                p.mkdir(parents=True, exist_ok=True)
            else:
                p.parent.mkdir(parents=True, exist_ok=True)
                if not p.exists():
                    p.write_text("x")
        # provenance recording the mode that "ran" (overwrites build_tree's stub)
        (self.d / "run_provenance.json").write_text(
            json.dumps({"config": {"cleanup_intermediates": "maximal"}}))

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def test_pass_after_cleanup(self):
        cleanup(self.d, "maximal", log=lambda m: None)
        mode, errs = check_cleanup(self.d)
        self.assertEqual(mode, "maximal")
        self.assertEqual(errs, [], f"clean tree flagged: {errs}")

    def test_fail_if_scratch_left(self):
        cleanup(self.d, "maximal", log=lambda m: None)
        (self.d / "DANTE_TIR").mkdir(exist_ok=True)
        (self.d / "DANTE_TIR/DANTE_TIR.RData").write_text("x")  # cleanup "missed" it
        _, errs = check_cleanup(self.d)
        self.assertTrue(any("left scratch behind" in e for e in errs), errs)

    def test_fail_if_output_deleted(self):
        cleanup(self.d, "maximal", log=lambda m: None)
        (self.d / "Repeat_Annotation_Unified.gff3").unlink()  # a deliverable lost
        _, errs = check_cleanup(self.d)
        self.assertTrue(any("manifest output missing" in e for e in errs), errs)


if __name__ == "__main__":
    r = unittest.main(exit=False, verbosity=0).result
    if r.wasSuccessful():
        print("test_cleanup_outputs: PASSED")
        sys.exit(0)
    sys.exit(1)

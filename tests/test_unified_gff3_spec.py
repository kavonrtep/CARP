#!/usr/bin/env python3
"""Format-drift guard for Repeat_Annotation_Unified.gff3.

Exercises scripts/validate_unified_gff3.py (the executable mirror of
docs/unified_annotation_gff3_spec.md): a synthetic, contract-conforming GFF3
covering every feature kind must pass, and each single-attribute drift — most
importantly a TideCluster satellite whose Name is not the bare TRC_<n> (the
regression that silently emptied the per-family BigWig outputs) — must be
caught. Runnable standalone: `python tests/test_unified_gff3_spec.py`
(exit 0 = pass).
"""
import copy
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from validate_unified_gff3 import validate  # noqa: E402


def feat(seqid, source, ftype, start, end, strand, **attrs):
    return {"cols": [seqid, source, ftype, str(start), str(end), ".", strand, "."],
            "attrs": attrs}


# One conforming feature per kind. IDs are referenced by the Level-2 Parent below.
GOOD = [
    feat("OZ1", "DANTE_LTR", "transposable_element", 100, 5000, "+",
         ID="UA_L1_00000001", Name="Class_I/LTR/Ty1_copia/Ale",
         classification="Class_I/LTR/Ty1_copia/Ale",
         source_tier="1", source_tool="DANTE_LTR", element_type="complete"),
    feat("OZ1", "DANTE_LTR", "transposable_element", 6000, 9000, "-",
         ID="UA_L1_00000002", Name="Class_I/LTR/Ty1_copia/Ale",
         classification="Class_I/LTR/Ty1_copia/Ale",
         source_tier="1", source_tool="DANTE_LTR", element_type="partial"),
    feat("OZ1", "DANTE_TIR", "transposable_element", 10000, 11000, "+",
         ID="UA_L1_00000003", Name="Class_II/Subclass_1/TIR/EnSpm_CACTA",
         classification="Class_II/Subclass_1/TIR/EnSpm_CACTA",
         source_tier="1", source_tool="DANTE_TIR"),
    feat("OZ1", "DANTE_LINE", "transposable_element", 12000, 13000, "+",
         ID="UA_L1_00000004", Name="Class_I/LINE", classification="Class_I/LINE",
         source_tier="1", source_tool="DANTE_LINE"),
    feat("OZ1", "DANTE", "transposable_element", 14000, 14500, "+",
         ID="UA_L1_00000005", Name="Class_I/LTR/Ty1_copia/Ale",
         classification="Class_I/LTR/Ty1_copia/Ale",
         source_tier="2", source_tool="DANTE"),
    feat("OZ1", "TideCluster_default", "repeat_region", 20000, 30000, ".",
         ID="UA_L1_00000006", Name="TRC_1", classification="Satellite/TideCluster/TRC_1",
         source_tier="3", source_tool="TideCluster_default"),
    feat("OZ1", "TideCluster_default", "repeat_region", 31000, 40000, ".",
         ID="UA_L1_00000007", Name="TRC_2", classification="rDNA/45S_rDNA",
         source_tier="3", source_tool="TideCluster_default"),
    feat("OZ1", "TideCluster_default", "repeat_region", 41000, 50000, ".",
         ID="UA_L1_00000008", Name="TRC_13", classification="Satellite/TideCluster/TRC_13",
         source_tier="3", source_tool="TideCluster_default",
         TE_origin="Class_I/LTR/Ty1_copia/Ale", TE_origin_structure="tandem_LTR_RT"),
    feat("OZ1", "TideCluster_short", "repeat_region", 51000, 52000, ".",
         ID="UA_L1_00000009", Name="TRC_5", classification="Satellite/TideCluster/TRC_5",
         source_tier="3", source_tool="TideCluster_short"),
    feat("OZ1", "TideCluster_RM", "repeat_region", 53000, 54000, "+",
         ID="UA_L1_00000010", Name="TRC_1", classification="Satellite/TideCluster/TRC_1",
         source_tier="4", source_tool="TideCluster_RM"),
    feat("OZ1", "RepeatMasker", "transposable_element", 60000, 61000, ".",
         ID="UA_L1_00000011",
         Name="Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tat/Retand",
         classification="Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tat/Retand",
         source_tier="5", source_tool="RepeatMasker"),
    feat("OZ1", "RepeatMasker", "repeat_region", 62000, 62100, ".",
         ID="UA_L1_00000012", Name="rDNA/45S_rDNA/18S", classification="rDNA/45S_rDNA/18S",
         source_tier="5", source_tool="RepeatMasker"),
    feat("OZ1", "RepeatMasker", "repeat_region", 63000, 63030, ".",
         ID="UA_L1_00000013", Name="Simple_repeat(AT)n", classification="Simple_repeat(AT)n",
         source_tier="5", source_tool="RepeatMasker"),
    feat("OZ1", "RepeatMasker", "repeat_region", 64000, 64050, ".",
         ID="UA_L1_00000014", Name="Low_complexity", classification="Low_complexity",
         source_tier="5", source_tool="RepeatMasker"),
    feat("OZ1", "TideHunter", "repeat_region", 70000, 71000, ".",
         ID="UA_L1_00000015", Name="Satellite/Unknown", classification="Satellite/Unknown",
         source_tier="6", source_tool="TideHunter"),
    # Level-2 nested simple repeat inside the TRC_1 satellite above
    feat("OZ1", "RepeatMasker", "repeat_region", 25000, 25030, ".",
         ID="UA_L2_00000001", Name="Simple_repeat(ACC)n", classification="Simple_repeat(ACC)n",
         source_tier="5", source_tool="RepeatMasker", Parent="UA_L1_00000006"),
    # Tandem LTR-RT: a container (Level 1; structure/copy_number, no element_type)
    # plus a member copy (Level 2, Parent=container).  -- indices 16, 17 --
    feat("OZ1", "DANTE_LTR", "transposable_element", 80000, 110000, "-",
         ID="UA_L1_00000016", Name="Class_I/LTR/Ty3_gypsy/chromovirus/Tekay",
         classification="Class_I/LTR/Ty3_gypsy/chromovirus/Tekay",
         source_tier="1", source_tool="DANTE_LTR", structure="LTR_RT_TR", copy_number="4"),
    feat("OZ1", "DANTE_LTR", "transposable_element", 80000, 90000, "-",
         ID="UA_L2_00000002", Name="Class_I/LTR/Ty3_gypsy/chromovirus/Tekay",
         classification="Class_I/LTR/Ty3_gypsy/chromovirus/Tekay",
         source_tier="1", source_tool="DANTE_LTR", element_type="complete",
         in_structure="LTR_RT_TR", member_of="UA_L1_00000016",
         Parent="UA_L1_00000016"),
    # Tier-4 (RM-on-TideCluster) rDNA array: keeps bare Name=TRC_<n>, classified
    # rDNA from TideCluster's authoritative <prefix>_rdna.tsv (not the GFF3
    # rDNA_type attribute, which tier 4 never carries). Guards that the contract
    # keeps accepting tier-4 rDNA — index 18.  -- the make_unified rDNA fix --
    feat("OZ1", "TideCluster_RM", "repeat_region", 120000, 121000, "+",
         ID="UA_L1_00000017", Name="TRC_2", classification="rDNA/45S_rDNA",
         source_tier="4", source_tool="TideCluster_RM"),
]


def render(features):
    lines = ["##gff-version 3"]
    for f in features:
        attrs = ";".join(f"{k}={v}" for k, v in f["attrs"].items()) + ";"
        lines.append("\t".join(f["cols"] + [attrs]))
    return "\n".join(lines) + "\n"


def violations_for(features):
    with tempfile.NamedTemporaryFile("w", suffix=".gff3", delete=False) as fh:
        fh.write(render(features))
        path = fh.name
    try:
        return validate(path)
    finally:
        os.unlink(path)


def drift(mutate):
    """Deep-copy GOOD, apply mutate(features), return violations."""
    f = copy.deepcopy(GOOD)
    mutate(f)
    return violations_for(f)


def _set(idx, key, val):
    def m(f):
        f[idx]["attrs"][key] = val
    return m


def _del(idx, key):
    def m(f):
        f[idx]["attrs"].pop(key, None)
    return m


# (label, mutate, substring expected in at least one violation)
BAD_CASES = [
    ("satellite Name not TRC_ (the regression)",
     _set(5, "Name", "Satellite/TideCluster/TRC_1"), "must be bare TRC_"),
    ("rDNA satellite Name not TRC_",
     _set(6, "Name", "rDNA_45S"), "must be bare TRC_"),
    ("bad source_tier value",
     _set(0, "source_tier", "9"), "source_tier"),
    ("tier/tool mismatch (DANTE_LTR tier 2)",
     _set(0, "source_tier", "2"), "expected for DANTE_LTR"),
    ("unknown source_tool",
     lambda f: (f[10]["attrs"].__setitem__("source_tool", "Bogus"),
                f[10]["cols"].__setitem__(1, "Bogus")), "not in"),
    ("missing classification",
     _del(0, "classification"), "missing required attribute 'classification'"),
    ("element_type on RepeatMasker",
     _set(10, "element_type", "complete"), "element_type only allowed"),
    ("DANTE_LTR missing element_type",
     _del(0, "element_type"), "missing element_type"),
    ("type/classification mismatch",
     lambda f: f[5]["cols"].__setitem__(2, "transposable_element"), "implies 'repeat_region'"),
    ("TE_origin on RepeatMasker",
     _set(10, "TE_origin", "Class_I/LTR/Ty1_copia/Ale"), "TE_origin only allowed"),
    ("Level-2 without Parent",
     _del(15, "Parent"), "has no Parent"),
    ("dangling Parent",
     _set(15, "Parent", "UA_L1_99999999"), "non-existent ID"),
    ("Level-1 with a Parent",
     _set(0, "Parent", "UA_L1_00000006"), "must not have a Parent"),
    ("malformed ID",
     _set(0, "ID", "L1_1"), "does not match"),
    ("bad structure value",
     _set(16, "structure", "FOO"), "structure 'FOO' not in"),
    ("structure on non-DANTE_LTR",
     _set(10, "structure", "LTR_RT_TR"), "structure only allowed"),
    ("copy_number without structure",
     _set(0, "copy_number", "3"), "copy_number requires"),
    ("tandem member missing element_type",
     _del(17, "element_type"), "missing element_type"),
    ("in_structure on the Level-1 container",
     _set(16, "in_structure", "LTR_RT_TR"), "in_structure only allowed on Level-2"),
    ("in_structure on a non-DANTE_LTR Level-2 feature",
     _set(15, "in_structure", "LTR_RT_TR"), "in_structure only allowed on Level-2"),
    ("bad in_structure value",
     _set(17, "in_structure", "FOO"), "in_structure 'FOO' not in"),
    ("member_of without in_structure",
     lambda f: (f[17]["attrs"].pop("in_structure", None),
                f[17]["attrs"].__setitem__("member_of", "UA_L1_00000016")),
     "member_of requires in_structure"),
    ("member_of not equal to Parent",
     _set(17, "member_of", "UA_L1_00000006"), "must equal Parent"),
    ("TE_origin_structure without TE_origin",
     _set(5, "TE_origin_structure", "tandem_LTR_RT"), "requires TE_origin"),
    ("bad TE_origin_structure value",
     _set(7, "TE_origin_structure", "FOO"), "not in"),
]


def main():
    failures = []

    good_v = violations_for(GOOD)
    if good_v:
        failures.append("conforming sample REJECTED:\n    " + "\n    ".join(good_v))
    else:
        print(f"OK: conforming sample ({len(GOOD)} features, all kinds) passes")

    for label, mutate, expect in BAD_CASES:
        v = drift(mutate)
        if not v:
            failures.append(f"drift NOT caught: {label}")
        elif not any(expect in msg for msg in v):
            failures.append(f"drift '{label}': expected substring {expect!r} not in {v}")
        else:
            print(f"OK: caught drift — {label}")

    if failures:
        print("\nFAIL:", file=sys.stderr)
        for f in failures:
            print("  - " + f, file=sys.stderr)
        return 1
    print(f"\nAll {1 + len(BAD_CASES)} unified-GFF3 contract checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

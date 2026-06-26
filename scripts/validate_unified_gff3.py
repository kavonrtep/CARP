#!/usr/bin/env python3
"""Validate Repeat_Annotation_Unified.gff3 against the output contract.

A schema/format-drift guard for the pipeline's headline output. The contract is
specified in docs/unified_annotation_gff3_spec.md; this module is the executable
mirror. Pure standard library (no rtracklayer / pandas) so it runs anywhere,
including the minimal CI env and as a fast pre-commit / release check.

Usage:
    validate_unified_gff3.py Repeat_Annotation_Unified.gff3
    validate_unified_gff3.py --quiet out/Repeat_Annotation_Unified.gff3

Exit status 0 = conforms, 1 = one or more violations (printed to stderr), 2 =
usage / unreadable file. Importable: ``validate(path) -> list[str]`` returns the
(possibly empty) list of human-readable violation strings.
"""
import re
import sys

# ── The contract (keep in lock-step with docs/unified_annotation_gff3_spec.md) ──

ALLOWED_TYPES = {"transposable_element", "repeat_region"}

# source_tool (GFF3 col 2) -> its fixed source_tier. The pairing is 1:1.
TOOL_TIER = {
    "DANTE_LTR": 1, "DANTE_TIR": 1, "DANTE_LINE": 1,
    "DANTE": 2,
    "TideCluster_default": 3, "TideCluster_short": 3,
    "TideCluster_RM": 4,
    "RepeatMasker": 5,
    "TideHunter": 6,
}
ALLOWED_SOURCES = set(TOOL_TIER)

# Satellites from these tools MUST keep Name = bare TRC_<n> (downstream apps and
# split_gff_by_name.R --name-prefix TRC_ key on it). This is the invariant whose
# breakage silently emptied the per-family BigWig outputs.
TRC_NAME_TOOLS = {"TideCluster_default", "TideCluster_short", "TideCluster_RM"}
TRC_NAME_RE = re.compile(r"^TRC_\d+$")

# A feature is a repeat_region iff its classification starts with one of these.
REPEAT_REGION_PREFIX_RE = re.compile(r"^(Satellite|Simple_repeat|Low_complexity|rDNA|Unknown)")

ID_RE = re.compile(r"^UA_L([12])_\d{8}$")
PARENT_RE = re.compile(r"^UA_L1_\d{8}$")
STRANDS = {"+", "-", "."}
REQUIRED_ATTRS = ("ID", "Name", "classification", "source_tier", "source_tool")
ELEMENT_TYPE_VALUES = {"complete", "partial"}
ELEMENT_TYPE_TOOL = "DANTE_LTR"          # element_type present iff this tool
TE_ORIGIN_TOOLS = {"TideCluster_default", "TideCluster_short"}  # only TE-derived satellites carry it


def _parse_attrs(col9):
    out = {}
    for field in col9.rstrip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        k, _, v = field.partition("=")
        out[k] = v
    return out


def validate(path):
    """Return a list of violation strings (empty == conforms)."""
    v = []
    try:
        with open(path) as fh:
            lines = fh.readlines()
    except OSError as e:
        return [f"cannot read {path}: {e}"]

    seen_ids = {}
    parents = []          # (lineno, feature_id, parent_id)
    n_feat = 0

    for lineno, raw in enumerate(lines, 1):
        if raw.startswith("#") or not raw.strip():
            continue
        cols = raw.rstrip("\n").split("\t")
        if len(cols) != 9:
            v.append(f"L{lineno}: expected 9 tab columns, got {len(cols)}")
            continue
        seqid, source, ftype, start, end, score, strand, phase, col9 = cols
        n_feat += 1
        a = _parse_attrs(col9)

        # column-level
        if not seqid:
            v.append(f"L{lineno}: empty seqid")
        if source not in ALLOWED_SOURCES:
            v.append(f"L{lineno}: source(col2) '{source}' not in {sorted(ALLOWED_SOURCES)}")
        if ftype not in ALLOWED_TYPES:
            v.append(f"L{lineno}: type(col3) '{ftype}' not in {sorted(ALLOWED_TYPES)}")
        try:
            s, e = int(start), int(end)
            if s < 1 or e < s:
                v.append(f"L{lineno}: bad coordinates start={start} end={end}")
        except ValueError:
            v.append(f"L{lineno}: non-integer start/end ({start},{end})")
        if strand not in STRANDS:
            v.append(f"L{lineno}: strand(col7) '{strand}' not in {sorted(STRANDS)}")

        # required attributes
        for k in REQUIRED_ATTRS:
            if k not in a:
                v.append(f"L{lineno}: missing required attribute '{k}'")
        if any(k not in a for k in REQUIRED_ATTRS):
            continue  # downstream checks need these

        fid, name, cls, tier, tool = (a["ID"], a["Name"], a["classification"],
                                      a["source_tier"], a["source_tool"])

        # ID format + uniqueness + level
        m = ID_RE.match(fid)
        if not m:
            v.append(f"L{lineno}: ID '{fid}' does not match UA_L[12]_<8 digits>")
            level = None
        else:
            level = m.group(1)
            if fid in seen_ids:
                v.append(f"L{lineno}: duplicate ID '{fid}' (also L{seen_ids[fid]})")
            seen_ids[fid] = lineno

        # source_tool == col2, and tier matches the tool
        if tool != source:
            v.append(f"L{lineno}: source_tool '{tool}' != col2 source '{source}'")
        if tool in TOOL_TIER and tier != str(TOOL_TIER[tool]):
            v.append(f"L{lineno}: source_tier '{tier}' != {TOOL_TIER[tool]} expected for {tool}")

        # type <-> classification consistency
        is_rr = bool(REPEAT_REGION_PREFIX_RE.match(cls))
        want = "repeat_region" if is_rr else "transposable_element"
        if ftype != want:
            v.append(f"L{lineno}: type '{ftype}' but classification '{cls}' implies '{want}'")

        # Name invariant for TideCluster satellites
        if tool in TRC_NAME_TOOLS and not TRC_NAME_RE.match(name):
            v.append(f"L{lineno}: {tool} Name '{name}' must be bare TRC_<n> "
                     f"(downstream apps / split_gff_by_name.R depend on it)")

        # element_type present iff DANTE_LTR; value constrained
        if "element_type" in a:
            if tool != ELEMENT_TYPE_TOOL:
                v.append(f"L{lineno}: element_type only allowed on {ELEMENT_TYPE_TOOL}, not {tool}")
            elif a["element_type"] not in ELEMENT_TYPE_VALUES:
                v.append(f"L{lineno}: element_type '{a['element_type']}' not in {sorted(ELEMENT_TYPE_VALUES)}")
        elif tool == ELEMENT_TYPE_TOOL:
            v.append(f"L{lineno}: {ELEMENT_TYPE_TOOL} feature missing element_type")

        # TE_origin only on TE-derived clustering satellites; value is a class path
        if "TE_origin" in a:
            if tool not in TE_ORIGIN_TOOLS:
                v.append(f"L{lineno}: TE_origin only allowed on {sorted(TE_ORIGIN_TOOLS)}, not {tool}")
            if "/" not in a["TE_origin"]:
                v.append(f"L{lineno}: TE_origin '{a['TE_origin']}' is not a slash classification path")

        # Parent <-> level
        if level == "2" and "Parent" not in a:
            v.append(f"L{lineno}: Level-2 feature {fid} has no Parent")
        if level == "1" and "Parent" in a:
            v.append(f"L{lineno}: Level-1 feature {fid} must not have a Parent")
        if "Parent" in a:
            if not PARENT_RE.match(a["Parent"]):
                v.append(f"L{lineno}: Parent '{a['Parent']}' does not match UA_L1_<8 digits>")
            else:
                parents.append((lineno, fid, a["Parent"]))

    # cross-feature: every Parent resolves to a Level-1 ID present in the file
    for lineno, fid, pid in parents:
        if pid not in seen_ids:
            v.append(f"L{lineno}: feature {fid} Parent '{pid}' references a non-existent ID")

    if n_feat == 0:
        v.append("file contains no feature rows")
    return v


def main(argv):
    quiet = False
    args = [x for x in argv[1:] if x != "--quiet"]
    if "--quiet" in argv[1:]:
        quiet = True
    if len(args) != 1:
        print("usage: validate_unified_gff3.py [--quiet] <Repeat_Annotation_Unified.gff3>",
              file=sys.stderr)
        return 2
    violations = validate(args[0])
    if violations:
        print(f"FAIL: {len(violations)} unified-GFF3 contract violation(s) in {args[0]}:",
              file=sys.stderr)
        for msg in violations:
            print("  - " + msg, file=sys.stderr)
        return 1
    if not quiet:
        print(f"OK: {args[0]} conforms to the unified-annotation GFF3 contract")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

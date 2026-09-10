#!/usr/bin/env python3
"""Guard: `[[name]]` lookups on a name->value map must test MEMBERSHIP.

In R, a named vector or list subscripted with a name it does not hold **aborts**:

    > period_map <- setNames(120L, "TRC_1"); period_map[["TRC_9"]]
    Error in period_map[["TRC_9"]] : subscript out of bounds

It does not return NULL. So neither of these guards works:

    if (length(m) > 0) m[[n]] else NULL          # size, not membership
    if (length(m) > 0 && !is.null(m[[n]])) ...   # `[[` throws before is.null()

Both shapes shipped in `make_unified_annotation.R` and killed a CARP 1.7.1 run
on a 1.04 Gbp genome at `make_unified_annotation` -- 70 % through, after four
hours of upstream work -- with `Error in period_map[[name]] : subscript out of
bounds` inside `identify_te_derived_trcs`. The intent was documented on the line
above the bug ("skipped only when no period is available for this TRC"), and it
is not an edge case: `read_trc_periods()` drops every TRC whose monomer cannot
be parsed, so a clustering TRC absent from the map is routine. That run had 68
periods in `trc_table.tsv` and a TE_origin candidate that was not one of them.

The correct idiom is already used elsewhere in the same file, for `rdna_map`:

    hit <- trc %in% names(rdna_map)
    cls[hit] <- unname(rdna_map[trc[hit]])

This test reads the real script as text rather than copying its functions, so it
cannot drift out of sync with the source the way a verbatim copy can.

Exit 0 on success, 1 with an actionable message on failure.
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
TARGET = ROOT / "scripts" / "make_unified_annotation.R"

# Any identifier that looks like a name->value map. Deliberately broad: a new
# map added later is covered without editing this test.
MAP_SUBSCRIPT = re.compile(r"\b(\w*(?:_map|map|per_map|periods?))\[\[")
# A membership test on the same line, in any of the accepted spellings.
MEMBERSHIP = re.compile(r"%in%\s*names\(|\bhasName\(|\bexists\(|match\([^)]*names\(")

FAILURES = []


def main():
    if not TARGET.exists():
        print(f"FAIL: {TARGET} not found")
        return 1

    lines = TARGET.read_text().split("\n")
    checked = 0
    for i, line in enumerate(lines, start=1):
        code = line.split("#", 1)[0]          # ignore comments (they cite the idiom)
        if not MAP_SUBSCRIPT.search(code):
            continue
        checked += 1
        if MEMBERSHIP.search(code):
            continue
        # A literal or a loop variable known to come from names(map) is fine only
        # if guarded; anything else is the bug shape.
        FAILURES.append((i, line.strip()))

    print(f"  scanned {TARGET.relative_to(ROOT)}: {checked} map [[...]] lookup(s)")
    if FAILURES:
        print("\nFAIL: map [[name]] lookup without a membership test:\n")
        for ln, txt in FAILURES:
            print(f"  {TARGET.relative_to(ROOT)}:{ln}:  {txt}")
        print(
            "\n  In R, m[[n]] ABORTS when n is not a name of m -- it does not\n"
            "  return NULL, and `!is.null(m[[n]])` cannot guard it because `[[`\n"
            "  throws first. Use the idiom already in this file:\n\n"
            "      if (n %in% names(m)) m[[n]] else <default>\n")
        return 1

    if checked == 0:
        print("FAIL: found no map [[...]] lookups at all — has the guard's\n"
              "      pattern gone stale against a refactor?")
        return 1

    print("OK: every map [[name]] lookup tests membership first.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

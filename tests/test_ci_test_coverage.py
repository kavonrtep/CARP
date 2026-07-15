#!/usr/bin/env python3
"""Guard: every tests/test_*.{py,R,sh} must be executed by a CI workflow, or be
explicitly listed in CI_EXEMPT with a reason.

Why this exists: a regression test that lives on disk but is never wired into a
workflow gives false confidence — it protects nothing. Five scaling-branch
regression tests sat unwired this way, and the 1.1.0 `make_unified` provenance
crash reached CI in a code path that had no running test at all. This guard
turns "silently orphaned" into "CI red until you either run the test or justify
skipping it".

Add a new test either by referencing it in a `.github/workflows/*.yml` step, or
by adding it to CI_EXEMPT below with a one-line reason (heavy env / benchmark /
integration-covered). Exemptions that become stale (the test now runs, or the
file was deleted) also fail this guard, so the list can't rot.
"""
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
TESTS = ROOT / "tests"
WORKFLOWS = ROOT / ".github" / "workflows"

# Tests intentionally NOT run in the lightweight unit job. Each MUST carry a
# reason. These need a heavy toolchain (bioconductor R stack / mmseqs / blast)
# or are manual benchmarks; the underlying code is exercised end-to-end by the
# small+medium fixture runs in pipeline.yml. Promote one out of this list by
# wiring it into a workflow (e.g. a dedicated bioconductor-enabled job).
CI_EXEMPT = {
    "test_density_utils.R":
        "needs bioconductor rtracklayer (not in carp-unit); density is run by the fixture jobs",
    "test_gff_cleanup.R":
        "needs bioconductor rtracklayer (not in carp-unit); gff_cleanup is run by the fixture jobs",
    "test_resolve_tier1_overlaps.R":
        "needs bioconductor GenomicRanges (not in carp-unit); the resolver is run by the fixture jobs",
    "test_reduce_dimer_parallel.sh":
        "needs mmseqs2; dimer reduction is run by the fixture jobs",
    "test_reduce_library_parity.sh":
        "needs mmseqs2 + blast; library reduction is run by the fixture jobs",
    "test_reduce_library_memprofile.sh":
        "manual memory-profiling benchmark, not a correctness gate",
}


def workflow_text() -> str:
    return "\n".join(p.read_text() for p in sorted(WORKFLOWS.glob("*.yml")))


def discover_tests() -> list[str]:
    names = []
    for pat in ("test_*.py", "test_*.R", "test_*.sh"):
        names += [p.name for p in TESTS.glob(pat)]
    return sorted(names)


def main() -> int:
    wf = workflow_text()
    tests = discover_tests()

    orphans, stale_exempt = [], []
    for t in tests:
        runs = t in wf                     # referenced by some workflow step
        exempt = t in CI_EXEMPT
        if runs and exempt:
            stale_exempt.append(t)
        elif not runs and not exempt:
            orphans.append(t)

    missing_exempt = [t for t in CI_EXEMPT if not (TESTS / t).exists()]

    ok = True
    if orphans:
        ok = False
        print("FAIL: these test files are never run by any CI workflow and are not exempted:")
        for t in orphans:
            print(f"  - {t}   -> add it to a workflow 'run:' step, or to CI_EXEMPT with a reason")
    if stale_exempt:
        ok = False
        print("FAIL: these tests are BOTH run by CI and listed in CI_EXEMPT — drop the stale exemption:")
        for t in stale_exempt:
            print(f"  - {t}")
    if missing_exempt:
        ok = False
        print("FAIL: CI_EXEMPT names test files that no longer exist — remove them:")
        for t in missing_exempt:
            print(f"  - {t}")

    if ok:
        n_run = sum(1 for t in tests if t in wf)
        print(f"OK: {n_run}/{len(tests)} test files are run by CI; "
              f"{len(CI_EXEMPT)} explicitly exempted (with reasons).")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())

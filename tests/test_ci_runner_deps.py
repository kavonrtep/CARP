#!/usr/bin/env python3
"""Guard: every test is wired into a CI job whose environment actually PROVIDES
the packages that test imports.

Why this exists: a test whose file lives on disk and is referenced by a workflow
still fails at *runtime* if the job's runner lacks its dependency — e.g. an R test
that does `library(GenomicRanges)` wired into the lightweight `carp-unit` env
(python/r-base/r-yaml only). That has bitten this repo more than once and only
surfaces as a red CI job, never at commit time. `test_ci_test_coverage.py` proves
a test is *run somewhere*; this proves it is run *where its dependencies exist*.

How it works — no drift, derived from the YAML itself:
- For each `.github/workflows/*.yml` job we read the conda/mamba `create-args`
  package list (plus any `pip/mamba/micromamba install` or R `install.packages` /
  `BiocManager::install` in run-steps), and whether it restores the full pipeline
  conda envs (`.snakemake/conda`). A job that restores those is treated as
  providing everything.
- For each R/Python test that job runs, we extract the packages the test imports
  and require the job to provide each one.

Scope (documented so a gap is a known gap, not a surprise):
- R: direct library()/require()/requireNamespace() in the test file. Deps pulled
  in transitively via source()'d scripts are not followed.
- Python: imports that are neither stdlib nor a repo-local module (scripts/*.py or
  tests/*.py — tests put scripts/ on sys.path). Local modules are always present.
- Shell (.sh) tests are skipped: they locate their own toolchain and skip
  gracefully when it is absent (e.g. test_unified_multibatch_determinism.sh).

What this CANNOT catch (only running the test does): a package that a test never
`library()`s but that its ops lazy-load at RUNTIME — e.g. GenomicRanges pulling
XVector via range()/revmap on CompressedIntegerList, when bioconda's
bioconductor-genomicranges did not install XVector. That is a conda-recipe
incompleteness, not a wiring mistake, and is invisible to static analysis; the
job must list such deps explicitly and the test run is what verifies the env.

To satisfy the guard: put the test in a job whose env installs its deps, or add
the dep to that job's create-args. If a package is provided by an install command
this parser doesn't recognise, extend it — don't weaken the check.
"""
import ast
import glob
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
TESTS = ROOT / "tests"
SCRIPTS = ROOT / "scripts"
WORKFLOWS = ROOT / ".github" / "workflows"

# Base + recommended R packages that ship with `r-base` (importable without any
# extra conda package). Stable across R versions.
BASE_R = {
    "base", "methods", "utils", "stats", "stats4", "graphics", "grDevices",
    "datasets", "grid", "splines", "tools", "tcltk", "parallel", "compiler",
    "MASS", "Matrix",
}
# Core Bioconductor infrastructure that IS a hard dependency of essentially every
# bioconductor-* conda package (verified: bioconductor-genomicranges pulls exactly
# these). Deliberately does NOT include XVector / GenomeInfoDb / Biobase: those are
# NOT reliably co-installed (bioconductor-genomicranges does not pull XVector on
# bioconda), so a test that library()'s one of them needs its package listed
# explicitly and the guard must flag its absence rather than assume it.
BIOC_CORE = {"BiocGenerics", "S4Vectors", "IRanges"}
# Python import name -> conda package name, where they differ.
PY_IMPORT_TO_CONDA = {
    "yaml": "pyyaml", "Bio": "biopython", "sklearn": "scikit-learn",
    "cv2": "opencv", "PIL": "pillow", "bs4": "beautifulsoup4",
}

# ── workflow parsing ──────────────────────────────────────────────────────────

# Repo directories that end up on sys.path when tests run (tests add scripts/ and
# run from the repo root, so root-level modules like version.py import too).
LOCAL_MODULE_DIRS = (ROOT, SCRIPTS, TESTS)

def _strip_version(tok: str) -> str:
    return re.split(r"[=<>~!]", tok, maxsplit=1)[0].strip()

def _is_local_module(mod: str) -> bool:
    return any((d / f"{mod}.py").exists() or (d / mod / "__init__.py").exists()
               for d in LOCAL_MODULE_DIRS)

def load_yaml(path):
    import yaml  # provided by pyyaml in the unit env
    return yaml.safe_load(path.read_text())

def parse_jobs():
    """Yield dicts: {label, tests:[...], pkgs:set, provides_all:bool}."""
    jobs = []
    for wf in sorted(WORKFLOWS.glob("*.yml")):
        doc = load_yaml(wf)
        for jname, job in (doc.get("jobs") or {}).items():
            pkgs, tests, provides_all = set(), [], False
            for st in job.get("steps", []) or []:
                with_ = st.get("with", {}) or {}
                ca = with_.get("create-args")
                if ca:
                    pkgs.update(_strip_version(t) for t in str(ca).split())
                uses = str(st.get("uses", ""))
                if "actions/cache/restore" in uses and ".snakemake/conda" in str(with_.get("path", "")):
                    provides_all = True
                run = st.get("run", "") or ""
                # install commands in run-steps
                for m in re.finditer(r"(?:pip|mamba|micromamba)\s+install\s+([^\n|&]+)", run):
                    pkgs.update(_strip_version(t) for t in m.group(1).split()
                                if not t.startswith("-"))
                for m in re.finditer(r"(?:install\.packages|BiocManager::install)\(\s*['\"]([^'\"]+)", run):
                    pkgs.add(m.group(1))
                for m in re.findall(r"tests/test_[A-Za-z0-9_]+\.(?:py|R|sh)", run):
                    tests.append(m)
            if tests:
                jobs.append({"label": f"{wf.name}::{jname}", "tests": sorted(set(tests)),
                             "pkgs": pkgs, "provides_all": provides_all})
    return jobs

# ── dependency extraction ─────────────────────────────────────────────────────

def r_libs(path: pathlib.Path) -> set:
    txt = path.read_text()
    libs = set()
    for m in re.finditer(r"\b(?:library|require|requireNamespace)\(\s*['\"]?([A-Za-z][A-Za-z0-9._]*)", txt):
        libs.add(m.group(1))
    return libs

def _stdlib() -> set:
    return set(getattr(sys, "stdlib_module_names", set())) | {
        "__future__",
    }

def py_third_party(path: pathlib.Path) -> set:
    try:
        tree = ast.parse(path.read_text())
    except SyntaxError:
        return set()
    mods = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            mods.update(a.name.split(".")[0] for a in n.names)
        elif isinstance(n, ast.ImportFrom) and n.level == 0 and n.module:
            mods.add(n.module.split(".")[0])
    std = _stdlib()
    out = set()
    for m in mods:
        if m in std or m == "tests":
            continue
        if _is_local_module(m):
            continue  # repo-local module (root / scripts / tests on sys.path)
        out.add(m)
    return out

# ── satisfaction checks ───────────────────────────────────────────────────────

def r_lib_ok(lib: str, job: dict) -> bool:
    if job["provides_all"] or lib in BASE_R:
        return True
    lc = lib.lower()
    if f"r-{lc}" in job["pkgs"] or f"bioconductor-{lc}" in job["pkgs"]:
        return True
    if lib in BIOC_CORE and any(p.startswith("bioconductor-") for p in job["pkgs"]):
        return True
    return False

def py_mod_ok(mod: str, job: dict) -> bool:
    if job["provides_all"]:
        return True
    conda = PY_IMPORT_TO_CONDA.get(mod, mod)
    return conda in job["pkgs"] or mod in job["pkgs"]

# ── main ──────────────────────────────────────────────────────────────────────

def main() -> int:
    jobs = parse_jobs()
    violations = []
    for job in jobs:
        for t in job["tests"]:
            p = ROOT / t
            if not p.exists() or t.endswith(".sh"):
                continue
            if t.endswith(".R"):
                for lib in sorted(r_libs(p)):
                    if not r_lib_ok(lib, job):
                        violations.append((job["label"], t, f"library({lib})"))
            elif t.endswith(".py"):
                for mod in sorted(py_third_party(p)):
                    if not py_mod_ok(mod, job):
                        violations.append((job["label"], t, f"import {mod}"))

    if violations:
        print("FAIL: tests wired into a job whose env does not provide their deps:")
        for label, t, dep in violations:
            print(f"  - {t}  needs  {dep}  but  {label}  does not install it")
        print("\nFix: move the test to a job whose env has the dependency, or add "
              "the package to that job's create-args.")
        return 1

    n = sum(len(j["tests"]) for j in jobs)
    print(f"OK: every wired R/Python test's deps are provided by its CI job "
          f"({n} test-in-job placements across {len(jobs)} jobs checked).")
    return 0

if __name__ == "__main__":
    sys.exit(main())

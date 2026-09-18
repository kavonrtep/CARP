#!/usr/bin/env python3
"""Guard: `dante_ltr_mode` defaults to lineage, and the library conflict policy
follows it.

Why this exists. `dante_ltr_to_library --annotation_conflict nested` changes
which clusters reach the RepeatMasker library, and upstream is explicit that
**nested is not a superset of strict** — a cluster mixing two sibling lineages
with their shared parent is kept by `strict` and dropped by `nested`. So the
policy must stay welded to `dante_ltr_mode`: a default (lineage) run has to get
`strict`, byte-identical to every release before the knob existed. A stray edit
that flips the default, or drops the `params` reference from the shell body,
would silently change the library on every genome.

The DAG is not available here — the lightweight unit job installs no snakemake —
so instead of dry-running the rule this evaluates the ACTUAL expression from the
Snakefile under both config values. That makes it a behavioural check rather
than a text match: rewriting the expression in any equivalent form still passes,
while deleting the wiring fails.

Coupled to the Snakefile's structure: it locates the `make_library_of_ltrs` rule
and the `conflict_policy` assignment inside it. If either is renamed, update the
names below.
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SNAKEFILE = (ROOT / "Snakefile").read_text()

RULE = "make_library_of_ltrs"
PARAM = "conflict_policy"
failures = []


def rule_body(name):
    """Text of `rule <name>:` up to the next top-level rule/checkpoint."""
    m = re.search(rf"^rule {name}:$", SNAKEFILE, re.M)
    if not m:
        return None
    rest = SNAKEFILE[m.end():]
    nxt = re.search(r"^(rule|checkpoint) [a-z_0-9]+:$", rest, re.M)
    return rest[: nxt.start()] if nxt else rest


def check(label, got, want):
    if got != want:
        failures.append(f"{label}: got {got!r}, want {want!r}")


# ── 1. the pipeline default is lineage ────────────────────────────────────
m = re.search(r"""config\[["']dante_ltr_mode["']\]\s*=\s*(["'][a-z]+["'])""", SNAKEFILE)
if not m:
    failures.append("no default assignment for config['dante_ltr_mode'] in the Snakefile")
else:
    check("default dante_ltr_mode", eval(m.group(1)), "lineage")

# ── 2. exactly two modes are accepted ─────────────────────────────────────
m = re.search(r"""config\[["']dante_ltr_mode["']\]\s+not\s+in\s+(\([^)]*\))""", SNAKEFILE)
if not m:
    failures.append("no membership validation for config['dante_ltr_mode']")
else:
    check("accepted modes", tuple(sorted(eval(m.group(1)))), ("core", "lineage"))

# ── 3. the policy expression maps mode -> flag value ──────────────────────
body = rule_body(RULE)
if body is None:
    failures.append(f"rule {RULE} not found")
else:
    m = re.search(rf"^\s*{PARAM}\s*=\s*(.+?)$", body, re.M)
    if not m:
        failures.append(f"no {PARAM} param in rule {RULE}")
    else:
        expr = m.group(1).rstrip(",")
        for mode, want in (("lineage", "strict"), ("core", "nested")):
            try:
                got = eval(expr, {}, {"config": {"dante_ltr_mode": mode}})
            except Exception as exc:                      # noqa: BLE001
                failures.append(f"{PARAM} expression failed for {mode}: {exc}")
                continue
            check(f"{PARAM} under dante_ltr_mode={mode}", got, want)

    # ── 4. the shell body actually uses it ────────────────────────────────
    if f"--annotation_conflict {{params.{PARAM}}}" not in body:
        failures.append(
            f"rule {RULE} does not pass --annotation_conflict {{params.{PARAM}}} "
            "to dante_ltr_to_library — the policy would fall back to the tool's "
            "default and the knob would be inert")

if failures:
    print("FAIL: dante_ltr_mode wiring")
    for f in failures:
        print(f"  - {f}")
    sys.exit(1)
print("OK: dante_ltr_mode defaults to lineage; library conflict policy is "
      "strict under lineage and nested under core, and is passed to the tool.")

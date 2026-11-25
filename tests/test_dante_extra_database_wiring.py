#!/usr/bin/env python3
"""Guard: `dante_extra_database` is off by default and reaches DANTE when set.

Why this exists. With the option off, the DANTE command and the rule graph must
not change; with it on, DANTE must get the CHECKED copy, never the user's file
directly. No snakemake here (the unit job installs none), so this evaluates the
actual expressions from the Snakefile, as test_dante_ltr_mode_wiring.py does.
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SNAKEFILE = (ROOT / "Snakefile").read_text()
KEY = "dante_extra_database"
failures = []


def check(label, ok, detail=""):
    if not ok:
        failures.append(f"{label}: {detail}")


def rule_body(name):
    m = re.search(rf"^\s*rule {name}:$", SNAKEFILE, re.M)
    if not m:
        failures.append(f"rule {name} not found")
        return ""
    rest = SNAKEFILE[m.end():]
    nxt = re.search(r"^\s*(rule|checkpoint) [a-z_0-9]+:$", rest, re.M)
    return rest[: nxt.start()] if nxt else rest


def evaluate(pattern, text, value, names):
    m = re.search(pattern, text, re.M | re.S)
    if not m:
        failures.append(f"pattern not found: {pattern}")
        return None
    return eval(m.group(1), {}, {"config": {KEY: value}, **names})


# 1. default handling: run the real block on missing / None / value / wrong type
block = re.search(r'^if config\.get\("dante_extra_database"\) is None:.*?'
                  r'raise ValueError\("Invalid value for dante_extra_database[^\n]*\n',
                  SNAKEFILE, re.M | re.S)
if not block:
    failures.append("no default/validation block for dante_extra_database")
else:
    for given, want in (({}, ""), ({KEY: None}, ""), ({KEY: ""}, ""), ({KEY: "x.fa"}, "x.fa")):
        cfg = dict(given)
        exec(block.group(0), {"config": cfg})
        check(f"default for {given}", cfg[KEY] == want, cfg[KEY])
    try:
        exec(block.group(0), {"config": {KEY: ["a"]}})
        failures.append("a list value was accepted")
    except ValueError:
        pass

# 2. dante: checked copy as input and -e, nothing when off
names = {"DANTE_EXTRA_DB": "OUT/extra.fa", "DANTE_EXTRA_DB_SUMMARY": "OUT/extra.tsv"}
dante = rule_body("dante")
for value, want_in, want_flag in (("", [], ""), ("user.fa", "OUT/extra.fa", "-e OUT/extra.fa")):
    got = evaluate(r"^\s*extra_db=(.+?)$", dante, value, names)
    check(f"dante input, {KEY}={value!r}", got == want_in, got)
    expr = re.search(r'^\s*extra_flag=F"(.+?)" if (.+?) else ""$', dante, re.M)
    if not expr:
        failures.append("dante extra_flag not found")
    else:
        got = (expr.group(1).format(**names)
               if eval(expr.group(2), {}, {"config": {KEY: value}}) else "")
        check(f"dante flag, {KEY}={value!r}", got == want_flag, got)
cmd = re.search(r"^\s*dante -q .*$", dante, re.M)
check("dante command uses the flag", cmd and "{params.extra_flag}" in cmd.group(0))
check("dante reads the genome only with -q", cmd and "-q {input.fasta}" in cmd.group(0))
check("dante states the REXdb version", cmd and "-D {params.rexdb}" in cmd.group(0))

# 3. the checking rule exists only when the option is set, same REXdb as dante
m = re.search(r'^if config\["dante_extra_database"\]:\n    rule prepare_dante_extra_database:$',
              SNAKEFILE, re.M)
check("prepare rule is conditional", bool(m))
prep = rule_body("prepare_dante_extra_database")
for text, label in ((prep, "prepare"), (dante, "dante")):
    check(f"{label} uses DANTE_REXDB", "rexdb=DANTE_REXDB" in text)
check("prepare env is dante's", '"envs/dante.yaml"' in prep)

# 4. the report gets the summary only when the option is set
rep = rule_body("make_repeat_report")
for value, want in (("", []), ("user.fa", "OUT/extra.tsv")):
    got = evaluate(r"^\s*dante_extra_db = (.+?)$", rep, value, names)
    check(f"report input, {KEY}={value!r}", got == want, got)
check("report command uses the flag", "{params.extra_db_flag}" in rep.split("shell:")[-1])

# 5. provenance keeps the key (the report falls back to it)
check("provenance key", f'"{KEY}"' in (ROOT / "scripts" / "record_provenance.py").read_text())

if failures:
    print("FAIL: dante_extra_database wiring")
    for f in failures:
        print(f"  - {f}")
    sys.exit(1)
print("OK: dante_extra_database is off by default, checked before DANTE, and reported.")

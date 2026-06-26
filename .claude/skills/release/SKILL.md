---
name: release
description: Cut a CARP pipeline release the project's way (version bump + annotated tag). Use whenever asked to "create/cut/make a release", "bump the version", or "tag X" (e.g. "create new release 1.0.0rc1"). Encodes version.py as the single source of truth, the unprefixed-tag convention, the CI gates, and the host-push step so a release is one command instead of an investigation.
---

# Cutting a CARP release

`version.py` is the **single source of truth** for the version — read by
`run_pipeline.py` (`--version`/banner), the `Singularity` `Version` label,
`scripts/record_provenance.py`, the per-run HTML report, and the CI gates. A
release is: **bump `version.py`** + a git **tag whose name equals `version.py`
exactly**.

Tags are **unprefixed** PEP-440-lite — `MAJOR.MINOR.PATCH[aN|bN|rcN]`
(`1.0.0rc1`, `1.0.0`, `1.1.0a1`). **No leading `v`** (both `version.yml` and
`release.yml` trigger on `tags: ['[0-9]*']`).

## One command — does everything reversible
```bash
.claude/skills/release/cut-release.sh <VERSION>      # e.g. 1.0.0rc1
```
The helper: validates the version shape, refuses a dirty working tree or an
already-existing tag, checks the new version is **> current** (monotonic),
bumps **both** `__version__` and `__version_info__` in lock-step, runs the cheap
CI gates (version parse + the Python unit tests; R mirror test best-effort),
commits `Release <VERSION>`, and creates the **annotated** tag. It does **not**
push.

## Push — from the HOST, not the sandbox
The remote is SSH (`git@github.com:kavonrtep/...`) and the container has **no
`ssh`**, so the agent cannot push. Ask the user to run, in their host terminal
(NOT via `!`, which runs inside this sandbox):
```bash
cd <repo>
git push origin main
git push origin <VERSION>
```

## What the tag push triggers (`.github/workflows/release.yml`)
unit tests → build Singularity SIF → **run the medium fixture inside the SIF**
(the real release gate) → push SIF to GHCR
(`ghcr.io/kavonrtep/carp/sif:<VERSION>` and `:latest`) → create the GitHub
Release → **Zenodo snapshot (permanent DOI)**. The stages are gated in series:
if the in-container fixture fails, **nothing is published** — fix, delete the
tag, re-tag.

## Guards / gotchas (so you don't re-discover them)
- `version.yml` enforces three checks the helper pre-runs: string **parses**,
  pushed **tag == version.py**, version **monotonic** (`>=` previous main).
- A `PreToolUse` hook in `.claude/settings.json` blocks history rewrites.
  **Never** amend or move a commit/tag already on `origin` — the user pushes
  from the host between turns, so a "local-only" commit may already be remote.
- A tooling-only commit (docs, this skill) needs **no** version bump — the
  monotonic gate allows an equal version.
- The release SIF does not include `.claude/` (not in `Singularity` `%files`),
  so adding a skill never affects the built image.
- If the bump is for a **final** release (no modifier), `__version_info__`'s 4th
  element is `None` (e.g. `(1, 0, 0, None)`) — the helper handles this.

## Manual fallback (if the helper can't run)
```bash
# 1. edit version.py: __version__ = "X" ; __version_info__ = (maj,min,pat,"mod" or None)
python3 -c "import sys;sys.path.insert(0,'.');from version import __version__,parse_version;print(parse_version(__version__))"
python3 tests/test_classification.py && python3 tests/test_manifest.py   # cheap gates
git add version.py && git commit -m "Release X"
git tag -a X -m "CARP X"
# then user pushes main + tag from the host (see above)
```

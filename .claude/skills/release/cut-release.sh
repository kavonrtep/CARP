#!/usr/bin/env bash
# cut-release.sh <VERSION>
#
# Cut a CARP release the project's way: validate, bump version.py, run the cheap
# CI gates, commit, and create the annotated tag. Does NOT push (the sandbox has
# no ssh; the maintainer pushes from the host). See SKILL.md beside this file.
#
# VERSION is unprefixed PEP-440-lite: MAJOR.MINOR.PATCH[aN|bN|rcN]  e.g. 1.0.0rc1
set -euo pipefail

VER="${1:-}"
if [ -z "$VER" ]; then
  echo "usage: $(basename "$0") <VERSION>   e.g. 1.0.0rc1" >&2
  exit 2
fi

cd "$(git rev-parse --show-toplevel)"

# ── 1. Version shape (same grammar as version.py:parse_version) ───────────────
python3 - "$VER" <<'PY'
import sys, re
v = sys.argv[1]
if not re.fullmatch(r"(\d+)\.(\d+)\.(\d+)((?:a|b|rc)\d+)?", v):
    sys.exit(f"ERROR: {v!r} is not PEP-440-lite (MAJOR.MINOR.PATCH with optional aN/bN/rcN).")
PY

# ── 2. Clean working tree (tracked files) ─────────────────────────────────────
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "ERROR: working tree has uncommitted tracked changes — commit/stash them first:" >&2
  git status --short --untracked-files=no >&2
  exit 1
fi

# ── 3. Tag must not already exist ─────────────────────────────────────────────
if git rev-parse -q --verify "refs/tags/$VER" >/dev/null 2>&1; then
  echo "ERROR: tag '$VER' already exists." >&2
  exit 1
fi

# ── 4. Monotonic: new version must exceed the current version.py ──────────────
PREV="$(python3 version.py)"
python3 - "$PREV" "$VER" <<'PY'
import sys
try:
    from packaging.version import Version
except ModuleNotFoundError:
    sys.exit(0)  # packaging absent locally; CI's version.yml enforces this authoritatively
prev, new = Version(sys.argv[1]), Version(sys.argv[2])
if new <= prev:
    sys.exit(f"ERROR: {new} is not greater than the current version {prev}.")
print(f"version: {prev} -> {new}")
PY

# ── 5. Bump version.py (__version__ and __version_info__ in lock-step) ─────────
python3 - "$VER" <<'PY'
import re, sys
v = sys.argv[1]
m = re.fullmatch(r"(\d+)\.(\d+)\.(\d+)((?:a|b|rc)\d+)?", v)
maj, mnr, pat, mod = m.group(1), m.group(2), m.group(3), m.group(4)
info = f'({maj}, {mnr}, {pat}, "{mod}")' if mod else f'({maj}, {mnr}, {pat}, None)'
s = open("version.py").read()
s, n1 = re.subn(r'__version__ = "[^"]*"', f'__version__ = "{v}"', s, count=1)
s, n2 = re.subn(r'__version_info__ = \([^)]*\)', f'__version_info__ = {info}', s, count=1)
if n1 != 1 or n2 != 1:
    sys.exit("ERROR: could not rewrite version.py (expected one __version__ and one __version_info__).")
open("version.py", "w").write(s)
PY

# ── 6. Cheap CI gates (the same ones release.yml's `unit` job runs) ───────────
echo "== gate: version parses =="
python3 -c "import sys; sys.path.insert(0,'.'); from version import __version__, parse_version; parse_version(__version__); print('  ok:', __version__)"
echo "== gate: tests/test_classification.py =="
python3 tests/test_classification.py
echo "== gate: tests/test_manifest.py =="
python3 tests/test_manifest.py
echo "== gate: tests/test_classification.R (best-effort) =="
if command -v Rscript >/dev/null 2>&1; then
  if Rscript tests/test_classification.R; then :; else
    echo "  WARN: local R mirror test failed/unavailable — CI runs it authoritatively." >&2
  fi
else
  echo "  WARN: Rscript not on PATH — skipping; CI runs the R mirror test." >&2
fi

# ── 7. Commit + annotated tag (no push) ───────────────────────────────────────
git add version.py
git commit -m "Release $VER"
git tag -a "$VER" -m "CARP $VER"

cat <<EOF

────────────────────────────────────────────────────────────────────────────
Release $VER prepared locally: commit $(git rev-parse --short HEAD) + tag $VER.
The sandbox has no ssh, so push from the HOST terminal:

    git push origin main
    git push origin $VER

The tag push runs .github/workflows/release.yml:
  unit → build SIF → in-container medium fixture → GHCR push
       → GitHub Release → Zenodo snapshot (permanent DOI).
Nothing is published unless the in-container fixture passes.
────────────────────────────────────────────────────────────────────────────
EOF

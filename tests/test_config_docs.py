#!/usr/bin/env python3
"""Docs-in-sync gate: config parameters and the CHANGELOG.

Two checks, run by ``.claude/skills/release/cut-release.sh`` and CI (unit.yml):

1. **Every config parameter the Snakefile reads is documented** in
   ``docs/configuration.md``. The Snakefile is the single source of truth —
   each ``config["X"]`` it reads is a user-facing knob. This is the guard that
   turns "added a config option but forgot to document it" into a build/release
   failure (it was added after exactly that slip for the culling knobs).

2. **CHANGELOG.md is maintained** — its top section is either ``## Unreleased``
   (pending changes between releases) or ``## <version.py>`` (labelled at
   release). Keeps release notes from silently lagging the code.

Exit 0 on success, 1 with an actionable message on failure.
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent

# Params intentionally NOT in the user-facing config reference. Keep empty;
# add a name here (with a comment justifying it) only for a genuinely internal
# key that must never appear in docs/configuration.md.
ALLOWLIST = set()


def snakefile_config_params():
    """Every distinct ``config["X"]`` / ``config['X']`` key read by the Snakefile."""
    txt = (ROOT / "Snakefile").read_text()
    return set(re.findall(r"""config\[["']([a-z_0-9]+)["']\]""", txt))


def documented_names():
    """Backtick-wrapped tokens in docs/configuration.md (param names appear as
    table-row keys ``| `name` | … |`` and in prose)."""
    txt = (ROOT / "docs" / "configuration.md").read_text()
    return set(re.findall(r"`([a-z_0-9]+)`", txt))


def check_config_params():
    read = snakefile_config_params()
    documented = documented_names()
    missing = sorted(p for p in read if p not in documented and p not in ALLOWLIST)
    if missing:
        print("FAIL: config parameters read by the Snakefile but NOT documented "
              "in docs/configuration.md:")
        for p in missing:
            print(f"  - {p}")
        print("\nFix: add each to docs/configuration.md (the right ## section) — "
              "and to config_full.yaml / the README table if user-facing — or, "
              "for a genuinely internal key, add it to ALLOWLIST in this test.")
        return False
    print(f"OK: all {len(read)} Snakefile config parameters are documented in "
          "docs/configuration.md.")
    return True


def check_changelog(release_version=None):
    sys.path.insert(0, str(ROOT))
    from version import __version__

    txt = (ROOT / "CHANGELOG.md").read_text()
    headers = re.findall(r"^##\s+(.+?)\s*$", txt, re.MULTILINE)
    if not headers:
        print("FAIL: CHANGELOG.md has no '## <version>' / '## Unreleased' sections.")
        return False

    if release_version is not None:
        # Release mode (cut-release.sh): the CHANGELOG must carry a labelled
        # section for the version being cut — forces renaming ## Unreleased ->
        # ## <VERSION> so the release notes never lag the tag.
        if release_version in headers:
            print(f"OK: CHANGELOG has a '## {release_version}' section for the release.")
            return True
        print(f"FAIL: cutting {release_version} but CHANGELOG.md has no "
              f"'## {release_version}' section (top is '## {headers[0]}').\n"
              f"Fix: rename the '## Unreleased' section to '## {release_version}' "
              f"(and make sure it lists this release's user-visible changes) "
              f"before cutting.")
        return False

    # Default (CI / dev): top section is Unreleased (pending) or the current
    # released version.
    top = headers[0]
    if top == "Unreleased" or top == __version__:
        print(f"OK: CHANGELOG top section is '## {top}' (version.py = {__version__}).")
        return True
    print(f"FAIL: CHANGELOG top section '## {top}' is neither 'Unreleased' nor "
          f"the current version '{__version__}'.\n"
          f"Fix: open an '## Unreleased' section for pending changes, or label it "
          f"'## {__version__}' when cutting the release.")
    return False


def main():
    release_version = None
    argv = sys.argv[1:]
    if argv and argv[0] == "--release":
        release_version = argv[1] if len(argv) > 1 else None
        if release_version is None:
            print("usage: test_config_docs.py [--release VERSION]")
            sys.exit(2)
    ok = True
    ok = check_config_params() and ok
    ok = check_changelog(release_version) and ok
    if ok:
        print("OK: test_config_docs — docs are in sync with the config + version.")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

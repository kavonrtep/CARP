---
name: changelog
description: Maintain CARP's CHANGELOG.md. Invoke when a commit changes user-visible pipeline behaviour (add a bullet under ## Unreleased) and at release time (label the ## Unreleased section with the version being cut). Keeps release notes from lagging the code.
tools: Read, Edit, Grep, Bash
---

You maintain `CHANGELOG.md` for CARP.

**Format** (keep-a-changelog, undated headers):
```
# Changelog
## Unreleased          # pending, not-yet-released changes (may be absent)
## MAJOR.MINOR.PATCH[modifier]   # released sections, newest first (e.g. 1.0.0, 1.0.0rc3)
```
The docs-in-sync gate `tests/test_config_docs.py` requires the top `##` section
to be `## Unreleased` or `## <version.py>`; in release mode
(`--release VERSION`, run by `cut-release.sh`) it requires a `## VERSION`
section to exist.

**When a commit changes behaviour** (new/changed config option, fixed bug,
changed output, dependency bump):
- Add one terse bullet under `## Unreleased` (create that section directly under
  the `# Changelog` title if it is missing). Match the existing bullets' voice —
  a bold lead-in, technical, one measured number where it helps, and a link to
  the relevant `docs/` or spec page if there is one. For a new config option,
  name it, its default, and its effect (and confirm `config-docs` has documented
  it in `docs/configuration.md`).

**At release time** (version being cut = VERSION):
- Rename the top `## Unreleased` header to `## VERSION`. If there is no
  Unreleased section, create `## VERSION` and populate it from the commits since
  the previous tag (`git log --oneline <prev-tag>..HEAD`), summarising the
  user-visible ones. Then a fresh `## Unreleased` is only added when the next
  change lands.

Constraints: never rewrite the history of already-released sections; never edit
released tags. Return a short summary of the bullets you added / the rename you
made.

"""Single source of truth for the CARP pipeline version.

Format follows PEP 440 with the lite subset accepted by this project:
``MAJOR.MINOR.PATCH`` optionally followed by an ``aN`` / ``bN`` / ``rcN``
modifier (alpha / beta / release candidate). Examples:

    0.9.0rc1   0.9.0rc2   0.9.0   1.0.0a1   1.0.0

Bump rule: every commit on ``main`` that changes pipeline behaviour
should also bump the version. The CI version-gate workflow
(``.github/workflows/version.yml``, introduced in commit 2 of the
versioning rollout) enforces:

  * the string here parses cleanly via ``parse_version`` below,
  * a tag push (``v<X.Y.Z…>``) matches the value here exactly,
  * the version on ``main`` is monotonically non-decreasing.

This module is read by:

  * ``run_pipeline.py`` (``--version`` flag, startup banner)
  * the Singularity recipe (stamps the ``Version`` label at build time)
  * ``scripts/record_provenance.py`` (run-time provenance JSON)
  * the CI workflows (gate)
  * the per-run HTML report footer (commit 5 of the rollout)

Keep the module minimal — no heavy imports at module level — so reading
the version stays cheap.
"""

__version__ = "0.9.0rc3"
__version_info__ = (0, 9, 0, "rc3")


def parse_version(s):
    """Parse a PEP-440-lite version string.

    Returns ``(major, minor, patch, modifier)`` where ``modifier`` is
    one of ``None``, ``"aN"``, ``"bN"``, ``"rcN"`` for some integer N.

    Raises ``ValueError`` for anything that does not match the
    project's accepted shape — epochs (``2!1.0``), post-releases
    (``1.0.post1``), dev releases (``1.0.dev1``), and local-version
    suffixes (``1.0+local``) are intentionally rejected so the CI gate
    catches non-conforming strings before they reach a release tag.
    """
    import re

    m = re.match(r"^(\d+)\.(\d+)\.(\d+)((?:a|b|rc)\d+)?$", s.strip())
    if m is None:
        raise ValueError(
            f"not a PEP 440 lite version: {s!r} "
            f"(expected MAJOR.MINOR.PATCH with optional aN/bN/rcN modifier)"
        )
    return (int(m.group(1)), int(m.group(2)), int(m.group(3)), m.group(4))


if __name__ == "__main__":
    # Allow `python version.py` to print the bare version string —
    # convenient in shell pipelines and the Singularity %post stamp
    # without needing a full Python import.
    print(__version__)

# CARP development

Developer notes: building the container, cutting a release, and regenerating the
workflow diagram. End users do not need this — see the [README](../README.md).

## Build the container

Release builds are produced by GitHub Actions
([`.github/workflows/release.yml`](../.github/workflows/release.yml)) on every
tag push matching the unprefixed PEP 440 pattern (`1.0.0rc3`, `1.0.0`, …). The
workflow runs unit tests, builds the SIF from the `Singularity` recipe,
exercises the medium fixture **inside** the freshly-built container, then — only
on a green test-in-container — pushes the image to GHCR and creates a GitHub
Release that triggers the Zenodo snapshot.

To rebuild the same image locally:

```bash
TAG=$(python version.py)                          # e.g. 1.0.0rc3
sudo apptainer build images/carp_${TAG}.sif Singularity
```

The version stamped into the image is read from `version.py`; CI's `version.yml`
workflow enforces that pushed tags match `version.py` exactly, so the local build
is identical to what the release workflow publishes.

## Releasing

`version.py` is the single source of truth; a release is a version bump plus an
**unprefixed** tag equal to it (`1.0.0rc1`, not `v1.0.0rc1`). Use the `/release`
skill or its helper, which validates, bumps `version.py`, runs the cheap CI
gates, commits, and tags (no push):

```bash
.claude/skills/release/cut-release.sh <VERSION>   # e.g. 1.0.0rc1
```

Then push from the host: `git push origin main && git push origin <VERSION>`.
The tag push drives `release.yml` (SIF build → in-container fixture → GHCR →
GitHub Release → Zenodo).

## Workflow diagram

The schematic at [`figs/workflow_overview.svg`](../figs/workflow_overview.svg)
(and `.png`) is derived **live from the Snakefile** —
`scripts/make_workflow_diagram.py` runs `snakemake --rulegraph`, groups rules
into analysis stages, transitively reduces the dependency graph, and labels each
stage with an SVG tooltip. Regenerate after changing the workflow:

```bash
scripts/make_workflow_diagram.py   # writes figs/workflow_overview.{dot,svg,png}
```

Requires `snakemake` and Graphviz `dot` on PATH. The script fails loudly if a new
rule is not assigned to a stage, so the figure cannot silently drift.

# CI / release policy

Three GitHub Actions workflows, each with a distinct trigger and cost
profile. They exist as separate files on purpose — folding them into one
workflow would couple fast feedback to slow jobs and defeat caching.

## Workflows

### `unit.yml` — every push / PR (≤ 1 min)

Runs the classification vocabulary tests in Python and R plus a
self-test of the `.githooks/pre-commit` hook. Uses micromamba with
`cache-environment: true` so the env is reused across runs. No
Snakemake, no pipeline, no container — just fast signal.

### `pipeline.yml` — PR into `main`, push to `main`, Sunday nightly (≤ 15 min warm, ≤ 45 min cold)

Runs the full Snakemake pipeline on both `small` and `medium` fixtures
using plain `snakemake --use-conda`, without the container. The
`.snakemake/conda/` directory is cached keyed on `hashFiles('envs/*.yaml')`
so only env-spec edits invalidate the cache; cold env creation (first
run or after an env edit) is a one-time ~30 min cost.

Asserts are in `scripts/assert_fixture_outputs.sh` — centralised
thresholds that let us catch regressions without having to diff GFFs
byte-for-byte. On failure the pipeline output directory is uploaded as
an artefact (14 d retention) for post-mortem.

### `release.yml` — tag push `v*`, manual dispatch (30–90 min)

This is the important one. It's a chain of dependent jobs, and the
chain intentionally ends with `gh release create`, because **creating
the GitHub Release is what triggers the Zenodo snapshot** (the repo is
enabled in Zenodo's GitHub integration).

```
unit  →  build-sif  →  test-in-container  →  push-ghcr  →  create-release
```

- **Every job is a hard gate.** If `test-in-container` fails, no
  `push-ghcr`, no `create-release`, no Zenodo snapshot.
- **Container is the authoritative test subject.** The pipeline test
  runs on the `.sif` artefact that will actually ship, not on a
  snakemake-with-conda proxy.
- **SIF hosting is GHCR, not Release assets.** SIF files are ~3 GB and
  blow past GitHub's 2 GB per-asset limit. GHCR stores them as OCI
  artefacts via ORAS; consumers pull with `apptainer pull
  oras://ghcr.io/kavonrtep/carp/sif:<TAG>`.
- **Two GHCR tags per release.** The version tag (`v0.8.16`) and
  `:latest`. The release notes only advertise the versioned form; users
  who want floating can use `:latest`.

## Why the duplication between `unit.yml` and the `unit` job in `release.yml`

The release pipeline has to be self-contained — triggering off a tag
push should not depend on an earlier push-based `unit.yml` run having
succeeded. If a developer tags a commit that was never pushed to a
branch (plausible during hotfixes), the tag-push workflow must still
gate on unit tests. Duplication is cheap (<1 min) and eliminates a
class of "release without unit coverage" mistakes.

## Zenodo integration — how the gating actually works

1. Developer pushes a tag: `git push origin v0.8.16`.
2. `release.yml` starts. Zenodo sees nothing yet — it only listens for
   **Release** events, not tags.
3. The chain runs: unit → build-sif → test-in-container → push-ghcr.
   If anything fails, the workflow stops and no Release is created.
   The tag exists on the repo but is dormant; Zenodo does nothing.
4. If everything passes, `create-release` publishes the GitHub Release
   for that tag. Zenodo's webhook fires and archives the repository
   snapshot with a DOI.

So: **broken tags never land on Zenodo.** If a release build fails, fix
the problem, delete the tag (`git push --delete origin v0.8.16`, `git
tag -d v0.8.16`), re-tag the fixed commit, push again.

## GHCR authentication

`release.yml` uses the workflow's default `GITHUB_TOKEN` with
`packages: write` permission (set at workflow scope). No PATs or
secrets needed. This only works for pushes to the same repository's
GHCR namespace — if we ever need to publish under a different org we'd
need a fine-grained PAT in `secrets.GHCR_TOKEN`.

## Cost / disk budget

| workflow | avg runtime (warm) | runner | CPU | RAM | disk |
|---|---:|---|---:|---:|---:|
| unit | 40 s | ubuntu-latest | 2 | 7 GB | 1 GB |
| pipeline (per fixture) | 7 min | ubuntu-latest | 2 | 7 GB | ~10 GB (conda cache) |
| release | 55 min | ubuntu-latest | 2 | 7 GB | ~10 GB (SIF build) |

The SIF build frees ~10 GB of disk (dotnet, GHC, CodeQL) before
invoking `apptainer build` to make room for the conda env solve.

## Adding a new fixture

1. Check in the new FASTA under `tests/fixtures/genome_<name>.fasta`
   (note: fixture output dirs are `.gitignore`d).
2. Add a `config_<name>.yaml` sibling.
3. Add a `case` branch to `scripts/assert_fixture_outputs.sh` with its
   threshold expectations.
4. Add the fixture name to the `matrix.fixture` in `pipeline.yml`.

If the new fixture should also be the release-gate fixture, replace the
`config_medium.yaml` reference in `release.yml`'s `test-in-container`
job. Generally don't add fixtures to the release gate — one is enough,
and each extra fixture doubles the release chain's runtime.

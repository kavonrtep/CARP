# DANTE-family & TideCluster: constrain runtime dependencies in the conda recipes

## Summary

The DANTE-family tools (**DANTE**, **DANTE_LTR**, **DANTE_TIR**) and
**TideCluster** are distributed as conda packages, and CARP installs them by
**exact version** (`dante=0.2.5`, `dante_ltr=0.4.0.4`, `dante_tir=0.2.6`,
`tidecluster=1.16.3`). But the recipes place **no version constraints on their
runtime dependencies**, so `conda install <tool>=<version>` still re-solves the
entire transitive dependency tree to whatever is *current* on
conda-forge/bioconda at build time.

The consequence: pinning the tool version does **not** pin the tool's behaviour.
A container (or HPC module, or fresh env) built from the same recipe today and
next month can contain different `mmseqs2`, `RepeatMasker`, `TideHunter`,
`r-base`, `bioconductor-*`, and `blast` versions — with no change to the tool
version. That makes downstream results non-reproducible, can silently change
annotations, and repeatedly breaks the build.

We would like the recipes to **constrain their result-affecting dependencies**
(version ranges tested per release, or at minimum upper bounds), so that a
consumer who pins the tool version gets a reproducible, tested dependency
closure.

## What we observe

CARP pins the four tool versions above and builds a released Singularity image
on every tag. Because the recipes leave the deps unconstrained, the deps float.
The clearest evidence is a **single CARP container** (release 1.0.7): because
the recipes don't constrain the deps, the *same* result-affecting dependency
resolved to **different versions across the three tools' environments — in one
image, built at one moment**:

| dependency | DANTE / DANTE_LTR env | DANTE_TIR env | TideCluster env | result-affecting? |
|------------|:---------------------:|:-------------:|:---------------:|-------------------|
| `mmseqs2` | `16.747c6` | **`18.8cc5c`** | `16.747c6` | **yes** — clustering (major-version split) |
| `bioconductor-rtracklayer` | `1.54.0` | **`1.62.0`** | **`1.62.0`** | yes — R I/O (8 minor versions apart) |
| `r-base` | `4.1.3` | **`4.3.3`** | **`4.3.3`** | yes — whole R stack |
| `blast` | `2.16.0` | **`2.17.0`** | `2.16.0` | yes |
| `repeatmasker` | `4.1.2.p1` | — | `4.1.2.p1` | yes — masking |
| `tidehunter` | `1.4.3` | — | `1.4.3` | yes — tandem detection |
| `cap3` | `10.2011` | `10.2011` | — | yes |

If a single build can't even resolve `mmseqs2` or `rtracklayer` consistently
across the three tools, then the *same recipe rebuilt next month* has no chance
of reproducing today's versions. (The version each tool happens to get is just
"latest compatible on build day".)

The same class of drift recently broke CARP's release build repeatedly (the
unpinned `continuumio/miniconda3` base drifted to a conda that enforced
Anaconda's channel Terms of Service, then to python 3.14, then to a conda
plugin conflict). We pinned the base image to stop that — but the **tool
recipes' unpinned deps are the same problem one level down**, and only the tool
authors can fix it at the source.

## Why it matters

1. **Reproducibility.** A build's scientific output becomes a function of the
   build *date*, not of the pinned tool versions. Two labs pinning
   `tidecluster=1.16.3` a few months apart can get different annotations.
2. **Silent behaviour drift — `mmseqs2` especially.** mmseqs clustering results
   depend on the version and its defaults. CARP already has to work around
   *"the mmseqs default k=15 segfaults on short tandem monomers"* in its own
   dimer reducer; an unconstrained `mmseqs2` bump inside DANTE_TIR
   (`--max_class_size` AA clustering) or TideCluster can change results or
   reintroduce/alter such edge cases with no version signal. RepeatMasker and
   the bioconductor stack are equally behaviour-defining.
3. **Build breakage.** An incompatible transitive dep can make the solve fail
   outright or the tool crash at runtime — discovered only at build/run time,
   often long after the tool release was cut.

## Requested change

In each tool's conda recipe (`meta.yaml`, `requirements: run:`):

1. **Constrain the result-affecting runtime dependencies** with version ranges
   that were tested for that tool release — at minimum an **upper bound** on the
   ones that break across majors. Priority list:
   - **DANTE_TIR:** `mmseqs2`, `blast`, `cap3`, the R stack it uses, `r-base`.
   - **DANTE_LTR:** `repeatmasker`, `blast`, `bioconductor-{rtracklayer,genomicranges,biostrings}`, `r-base`, `r-igraph`.
   - **DANTE:** `blast`/`last`, `hmmer` (or whatever the domain search uses), `r-base` where applicable.
   - **TideCluster:** `tidehunter`, `mmseqs2`, `repeatmasker`, `bioconductor-rtracklayer`, `r-base`, `r-igraph`.
2. **Or** publish a tested, fully-resolved spec per release (an exported
   `conda list --explicit` / `conda-lock` alongside the package) so consumers
   can build the exact dependency closure the release was validated against.
3. **At minimum**, upper-bound the handful of deps known to break across major
   versions, so pinning the tool version yields a reproducible, tested closure
   rather than "latest wins".

We understand fully-exact pins can make recipes brittle over time (old builds
get garbage-collected from channels). **Version ranges** (`>=x,<y`) tested per
release are the pragmatic middle ground and are what we're really asking for.

## Impact / who hits it

Anyone who installs these tools by version and expects reproducible behaviour —
pipelines, containers, HPC modules, Galaxy wrappers. Fixing it in the recipes
solves it **once, upstream, for every consumer**, instead of each downstream
pipeline hand-pinning and patching its environments (which CARP currently does
only partially, e.g. adding `jq` and env-level pins in its own env files).

For CARP specifically: we pin the four tool versions and ship a released
container per tag; because the recipe deps are unconstrained, that container's
contents — and therefore the repeat annotation it produces — are not fully
determined by the versions we pin. Constraining the recipe deps would let CARP
(and everyone else) treat a pinned tool version as a genuine, reproducible
contract.

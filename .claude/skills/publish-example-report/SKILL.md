# Publish / refresh the CARP example-output website

A **separate** GitHub repo, served by GitHub Pages, hosts one real, fully
explorable CARP report so users (biologists) can see what CARP produces without
running it. It is refreshed **on demand** — not on every release — because the
report format changes rarely and the example genome is large. Keeping it in its
own repo leaves CARP's git history untouched (the bundle is ~30 MB, mostly the
TideCluster per-family drill-down).

`build-example-site.sh` assembles a portable, self-contained site from any CARP
output directory: the main `repeat_annotation_report.html` + `plotly.min.js`, the
`DANTE_LTR` structural summary, and the whole self-contained `TideCluster/default`
report tree, plus a biologist-facing `index.html` landing page whose metadata
(CARP version, run date, genome accession, size) is read from the run's
`run_provenance.json`. No heavy raw data (masked FASTA, bigwigs, GFF3s) is
included. It finishes as a single-commit git repo, ready to force-push.

## One-time setup
1. Create a new **empty, public** GitHub repo, e.g. `kavonrtep/carp-example-report`.
2. Pages is enabled after the first push (step below).

## Build the site (in the sandbox or host)
The output dir doubles as the local clone of the example repo. `carp-example-report/`
(a sibling-style subdir of the CARP checkout) is gitignored in the CARP repo, so
it never pollutes CARP's tree. Only that dir is writable in the sandbox; on the
host you may prefer a true sibling path.
```bash
.claude/skills/publish-example-report/build-example-site.sh \
    -i <carp_output_dir> \
    -o carp-example-report \
    -s "Boechera stricta"        # species label for the landing page
# -a GCA_XXXX.N  overrides the accession (else auto-detected from provenance;
#                a GCA/GCF accession is auto-linked to the ENA browser)
```
Preview locally by opening `carp-example-report/index.html` in a browser. Each
build is a fresh single commit but an existing `origin` remote is preserved, so
refresh is just rebuild + `git push -f`.

## Publish — from the HOST (the sandbox has no ssh)
```bash
cd carp-example-report
git remote add origin git@github.com:kavonrtep/carp-example-report.git   # first time only
git push -f origin main
```
Then, **first time only**, enable Pages: repo → Settings → Pages → Source
"Deploy from a branch" → Branch `main` / `/ (root)` → Save.
Live at: `https://kavonrtep.github.io/carp-example-report/`
(consider linking that URL from the main CARP `README.md`).

## Refresh on demand
Run a fresh CARP job on the chosen genome (e.g. after a change that alters the
report), then rebuild + force-push:
```bash
.claude/skills/publish-example-report/build-example-site.sh -i <new_output> -o carp-example-report -s "…"
cd carp-example-report && git push -f origin main
```
Each build is a fresh single commit (no history accumulation), so the example
repo stays flat at ~30 MB and Pages redeploys automatically on push.

## Notes / gotchas
- `.nojekyll` is written so Pages serves the asset dirs verbatim (no Jekyll).
- The site is self-contained: every `src`/`href` stays within the site root
  (only outbound links are to `github.com/kavonrtep/carp`). If a future report
  adds a new external asset, extend the copy list in `build-example-site.sh`.
- `cp -rL` dereferences the top-level report symlinks so the bundle has real files.
- This lives under `.claude/` so it is **not** copied into the release SIF
  (`Singularity` `%files`), and it is not a pipeline rule — run it by hand.

# TideCluster report v2: make `<prefix>_report/` self-contained (vendor the PNGs)

**Status:** ✅ RESOLVED in **TideCluster 1.17.0** (petrnovak channel); CARP's
`envs/tidecluster_run.yaml` pin bumped 1.16.6 → 1.17.0. The report now vendors
its per-TRC PNGs into `<prefix>_report/img/`. Verified end-to-end on the small
fixture: after `cleanup_intermediates: maximal` purges `TideCluster_kite` /
`TideCluster_tarean` / `TideCluster_consensus`, 0 of 9 report PNG references are
missing — the report is fully self-contained.

## Summary

The modernised HTML report (`<prefix>_report/`, built by `tc_rerender_report.py`)
**vendors its CSS/JS assets** into `<prefix>_report/assets/` and is self-contained
*for styling* — but it **references its per-TRC images by reaching OUTSIDE the
report directory** with `../`, into the sibling `TideCluster_kite/` and
`TideCluster_tarean/` working trees:

- `kite.html` → `../TideCluster_kite/profile_plots/profile_top3_<TRC>.png`
- `tarean.html` (and the legacy `<prefix>_report_legacy/TideCluster_tarean_report.html`)
  → `../TideCluster_tarean/<TRC>.fasta_tarean/img/{graph,logo}_<k>mer_<n>.png`

Those two trees are large intermediate/working directories (kite ≈ 3 GB, tarean
≈ 1.6 GB on a real plant genome), so a downstream consumer that reclaims disk by
deleting them **breaks the report** — every image renders as a broken link even
though the HTML pages themselves survive.

**Request:** vendor the referenced PNGs into the report directory the same way
the CSS/JS assets already are (e.g. `<prefix>_report/img/`) and rewrite the
`<img src>` / `<a href>` to point there. Then `<prefix>_report/` is a fully
portable, self-contained deliverable, and the `TideCluster_kite` /
`TideCluster_tarean` working trees can be deleted independently without
affecting the report.

## Why this matters (downstream context)

CARP (the Assembly Repeat Annotation Pipeline) treats `TideCluster_kite/` and
`TideCluster_tarean/` as reclaimable working data and offers an aggressive
post-run cleanup (`cleanup_intermediates: maximal`) that deletes them — they are
multi-GB and hold `.kmers`, dimer FASTAs, periodograms, etc. that nothing
downstream consumes. The **report images are the only reason those trees cannot
be deleted**, and the images are a tiny fraction of the trees:

| tree | total | report-referenced images |
|------|-------|--------------------------|
| `TideCluster_kite` | ~3.0 GB | `profile_plots/*.png` — **~2.1 MB** |
| `TideCluster_tarean` | ~1.6 GB | all `<TRC>.fasta_tarean/img/*.png` — **~57 MB** |

(Real genome, ~2.8 Gb assembly: 254 kite profile PNGs + 564 tarean img PNGs.)

So vendoring the images into the report costs ≈ 59 MB, after which the ≈ 4.5 GB
of genuine scratch is freely deletable. This is a general portability win too:
today the report cannot be moved/zipped/served on its own because its images
live in sibling directories.

## Evidence

Report page → image reference (verbatim from a real run):

```
# <prefix>_report/kite.html
href="../TideCluster_kite/profile_plots/profile_top3_TRC_1.png"

# <prefix>_report/tarean.html   (and <prefix>_report_legacy/TideCluster_tarean_report.html)
href="../TideCluster_tarean/TRC_1.fasta_tarean/img/graph_11mer_1.png"
href="../TideCluster_tarean/TRC_1.fasta_tarean/img/logo_11mer_1.png"
```

Directory layout (per clustering run, i.e. both `default/` and `short_monomer/`):

```
<outdir>/TideCluster/<default|short_monomer>/
  TideCluster_index.html          # landing page (redirect target from the top level)
  TideCluster_report/             # report v2 — self-contained for CSS/JS…
    assets/                       #   …assets ARE vendored here
    index.html trcs.html kite.html tarean.html superfamilies.html
    trc/TRC_<N>.html              #   per-TRC dashboards (one dir deeper)
  TideCluster_kite/               # ← report reaches in here via ../
    profile_plots/*.png
  TideCluster_tarean/             # ← and here via ../
    <TRC>.fasta_tarean/img/*.png
```

## Where in the code (TideCluster 1.16.x, `share/tidecluster/tc_rerender_report.py`)

The renderer already copies `assets/` into the report dir — the same mechanism
just needs to cover the PNGs:

- **Assets are vendored** (the pattern to extend): `shutil.copytree(src, dst)`
  around **line 345**, referenced with `assets_href="assets/"` (self-contained).
- **KITE image paths** are built relative to the *sibling* kite dir, not the
  report dir:
  - **`kite_dir_rel = paths["kite_dir"].name`** (line ~1313) — the bare dir name,
    so it resolves as a sibling.
  - **line ~1396-1397**:
    `profile_png = f"{kite_dir_rel}/profile_plots/profile_{trc_id}.png"` and
    `profile_top3_png = f"{kite_dir_rel}/profile_plots/profile_top3_{trc_id}.png"`.
- **TAREAN image paths**: `graph_link` / `logo_link` (line ~1233-1234), rendered
  by `_rel_img(...)` / `_thumb(...)` with `src_prefix="../"` (top-level pages) or
  `"../../"` (per-TRC dashboards, one dir deeper) — **lines ~1469, ~2017**. That
  `../` prefix is exactly what points outside the report dir.

## Requested change

1. When building `<prefix>_report/`, **copy the referenced PNGs into the report
   dir** — e.g. `<prefix>_report/img/kite/profile_top3_<TRC>.png` and
   `<prefix>_report/img/tarean/<TRC>/{graph,logo}_<k>mer_<n>.png` — mirroring the
   existing `assets/` `copytree`. (Copy only the PNGs actually referenced, not the
   whole trees, so the report stays ~59 MB, not multi-GB.)
2. **Rewrite the references** to the vendored copies: drop the `../<tree>/` prefix
   and the `src_prefix="../"`/`"../../"` reach-out, pointing at the local `img/…`
   path (with the correct relative depth for top-level pages vs `trc/` dashboards).
3. Result: `<prefix>_report/` opens correctly with **nothing outside it required**.

Nice-to-haves / questions for you:
- **Legacy report** (`<prefix>_report_legacy/TideCluster_tarean_report.html`)
  also references `../TideCluster_tarean/*/img/*.png`. If the legacy report is
  still shipped, it needs the same treatment (or its own vendored `img/`);
  if it is deprecated, that's fine — CARP consumes the v2 report.
- A `--self-contained-report` toggle is not needed from our side; making it the
  default behaviour is strictly better (the ~59 MB duplication is negligible next
  to the trees, and portability is a clear win). But a flag is acceptable if you
  prefer to keep the reach-out layout as an option.
- Please keep the on-disk image files inside the trees unchanged (TideCluster's
  own re-runs may still use them); this request is only about the report having
  its *own* copies so it no longer depends on them.

## Benefit

- The report becomes a portable, self-contained deliverable (movable / zippable /
  servable on its own) — independent of TideCluster's working directories.
- Downstream tools (CARP's `maximal` cleanup) can reclaim the ~4.5 GB of genuine
  scratch in `TideCluster_kite` / `TideCluster_tarean` **without breaking the
  report**. Current state: CARP's `maximal` mode deletes those trees, so the
  report's images break under `maximal` (it is opt-in; CARP's default `minimal`
  mode never touches them). Once the report is self-contained, `maximal` becomes
  safe and CARP will let it purge the trees whole.

(Author == platform operator, so this is directly actionable in TideCluster.)

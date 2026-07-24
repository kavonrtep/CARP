#!/usr/bin/env bash
# Assemble a portable, self-contained CARP example-report website from a CARP
# output directory, ready to publish to GitHub Pages (a SEPARATE repo).
#
# It copies only what the reports actually reference (the main report + plotly +
# the DANTE_LTR summary subtree + the whole TideCluster/default subtree, which is
# self-contained), and generates a biologist-facing index.html landing page whose
# metadata (CARP version, run dates, genome accession) is read from the run's
# run_provenance.json. No heavy raw data (masked FASTA, bigwigs) is included.
#
# Usage:
#   build-example-site.sh -i <carp_output_dir> -o <site_dir> \
#       [-s "Species name"] [-a GCA_XXXXXXXXX.N] [--no-git]
#
# Then publish from the HOST (the sandbox has no ssh):
#   see SKILL.md — create the example repo once, enable Pages, force-push <site_dir>.
set -euo pipefail

IN="" ; OUT="" ; SPECIES="" ; ACCESSION="" ; DO_GIT=1
while [ $# -gt 0 ]; do
  case "$1" in
    -i|--input)     IN="$2"; shift 2;;
    -o|--output)    OUT="$2"; shift 2;;
    -s|--species)   SPECIES="$2"; shift 2;;
    -a|--accession) ACCESSION="$2"; shift 2;;
    --no-git)       DO_GIT=0; shift;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$IN" ] && [ -n "$OUT" ] || { echo "usage: build-example-site.sh -i <carp_output_dir> -o <site_dir> [-s species] [-a accession] [--no-git]" >&2; exit 2; }
[ -d "$IN" ] || { echo "input dir not found: $IN" >&2; exit 1; }
[ -f "$IN/repeat_annotation_report.html" ] || { echo "no repeat_annotation_report.html in $IN" >&2; exit 1; }

# ---- metadata from run_provenance.json (best-effort; falls back gracefully) ----
prov="$IN/run_provenance.json"
jget() { grep -oE "\"$1\" *: *\"[^\"]*\"" "$prov" 2>/dev/null | head -1 | sed 's/.*: *"//; s/"$//'; }
VERSION="$(jget pipeline_version)"; [ -n "$VERSION" ] || VERSION="unknown"
RUN_FIN="$(jget run_finished)"; RUN_DATE="${RUN_FIN%%T*}"; [ -n "$RUN_DATE" ] || RUN_DATE="$(jget run_started | cut -dT -f1)"
GPATH="$(jget genome_fasta)"; GBASE="$(basename "${GPATH:-}")"
[ -n "$ACCESSION" ] || ACCESSION="$(printf '%s' "$GBASE" | grep -oE 'GC[AF]_[0-9]+\.[0-9]+' | head -1)"
[ -n "$ACCESSION" ] || ACCESSION="(unknown accession)"
[ -n "$SPECIES" ] || SPECIES="(species — set with -s)"

CARP_URL="https://github.com/kavonrtep/assembly_repeat_annotation_pipeline"
# link the accession to the ENA browser when it is a real GCA/GCF accession
if printf '%s' "$ACCESSION" | grep -qE '^GC[AF]_[0-9]+\.[0-9]+$'; then
  ACCESSION_HTML="<a href=\"https://www.ebi.ac.uk/ena/browser/view/${ACCESSION}\">${ACCESSION}</a>"
else
  ACCESSION_HTML="${ACCESSION}"
fi

# genome facts from a .fai if present
FAI="$IN/genome_cleaned.fasta.fai"
GSTAT=""
if [ -f "$FAI" ]; then
  GSTAT="$(awk '{n++; bp+=$2; if($2>1000000)big++} END{printf "%.0f Mb across %d sequences (%d chromosome-scale)", bp/1e6, n, big}' "$FAI")"
fi

echo "Assembling example site:"
echo "  species=$SPECIES  accession=$ACCESSION  CARP=$VERSION  run=$RUN_DATE"
echo "  genome: ${GSTAT:-n/a}"

# ---- assemble site (preserve the relative paths the main report links) ----
# Preserve an existing 'origin' remote across rebuilds so refresh is just
# build + `git push -f` (the fresh single commit keeps the repo flat).
PREV_REMOTE=""
[ -d "$OUT/.git" ] && PREV_REMOTE="$(git -C "$OUT" remote get-url origin 2>/dev/null || true)"
rm -rf "$OUT"; mkdir -p "$OUT"
cp -rL "$IN/repeat_annotation_report.html" "$OUT/"
cp -rL "$IN/plotly.min.js" "$OUT/"
# DANTE_LTR structural summary (linked as DANTE_LTR/DANTE_LTR_summary.html)
if [ -f "$IN/DANTE_LTR/DANTE_LTR_summary.html" ]; then
  mkdir -p "$OUT/DANTE_LTR"
  cp -rL "$IN/DANTE_LTR/DANTE_LTR_summary.html" "$OUT/DANTE_LTR/"
  [ -d "$IN/DANTE_LTR/DANTE_LTR_summary_plots" ] && cp -rL "$IN/DANTE_LTR/DANTE_LTR_summary_plots" "$OUT/DANTE_LTR/"
  [ -f "$IN/DANTE_LTR/DANTE_LTR_summary.csv" ] && cp -rL "$IN/DANTE_LTR/DANTE_LTR_summary.csv" "$OUT/DANTE_LTR/"
fi
# DANTE_TIR DNA-transposon summary (linked as DANTE_TIR/report.html; refs img/*.png)
if [ -f "$IN/DANTE_TIR/report.html" ]; then
  mkdir -p "$OUT/DANTE_TIR"
  cp -rL "$IN/DANTE_TIR/report.html" "$OUT/DANTE_TIR/"
  [ -d "$IN/DANTE_TIR/img" ] && cp -rL "$IN/DANTE_TIR/img" "$OUT/DANTE_TIR/"
fi
# TideCluster satellite/tandem report (linked as TideCluster/default/TideCluster_index.html);
# the whole default/ subtree is self-contained (its ../ links stay inside it).
if [ -f "$IN/TideCluster/default/TideCluster_index.html" ]; then
  mkdir -p "$OUT/TideCluster"
  cp -rL "$IN/TideCluster/default" "$OUT/TideCluster/"
  # drop bulky non-report intermediates that no report page links
  rm -f "$OUT"/TideCluster/default/*.log "$OUT"/TideCluster/default/*consensus_dimer_library*.fasta 2>/dev/null || true
fi
# provenance JSON (linked from the main report) — sanitise internal infra details
# (host, user, and absolute /nfsroot cluster paths) before it goes on a public page.
if [ -f "$IN/run_provenance.json" ]; then
  cp -rL "$IN/run_provenance.json" "$OUT/run_provenance.json"
  sed -i -E \
    -e 's#("genome_fasta": ")[^"]*/([^"/]+)"#\1\2"#' \
    -e 's#("container_sif": ")[^"]*/([^"/]+)"#\1\2"#' \
    -e 's#("output_dir": ")[^"]*"#\1(redacted)"#' \
    -e 's#("host": ")[^"]*"#\1(redacted)"#' \
    -e 's#("user": ")[^"]*"#\1(redacted)"#' \
    "$OUT/run_provenance.json"
fi

# ---- biologist-facing landing page ----
cat > "$OUT/index.html" <<HTML
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CARP example output — ${SPECIES} (${ACCESSION})</title>
<style>
  :root { color-scheme: light dark; }
  body { font: 16px/1.6 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
         max-width: 820px; margin: 0 auto; padding: 2.5rem 1.25rem; }
  h1 { font-size: 1.7rem; margin: 0 0 .25rem; }
  .sub { color: #666; margin: 0 0 1.5rem; }
  .meta { font-size: .9rem; color: #666; border:1px solid #8884; border-radius:8px;
          padding:.75rem 1rem; margin:1.25rem 0; }
  .cards { display:grid; gap:.9rem; margin:1.5rem 0; }
  a.card { display:block; border:1px solid #8886; border-radius:10px; padding:1rem 1.2rem;
           text-decoration:none; color:inherit; }
  a.card:hover { border-color:#3b82f6; }
  a.card.primary { border-width:2px; border-color:#3b82f6; }
  a.card b { font-size:1.05rem; }
  a.card span { display:block; color:#666; font-size:.9rem; margin-top:.2rem; }
  footer { margin-top:2rem; font-size:.85rem; color:#888; border-top:1px solid #8883; padding-top:1rem; }
  code { background:#8881; padding:.1em .35em; border-radius:4px; }
  .note { background:#f59e0b1a; border:1px solid #f59e0b66; border-radius:8px;
          padding:.8rem 1rem; font-size:.92rem; margin:1.25rem 0; }
</style>
</head>
<body>
  <h1>CARP — example repeat annotation</h1>
  <p class="sub"><em>${SPECIES}</em> &middot; ${ACCESSION_HTML}</p>

  <p>This is a real example of the interactive HTML report produced by
  <a href="${CARP_URL}">CARP</a> (Comprehensive Annotation of Repeats Pipeline)
  on a whole genome. It shows how CARP classifies and quantifies the repetitive
  fraction — LTR retrotransposons, DNA transposons, LINEs, tandem repeats /
  satellites and rDNA — combining structural (protein-domain) detection with
  similarity-based annotation.</p>

  <p class="note"><strong>This is an example report only, not a full analysis.</strong>
  It contains the HTML reports below. A complete CARP run also produces the full
  annotation as <strong>GFF3 files</strong> (per class and unified), genome-wide
  and per-family repeat-density tracks (BigWig), summary tables and soft-masked
  sequence — those data files are <strong>not</strong> included here.</p>

  <div class="cards">
    <a class="card primary" href="repeat_annotation_report.html">
      <b>▶ Repeat annotation report</b>
      <span>Genome-wide repeat composition, interactive class breakdown, genomic
      distribution, and links into the structural and tandem-repeat sub-reports.</span>
    </a>
    <a class="card" href="DANTE_LTR/DANTE_LTR_summary.html">
      <b>LTR retrotransposon structural summary (DANTE_LTR)</b>
      <span>Intact LTR-RTs identified by protein-domain structure, by lineage.</span>
    </a>
    <a class="card" href="DANTE_TIR/report.html">
      <b>DNA transposon (TIR) structural summary (DANTE_TIR)</b>
      <span>Terminal-inverted-repeat elements (EnSpm/CACTA, hAT, PIF/Harbinger…):
      counts, TIR lengths and consensus logos.</span>
    </a>
    <a class="card" href="TideCluster/default/TideCluster_index.html">
      <b>Tandem repeat &amp; satellite report (TideCluster)</b>
      <span>Per-family (TRC) drill-down: monomer structure, dot-plots, TAREAN.</span>
    </a>
  </div>

  <div class="meta">
    Genome: ${GSTAT:-n/a}<br>
    Generated by CARP <b>${VERSION}</b>${RUN_DATE:+ on ${RUN_DATE}}.
  </div>

  <footer>
    Generated with the CARP pipeline. Source &amp; documentation:
    <a href="${CARP_URL}">the CARP repository on GitHub</a>.
    This page is a static example output; it is refreshed on demand and may lag
    the latest CARP release.
  </footer>
</body>
</html>
HTML

SIZE="$(du -sh "$OUT" | cut -f1)"
echo "Site assembled at $OUT  (total $SIZE)"

if [ "$DO_GIT" = 1 ]; then
  # fresh single-commit repo (no history accumulation); user force-pushes from host
  ( cd "$OUT"
    git init -q
    git symbolic-ref HEAD refs/heads/main   # portable across git versions (no `init -b`)
    printf '%s\n' '# CARP example output' \
      "Static GitHub Pages site — CARP ${VERSION} example on ${SPECIES} (${ACCESSION})." \
      'Rebuilt by .claude/skills/publish-example-report/build-example-site.sh; force-pushed.' > README.md
    # Pages must not run Jekyll (it would hide the _-prefixed/asset dirs)
    : > .nojekyll
    git add -A
    GIT_AUTHOR_NAME=CARP GIT_AUTHOR_EMAIL=noreply@example.invalid \
    GIT_COMMITTER_NAME=CARP GIT_COMMITTER_EMAIL=noreply@example.invalid \
      git commit -q -m "CARP ${VERSION} example — ${ACCESSION} (${RUN_DATE:-n/a})"
    if [ -n "${PREV_REMOTE:-}" ]; then git remote add origin "$PREV_REMOTE"; echo "restored origin: $PREV_REMOTE"; fi
  )
  echo "git repo initialized in $OUT (single commit)."
fi
echo "DONE"

#!/usr/bin/env bash
# Post-deploy hook for the tidecluster env. Snakemake runs this file
# automatically after `conda env create` if it exists alongside the
# matching <env>.yaml file.
#
# RepeatMasker (pulled in transitively via tidecluster → dante_ltr)
# ships its bundled reference libraries (is.lib, RepeatPeps.lib,
# Dfam pieces) as plain FASTA that must be BLAST-indexed before use.
# Bioconda's post-link.sh does this when the package is installed
# directly, but when snakemake creates the env through micromamba the
# post-link scripts can be silently skipped — leaving the `.nsq` /
# `.phr` companions missing and the first `repeatmasker` rule failing
# with:
#   RepeatMasker::createLib(): Error invoking makeblastdb on file ...
#
# We handle both is.lib (nucleotide) and RepeatPeps.lib (protein) plus
# any other *.lib files under Libraries/. Treat all errors as
# non-fatal: if a particular lib is already indexed, makeblastdb will
# say so and we continue.
set -u

conda_prefix="${CONDA_PREFIX:-}"
if [[ -z "$conda_prefix" ]]; then
  echo "tidecluster.post-deploy.sh: CONDA_PREFIX not set; skipping" >&2
  exit 0
fi

libs_dir="$conda_prefix/share/RepeatMasker/Libraries"
if [[ ! -d "$libs_dir" ]]; then
  echo "tidecluster.post-deploy.sh: no RepeatMasker/Libraries at $libs_dir; skipping" >&2
  exit 0
fi

makeblastdb="$conda_prefix/bin/makeblastdb"
if [[ ! -x "$makeblastdb" ]]; then
  echo "tidecluster.post-deploy.sh: makeblastdb missing; skipping" >&2
  exit 0
fi

index_nucl() {
  local lib="$1"
  # Already indexed? skip.
  [[ -e "${lib}.nsq" || -e "${lib}.nal" ]] && return 0
  "$makeblastdb" -in "$lib" -dbtype nucl -hash_index >/dev/null 2>&1 \
    && echo "  indexed $lib (nucl)" \
    || echo "  warn: failed to index $lib (nucl) — RepeatMasker may still work" >&2
}

index_prot() {
  local lib="$1"
  [[ -e "${lib}.psq" || -e "${lib}.pal" ]] && return 0
  "$makeblastdb" -in "$lib" -dbtype prot -hash_index >/dev/null 2>&1 \
    && echo "  indexed $lib (prot)" \
    || echo "  warn: failed to index $lib (prot) — RepeatMasker may still work" >&2
}

echo "tidecluster.post-deploy.sh: priming RepeatMasker library BLAST indices..."

# is.lib (insertion sequences) — required even when the RepeatMasker
# run passes -no_is; RepeatMasker::createLib still validates it.
[[ -f "$libs_dir/general/is.lib" ]] && index_nucl "$libs_dir/general/is.lib"

# Any other nucleotide libs under general/
while IFS= read -r -d '' lib; do
  [[ "$(basename "$lib")" == "is.lib" ]] && continue
  index_nucl "$lib"
done < <(find "$libs_dir/general" -maxdepth 1 -name '*.lib' -print0 2>/dev/null)

# Protein libs: RepeatPeps.lib and any *.phr-less *.lib at the top level
[[ -f "$libs_dir/RepeatPeps.lib" ]] && index_prot "$libs_dir/RepeatPeps.lib"

echo "tidecluster.post-deploy.sh: done"

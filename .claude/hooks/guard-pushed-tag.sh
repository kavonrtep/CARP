#!/usr/bin/env bash
# PreToolUse (Bash) guard: block creating / moving / deleting a git tag — or
# re-cutting a release — when that tag ALREADY EXISTS on the remote origin.
#
# Why: this sandbox has no ssh, so `git ls-remote` / `git push` can't reach
# origin and we must not ASSUME a local tag is unpushed. A published release tag
# is immutable; re-creating or moving it diverges from what CI/GHCR/Zenodo
# already built. This hook checks origin over the GitHub HTTPS API and denies the
# tool call if the target tag is already there.
#
# Complements (does not replace) the existing history-rewrite guard in
# .claude/settings.json, which only `ask`s on amend/rebase/reset/tag -d/-f/push
# --force and never consults the remote.
#
# Behaviour:
#   - tag already on origin  -> deny (permissionDecision: deny) with a reason
#   - remote can't be verified (offline / non-JSON) -> warn on stderr, ALLOW
#   - not a tag/release op, or tag not on origin -> allow (no output)
#
# Deps: jq, curl, git, sed, grep — all already used by this repo's tooling.
set -uo pipefail

cmd=$(jq -r '.tool_input.command // ""' 2>/dev/null) || exit 0
[ -n "$cmd" ] || exit 0

# ── Extract a version-shaped tag from a tag-affecting / release command ──────
VER='v?[0-9]+\.[0-9]+(\.[0-9]+)?(a|b|rc)?[0-9]*'   # PEP-440-lite-ish, matches this repo's tags
tag=""

# .claude/skills/release/cut-release.sh <VERSION>  (creates a tag)
if printf '%s' "$cmd" | grep -qE 'cut-release\.sh'; then
  tag=$(printf '%s' "$cmd" | grep -oE "cut-release\.sh[[:space:]]+[^[:space:];&|]+" | head -1 | awk '{print $2}')
fi

# git tag <name>            (create)
# git tag -a/-f/-d <name>   (annotate / move / delete)
if [ -z "$tag" ] && printf '%s' "$cmd" | grep -qE 'git[[:space:]]+tag([[:space:]]|$)'; then
  seg=$(printf '%s' "$cmd" | grep -oE 'git[[:space:]]+tag[^;&|]*' | head -1)
  # ignore pure listing (git tag -l / --list with no name still yields no VER match)
  tag=$(printf '%s' "$seg" | grep -oE "\b${VER}\b" | head -1)
fi

# git push ... --delete <tag> | --force | -f     (move/delete a pushed tag)
if [ -z "$tag" ] && printf '%s' "$cmd" | grep -qE 'git[[:space:]]+push([^;&|]*)(--delete|--force|--force-with-lease|[[:space:]]-f([[:space:]]|$))'; then
  seg=$(printf '%s' "$cmd" | grep -oE 'git[[:space:]]+push[^;&|]*' | head -1)
  tag=$(printf '%s' "$seg" | grep -oE "\b${VER}\b" | head -1)
fi

[ -n "$tag" ] || exit 0                 # not a tag/release op -> allow
tag=$(printf '%s' "$tag" | tr -d "\"'")

# ── Resolve owner/repo from origin (not hard-coded) ─────────────────────────
url=$(git remote get-url origin 2>/dev/null) || {
  echo "[tag-guard] no 'origin' remote; skipping remote-tag check for '$tag'" >&2; exit 0; }
slug=$(printf '%s' "$url" | sed -E 's#\.git$##' | sed -E 's#^git@[^:]+:##; s#^[a-zA-Z]+://[^/]+/##')
case "$slug" in
  */*) : ;;
  *) echo "[tag-guard] could not parse owner/repo from '$url'; skipping check for '$tag'" >&2; exit 0 ;;
esac

# ── Query origin's tags over HTTPS (no ssh). -L: the API 301-redirects. ─────
remote_tags=$(curl -sL --max-time 15 \
  -H "Accept: application/vnd.github+json" \
  "https://api.github.com/repos/${slug}/tags?per_page=100" 2>/dev/null \
  | jq -r '.[]?.name' 2>/dev/null)

if [ -z "$remote_tags" ]; then
  echo "[tag-guard] could not verify remote tags for '$slug' (offline / API error); ALLOWING '$tag' — verify manually that it is not already pushed before you push." >&2
  exit 0
fi

if printf '%s\n' "$remote_tags" | grep -qxF "$tag"; then
  reason="Tag '$tag' already exists on origin (already pushed/published) — a released tag is immutable and CI/GHCR/Zenodo already built it. Do NOT re-create, move, or delete it. Cut a NEW version instead (bump the number). If you believed this tag was unpushed, it is not: origin already has it."
  printf '{"hookSpecificOutput":{"hookEventName":"PreToolUse","permissionDecision":"deny","permissionDecisionReason":"%s"}}' "$reason"
fi
exit 0

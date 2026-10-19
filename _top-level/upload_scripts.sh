#!/usr/bin/env bash
# upload_scripts.sh - find scripts under a directory, copy them into a git repo
# (organized by project or language), and push to GitHub. Only new or changed
# scripts are copied/committed; duplicates are skipped.
#
# Usage:
#   ./upload_scripts.sh [-s SOURCE_DIR] [-r REPO_DIR] [-u REMOTE_URL]
#                       [-l project|language] [-x EXCLUDE_DIR]... [-p] [-D]
#                       [-n] [-f] [-m MAX_KB] [-b BRANCH]
#
#   -s  directory to search            (default: $HOME)
#   -r  local repo to build/update     (default: $HOME/my-scripts)
#   -u  GitHub remote URL              (e.g. git@github.com:you/my-scripts.git)
#   -l  layout: project -> <project>/<path>   language -> <lang>/<project>/<path>
#   -x  exclude a directory from the search (repeatable), e.g. -x ~/my-scripts.bak
#   -D  keep duplicates (default: files with identical content are copied once;
#       the first path in sorted order wins)
#   -p  prune: also delete files from the repo that are no longer selected
#       (deleted/renamed sources, duplicates, newly flagged files)
#   -n  dry run: show what would happen, change nothing
#   -f  also include files flagged as possibly containing secrets (not advised)
#   -m  skip files larger than this many KB (default 1024)
#   -b  branch name (default main)
#
# Always start with:  ./upload_scripts.sh -n -s ~ -r ~/SV_analysis

set -euo pipefail

SRC="$HOME"
REPO="$HOME/my-scripts"
REMOTE=""
LAYOUT="project"
DRY=0
INCLUDE_FLAGGED=0
MAXKB=1024
BRANCH="main"
DEDUPE=1
PRUNE=0
EXCLUDES=()

while getopts "s:r:u:l:x:Dpnfm:b:h" opt; do
  case $opt in
    s) SRC="$OPTARG" ;;
    r) REPO="$OPTARG" ;;
    u) REMOTE="$OPTARG" ;;
    l) LAYOUT="$OPTARG" ;;
    x) EXCLUDES+=("$OPTARG") ;;
    D) DEDUPE=0 ;;
    p) PRUNE=1 ;;
    n) DRY=1 ;;
    f) INCLUDE_FLAGGED=1 ;;
    m) MAXKB="$OPTARG" ;;
    b) BRANCH="$OPTARG" ;;
    h|*) sed -n '2,/^$/p' "$0"; exit 0 ;;
  esac
done

[[ "$LAYOUT" == "project" || "$LAYOUT" == "language" ]] || { echo "Layout must be project or language" >&2; exit 1; }
[[ -d "$SRC" ]] || { echo "Source dir not found: $SRC" >&2; exit 1; }
SRC="$(cd "$SRC" && pwd)"
mkdir -p "$REPO"
REPO="$(cd "$REPO" && pwd)"

# ---- extra excludes -> find prune arguments ---------------------------------
PRUNE_ARGS=()
for x in ${EXCLUDES[@]+"${EXCLUDES[@]}"}; do
  [[ -d "$x" ]] && x="$(cd "$x" && pwd)"
  PRUNE_ARGS+=(-o -path "$x")
done

# ---- sync: pull latest from GitHub before copying ---------------------------
# Only when the repo already has an origin (e.g. a clone). Failures (empty
# remote, no network, etc.) are ignored so the script can still continue.
if [[ $DRY -eq 0 && -d "$REPO/.git" ]] && git -C "$REPO" remote get-url origin >/dev/null 2>&1; then
  git -C "$REPO" pull --rebase -q origin "$BRANCH" || true
fi

# ---- extension -> language folder -------------------------------------------
lang_for() {
  local f base ext
  f="$1"; base="$(basename "$f")"
  [[ "$base" == "Snakefile" ]] && { echo "snakemake"; return; }
  ext="${base##*.}"; ext="${ext,,}"
  case "$ext" in
    py)                 echo "python" ;;
    sh|bash|zsh)        echo "shell" ;;
    r|rmd)              echo "r" ;;
    pl|pm)              echo "perl" ;;
    smk)                echo "snakemake" ;;
    nf)                 echo "nextflow" ;;
    pbs|slurm|sbatch)   echo "hpc-jobs" ;;
    ipynb)              echo "notebooks" ;;
    sql)                echo "sql" ;;
    js|ts)              echo "javascript" ;;
    *)                  echo "other" ;;
  esac
}

# ---- find candidate scripts -------------------------------------------------
# Skips hidden dirs, env/cache/output dirs, *.bak dirs, -x excludes, and the
# target repo itself.
mapfile -d '' FILES < <(
  find "$SRC" \
    \( -path "$REPO" ${PRUNE_ARGS[@]+"${PRUNE_ARGS[@]}"} \
       -o \( -type d ! -path "$SRC" -name '.*' \) \
       -o -type d \( -name node_modules -o -name __pycache__ -o -name .snakemake \
                     -o -name venv -o -name .venv -o -name env -o -name envs \
                     -o -name miniconda3 -o -name anaconda3 -o -name site-packages \
                     -o -name results -o -name work -o -name '*.bak' \) \
    \) -prune -o \
    -type f -size -"${MAXKB}k" \
    \( -name '*.py' -o -name '*.sh' -o -name '*.bash' -o -name '*.zsh' \
       -o -name '*.R' -o -name '*.r' -o -name '*.Rmd' -o -name '*.pl' -o -name '*.pm' \
       -o -name '*.smk' -o -name 'Snakefile' -o -name '*.nf' \
       -o -name '*.pbs' -o -name '*.slurm' -o -name '*.sbatch' \
       -o -name '*.ipynb' -o -name '*.sql' -o -name '*.js' -o -name '*.ts' \) \
    -print0 | sort -z
)

echo "Found ${#FILES[@]} candidate script(s) under $SRC"

# ---- very basic secret scan -------------------------------------------------
SECRET_RE='(AKIA[0-9A-Z]{16}|ghp_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{20,}|xox[baprs]-[A-Za-z0-9-]{10,}|-----BEGIN [A-Z ]*PRIVATE KEY-----|(password|passwd|secret|api[_-]?key|token)[[:space:]]*[=:][[:space:]]*["'"'"'][^"'"'"' ]{6,})'

declare -a MANIFEST=() FLAGGED=() DUPES=()
declare -A SEEN=()
COPIED=0
UNCHANGED=0
REMOVED=0

for f in "${FILES[@]}"; do
  rel="${f#"$SRC"/}"
  if [[ "$rel" == */* ]]; then project="${rel%%/*}"; inner="${rel#*/}"; else project="_top-level"; inner="$rel"; fi
  lang="$(lang_for "$f")"

  if grep -IqiE "$SECRET_RE" "$f" 2>/dev/null && [[ $INCLUDE_FLAGGED -eq 0 ]]; then
    FLAGGED+=("$rel"); continue
  fi

  # Skip files whose content is identical to one already selected (empty files excepted)
  if [[ $DEDUPE -eq 1 && -s "$f" ]]; then
    h="$(md5sum "$f" | cut -d' ' -f1)"
    if [[ -n "${SEEN[$h]:-}" ]]; then
      DUPES+=("$rel  (same as ${SEEN[$h]})"); continue
    fi
    SEEN[$h]="$rel"
  fi

  if [[ "$LAYOUT" == "language" ]]; then dest="$lang/$project/$inner"; else dest="$project/$inner"; fi
  MANIFEST+=("$project|$lang|$dest")

  status="new"
  if [[ -f "$REPO/$dest" ]]; then
    if cmp -s "$f" "$REPO/$dest"; then status="same"; else status="changed"; fi
  fi

  if [[ $DRY -eq 1 ]]; then
    if [[ "$status" == "same" ]]; then UNCHANGED=$((UNCHANGED + 1)); else
      echo "  [dry-run: $status] $rel  ->  $dest"; COPIED=$((COPIED + 1))
    fi
  else
    if [[ "$status" == "same" ]]; then UNCHANGED=$((UNCHANGED + 1)); continue; fi
    mkdir -p "$REPO/$(dirname "$dest")"
    cp -p "$f" "$REPO/$dest"
    COPIED=$((COPIED + 1))
  fi
done

if [[ ${#DUPES[@]} -gt 0 ]]; then
  echo
  echo "Skipped ${#DUPES[@]} duplicate file(s) (identical content to an earlier file):"
  printf '  - %s\n' "${DUPES[@]}"
fi

if [[ ${#FLAGGED[@]} -gt 0 ]]; then
  echo
  echo "Skipped ${#FLAGGED[@]} file(s) that look like they contain credentials (review, remove secrets, re-run):"
  printf '  - %s\n' "${FLAGGED[@]}"
fi

# ---- prune: remove repo files that are no longer selected -------------------
if [[ $PRUNE -eq 1 ]]; then
  declare -A KEEP=()
  for m in "${MANIFEST[@]}"; do KEEP["${m##*|}"]=1; done
  while IFS= read -r -d '' g; do
    r="${g#"$REPO"/}"
    [[ -n "${KEEP[$r]:-}" ]] && continue
    if [[ $DRY -eq 1 ]]; then
      echo "  [dry-run: remove from repo] $r"
    else
      rm -f "$g"
    fi
    REMOVED=$((REMOVED + 1))
  done < <(find "$REPO" -path "$REPO/.git" -prune -o -type f \
             ! -path "$REPO/README.md" ! -path "$REPO/.gitignore" -print0)
  if [[ $DRY -eq 0 ]]; then
    find "$REPO" -mindepth 1 -type d -empty -not -path "$REPO/.git" -not -path "$REPO/.git/*" -delete
  fi
fi

if [[ $DRY -eq 1 ]]; then
  echo
  echo "Dry run complete: $COPIED new/changed file(s) would be copied, $UNCHANGED already up to date, ${#DUPES[@]} duplicate(s) skipped, $REMOVED would be pruned. Repo: $REPO"
  echo "(Compared against the local repo as it is now; run 'git -C $REPO pull' first to compare against GitHub.)"
  exit 0
fi

# ---- README index + .gitignore ---------------------------------------------
if [[ ! -f "$REPO/.gitignore" ]]; then
  cat > "$REPO/.gitignore" <<'EOF'
__pycache__/
*.pyc
.snakemake/
.ipynb_checkpoints/
*.bam
*.bai
*.vcf.gz
*.fastq*
*.fq*
.env
*.log
EOF
fi

{
  echo "# My scripts"
  echo
  echo "Auto-organized collection of scripts (layout: \`$LAYOUT\`)."
  echo
  echo "| Project | Language | File |"
  echo "|---|---|---|"
  printf '%s\n' "${MANIFEST[@]}" | sort | while IFS='|' read -r p l d; do
    echo "| $p | $l | [\`$d\`](./$(printf '%s' "$d" | sed 's/ /%20/g')) |"
  done
} > "$REPO/README.md"

# ---- git --------------------------------------------------------------------
cd "$REPO"
if [[ ! -d .git ]]; then
  git init -q -b "$BRANCH" 2>/dev/null || { git init -q; git checkout -q -b "$BRANCH"; }
fi

git add -A
if git diff --cached --quiet; then
  echo "No changes to commit."
else
  git commit -q -m "Sync scripts ($(date +%F)): $COPIED new/changed, $REMOVED removed"
  echo "Committed: $COPIED new/changed file(s), $REMOVED removed."
fi

if [[ -n "$REMOTE" ]]; then
  if git remote get-url origin >/dev/null 2>&1; then git remote set-url origin "$REMOTE"; else git remote add origin "$REMOTE"; fi
fi

if git remote get-url origin >/dev/null 2>&1; then
  git push -u origin "$BRANCH"
  echo "Pushed to $(git remote get-url origin)"
else
  echo
  echo "No remote configured. Create an empty repo on GitHub, then either re-run with -u <url>"
  echo "or, if you have the GitHub CLI:  cd $REPO && gh repo create my-scripts --private --source=. --push"
fi

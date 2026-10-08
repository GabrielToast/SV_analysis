#!/usr/bin/env bash
# upload_scripts.sh - find scripts under a directory, copy them into a clean
# git repo (organized by project or language), and push to GitHub.
#
# Usage:
#   ./upload_scripts.sh [-s SOURCE_DIR] [-r REPO_DIR] [-u REMOTE_URL]
#                       [-l project|language] [-n] [-f] [-m MAX_KB] [-b BRANCH]
#
#   -s  directory to search            (default: $HOME)
#   -r  local repo to build/update     (default: $HOME/my-scripts)
#   -u  GitHub remote URL              (e.g. git@github.com:you/my-scripts.git)
#   -l  layout: project -> <project>/<path>   language -> <lang>/<project>/<path>
#   -n  dry run: show what would happen, change nothing
#   -f  also include files flagged as possibly containing secrets (not advised)
#   -m  skip files larger than this many KB (default 1024)
#   -b  branch name (default main)
#
# Always start with:  ./upload_scripts.sh -n -s ~/projects

set -euo pipefail

SRC="$HOME"
REPO="$HOME/my-scripts"
REMOTE=""
LAYOUT="project"
DRY=0
INCLUDE_FLAGGED=0
MAXKB=1024
BRANCH="main"

while getopts "s:r:u:l:nfm:b:h" opt; do
  case $opt in
    s) SRC="$OPTARG" ;;
    r) REPO="$OPTARG" ;;
    u) REMOTE="$OPTARG" ;;
    l) LAYOUT="$OPTARG" ;;
    n) DRY=1 ;;
    f) INCLUDE_FLAGGED=1 ;;
    m) MAXKB="$OPTARG" ;;
    b) BRANCH="$OPTARG" ;;
    h|*) sed -n '2,20p' "$0"; exit 0 ;;
  esac
done

[[ "$LAYOUT" == "project" || "$LAYOUT" == "language" ]] || { echo "Layout must be project or language" >&2; exit 1; }
[[ -d "$SRC" ]] || { echo "Source dir not found: $SRC" >&2; exit 1; }
SRC="$(cd "$SRC" && pwd)"
mkdir -p "$REPO"
REPO="$(cd "$REPO" && pwd)"

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
# Skips hidden dirs, env/cache/output dirs, and the target repo itself.
mapfile -d '' FILES < <(
  find "$SRC" \
    \( -path "$REPO" \
       -o \( -type d ! -path "$SRC" -name '.*' \) \
       -o -type d \( -name node_modules -o -name __pycache__ -o -name .snakemake \
                     -o -name venv -o -name .venv -o -name env -o -name envs \
                     -o -name miniconda3 -o -name anaconda3 -o -name site-packages \
                     -o -name results -o -name work \) \
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

declare -a MANIFEST=() FLAGGED=()
COPIED=0

for f in "${FILES[@]}"; do
  rel="${f#"$SRC"/}"
  if [[ "$rel" == */* ]]; then project="${rel%%/*}"; inner="${rel#*/}"; else project="_top-level"; inner="$rel"; fi
  lang="$(lang_for "$f")"

  if grep -IqiE "$SECRET_RE" "$f" 2>/dev/null && [[ $INCLUDE_FLAGGED -eq 0 ]]; then
    FLAGGED+=("$rel"); continue
  fi

  if [[ "$LAYOUT" == "language" ]]; then dest="$lang/$project/$inner"; else dest="$project/$inner"; fi
  MANIFEST+=("$project|$lang|$dest")

  if [[ $DRY -eq 1 ]]; then
    echo "  [dry-run] $rel  ->  $dest"
  else
    mkdir -p "$REPO/$(dirname "$dest")"
    cp -p "$f" "$REPO/$dest"
    COPIED=$((COPIED + 1))
  fi
done

if [[ ${#FLAGGED[@]} -gt 0 ]]; then
  echo
  echo "Skipped ${#FLAGGED[@]} file(s) that look like they contain credentials (review, remove secrets, re-run):"
  printf '  - %s\n' "${FLAGGED[@]}"
fi

if [[ $DRY -eq 1 ]]; then
  echo; echo "Dry run complete: ${#MANIFEST[@]} file(s) would be copied to $REPO"; exit 0
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
  echo "Auto-organized collection of scripts (layout: \`$LAYOUT\`). Generated $(date +%F)."
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
  git commit -q -m "Sync scripts ($(date +%F)): $COPIED file(s)"
  echo "Committed $COPIED file(s)."
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
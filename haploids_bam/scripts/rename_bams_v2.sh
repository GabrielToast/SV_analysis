#!/usr/bin/env bash
# Rename BAM files to their NGI ID using an exact old-filename -> NGI-ID table.
#
# Usage: ./rename_bams_v2.sh new_old_name_bam.txt BAM_DIR [--apply]
#   Mapping file (tab-separated):  NGI_ID <TAB> old_filename.bam
#   Lines without a second column (e.g. the c_80_SRR... rows) are ignored.
#   Without --apply this is a DRY RUN: nothing is renamed, only printed.
#
# Example: AK6:1.sorted.marked_dupes.rg.bam -> P11604_839_sorted_marked_dupes_rg.bam
# Rule:    NGI ID + "_" + rest of the old name with every "." turned into "_".
# Existing .bai files (name.bam.bai or name.bai) are renamed along with their BAM.

set -euo pipefail

MAP="${1:?usage: $0 mapping.txt BAM_DIR [--apply]}"
DIR="${2:?usage: $0 mapping.txt BAM_DIR [--apply]}"
APPLY="${3:-}"

declare -A NGI_OF        # old filename -> NGI ID
declare -A OLD_OF        # NGI ID -> old filename (to catch duplicated IDs)
n_lines=0
while IFS=$'\t' read -r ngi old _; do
    [[ -z "${ngi:-}" || -z "${old:-}" ]] && continue      # skips 1-column rows
    n_lines=$((n_lines + 1))
    if [[ -n "${NGI_OF[$old]:-}" && "${NGI_OF[$old]}" != "$ngi" ]]; then
        echo "WARNING: $old listed with two IDs: ${NGI_OF[$old]} and $ngi" >&2
    fi
    if [[ -n "${OLD_OF[$ngi]:-}" && "${OLD_OF[$ngi]}" != "$old" ]]; then
        echo "WARNING: $ngi assigned to two files: ${OLD_OF[$ngi]} and $old" >&2
    fi
    NGI_OF["$old"]="$ngi"
    OLD_OF["$ngi"]="$old"
done < <(tr -d '\r' < "$MAP")                              # strip Windows CRs

echo "Loaded $n_lines mappings"

declare -A TARGETS
LOG="$DIR/rename_log_$(date +%Y%m%d_%H%M%S).tsv"
n_ok=0; n_missing=0; n_conflict=0

for f in "$DIR"/*.bam; do
    [[ -e "$f" ]] || { echo "No .bam files in $DIR"; exit 1; }
    base=$(basename "$f")

    if [[ -z "${NGI_OF[$base]:-}" ]]; then
        echo "NO MATCH   $base (not in mapping, left untouched)" >&2
        n_missing=$((n_missing + 1)); continue
    fi

    rest="${base#*.}"; rest="${rest%.bam}"; rest="${rest//./_}"
    new="${NGI_OF[$base]}_${rest}.bam"

    if [[ -n "${TARGETS[$new]:-}" || -e "$DIR/$new" ]]; then
        echo "CONFLICT   $base -> $new (target already exists or is duplicated)" >&2
        n_conflict=$((n_conflict + 1)); continue
    fi
    TARGETS["$new"]=1

    echo "$base -> $new"
    n_ok=$((n_ok + 1))

    if [[ "$APPLY" == "--apply" ]]; then
        mv -n -- "$f" "$DIR/$new"
        printf '%s\t%s\n' "$base" "$new" >> "$LOG"
        [[ -e "$f.bai" ]]        && mv -n -- "$f.bai"        "$DIR/$new.bai"
        [[ -e "${f%.bam}.bai" ]] && mv -n -- "${f%.bam}.bai" "$DIR/${new%.bam}.bai"
    fi
done

# Dry run only: mapping entries whose BAM is not in the directory
if [[ "$APPLY" != "--apply" ]]; then
    for old in "${!NGI_OF[@]}"; do
        [[ -e "$DIR/$old" ]] || echo "MAPPED BUT NO FILE   $old (${NGI_OF[$old]})" >&2
    done
fi

echo
echo "Renamable: $n_ok | No match: $n_missing | Conflicts: $n_conflict"
if [[ "$APPLY" == "--apply" ]]; then
    echo "Renamed. Log of old -> new names: $LOG"
else
    echo "DRY RUN only. Re-run with --apply to rename."
fi

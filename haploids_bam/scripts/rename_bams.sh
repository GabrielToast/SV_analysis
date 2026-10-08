#!/usr/bin/env bash
# Rename BAM files from user ID to NGI ID using a 2-column TSV mapping.
#
# Usage: ./rename_bams.sh mapping.tsv BAM_DIR [--apply]
#   mapping.tsv : column 1 = user ID (e.g. AK6:1), column 2 = NGI ID (e.g. P11604_839)
#   Without --apply this is a DRY RUN: nothing is renamed, only printed.
#
# Example:  AK6:1.sorted.marked_dupes.rg.bam -> P11604_839_sorted_marked_dupes_rg.bam
# Rule:     sample = text before the first "."; the remaining dots become "_".
# Any existing .bai files (name.bam.bai or name.bai) are renamed too.

set -euo pipefail

MAP="${1:?usage: $0 mapping.tsv BAM_DIR [--apply]}"
DIR="${2:?usage: $0 mapping.tsv BAM_DIR [--apply]}"
APPLY="${3:-}"

# norm: build a comparable key from an ID. The text after the LAST separator
# (-, _ or :) is the replicate number; everything before it (the prefix) is
# reduced to lowercase letters/digits only. Result goes in $NORM as "prefix|rep".
#   AK6:1 / AK_6_1            -> ak6|1
#   VSH-014:7 / VSH_014_7     -> vsh014|7
#   RIC-SG-NL3:4 / RICSG_NL3_4 -> ricsgnl3|4
# Keeping the replicate separate avoids false matches such as C1:11 vs C11:1.
# Applied to both mapping IDs and filename IDs.
norm() {
    local s="$1" pre rep
    if [[ $s =~ ^(.*)[-_:]([^-_:]*)$ ]]; then
        pre="${BASH_REMATCH[1]}"; rep="${BASH_REMATCH[2]}"
    else
        pre="$s"; rep=""
    fi
    pre="${pre//[^[:alnum:]]/}"
    NORM="${pre,,}|${rep,,}"
}

declare -A NGI
while IFS=$'\t' read -r user ngi _; do
    [[ -z "${user:-}" || -z "${ngi:-}" ]] && continue
    norm "$user"; key="$NORM"
    if [[ -n "${NGI[$key]:-}" && "${NGI[$key]}" != "$ngi" ]]; then
        echo "WARNING: '$user' normalizes to '$key' which maps to both ${NGI[$key]} and $ngi" >&2
    fi
    NGI["$key"]="$ngi"
done < <(tr -d '\r' < "$MAP")      # tr strips Windows line endings

echo "Loaded ${#NGI[@]} ID mappings"

declare -A TARGETS
LOG="$DIR/rename_log_$(date +%Y%m%d_%H%M%S).tsv"
n_ok=0; n_missing=0; n_conflict=0

for f in "$DIR"/*.bam; do
    [[ -e "$f" ]] || { echo "No .bam files in $DIR"; exit 1; }
    base=$(basename "$f")
    sample="${base%%.*}"                 # AK6:1
    rest="${base#*.}"                    # sorted.marked_dupes.rg.bam
    rest="${rest%.bam}"                  # sorted.marked_dupes.rg
    rest="${rest//./_}"                  # sorted_marked_dupes_rg

    norm "$sample"; key="$NORM"

    if [[ -z "${NGI[$key]:-}" ]]; then
        echo "NO MATCH   $base  (looked up '$key' in mapping)" >&2
        n_missing=$((n_missing + 1)); continue
    fi

    new="${NGI[$key]}_${rest}.bam"

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

echo
echo "Renamable: $n_ok | No match: $n_missing | Conflicts: $n_conflict"
if [[ "$APPLY" == "--apply" ]]; then
    echo "Renamed. Log of old -> new names: $LOG"
else
    echo "DRY RUN only. Re-run with --apply to rename."
fi

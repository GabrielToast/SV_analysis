#!/bin/bash
set -u

# ---------- settings ----------
BAMDIR=/home/gato/haploids_bam/bam
SIBFILE=siblings.txt
OUTDIR=samplot_test
SAMPLE=P27152_182
CHROM=NC_037639.1
START=14695130
END=14695963
TYPE=DEL
WINDOW=500
# ------------------------------

mkdir -p "$OUTDIR"

find_bam() { ls "$BAMDIR"/"$1"*.bam 2>/dev/null | head -n 1; }

names=("$SAMPLE")
bams=("$(find_bam "$SAMPLE")")

# every other ID in siblings.txt becomes a control
while read -r c; do
  [ -z "$c" ] && continue
  [ "$c" = "$SAMPLE" ] && continue
  cbam=$(find_bam "$c")
  if [ -z "$cbam" ]; then echo "No BAM for $c, skipping" >&2; continue; fi
  names+=("$c")
  bams+=("$cbam")
done < "$SIBFILE"

echo "Plotting ${#names[@]} samples: ${names[*]}"

samplot plot \
  -n "${names[@]}" \
  -b "${bams[@]}" \
  -o "$OUTDIR/${SAMPLE}_${TYPE}_${CHROM}_${START}_${END}_with_sibs.png" \
  -c "$CHROM" -s "$START" -e "$END" \
  -t "$TYPE" \
  -w "$WINDOW" \
  -H 3 -W 10 \
  --same_yaxis_scales

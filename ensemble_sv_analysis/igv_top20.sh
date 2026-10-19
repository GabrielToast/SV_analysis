#!/bin/bash
set -u

# ---------- settings ----------
REF=/path/to/Amel_HAv3.1.fna          # EDIT: reference FASTA (needs .fai)
BAMDIR=/home/gato/haploids_bam/bam
OUTDIR=samplot_top20                  # same folder as plot_top20.sh
IGV=igv
PAD=40                                # bp shown each side of each breakpoint
# ------------------------------

DT="$OUTDIR/depth_table.tsv"
SNAPDIR="$(realpath "$OUTDIR")/igv"
BATCH="$SNAPDIR/igv_batch.txt"
mkdir -p "$SNAPDIR"
[ -f "$DT" ] || { echo "ERROR: $DT not found, run plot_top20.sh first" >&2; exit 1; }
[ -f "$REF.fai" ] || [ -f "${REF%.*}.fai" ] || echo "WARNING: no .fai index for $REF" >&2

{
  echo "preference SAM.SHOW_SOFT_CLIPPED TRUE"
  echo "preference SAM.DOWNSAMPLE_READS FALSE"
  echo "snapshotDirectory $SNAPDIR"
} > "$BATCH"

for cand in $(tail -n +2 "$DT" | cut -f2 | sort -u); do
  # chrom, start, end are the same on every row of a candidate
  read -r chrom start end < <(awk -F'\t' -v c="$cand" '$2==c{print $4, $5, $6; exit}' "$DT")

  echo "new" >> "$BATCH"
  echo "genome $REF" >> "$BATCH"

  # focal first (rows are already written focal-first), then controls
  n=0
  for id in $(awk -F'\t' -v c="$cand" '$2==c{print $8}' "$DT"); do
    bam=$(ls "$BAMDIR"/${id}_*.bam 2>/dev/null | head -n 1)
    [ -z "$bam" ] && { echo "no BAM for $id" >&2; continue; }
    echo "load $(realpath "$bam")" >> "$BATCH"
    n=$((n+1))
  done

  # two panels side by side: left breakpoint | right breakpoint
  l1=$((start-PAD)); [ $l1 -lt 1 ] && l1=1
  l2=$((start+PAD)); r1=$((end-PAD)); r2=$((end+PAD))
  echo "goto $chrom:$l1-$l2 $chrom:$r1-$r2" >> "$BATCH"
  echo "expand" >> "$BATCH"
  echo "snapshot ${cand}_sequence.png" >> "$BATCH"
  echo "[$cand] $n tracks, breakpoints $chrom:$start and $chrom:$end"
done
echo "exit" >> "$BATCH"

echo "Running IGV ..."
if command -v xvfb-run >/dev/null; then
  xvfb-run -a "$IGV" -b "$BATCH"
else
  "$IGV" -b "$BATCH"
fi
echo "Sequence views: $SNAPDIR/*_sequence.png"

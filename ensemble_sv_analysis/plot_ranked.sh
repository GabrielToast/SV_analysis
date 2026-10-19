#!/bin/bash
set -u

# ---------- settings ----------
TSV=ranked_to_plot.tsv
ALL=ranked_all.tsv      # used only to look up each call's score
NTOP=200                # first NTOP rows are top-ranked, the rest are random
BAMDIR=/home/gato/haploids_bam/bam
MAPFILE=colonies.tsv          # NGI_ID<TAB>Colony_ID, no header
OUTDIR=samplot_ranked
TOPN=all
WINDOW=500                    # samplot flanking window
FLANK=1000                    # flank on each side for baseline depth
CARRIER_RATIO=0.1             # DEL carrier if interval/flank depth < this
DUP_RATIO=1.5                 # DUP carrier if ratio > this
MIN_FLANK_DEPTH=5             # below this the sample is "lowcov" (no call)
ALLOW_PREFIX_FALLBACK=0       # 1 = use ID prefix as colony if no map file
# ------------------------------

mkdir -p "$OUTDIR"
declare -A BAMOF COLONY MEMBERS OK

# 1. index BAM paths by sample ID (first two underscore fields of filename)
for f in "$BAMDIR"/*.bam; do
  b=${f##*/}
  rest=${b#*_}
  id="${b%%_*}_${rest%%_*}"
  BAMOF[$id]=$f
done
echo "BAM files found: ${#BAMOF[@]}"

# 2. colony assignment
if [ -f "$MAPFILE" ]; then
  while IFS=$'\t' read -r id col || [ -n "$id" ]; do
    id=${id//$'\r'/}; col=${col//$'\r'/}
    id=${id// /}; col=${col// /}
    [ -z "$id" ] && continue
    COLONY[$id]=$col
  done < "$MAPFILE"
  echo "Colony map: $MAPFILE (${#COLONY[@]} entries)"
elif [ "$ALLOW_PREFIX_FALLBACK" = 1 ]; then
  echo "No $MAPFILE: using ID prefix as colony" >&2
  for id in "${!BAMOF[@]}"; do COLONY[$id]=${id%%_*}; done
else
  echo "ERROR: $MAPFILE not found. Create it first (make_colonies.py)." >&2
  exit 1
fi

nomap=0
for id in "${!BAMOF[@]}"; do
  col=${COLONY[$id]:-}
  if [ -n "$col" ]; then MEMBERS[$col]+="$id "; else nomap=$((nomap+1)); fi
done
[ $nomap -gt 0 ] && echo "WARNING: $nomap BAMs have no colony entry" >&2

# lazy BAM validation with cache
usable() {
  local id=$1 f=${BAMOF[$1]:-}
  [ -z "$f" ] && return 1
  if [ -z "${OK[$id]:-}" ]; then
    if samtools quickcheck "$f" 2>/dev/null \
       && [ "$(samtools view -H "$f" 2>/dev/null | grep -c '^@SQ')" -gt 0 ] \
       && { [ -f "$f.bai" ] || [ -f "${f%.bam}.bai" ]; }; then
      OK[$id]=1
    else
      OK[$id]=0
      echo "UNUSABLE BAM (quickcheck/@SQ/index): $f" >&2
    fi
  fi
  [ "${OK[$id]}" = 1 ]
}

mean_depth() {   # bam region [region ...]
  local bam=$1; shift
  local r
  for r in "$@"; do samtools depth -a -r "$r" "$bam"; done |
    awk '{s+=$3;n++} END{if(n>0) printf "%.2f", s/n; else print "NA"}'
}

# 3. top N candidates (score desc, qual desc), header skipped
cp "$TSV" top${TOPN}.tsv

DT="$OUTDIR/depth_table.tsv"
SM="$OUTDIR/summary.tsv"
printf "rank\tcandidate\tcolony\tchrom\tstart\tend\ttype\tsample\trole\tinterval_depth\tflank_depth\tratio\tcarrier\n" > "$DT"
printf "rank\tcandidate\tcolony\tfocal_carrier\tn_controls\tn_control_carriers\tcarrier_ids\tverdict\n" > "$SM"

n=0
while IFS=$'\t' read -r -u 3 chrom start end type svlen sample gt; do
  n=$((n+1)); rank=$(printf "%03d" "$n")
  score=$(awk -F'\t' -v c="$chrom" -v s="$start" -v e="$end" -v id="$sample" '$1==c && $2==s && $3==e && $6==id {print $15; exit}' "$ALL")
  score=${score:-NA}
  grp=top; [ "$n" -gt "$NTOP" ] && grp=rnd
  cand="${rank}_${grp}_s${score}_${sample}_${type}_${chrom}_${start}_${end}"

  if ! usable "$sample"; then echo "[$rank] no usable BAM for $sample, skipping" >&2; continue; fi
  col=${COLONY[$sample]:-}
  if [ -z "$col" ]; then echo "[$rank] no colony for $sample, skipping" >&2; continue; fi

  # focal first, then usable siblings (sorted)
  names=("$sample"); bams=("${BAMOF[$sample]}")
  for id in $(echo ${MEMBERS[$col]} | tr ' ' '\n' | sort); do
    [ "$id" = "$sample" ] && continue
    usable "$id" || continue
    names+=("$id"); bams+=("${BAMOF[$id]}")
  done
  echo "[$rank] $cand colony=$col samples=${#names[@]} (${names[*]})"

  # depth table
  lstart=$((start-FLANK)); [ $lstart -lt 1 ] && lstart=1
  lend=$((start-1)); rstart=$((end+1)); rend=$((end+FLANK))
  focal_carrier="NA"; nctrl=0; nctrl_car=0; carriers=""
  for i in "${!names[@]}"; do
    id=${names[$i]}; bam=${bams[$i]}
    role=control; [ $i -eq 0 ] && role=focal
    idepth=$(mean_depth "$bam" "$chrom:$start-$end")
    fdepth=$(mean_depth "$bam" "$chrom:$lstart-$lend" "$chrom:$rstart-$rend")
    read -r ratio carrier < <(awk -v i="$idepth" -v f="$fdepth" -v t="$type" \
      -v cr="$CARRIER_RATIO" -v dr="$DUP_RATIO" -v mf="$MIN_FLANK_DEPTH" 'BEGIN{
        if(i=="NA"||f=="NA"||f+0<mf+0){print "NA lowcov"; exit}
        r=i/f
        if(t=="DUP"){c=(r>dr)?"yes":"no"} else {c=(r<cr)?"yes":"no"}
        printf "%.3f %s\n", r, c}')
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$rank" "$cand" "$col" "$chrom" "$start" "$end" "$type" "$id" "$role" "$idepth" "$fdepth" "$ratio" "$carrier" >> "$DT"
    if [ "$role" = focal ]; then
      focal_carrier=$carrier
    else
      nctrl=$((nctrl+1))
      if [ "$carrier" = yes ]; then nctrl_car=$((nctrl_car+1)); carriers+="$id,"; fi
    fi
  done

  if [ "$focal_carrier" != yes ]; then verdict="FOCAL_NOT_SUPPORTED_BY_DEPTH"
  elif [ $nctrl_car -eq 0 ]; then verdict="DENOVO_SUPPORTED_BY_DEPTH"
  else verdict="SHARED_WITH_SIBLINGS"; fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$rank" "$cand" "$col" "$focal_carrier" "$nctrl" "$nctrl_car" "${carriers%,}" "$verdict" >> "$SM"

  # plot
  samplot plot \
    -n "${names[@]}" -b "${bams[@]}" \
    -o "$OUTDIR/${cand}.png" \
    -c "$chrom" -s "$start" -e "$end" -t "$type" \
    -w "$WINDOW" -H 3 -W 10 --same_yaxis_scales
done 3< top${TOPN}.tsv

echo
echo "Plots:       $OUTDIR/*.png"
echo "Depth table: $DT"
echo "Summary:     $SM"
column -t -s$'\t' "$SM"

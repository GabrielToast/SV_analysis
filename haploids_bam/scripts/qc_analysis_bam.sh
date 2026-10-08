#!/usr/bin/env bash
# QC for all .bam files in a directory: mosdepth + samtools stats/flagstat + MultiQC.
#
# Usage: ./qc_bams.sh [BAM_DIR] [OUT_DIR] [TOTAL_CORES] [THREADS_PER_JOB]
#   BAM_DIR          directory with .bam files          (default: all_bams)
#   OUT_DIR          where results go                   (default: qc)
#   TOTAL_CORES      cores to use in total              (default: 32)
#   THREADS_PER_JOB  threads per sample                 (default: 4)
#
# Concurrent samples = TOTAL_CORES / THREADS_PER_JOB (default 8 x 4).
# Per sample it runs:
#   mosdepth  -Q 20, 10 kb windows, no per-base output  -> OUT_DIR/mosdepth/
#   samtools stats + flagstat                           -> OUT_DIR/samtools/
# Then MultiQC aggregates everything                    -> OUT_DIR/multiqc/
# Finished samples are skipped, so you can safely rerun after an interruption.
# BAMs without an index are indexed automatically.

set -euo pipefail

BAM_DIR="${1:-all_bams}"
OUT="${2:-qc}"
CORES="${3:-32}"
THREADS="${4:-4}"
MAPQ=20
WINDOW=10000

JOBS=$(( CORES / THREADS ))
(( JOBS < 1 )) && JOBS=1

for tool in samtools mosdepth; do
    command -v "$tool" >/dev/null || { echo "$tool not found in PATH" >&2; exit 1; }
done

mkdir -p "$OUT"/{mosdepth,samtools,logs}

qc_one() {
    local bam="$1" threads="$2" out="$3" mapq="$4" window="$5"
    local name; name=$(basename "$bam" .bam)
    local log="$out/logs/$name.log"

    # Done already? (mosdepth summary and samtools stats both present)
    if [[ -s "$out/mosdepth/$name.mosdepth.summary.txt" && -s "$out/samtools/$name.flagstat" ]]; then
        echo "SKIP  $name"
        return 0
    fi

    {
        # Index if needed (mosdepth requires it)
        if [[ ! -e "$bam.bai" && ! -e "${bam%.bam}.bai" ]]; then
            samtools index -@ "$threads" "$bam"
        fi

        mosdepth -t "$threads" -Q "$mapq" --no-per-base --by "$window" \
            "$out/mosdepth/$name" "$bam"

        samtools stats -@ "$threads" "$bam" > "$out/samtools/$name.stats"
        samtools flagstat -@ "$threads" "$bam" > "$out/samtools/$name.flagstat"
    } > "$log" 2>&1 || {
        echo "FAIL  $name (see $log)" >&2
        # remove partial outputs so a rerun redoes this sample
        rm -f "$out/mosdepth/$name".* "$out/samtools/$name".stats "$out/samtools/$name".flagstat
        return 1
    }
    echo "OK    $name"
}
export -f qc_one

n=$(find "$BAM_DIR" -maxdepth 1 \( -type f -o -type l \) -name '*.bam' | wc -l)
if (( n == 0 )); then
    echo "No .bam files found directly in '$BAM_DIR' (symlinks are included, subfolders are not)." >&2
    echo "Check the path, and look in subfolders with: find '$BAM_DIR' -name '*.bam' | head" >&2
    exit 1
fi
echo "QC on $n BAMs: $JOBS parallel samples x $THREADS threads (MAPQ >= $MAPQ, ${WINDOW} bp windows)"

set +e
find "$BAM_DIR" -maxdepth 1 \( -type f -o -type l \) -name '*.bam' -print0 \
    | xargs -r -0 -n 1 -P "$JOBS" bash -c 'qc_one "$1" '"$THREADS $OUT $MAPQ $WINDOW" _
xargs_status=$?
set -e

# Aggregate report
if command -v multiqc >/dev/null; then
    multiqc "$OUT/mosdepth" "$OUT/samtools" -o "$OUT/multiqc" --force
else
    echo "multiqc not found in PATH; skipping the aggregate report." >&2
    echo "Install it with: pip install multiqc  (or: conda install -c bioconda multiqc)" >&2
fi

# Simple cohort table: mean depth (total row) and % mapped per sample
{
    printf 'sample\tmean_depth\tpct_mapped\n'
    for s in "$OUT"/mosdepth/*.mosdepth.summary.txt; do
        [[ -e "$s" ]] || continue
        name=$(basename "$s" .mosdepth.summary.txt)
        depth=$(awk '$1=="total"{print $4}' "$s")
        mapped=$(grep -m1 ' mapped (' "$OUT/samtools/$name.flagstat" 2>/dev/null \
                 | sed -E 's/.*\(([0-9.]+)%.*/\1/' || true)
        printf '%s\t%s\t%s\n' "$name" "$depth" "${mapped:-NA}"
    done
} > "$OUT/qc_summary.tsv"

echo
echo "Summary table: $OUT/qc_summary.tsv"
[[ -e "$OUT/multiqc/multiqc_report.html" ]] && echo "MultiQC report: $OUT/multiqc/multiqc_report.html"
if (( xargs_status != 0 )); then
    echo "Some samples failed; check messages above and $OUT/logs/. Rerun to retry only those." >&2
    exit 1
fi
echo "Done."
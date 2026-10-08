#!/usr/bin/env bash
# Index all .bam files in a directory in parallel using samtools.
#
# Usage: ./index_bams.sh [BAM_DIR] [TOTAL_CORES] [THREADS_PER_JOB]
#   BAM_DIR          directory with .bam files      (default: .)
#   TOTAL_CORES      cores to use in total          (default: 32)
#   THREADS_PER_JOB  samtools -@ threads per file   (default: 1)
#
# Concurrent jobs = TOTAL_CORES / THREADS_PER_JOB.
# With ~200 files, 32 jobs x 1 thread keeps all cores busy and avoids
# idle cores at the tail. If your storage is slow, try 16 jobs x 2 threads.

set -euo pipefail

BAM_DIR="${1:-.}"
CORES="${2:-32}"
THREADS="${3:-1}"
JOBS=$(( CORES / THREADS ))
(( JOBS < 1 )) && JOBS=1

command -v samtools >/dev/null || { echo "samtools not found in PATH" >&2; exit 1; }

index_one() {
    local bam="$1" threads="$2"
    local bai="${bam}.bai"

    # Skip if an index exists and is newer than the BAM
    if [[ -s "$bai" && "$bai" -nt "$bam" ]]; then
        echo "SKIP  $bam"
        return 0
    fi

    if samtools index -@ "$threads" "$bam" "$bai" 2> "${bai}.err"; then
        rm -f "${bai}.err"
        echo "OK    $bam"
    else
        echo "FAIL  $bam (see ${bai}.err)" >&2
        rm -f "$bai"
        return 1
    fi
}
export -f index_one

n=$(find "$BAM_DIR" -maxdepth 1 -type f -name '*.bam' | wc -l)
echo "Indexing $n BAM files: $JOBS parallel jobs x $THREADS thread(s) each"

# -print0/-0 handles odd filenames; xargs exits non-zero if any job failed
find "$BAM_DIR" -maxdepth 1 -type f -name '*.bam' -print0 \
    | xargs -0 -n 1 -P "$JOBS" bash -c 'index_one "$1" '"$THREADS" _

echo "Done."

#!/bin/bash
set -uo pipefail

BAM_DIR="all_bams"
OUTDIR="mosdepth_results_all_bams"
TOTAL_CORES=32
THREADS_PER_JOB=4
JOBS=$(( TOTAL_CORES / THREADS_PER_JOB ))   # 8 samples at once

mkdir -p "$OUTDIR"

run_mosdepth() {
    bam="$1"
    sample=$(basename "$bam" .bam)
    summary="$OUTDIR/$sample.mosdepth.summary.txt"

    # Skip samples that already finished
    if [[ ! -s "$summary" ]]; then
        echo "Processing $sample"
        if ! mosdepth -n -t "$THREADS_PER_JOB" "$OUTDIR/$sample" "$bam"; then
            echo "$sample" >> "$OUTDIR/failed_samples.txt"
            return 0
        fi
    fi

    depth=$(awk '$1=="total" {print $4}' "$summary")
    echo -e "${sample}\t${depth}" > "$OUTDIR/$sample.depth.tmp"
}
export -f run_mosdepth
export OUTDIR THREADS_PER_JOB

ls "$BAM_DIR"/*.bam | xargs -P "$JOBS" -I{} bash -c 'run_mosdepth "$@"' _ {}

echo -e "sample\tmean_depth" > "$OUTDIR/depths.tsv"
cat "$OUTDIR"/*.depth.tmp | sort -k1,1 >> "$OUTDIR/depths.tsv"
rm -f "$OUTDIR"/*.depth.tmp

[[ -s "$OUTDIR/failed_samples.txt" ]] && echo "WARNING: some samples failed, see $OUTDIR/failed_samples.txt"

python3 << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

df = pd.read_csv("mosdepth_results/depths.tsv", sep="\t")
depths = df["mean_depth"]

print("\nSummary statistics:")
print(depths.describe())

XMAX = 46
ticks = np.arange(0, XMAX + 2, 2)  # 0, 2, 4, ..., 46


def set_xaxis():
    plt.xlim(0, XMAX)
    plt.xticks(ticks, [str(t) for t in ticks])


# Figure 1: histogram + KDE
plt.figure(figsize=(10, 6))
bins = np.arange(0, XMAX + 0.25, 0.25)
plt.hist(depths, bins=bins, density=True, alpha=0.5, label="Histogram")

kde = gaussian_kde(depths)
x = np.linspace(depths.min(), depths.max(), 5000)
plt.plot(x, kde(x), linewidth=2, color="red", label="Density")

plt.xlabel("Mean sequencing depth (X)")
plt.ylabel("Density")
plt.title("Sequencing Depth Distribution")
plt.legend()
set_xaxis()
plt.tight_layout()
plt.savefig("mosdepth_results/depth_combined.png", dpi=300)

# Figure 2: density curve
plt.figure(figsize=(8, 5))
depths.plot(kind="density")
plt.xlabel("Mean sequencing depth (X)")
plt.title("Sequencing Depth Density")
set_xaxis()
plt.tight_layout()
plt.savefig("mosdepth_results/depth_density.png", dpi=300)

# Figure 3: histogram + density
plt.figure(figsize=(8, 5))
depths.hist(bins=300, density=True, alpha=0.6)
depths.plot(kind="density")
plt.xlabel("Mean sequencing depth (X)")
plt.ylabel("Density")
plt.title("Sequencing Depth Distribution")
set_xaxis()
plt.tight_layout()
plt.savefig("mosdepth_results/depth_hist_density.png", dpi=300)
EOF

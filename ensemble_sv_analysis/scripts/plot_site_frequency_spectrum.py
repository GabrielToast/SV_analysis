#!/usr/bin/env python3
"""
Site frequency spectrum (SFS) of structural variants.

For each SV in the cohort-merged VCF, SURVIVOR's SUPP tag gives the number of
samples that carry it. The SFS is the histogram of that number: how many SVs
are found in exactly 1 sample, exactly 2 samples, ... up to all samples.

Use the UNFILTERED cohort_merged.vcf. The filtered file has the recurrent SVs
removed, so its spectrum is cut off.

Example
  python3 plot_sfs.py
  python3 plot_sfs.py --by-type --log-y
  python3 plot_sfs.py --vcf results/cohort/cohort_merged.vcf --out sfs.png

Outputs
  <out>.png        the plot
  <out>.tsv        the counts behind the plot (overall and per SV type)
"""

import argparse
import os
import re
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")  # no display needed on the cluster
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

SUPP_RE = re.compile(r"(?:^|;)SUPP=(\d+)")
SUPP_VEC_RE = re.compile(r"(?:^|;)SUPP_VEC=([0-9]+)")
SVTYPE_RE = re.compile(r"(?:^|;)SVTYPE=([A-Za-z]+)")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf",
                   default="/home/gato/ensemble_sv_analysis/results/cohort/cohort_merged.vcf",
                   help="cohort-merged VCF from SURVIVOR (unfiltered)")
    p.add_argument("--out",
                   default="/home/gato/ensemble_sv_analysis/results/cohort/sv_sfs.png",
                   help="output PNG path (a .tsv is written next to it)")
    p.add_argument("--n-samples", type=int, default=None,
                   help="total number of samples (default: inferred from SUPP_VEC)")
    p.add_argument("--by-type", action="store_true",
                   help="stack the bars by SV type (DEL, DUP, INS, INV, TRA)")
    p.add_argument("--log-y", action="store_true",
                   help="log-scale y axis")
    p.add_argument("--title", default="SV site frequency spectrum")
    return p.parse_args()


def read_spectrum(vcf_path):
    overall = Counter()
    by_type = defaultdict(Counter)
    n_vec = 0
    n_records = n_no_supp = 0
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            n_records += 1
            info = fields[7]
            m = SUPP_RE.search(info)
            if m is None:
                n_no_supp += 1
                continue
            k = int(m.group(1))
            overall[k] += 1
            t = SVTYPE_RE.search(info)
            by_type[t.group(1) if t else "UNKNOWN"][k] += 1
            v = SUPP_VEC_RE.search(info)
            if v:
                n_vec = max(n_vec, len(v.group(1)))
    return overall, by_type, n_vec, n_records, n_no_supp


def main():
    a = parse_args()
    if not os.path.isfile(a.vcf):
        sys.exit(f"ERROR: VCF not found: {a.vcf}")

    overall, by_type, n_vec, n_records, n_no_supp = read_spectrum(a.vcf)
    if not overall:
        sys.exit("ERROR: no records with a SUPP tag found; is this a SURVIVOR-merged VCF?")

    n = a.n_samples or max(n_vec, max(overall))
    ks = list(range(1, n + 1))

    dropped = sum(c for k, c in overall.items() if k > n or k < 1)
    if dropped:
        print(f"WARNING: {dropped} SVs have SUPP outside 1..{n} and are not plotted "
              f"(check --n-samples)", file=sys.stderr)

    # --- table ---
    tsv = os.path.splitext(a.out)[0] + ".tsv"
    types = sorted(by_type)
    with open(tsv, "w") as fh:
        fh.write("n_samples\tall\t" + "\t".join(types) + "\n")
        for k in ks:
            fh.write(f"{k}\t{overall.get(k, 0)}\t" +
                     "\t".join(str(by_type[t].get(k, 0)) for t in types) + "\n")

    # --- plot ---
    many = n > 40                      # switch to a compact style for big cohorts
    width_in = 10 if many else max(6, 0.6 * n + 3)
    fig, ax = plt.subplots(figsize=(width_in, 4.5))
    bar_kw = dict(width=1.0, linewidth=0) if many else dict(edgecolor="black", linewidth=0.4)

    if a.by_type:
        bottom = [0] * len(ks)
        for t in types:
            vals = [by_type[t].get(k, 0) for k in ks]
            ax.bar(ks, vals, bottom=bottom, label=t, **bar_kw)
            bottom = [b + v for b, v in zip(bottom, vals)]
        ax.legend(title="SV type", frameon=False)
        totals = bottom
    else:
        totals = [overall.get(k, 0) for k in ks]
        ax.bar(ks, totals, color="#4c72b0", **bar_kw)

    if not many:                       # per-bar count labels only when readable
        top = max(totals) if totals else 1
        for k, tot in zip(ks, totals):
            if tot:
                ax.text(k, tot + (top * 0.01 if not a.log_y else 0), str(tot),
                        ha="center", va="bottom", fontsize=8)

    if a.log_y:
        ax.set_yscale("log")
    if many:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=12))
        ax.set_xlim(0, n + 1)
    else:
        ax.set_xticks(ks)
    ax.set_xlabel("Number of samples carrying the SV")
    ax.set_ylabel("Number of SVs")
    ax.set_title(f"{a.title} (n = {n} samples, {sum(overall.values())} SVs)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(a.out, dpi=200)

    print(f"SVs plotted: {sum(totals)}  (records without SUPP skipped: {n_no_supp})")
    print(f"Samples (x-axis max): {n}")
    print(f"Plot:  {a.out}")
    print(f"Table: {tsv}")


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
Summarise cohort_filtered.vcf: SV length distribution and which sample each
SV came from.

Sample of origin
  SURVIVOR's merged VCF has one column per input VCF and a SUPP_VEC INFO tag
  (e.g. 0010000) saying which inputs carry the SV. Position i in SUPP_VEC is
  input i. Sample names are taken from cohort_vcf_list.txt (written by
  cohort_merge_filter.py; the filename is stripped of its suffix) when it
  matches the number of columns, otherwise from the VCF header columns.
  If an SV is carried by more than one sample (e.g. --max-samples 2), the
  per-sample tables count it once for every sample that carries it.

Example
  python3 cohort_stats.py

  python3 cohort_stats.py \
      --vcf /home/gato/ensemble_sv_analysis/results/cohort/cohort_filtered.vcf \
      --out-dir /home/gato/ensemble_sv_analysis/results/cohort/stats

Outputs (in --out-dir)
  sv_table.tsv               one row per SV and carrying sample
  length_summary.tsv         length stats (n, min, median, mean, max) overall and per SVTYPE
  length_bins.tsv            counts per length bin and SVTYPE
  sample_by_type.tsv         SV counts per sample and SVTYPE
  sample_length_summary.tsv  length stats per sample
  length_distribution.png    histogram of SV lengths (log x), stacked by SVTYPE
  svs_per_sample.png         SVs per sample, stacked by SVTYPE
  (plots are skipped if matplotlib is not installed)
"""

import argparse
import os
import re
import statistics
import sys
from collections import Counter, defaultdict

SVTYPE_RE = re.compile(r"(?:^|;)SVTYPE=([^;]+)")
SVLEN_RE = re.compile(r"(?:^|;)SVLEN=(-?\d+)")
END_RE = re.compile(r"(?:^|;)END=(\d+)")
SUPP_VEC_RE = re.compile(r"(?:^|;)SUPP_VEC=([01]+)")

BIN_EDGES = [50, 100, 250, 500, 1_000, 2_500, 5_000, 10_000,
             50_000, 100_000, 1_000_000]


def bin_label(length):
    def fmt(x):
        if x >= 1_000_000:
            return f"{x // 1_000_000}M"
        if x >= 1_000:
            return f"{x / 1_000:g}k"
        return str(x)

    if length < BIN_EDGES[0]:
        return f"<{fmt(BIN_EDGES[0])}"
    for lo, hi in zip(BIN_EDGES, BIN_EDGES[1:]):
        if lo <= length < hi:
            return f"{fmt(lo)}-{fmt(hi)}"
    return f">={fmt(BIN_EDGES[-1])}"


def bin_order():
    labels = [bin_label(BIN_EDGES[0] - 1)]
    labels += [bin_label(lo) for lo in BIN_EDGES[:-1]]
    labels.append(bin_label(BIN_EDGES[-1]))
    return labels


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf",
                   default="/home/gato/ensemble_sv_analysis/results/cohort/cohort_filtered.vcf",
                   help="filtered cohort VCF")
    p.add_argument("--list-file", default=None,
                   help="VCF list used for the merge (default: cohort_vcf_list.txt "
                        "next to --vcf)")
    p.add_argument("--out-dir", default=None,
                   help="output directory (default: <vcf dir>/stats)")
    return p.parse_args()


def sample_names_from_list(list_file, n):
    """Names from the merge input list, or None if unusable."""
    if not os.path.isfile(list_file):
        return None
    with open(list_file) as fh:
        paths = [l.strip() for l in fh if l.strip()]
    if len(paths) != n:
        return None
    names = []
    for p in paths:
        base = os.path.basename(p)
        base = re.sub(r"\.vcf(\.gz)?$", "", base)
        base = re.sub(r"_merged(_nohet)?$", "", base)
        names.append(base)
    return names


def sv_length(info, pos):
    m = SVLEN_RE.search(info)
    if m and int(m.group(1)) != 0:
        return abs(int(m.group(1)))
    m = END_RE.search(info)
    if m:
        return abs(int(m.group(1)) - pos)
    return None


def summarise(values):
    if not values:
        return ["0", "NA", "NA", "NA", "NA"]
    return [str(len(values)), str(min(values)),
            f"{statistics.median(values):g}",
            f"{statistics.mean(values):.1f}", str(max(values))]


def main():
    a = parse_args()
    vcf = os.path.expanduser(a.vcf)
    if not os.path.isfile(vcf):
        sys.exit(f"ERROR: VCF not found: {vcf}")
    out_dir = os.path.expanduser(a.out_dir or os.path.join(os.path.dirname(vcf), "stats"))
    os.makedirs(out_dir, exist_ok=True)
    list_file = os.path.expanduser(
        a.list_file or os.path.join(os.path.dirname(vcf), "cohort_vcf_list.txt"))

    header_samples = []
    records = []        # (chrom, pos, end, svtype, length, [samples])
    skipped_nolen = skipped_nosupp = 0
    names = None

    with open(vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_samples = line.rstrip("\n").split("\t")[9:]
                names = (sample_names_from_list(list_file, len(header_samples))
                         or header_samples)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, info = f[0], int(f[1]), f[7]
            m = SVTYPE_RE.search(info)
            svtype = m.group(1) if m else "NA"
            length = sv_length(info, pos)
            if length is None:
                skipped_nolen += 1
                continue
            m = END_RE.search(info)
            end = int(m.group(1)) if m else pos + length

            m = SUPP_VEC_RE.search(info)
            if m is None or names is None or len(m.group(1)) != len(names):
                skipped_nosupp += 1
                samples = ["unknown"]
            else:
                samples = [names[i] for i, c in enumerate(m.group(1)) if c == "1"]
                if not samples:
                    samples = ["unknown"]
            records.append((chrom, pos, end, svtype, length, samples))

    if not records:
        sys.exit("ERROR: no SV records with a usable length found")

    # ---- per-SV table
    with open(os.path.join(out_dir, "sv_table.tsv"), "w") as out:
        out.write("chrom\tpos\tend\tsvtype\tsvlen\tn_samples\tsample\n")
        for chrom, pos, end, svtype, length, samples in records:
            for s in samples:
                out.write(f"{chrom}\t{pos}\t{end}\t{svtype}\t{length}\t{len(samples)}\t{s}\n")

    types = sorted({r[3] for r in records})
    all_samples = sorted({s for r in records for s in r[5]})

    # ---- length summary overall and by type
    by_type = defaultdict(list)
    for r in records:
        by_type[r[3]].append(r[4])
    with open(os.path.join(out_dir, "length_summary.tsv"), "w") as out:
        out.write("group\tn\tmin\tmedian\tmean\tmax\n")
        out.write("\t".join(["ALL"] + summarise([r[4] for r in records])) + "\n")
        for t in types:
            out.write("\t".join([t] + summarise(by_type[t])) + "\n")

    # ---- length bins x type
    order = bin_order()
    bins = defaultdict(Counter)
    for r in records:
        bins[bin_label(r[4])][r[3]] += 1
    with open(os.path.join(out_dir, "length_bins.tsv"), "w") as out:
        out.write("bin\t" + "\t".join(types) + "\ttotal\n")
        for b in order:
            counts = [bins[b][t] for t in types]
            out.write(f"{b}\t" + "\t".join(map(str, counts)) + f"\t{sum(counts)}\n")

    # ---- per-sample tables
    sample_type = defaultdict(Counter)
    sample_len = defaultdict(list)
    for r in records:
        for s in r[5]:
            sample_type[s][r[3]] += 1
            sample_len[s].append(r[4])
    with open(os.path.join(out_dir, "sample_by_type.tsv"), "w") as out:
        out.write("sample\t" + "\t".join(types) + "\ttotal\n")
        for s in sorted(all_samples, key=lambda x: -sum(sample_type[x].values())):
            counts = [sample_type[s][t] for t in types]
            out.write(f"{s}\t" + "\t".join(map(str, counts)) + f"\t{sum(counts)}\n")
    with open(os.path.join(out_dir, "sample_length_summary.tsv"), "w") as out:
        out.write("sample\tn\tmin\tmedian\tmean\tmax\n")
        for s in sorted(all_samples):
            out.write("\t".join([s] + summarise(sample_len[s])) + "\n")

    # ---- console summary
    print(f"SVs read: {len(records)}  (samples seen: {len(all_samples)})")
    if skipped_nolen:
        print(f"  skipped (no SVLEN/END): {skipped_nolen}")
    if skipped_nosupp:
        print(f"  assigned 'unknown' sample (no usable SUPP_VEC): {skipped_nosupp}")
    print("\nBy type:")
    for t in types:
        print(f"  {t:<6} n={len(by_type[t]):<8} median length={statistics.median(by_type[t]):g}")
    top = sorted(all_samples, key=lambda x: -len(sample_len[x]))[:10]
    print("\nTop samples by SV count:")
    for s in top:
        print(f"  {s:<30} {len(sample_len[s])}")

    # ---- plots
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import math
    except ImportError:
        print("\nmatplotlib not installed: skipping plots")
        print(f"\nTables written to {out_dir}")
        return

    # length histogram, log x, stacked by type
    lo = max(1, min(r[4] for r in records))
    hi = max(r[4] for r in records)
    nb = 40
    edges = [10 ** (math.log10(lo) + i * (math.log10(hi + 1) - math.log10(lo)) / nb)
             for i in range(nb + 1)]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist([by_type[t] for t in types], bins=edges, stacked=True, label=types)
    ax.set_xscale("log")
    ax.set_xlabel("SV length (bp)")
    ax.set_ylabel("number of SVs")
    ax.set_title("SV length distribution")
    ax.legend(title="SVTYPE")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "length_distribution.png"), dpi=150)
    plt.close(fig)

    # SVs per sample, stacked by type (top 40 samples if the cohort is large)
    ordered = sorted(all_samples, key=lambda x: -len(sample_len[x]))
    shown = ordered[:40]
    fig, ax = plt.subplots(figsize=(max(8, 0.35 * len(shown) + 3), 5))
    bottoms = [0] * len(shown)
    for t in types:
        vals = [sample_type[s][t] for s in shown]
        ax.bar(range(len(shown)), vals, bottom=bottoms, label=t)
        bottoms = [b + v for b, v in zip(bottoms, vals)]
    ax.set_xticks(range(len(shown)))
    ax.set_xticklabels(shown, rotation=90, fontsize=7)
    ax.set_ylabel("number of SVs")
    title = "SVs per sample"
    if len(ordered) > len(shown):
        title += f" (top {len(shown)} of {len(ordered)})"
    ax.set_title(title)
    ax.legend(title="SVTYPE")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "svs_per_sample.png"), dpi=150)
    plt.close(fig)

    print(f"\nTables and plots written to {out_dir}")


if __name__ == "__main__":
    main()
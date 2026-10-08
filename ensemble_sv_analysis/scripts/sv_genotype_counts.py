#!/usr/bin/env python3
"""
Count genotypes (hom-ref / het / hom-alt / missing) for the SVs in the
Snakemake pipeline's VCFs.

It reads two kinds of files:
  callers  the per-caller plain VCFs (results/{manta,delly,smoove}/<sample>/<sample>_<caller>.vcf)
  merged   the per-sample SURVIVOR-merged VCFs (results/merged/<sample>_merged.vcf)

Counting from the per-caller VCFs is the reliable way to get each caller's own
genotype. The merged VCFs have one sample column per input caller, but SURVIVOR
rewrites the FORMAT fields, so this script also checks whether the merged
genotypes look usable (e.g. not all missing) and shows how the callers'
genotypes agree for the same SV.

Genotype classes (haploid calls "0" and "1" count as hom_ref and hom_alt):
  hom_ref  0/0      het  0/1 or 1/0      hom_alt  1/1
  missing  ./. or . other  anything else (e.g. 1/2)

Example
  python3 sv_genotype_counts.py
  python3 sv_genotype_counts.py --pass-only --min-size 50
  python3 sv_genotype_counts.py --vcf some_file.vcf other_file.vcf

Outputs (in --out-dir)
  genotype_counts.tsv          counts per file, sample column and SV type
  genotype_patterns_merged.tsv how the callers' genotypes combine per merged SV
"""

import argparse
import glob
import os
import re
import sys
from collections import Counter, defaultdict

CLASSES = ["hom_ref", "het", "hom_alt", "missing", "other"]
SVTYPE_RE = re.compile(r"(?:^|;)SVTYPE=([A-Za-z]+)")
SVLEN_RE = re.compile(r"(?:^|;)SVLEN=(-?\d+)")
END_RE = re.compile(r"(?:^|;)END=(\d+)")
SPLIT_RE = re.compile(r"[/|]")


def classify(gt):
    """Map a GT string like 0/1, 1|1, 1 or ./. to one of CLASSES."""
    alleles = SPLIT_RE.split(gt.strip())
    if not alleles or all(a in (".", "") for a in alleles):
        return "missing"
    if any(a in (".", "") for a in alleles):
        return "missing"
    try:
        nums = {int(a) for a in alleles}
    except ValueError:
        return "other"
    if nums == {0}:
        return "hom_ref"
    if nums == {1}:
        return "hom_alt"
    if nums == {0, 1}:
        return "het"
    return "other"


def sv_size(info):
    m = SVLEN_RE.search(info)
    if m:
        return abs(int(m.group(1)))
    return None


def parse_vcf(path, pass_only, min_size):
    """Yield (svtype, [class per sample column]) for each record, plus the
    sample column names. Returns (sample_names, records)."""
    samples = []
    records = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            if pass_only and f[6] not in ("PASS", "."):
                continue
            info = f[7]
            if min_size:
                s = sv_size(info)
                if s is not None and s < min_size:
                    continue
            t = SVTYPE_RE.search(info)
            svtype = t.group(1) if t else "UNKNOWN"
            classes = []
            if len(f) > 9 and len(f) >= 9 + len(samples):
                fmt = f[8].split(":")
                gt_idx = fmt.index("GT") if "GT" in fmt else None
                for col in f[9:9 + len(samples)]:
                    if gt_idx is None:
                        classes.append("missing")
                    else:
                        parts = col.split(":")
                        classes.append(classify(parts[gt_idx]) if gt_idx < len(parts) else "missing")
            else:
                classes = ["missing"] * max(len(samples), 1)
                if not samples:
                    samples = ["(no sample column)"]
            records.append((svtype, classes))
    return samples, records


def sample_from_name(path):
    base = os.path.basename(path)
    for suffix in ("_merged.vcf", "_manta.vcf", "_delly.vcf", "_smoove.vcf"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base.rsplit(".", 1)[0]


def discover(results_dir):
    files = []  # (kind, label, path)
    for caller in ("manta", "delly", "smoove"):
        pat = os.path.join(results_dir, caller, "*", f"*_{caller}.vcf")
        for p in sorted(glob.glob(pat)):
            files.append(("callers", caller, p))
    for p in sorted(glob.glob(os.path.join(results_dir, "merged", "*_merged.vcf"))):
        files.append(("merged", "merged", p))
    return files


def pct(n, d):
    return f"{100.0 * n / d:.1f}%" if d else "NA"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="/home/gato/ensemble_sv_analysis/results")
    ap.add_argument("--out-dir", default=None,
                    help="default: <results-dir>/genotype_qc")
    ap.add_argument("--vcf", nargs="+", default=None,
                    help="count these VCFs instead of auto-discovering the pipeline outputs")
    ap.add_argument("--pass-only", action="store_true",
                    help="only count records with FILTER = PASS (or '.')")
    ap.add_argument("--min-size", type=int, default=0,
                    help="skip SVs smaller than this (only when SVLEN is present)")
    a = ap.parse_args()

    if a.vcf:
        files = [("custom", "custom", os.path.expanduser(p)) for p in a.vcf]
    else:
        files = discover(os.path.expanduser(a.results_dir))
    if not files:
        sys.exit("ERROR: no VCFs found. Check --results-dir or pass --vcf.")

    out_dir = os.path.expanduser(a.out_dir or os.path.join(a.results_dir, "genotype_qc"))
    os.makedirs(out_dir, exist_ok=True)

    rows = []  # (kind, label, sample, column, svtype, counts)
    patterns = defaultdict(Counter)  # sample -> Counter of tuple(classes)
    merged_columns = {}

    for kind, label, path in files:
        sample = sample_from_name(path)
        try:
            cols, records = parse_vcf(path, a.pass_only, a.min_size)
        except OSError as e:
            print(f"WARNING: could not read {path}: {e}", file=sys.stderr)
            continue
        if kind == "merged":
            merged_columns[sample] = cols
        by = defaultdict(lambda: [Counter() for _ in cols])
        for svtype, classes in records:
            for i, c in enumerate(classes[: len(cols)]):
                by[svtype][i][c] += 1
            if kind == "merged":
                patterns[sample][tuple(classes)] += 1
        for svtype, per_col in by.items():
            for i, cnt in enumerate(per_col):
                rows.append((kind, label, sample, cols[i], svtype, cnt))

    # ---- write the long table ----
    tsv = os.path.join(out_dir, "genotype_counts.tsv")
    with open(tsv, "w") as fh:
        fh.write("source\tlabel\tsample\tsample_column\tsvtype\t" + "\t".join(CLASSES) + "\ttotal\n")
        for kind, label, sample, col, svtype, cnt in sorted(rows):
            total = sum(cnt.values())
            fh.write(f"{kind}\t{label}\t{sample}\t{col}\t{svtype}\t" +
                     "\t".join(str(cnt.get(c, 0)) for c in CLASSES) + f"\t{total}\n")

    # ---- console summary: per source/column (all SV types pooled) ----
    pooled = defaultdict(Counter)
    for kind, label, sample, col, svtype, cnt in rows:
        key = (kind, label, sample, col)
        pooled[key].update(cnt)

    print(f"{'source':<8} {'sample':<14} {'column':<24} {'hom_ref':>8} {'het':>8} {'hom_alt':>8} "
          f"{'missing':>8} {'other':>6} {'het/(het+alt)':>14}")
    for key in sorted(pooled):
        kind, label, sample, col = key
        c = pooled[key]
        called = c["het"] + c["hom_alt"]
        name = label if kind == "callers" else col
        print(f"{kind:<8} {sample:<14} {name[:24]:<24} {c['hom_ref']:>8} {c['het']:>8} {c['hom_alt']:>8} "
              f"{c['missing']:>8} {c['other']:>6} {pct(c['het'], called):>14}")

    # ---- did the SURVIVOR merge keep usable genotypes? ----
    if merged_columns:
        print("\nMerged VCF genotype check")
        for sample in sorted(merged_columns):
            tot = Counter()
            for key, c in pooled.items():
                if key[0] == "merged" and key[2] == sample:
                    tot.update(c)
            n = sum(tot.values())
            if n and tot["missing"] == n:
                print(f"  {sample}: ALL genotypes are missing -> the merge did not keep GT; "
                      f"use the per-caller counts instead.")
            elif n and (tot["hom_ref"] + tot["missing"]) == n:
                print(f"  {sample}: no het/hom_alt calls at all -> the merge probably did not "
                      f"keep GT; use the per-caller counts instead.")
            else:
                print(f"  {sample}: genotypes present ({n} calls across {len(merged_columns[sample])} column(s)); "
                      f"compare with the per-caller counts before trusting them.")

        # pattern table
        ptsv = os.path.join(out_dir, "genotype_patterns_merged.tsv")
        with open(ptsv, "w") as fh:
            fh.write("sample\tcolumns\tpattern\tn_svs\n")
            for sample in sorted(patterns):
                cols = ",".join(merged_columns[sample])
                for pat, n in patterns[sample].most_common():
                    fh.write(f"{sample}\t{cols}\t{'|'.join(pat)}\t{n}\n")
        print(f"Patterns: {ptsv}")

    print(f"\nTable: {tsv}")


if __name__ == "__main__":
    main()

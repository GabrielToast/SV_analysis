#!/usr/bin/env python3
"""
Merge all per-sample SV VCFs into one cohort VCF, keep only SVs that appear in
a single sample, and remove unresolved breakpoints (BNDs/breakends).

How it works
  1. Writes a list of the per-sample VCFs.
  2. Runs `SURVIVOR merge` on that list with min_callers=1 (union). Each input
     VCF counts as one "sample", so every distinct SV appears once and gets a
     SUPP INFO tag = number of samples carrying it (SUPP_VEC shows which).
  3. Keeps only records with SUPP <= --max-samples (default 1, i.e. private
     to one sample) and drops BND records (SVTYPE=BND, or an ALT allele
     written in breakend notation such as N[chr2:100[ ).

Example
  python3 cohort_merge_filter.py

  python3 cohort_merge_filter.py \
      --vcf-glob "/home/gato/ensemble_sv_analysis/results/merged_nohet/*_merged_nohet.vcf" \
      --out-dir  /home/gato/ensemble_sv_analysis/results/cohort \
      --max-samples 1

  # keep BNDs as well
  python3 cohort_merge_filter.py --keep-bnd

Outputs (in --out-dir)
  cohort_vcf_list.txt     list of input VCFs
  cohort_merged.vcf       all SVs, before filtering
  cohort_filtered.vcf     SVs present in <= --max-samples samples, no BNDs
  cohort_survivor.log     SURVIVOR's stdout/stderr
"""

import argparse
import glob
import os
import re
import subprocess
import sys

SUPP_RE = re.compile(r"(?:^|;)SUPP=(\d+)")
SVTYPE_RE = re.compile(r"(?:^|;)SVTYPE=([^;]+)")
BND_ALT_RE = re.compile(r"[\[\]]")  # breakend notation in ALT, e.g. A]chr2:123]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf-glob",
                   default="/home/gato/ensemble_sv_analysis/results/merged_nohet/*_merged_nohet.vcf",
                   help="glob matching the per-sample VCFs to combine")
    p.add_argument("--out-dir",
                   default="/home/gato/ensemble_sv_analysis/results/cohort",
                   help="output directory (created if missing)")
    p.add_argument("--survivor",
                   default="/home/gato/miniconda3/envs/survivor/bin/SURVIVOR",
                   help="path to the SURVIVOR binary")
    p.add_argument("--max-samples", type=int, default=1,
                   help="keep SVs present in at most this many samples (default 1)")
    p.add_argument("--keep-bnd", action="store_true",
                   help="do not remove BND (breakend) records")
    p.add_argument("--max-dist", type=int, default=1000,
                   help="SURVIVOR: max distance between breakpoints (default 1000)")
    p.add_argument("--type-match", type=int, default=1, choices=[0, 1],
                   help="SURVIVOR: require same SV type (default 1)")
    p.add_argument("--strand-match", type=int, default=0, choices=[0, 1],
                   help="SURVIVOR: require same strand (default 0)")
    p.add_argument("--dist-est", type=int, default=0, choices=[0, 1],
                   help="SURVIVOR: estimate distance based on SV size (default 0)")
    p.add_argument("--min-size", type=int, default=50,
                   help="SURVIVOR: minimum SV size in bp (default 50)")
    return p.parse_args()


def is_bnd(fields):
    """True if the record is a breakend (SVTYPE=BND or breakend-style ALT)."""
    m = SVTYPE_RE.search(fields[7])
    if m is not None and m.group(1) == "BND":
        return True
    alt = fields[4]
    return alt == "<BND>" or BND_ALT_RE.search(alt) is not None


def filter_records(merged_vcf, filtered_vcf, max_samples, drop_bnd):
    kept = dropped_recurrent = dropped_bnd = no_supp = 0
    with open(merged_vcf) as fin, open(filtered_vcf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            if drop_bnd and is_bnd(fields):
                dropped_bnd += 1
                continue
            m = SUPP_RE.search(fields[7])
            if m is None:
                # No SUPP tag, so recurrence can't be judged: keep, but report
                no_supp += 1
                fout.write(line)
            elif int(m.group(1)) <= max_samples:
                kept += 1
                fout.write(line)
            else:
                dropped_recurrent += 1
    return kept, dropped_recurrent, dropped_bnd, no_supp


def main():
    a = parse_args()

    vcfs = sorted(os.path.realpath(f) for f in glob.glob(os.path.expanduser(a.vcf_glob)))
    if len(vcfs) < 2:
        sys.exit(f"ERROR: need at least 2 VCFs, found {len(vcfs)} for glob: {a.vcf_glob}")
    if not os.path.isfile(a.survivor):
        sys.exit(f"ERROR: SURVIVOR not found at {a.survivor} (use --survivor)")

    out_dir = os.path.expanduser(a.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    list_file = os.path.join(out_dir, "cohort_vcf_list.txt")
    merged = os.path.join(out_dir, "cohort_merged.vcf")
    filtered = os.path.join(out_dir, "cohort_filtered.vcf")
    log_file = os.path.join(out_dir, "cohort_survivor.log")

    with open(list_file, "w") as fh:
        fh.write("\n".join(vcfs) + "\n")
    print(f"Merging {len(vcfs)} VCFs:")
    for v in vcfs:
        print(f"  {v}")

    # SURVIVOR merge <list> <max_dist> <min_callers> <type> <strand> <dist_est> <min_size> <out>
    # min_callers=1 -> union, so SUPP can then be used to count samples per SV
    cmd = [a.survivor, "merge", list_file, str(a.max_dist), "1",
           str(a.type_match), str(a.strand_match), str(a.dist_est),
           str(a.min_size), merged]
    with open(log_file, "w") as lf:
        rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
    if rc != 0 or not os.path.isfile(merged):
        sys.exit(f"ERROR: SURVIVOR failed (exit {rc}); see {log_file}")

    kept, dropped_rec, dropped_bnd, no_supp = filter_records(
        merged, filtered, a.max_samples, drop_bnd=not a.keep_bnd)

    print(f"\nFilter: keep SVs present in <= {a.max_samples} sample(s)"
          f"{'' if a.keep_bnd else ', BNDs removed'}")
    print(f"  kept:                  {kept}")
    print(f"  dropped (recurrent):   {dropped_rec}")
    print(f"  dropped (BND):         {dropped_bnd}")
    print(f"  kept, no SUPP tag:     {no_supp}")
    print(f"\nUnfiltered: {merged}")
    print(f"Filtered:   {filtered}")


if __name__ == "__main__":
    main()
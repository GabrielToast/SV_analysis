#!/usr/bin/env python3
"""
Rank de novo (single-sample) SVs so the most promising ones get inspected first.

Input table (denovo_table.tsv), tab-separated, no header:
    chr  start  end  svtype  svlen  sample  gt  [qual]
The 8th column (QUAL, phred-scaled, from the cohort VCF) is optional; without
it the QUAL term is skipped.

Scoring (higher = better candidate):
    +2  all three callers support the call in the carrier's per-sample VCF (SUPP=3)
    +1  two callers support it (SUPP=2)
    +2  length 100 bp - 50 kb
    +1  length 50-99 bp or 50 kb - 200 kb
    +2  overlaps a CDS
    +1  overlaps a gene (but no CDS)
    +2  QUAL in the top third of all calls in the table
    +1  QUAL in the middle third
    -2  a breakpoint lies within 100 bp of a repeat
    -2  within 10 kb of a contig end
    -3  fewer than --low-sibs siblings (default 4): weak control set
    +2  more than --high-sibs siblings (default 5, i.e. 6 or more)
    -3  sample has no usable BAM in bam_map.tsv (--nobam-penalty); these calls
        cannot be checked against siblings and are left out of <out>_to_plot.tsv

Siblings = the other samples of the same colony listed in bam_map.tsv
(made by make_bam_map.py). Samples with several BAMs, or without an index, are
not in the map, so they are never counted as siblings.

Outputs:
    <out>_all.tsv      every call with its annotations and score, sorted by
                       score (ties: higher QUAL first). Columns 1-15 as before
                       (score = 15), then n_siblings, colony, has_bam.
                       Column 6 is the cleaned sample ID (same as bam_map.tsv).
    <out>_to_plot.tsv  top N calls + a random set of M others, taken only from
                       calls with a usable BAM (7 columns, no header)

Example:
    python rank_svs_v3.py \
      --table denovo_table_remapped.tsv \
      --merged-dir results/merged_nohet \
      --merged-suffix _merged_nohet.vcf \
      --gff genomic.gff \
      --fai GCF_003254395.2_Amel_HAv3.1_genomic.fna.fai \
      --repeats repeats.bed \
      --bam-map bam_map.tsv \
      --top 200 --random 50 --out ranked
"""
import argparse
import bisect
import os
import random
import re
from collections import defaultdict

p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--table", required=True, help="denovo_table.tsv (7-8 columns)")
p.add_argument("--merged-dir", required=True,
               help="folder with per-sample merged VCFs (for caller support)")
p.add_argument("--merged-suffix", default="_merged.vcf",
               help="filename suffix of per-sample merged VCFs "
                    "(e.g. _merged_nohet.vcf)")
p.add_argument("--gff", required=True, help="genomic.gff (NCBI GFF3)")
p.add_argument("--fai", required=True, help="reference .fai index")
p.add_argument("--repeats", default=None, help="optional BED of repeat regions")
p.add_argument("--bam-map", default="bam_map.tsv",
               help="bam_map.tsv from make_bam_map.py (sample_id, bam, colony)")
p.add_argument("--top", type=int, default=200)
p.add_argument("--random", type=int, default=50)
p.add_argument("--window", type=int, default=1000,
               help="bp window for matching a call in the per-sample VCF")
p.add_argument("--out", default="ranked")
p.add_argument("--low-sibs", type=int, default=4,
               help="fewer siblings than this is penalised (default 4)")
p.add_argument("--high-sibs", type=int, default=5,
               help="more siblings than this earns a bonus (default 5, i.e. 6+)")
p.add_argument("--sib-penalty", type=int, default=3)
p.add_argument("--sib-bonus", type=int, default=2)
p.add_argument("--nobam-penalty", type=int, default=3,
               help="points subtracted when the sample has no usable BAM")
a = p.parse_args()


def canon(x):
    """Same cleaning rule as make_bam_map.py / check_ids.py."""
    x = os.path.basename(x.strip())
    x = re.sub(r"\.bam$", "", x)
    x = x.replace("-", "_")
    x = re.sub(r"_S\d+(_.*)?$", "", x)
    x = re.sub(r"_(Amel|rg|rmd|sorted|marked)(_.*)?$", "", x)
    return x


class Intervals:
    """Overlap queries against a set of intervals per contig."""

    def __init__(self):
        self.d = defaultdict(list)
        self.starts = {}
        self.pmax = {}

    def add(self, chrom, s, e):
        self.d[chrom].append((s, e))

    def build(self):
        for c, lst in self.d.items():
            lst.sort()
            self.starts[c] = [s for s, _ in lst]
            m, pm = -1, []
            for _, e in lst:
                m = max(m, e)
                pm.append(m)
            self.pmax[c] = pm

    def hits(self, chrom, s, e):
        if chrom not in self.starts:
            return False
        i = bisect.bisect_right(self.starts[chrom], e)
        return i > 0 and self.pmax[chrom][i - 1] >= s


# ---- annotations ----------------------------------------------------------
cds, genes, reps = Intervals(), Intervals(), Intervals()
with open(a.gff) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        if f[2] == "CDS":
            cds.add(f[0], int(f[3]), int(f[4]))
        elif f[2] == "gene":
            genes.add(f[0], int(f[3]), int(f[4]))
cds.build()
genes.build()

if a.repeats:
    with open(a.repeats) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) >= 3:
                reps.add(f[0], int(f[1]), int(f[2]))
    reps.build()

contig_len = {}
with open(a.fai) as fh:
    for line in fh:
        f = line.split("\t")
        contig_len[f[0]] = int(f[1])

# ---- bam_map: usable BAMs, colonies, siblings -------------------------------
in_map = set()
colony_of = {}
members = defaultdict(set)
with open(a.bam_map) as fh:
    for line in fh:
        f = line.rstrip("\r\n").split("\t")
        if len(f) < 3 or not f[0] or f[0] == "sample_id":
            continue
        sid, col = f[0].strip(), f[2].strip()
        in_map.add(sid)
        if col and col != "NA":
            colony_of[sid] = col
            members[col].add(sid)
print(f"bam_map: {len(in_map)} usable samples in {len(members)} colonies")


def n_siblings(sid):
    col = colony_of.get(sid)
    if col is None:
        return None
    return len(members[col] - {sid})


# ---- caller support from the carrier's per-sample merged VCF ---------------
cache = {}


def load_merged(sample):
    """Return {chrom: [(pos, SUPP), ...]} for the sample's merged VCF."""
    if sample in cache:
        return cache[sample]
    short = re.sub(r"_S\d+$", "", sample)  # P27152_283_S171 -> P27152_283
    path = None
    for name in (sample, short):
        cand = os.path.join(a.merged_dir, f"{name}{a.merged_suffix}")
        if os.path.exists(cand):
            path = cand
            break
    recs = defaultdict(list)
    if path:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 8)
                supp = None
                for kv in f[7].split(";"):
                    if kv.startswith("SUPP="):
                        supp = int(kv[5:])
                if supp is not None:
                    recs[f[0]].append((int(f[1]), supp))
    cache[sample] = recs
    return recs


def callers(sample, chrom, pos):
    recs = load_merged(sample)
    if not recs:
        return None
    hits = [s for p_, s in recs.get(chrom, []) if abs(p_ - pos) <= a.window]
    return max(hits) if hits else None


# ---- read the table and work out QUAL thresholds ---------------------------
table = []
with open(a.table) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 7 and f[5]:
            table.append(f)


def parse_qual(f):
    if len(f) < 8 or f[7] in (".", ""):
        return None
    try:
        return float(f[7])
    except ValueError:
        return None


quals = sorted(q for q in (parse_qual(f) for f in table) if q is not None)
if quals:
    q_lo = quals[len(quals) // 3]
    q_hi = quals[(2 * len(quals)) // 3]
else:
    q_lo = q_hi = None

# ---- score every call ------------------------------------------------------
rows = []
for f in table:
    chrom, start, end, svtype, svlen, sample, gt = f[:7]
    sid = canon(sample)
    start = int(start)
    end = int(end) if end not in (".", "") else start
    length = abs(end - start)
    if svlen not in (".", ""):
        try:
            length = abs(int(svlen))
        except ValueError:
            pass

    qual = parse_qual(f)
    n_call = callers(sample, chrom, start)
    in_cds = cds.hits(chrom, start, end)
    in_gene = genes.hits(chrom, start, end)
    near_rep = bool(a.repeats) and (
        reps.hits(chrom, start - 100, start + 100)
        or reps.hits(chrom, end - 100, end + 100))
    clen = contig_len.get(chrom, 0)
    near_end = start < 10000 or (clen and end > clen - 10000)
    has_bam = sid in in_map
    nsib = n_siblings(sid) if has_bam else None

    score = 0
    if n_call == 3:
        score += 2
    elif n_call == 2:
        score += 1
    if 100 <= length <= 50000:
        score += 2
    elif 50 <= length < 100 or 50000 < length <= 200000:
        score += 1
    if in_cds:
        score += 2
    elif in_gene:
        score += 1
    if qual is not None and q_hi is not None:
        if qual >= q_hi:
            score += 2
        elif qual >= q_lo:
            score += 1
    if near_rep:
        score -= 2
    if near_end:
        score -= 2
    if not has_bam:
        score -= a.nobam_penalty
    elif nsib is not None:
        if nsib < a.low_sibs:
            score -= a.sib_penalty
        elif nsib > a.high_sibs:
            score += a.sib_bonus

    rows.append([chrom, start, end, svtype, svlen, sid, gt, length,
                 n_call if n_call is not None else "NA",
                 int(bool(in_cds)), int(bool(in_gene)),
                 int(bool(near_rep)), int(bool(near_end)),
                 qual if qual is not None else "NA",
                 score,
                 nsib if nsib is not None else "NA",
                 colony_of.get(sid, "NA"),
                 int(has_bam)])

# score descending; ties broken by QUAL descending (stable otherwise)
rows.sort(key=lambda r: (-r[14], -(r[13] if r[13] != "NA" else -1)))

hdr = ["chr", "start", "end", "type", "svlen", "sample", "gt", "length",
       "callers", "in_cds", "in_gene", "repeat_bp", "near_contig_end",
       "qual", "score", "n_siblings", "colony", "has_bam"]
with open(f"{a.out}_all.tsv", "w") as o:
    o.write("\t".join(hdr) + "\n")
    for r in rows:
        o.write("\t".join(map(str, r)) + "\n")

eligible = [r for r in rows if r[17] == 1]
top = eligible[:a.top]
rest = eligible[a.top:]
random.seed(1)
rnd = random.sample(rest, min(a.random, len(rest)))
with open(f"{a.out}_to_plot.tsv", "w") as o:
    for r in top + rnd:
        o.write("\t".join(map(str, r[:7])) + "\n")

n_na = sum(1 for r in rows if r[8] == "NA")
n_nobam = sum(1 for r in rows if r[17] == 0)
n_low = sum(1 for r in rows if r[15] != "NA" and r[15] < a.low_sibs)
n_high = sum(1 for r in rows if r[15] != "NA" and r[15] > a.high_sibs)
n_nocol = sum(1 for r in rows if r[17] == 1 and r[15] == "NA")
print(f"{len(rows)} calls scored; wrote {a.out}_all.tsv and "
      f"{a.out}_to_plot.tsv ({len(top)} top + {len(rnd)} random)")
if quals:
    print(f"QUAL thresholds: middle third >= {q_lo}, top third >= {q_hi}")
else:
    print("No QUAL values found (8th column missing or all '.'); "
          "QUAL term skipped.")
if rows and n_na == len(rows):
    print("WARNING: caller support is NA for every call - check --merged-dir, "
          "--merged-suffix and the sample names.")
elif n_na:
    print(f"Note: {n_na} calls had no matching record in the per-sample VCFs "
          f"(callers = NA).")
print(f"Siblings: {n_low} calls penalised (<{a.low_sibs}), {n_high} rewarded "
      f"(>{a.high_sibs}), {n_nocol} with a BAM but no colony.")
print(f"No usable BAM: {n_nobam} calls (penalised by {a.nobam_penalty}, "
      f"kept in {a.out}_all.tsv, excluded from {a.out}_to_plot.tsv).")
if rows and n_nobam == len(rows):
    print("WARNING: no sample matched bam_map.tsv - check the sample IDs.")
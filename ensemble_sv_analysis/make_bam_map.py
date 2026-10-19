#!/usr/bin/env python3
import glob, os, re, sys
from collections import defaultdict

BAMDIR   = "/home/gato/haploids_bam/bam"
COLONIES = "colonies.tsv"
TABLE    = sys.argv[1] if len(sys.argv) > 1 else "denovo_table_remapped.tsv"

def canon(x):
    x = os.path.basename(x.strip())
    x = re.sub(r"\.bam$", "", x)
    x = x.replace("-", "_")
    x = re.sub(r"_S\d+(_.*)?$", "", x)
    x = re.sub(r"_(Amel|rg|rmd|sorted|marked)(_.*)?$", "", x)
    return x

bams = defaultdict(list)
for f in glob.glob(f"{BAMDIR}/*.bam"):
    bams[canon(f)].append(f)

dups   = {k: v for k, v in bams.items() if len(v) > 1}
unique = {k: v[0] for k, v in bams.items() if len(v) == 1}

noidx = [k for k, f in unique.items()
         if not (os.path.exists(f + ".bai") or os.path.exists(f[:-4] + ".bai"))]
for k in noidx:
    del unique[k]

col = {}
with open(COLONIES) as fh:
    for line in fh:
        f = line.rstrip("\r\n").split("\t")
        if len(f) >= 2 and f[0].strip() and f[1].strip():
            col[canon(f[0])] = f[1].strip()

with open("bam_map.tsv", "w") as o:
    o.write("sample_id\tbam\tcolony\n")
    for k in sorted(unique):
        o.write(f"{k}\t{unique[k]}\t{col.get(k, 'NA')}\n")

with open("excluded_samples.txt", "w") as o:
    for k in sorted(dups):
        o.write(f"{k}\tmultiple_bams\t" + ",".join(os.path.basename(x) for x in dups[k]) + "\n")
    for k in sorted(noidx):
        o.write(f"{k}\tno_index\t{os.path.basename(bams[k][0])}\n")

n_nocol = sum(1 for k in unique if k not in col)
print(f"usable BAMs: {len(unique)}  (excluded: {len(dups)} with multiple BAMs, {len(noidx)} without index)")
print(f"usable BAMs with no colony: {n_nocol}")
for k in sorted(dups): print("   excluded duplicate:", k)

tab = set()
with open(TABLE) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 6 and f[5]:
            tab.add(canon(f[5]))
t_dup  = sorted(t for t in tab if t in dups)
t_nobam = sorted(t for t in tab if t not in unique and t not in dups)
t_nocol = sorted(t for t in tab if t in unique and t not in col)
print(f"\ntable ({TABLE}): {len(tab)} samples")
print(f"  in the table but excluded as duplicate BAM: {len(t_dup)} {t_dup[:6]}")
print(f"  in the table with no BAM at all:            {len(t_nobam)} {t_nobam[:6]}")
print(f"  in the table with a BAM but no colony:      {len(t_nocol)} {t_nocol[:6]}")

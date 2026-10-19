#!/usr/bin/env python3
import re, os, sys, glob
from collections import defaultdict

BAMDIR = "/home/gato/haploids_bam/bam"
TABLE = sys.argv[1] if len(sys.argv) > 1 else "denovo_table_remapped.tsv"
COLONIES = "colonies.tsv"

def canon(x):
    x = os.path.basename(x.strip())
    x = re.sub(r"\.bam$", "", x)
    x = x.replace("-", "_")
    x = re.sub(r"_S\d+(_.*)?$", "", x)
    x = re.sub(r"_(Amel|rg|rmd|sorted|marked)(_.*)?$", "", x)
    return x

bam = defaultdict(list)
for f in glob.glob(f"{BAMDIR}/*.bam"):
    bam[canon(f)].append(os.path.basename(f))
dups = {k: v for k, v in bam.items() if len(v) > 1}

col = {}
with open(COLONIES) as fh:
    for line in fh:
        f = line.rstrip("\r\n").split("\t")
        if len(f) >= 2 and f[0].strip():
            col[canon(f[0])] = f[1].strip()

tab = defaultdict(set)
with open(TABLE) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 6 and f[5]:
            tab[canon(f[5])].add(f[5])

print(f"BAMs: {len(bam)} unique IDs; {len(dups)} IDs shared by >1 file")
for k, v in dups.items(): print("   duplicate:", k, v)
print(f"colonies.tsv: {len(col)} IDs")
print(f"table ({TABLE}): {len(tab)} unique samples")

no_bam = [k for k in tab if k not in bam]
no_col = [k for k in tab if k not in col]
print(f"\ntable samples with no BAM:    {len(no_bam)}")
for k in no_bam[:8]: print("   ", k, "<-", sorted(tab[k]))
print(f"table samples with no colony: {len(no_col)}")
for k in no_col[:8]: print("   ", k, "<-", sorted(tab[k]))
bam_nocol = [k for k in bam if k not in col]
print(f"BAMs with no colony:          {len(bam_nocol)}")
for k in bam_nocol[:8]: print("   ", k, bam[k])

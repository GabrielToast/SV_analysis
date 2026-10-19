#!/usr/bin/env python
import argparse, csv, glob, os
from collections import Counter
import pysam

p = argparse.ArgumentParser()
p.add_argument("--depth-table", default="samplot_top20/depth_table.tsv")
p.add_argument("--bam-map", default="bam_map.tsv")
p.add_argument("--ref", default=None, help="reference FASTA (optional, enables junction check)")
p.add_argument("--outdir", default="samplot_top20/breakpoints")
p.add_argument("--tol", type=int, default=10, help="bp tolerance around breakpoint for clips")
p.add_argument("--min-clip", type=int, default=8, help="min clipped bases to count")
p.add_argument("--k", type=int, default=12, help="junction-adjacent bases used for consensus")
p.add_argument("--min-mapq", type=int, default=10)
p.add_argument("--min-support", type=int, default=3, help="reads needed to call a carrier")
p.add_argument("--pair-win", type=int, default=500)
p.add_argument("--sa-tol", type=int, default=200)
a = p.parse_args()
os.makedirs(a.outdir, exist_ok=True)

ref = pysam.FastaFile(a.ref) if a.ref else None
_bams = {}

_map = {}
with open(a.bam_map) as _fh:
    for _line in _fh:
        _f = _line.rstrip("\r\n").split("\t")
        if len(_f) >= 2 and _f[0] and _f[0] != "sample_id":
            _map[_f[0]] = _f[1]

def get_bam(sample):
    if sample not in _bams:
        path = _map.get(sample)
        _bams[sample] = pysam.AlignmentFile(path, "rb") if path else None
    return _bams[sample]

def sa_near(tag, chrom, pos):
    for e in tag.split(";"):
        if not e: continue
        f = e.split(",")
        if f[0] == chrom and abs(int(f[1]) - pos) <= a.sa_tol: return True
    return False

def consensus(seqs):
    if not seqs: return "", 0.0
    L = len(seqs[0]); out = []; agree = 0
    for i in range(L):
        base, n = Counter(s[i] for s in seqs).most_common(1)[0]
        out.append(base); agree += n
    return "".join(out), agree / (L * len(seqs))

def match_other(cons, side, chrom, other):
    """Does the clipped consensus match the reference near the other breakpoint?"""
    if ref is None or not cons: return "NA"
    lo = max(0, other - 30 - 1)
    seg = ref.fetch(chrom, lo, other + 30 + len(cons)).upper()
    best = None
    for i in range(len(seg) - len(cons) + 1):
        mm = sum(1 for x, y in zip(cons, seg[i:i + len(cons)]) if x != y)
        if mm <= 1 and (best is None or mm < best[0]):
            edge = (lo + i + 1) if side == "R" else (lo + i + len(cons))
            best = (mm, edge - other)
    return "no" if best is None else f"yes(mm={best[0]},off={best[1]:+d})"

def analyse_bp(bam, chrom, bp, other):
    R, L = [], []
    nR = nL = span = sa = tot = 0
    for r in bam.fetch(chrom, max(0, bp - a.tol - 1), bp + a.tol + 1):
        if r.is_unmapped or r.is_secondary or r.is_duplicate or r.mapping_quality < a.min_mapq:
            continue
        ct = r.cigartuples
        if not ct: continue
        tot += 1
        left = ct[0][1] if ct[0][0] in (4, 5) else 0
        right = ct[-1][1] if ct[-1][0] in (4, 5) else 0
        rs, re_ = r.reference_start, r.reference_end
        qs = r.query_sequence
        is_sa = r.has_tag("SA") and sa_near(r.get_tag("SA"), chrom, other)
        if right >= a.min_clip and abs(re_ - bp) <= a.tol:
            nR += 1; sa += is_sa
            if ct[-1][0] == 4 and qs:
                s = qs[len(qs) - right: len(qs) - right + a.k]
                if len(s) == a.k: R.append(s)
        elif left >= a.min_clip and abs(rs - bp) <= a.tol:
            nL += 1; sa += is_sa
            if ct[0][0] == 4 and qs and left >= a.k:
                L.append(qs[left - a.k:left])
        elif left < a.min_clip and right < a.min_clip and rs < bp - a.tol and re_ > bp + a.tol:
            span += 1
    side = "R" if nR >= nL else "L"
    cons, agree = consensus(R if side == "R" else L)
    return dict(nR=nR, nL=nL, span=span, sa=sa, tot=tot, side=side, ndom=max(nR, nL),
                cons=cons, agree=agree, match=match_other(cons, side, chrom, other))

def count_disc(bam, chrom, bp1, bp2):
    pw = min(a.pair_win, max(50, (bp2 - bp1) // 3))
    n = 0
    for r in bam.fetch(chrom, max(0, bp1 - pw), bp1 + pw):
        if r.is_unmapped or r.is_secondary or r.is_duplicate or r.mapping_quality < a.min_mapq: continue
        if r.is_proper_pair or r.mate_is_unmapped or r.next_reference_id != r.reference_id: continue
        if r.reference_start <= r.next_reference_start and bp2 - pw <= r.next_reference_start <= bp2 + pw:
            n += 1
    return n

# ---- read candidates and samples from the depth table ----
cands = {}
with open(a.depth_table) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        cands.setdefault(row["candidate"], []).append(row)

ev_cols = ["candidate", "sample", "role", "chrom", "bp1", "bp2", "type",
           "bp1_clipR", "bp1_clipL", "bp1_span", "bp1_SA", "bp1_side", "bp1_consensus", "bp1_agree", "bp1_match_other",
           "bp2_clipR", "bp2_clipL", "bp2_span", "bp2_SA", "bp2_side", "bp2_consensus", "bp2_agree", "bp2_match_other",
           "disc_pairs", "support", "span_frac", "read_call", "depth_call"]
evf = open(f"{a.outdir}/breakpoint_evidence.tsv", "w")
evf.write("\t".join(ev_cols) + "\n")
smf = open(f"{a.outdir}/candidate_summary.tsv", "w")
smf.write("candidate\tcolony\tfocal_read_call\tn_controls\tcontrol_carriers_reads\tcontrol_carriers_depth\tcontrols_weak_or_mixed\tverdict\n")
rep = open(f"{a.outdir}/junction_comparison.txt", "w")

for cand, rows in cands.items():
    chrom, bp1, bp2, typ = rows[0]["chrom"], int(rows[0]["start"]), int(rows[0]["end"]), rows[0]["type"]
    colony = rows[0].get("colony", "")
    print(f"{cand}: {len(rows)} samples")
    rep.write(f"\n{'=' * 100}\n{cand}  colony={colony}  {typ} {chrom}:{bp1}-{bp2}\n")
    if ref is not None:
        rep.write(f"ref @bp1  {ref.fetch(chrom, bp1 - 20, bp1).upper()} | {ref.fetch(chrom, bp1, bp1 + 20).upper()}\n")
        rep.write(f"ref @bp2  {ref.fetch(chrom, bp2 - 20, bp2).upper()} | {ref.fetch(chrom, bp2, bp2 + 20).upper()}\n")
        if typ == "DEL":
            rep.write(f"expected  {ref.fetch(chrom, bp1 - a.k, bp1).upper()} | {ref.fetch(chrom, bp2, bp2 + a.k).upper()}"
                      "   (left-of-bp1 | right-of-bp2, may shift by a few bp)\n")
    rep.write(f"{'sample':<14}{'role':<8}{'call':<11}{'dep':<8}| bp1 R/L/span/SA  seq(side)          match"
              f"{'':<18}| bp2 R/L/span/SA  seq(side)          match{'':<18}| disc\n")

    focal_call = "NA"; ctrl_read = []; ctrl_depth = []; ctrl_weak = []; nctrl = 0
    for row in rows:
        sample, role, depth_call = row["sample"], row["role"], row["carrier"]
        bam = get_bam(sample)
        if bam is None:
            print(f"  no BAM for {sample}"); continue
        e1 = analyse_bp(bam, chrom, bp1, bp2)
        e2 = analyse_bp(bam, chrom, bp2, bp1)
        disc = count_disc(bam, chrom, bp1, bp2)
        support = disc + max(e1["ndom"], e2["ndom"])
        fr = []
        for e in (e1, e2):
            d = e["span"] + e["ndom"]
            if d > 0: fr.append(e["span"] / d)
        span_frac = sum(fr) / len(fr) if fr else None
        if e1["tot"] < 3 and e2["tot"] < 3: call = "lowcov"
        elif support >= a.min_support: call = "carrier" if (span_frac is None or span_frac < 0.5) else "mixed"
        elif support > 0: call = "weak"
        else: call = "none"
        sf = "NA" if span_frac is None else f"{span_frac:.2f}"
        evf.write("\t".join(map(str, [cand, sample, role, chrom, bp1, bp2, typ,
            e1["nR"], e1["nL"], e1["span"], e1["sa"], e1["side"], e1["cons"], f'{e1["agree"]:.2f}', e1["match"],
            e2["nR"], e2["nL"], e2["span"], e2["sa"], e2["side"], e2["cons"], f'{e2["agree"]:.2f}', e2["match"],
            disc, support, sf, call, depth_call])) + "\n")
        def cell(e):
            seq = f'{e["cons"] or "-"}({e["side"]})'
            return f'{e["nR"]}/{e["nL"]}/{e["span"]}/{e["sa"]}'.ljust(17) + seq.ljust(20) + e["match"].ljust(23)
        rep.write(f"{sample:<14}{role:<8}{call:<11}{depth_call:<8}| {cell(e1)}| {cell(e2)}| {disc}\n")
        if role == "focal": focal_call = call
        else:
            nctrl += 1
            if call == "carrier": ctrl_read.append(sample)
            elif call in ("weak", "mixed"): ctrl_weak.append(sample)
            if depth_call == "yes": ctrl_depth.append(sample)

    if focal_call != "carrier": verdict = "FOCAL_NOT_SUPPORTED_BY_READS"
    elif ctrl_read or ctrl_depth: verdict = "SHARED_WITH_SIBLINGS"
    elif ctrl_weak: verdict = "POSSIBLY_SHARED_WEAK_CONTROL_EVIDENCE"
    else: verdict = "DENOVO_SUPPORTED_BY_READS"
    smf.write(f"{cand}\t{colony}\t{focal_call}\t{nctrl}\t{','.join(ctrl_read) or '-'}\t{','.join(ctrl_depth) or '-'}\t{','.join(ctrl_weak) or '-'}\t{verdict}\n")

evf.close(); smf.close(); rep.close()
print(f"\nWrote: {a.outdir}/candidate_summary.tsv, breakpoint_evidence.tsv, junction_comparison.txt")

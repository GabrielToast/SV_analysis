#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")   # headless cluster, no display
import pandas as pd
import matplotlib.pyplot as plt

# SV table (supp_vec is a 773-digit string: read it as text)
sv = pd.read_csv("svs.bed", sep="\t",
                 names=["chrom", "start", "end", "sv_id", "svtype", "svlen", "supp_vec"],
                 dtype={"supp_vec": str, "sv_id": str})
samples = open("samples_vcf_order.txt").read().split()
print(sv.supp_vec.str.len().unique(), len(samples))
print(sv.supp_vec.str.count("1").value_counts())

sv["sample"] = sv.supp_vec.apply(
    lambda v: ",".join(samples[i] for i, c in enumerate(v) if c == "1"))


def load(path, label):
    """Read a bedtools intersect -wo file: gene is column 10, overlap bp is the last column."""
    try:
        d = pd.read_csv(path, sep="\t", header=None, dtype=str)   # everything as text
    except pd.errors.EmptyDataError:                              # no overlaps for this feature
        return pd.DataFrame(columns=["sv_id", "gene", "ov_bp", "region"])
    d = d[[3, 10]].assign(ov_bp=pd.to_numeric(d.iloc[:, -1]))
    d.columns = ["sv_id", "gene", "ov_bp"]
    d = d.groupby("sv_id").agg(gene=("gene", lambda x: ",".join(sorted(set(x)))),
                               ov_bp=("ov_bp", "sum")).reset_index()
    d["region"] = label
    return d


hits = pd.concat([load("sv_cds.tsv", "CDS"),
                  load("sv_exon.tsv", "UTR/non-coding exon"),
                  load("sv_gene.tsv", "intronic")]
                 ).drop_duplicates("sv_id", keep="first")   # keeps highest priority

res = sv.merge(hits, on="sv_id", how="left")
res["region"] = res["region"].fillna("intergenic")
res.drop(columns="supp_vec").to_csv("sv_gene_annotation.tsv", sep="\t", index=False)
print(res.region.value_counts())

# Overview plot: SV type x region
order = ["CDS", "UTR/non-coding exon", "intronic", "intergenic"]
ct = res.groupby(["region", "svtype"]).size().unstack(fill_value=0).reindex(order).fillna(0)
ct.plot(kind="bar", stacked=True, figsize=(7, 4), edgecolor="black", linewidth=0.4)
plt.ylabel("Number of de novo SVs"); plt.xlabel(""); plt.xticks(rotation=20, ha="right")
plt.tight_layout(); plt.savefig("sv_region_summary.png", dpi=200)
plt.close()

# Per-sample burden of coding SVs
coding = res[res.region == "CDS"]
if len(coding):
    coding.groupby("sample").size().plot(kind="hist", bins=30, figsize=(6, 4), edgecolor="black")
    plt.xlabel("Coding SVs per sample"); plt.tight_layout()
    plt.savefig("coding_sv_per_sample.png", dpi=200)
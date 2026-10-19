import pandas as pd

XLSX  = "Drones_GWAS_updated_10032025.xlsx"
SHEET = "20240409_Drone_Metadata_TE_rm_w"# sheet name or index
OUT   = "colonies.tsv"

df = pd.read_excel(XLSX, sheet_name=SHEET, dtype=str)
df.columns = df.columns.str.strip()
df = df[["NGI_ID", "Colony_ID"]].copy()

for c in df.columns:
    df[c] = df[c].str.strip()
df = df.dropna(subset=["NGI_ID", "Colony_ID"])
df = df[(df.NGI_ID != "") & (df.Colony_ID != "")]

# reduce NGI_ID to the first two underscore fields (e.g. P27152_182_S82 -> P27152_182),
# which is how plot_top20.sh names each BAM
df["NGI_ID"] = df["NGI_ID"].str.split("_").str[:2].str.join("_")

# sanity checks
dup = df.groupby("NGI_ID")["Colony_ID"].nunique()
bad = dup[dup > 1]
if len(bad):
    print("WARNING: IDs assigned to more than one colony:\n", bad)
df = df.drop_duplicates()

df.to_csv(OUT, sep="\t", header=False, index=False)
print(f"Wrote {len(df)} rows to {OUT}")
print(df.groupby("Colony_ID").size().describe())

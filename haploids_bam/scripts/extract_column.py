import pandas as pd

df = pd.read_excel("Drones_GWAS_updated_10032025.xlsx")

# Extract columns 2 and 3 (Excel columns B and C)
df.iloc[:, [1, 2]].to_csv("columns_2_3.tsv", sep="\t", index=False)

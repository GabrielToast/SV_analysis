import pandas as pd
import os

# Excel metadata file
metadata_file = "metadata.xlsx"

# Directory containing BAM files
bam_dir = "/path/to/bam_directory"

# Read Excel file
df = pd.read_excel(metadata_file)

# Change 'SampleID' to the actual column name in your Excel file
metadata_samples = set(df["SampleID"].astype(str))

# Get sample names from BAM files
bam_samples = set()

for f in os.listdir(bam_dir):
    if f.endswith(".bam"):
        sample = f.replace(".bam", "")
        bam_samples.add(sample)

# Compare
missing_bams = metadata_samples - bam_samples
extra_bams = bam_samples - metadata_samples

print(f"Samples in metadata: {len(metadata_samples)}")
print(f"BAM files found: {len(bam_samples)}")
print(f"Missing BAMs: {len(missing_bams)}")
print(f"Extra BAMs: {len(extra_bams)}")

print("\nMissing BAM files:")
for sample in sorted(missing_bams):
    print(sample)

print("\nExtra BAM files:")
for sample in sorted(extra_bams):
    print(sample)

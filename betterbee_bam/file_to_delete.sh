# Extract metadata IDs (first 10 characters of column 2)
awk 'NR>1 {print substr($2,1,10)}' metadata_simple_text.txt | sort -u > metadata_ids.txt

# Find BAMs whose sample ID (first 10 chars) is not in metadata
for bam in *.bam; do
    sample=$(echo "$bam" | sed 's/_S[0-9].*\.bam$//')
    sample10=${sample:0:10}

    if ! grep -qxF "$sample10" metadata_ids.txt; then
        echo "$bam"
    fi
done > bam_files_without_metadata.txt


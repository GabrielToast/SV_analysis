for bam in all_bams/*.bam; do
    n=$(samtools view -H "$bam" | grep "^@RG" | \
        sed 's/\t/\n/g' | grep "^SM:" | sort -u | wc -l)

    if [ "$n" -gt 1 ]; then
        echo "$bam : $n samples"
    fi
done

#!/bin/bash
bam_list=($(ls *_rg_rmd.bam))
for bam in "${bam_list[@]}"
do
    echo "Processing BAM file: $bam" >> log.txt
        # Run GATK HaplotypeCaller for each BAM file
        gatk HaplotypeCaller \
            -I "$bam" \
            -O "${bam%.bam}_chr$chr.gvcf.gz" \
            -R /home/hieu/try_hieu/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
            -ERC GVCF \
            -ploidy 1 &
    wait
    echo "Finished processing BAM file: $bam" >> log.txt
done

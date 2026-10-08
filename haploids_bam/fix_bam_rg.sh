#!/usr/bin/env bash
set -euo pipefail

SRC=/home/gato/haploids_bam/bam

fix_one() {
    bam=$1

    # Safety guard: never touch symlinks
    if [ -L "$bam" ]; then
        echo "SYMLINK $(basename "$bam")"
        return 0
    fi

    # Only touch BAMs whose SM tag contains a colon
    if ! samtools view -H "$bam" | grep '^@RG' | grep -qP '\tSM:[^\t]*:'; then
        echo "SKIP    $(basename "$bam")"
        return 0
    fi

    tmp="${bam%.bam}.tmpfix.bam"
    hdr="${bam%.bam}.tmpfix.header.sam"

    samtools view -H "$bam" | awk 'BEGIN{FS=OFS="\t"}
        /^@RG/ { for(i=2;i<=NF;i++) if($i ~ /^SM:/){ v=substr($i,4); gsub(/:/,"_",v); $i="SM:" v } }
        {print}' > "$hdr"

    samtools reheader "$hdr" "$bam" > "$tmp"

    # Verify before replacing: file intact and same read count
    samtools quickcheck "$tmp"
    [ "$(samtools idxstats "$bam" | awk '{s+=$3+$4} END{print s}')" = \
      "$(samtools view -c "$tmp")" ] || { echo "COUNT MISMATCH $bam"; rm -f "$tmp" "$hdr"; return 1; }

    # Replace original, keep timestamp, rebuild index
    touch -r "$bam" "$tmp"
    mv -f "$tmp" "$bam"
    samtools index -@ 2 "$bam"
    touch -r "$bam" "$bam.bai"
    rm -f "$hdr"
    echo "FIXED   $(basename "$bam")"
}
export -f fix_one

# -type f = regular files only (excludes symlinks); also skip leftover temp files
find "$SRC" -maxdepth 1 -type f -name '*.bam' ! -name '*.tmpfix.bam' -print0 \
  | xargs -0 -P 4 -I{} bash -c 'fix_one "$1"' _ {} 2>&1 | tee fix_sm_inplace.log
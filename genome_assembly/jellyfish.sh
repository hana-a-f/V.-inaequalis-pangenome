#!/bin/bash
export PATH=/programs/jellyfish-2.3.0/bin:$PATH

# Loop through all .fastq files
for fq in VI_*.fastq
do
    # Extract sample name (e.g., VI_19_004 from VI_19_004.fastq)
    sample=$(basename "$fq" .fastq)
    echo "Processing $sample..."

    # Run jellyfish count
    jellyfish count -m 21 -s 100M -t 10 -C -o ${sample}_kmers_21.jf "$fq"

    # Run jellyfish histo
    jellyfish histo ${sample}_kmers_21.jf > ${sample}_kmers_21.histo
done

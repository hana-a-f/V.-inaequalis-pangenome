#!/bin/bash

# Number of threads for hifiasm
THREADS=32

# Loop through all FASTQ files
for fq in VI_*.fastq; do
    # Extract base name (e.g., VI_1797_2)
    sample=$(basename "$fq" .fastq)

    echo "Starting HiFiASM for $sample..."

    # Make output directory
    mkdir -p "$sample"
    cd "$sample" || continue
    export PATH=/programs/hifiasm-0.19.9:$PATH
    # Run hifiasm with full path to FASTQ (since we're now in a subfolder)
    hifiasm -o "${sample}_hifi.asm" -t $THREADS --telo-m GGGTTA "../$fq" > "${sample}.log" 2>&1

    # Check success
    if [ $? -ne 0 ]; then
        echo "$sample failed. See $sample/$sample.log"
    else
        echo "$sample completed."
    fi

    # Return to parent directory
    cd ..
done

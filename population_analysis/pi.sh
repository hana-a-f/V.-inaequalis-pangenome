#!/bin/bash
set -euo pipefail

VCF_SNP="pangenie.snps.miss10.maf01.vcf.gz"
VCF_INS="pangenie.ins.NONREDUNDANT.vcf.gz"
VCF_DEL="pangenie.del.NONREDUNDANT.vcf.gz"

declare -A VCF_MAP
VCF_MAP["SNP"]=$VCF_SNP
VCF_MAP["INS"]=$VCF_INS
VCF_MAP["DEL"]=$VCF_DEL

pops=("EU" "US" "CAM" "CAP" "Rvi6YES" "Rvi6NO" "NonMalus")

echo -e "Population\tVariant\tMean_Pi" > Pi_summary.txt

run_pi() {

    POP=$1
    TYPE=$2
    INPUT_VCF=${VCF_MAP[$TYPE]}

    PREFIX="${POP}.${TYPE}"

    echo "Calculating π for ${POP} (${TYPE})"

    vcftools \
        --gzvcf ${INPUT_VCF} \
        --keep ${POP}.samples.txt \
        --window-pi 100000 \
        --out ${PREFIX}

    MEAN=$(awk 'NR>1 {sum+=$5; n++} END {if(n>0) print sum/n; else print "NA"}' ${PREFIX}.windowed.pi)

    echo -e "${POP}\t${TYPE}\t${MEAN}" >> Pi_summary.txt
}

for POP in "${pops[@]}"; do
    for TYPE in SNP INS DEL; do
        run_pi $POP $TYPE
    done
done

echo "Done."

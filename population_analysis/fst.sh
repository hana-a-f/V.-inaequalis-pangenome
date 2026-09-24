#!/bin/bash
set -euo pipefail

############################################
# INPUT VCF FILES
############################################

VCF_SNP="pangenie.snps.miss10.maf01.vcf.gz"
VCF_INS="pangenie.ins.gt50bp.maf01.miss20.vcf.gz"
VCF_DEL="pangenie.del.gt50bp.maf01.miss20.vcf.gz"

############################################
# Map variant types
############################################

declare -A VCF_MAP
VCF_MAP["SNP"]=$VCF_SNP
VCF_MAP["INS"]=$VCF_INS
VCF_MAP["DEL"]=$VCF_DEL

############################################
# Comparisons
############################################

comparisons=(
"EU US"
"EU CAM"
"EU CAP"
"US CAM"
"US CAP"
"CAM CAP"
"Malus NonMalus"
"Rvi6YES Rvi6NO"
)

############################################
# Output summary file
############################################

echo -e "Pop1\tPop2\tVariant\tMean_FST" > FST_summary.txt

############################################
# Function
############################################

run_fst() {

    POP1=$1
    POP2=$2
    TYPE=$3
    INPUT_VCF=${VCF_MAP[$TYPE]}

    PREFIX="${POP1}.${POP2}.${TYPE}"

    echo "Running ${POP1} vs ${POP2} (${TYPE})"

    # 1. Subset samples
    bcftools view \
        -S ${POP1}.v.${POP2}.samples.txt \
        "$INPUT_VCF" \
        -Oz -o ${PREFIX}.raw.vcf.gz

    tabix -p vcf ${PREFIX}.raw.vcf.gz

    # 2. Recalculate MAF + missingness
    bcftools +fill-tags \
        ${PREFIX}.raw.vcf.gz \
        -Oz -o ${PREFIX}.filled.vcf.gz \
        -- -t MAF,F_MISSING

    tabix -p vcf ${PREFIX}.filled.vcf.gz

    # 3. Filter within subset
    bcftools view \
        -i 'INFO/MAF >= 0.01 && INFO/F_MISSING <= 0.2' \
        ${PREFIX}.filled.vcf.gz \
        -Oz -o ${PREFIX}.filtered.vcf.gz

    tabix -p vcf ${PREFIX}.filtered.vcf.gz

    # 4. Site-level FST
    vcftools \
        --gzvcf ${PREFIX}.filtered.vcf.gz \
        --weir-fst-pop ${POP1}.samples.txt \
        --weir-fst-pop ${POP2}.samples.txt \
        --out ${PREFIX}

    # 5. Compute mean FST
    MEAN=$(awk 'NR>1 && $3!="nan" {sum+=$3; n++} END {if(n>0) print sum/n; else print "NA"}' ${PREFIX}.weir.fst)

    echo -e "${POP1}\t${POP2}\t${TYPE}\t${MEAN}" >> FST_summary.txt

    echo "Mean FST = ${MEAN}"
    echo
}

############################################
# Run all comparisons
############################################

for comp in "${comparisons[@]}"; do

    POP1=$(echo $comp | awk '{print $1}')
    POP2=$(echo $comp | awk '{print $2}')

    for TYPE in SNP INS DEL; do
        run_fst $POP1 $POP2 $TYPE
    done

done

echo "All mean FST analyses complete."



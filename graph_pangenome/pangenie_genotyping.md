# PanGenie Genotyping

This folder documents the **PanGenie** pipeline: genotyping short-read sequencing data from  *V. inaequalis* populations against the pangenome graph built in `pangenome_graph/`, then filtering to a population-level SV and SNP set for selection / PCA / GWAS-style downstream work.

## Pipeline overview

1. Install PanGenie (Singularity build from EBLERJANA's `pangenie.def`)
2. Prepare the graph VCF: decompose bubbles (`vcfbub`), filter by missingness, decompose to biallelic (PanGenie's `annotate_vcf.py`)
3. Force all genotypes to phased diploid (PanGenie requires this, even for haploid organisms)
4. Build a PanGenie index against the reference (VI_19_031)
5. Interleave paired-end FASTQs for every population sample
6. Genotype every sample against the graph
7. Sort, restore contig headers (PanGenie strips them), and index each per-sample VCF
8. Merge into a single population VCF
9. Filter: mask het calls (Venturia is haploid → het = genotyping error), mask low-GQ
10. Fill INFO tags (MAF, F_MISSING) and filter by MAF ≥ 1% and missingness ≤ 20%
11. Split by variant type: SNPs / insertions / deletions
12.run PLINK PCA

## Requirements

- The graph VCF must be:
  - **Phased** , haploid organism here
  - **Sequence-resolved** (not just structural like a raw MC bubble VCF)
  - **Non-overlapping** — no two variants at the same position
- PanGenie treats genotypes as diploid regardless of the underlying organism, so haploid samples (like *V. inaequalis*) we made diploid, this shouldn't be a problem since highly homozygous plant species have been used with this tool

## Tools used

| Tool | Version / source |
| --- | --- |
| PanGenie | Built from https://raw.githubusercontent.com/eblerjana/pangenie/master/container/pangenie.def |
| vcfbub | Bundled with PGGB (`pggb_latest.sif`) |
| bcftools | (system) |
| SeqFu | 1.18.0 (`/programs/SeqFu-1.18.0/seqfu`) |
| PLINK | 1.9 beta 7 (`/programs/plink-1.9-x86_64-beta7`) |
| tabix / bgzip | (system) |
| GNU parallel | (system) |
| External: `annotate_vcf.py` | PanGenie repo |
| External: `convert-to-biallelic.py` | PanGenie repo |


## Files in this folder

| File | Purpose |
| --- | --- |
| `README.md` | This file |
| `pangenie_genotyping.sh` | Loop over interleaved FASTQs, genotype each with PanGenie, convert to biallelic |
| `headers_sort_index.sh` | Add contig headers back to each per-sample PanGenie output, sort, index |

---

## Step 1 — Install PanGenie

```bash
/programs/bin/labutils/fakeroot
wget https://raw.githubusercontent.com/eblerjana/pangenie/master/container/pangenie.def
singularity build --fakeroot pangenie.sif pangenie.def
```

## Step 2 — Prepare the graph VCF

### 2a. vcfbub (from PGGB) — expand bubbles

```bash
singularity run --bind $PWD --pwd $PWD pggb_latest.sif vcfbub \
  -l 0 -r 100000 \
  -i venturia-pangenome.raw.vcf.gz \
  > venturia-pangenome.raw.vcf.bub.r100k.gz
```

### 2b. Filter by missingness (>20%)

```bash
bcftools filter -e 'F_MISSING > 0.20' \
  venturia-pangenome.raw.vcf.bub.r100k.gz \
  -Oz -o venturia-pangenome.raw.vcf.bub.r100k.f20.gz

bcftools index venturia-pangenome.raw.vcf.bub.r100k.f20.gz
gunzip venturia-pangenome.raw.vcf.bub.r100k.f20.gz
```

### 2c. PanGenie decomposition (multi- to biallelic)

```bash
python3 annotate_vcf.py \
  -vcf venturia-pangenome.raw.vcf.bub.r100k.f20 \
  -gfa venturia-pangenome.gfa.gz \
  -o mc.pangenie
```

Result:
```
Wrote 969672 multi-allelic records.
Wrote 1623227 bi-allelic records.
```

### 2d. Sort + normalize

```bash
# Sort graph VCF
bcftools sort --write-index -Oz -o mc.pangenie.sort.vcf.gz mc.pangenie.vcf

# Sort and normalize biallelic VCF
bcftools norm -f VI_19_031_pangenie.fasta -Ou mc.pangenie_biallelic.vcf \
  | bcftools sort --write-index -Oz -o mc.pangenie.biallelic.sort.vcf.gz
```

## Step 3 — Force phased diploid genotypes

PanGenie assumes diploid. *V. inaequalis* is haploid, so we duplicate each haploid call to a homozygous diploid; unphased hets become missing (real hets shouldn't exist in a haploid).

**The same awk block is applied to both `mc.pangenie.sort.vcf.gz` and `mc.pangenie.biallelic.sort.vcf.gz`:**

```bash
bcftools view <INPUT.vcf.gz> | \
awk '
BEGIN { FS=OFS="\t" }

/^#/ { print; next }

{
    split($9, fmt, ":")
    fmt_n = length(fmt)

    for (i=10; i<=NF; i++) {
        split($i, fields, ":")
        gt = fields[1]

        # Rule 1 — missing
        if (gt == "." || gt == "./." || gt == ".|." || gt == "") {
            fields[1] = ".|."
        }
        # Rule 2 — haploid numeric → duplicated
        else if (gt ~ /^[0-9]+$/) {
            fields[1] = gt "|" gt
        }
        # Rule 3 — unphased diploid: only keep if homozygous
        else if (gt ~ /^[0-9]+\/[0-9]+$/) {
            split(gt, g2, "/")
            if (g2[1] == g2[2]) fields[1] = g2[1] "|" g2[1]
            else                fields[1] = ".|."
        }
        # Rule 4 — already phased
        else if (gt ~ /^[0-9]+\|[0-9]+$/) {
            fields[1] = gt
        }

        sample = fields[1]
        for (j=2; j<=fmt_n; j++) sample = sample ":" fields[j]
        $i = sample
    }

    print
}' | bgzip > <OUTPUT.vcf.gz>
```

Apply to both files:

```bash
# Full graph
… > mc.pangenie.sort.diploid.vcf.gz

# Biallelic graph
… > mc.pangenie.biallelic.sort.diplod.vcf.gz
```

## Step 4 — Build PanGenie index

```bash
gunzip mc.pangenie.sort.diploid.vcf.gz

singularity run --bind $PWD --pwd $PWD pangenie.sif PanGenie-index \
  -v mc.pangenie.sort.diploid.vcf \
  -r VI_19_031_pangenie.fasta \
  -o mc.pangenie.index

# Prepare the biallelic reference for later biallelic conversion
gunzip mc.pangenie.biallelic.sort.diplod.vcf.gz
bgzip  mc.pangenie.biallelic.sort.diplod.vcf
```

## Step 5 — Interleave paired-end FASTQs

PanGenie takes interleaved reads. SeqFu handles it:

```bash
#!/bin/bash
set -euo pipefail

ILV="/programs/SeqFu-1.18.0/seqfu ilv"

interleave_one() {
    r1="$1"
    prefix="${r1%%_1P.fastq.gz}"
    r2="${prefix}_2P.fastq.gz"
    out="${prefix}.interleaved.fastq.gz"

    [[ ! -f "$r2" ]] && { echo "Missing R2 for $prefix — skipping"; return; }
    [[ -f "$out" ]] && { echo "Already exists: $out — skipping"; return; }

    echo "Interleaving $prefix → $out"
    $ILV -1 "$r1" -2 "$r2" | bgzip > "$out"
}

export -f interleave_one
export ILV

parallel -j 12 interleave_one ::: *_1P.fastq.gz
```

## Step 6 — Genotype all samples (`pangenie_genotyping.sh`)

Loops over every interleaved FASTQ, runs PanGenie against the index, converts each per-sample VCF to biallelic:

```bash
#!/bin/bash

cp /workdir/hf332/pangenie/short_reads/*.interleaved.fastq.gz .
gunzip -k *.interleaved.fastq.gz || true

REF_DIPLOID_VCF="mc.pangenie.biallelic.sort.diplod.vcf.gz"

for f in *.interleaved.fastq; do
    sample="${f%.interleaved.fastq}"
    echo "=== Running PanGenie for sample: $sample ==="

    singularity run --bind $PWD --pwd $PWD pangenie.sif \
        PanGenie \
            -f mc.pangenie.index \
            -i "$f" \
            -o "$sample" \
            -s "$sample" \
            -j 4 -t 12

    INVCF="${sample}_genotyping.vcf"
    [[ ! -f "$INVCF" ]] && { echo "ERROR: PanGenie output not found: $INVCF"; continue; }

    echo "Converting to biallelic diploid VCF…"
    cat "$INVCF" | \
        python3 convert-to-biallelic.py "$REF_DIPLOID_VCF" | \
        bgzip > "${sample}.pangenie.biallelic.raw.vcf.gz"

    echo "Done: ${sample}.pangenie.biallelic.raw.vcf.gz"
done

echo "=== ALL SAMPLES COMPLETED ==="
```

## Step 7 — Restore contig headers, sort, index (`headers_sort_index.sh`)

PanGenie's output drops `##contig` headers, which breaks bcftools sort and downstream tools. This script grabs them from the reference VCF and prepends them:

```bash
#!/bin/bash

# Grab contig headers from the reference biallelic VCF (only need to do this once)
bcftools view -h mc.pangenie.biallelic.sort.diplod.vcf.gz | grep "^##contig" > contigs.hdr

for f in *.pangenie.biallelic.raw.vcf.gz; do
    echo "Processing: $f"
    sample="${f%.pangenie.biallelic.raw.vcf.gz}"

    # 1) Add contig headers to VCF
    zcat "$f" \
      | awk '
          BEGIN {added=0}
          /^#CHROM/ && !added {
              while ((getline line < "contigs.hdr") > 0) print line
              close("contigs.hdr")
              added=1
          }
          {print}
        ' \
      | bgzip > "${sample}.withcontigs.vcf.gz"

    # 2) Sort
    bcftools sort "${sample}.withcontigs.vcf.gz" -Oz -o "${sample}.sorted.vcf.gz"

    # 3) Index
    tabix -p vcf "${sample}.sorted.vcf.gz"
done
```

## Step 8 — Merge per-sample VCFs into a population VCF

```bash
bcftools merge SRR*.sorted.vcf.gz -Oz -o population.pangenie.sorted.vcf.gz
tabix -p vcf population.pangenie.sorted.vcf.gz

bcftools stats population.pangenie.sorted.vcf.gz > population.stats.txt
```

## Step 9 — Mask spurious het calls, mask low-GQ, drop problem sample

*V. inaequalis* is haploid, so any heterozygous PanGenie call is a genotyping error — set to missing:

```bash
bcftools +setGT \
  -Oz -o population.pangenie.sorted.noHET.vcf.gz \
  population.pangenie.sorted.vcf.gz \
  -- -t q -n . -i 'GT="het"'
tabix -p vcf population.pangenie.sorted.noHET.vcf.gz
```

Mask low-support calls (GQ < 15):

```bash
bcftools +setGT \
  -Oz -o population.pangenie.nohetGT.noLowGQ.vcf.gz \
  population.pangenie.sorted.noHET.vcf.gz \
  -- -t q -n . -i 'FMT/GQ < 15'
tabix -p vcf population.pangenie.nohetGT.noLowGQ.vcf.gz
```

## Step 10 — Fill INFO tags and apply MAF / missingness filters

Fill `MAF`, `F_MISSING`, `END`:

```bash
bcftools +fill-tags \
  population.pangenie.nohetGT.noLowGQ.f.vcf.gz \
  -Oz -o population.pangenie.filled.vcf.gz \
  -- -t MAF,F_MISSING,END
tabix -p vcf population.pangenie.filled.vcf.gz

# Sanity check
bcftools query -f '%INFO/MAF\t%INFO/F_MISSING\n' population.pangenie.filled.vcf.gz | head
```

Filter: MAF ≥ 1%, missingness ≤ 20%:

```bash
bcftools view \
  -i 'INFO/MAF >= 0.01 && INFO/F_MISSING <= 0.2' \
  population.pangenie.filled.vcf.gz \
  -Oz -o population.pangenie.filtered.vcf.gz
tabix -p vcf population.pangenie.filtered.vcf.gz
```

## Step 11 — SNPs

```bash
bcftools view -v snps -Oz -o pangenie.snps.raw.vcf.gz population.pangenie.filtered.vcf.gz
tabix -p vcf pangenie.snps.raw.vcf.gz

bcftools view \
  -i 'MAF >= 0.01 && F_MISSING <= 0.1' \
  -Oz -o pangenie.snps.miss10.maf01.vcf.gz \
  pangenie.snps.raw.vcf.gz
tabix -p vcf pangenie.snps.miss10.maf01.vcf.gz
```

### PLINK PCA for SNPs

```bash
/programs/plink-1.9-x86_64-beta7/plink \
  --vcf pangenie.snps.miss10.maf01.vcf.gz \
  --double-id --allow-extra-chr --make-bed \
  --out snps

# LD prune
/programs/plink-1.9-x86_64-beta7/plink \
  --bfile snps --indep-pairwise 50 10 0.2 \
  --out snps.ld

# PCA on LD-pruned SNPs
/programs/plink-1.9-x86_64-beta7/plink \
  --bfile snps --extract snps.ld.prune.in \
  --pca 20 --out snps.pca
```

## Step 12 — SV pipeline

```bash
# Extract all SVs >= 50 bp
bcftools view \
  -i 'abs(strlen(REF) - strlen(ALT)) >= 50' \
  -Oz -o pangenie.sv.gt50bp.vcf.gz \
  population.pangenie.filled.vcf.gz
tabix -p vcf pangenie.sv.gt50bp.vcf.gz

# Apply MAF/missingness
bcftools view \
  -i 'MAF >= 0.01 && F_MISSING <= 0.2' \
  -Oz -o pangenie.sv.gt50bp.maf01.miss20.vcf.gz \
  pangenie.sv.gt50bp.vcf.gz
tabix -p vcf pangenie.sv.gt50bp.maf01.miss20.vcf.gz

# Split by direction
bcftools view -i 'strlen(REF) - strlen(ALT) >= 50' \
  -Oz -o pangenie.del.gt50bp.maf01.miss20.vcf.gz pangenie.sv.gt50bp.maf01.miss20.vcf.gz
tabix -p vcf pangenie.del.gt50bp.maf01.miss20.vcf.gz

bcftools view -i 'strlen(ALT) - strlen(REF) >= 50' \
  -Oz -o pangenie.ins.gt50bp.maf01.miss20.vcf.gz pangenie.sv.gt50bp.maf01.miss20.vcf.gz
tabix -p vcf pangenie.ins.gt50bp.maf01.miss20.vcf.gz
```

### Shorten IDs for PLINK

PLINK barfs on the long PanGenie SV IDs. Replace them with `SV_<row>`:

```bash
for vcf in \
  pangenie.del.gt50bp.maf01.miss20.vcf.gz \
  pangenie.ins.gt50bp.maf01.miss20.vcf.gz
do
  base=$(basename $vcf .vcf.gz)

  bcftools view "$vcf" \
    | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$3="SV_" NR; print}' \
    | bgzip > ${base}.shortID.vcf.gz

  tabix -p vcf ${base}.shortID.vcf.gz
done
```

### PLINK PCA — deletions

```bash
/programs/plink-1.9-x86_64-beta7/plink \
  --vcf pangenie.del.gt50bp.maf01.miss20.shortID.vcf.gz \
  --double-id --allow-extra-chr --make-bed \
  --out del50

/programs/plink-1.9-x86_64-beta7/plink \
  --bfile del50 --pca 20 --out del50.pca
```

### PLINK PCA — insertions

```bash
/programs/plink-1.9-x86_64-beta7/plink \
  --vcf pangenie.ins.gt50bp.maf01.miss20.shortID.vcf.gz \
  --double-id --allow-extra-chr --make-bed \
  --out ins50

/programs/plink-1.9-x86_64-beta7/plink \
  --bfile ins50 --pca 20 --out ins50.pca
```

### PLINK PCA — all SVs combined

```bash
/programs/plink-1.9-x86_64-beta7/plink \
  --vcf pangenie.SV.maf01.miss20.shortID.vcf.gz \
  --double-id --allow-extra-chr --make-bed \
  --out sv.50

/programs/plink-1.9-x86_64-beta7/plink \
  --bfile sv.50 --pca 20 --out sv.50.pca
```

## Outputs used downstream

- `pangenie.snps.miss10.maf01.vcf.gz` — filtered SNP set → `sv_analysis/`
- `pangenie.ins.gt50bp.maf01.miss20.vcf.gz` — filtered insertions → `sv_analysis/`
- `pangenie.del.gt50bp.maf01.miss20.vcf.gz` — filtered deletions → `sv_analysis/`
- `snps.pca.eigenvec`, `del50.pca.eigenvec`, `ins50.pca.eigenvec`, `sv.50.pca.eigenvec` — PCA coordinates for downstream plots

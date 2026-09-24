# Population Structure

## Pipeline overview

1. LD-prune the filtered SNP VCF
2. Run **ADMIXTURE** with cross-validation on SNPs (K = 1–6)
3. LD-prune the filtered SV VCF (500-SV window)
4. Run **ADMIXTURE** with cross-validation on SVs (K = 1–6)
5. Pick the best K per marker set (minimum CV error)
6. Plot in R

## Inputs

- `pangenie.snps.miss10.maf01.vcf.gz` — filtered SNP VCF (from `pangenie_genotyping/`)
- `pangenie.SV.maf01.miss20.shortID.vcf` — filtered SV VCF with short IDs (from `pangenie_genotyping/`)

## Tools used

| Tool | Version / source |
| --- | --- |
| ADMIXTURE | 1.3.0 (`/programs/admixture_linux-1.3.0`) — Alexander and Lange, 2011 |
| PLINK | 1.9 beta 7 (`/programs/plink-1.9-x86_64-beta7`) — Purcell et al., 2007 |


## Step 1 — ADMIXTURE on SNPs

Convert SNPs to PLINK format:

```bash
export PATH=/programs/admixture_linux-1.3.0:$PATH

/programs/plink-1.9-x86_64-beta7/plink \
  --vcf pangenie.snps.miss10.maf01.vcf.gz \
  --double-id \
  --allow-extra-chr \
  --make-bed \
  --out snps
```

LD-prune (200-SNP window, 50-SNP step, r² < 0.5):

```bash
/programs/plink-1.9-x86_64-beta7/plink \
  --bfile snps \
  --indep-pairwise 200 50 0.5 \
  --out snps.ld
```

Extract the pruned set:

```bash
/programs/plink-1.9-x86_64-beta7/plink \
  --bfile snps \
  --extract snps.ld.prune.in \
  --make-bed \
  --out snps.ld.pruned
```

Run ADMIXTURE for K = 1–6 with cross-validation:

```bash
for K in 1 2 3 4 5 6; do
  admixture --cv snps.ld.pruned.bed $K | tee log${K}.out
done
```

## Step 2 — ADMIXTURE on SVs

Convert SVs to PLINK format:

```bash
/programs/plink-1.9-x86_64-beta7/plink \
  --vcf pangenie.SV.maf01.miss20.shortID.vcf \
  --double-id \
  --allow-extra-chr \
  --make-bed \
  --out svs
```

LD-prune:

```bash
/programs/plink-1.9-x86_64-beta7/plink \
  --bfile svs \
  --indep-pairwise 500 50 0.5 \
  --out svs.ld
```

Extract the pruned set:

```bash
/programs/plink-1.9-x86_64-beta7/plink \
  --bfile svs \
  --extract svs.ld.prune.in \
  --make-bed \
  --out svs.ld.pruned
```

Run ADMIXTURE for K = 1–6 with cross-validation:

```bash
for K in 1 2 3 4 5 6; do
  admixture --cv svs.ld.pruned.bed $K | tee sv_log${K}.out
done
```

## Step 3 — Pick the best K

The K with the smallest cross-validation error is ADMIXTURE's estimate of the true number of ancestral populations:

```bash
# SNP-based
grep "CV error" log*.out

# SV-based
grep "CV error" sv_log*.out
```

Outputs:
- `snps.ld.pruned.<K>.Q` — ancestry proportions per isolate at each K
- `snps.ld.pruned.<K>.P` — allele frequencies per ancestral population
- Same for SVs with `svs.ld.pruned.<K>.*`

## Step 4 — Plot in R


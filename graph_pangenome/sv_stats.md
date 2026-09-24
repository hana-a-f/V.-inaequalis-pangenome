# Structural Variant Analysis

This folder documents structural-variant (SV) discovery, characterisation, and preparation for downstream population-level analysis on the pangenome graph built in `pangenome_graph/`.

## Two VCFs, two purposes

**This distinction matters** — reviewers will ask, and it's easy to lose track:

| VCF | Source | Use |
| --- | --- | --- |
| **Raw graph VCF** `venturia-pangenome.raw.vcf.bub.r100k.vcf` | `cactus-pangenome` output (snarl bubbles decomposed by PanGenie's `annotate_vcf.py`) | **SV discovery** — every possible SV the graph found across all 18 assemblies. Answers "what SVs exist in this species?" |
| **PanGenie-filtered VCF** `pangenie.<type>.gt50bp.maf01.miss20.vcf.gz` | PanGenie short-read genotyping → filtered to MAF ≥1% and ≤20% missingness | **Population genetics** — SVs common enough in the surveyed population to test for selection. Answers "which SVs vary in the population and might be under selection?" |

The PanGenie genotyping step itself is not documented here — that's its own pipeline (short-read isolates → PanGenie against the graph → filtered VCF). Note it as a separate folder when you get to it.

## Pipeline overview

1. Decompose the raw graph VCF into biallelic records with PanGenie's `annotate_vcf.py`
2. Normalize + sort with bcftools
3. Filter to SVs ≥50 bp
4. Set canonical IDs (critical for downstream matching)
5. Split into insertions and deletions
6. Extract length + genotypes with `bcftools query` + awk
7. Annotate SVs against the gene GFF with **Ensembl VEP**
8. Extract SV sequences as FASTA
9. Classify SVs by TE content with **RepeatMasker**
10. Extract PanGenie-filtered genotypes in long format (for population analyses)
11. Prepare population-filtered VCFs with canonical IDs

## Length cutoff for SVs

**≥50 bp** 
## Tools used

| Tool | Version / source |
| --- | --- |
| `annotate_vcf.py` | PanGenie toolkit (from PanGenie repo) |
| bcftools | (system) |
| Ensembl VEP | 110.1 (Singularity: `/programs/ensembl-vep-110.1/vep.sif`) |
| RepeatMasker | Via `tetools.sif` (see `repeat_masking/`) |

## Files in this folder

| File | Purpose |
| --- | --- |
| `README.md` | This file |

External scripts referenced (not custom to this project — from other tools' repos):
- `annotate_vcf.py` — from the PanGenie GitHub repo

---

## Step 1 — Decompose the raw graph VCF

The raw `cactus-pangenome` VCF is in "bubble" form — variants are represented as snarls (multi-allelic bubbles in the graph). PanGenie's decomposition script splits these into standard biallelic records that bcftools and other downstream tools can work with:

```bash
python3 annotate_vcf.py \
  -vcf venturia-pangenome.raw.vcf.bub.r100k.vcf \
  -gfa venturia-pangenome.gfa.gz \
  -o mc.pangenie.unfiltered
```

Result:
```
Wrote 1461963 multi-allelic records.
Wrote 2332458 bi-allelic records.
```

## Step 2 — Normalize and sort

Left-align indels against the VI_19_031 reference and produce a sorted, indexed VCF:

```bash
bcftools norm -f VI_19_031_pangenie.fasta -Ou mc.pangenie.unfiltered_biallelic.vcf | \
bcftools sort -Oz --write-index -o mc.pangenie.unfiltered_biallelic.sort.vcf.gz
```

This is the master biallelic-SV file — everything below derives from it.

## Step 3 — Filter to SVs ≥50 bp

```bash
bcftools view \
  -i 'abs(strlen(REF) - strlen(ALT)) >= 50' \
  -Oz \
  -o graph.sv.gt50bp.vcf.gz \
  mc.pangenie.unfiltered_biallelic.sort.vcf.gz

bcftools index graph.sv.gt50bp.vcf.gz
bcftools view -H graph.sv.gt50bp.vcf.gz | wc -l
# 38799
```

## Step 4 — Set IDs (**critical — do not skip**)

Copies the `INFO/ID` field into the VCF's `ID` column so every SV has a stable unique identifier. Every downstream step (VEP annotation, ID-mapping for the TE classification, population selection tests) depends on this:

```bash
bcftools annotate \
  --set-id '%INFO/ID' \
  graph.sv.gt50bp.vcf.gz \
  -Oz -o graph.sv.gt50bp.withID.vcf.gz
```

## Step 5 — Split into insertions and deletions

### Deletions

```bash
bcftools view \
  -i 'strlen(REF) - strlen(ALT) >= 50' \
  -Oz \
  -o graph.del.gt50bp.vcf.gz \
  graph.sv.gt50bp.withID.vcf.gz

bcftools index graph.del.gt50bp.vcf.gz
bcftools view -H graph.del.gt50bp.vcf.gz | wc -l
# 13418
```

### Insertions

```bash
bcftools view \
  -i 'strlen(ALT) - strlen(REF) >= 50' \
  -Oz \
  -o graph.ins.gt50bp.vcf.gz \
  graph.sv.gt50bp.withID.vcf.gz

bcftools index graph.ins.gt50bp.vcf.gz
bcftools view -H graph.ins.gt50bp.vcf.gz | wc -l
# 25381
```

## Step 6 — Extract length and genotypes (raw graph SVs)

`bcftools query` pulls chrom/pos/ref/alt/GT; awk appends the SV length:

```bash
# Deletions
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
  graph.del.gt50bp.vcf.gz | \
awk 'BEGIN{OFS="\t"} {
  ref_len=length($3); alt_len=length($4);
  svlen=(alt_len-1)-(ref_len-1);
  printf "%s\t%s\t%s\t%s", $1,$2,$3,$4;
  for(i=5;i<=NF;i++) printf "\t%s",$i;
  printf "\t%s\n", svlen;
}' > del_sv_withGT.tsv

# Insertions (same shape)
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
  graph.ins.gt50bp.vcf.gz | \
awk 'BEGIN{OFS="\t"} {
  ref_len=length($3); alt_len=length($4);
  svlen=(alt_len-1)-(ref_len-1);
  printf "%s\t%s\t%s\t%s", $1,$2,$3,$4;
  for(i=5;i<=NF;i++) printf "\t%s",$i;
  printf "\t%s\n", svlen;
}' > ins_sv_withGT.tsv
```

## Step 7 — Annotate SVs with VEP against the gene GFF

Cross-references each SV against the merged annotation from `genome_annotation/`; `--distance 2000` includes genes within 2 kb of the SV, `--pick` reports the single best transcript per SV:

```bash
tabix -p gff myAnnot.gff.gz

singularity run -C --bind $PWD --pwd $PWD /programs/ensembl-vep-110.1/vep.sif vep \
  -i graph.del.gt50bp.vcf.gz \
  -o myoutput.SV.del \
  --fork 12 \
  -gff myAnnot.gff.gz \
  -fasta sequences.fa \
  --stats_text --distance 2000 --pick \
  --stats_file mysummary.SV.del.txt

singularity run -C --bind $PWD --pwd $PWD /programs/ensembl-vep-110.1/vep.sif vep \
  -i graph.ins.gt50bp.vcf.gz \
  -o myoutput.SV.ins \
  --fork 12 \
  -gff myAnnot.gff.gz \
  -fasta sequences.fa \
  --stats_text --distance 2000 --pick \
  --stats_file mysummary.SV.ins.txt
```

## Step 8 — Extract graph genotypes (per-sample GT), population genetics

```bash
bcftools query -l graph.ins.gt50bp.vcf.gz > samples.txt

bcftools query \
  -f '%ID\t%CHROM\t%POS[\t%GT]\n' \
  graph.ins.gt50bp.vcf.gz | \
awk 'BEGIN{OFS="\t"} {print $0, "row_"NR}' \
> ins.genotypes.tsv

bcftools query \
  -f '%ID\t%CHROM\t%POS[\t%GT]\n' \
  graph.del.gt50bp.vcf.gz | \
awk 'BEGIN{OFS="\t"} {print $0, "row_"NR}' \
> del.genotypes.tsv
```

## Step 9 — Extract SV sequences as FASTA (for TE classification)

Two versions were produced. The **first** used per-locus allele counting; the **second** uses the canonical `INFO/ID` (Step 4) so every SV sequence matches back to the VCF row unambiguously. **Use the second version for anything that needs matching back.**

### Version 1 — per-allele naming (superseded)

```bash
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\n' \
  graph.sv.gt50bp.vcf.gz | \
awk 'BEGIN{OFS="\t"}{
    ref = $3; alt = $4;
    key = $1"_"$2;
    count[key]++;
    id = key"_"count[key];

    if (length(alt) > length(ref)) {
        seq = substr(alt,2);
        print ">INS_"id"\n"seq
    }
    else if (length(ref) > length(alt)) {
        seq = substr(ref,2);
        print ">DEL_"id"\n"seq
    }
}' > sv_sequences.fa
```

### Version 2 — matched to canonical IDs (**use this one**)

```bash
bcftools query \
  -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
  graph.sv.gt50bp.withID.vcf.gz | \
awk 'BEGIN{OFS="\t"}{
    chrom=$1; pos=$2; id=$3; ref=$4; alt=$5;
    if (length(alt) > length(ref)) {
        n++;
        seq = substr(alt,2);
        print ">SV_"n > "sv_sequences_allmatched.fa";
        print seq   > "sv_sequences_allmatched.fa";
        print "SV_"n"\t"id"\t"chrom"\t"pos"\tINS" > "sv_id_map.tsv";
    }
    else if (length(ref) > length(alt)) {
        n++;
        seq = substr(ref,2);
        print ">SV_"n > "sv_sequences_allmatched.fa";
        print seq   > "sv_sequences_allmatched.fa";
        print "SV_"n"\t"id"\t"chrom"\t"pos"\tDEL" > "sv_id_map.tsv";
    }
}'

# Sanity checks
grep -c ">" sv_sequences_allmatched.fa
wc -l sv_id_map.tsv                       # should match the grep -c above
cut -f1 sv_id_map.tsv | sort | uniq -d    # should print nothing (no duplicate SV_n)
```

## Step 10 — Classify SVs by TE content (RepeatMasker)

Uses the VI_19_031 TE library from `repeat_masking/`. Note the doubled `.classified.classified` suffix on the library — that's real (RepeatModeler classifies the raw output, so both suffixes are present).

**Note:** This is a more compleet classification, used all entries from  RepBase (fungi)

```bash
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker \
  -a -gff -xsmall -pa 10 \
  -lib VI_19_031_consensi.fa.classified.classified \
  sv_sequences.fa \
  > repeatmasker.log 2>&1

singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker \
  -a -gff -xsmall -pa 10 \
  -lib VI_19_031_consensi.fa.classified.classified \
  sv_sequences_allmatched.fa \
  > repeatmasker.log 2>&1
```

Clean the RepeatMasker `.out` for downstream parsing (removes the 3-line header and the trailing `*`):

```bash
tail -n +4 sv_sequences.fa.out \
  | sed 's/^ *//' \
  | tr -s ' ' '\t' \
  | sed 's/\t\*$//' \
  > rm.clean.tsv
```

## Step 11 — Extract PanGenie-filtered SV lengths and per-sample genotypes (population)

Now switches to the **PanGenie-genotyped** VCF (filtered to MAF ≥1%, missingness ≤20%). BED format for downstream intersections:

### Insertions BED

```bash
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\tINS\n' \
  pangenie.ins.gt50bp.maf01.miss20.vcf.gz \
| awk 'BEGIN{OFS="\t"} {
    svlen = length($4) - length($3);
    print $1, $2, $2+1, $5, $6, svlen
}' > ins.bed
```

### Deletions BED

```bash
bcftools query \
  -f '%CHROM\t%POS\t%INFO/END\t%ID\tDEL\n' \
  pangenie.del.gt50bp.maf01.miss20.vcf.gz \
| awk 'BEGIN{OFS="\t"} {
    svlen = $3 - $2;
    print $1, $2, $3, $4, $5, svlen
}' > del.bed
```

### Combined sorted BED

```bash
cat del.bed ins.bed | sort -k1,1 -k2,2n > all_sv.sorted.bed
```

### Long-format sample × SV genotype tables

```bash
bcftools query -l pangenie.del.gt50bp.maf01.miss20.vcf.gz > samples.txt

# Deletions
bcftools query \
  -f '%ID\t%CHROM\t%POS\t%INFO/END[\t%SAMPLE\t%GT]\n' \
  pangenie.del.gt50bp.maf01.miss20.vcf.gz \
| awk 'BEGIN{OFS="\t"} {
    svlen = $4 - $3;
    for (i=5; i<=NF; i+=2) {
        print $1, $2, $3, $4, $(i), $(i+1), svlen
    }
}' > del.genotypes.long.fixed.tsv

# Insertions
bcftools query \
  -f '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE\t%GT]\n' \
  pangenie.ins.gt50bp.maf01.miss20.vcf.gz \
| awk 'BEGIN{OFS="\t"} {
    svlen = length($5) - length($4);
    end = $3 + 1;
    for (i=6; i<=NF; i+=2) {
        print $1, $2, $3, end, $(i), $(i+1), svlen
    }
}' > ins.genotypes.long.fixed.tsv
```

## Step 12 — Prep population-filtered VCFs with canonical IDs (for selection tests)

Same `--set-id '%INFO/ID'` step as Step 4 but on the PanGenie-filtered VCFs. Needed before any selection scan / GWAS-style workflow keyed on SV ID.

```bash
bcftools annotate --set-id '%INFO/ID' pangenie.del.gt50bp.maf01.miss20.vcf.gz    -Oz -o pangenie.del.gt50bp.maf01.miss20.ID.vcf.gz
bcftools annotate --set-id '%INFO/ID' pangenie.ins.gt50bp.maf01.miss20.vcf.gz    -Oz -o pangenie.ins.gt50bp.maf01.miss20.ID.vcf.gz
bcftools annotate --set-id '%INFO/ID' pangenie.snps.miss10.maf01.vcf.gz          -Oz -o pangenie.snps.miss10.maf01.ID.vcf.gz
```

Wide-format sample × GT tables for downstream analyses:

```bash
bcftools query -l pangenie.ins.gt50bp.maf01.miss20.ID.vcf.gz > samples.pop.txt

bcftools query \
  -f '%ID\t%CHROM\t%POS[\t%GT]\n' \
  pangenie.ins.gt50bp.maf01.miss20.ID.vcf.gz \
> ins.pop.genotypes.tsv

bcftools query \
  -f '%ID\t%CHROM\t%POS[\t%GT]\n' \
  pangenie.del.gt50bp.maf01.miss20.ID.vcf.gz \
> del.pop.genotypes.tsv
```

## Step 13 — Overall stats

Count SVs and SNPs in the biallelic + diploid-genotype-called VCF that fed into all downstream analyses:

```bash
bcftools view \
  -i 'abs(strlen(REF) - strlen(ALT)) >= 50' \
  -Oz \
  -o graph.sv.gt50bp.filt.vcf.gz \
  mc.pangenie.biallelic.sort.diplod.vcf.gz

bcftools index graph.sv.gt50bp.filt.vcf.gz

# SV count
bcftools view -H graph.sv.gt50bp.filt.vcf.gz | wc -l

# SNP count
bcftools view -H -v snps mc.pangenie.biallelic.sort.diplod.vcf.gz | wc -l
```

## Headline numbers

| Metric | Count |
| ---: | ---: |
| Total biallelic records (Step 1 decomposition) | 2,332,458 |
| SVs ≥50 bp (Step 3) | 38,799 |
| Deletions ≥50 bp (Step 5) | 13,418 |
| Insertions ≥50 bp (Step 5) | 25,381 |

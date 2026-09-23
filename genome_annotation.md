# Gene Annotation

This folder documents the full gene-annotation pipeline. Each PanSeq-scaffolded, soft-masked assembly (from `ragtag_scaffolding/` + `repeat_masking/` round 2) is annotated with **BRAKER3** using public RNA-Seq and a fungal protein reference; **miniprot** is used in parallel to align known AVR effector homologs so effector candidates are not missed. The two annotation sets are then QC-filtered, merged with **gffread**, reduced to one representative isoform per locus, and formatted for downstream orthology (**OrthoFinder**) and synteny work.

## Pipeline overview

1. Align RNA-Seq reads with **HISAT2** (see `align_bam.sh`)
2. Set up the **BRAKER3** Singularity container
3. Run **BRAKER3** per isolate (RNA-Seq BAMs + fungal protein evidence)
4. Standardise each annotation to the longest isoform per gene with **AGAT**
6. QC — BUSCO on the protein sets
7. Effector-focused annotation with **miniprot** using Rocafort et al. 2023 effector sequences
8. Extract CDS from both BRAKER and miniprot with **gffread**
9. QC the CDS with the custom `gff_qc.mod.py` script
10. Filter by QC flags (frameshifts, internal stops, missing start codons, short ORFs)
11. Merge BRAKER + miniprot annotations per isolate with **gffread -M**
12. Reduce merged annotations to longest isoform per gene (AGAT)
13. Collapse overlapping transcripts to one locus representative with `collapse_longest_isoform_by_locus.py`
14. Final naming — proteins and BEDs renamed to `ISOLATE.gene` format for GENESPACE
15. Final BUSCO on the collapsed protein sets
16. Compute annotation stats — genes per isolate, per chromosome, per Mb

## Tools used

| Tool | Version / source |
| --- | --- |
| HISAT2 | (from `/programs/HISAT2/hisat2.sh`) |
| samtools | (system) |
| BRAKER3 | `docker://teambraker/braker3:latest` (Singularity) |
| Augustus | Bundled inside BRAKER3 container; config copied out to `$PWD/config_${iso}` |
| AGAT | 1.2.0 (Singularity, `/programs/agat-1.2.0/agat.sif`) |
| BUSCO | 5.5.0 (`ascomycota_odb10`) |
| miniprot | 0.13 |
| gffread | 0.12.7 |
| jcvi | (conda env `jcvi`) — deprecated pathway, see end of README |
| GNU parallel | (system) |
| Custom: `gff_qc.mod.py` | In this folder |
| Custom: `collapse_longest_isoform_by_locus.py` | In this folder |


## Files in this folder

| File | Purpose |
| --- | --- |
| `README.md` | This file |
| `align_bam.sh` | HISAT2 alignment of RNA-Seq to every soft-masked assembly |
| `braker.cmds` | Auto-generated BRAKER command list (one per isolate) |
| `busco.prot.cmds` | Auto-generated BUSCO command list (on protein sets) |
| `gff_qc.mod.py` | Custom CDS-QC script (flags frameshifts, internal stops, short ORFs) |
| `collapse_longest_isoform_by_locus.py` | Collapses overlapping longest-isoform transcripts to one representative per locus; also writes a per-locus BED |

## RNA-Seq samples (SRA accessions)

The following 7 paired-end RNA-Seq accessions were downloaded and trimmed prior to alignment (add BioProject in methods) from Rocafort et al 2023

`SRR18280441`, `SRR18280438`, `SRR18280428`, `SRR18280424`, `SRR18280420`, `SRR18280436`, `SRR18280431`

Trimmed reads are expected as `${SRR}_output_forward_paired.fastq` and `${SRR}_output_reverse_paired.fastq`.

---

## Step 1 — HISAT2 alignment (`align_bam.sh`)

Builds a HISAT2 index per genome (if not already built), then aligns each of the 7 RNA-Seq pairs to each genome. Q≥2 filter with `samtools view -q 2`, then sort. Up to 4 genome jobs in parallel; each HISAT2 uses 14 threads, samtools sort uses 6. Idempotent — skips already-completed BAMs.

Full script: `align_bam.sh` (in this folder).

## Step 2 — BRAKER3 setup

```bash
mkdir /workdir/$USER
cd /workdir/$USER
singularity build braker3.sif docker://teambraker/braker3:latest

# Augustus needs a writable config folder outside the container
singularity run --bind $PWD braker3.sif cp -r /opt/Augustus/config $PWD

# Copy protein evidence into the working directory
cp -r /home/khanlab/Hana/braker_102725/evidence/ .
```

Protein evidence: `evidence/Fungi.fa` (add source — OrthoDB v11 fungi or similar).

## Step 3 — Run BRAKER3

Auto-generate one command per isolate (using all 7 RNA-Seq BAMs for that genome) and run 6 in parallel:

```bash
#!/bin/bash
set -euo pipefail

CMD_FILE="braker.cmds"
> "$CMD_FILE"

for f in *_panseq.fasta.masked; do
    ISO=$(basename "$f" _panseq.fasta.masked)
    BAMLIST=$(ls align/SRR*_${ISO}_*.bam | tr $'\n' ',' | sed 's/,$//')

    echo "( ISO=${ISO} && \
mkdir -p braker_${ISO} && \
singularity run --bind \$PWD braker3.sif cp -r /opt/Augustus/config \$PWD/config_${ISO} && \
singularity exec -C --env AUGUSTUS_CONFIG_PATH=\$PWD/config_${ISO} --bind \$PWD --pwd \$PWD braker3.sif braker.pl \
  --genome=\$PWD/${ISO}_panseq.fasta.masked \
  --prot_seq=\$PWD/evidence/Fungi.fa \
  --workingdir=\$PWD/braker_${ISO} \
  --threads 16 \
  --gff3 \
  --fungus \
  --busco_lineage=ascomycota_odb10 \
  --bam=${BAMLIST} \
  > braker_${ISO}.log 2>&1 )" >> braker.cmds
done

parallel -j 6 < braker.cmds
```

Each isolate gets its own `AUGUSTUS_CONFIG_PATH` (`config_${ISO}`) — Augustus writes new species-specific parameters during training, and shared config paths cause races when running in parallel.

## Step 4 — Rename outputs and gather

```bash
# Prefix every file in each braker_* folder with the isolate name
for d in braker_*; do
    iso=${d#braker_}
    for f in "$d"/*; do
        base=$(basename "$f")
        cp "$f" "$d/${iso}_${base}"
    done
done

# Gather all annotations in one folder
mkdir -p isolate_annotations
for d in braker_*; do
    iso=${d#braker_}
    for f in "$d"/*.{gff3,gtf,aa,codingseq}; do
        [ -e "$f" ] || continue
        base=$(basename "$f")
        cp "$f" isolate_annotations/"$base"
    done
done
```

## Step 5 — AGAT: longest isoform + extract proteins (first pass)

For each isolate, keep the longest isoform per gene and extract its protein sequence:

```bash
mkdir -p agat.proteins

for d in braker_*; do
    iso=${d#braker_}
    echo "Processing $iso ..."

    gff=$d/braker.gff3
    genome=${iso}_panseq.fasta.masked
    [[ -f "$gff" && -f "$genome" ]] || { echo "  Skipping $iso"; continue; }

    longest_gff=${iso}_braker_longest.gff3
    longest_prot=${iso}_braker_longest_proteins.fa

    # Keep longest isoform
    singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
        agat_sp_keep_longest_isoform.pl \
        --gff "$gff" \
        -o "$longest_gff"

    # Extract protein sequences
    singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
        agat_sp_extract_sequences.pl \
        -g "$longest_gff" \
        -f "$genome" \
        -p \
        -o "$longest_prot"

    cp "$longest_gff" agat.proteins/
    cp "$longest_prot" agat.proteins/
done
```

## Step 6 — Augustus fallback for VI_18_019, VI_18_033, VI_EUB04

For these three isolates the combined BRAKER output had issues (short/broken gene models), so the Augustus-hints prediction from the same BRAKER run (`*_augustus.hints.gff3`) was used instead — it produced a cleaner, more complete annotation:

```bash
# Longest isoform on Augustus hints GFF
singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
  agat_sp_keep_longest_isoform.pl --gff VI_18_019_augustus.hints.gff3 -o VI_18_019_augustus_longest.gff3
singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
  agat_sp_keep_longest_isoform.pl --gff VI_18_033_augustus.hints.gff3 -o VI_18_033_augustus_longest.gff3
singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
  agat_sp_keep_longest_isoform.pl --gff VI_EUB04_augustus.hints.gff3  -o VI_EUB04_augustus_longest.gff3

# Extract proteins (VI_EUB04 needed wrapping first — long single-line sequences)
singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
  agat_sp_extract_sequences.pl -g VI_18_019_augustus_longest.gff3 -f VI_18_019_panseq.fasta.masked -p -o VI_18_019_augustus_longest_proteins.fa
singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
  agat_sp_extract_sequences.pl -g VI_18_033_augustus_longest.gff3 -f VI_18_033_panseq.fasta.masked -p -o VI_18_033_augustus_longest_proteins.fa

export PATH=/programs/seqkit-0.15.0:$PATH
seqkit seq -w 60 VI_EUB04_panseq.fasta.masked > VI_EUB04_panseq.fasta.masked.wrapped.fasta

singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
  agat_sp_extract_sequences.pl -g VI_EUB04_augustus.hints.gff3      -f VI_EUB04_panseq.fasta.masked.wrapped.fasta -p -o VI_EUB04_augustus_longest_proteins.fa
singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
  agat_sp_extract_sequences.pl -g VI_EUB04_braker_longest.gff3      -f VI_EUB04_panseq.fasta.masked.wrapped.fasta -p -o VI_EUB04_braker_longest_proteins.fa

cp *augustus_longest*                agat.proteins/
cp VI_EUB04_braker_longest_proteins.fa agat.proteins/
```

## Step 7 — BUSCO on protein sets

Auto-generate a BUSCO command per isolate:

```bash
source /programs/miniconda3/bin/activate busco-5.5.0
touch busco.prot.cmds

for f in *_braker_longest_proteins.fa; do
    iso=${f%_braker_longest_proteins.fa}
    echo "busco -f -i $f -l ascomycota_odb10 -o ${iso}.protein.busco -m prot -c 10" >> busco.prot.cmds
done
```

Individual BUSCO runs (Augustus fallback samples):

```bash
busco -f -i VI_18_019_augustus_longest_proteins.fa -l ascomycota_odb10 -o VI_18_019.protein.busco -m prot -c 10
busco -f -i VI_18_033_augustus_longest_proteins.fa -l ascomycota_odb10 -o VI_18_033.protein.busco -m prot -c 10
busco -f -i VI_EUB04_braker_longest_proteins.fa    -l ascomycota_odb10 -o VI_EUB04.protein.busco  -m prot -c 10
busco -f -i VI_1771_2_braker_longest_proteins.fa   -l ascomycota_odb10 -o VI_1771_2.protein.busco -m prot -c 10
```

## Step 8 — miniprot: align AVR effector homologs

Targeted annotation of avirulence-gene homologs (AVR6 family from Rocafort et al. 2024), which are prone to being missed by general-purpose gene predictors.

Refs:
- Rocafort et al., MBE 2024: https://academic.oup.com/mbe/article/41/8/msae164/7728407
- Barragan & Latorre et al., 2024 pipeline: https://github.com/YuSugihara/Barragan_and_Latorre_et_al_2024/tree/main

```bash
export PATH=/programs/miniprot-0.13:$PATH

OUTDIR=/workdir/hf332/miniprot_braker_annotation/miniprot_gff
GENOMEDIR=/workdir/hf332/miniprot_braker_annotation/genomes

mkdir -p "$OUTDIR"

for f in ${GENOMEDIR}/*masked*; do
    base=$(basename "$f")
    isolate=${base%%_panseq*}

    echo "Running miniprot for $isolate"

    miniprot -t 2 \
        -G 3k \
        -P SEC \
        -p 0.3 \
        --outs=0.5 \
        --gff \
        "$f" \
        all.Rocafort.AVR6homologs.fasta \
    | grep -v "^##PAF" \
    > "${OUTDIR}/${isolate}.miniprot.gff3"
done
```

Key flags:
- `-G 3k` — maximum intron length 3 kb (fungal introns are short)
- `-P SEC` — output selenocysteine-safe format
- `-p 0.3 --outs=0.5` — sensitivity for divergent effector homologs

## Step 9 — Extract CDS from both BRAKER and miniprot GFFs

```bash
GFFREAD=/workdir/hf332/gffread-0.12.7.Linux_x86_64/gffread
GENOME_DIR=/workdir/hf332/miniprot_braker_annotation/genomes
BRAKER_GFF_DIR=/workdir/hf332/miniprot_braker_annotation/braker.gff
MINIPROT_GFF_DIR=/workdir/hf332/miniprot_braker_annotation/miniprot_gff

mkdir -p cds

for f in ${GENOME_DIR}/*masked*; do
    base=$(basename "$f")
    isolate=${base%%_panseq*}

    # Find the BRAKER GFF (naming was inconsistent)
    if   [[ -f ${BRAKER_GFF_DIR}/${isolate}.braker.gff3 ]]; then
        BRAKER_GFF=${BRAKER_GFF_DIR}/${isolate}.braker.gff3
    elif [[ -f ${BRAKER_GFF_DIR}/${isolate}.gff3 ]]; then
        BRAKER_GFF=${BRAKER_GFF_DIR}/${isolate}.gff3
    else
        echo "WARNING: No BRAKER GFF for $isolate"
        BRAKER_GFF=""
    fi

    # BRAKER CDS
    [[ -n "$BRAKER_GFF" ]] && \
      ${GFFREAD} -g "$f" -x cds/${isolate}.braker.cds.fa   "$BRAKER_GFF"

    # miniprot CDS
    ${GFFREAD} -g "$f" -x cds/${isolate}.miniprot.cds.fa  ${MINIPROT_GFF_DIR}/${isolate}.miniprot.gff3
done
```

## Step 10 — QC (`gff_qc.mod.py`)

Custom script checks each CDS against its GFF and flags:
- CDS length not a multiple of 3 (`not_multiple_of_3`)
- Internal stop codons (`stop_codon_in_cds`)
- Missing start codon (`no_start_codon`)
- Shorter than 150 / 180 / 195 nt thresholds

```bash
QC_DIR=/workdir/hf332/miniprot_braker_annotation/qc
CDS_DIR=/workdir/hf332/miniprot_braker_annotation/cds
BRAKER_GFF_DIR=/workdir/hf332/miniprot_braker_annotation/braker.gff
MINIPROT_GFF_DIR=/workdir/hf332/miniprot_braker_annotation/miniprot_gff

mkdir -p "$QC_DIR"

for cds in ${CDS_DIR}/*.braker.cds.fa; do
    isolate=$(basename "$cds" .braker.cds.fa)

    # find the right BRAKER GFF
    if   [[ -f ${BRAKER_GFF_DIR}/${isolate}.braker.gff3 ]]; then BRAKER_GFF=${BRAKER_GFF_DIR}/${isolate}.braker.gff3
    elif [[ -f ${BRAKER_GFF_DIR}/${isolate}.gff3 ]]; then        BRAKER_GFF=${BRAKER_GFF_DIR}/${isolate}.gff3
    else continue
    fi

    # BRAKER QC
    python gff_qc.mod.py \
        "$BRAKER_GFF" \
        "${CDS_DIR}/${isolate}.braker.cds.fa" \
        150,180,195 \
        10,25,50 \
        1> ${QC_DIR}/${isolate}.braker_qc.gff3 \
        2> ${QC_DIR}/${isolate}.braker_qc.txt

    # miniprot QC — strip header comments first
    grep -v '^#' ${MINIPROT_GFF_DIR}/${isolate}.miniprot.gff3 \
      > ${QC_DIR}/${isolate}.miniprot.noheader.gff3

    python gff_qc.mod.py \
        ${QC_DIR}/${isolate}.miniprot.noheader.gff3 \
        "${CDS_DIR}/${isolate}.miniprot.cds.fa" \
        150,180,195 \
        10,25,50 \
        1> ${QC_DIR}/${isolate}.miniprot_qc.gff3 \
        2> ${QC_DIR}/${isolate}.miniprot_qc.txt
done
```

## Step 11 — Filter by QC flags

Drop transcripts with frame issues, internal stops, missing start codons, or under 150 nt:

```bash
FILTERED_DIR=/workdir/hf332/miniprot_braker_annotation/qc_filtered
mkdir -p "$FILTERED_DIR"

for f in ${QC_DIR}/*.braker_qc.gff3; do
    isolate=$(basename "$f" .braker_qc.gff3)
    grep -v 'not_multiple_of_3' "$f" | \
    grep -v 'stop_codon_in_cds' | \
    grep -v 'no_start_codon'   | \
    grep -v 'shorter_than_150nt' | \
    cut -f 1-9 \
    > ${FILTERED_DIR}/${isolate}.braker_qc.filtered.gff3
done

for f in ${QC_DIR}/*.miniprot_qc.gff3; do
    isolate=$(basename "$f" .miniprot_qc.gff3)
    grep -v 'not_multiple_of_3' "$f" | \
    grep -v 'stop_codon_in_cds' | \
    grep -v 'no_start_codon'   | \
    grep -v 'shorter_than_150nt' | \
    cut -f 1-9 \
    > ${FILTERED_DIR}/${isolate}.miniprot_qc.filtered.gff3
done
```

## Step 12 — Merge BRAKER + miniprot, extract CDS + protein

`gffread -M -K` merges overlapping features and keeps unique isoforms; `--force-exons` adds explicit exon lines where missing:

```bash
MERGED_DIR=/workdir/hf332/miniprot_braker_annotation/merged
mkdir -p "$MERGED_DIR"

for gff in ${QC_FILT_DIR}/*.braker_qc.filtered.gff3; do
    isolate=$(basename "$gff" .braker_qc.filtered.gff3)
    GENOME=${GENOME_DIR}/${isolate}_panseq.fasta.masked
    BRAKER_GFF=${QC_FILT_DIR}/${isolate}.braker_qc.filtered.gff3
    MINIPROT_GFF=${QC_FILT_DIR}/${isolate}.miniprot_qc.filtered.gff3

    ${GFFREAD} \
        --sort-alpha \
        --force-exons \
        -M \
        -K \
        "$BRAKER_GFF" \
        "$MINIPROT_GFF" \
        > ${MERGED_DIR}/${isolate}.merged.gff3

    ${GFFREAD} -g "$GENOME" -x ${MERGED_DIR}/${isolate}.merged.cds.fa      ${MERGED_DIR}/${isolate}.merged.gff3
    ${GFFREAD} -g "$GENOME" -y ${MERGED_DIR}/${isolate}.merged.protein.fa  ${MERGED_DIR}/${isolate}.merged.gff3
done
```

## Step 13 — Longest isoform per gene on the merged GFF (AGAT)

```bash
OUT_DIR=/workdir/hf332/miniprot_braker_annotation/longest_isoform
mkdir -p $OUT_DIR

for gff in $GFF_DIR/*.merged.gff3; do
    isolate=$(basename "$gff" .merged.gff3)
    genome=$GENOME_DIR/${isolate}_panseq.fasta.masked

    # Convert merged GFF3 to GTF
    /workdir/hf332/gffread-0.12.7.Linux_x86_64/gffread \
        -E -T "$gff" -o $OUT_DIR/${isolate}.gtf

    # Longest isoform
    singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
        agat_sp_keep_longest_isoform.pl \
        --gff $OUT_DIR/${isolate}.gtf \
        -o $OUT_DIR/${isolate}.longest.gtf

    # Proteins
    singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
        agat_sp_extract_sequences.pl \
        -g $OUT_DIR/${isolate}.longest.gtf \
        -f "$genome" \
        -p \
        -o $OUT_DIR/${isolate}.longest.prot.fa

    # CDS
    singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
        agat_sp_extract_sequences.pl \
        -g $OUT_DIR/${isolate}.longest.gtf \
        -f "$genome" \
        -t cds \
        -o $OUT_DIR/${isolate}.longest.cds.fa
done
```

## Step 14 — VI_EUB04 special-case: Augustus annotation

The Augustus-hints prediction was used instead of the merged BRAKER/miniprot output:

```bash
singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
    agat_sp_keep_longest_isoform.pl --gff VI_EUB04.augustus.gtf -o VI_EUB04.longest.gtf

fold VI_EUB04_panseq.fasta.masked > VI_EUB04_panseq.fasta.fold.masked

singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
    agat_sp_extract_sequences.pl -g VI_EUB04.longest.gtf -f VI_EUB04_panseq.fasta.fold.masked -p -o VI_EUB04.longest.prot.fa
singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif \
    agat_sp_extract_sequences.pl -g VI_EUB04.longest.gtf -f VI_EUB04_panseq.fasta.fold.masked -t cds -o VI_EUB04.longest.cds.fa

# BUSCO check
busco -f -i VI_EUB04.longest.prot.fa -l ascomycota_odb10 -o VI_EUB04.FOLD -m prot -c 12
```

## Step 15 — Collapse by locus

Longest-isoform per gene still left nested/overlapping transcripts. `collapse_longest_isoform_by_locus.py` merges these into one representative per locus and writes a BED file for each isolate (the BED is what feeds OrthoFinder / synteny):

```bash
OUTDIR=/workdir/hf332/miniprot_braker_annotation/collapse
mkdir -p $OUTDIR
cd $OUTDIR

for gff in $GFF_DIR/*.merged.gff3; do
    iso=$(basename "$gff" .merged.gff3)
    prot=$LONGEST/${iso}.longest.prot.fa

    python collapse_longest_isoform_by_locus.py "$gff" "$prot" "$iso"
done
```

## Step 16 — Final naming: `ISOLATE.gene` format

Rename protein and BED headers to the format OrthoFinder/synteny expect (`VI_1771_2.gene0001` etc.):

```bash
mkdir -p named.collapsed.proteins

for f in *.locus_rep.prot.fa; do
    iso=${f%.locus_rep.prot.fa}
    out="named.collapsed.proteins/${iso}.fa"

    awk -v ISO="$iso" '
        /^>/ {
            hdr = substr($0, 2)
            split(hdr, a, /[ \t]/)
            gene = a[1]
            sub(/\.t[0-9]+$/, "", gene)
            print ">" ISO "." gene
            next
        }
        { print }
    ' "$f" > "$out"
done

mkdir -p named.collapsed.bed

for bed in *.locus_rep.bed; do
    iso=$(basename "$bed" .locus_rep.bed)
    out="named.collapsed.bed/${iso}.bed"

    awk -v ISO="$iso" 'BEGIN{OFS="\t"}
    {
        gene = $4
        sub(/\.t[0-9]+$/, "", gene)
        $4 = ISO "." gene
        print
    }' "$bed" > "$out"
done

# Strip the PanSeq `SAMPLE#H#` prefix from chromosome names in the BEDs
for bed in *.bed; do
    awk 'BEGIN{OFS="\t"}
    {
        sub(/^.*#.*#/, "", $1)
        print $1, $2, $3, $4
    }' "$bed" > tmp && mv tmp "$bed"
done
```

## Step 17 — Final BUSCO on collapsed protein sets

```bash
LONG=/workdir/hf332/miniprot_braker_annotation/longest_isoform
OUTDIR=/workdir/hf332/miniprot_braker_annotation/busco
source /programs/miniconda3/bin/activate busco-5.5.0
mkdir -p "$OUTDIR"

for long in "$LONG"/*.longest.prot.fa; do
    isolate=$(basename "$long" .longest.prot.fa)
    busco -f -i "$long" -l ascomycota_odb10 -o "$isolate" --out_path "$OUTDIR" -m prot -c 12
done
```

## Step 18 — Annotation stats

### Genes per isolate

```bash
echo -e "isolate\ttotal_genes" > genes_per_isolate.tsv

for gtf in *.longest.gtf; do
    isolate=${gtf%.longest.gtf}

    n=$(awk '
        $3=="gene" {
            for(i=1;i<=NF;i++)
                if($i ~ /gene_id/) {
                    gsub(/"|;/,"",$i)
                    split($i,a,"=")
                    print a[2]
                }
        }' "$gtf" | sort -u | wc -l)

    # Fallback: some BRAKER-only GTFs don't have `gene` features
    if [[ "$n" -eq 0 ]]; then
        n=$(awk '
            $3~/transcript|mRNA/ {
                for(i=1;i<=NF;i++)
                    if($i ~ /gene_id/) {
                        gsub(/"|;/,"",$i)
                        split($i,a,"=")
                        print a[2]
                    }
            }' "$gtf" | sort -u | wc -l)
    fi

    echo -e "${isolate}\t${n}" >> genes_per_isolate.tsv
done
```

### Genes per chromosome (example, single sample)

```bash
awk '$3=="gene"{count[$1]++}
     END{for (c in count) print c "\t" count[c]}' \
    VI_EUB04.longest.gtf | sort
```

### Chromosome sizes + genes per Mb

```bash
echo -e "isolate\tchromosome\tlength_bp" > chrom_sizes.tsv

for fa in *_panseq.fasta.masked; do
    isolate=${fa%_panseq.fasta.masked}
    samtools faidx "$fa"
    awk -v iso="$isolate" '{print iso "\t" $1 "\t" $2}' "${fa}.fai" >> chrom_sizes.tsv
done

echo -e "isolate\tchromosome\tgenes\tlength_mb\tgenes_per_mb" > genes_per_mb_per_chromosome.tsv

awk '
    NR==FNR { g[$1,$2]=$3; next }
    {
        iso=$1; chr=$2; len=$3
        key=iso SUBSEP chr
        if (key in g) {
            mb=len/1e6
            printf "%s\t%s\t%d\t%.3f\t%.2f\n", iso, chr, g[key], mb, g[key]/mb
        }
    }
' genes_per_chromosome_all.tsv chrom_sizes.tsv >> genes_per_mb_per_chromosome.tsv
```

## Sanity checks

```bash
grep -c "^>" named.collapsed.proteins/*.fa   # protein counts per isolate
grep -c "^>" named.cds/*.fa                  # CDS counts per isolate
wc -l bed/*.bed                              # BED line counts
```

---

## Deprecated pathway (kept for reference)

An earlier version of the naming step used `jcvi` to make BEDs from the *uncollapsed* longest-isoform GTFs. This produced correct BEDs but did not collapse nested/overlapping transcripts to a single locus, which caused downstream problems in OrthoFinder. Superseded by Step 15 (`collapse_longest_isoform_by_locus.py`).

```bash
source /programs/miniconda3/bin/activate jcvi

mkdir -p bed
for f in *.longest.gtf; do
    isolate=${f%.longest.gtf}
    python -m jcvi.formats.gff bed --type=transcript --key=ID "$f" -o "bed/${isolate}.bed"
done

# Normalize chromosome names (strip PanSeq prefix)
for bed in bed/*.bed; do
    awk 'BEGIN{OFS="\t"}
    { sub(/^.*#.*#/, "", $1); print $1, $2, $3, $4 }' "$bed" > tmp && mv tmp "$bed"
done

# Rename proteins / CDS / BEDs to ISOLATE.gene format (same awk block as Step 16)
```

# Hi-C Scaffolding — VI_19_031

This folder documents the Hi-C scaffolding pipeline used to produce the **chromosome-level reference assembly for VI_19_031**. This is the only isolate in the project with Hi-C data, so it is also the only one scaffolded to chromosome level; every other isolate remains at contig level (see `genome_assembly/`). The finished VI_19_031 assembly is used downstream as the anchor for RagTag scaffolding of the other genomes.

## Pipeline overview

1. QC the Hi-C reads (fastp + FastQC + assembly-stats)
2. Clean the HiFi reads (fastp)
3. Estimate genome size / heterozygosity (jellyfish + GenomeScope)
4. Assemble contigs with hifiasm
5. Remove mitochondrial contigs (MitoHiFi)
6. Purge duplicate haplotype copies
7. Remove small unplaced contigs (BLAST-identified as mito)
8. Align Hi-C reads with Juicer
9. Scaffold with YaHS (tested multiple `-r` bin sizes; used 50000)
10. Manual curation in Juicebox → post-review with `juicer post`
11. Rename to panSeq header format and final QC (Merqury, BUSCO, telomere check)

## Sample notes

- **Only VI_19_031 has Hi-C data**, so only this sample was scaffolded to chromosome level.
- **GenomeScope report:** http://genomescope.org/analysis.php?code=AyxYQC3Q1hS7DqT4Ujhw
- **BUSCO after purge_dups:** C:98.0% [S:97.9%, D:0.1%], F:0.6%, M:1.4%, n:1706 (ascomycota_odb10)
- The Hi-C reads were provided by Delaware already demultiplexed and filtered.

## Dead ends / things tried but not used

- **Initial MitoHiFi reference `CP118466.1`** was auto-selected by `findMitoReference.py` as *Ascochyta lentis* — wrong species. Replaced with the correct *V. inaequalis* reference `PV785816.1` (see `genome_assembly/`).
- **HapHic** scaffolding — produced a more fragmented assembly than YaHS.
- **3D-DNA** scaffolding — also produced a more fragmented assembly than YaHS.
- **YaHS `-r` (bin size) parameter sweep:** tested 25000, 50000, 100000, 150000 — settled on **50000**.

## Tools used

| Tool | Version |
| --- | --- |
| fastp | 0.23.4 |
| FastQC | (system) |
| assembly-stats | `/programs/assembly-stats` |
| jellyfish | 2.3.0 |
| GenomeScope | (web version) |
| hifiasm | 0.19.9 |
| MitoHiFi | 3.0.0 (Singularity) |
| seqkit | 0.15.0 |
| purge_dups | 1.2.6 |
| minimap2 | 2.17 |
| BUSCO | 5.5.0 (`ascomycota_odb10`) |
| BWA | (system) |
| samtools | (system) |
| Juicer | 2 (`/programs/juicer2`) |
| juicer_tools | 1.9.9 |
| YaHS | 1.2.2 |
| Merqury | 1.3 (Singularity) |
| RagTag | (add version) |

## Files in this folder

| File | What it does |
| --- | --- |
| `README.md` | This file |
| `telomeres_4.0.py` | Telomere detection script (referenced in QC step) |
| `re-check_commands.sh` | Full commands for the second-iteration re-alignment (referenced) |

---

## Step 1 — QC the Hi-C reads

Remove duplicates and adapters, then check quality:

```bash
/programs/fastp-0.23.4/fastp \
  -i VI-19-031-1_L3_137A37.R1.fastq.gz \
  -I VI-19-031-1_L3_137A37.R2.fastq.gz \
  -o VI_19_031_R1.fastq.gz \
  -O VI_19_031_R2.fastq.gz \
  --detect_adapter_for_pe

fastqc VI_19_031_R1.fastq.gz VI_19_031_R2.fastq.gz

export PATH=/programs/assembly-stats:$PATH
assembly-stats VI_19_031_R1.fastq
assembly-stats VI_19_031_R2.fastq
```

## Step 2 — Clean the HiFi reads

Delaware pre-demultiplexed and pre-filtered the PacBio reads, but a fastp pass was still run:

```bash
/programs/fastp-0.23.4/fastp -i VI_19_031.fastq.gz -o filtered_VI_19_031.hifi.fastq.gz
zcat -dk filtered_VI_19_031.hifi.fastq.gz > filtered_VI_19_031.hifi.fastq
```

## Step 3 — Jellyfish + GenomeScope

```bash
export PATH=/programs/jellyfish-2.3.0/bin:$PATH
jellyfish count -m 21 -s 100M -t 10 -C -o VI_19_031_hifi_kmers_21.jf filtered_VI_19_031.hifi.fastq
jellyfish histo VI_19_031_hifi_kmers_21.jf > VI_19_031_hifi_kmers_21.histo
```

GenomeScope output: http://genomescope.org/analysis.php?code=AyxYQC3Q1hS7DqT4Ujhw

## Step 4 — Assemble contigs with hifiasm

```bash
hifiasm -o VI_19_031_hifi.asm -t 32 --telo-m GGGTTA filtered_VI_19_031.hifi.fastq.gz

# Convert primary contig GFA to FASTA
awk '/^S/{print ">"$2;print $3}' VI_19_031_hifi.asm.bp.p_ctg.gfa \
  > VI_19_031_hifi.telo.asm.bp.p_ctg.fa

# Assembly stats
export PATH=/programs/assembly-stats:$PATH
export PATH=/programs/seqkit-0.15.0:$PATH
assembly-stats VI_19_031_hifi.telo.asm.bp.p_ctg.fa
seqkit stats VI_19_031_hifi.telo.asm.bp.p_ctg.fa
```

## Step 5 — Remove mitochondrial contigs

Correct *V. inaequalis* mito reference (`PV785816.1`) is used in `genome_assembly/`; that same reference should be used here for consistency. The commands below reflect what was run.

```bash
# NOTE: the auto-picked reference was wrong (Ascochyta lentis) — see Dead Ends
singularity run --bind $PWD --pwd $PWD /programs/mitohifi-3.0.0/mitohifi.sif \
  mitohifi.py -c VI_19_031_hifi.telo.asm.bp.p_ctg.fa \
  -f CP118466.1.fasta -g CP118466.1.gb \
  -t 24 -a fungi

# Build a list of contigs to remove
find . -type d -name "ptg*.annotation" | sed 's|.*/||' | sed 's/\.annotation$//' > mt_contigs.txt

# Remove them (started ~380 contigs → 25 after)
export PATH=/programs/seqkit-0.15.0:$PATH
seqkit grep -v -f mt_contigs.txt VI_19_031_hifi.telo.asm.bp.p_ctg.fa \
  > VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.fa
```

## Step 6 — Purge duplicate haplotypes 

```bash
export PYTHONPATH=/programs/purge_dups-1.2.6/lib/python3.9/site-packages:/programs/purge_dups-1.2.6/lib/python3.9/site-packages/runner-0.0.0-py3.9.egg
export PATH=/programs/purge_dups-1.2.6/bin:/programs/purge_dups-1.2.6/scripts/:$PATH
export PATH=/programs/minimap2-2.17:$PATH

ls /home/khanlab/venturia_pangenome/VI_19_031.fastq.gz > fastq.fofn
pd_config.py -l . -n config.json VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.fa fastq.fofn
run_purge_dups.py -p bash config.json /programs/purge_dups-1.2.6/bin VI_19_031

# BUSCO check
source /programs/miniconda3/bin/activate busco-5.5.0
busco -f -i VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.fa \
  -l ascomycota_odb10 -o purged_busco -m geno -c 24
# → C:98.0% [S:97.9%, D:0.1%], F:0.6%, M:1.4%, n:1706
```

## Step 7 — Remove small mito-hit contigs

Small contigs kept producing Hi-C hits; BLAST (NCBI GUI) confirmed they were all mitochondrial.

```bash
export PATH=/programs/seqkit-0.15.0:$PATH

# Extract contigs < 50 kb
seqkit seq -m 1 -M 50000 VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.fa > small_contigs.fasta

# (BLASTed via NCBI GUI — all confirmed mito)

# Build the list and remove
grep '^>' small_contigs.fasta | sed 's/^>//' > mito_contigs.txt
seqkit grep -v -f mito_contigs.txt VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.fa \
  > VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.noMT.fa
```

## Step 8 — Set up Juicer inputs

Juicer needs four inputs — genome (`-z`), restriction sites (`-y`), chromosome sizes (`-p`), and the paired fastqs.

```bash
# -z: wrapped reference + BWA index
mkdir references/
export PATH=/programs/seqkit-0.15.0:$PATH
seqkit seq -w 60 VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.noMT.fa \
  > VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.noMT.wrapped.fa
bwa index VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.noMT.wrapped.fa

# -y: MboI restriction site positions
mkdir restriction_sites/
python /programs/juicer2/misc/generate_site_positions.py MboI \
  VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.noMT.wrapped.fa \
  ../references/VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.noMT.wrapped.fa

# -p: chromosome sizes
mkdir chr_size/
samtools faidx VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.noMT.wrapped.fa
cut -f1,2 VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.noMT.wrapped.fa.fai \
  > VI_19_031.chrom.sizes.txt

# Paired Hi-C fastqs
mkdir fastq/
cp /home/hf332/hc_wrk/fastq/*.fastq .
```

## Step 9 — Run Juicer

```bash
manage_slurm new `hostname -s`

# Copy juicer to a workdir and set it up for SLURM
cp -r /programs/juicer2 /workdir/hf332/
cd /workdir/hf332/juicer2/
ln -s SLURM/scripts ./
cp collisions.awk scripts/

# Move all Juicer inputs into the juicer2 dir (must be clean — do not reuse from prior runs)
cp -r ../references/ .
cp -r ../chr_size/ .
cp -r ../restriction_sites/ .

screen -S juicer_run
export JuicerDir=/workdir/hf332/juicer2
mkdir logs

/workdir/hf332/juicer2/scripts/juicer.sh \
  -D $JuicerDir \
  -g VI_19_031 \
  -s MboI \
  -z $JuicerDir/references/VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.noMT.wrapped.fa \
  -y $JuicerDir/restriction_sites/VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.noMT.wrapped.fa_MboI.txt \
  -p $JuicerDir/chr_size/VI_19_031.chrom.sizes.txt \
  -q regular \
  -l regular \
  -C 500000 \
  -t 64 \
  --assembly \
  > $JuicerDir/logs/juicer_run_61024.log 2>&1

# Build the .hic file for Juicebox
java -Xmx16g -jar /workdir/hf332/juicer2/scripts/juicer_tools.jar pre \
  -q 30 \
  /workdir/hf332/juicer2/aligned/merged_nodups.txt \
  /workdir/hf332/juicer2/aligned/VI_19_031.hic \
  /workdir/hf332/juicer2/chr_size/VI_19_031.chrom.sizes.txt
```

## Step 10 — Scaffold with YaHS (first iteration)

Tested several `-r` bin sizes; **50000** was chosen.

```bash
export PATH=/programs/yahs-1.2.2:$PATH

# Name-sort the Juicer BAM for YaHS
samtools sort -n -@ 16 hic.sorted.bam -o hic.namesort.bam

# Bin-size sweep
for r in 25000 50000 100000 150000; do
  out="yahs_r${r}"
  yahs \
    -r ${r} \
    -e GATC \
    --telo-motif TTAGGG \
    -o ${out} \
    VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.noMT.fa \
    hic.namesort.bam \
    > ${out}.log 2>&1
done
```

Build the `.hic` for Juicebox review (using r=50000):

```bash
wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar

juicer pre -a -o out_JBAT_50000 yahs_r50000.bin yahs_r50000_scaffolds_final.agp \
  VI_19_031_hifi.telo.nuc.asm.bp.p_ctg.purged.noMT.fa.fai \
  > out_JBAT.50000.log 2>&1

grep "PRE_C_SIZE" out_JBAT.50000.log | awk '{print $2"\t"$3}' > out_JBAT_50000.chrom.sizes

java -jar -Xmx32G juicer_tools.1.9.9_jcuda.0.8.jar pre \
  out_JBAT_50000.txt \
  out_JBAT_50000.hic.part \
  out_JBAT_50000.chrom.sizes \
  && mv out_JBAT_50000.hic.part out_JBAT_50000.part.hic
```

## Step 11 — Manual curation in Juicebox

The `.hic` was opened in Juicebox for manual review; the review-corrected `.assembly` file was then fed to `juicer post` to apply the edits (primarily fixing telomere sequences).


```bash
samtools sort -n -@ 16 merged_dedup.sorted.bam -o merged_dedup.sorted.bam

yahs \
  -r 50000 \
  -e GATC \
  --telo-motif TTAGGG \
  -o yahs.r.50000 \
  ragtag.scaffold.wrapped.fasta \
  merged_dedup.sorted.bam \
  > yahs.r.50000.log 2>&1

# Rebuild .hic for the second-iteration Juicebox review
juicer pre -a -o out_JBAT.r_50000 yahs.r.50000.bin yahs.r.50000_scaffolds_final.agp \
  ragtag.scaffold.wrapped.fasta.fai \
  > out_JBAT.r.50000.log 2>&1

grep "PRE_C_SIZE" out_JBAT.r.50000.log | awk '{print $2"\t"$3}' > out_JBAT.r.50000.log.chrom.sizes

java -jar -Xmx32G juicer_tools.1.9.9_jcuda.0.8.jar pre \
  out_JBAT.r_50000.txt \
  out_JBAT.r_50000.hic.part \
  out_JBAT.r.50000.log.chrom.sizes \
  && mv out_JBAT.r_50000.hic.part out_JBAT.r._50000.hic
```

## Step 12 — Rename and final QC

Rename contigs to the PanSeq header format (`SAMPLE#HAPLOTYPE#chrNN`) for downstream pangenome work:

```bash
awk '/^>/{i++; printf(">VI_19_031#1#chr%02d\n", i); next} {print}' \
  VI_19_031_assembly.fa > VI_19_031_assembly.panseq.fa

export PATH=/programs/seqkit-0.15.0:$PATH
seqkit stats -a VI_19_031_assembly.panseq.fa
```

Merqury k-mer QC:

```bash
export MERQURY=/merqury
singularity run -B /workdir/$USER --pwd /workdir/$USER /programs/merqury-1.3/merqury.sif \
  meryl k=21 threads=16 count VI_19_031.fastq output VI_19_031.meryl
singularity run -B /workdir/$USER --pwd /workdir/$USER /programs/merqury-1.3/merqury.sif \
  merqury.sh VI_19_031.meryl VI_19_031_assembly.panseq.fa VI_19_031.FINAL
```

Final BUSCO:

```bash
source /programs/miniconda3/bin/activate busco-5.5.0
busco -f -i VI_19_031_assembly.panseq.fa -l ascomycota_odb10 -o Final_busco -m geno -c 24
```

Telomere check: see `telomeres_4.0.py`.

Chromosome sizes:

```bash
samtools faidx VI_19_031_assembly.panseq.fa
cut -f1,2 VI_19_031_assembly.panseq.fa.fai > chrom.sizes
```

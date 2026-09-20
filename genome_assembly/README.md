# Genome Assembly

This step takes raw PacBio HiFi reads for each *Venturia inaequalis* isolate and produces cleaned, haplotig-purged nuclear assemblies ready for downstream analysis (repeat masking, annotation, pangenome build).

## Overview of the pipeline

1. Summarise raw reads with **seqkit**
2. Estimate genome size and heterozygosity with **jellyfish** (k=21)
3. Assemble each genome with **hifiasm** (fungal telomere motif `GGGTTA`)
4. Convert the primary contig graph (`.gfa`) to FASTA
5. Assemble the mitochondrial genome and identify mito contigs with **MitoHiFi**
6. Remove mitochondrial contigs from the nuclear assemblies with **seqkit grep**
7. Remove duplicate haplotype copies with **purge_haplotigs**
8. Run a final QC with **seqkit stats**

## Tools used

- [seqkit](https://bioinf.shenwei.me/seqkit/) 0.15.0
- [jellyfish](https://github.com/gmarcais/Jellyfish) 2.3.0
- [hifiasm](https://github.com/chhoward/hifiasm) 0.19.9
- [MitoHiFi](https://github.com/marcelauliano/MitoHiFi) 3.0.0 (via Singularity)
- [minimap2](https://github.com/lh3/minimap2) + [samtools](https://www.htslib.org/)
- [purge_haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) (dev branch)

## Sample notes

- **16 isolates** were assembled fresh from HiFi reads with hifiasm.
- **VI_19_031** was pre-assembled to chromosome level in prior work; it joins the sample set from downstream steps onwards (no hifiasm output in this folder).
- **Four samples had contamination** and their NCBI pre-screened assemblies were used for further analysis instead of resubmitting:
  - VI_19_004
  - VI_18_019
  - VI_18_033
  - VI_1797_2
- **VI_19_011** returned no MitoHiFi hits.

## Files in this folder

| File | What it does |
| --- | --- |
| `jellyfish.sh` | k-mer counting + histogram for each sample (input for genome-size estimation) |
| `hifiasm.sh` | Loop over all `VI_*.fastq` files and assemble each with hifiasm |
| `run_all_mitohifi.sh` | Run MitoHiFi in contigs mode (`-c`) on every hifiasm assembly |
| `contigs_ids.txt` | Curated list of mitochondrial contig IDs to filter out |
| `README.md` | This file |

---

## Step 1 — Read stats

```bash
export PATH=/programs/seqkit-0.15.0:$PATH
seqkit stats -a *.fastq > summary.tsv
```

## Step 2 — Jellyfish (k-mer counting)

See `jellyfish.sh`. Runs `jellyfish count` and `jellyfish histo` per sample; use the `.histo` files with GenomeScope to estimate genome size and heterozygosity.

## Step 3 — Assemble with hifiasm

See `hifiasm.sh`. Each sample gets its own subfolder; hifiasm is run with `--telo-m GGGTTA` (fungal telomere motif) and 32 threads.

## Step 4 — Convert primary GFA to FASTA

```bash
awk '/^S/{print ">"$2;print $3}' VI_1771_2_hifi.asm.bp.p_ctg.gfa > VI_1771_2_hifi.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' VI_1797_9_hifi.asm.bp.p_ctg.gfa > VI_1797_9_hifi.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' VI_18_030_hifi.asm.bp.p_ctg.gfa > VI_18_030_hifi.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' VI_18_037_hifi.asm.bp.p_ctg.gfa > VI_18_037_hifi.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' VI_18_043_hifi.asm.bp.p_ctg.gfa > VI_18_043_hifi.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' VI_19_011_hifi.asm.bp.p_ctg.gfa > VI_19_011_hifi.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' VI_EU104_hifi.asm.bp.p_ctg.gfa > VI_EU104_hifi.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' VI_EU160_hifi.asm.bp.p_ctg.gfa > VI_EU160_hifi.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' VI_EU301_hifi.asm.bp.p_ctg.gfa > VI_EU301_hifi.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' VI_EU302_hifi.asm.bp.p_ctg.gfa > VI_EU302_hifi.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' VI_EU413_hifi.asm.bp.p_ctg.gfa > VI_EU413_hifi.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' VI_EUNL19_hifi.asm.bp.p_ctg.gfa > VI_EUNL19_hifi.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' VI_EUNL24_hifi.asm.bp.p_ctg.gfa > VI_EUNL24_hifi.asm.bp.p_ctg.fa
```

## Step 5 — MitoHiFi

See `run_all_mitohifi.sh`. Uses MitoHiFi in contigs mode (`-c`) with the *V. inaequalis* mitochondrial reference `PV785816.1`. Each sample gets its own output subfolder.

For VI_19_011 (which had no MitoHiFi hits from the primary reference), the alternative reference `PP484585.4` was used:

```bash
singularity run --bind $PWD --pwd $PWD /programs/mitohifi-3.0.0/mitohifi.sif \
  mitohifi.py -c VI_19_011_curated.fasta \
  -f mito_reference/PP484585.4.fasta \
  -g mito_reference/PP484585.4.gb \
  -t 24 -a fungi
```

## Step 6 — Remove mitochondrial contigs from nuclear assemblies

Build the mito-contig ID list (per-sample version, uses MitoHiFi output folders named `ptg*.annotation`):

```bash
find . -type d -name "ptg*.annotation" | sed 's|.*/||' | sed 's/\.annotation$//' > mt_contigs.txt
```

Filter out mito contigs with seqkit grep. The curated `contigs_ids.txt` was used across all samples:

```bash
export PATH=/programs/seqkit-0.15.0:$PATH

seqkit grep -v -f contigs_filtering/contigs_ids.txt VI_1771_2_hifi.asm.bp.p_ctg.fa > VI_1771_2_hifi.asm.bp.p_ctg.nuc.fa
seqkit grep -v -f contigs_filtering/contigs_ids.txt VI_1797_2_hifi.asm.bp.p_ctg.fa > VI_1797_2_hifi.asm.bp.p_ctg.nuc.fa
seqkit grep -v -f contigs_filtering/contigs_ids.txt VI_1797_9_hifi.asm.bp.p_ctg.fa > VI_1797_9_hifi.asm.bp.p_ctg.nuc.fa
seqkit grep -v -f contigs_filtering/contigs_ids.txt VI_18_019_hifi.asm.bp.p_ctg.fa > VI_18_019_hifi.asm.bp.p_ctg.nuc.fa
seqkit grep -v -f contigs_filtering/contigs_ids.txt VI_18_030_hifi.asm.bp.p_ctg.fa > VI_18_030_hifi.asm.bp.p_ctg.nuc.fa
seqkit grep -v -f contigs_filtering/contigs_ids.txt VI_18_037_hifi.asm.bp.p_ctg.fa > VI_18_037_hifi.asm.bp.p_ctg.nuc.fa
seqkit grep -v -f contigs_filtering/contigs_ids.txt VI_18_043_hifi.asm.bp.p_ctg.fa > VI_18_043_hifi.asm.bp.p_ctg.nuc.fa
seqkit grep -v -f contigs_filtering/contigs_ids.txt VI_EU104_hifi.asm.bp.p_ctg.fa > VI_EU104_hifi.asm.bp.p_ctg.nuc.fa
seqkit grep -v -f contigs_filtering/contigs_ids.txt VI_EU160_hifi.asm.bp.p_ctg.fa > VI_EU160_hifi.asm.bp.p_ctg.nuc.fa
seqkit grep -v -f contigs_filtering/contigs_ids.txt VI_EU302_hifi.asm.bp.p_ctg.fa > VI_EU302_hifi.asm.bp.p_ctg.nuc.fa
seqkit grep -v -f contigs_filtering/contigs_ids.txt VI_EU413_hifi.asm.bp.p_ctg.fa > VI_EU413_hifi.asm.bp.p_ctg.nuc.fa
seqkit grep -v -f contigs_filtering/contigs_ids.txt VI_EUNL19_hifi.asm.bp.p_ctg.fa > VI_EUNL19_hifi.asm.bp.p_ctg.nuc.fa
seqkit grep -v -f contigs_filtering/contigs_ids.txt VI_EUNL24_hifi.asm.bp.p_ctg.fa > VI_EUNL24_hifi.asm.bp.p_ctg.nuc.fa
```

Same treatment for the pre-assembled chromosome-level sample:

```bash
seqkit grep -v -f contigs_ids.txt VI_19_031.chr.ref.fasta > VI_19_031.chr.ref.nuc.fasta
```

## Step 7 — Purge haplotigs

One-time install (only needs to be done once, ever):

```bash
git clone https://mroachawri@bitbucket.org/mroachawri/purge_haplotigs.git
cd purge_haplotigs/
git checkout dev
mkdir -p ~/bin
ln -s "$(pwd)/bin/purge_haplotigs" ~/bin/purge_haplotigs
which purge_haplotigs
```

### 7a. Map long reads to each assembly (minimap2 + samtools)

```bash
minimap2 -x map-pb -t 24 -a your_assembly.fa your_reads.fastq | \
  samtools view -b -o aligned.bam -
samtools sort -@ 8 -o aligned.sorted.bam aligned.bam
samtools index aligned.sorted.bam
```

### 7b. Coverage histogram — pick cutoffs from the PNGs

```bash
for asm in assemblies/*_hifi.asm.bp.p_ctg.nuc.fa; do
  sample=$(basename "$asm" _hifi.asm.bp.p_ctg.nuc.fa)
  bam="${sample}.sorted.bam"

  if [ ! -f "$bam" ]; then
    echo "Missing BAM for $sample — expected: $bam"
    continue
  fi

  echo "Running histogram step for $sample"
  purge_haplotigs hist \
    -b "$bam" \
    -g "$asm" \
    -t 8
done
```

### 7c. Apply per-sample coverage cutoffs (`-l`, `-m`, `-h` read from the histograms above)

```bash
purge_haplotigs cov -i VI_1771_2.sorted.bam.200.gencov -l 1  -m 1   -h 25  -o VI_1771_2_coverage.csv
purge_haplotigs cov -i VI_1797_2.sorted.bam.200.gencov -l 2  -m 13  -h 25  -o VI_1797_2_coverage.csv
purge_haplotigs cov -i VI_1797_9.sorted.bam.200.gencov -l 5  -m 20  -h 65  -o VI_1797_9_coverage.csv
purge_haplotigs cov -i VI_18_019.sorted.bam.200.gencov -l 5  -m 37  -h 65  -o VI_18_019_coverage.csv
purge_haplotigs cov -i VI_18_030.sorted.bam.200.gencov -l 80 -m 105 -h 175 -o VI_18_030_coverage.csv
purge_haplotigs cov -i VI_18_033.sorted.bam.200.gencov -l 1  -m 2   -h 20  -o VI_18_033_coverage.csv
purge_haplotigs cov -i VI_18_037.sorted.bam.200.gencov -l 20 -m 25  -h 70  -o VI_18_037_coverage.csv
purge_haplotigs cov -i VI_18_043.sorted.bam.200.gencov -l 35 -m 45  -h 105 -o VI_18_043_coverage.csv
purge_haplotigs cov -i VI_19_004.sorted.bam.200.gencov -l 5  -m 10  -h 50  -o VI_19_004_coverage.csv
purge_haplotigs cov -i VI_19_011.sorted.bam.200.gencov -l 10 -m 35  -h 80  -o VI_19_011_coverage.csv
purge_haplotigs cov -i VI_EU104.sorted.bam.200.gencov  -l 5  -m 15  -h 55  -o VI_EU104_coverage.csv
purge_haplotigs cov -i VI_EU160.sorted.bam.200.gencov  -l 2  -m 5   -h 40  -o VI_EU160_coverage.csv
purge_haplotigs cov -i VI_EU301.sorted.bam.200.gencov  -l 25 -m 35  -h 100 -o VI_EU301_coverage.csv
purge_haplotigs cov -i VI_EU302.sorted.bam.200.gencov  -l 10 -m 25  -h 100 -o VI_EU302_coverage.csv
purge_haplotigs cov -i VI_EU413.sorted.bam.200.gencov  -l 10 -m 20  -h 75  -o VI_EU413_coverage.csv
purge_haplotigs cov -i VI_EUNL19.sorted.bam.200.gencov -l 5  -m 25  -h 75  -o VI_EUNL19_coverage.csv
purge_haplotigs cov -i VI_EUNL24.sorted.bam.200.gencov -l 5  -m 20  -h 90  -o VI_EUNL24_coverage.csv
```

### 7d. Run the purging pipeline

```bash
purge_haplotigs purge -g assemblies/VI_1771_2_hifi.asm.bp.p_ctg.nuc.fa -c VI_1771_2_coverage.csv -t 8 -d -b VI_1771_2.sorted.bam -o VI_1771_2_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_1797_2_hifi.asm.bp.p_ctg.nuc.fa -c VI_1797_2_coverage.csv -t 8 -d -b VI_1797_2.sorted.bam -o VI_1797_2_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_1797_9_hifi.asm.bp.p_ctg.nuc.fa -c VI_1797_9_coverage.csv -t 8 -d -b VI_1797_9.sorted.bam -o VI_1797_9_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_18_019_hifi.asm.bp.p_ctg.nuc.fa -c VI_18_019_coverage.csv -t 8 -d -b VI_18_019.sorted.bam -o VI_18_019_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_18_030_hifi.asm.bp.p_ctg.nuc.fa -c VI_18_030_coverage.csv -t 8 -d -b VI_18_030.sorted.bam -o VI_18_030_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_18_033_hifi.asm.bp.p_ctg.nuc.fa -c VI_18_033_coverage.csv -t 8 -d -b VI_18_033.sorted.bam -o VI_18_033_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_18_037_hifi.asm.bp.p_ctg.nuc.fa -c VI_18_037_coverage.csv -t 8 -d -b VI_18_037.sorted.bam -o VI_18_037_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_18_043_hifi.asm.bp.p_ctg.nuc.fa -c VI_18_043_coverage.csv -t 8 -d -b VI_18_043.sorted.bam -o VI_18_043_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_19_004_hifi.asm.bp.p_ctg.nuc.fa -c VI_19_004_coverage.csv -t 8 -d -b VI_19_004.sorted.bam -o VI_19_004_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_19_011_hifi.asm.bp.p_ctg.nuc.fa -c VI_19_011_coverage.csv -t 8 -d -b VI_19_011.sorted.bam -o VI_19_011_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_EU104_hifi.asm.bp.p_ctg.nuc.fa -c VI_EU104_coverage.csv -t 8 -d -b VI_EU104.sorted.bam -o VI_EU104_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_EU160_hifi.asm.bp.p_ctg.nuc.fa -c VI_EU160_coverage.csv -t 8 -d -b VI_EU160.sorted.bam -o VI_EU160_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_EU301_hifi.asm.bp.p_ctg.nuc.fa -c VI_EU301_coverage.csv -t 8 -d -b VI_EU301.sorted.bam -o VI_EU301_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_EU302_hifi.asm.bp.p_ctg.nuc.fa -c VI_EU302_coverage.csv -t 8 -d -b VI_EU302.sorted.bam -o VI_EU302_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_EU413_hifi.asm.bp.p_ctg.nuc.fa -c VI_EU413_coverage.csv -t 8 -d -b VI_EU413.sorted.bam -o VI_EU413_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_EUNL19_hifi.asm.bp.p_ctg.nuc.fa -c VI_EUNL19_coverage.csv -t 8 -d -b VI_EUNL19.sorted.bam -o VI_EUNL19_hifi.asm.bp.p_ctg.nuc.purge
purge_haplotigs purge -g assemblies/VI_EUNL24_hifi.asm.bp.p_ctg.nuc.fa -c VI_EUNL24_coverage.csv -t 8 -d -b VI_EUNL24.sorted.bam -o VI_EUNL24_hifi.asm.bp.p_ctg.nuc.purge
```

## Step 8 — Final QC

```bash
for f in *.fasta; do
  seqkit stats "$f"
done

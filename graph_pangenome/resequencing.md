# Population Short-Read Data

This folder documents acquisition and QC of the public short-read sequencing data used for PanGenie genotyping against the pangenome graph. Reads come from four BioProjects on the SRA, are downloaded via ENA, filtered to the final set of **136 isolates** used in this study, and then adapter- and quality-trimmed with Trimmomatic. The trimmed paired reads feed directly into `pangenie_genotyping/`.

## BioProjects used

| BioProject | # SRRs kept |
| --- | ---: | --- |
| PRJNA354841 | 24 | 
| PRJNA407103 | 78 | 
| PRJNA817384 | 2 | 
| PRJNA962118 | 32 | 


## Pipeline overview

1. Download every SRR run's paired FASTQs from ENA via `wget`
2. Filter to the final list (`keep_srrs.txt`) — copy only the kept SRRs into a working directory
3. Quality/adapter trim with **Trimmomatic** — produces the `_1P.fastq.gz` / `_2P.fastq.gz` paired outputs

## Tools used

| Tool | Version / source |
| --- | --- |
| wget | (system) |
| Trimmomatic | 0.39 (`/programs/trimmomatic/trimmomatic-0.39.jar`) |
| Adapters | TruSeq3-PE.fa (`/programs/trimmomatic/adapters/`) |
| GNU parallel | (system) |

## Files in this folder

| File | Purpose |
| --- | --- |
| `README.md` | This file |
| `keep_srrs.txt` | The 105 SRR IDs used in the analysis |
| `download_urls.txt` | All raw ENA download URLs (one `wget -nc <url>` line per FASTQ) |
| `filter_reads.sh` | Copy only the SRRs in `keep_srrs.txt` from the download directory into the working directory |
| `trim_reads.sh` | Trimmomatic paired-end quality/adapter trimming in parallel |

---

## Step 1 — Download from ENA

All reads were downloaded via `wget -nc` (no-clobber, so re-running is idempotent) from ENA's FTP mirror of the SRA. Example lines:

```bash
# PRJNA354841 (25 runs — kept 24)
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR518/003/SRR5183043/SRR5183043_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR518/003/SRR5183043/SRR5183043_2.fastq.gz
# ...

# North-American isolates (SRR24315*)
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR243/046/SRR24315646/SRR24315646_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR243/046/SRR24315646/SRR24315646_2.fastq.gz
# ...
```

The complete list of URLs is in `download_urls.txt` in this folder. For BioProjects PRJNA407103 and the NA batch, all runs were downloaded and then filtered (Step 2) — this was simpler than trying to filter at download time.

Run everything:

```bash
bash download_urls.txt
```

## Step 2 — Filter to the keep list

The keep list (`keep_srrs.txt`) is the definitive record of which isolates are in the study, listed in supplemental table 16.
```bash
#!/bin/bash
FASTQ_DIR="/workdir/hf332/fastq"
OUT_DIR="/workdir/hf332/filtered_fastq"
mkdir -p "$OUT_DIR"

while read -r srr; do
  for pair in 1 2; do
    file=$(find "$FASTQ_DIR" -type f -name "${srr}_${pair}.fastq.gz")
    if [[ -n "$file" ]]; then
      cp "$file" "$OUT_DIR/"
    else
      echo "Missing: ${srr}_${pair}.fastq.gz"
    fi
  done
done < keep_srrs.txt
```

The "Missing:" warnings surface any SRR listed in `keep_srrs.txt` but not present in the download directory — investigate those before proceeding.

## Step 3 — Quality and adapter trim (Trimmomatic)

Adapter clipping (TruSeq3-PE) + leading/trailing Q20 trim + 4-bp sliding-window Q15 + average-quality Q20 + minimum length 25:

```bash
#!/bin/bash
set -euo pipefail

RAW_DIR="/workdir/hf332/sra/filtered_fastq"
OUT_DIR="/workdir/hf332/sra/trimmed_fastq"
LOG_DIR="/workdir/hf332/sra/logs"
THREADS=8    # threads per sample
JOBS=8       # number of samples in parallel

mkdir -p "$OUT_DIR" "$LOG_DIR"

TRIMMOMATIC_JAR="/programs/trimmomatic/trimmomatic-0.39.jar"
ADAPTERS="/programs/trimmomatic/adapters/TruSeq3-PE.fa"

export RAW_DIR OUT_DIR LOG_DIR THREADS TRIMMOMATIC_JAR ADAPTERS

# Auto-build sample list from R1 filenames
ls ${RAW_DIR}/*_1.fastq.gz | sed 's/.*\///;s/_1.fastq.gz//' > sample_list.txt

cat sample_list.txt | parallel -j $JOBS "
  echo 'Processing {1}...'
  java -jar $TRIMMOMATIC_JAR PE -threads $THREADS -phred33 \
    ${RAW_DIR}/{1}_1.fastq.gz ${RAW_DIR}/{1}_2.fastq.gz \
    ${OUT_DIR}/{1}_1P.fastq.gz ${OUT_DIR}/{1}_1U.fastq.gz \
    ${OUT_DIR}/{1}_2P.fastq.gz ${OUT_DIR}/{1}_2U.fastq.gz \
    ILLUMINACLIP:$ADAPTERS:2:30:10 \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:25 \
    &> ${LOG_DIR}/{1}_trimmomatic.log
  echo 'Finished {1}'
"
```

Trimmomatic output files per sample:
- `{SRR}_1P.fastq.gz` — R1 paired
- `{SRR}_1U.fastq.gz` — R1 unpaired (mate was dropped)
- `{SRR}_2P.fastq.gz` — R2 paired
- `{SRR}_2U.fastq.gz` — R2 unpaired

The `_1P` and `_2P` paired files are what `pangenie_genotyping/` uses. 

## What goes to PanGenie

After trimming, `pangenie_genotyping/` picks up from here:
1. Interleaves each isolate's `_1P`/`_2P` with SeqFu → `{SRR}.interleaved.fastq.gz`
2. Runs PanGenie against the graph 

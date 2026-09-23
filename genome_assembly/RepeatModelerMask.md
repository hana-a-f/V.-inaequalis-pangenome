#Repeat Modeling & Masking

This folder documents the repeat-modeling and soft-masking pipeline. Each assembly gets its own *de novo* repeat library from **RepeatModeler**, which is then used by **RepeatMasker** to soft-mask that assembly (repeats lower-cased in place with `-xsmall`, so downstream tools can still see the sequence but distinguish it as repetitive).

**Repeat masking is run twice in this project.** RagTag ignores soft-masked bases during scaffolding, so it was fine to feed it soft-masked FASTAs — but I hadn't planned in advance that I'd want the *final* repeat GFF3s to line up with the *scaffolded* PanSeq-header assemblies that feed into the pangenome graph. That's why round 2 exists.

1. **First round** — on the *purged, mito-removed contig-level assemblies* (this folder). The masked FASTAs then feed into RagTag scaffolding.
2. **Second round** — on the *scaffolded PanSeq-formatted FASTAs* (see `ragtag_scaffolding/` — a note on the second-round commands). The second-round GFF3s stay in sync with the exact FASTA headers used to build the pangenome graph.

## Tools used

| Tool | Version |
| --- | --- |
| Dfam TE Tools container | `docker://dfam/tetools:latest` (bundles RepeatModeler2, RepeatMasker, RECON, RepeatScout, TRF, LTR_retriever, ninja, etc.) |
| Dfam database | 3.8 |
| GNU parallel | (system) |
| Parsing-RepeatMasker-Outputs (`parseRM.pl`) | From https://github.com/4ureliek/Parsing-RepeatMasker-Outputs |

Key flags used:
- **RepeatModeler**: `-LTRStruct` (adds structural LTR-retrotransposon discovery)
- **RepeatMasker**: `-a -gff -xsmall` (write `.align` file, write GFF3, **soft-mask** with lowercase rather than replace with N)

## Samples (17 total)

VI_18_043 was originally masked but was later dropped from the analysis, so it is not listed here.

| Isolate | Input assembly (round 1) |
| --- | --- |
| VI_19_031 | `VI_19_031_CHR_final.fa` (chromosome-level, from `hic_assembly/`) |
| VI_1771_2 | `VI_1771_2_hifi.asm.bp.p_ctg.nuc.purge.fasta` |
| VI_1797_2 | `VI_1797_2_hifi.asm.bp.p_ctg.nuc.purge.fasta` |
| VI_1797_9 | `VI_1797_9_hifi.asm.bp.p_ctg.nuc.purge.fasta` |
| VI_18_019 | `VI_18_019_hifi.asm.bp.p_ctg.nuc.purge.fasta` |
| VI_18_030 | `VI_18_030_hifi.asm.bp.p_ctg.nuc.purge.fasta` |
| VI_18_033 | `VI_18_033_curated.fasta` |
| VI_18_037 | `VI_18_037_hifi.asm.bp.p_ctg.nuc.purge.fasta` |
| VI_19_004 | `VI_19_004_hifi.asm.bp.p_ctg.nuc.purge.fasta` |
| VI_19_011 | `VI_19_011_curated.fasta` |
| VI_EU104 | `VI_EU104_hifi.asm.bp.p_ctg.nuc.purge.fasta` |
| VI_EU160 | `VI_EU160_hifi.asm.bp.p_ctg.nuc.purge.fasta` |
| VI_EU301 | `VI_EU301_hifi.asm.bp.p_ctg.nuc.purge.fasta` |
| VI_EU302 | `VI_EU302_hifi.asm.bp.p_ctg.nuc.purge.fasta` |
| VI_EU413 | `VI_EU413_hifi.asm.bp.p_ctg.nuc.purge.fasta` |
| VI_EUNL19 | `VI_EUNL19_hifi.asm.bp.p_ctg.nuc.purge.fasta` |
| VI_EUNL24 | `VI_EUNL24_hifi.asm.bp.p_ctg.nuc.purge.fasta` |

## Files in this folder

| File | Purpose |
| --- | --- |
| `README.md` | This file |
| `database.txt` | Auto-generated RepeatModeler command list (one line per isolate) |
| `mask.txt` | Auto-generated RepeatMasker command list |

---

## Step 1 — Pull the Dfam TE Tools container

```bash
singularity pull tetools.sif docker://dfam/tetools:latest
```

## Step 2 — Build a BLAST database per isolate

```bash
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_19_031  VI_19_031_CHR_final.fa
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_1771_2  VI_1771_2_hifi.asm.bp.p_ctg.nuc.purge.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_1797_2  VI_1797_2_hifi.asm.bp.p_ctg.nuc.purge.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_1797_9  VI_1797_9_hifi.asm.bp.p_ctg.nuc.purge.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_18_019  VI_18_019_hifi.asm.bp.p_ctg.nuc.purge.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_18_030  VI_18_030_hifi.asm.bp.p_ctg.nuc.purge.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_18_033  VI_18_033_curated.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_18_037  VI_18_037_hifi.asm.bp.p_ctg.nuc.purge.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_19_004  VI_19_004_hifi.asm.bp.p_ctg.nuc.purge.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_19_011  VI_19_011_hifi.asm.bp.p_ctg.nuc.purge.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_EU104   VI_EU104_hifi.asm.bp.p_ctg.nuc.purge.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_EU160   VI_EU160_hifi.asm.bp.p_ctg.nuc.purge.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_EU301   VI_EU301_hifi.asm.bp.p_ctg.nuc.purge.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_EU302   VI_EU302_hifi.asm.bp.p_ctg.nuc.purge.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_EU413   VI_EU413_hifi.asm.bp.p_ctg.nuc.purge.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_EUNL19  VI_EUNL19_hifi.asm.bp.p_ctg.nuc.purge.fasta
singularity run --bind $PWD --pwd $PWD ./tetools.sif BuildDatabase -name VI_EUNL24  VI_EUNL24_hifi.asm.bp.p_ctg.nuc.purge.fasta
```

## Step 3 — RepeatModeler (build a *de novo* TE library per isolate)

Make an output directory per isolate:

```bash
for iso in VI_19_031 VI_1771_2 VI_1797_2 VI_1797_9 VI_18_019 VI_18_030 VI_18_033 \
           VI_18_037 VI_19_004 VI_19_011 VI_EU104 VI_EU160 VI_EU301 \
           VI_EU302 VI_EU413 VI_EUNL19 VI_EUNL24; do
    mkdir -p "${iso}_out"
done
```

Run RepeatModeler per isolate — added to `database.txt` and executed 3 in parallel (`parallel -j 3 < database.txt`):

```bash
singularity run --bind $PWD --pwd $PWD/VI_EU302_out  ./tetools.sif RepeatModeler -database ../VI_EU302  -threads 10 -LTRStruct 2>&1 | tee VI_EU302_out/VI_EU302_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_EU104_out  ./tetools.sif RepeatModeler -database ../VI_EU104  -threads 10 -LTRStruct 2>&1 | tee VI_EU104_out/VI_EU104_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_EU160_out  ./tetools.sif RepeatModeler -database ../VI_EU160  -threads 20 -LTRStruct 2>&1 | tee VI_EU160_out/VI_EU160_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_EU413_out  ./tetools.sif RepeatModeler -database ../VI_EU413  -threads 20 -LTRStruct 2>&1 | tee VI_EU413_out/VI_EU413_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_EUNL24_out ./tetools.sif RepeatModeler -database ../VI_EUNL24 -threads 20 -LTRStruct 2>&1 | tee VI_EUNL24_out/VI_EUNL24_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_EUNL19_out ./tetools.sif RepeatModeler -database ../VI_EUNL19 -threads 10 -LTRStruct 2>&1 | tee VI_EUNL19_out/VI_EUNL19_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_EU301_out  ./tetools.sif RepeatModeler -database ../VI_EU301  -threads 10 -LTRStruct 2>&1 | tee VI_EU301_out/VI_EU301_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_18_030_out ./tetools.sif RepeatModeler -database ../VI_18_030 -threads 20 -LTRStruct 2>&1 | tee VI_18_030_out/VI_18_030_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_19_031_out ./tetools.sif RepeatModeler -database ../VI_19_031 -threads 20 -LTRStruct 2>&1 | tee VI_19_031_out/VI_19_031_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_19_004_out ./tetools.sif RepeatModeler -database ../VI_19_004 -threads 20 -LTRStruct 2>&1 | tee VI_19_004_out/VI_19_004_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_19_011_out ./tetools.sif RepeatModeler -database ../VI_19_011 -threads 20 -LTRStruct 2>&1 | tee VI_19_011_out/VI_19_011_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_18_019_out ./tetools.sif RepeatModeler -database ../VI_18_019 -threads 10 -LTRStruct 2>&1 | tee VI_18_019_out/VI_18_019_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_18_033_out ./tetools.sif RepeatModeler -database ../VI_18_033 -threads 10 -LTRStruct 2>&1 | tee VI_18_033_out/VI_18_033_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_18_037_out ./tetools.sif RepeatModeler -database ../VI_18_037 -threads 10 -LTRStruct 2>&1 | tee VI_18_037_out/VI_18_037_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_1771_2_out ./tetools.sif RepeatModeler -database ../VI_1771_2 -threads 10 -LTRStruct 2>&1 | tee VI_1771_2_out/VI_1771_2_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_1797_2_out ./tetools.sif RepeatModeler -database ../VI_1797_2 -threads 10 -LTRStruct 2>&1 | tee VI_1797_2_out/VI_1797_2_repeatmodeler.log
singularity run --bind $PWD --pwd $PWD/VI_1797_9_out ./tetools.sif RepeatModeler -database ../VI_1797_9 -threads 10 -LTRStruct 2>&1 | tee VI_1797_9_out/VI_1797_9_repeatmodeler.log

parallel -j 3 < database.txt
```

## Step 4 — Rename and gather RepeatModeler consensus libraries

Each RepeatModeler run outputs a generic `consensi.fa.classified`. Rename to `<isolate>_consensi.fa.classified` and gather all 17 in one directory:

```bash
# Inside each *_out folder, rename the consensus file with the isolate name.
# Example for one sample:
cp consensi.fa.classified VI_EU302_consensi.fa.classified

# Then move them all to a single working directory:
cp VI_EU302_consensi.fa.classified /workdir/hf332/consensi
# ...repeat for every isolate
```

## Step 5 — RepeatMasker (round 1 — on contig-level assemblies)

`-xsmall` is essential — it soft-masks (lowercase) rather than replacing repeats with `N`, so downstream tools can still read the sequence. `-a` writes the alignment file needed for the TE landscape step. `-gff` writes a GFF3 annotation of repeats.

Copy `tetools.sif` into the consensi directory and run from there (RepeatMasker expects the library and the genome in the working directory):

```bash
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_EU302_consensi.fa.classified   VI_EU302_hifi.asm.bp.p_ctg.nuc.purge.fasta   > VI_EU302_repeatmasker.log   2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_EU104_consensi.fa.classified   VI_EU104_hifi.asm.bp.p_ctg.nuc.purge.fasta   > VI_EU104_repeatmasker.log   2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_EU160_consensi.fa.classified   VI_EU160_hifi.asm.bp.p_ctg.nuc.purge.fasta   > VI_EU160_repeatmasker.log   2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_EU413_consensi.fa.classified   VI_EU413_hifi.asm.bp.p_ctg.nuc.purge.fasta   > VI_EU413_repeatmasker.log   2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_EUNL24_consensi.fa.classified  VI_EUNL24_hifi.asm.bp.p_ctg.nuc.purge.fasta  > VI_EUNL24_repeatmasker.log  2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_EUNL19_consensi.fa.classified  VI_EUNL19_hifi.asm.bp.p_ctg.nuc.purge.fasta  > VI_EUNL19_repeatmasker.log  2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_EU301_consensi.fa.classified   VI_EU301_hifi.asm.bp.p_ctg.nuc.purge.fasta   > VI_EU301_repeatmasker.log   2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_18_030_consensi.fa.classified  VI_18_030_hifi.asm.bp.p_ctg.nuc.purge.fasta  > VI_18_030_repeatmasker.log  2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 5  -lib VI_19_031_consensi.fa.classified  VI_19_031_CHR_final.fa                       > VI_19_031_repeatmasker.log  2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_18_019_consensi.fa.classified  VI_18_019_hifi.asm.bp.p_ctg.nuc.purge.fasta  > VI_18_019_repeatmasker.log  2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_1771_2_consensi.fa.classified  VI_1771_2_hifi.asm.bp.p_ctg.nuc.purge.fasta  > VI_1771_2_repeatmasker.log  2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_1797_2_consensi.fa.classified  VI_1797_2_hifi.asm.bp.p_ctg.nuc.purge.fasta  > VI_1797_2_repeatmasker.log  2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_1797_9_consensi.fa.classified  VI_1797_9_hifi.asm.bp.p_ctg.nuc.purge.fasta  > VI_1797_9_repeatmasker.log  2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_19_004_consensi.fa.classified  VI_19_004_hifi.asm.bp.p_ctg.nuc.purge.fasta  > VI_19_004_repeatmasker.log  2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_18_033_consensi.fa.classified  VI_18_033_curated.fasta                      > VI_18_033_repeatmasker.log  2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_18_037_consensi.fa.classified  VI_18_037_hifi.asm.bp.p_ctg.nuc.purge.fasta  > VI_18_037_repeatmasker.log  2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_19_011_consensi.fa              VI_19_011_curated.fasta                      > VI_19_011_repeatmasker.log  2>&1
```

Written to `mask.txt` and dispatched:

```bash
parallel -j 1 < mask.txt
```

### Reruns

- **VI_19_031** was re-masked using `VI_19_031.chr.ref.nuc.fasta` (the canonical PanSeq-header version) to keep GFF3 coordinates in sync with the reference used downstream.
- **VI_EUB04** was masked using VI_EU302's library (question: is this the correct pairing? Flagged below).

```bash
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_19_031_consensi.fa            VI_19_031.chr.ref.nuc.fasta      > VI_19_031_repeatmasker.log  2>&1
singularity run --bind $PWD --pwd $PWD ./tetools.sif RepeatMasker -a -gff -xsmall -pa 10 -lib VI_EU302_consensi.fa.classified  VI_EUB04_panseq.fasta.masked     > VI_EUB04_repeatmasker.log   2>&1
```

## Step 6 — RepeatMasker (round 2 — on scaffolded PanSeq assemblies)

After RagTag scaffolding, RepeatMasker is re-run on the `*_panseq.fasta` files so that the resulting GFF3 coordinates match the exact headers of the FASTAs fed into the pangenome graph build. See `ragtag_scaffolding/README.md` for the auto-generated `repeat_masker.cmds` loop that dispatches these round-2 jobs.

## Step 7 — TE landscape plots (parseRM.pl)

### Install parseRM

```bash
git clone https://github.com/4ureliek/Parsing-RepeatMasker-Outputs.git
cd Parsing-RepeatMasker-Outputs
chmod +x parseRM.pl
perl parseRM.pl --help
```

### Flags used

| Flag | Meaning |
| --- | --- |
| `--land 50,1` | Landscape plot: 50 bins × 1% divergence increments |
| `--parse` | Parse `.align` file to extract data |
| `--fa` | Output consensus sequences in FASTA |
| `-n` (`--nrem`) | Remove nested repeats from the analysis |
| `-p` | (parse — same behaviour as `--parse`) |
| `-i` | Input `.align` file |

Repeat divergence is measured as the percentage of substitutions from the consensus sequence. Higher divergence generally indicates older repeats; low divergence suggests recent insertions.

### Run on every `.align` file

```bash
touch parseRM.cmds

for f in *.align; do
    echo "/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl \
      --land 50,1 --parse --fa -n -p -i $f" >> parseRM.cmds
done

parallel -j 8 < parseRM.cmds
```

Or, listed explicitly per sample (15 isolates — VI_19_011 and VI_19_031 to be added when their `.align` files are updated):

```bash
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_1771_2_hifi.asm.bp.p_ctg.nuc.purge.fasta.align
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_1797_2_hifi.asm.bp.p_ctg.nuc.purge.fasta.align
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_1797_9_hifi.asm.bp.p_ctg.nuc.purge.fasta.align
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_18_019_hifi.asm.bp.p_ctg.nuc.purge.fasta.align
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_18_030_hifi.asm.bp.p_ctg.nuc.purge.fasta.align
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_18_033_curated.fasta.align
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_18_037_hifi.asm.bp.p_ctg.nuc.purge.fasta.align
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_19_004_hifi.asm.bp.p_ctg.nuc.purge.fasta.align
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_EU104_hifi.asm.bp.p_ctg.nuc.purge.fasta.align
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_EU160_hifi.asm.bp.p_ctg.nuc.purge.fasta.align
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_EU301_hifi.asm.bp.p_ctg.nuc.purge.fasta.align
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_EU302_hifi.asm.bp.p_ctg.nuc.purge.fasta.align
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_EU413_hifi.asm.bp.p_ctg.nuc.purge.fasta.align
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_EUNL19_hifi.asm.bp.p_ctg.nuc.purge.fasta.align
/workdir/hf332/parse_RM/Parsing-RepeatMasker-Outputs/parseRM.pl --land 50,1 --parse --fa -n -p -i ../VI_EUNL24_hifi.asm.bp.p_ctg.nuc.purge.fasta.align
```

### Downstream figures

The parseRM outputs (`*.landscape*.tab`) feed into the **stacked bar plot** and **repeat landscape** figures in R.


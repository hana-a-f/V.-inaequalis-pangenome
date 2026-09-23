# Functional Annotation (InterProScan + DIAMOND + BLAST2GO)

This folder documents the *functional* annotation pipeline. Each isolate's representative protein set (from `genome_annotation/`) is scanned for protein-domain signatures with **InterProScan**, aligned against **UniRef90** with **DIAMOND**, and both result files are then fed into **BLAST2GO** to assign Gene Ontology (GO) terms and produce annotation summary reports.

This is the answer to *"what do the genes we predicted actually do?"* — as opposed to the structural annotation in `genome_annotation/`, which is *"where are the genes and what are their exon boundaries?"*.

## Pipeline overview

1. Set up InterProScan (one-time, per machine)
2. Run **InterProScan** on each isolate's protein FASTA → XML + GFF3 of domain hits
3. Run **DIAMOND blastp** against UniRef90 → XML BLAST result
4. Run **BLAST2GO CLI** using both files → `.annot`, `.txt`, and `.pdf` reports

## Where to run this

Run on a **Cornell BioHPC medium-memory gen2 machine** (using all 40 threads).

BLAST2GO specifically must run on **cbsumm10**, which has the `go_db` MongoDB database mounted. Verify with:

```bash
echo "show dbs" | mongo
```

Expected output includes:

```
MongoDB server version: 3.6.12
admin    0.000GB
config   0.000GB
go_db   45.877GB
local    0.000GB
bye
```

If `go_db` isn't listed, you're on the wrong machine.

## Inputs

Per isolate, the representative protein FASTA from the structural annotation:

- `<isolate>.locus_rep.prot.fa` — one entry per locus, from `genome_annotation/named.collapsed.proteins/`

The commands below are shown for **VI_19_031** as an example. Loop them over every isolate for a full run (see loops at the end of each section).

## Tools used

| Tool | Version / source |
| --- | --- |
| InterProScan | 5.71-102.0 |
| DIAMOND | Cornell BioHPC (`/programs/diamond/diamond`) |
| UniRef90 database | Cornell BioHPC (`/home2/shared/genome_db/uniref90`) |
| BLAST2GO CLI | Cornell BioHPC (`/usr/local/blast2go/blast2go_cli.run`) |
| MongoDB `go_db` | Cornell BioHPC (`cbsumm10`) |
| GO obo file | `/shared_data/blast2go/go.obo` |

## Files in this folder

| File | Purpose |
| --- | --- |
| `README.md` | This file |
| `annotation.prop` | BLAST2GO parameter file (copied from `/shared_data/blast2go/` and edited — see Step 3) |

---

## Step 1 — Install and configure InterProScan (one-time)

```bash
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.71-102.0/interproscan-5.71-102.0-64-bit.tar.gz
tar xvfz interproscan-5.71-102.0-64-bit.tar.gz
```

Edit `interproscan.properties` to match the CPU count of the machine. For a 40-thread BioHPC medium-memory gen2 machine, change these two lines:

```properties
number.of.embedded.workers=12
maxnumber.of.embedded.workers=16
```

Run the one-time setup:

```bash
cd /workdir/$USER/interproscan-5.71-102.0
python3 setup.py -f interproscan.properties
```

## Step 2 — Run InterProScan per isolate

```bash
./interproscan.sh \
  -b ipr_out_VI_19_031 \
  -f XML,GFF3 \
  -i VI_19_031.fa \
  --goterms \
  --pathways \
  --iprlookup \
  --disable-precalc \
  -t p \
  -T ./
```

Flag reference:

| Flag | Meaning |
| --- | --- |
| `-b` | Output prefix (produces `ipr_out_<iso>.xml` and `.gff3`) |
| `-f XML,GFF3` | Output formats needed downstream |
| `-i` | Input protein FASTA |
| `--goterms` | Include GO term annotations |
| `--pathways` | Include pathway (Reactome/KEGG) annotations |
| `--iprlookup` | Look up InterPro IDs |
| `--disable-precalc` | Don't use the EBI precalculated match lookup (needed for non-model organisms) |
| `-t p` | Sequence type: protein |
| `-T ./` | Temp folder in current directory |

To loop over all isolates:

```bash
for f in /path/to/genome_annotation/named.collapsed.proteins/*.fa; do
    iso=$(basename "$f" .fa)
    ./interproscan.sh -b ipr_out_${iso} -f XML,GFF3 -i "$f" \
        --goterms --pathways --iprlookup --disable-precalc -t p -T ./
done
```

## Step 3 — DIAMOND blastp vs. UniRef90

```bash
/programs/diamond/diamond blastp \
  --db /home2/shared/genome_db/uniref90 \
  --query VI_19_031.fa \
  --outfmt 5 \
  --max-target-seqs 100 \
  --max-hsps 1 \
  --evalue 1e-10 \
  -t ./ \
  --block-size 10 \
  --index-chunks 1 \
  -o VI_19_031_blastresult.xml
```

Key flags:

| Flag | Meaning |
| --- | --- |
| `--outfmt 5` | Output XML (required by BLAST2GO) |
| `--max-target-seqs 100` | Top 100 hits per query |
| `--max-hsps 1` | One HSP per subject |
| `--evalue 1e-10` | Strict E-value cutoff |
| `--block-size 10` `--index-chunks 1` | Memory/speed tuning for a machine with plenty of RAM |

## Step 4 — BLAST2GO (on cbsumm10)

Set up the working directory:

```bash
mkdir /workdir/$USER
cd /workdir/$USER
cp /shared_data/blast2go/annotation.prop ./
cp /shared_data/blast2go/go.obo ./
```

Copy the DIAMOND XML and InterProScan XML for the isolate into `/workdir/$USER`.

### Edit `annotation.prop`

Key parameters to review before running (defaults are usually fine, but tune the top block for divergent / non-model organisms):

```properties
ImportBlastResultsAlgoParameters.numberOfHits=100
ImportBlastResultsAlgoParameters.blastMinHSPLength=33
AnnotationAlgoParameters.eValueHitFilter=1.0E-10
AnnotationAlgoParameters.hspHitCoverageCutoff=0
```

If many queries return no GO term (unmapped against the BLAST2GO mapping database), raising `numberOfHits` can help — but going too high over-annotates. 100 is a reasonable ceiling.

If you're feeding in InterProScan XML (as we are), also make sure `InterProScanImportParameters.inputFormat` in `annotation.prop` is set to match the format we produced (`XML`).

### Run BLAST2GO CLI (with both inputs)

```bash
/usr/local/blast2go/blast2go_cli.run \
  -properties annotation.prop \
  -useobo go.obo \
  -loadblast VI_19_031_blastresult.xml \
  -loadips50 ipr_out_VI_19_031.xml \
  -mapping \
  -annotation \
  -statistics all \
  -saveannot myresult \
  -saveseqtable myresult \
  -savereport myresult \
  -tempfolder ./ \
  >& annotatelogfile &
```

Outputs (per isolate):

- `myresult.annot` — tab-delimited GO annotations
- `myresult.txt` — sequence table with GO IDs and descriptions
- `myresult.pdf` — summary report with charts

## Downstream

The `.annot` files feed into any downstream GO enrichment analysis (e.g. topGO, clusterProfiler in R). The `.pdf` reports are useful sanity-checks and figure material for the manuscript.

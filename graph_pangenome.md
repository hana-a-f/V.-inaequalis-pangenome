# Pangenome Graph (Minigraph-Cactus)

This folder documents the pangenome graph build. All 18 chromosome-scaffolded, PanSeq-formatted, soft-masked assemblies are combined into a single graph using **`cactus-pangenome`** (the Minigraph-Cactus pipeline), with **VI_19_031** as the reference haplotype. Graph statistics and pangenome growth (core / accessory / private) curves are then produced with **odgi** and **panacus**.

## Tool nomenclature

The tool is called several things depending on where you look — worth pinning down for the manuscript:

- **`cactus-pangenome`** — the executable
- **Minigraph-Cactus** — the pipeline name in the Cactus documentation and papers
- Not to be confused with plain **minigraph** (Li et al., which is a different tool that only produces a structural graph without base-level alignment). Minigraph-Cactus *uses* minigraph internally as one step, then adds base-level Cactus alignment on top.

Cite: Hickey et al., *Nature Biotechnology*, 2024, "Pangenome graph construction from genome alignments with Minigraph-Cactus."

## Pipeline overview

1. Simplify contig headers
2. Write the `venturia.seqfile` (sample-to-FASTA mapping)
3. Pull the Cactus Docker image
4. Run `cactus-pangenome` — produces `.full.og`, `.gfa`, `.vcf`, and viz outputs
5. Graph stats with **odgi**
6. Pangenome growth analysis with **panacus** (hist + histgrowth)
7. Plot growth curves with `panacus-visualize.py`

## Machine requirements

- **AVX2 CPU support** is required. Some Cactus steps silently fail without it.
  Check with:
  ```bash
  grep -l avx2 /proc/cpuinfo
  ```
  If nothing prints, don't run on that machine.
- Cornell BioHPC medium/large-memory machine. The graph produced here has ~11.3M nodes and ~68M path steps — heavy on memory during the build.

## Inputs (18 assemblies)

VI_19_031 is the reference. Every other isolate is a query haplotype.

## Tools used

| Tool | Version / source |
| --- | --- |
| cactus-pangenome (Minigraph-Cactus) | `quay.io/comparative-genomics-toolkit/cactus:latest` (Docker) |
| Docker wrapper | `docker1` (Cornell BioHPC rootless docker) |
| panacus | 0.2.3 (`/programs/panacus-0.2.3/bin`) |
| panacus-visualize | Python plotting script bundled with panacus |

## Files in this folder

| File | Purpose |
| --- | --- |
| `README.md` | This file |
| `venturia.seqfile` | Sample → FASTA path mapping for Cactus |
| `isolates.txt` | Auto-generated isolate list (used by panacus) |

---

## Step 1 — Simplify contig headers (**do this before anything else**)

Cactus is fussy about complex headers. Before running, check every input FASTA:

```bash
grep "^>" *_panseq.fasta.masked | head
```

Anything with spaces, unusual characters, or long descriptions should be trimmed to the bare identifier (e.g. `VI_1771_2#1#chr01`).

## Step 2 — Write the seqfile

`venturia.seqfile` is a two-column, tab-separated (or whitespace-separated) file mapping sample IDs to FASTA paths (gzipped FASTAs are supported). No tree is specified — Minigraph-Cactus infers structure from alignments.

```
VI_19_031   /workdir/cactus/VI_19_031_panseq.fasta.masked
VI_EU302    /workdir/cactus/VI_EU302_panseq.fasta.masked
VI_EU104    /workdir/cactus/VI_EU104_panseq.fasta.masked
VI_EU160    /workdir/cactus/VI_EU160_panseq.fasta.masked
VI_EU413    /workdir/cactus/VI_EU413_panseq.fasta.masked
VI_EUNL24   /workdir/cactus/VI_EUNL24_panseq.fasta.masked
VI_EUNL19   /workdir/cactus/VI_EUNL19_panseq.fasta.masked
VI_EU301    /workdir/cactus/VI_EU301_panseq.fasta.masked
VI_EUB04    /workdir/cactus/VI_EUB04_panseq.fasta.masked
VI_18_030   /workdir/cactus/VI_18_030_panseq.fasta.masked
VI_19_004   /workdir/cactus/VI_19_004_panseq.fasta.masked
VI_19_011   /workdir/cactus/VI_19_011_panseq.fasta.masked
VI_18_019   /workdir/cactus/VI_18_019_panseq.fasta.masked
VI_18_033   /workdir/cactus/VI_18_033_panseq.fasta.masked
VI_18_037   /workdir/cactus/VI_18_037_panseq.fasta.masked
VI_1771_2   /workdir/cactus/VI_1771_2_panseq.fasta.masked
VI_1797_2   /workdir/cactus/VI_1797_2_panseq.fasta.masked
VI_1797_9   /workdir/cactus/VI_1797_9_panseq.fasta.masked
```

## Step 3 — Pull the Cactus Docker image

```bash
docker1 pull quay.io/comparative-genomics-toolkit/cactus:latest
```

## Step 4 — Run `cactus-pangenome`

Basic command syntax:

```
cactus-pangenome <jobStorePath> <seqFile> --outDir <output directory> --outName <output file prefix> --reference <reference sample name>
```

Where:
- `<jobStorePath>` — intermediate-file directory, must be accessible to all worker processes
- `<seqFile>` — the sample-to-FASTA mapping above
- `<outDir>` — output directory
- `<outName>` — output file prefix
- `<reference>` — the reference sample name (`VI_19_031`)

Actual command run:

```bash
docker1 run --rm -it \
  -v /workdir/hf332/cactus:/workdir/cactus \
  quay.io/comparative-genomics-toolkit/cactus:latest \
  cactus-pangenome \
  /workdir/cactus/jobStore \
  /workdir/cactus/venturia.seqfile \
  --outName venturia-pangenome \
  --outDir /workdir/cactus/venturia-output \
  --reference VI_19_031 \
  --viz \
  --odgi \
  --gfa \
  --vcf \
  --mgCores 8 \
  --mapCores 4 \
  --consCores 8 \
  --indexCores 4 \
  --logFile /workdir/cactus/venturia-pangenome.log
```


## Step 5 — Pangenome growth analysis (panacus)

Panacus estimates core / accessory / private-genome sizes by treating each isolate as a "sample" and counting how much sequence is shared among varying numbers of them.

Get the list of isolates from the graph paths:

```bash
export PATH=/programs/panacus-0.2.3/bin:$PATH

odgi paths -i venturia-pangenome.full.og -L > all.paths.txt
cut -d'#' -f1 all.paths.txt | sort -u > isolates.txt
```

### Histogram

```bash
RUST_LOG=info panacus hist \
  -t 16 \
  -c bp \
  -s isolates.txt \
  -S \
  venturia-pangenome.gfa \
  > v_inaequalis.hist.tsv
```

### Growth curves — main quorum settings

`-l 1,18,17,2` and `-q 0,1,1,0` give the classic core (present in all), accessory, and private-genome curves:

```bash
RUST_LOG=info panacus histgrowth \
  -t 8 \
  -l 1,18,17,2 \
  -q 0,1,1,0 \
  -c bp \
  -S -a \
  venturia-pangenome.gfa \
  > venturia.histgrowth.bp.tsv
```

### Growth curves — alternative quorums (example run)

`-q 0,1,0.5,0.1` explores relaxed quorum thresholds (50% and 10%):

```bash
RUST_LOG=info panacus histgrowth \
  -t 8 \
  -q 0,1,0.5,0.1 \
  -c bp \
  -S -a \
  venturia-pangenome.gfa \
  > venturia.histgrowth.bp.example.tsv
```

## Step 6 — Plot growth curves

```bash
python panacus-visualize.py venturia.histgrowth.bp.tsv > v_inaequalis.histgrowth.bp.pdf
python panacus-visualize.py venturia.histgrowth.bp.example.tsv > venturia.histgrowth.bp.example.pdf
```

`panacus-visualize.py` ships with panacus.

## Outputs

Key files under `venturia-output/`:

| File | Description |
| --- | --- |
| `venturia-pangenome.full.og` | odgi-format graph — input to odgi stats and panacus |
| `venturia-pangenome.gfa` | GFA graph — universal format, input to panacus |
| `venturia-pangenome.vcf.gz` | VCF of variants relative to VI_19_031 |
| `venturia-pangenome.viz.*` | 1D visualisations |
| `v_inaequalis.hist.tsv` | Panacus histogram of shared sequence content |
| `venturia.histgrowth.bp.tsv` | Growth curve, main quorum settings |
| `venturia.histgrowth.bp.example.tsv` | Growth curve, alternative quorums |
| `v_inaequalis.histgrowth.bp.pdf` | Growth curve figure |
| `venturia.histgrowth.bp.example.pdf` | Alternative growth curve figure |

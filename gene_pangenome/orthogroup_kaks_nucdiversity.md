# Orthogroup Analysis (Ka/Ks, π)

This folder documents per-orthogroup evolutionary analyses starting from OrthoFinder output. For every orthogroup:

- **Codon-level alignment** with **MAFFT**
- **Ka/Ks** (nonsynonymous vs. synonymous substitution ratio) for every pair of genes with **KaKs_Calculator** — summarised as the median per orthogroup
- **Nucleotide diversity (π)** with **EggLib** 


## Pipeline overview

1. **Extract CDS per orthogroup** from OrthoFinder's `Orthogroups.tsv` + your combined `all.cds.fa`
2. **Align** each orthogroup's CDS with MAFFT (`--auto`)
3. **Convert** each alignment to pairwise AXT format for KaKs_Calculator
4. **Compute Ka/Ks** for every gene pair (Nei-Gojobori method, `-m NG`)
5. **Aggregate** to median Ka/Ks per orthogroup
6. **Compute π** per orthogroup with EggLib

## Inputs

- `Orthogroups.tsv` — from OrthoFinder 
- `all.cds.fa` — combined CDS FASTA across all isolates

## Tools used

| Tool | Version / source |
| --- | --- |
| MAFFT | 7.520 (`/programs/mafft/bin`) — Katoh & Standley, 2013 |
| KaKs_Calculator | 2.0 — https://github.com/kullrich/kakscalculator2 — Wang et al., 2010 |
| EggLib | 3.6.0 — Siol et al., 2022 |
| topGO | Alexa & Rahnenführer, 2024 (R Bioconductor) |
| Biopython | For CDS extraction and AXT conversion |

## Files in this folder

| File | Purpose |
| --- | --- |
| `README.md` | This file |
| `extract_orthogroup_cds.py` | Read `Orthogroups.tsv`, write one CDS FASTA per orthogroup |
| `mafft_align.sh` | Loop MAFFT over every orthogroup FASTA |
| `mafft_to_axt.py` | Convert MAFFT alignments to per-pair AXT for KaKs_Calculator |
| `run_kaks.sh` | Loop KaKs_Calculator over every AXT file |
| `compute_pi_egglib.py` | Compute per-orthogroup π with EggLib 

---

## Step 1 — Extract CDS per orthogroup

Reads OrthoFinder's `Orthogroups.tsv` (rows = orthogroups, columns = species, cells = comma-separated gene IDs) and writes one FASTA per orthogroup containing all its CDS sequences:

```python
import sys
from Bio import SeqIO

orthogroups = "Orthogroups.tsv"
cds_fasta = "all.cds.fa"

seqs = SeqIO.to_dict(SeqIO.parse(cds_fasta, "fasta"))

with open(orthogroups) as f:
    for line in f:
        parts = line.strip().split("\t")
        og = parts[0]
        genes = []

        for col in parts[1:]:
            if col != "":
                genes.extend([g.strip() for g in col.split(",")])

        with open(f"orthogroups/{og}.fa", "w") as out:
            for g in genes:
                if g in seqs:
                    SeqIO.write(seqs[g], out, "fasta")
```

Output: `orthogroups/OG*.fa` (one per orthogroup).

## Step 2 — MAFFT alignment

```bash
export PATH=/programs/mafft/bin:$PATH
mkdir alignments

for f in orthogroups/*.fa; do
    base=$(basename $f .fa)
    mafft --auto $f > alignments/${base}.aln.fa
done
```

Output: `alignments/OG*.aln.fa` — aligned CDS per orthogroup.

## Step 3 — Convert MAFFT alignments to AXT

KaKs_Calculator wants **pairwise** AXT format (a name line + two aligned sequences per pair). This script explodes each orthogroup's multi-sequence alignment into every gene pair:

```python
from Bio import SeqIO
import itertools
import os

indir = "alignments"
outdir = "axt"
os.makedirs(outdir, exist_ok=True)

for file in os.listdir(indir):
    if not file.endswith(".aln.fa"):
        continue

    og = file.replace(".aln.fa","")
    seqs = list(SeqIO.parse(f"{indir}/{file}", "fasta"))

    with open(f"{outdir}/{og}.axt","w") as out:
        for a, b in itertools.combinations(seqs, 2):
            out.write(f"{a.id}_{b.id}\n")
            out.write(str(a.seq) + "\n")
            out.write(str(b.seq) + "\n\n")
```

Output: `axt/OG*.axt` — one AXT per orthogroup with all pairs stacked inside.

## Step 4 — Run KaKs_Calculator

Install and build once:

```bash
git clone https://github.com/kullrich/kakscalculator2.git
cd kakscalculator2
export PREFIX="."
make clean && make

export PATH=$PATH:/workdir/hf332/dnds_kaks/kakscalculator2/bin
```

Run over every orthogroup's AXT file, using the **Nei-Gojobori** method (`-m NG`):

```bash
mkdir kaks_results

for f in axt/*.axt; do
    base=$(basename $f .axt)
    KaKs_Calculator \
        -i $f \
        -o kaks_results/${base}.kaks \
        -m NG
done
```

Output: `kaks_results/OG*.kaks` — one row per gene pair with Ka, Ks, and Ka/Ks columns.


## Step 5 — Compute π per orthogroup (EggLib)

```python
import os
import egglib

indir = "alignments"

cs = egglib.stats.ComputeStats()
cs.add_stats("Pi")
cs.add_stats("lseff")

print("Orthogroup\tPi")

for file in os.listdir(indir):
    if not file.endswith(".aln.fa"):
        continue

    og = file.replace(".aln.fa", "")

    aln = egglib.io.from_fasta(
        os.path.join(indir, file),
        alphabet=egglib.alphabets.DNA
    )

    stats = cs.process_align(aln)
    lseff = stats["lseff"]
    if lseff == 0:
        continue

    pi = stats["Pi"] / lseff
    print(f"{og}\t{pi}")
```

Output: TSV with `Orthogroup` and `Pi` columns.

normalise differently — EggLib uses effective aligned length (`lseff`) which excludes gappy positions, while the manual version normalises by full alignment length. That's why EggLib is preferred: it handles missing/gapped positions correctly.

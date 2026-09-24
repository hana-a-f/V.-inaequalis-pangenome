# Phylogenomics (BUSCO-based species tree)

This folder documents the species-level phylogenetic analysis. **BUSCO** is run on each isolate's proteome to identify single-copy orthologs; **BUSCO_phylogenomics** collects the shared single-copies across taxa into concatenated alignments; and **IQ-TREE** builds the final species tree from the resulting supermatrix with model selection and bootstrap support.

 **Venturia nashicola** (`v_nash`) 

## Pipeline overview

1. Build the BUSCO Apptainer image (once)
2. Run BUSCO on each proteome (outgroups + all in-group *V. inaequalis* isolates)
3. Run **BUSCO_phylogenomics.py** — aligns shared single-copy BUSCOs, trims, builds gene trees, concatenates into a supermatrix
4. Run **IQ-TREE** on the supermatrix with ModelFinder Plus and SH-aLRT + ultrafast bootstrap support

## Tools used

| Tool | Version / source |
| --- | --- |
| BUSCO | 6.1.0 (`docker://quay.io/biocontainers/busco:6.1.0--pyhdfd78af_2`, Apptainer) |
| BUSCO lineage | `ascomycota_odb10` |
| BUSCO_phylogenomics | `/programs/BUSCO_phylogenomics/BUSCO_phylogenomics.py` (conda env `BUSCO_phylogenomics`) |
| trimAl | Via BUSCO_phylogenomics (strategy `automated1`) |
| IQ-TREE | 2.2.2.6 (`/programs/iqtree-2.2.2.6-Linux/bin/iqtree2`) |

## Files in this folder

| File | Purpose |
| --- | --- |
| `README.md` | This file |
| `busco.sif` | Built BUSCO Apptainer image (or rebuild — see Step 1) |

---

## Step 1 — Build the BUSCO container

Replace `6.1.0--pyhdfd78af_2` with the latest tag if needed:

```bash
apptainer build busco.sif docker://quay.io/biocontainers/busco:6.1.0--pyhdfd78af_2
```

Sanity checks:

```bash
# Print software menu
apptainer exec --bind $PWD --pwd $PWD busco.sif busco -h

# List available lineage datasets
apptainer exec --bind $PWD --pwd $PWD busco.sif busco --list-datasets
```

Pick the correct lineage for the taxonomic level of interest — here **`ascomycota_odb10`** for fungi.

## Step 2 — Run BUSCO per proteome

Parameters:
- `-i` — input protein FASTA
- `-m protein` — protein mode (proteomes; use `-m geno` for raw genomes)
- `-l ascomycota_odb10` — lineage dataset
- `-o` — output prefix
- `-c` — CPU cores


Outgroup 1 — *Venturia nashicola* (run in parallel with `&`):

```bash
apptainer exec --bind $PWD --pwd $PWD busco.sif \
    busco -i v_nash.fa -m protein -l ascomycota_odb10 -o V.nash -c 12 &

wait
echo "Both BUSCO runs complete"
```

Repeat this same pattern for each of the 18 *V. inaequalis* isolate proteomes (from `genome_annotation/`). The BUSCO output for each isolate goes in its own directory inside the working folder.

## Step 3 — BUSCO_phylogenomics: align, trim, gene trees

Activate the environment and run:

```bash
source /programs/miniconda3/bin/activate BUSCO_phylogenomics

/programs/BUSCO_phylogenomics/BUSCO_phylogenomics.py \
    -i /workdir/hf332/busco_phy/busco \
    -o BUSCO_phylogeny_90 \
    -t 12 \
    -psc 90 \
    --trimal_strategy automated1 \
    --gene_tree_program iqtree
```

Parameters:
- `-i` — parent directory containing one BUSCO output subfolder per taxon
- `-o` — output directory
- `-t` — threads
- `-psc 90` — **percent single-copy threshold**: only include BUSCOs present as complete single-copy in ≥90% of taxa
- `--trimal_strategy automated1` — trimAl's automated1 heuristic for alignment trimming
- `--gene_tree_program iqtree` — use IQ-TREE for per-gene tree inference

Output includes per-gene alignments, trimmed alignments, gene trees, and a concatenated **`SUPERMATRIX.phylip`** (or similar).

## Step 4 — Build the final species tree with IQ-TREE

ModelFinder Plus (`-m MFP`) picks the best-fit substitution model automatically; 1000 SH-aLRT replicates + 1000 ultrafast bootstrap replicates give branch support:

```bash
/programs/iqtree-2.2.2.6-Linux/bin/iqtree2 \
    -s SUPERMATRIX.phylip \
    -m MFP \
    -alrt 1000 \
    -B 1000 \
    -T AUTO
```

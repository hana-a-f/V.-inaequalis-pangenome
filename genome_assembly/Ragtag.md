# RagTag Scaffolding

This folder documents the reference-guided scaffolding of the 16 contig-level assemblies against the chromosome-level VI_19_031 reference. The output is a set of pseudochromosome-level assemblies for every isolate, with headers in the PanSeq format.

## Inputs

Every input assembly is already:
- **Purged** (purge_haplotigs for 16 samples)
- **Mitochondrial-DNA-removed**
- **Repeat-masked**
- Named consistently: `{isolate}.nuc.purge.masked`

Reference: `VI_19_031.chr.ref.nuc.fasta.masked` (from `hic_assembly/`).


## Tools used

| Tool | Version |
| --- | --- |
| RagTag | (add version) |
| MUMmer / nucmer | (system) |
| minimap2 | (system) |
| seqkit | 0.15.0 |
| samtools | (system) |
| GNU parallel | (system) |
| BUSCO | 5.5.0 (`ascomycota_odb10`) |
| Merqury | 1.3 (Singularity) |

## Files in this folder

| File | Purpose |
| --- | --- |
| `README.md` | This file |
| `ragtag.cmds.txt` | Auto-generated RagTag scaffold command list |

---

## Step 1 — Preprocess: standardise names + fix headers

Activate environment and set up seqkit:

```bash
export PATH=/programs/seqkit-0.15.0:$PATH
source /programs/miniconda3/bin/activate ragtag
```

Rename every masked assembly to a consistent `{isolate}.nuc.purge.masked` filename (ran before to make sure works before adding in echo):

```bash
for f in *.fasta.masked *.fa.masked; do
    [[ "$f" == "VI_19_031.chr.ref.nuc.fasta.masked" ]] && continue
    iso=$(echo "$f" | sed -E 's/(_hifi\.asm\.bp\.p_ctg\.nuc\.purge|_curated|\.2\.0)?(\.fasta|\.fa)?\.masked//')
    new="${iso}.nuc.purge.masked"
    mv -n "$f" "$new" && echo "renamed: $f → $new"
done
```

Fix the NCBI-style `lcl|...` headers on the four pre-screened samples (VI_1797_2, VI_18_019, VI_18_033, VI_19_004 — same four flagged for contamination in `genome_assembly/`):

```bash
sed -i 's/^>lcl|/>/'  VI_1797_2.nuc.purge.masked
sed -i 's/ .*//'      VI_1797_2.nuc.purge.masked

sed -i 's/^>lcl|/>/'  VI_18_019.nuc.purge.masked
sed -i 's/ .*//'      VI_18_019.nuc.purge.masked

sed -i 's/^>lcl|/>/'  VI_18_033.nuc.purge.masked
sed -i 's/ .*//'      VI_18_033.nuc.purge.masked

sed -i 's/^>lcl|/>/'  VI_19_004.nuc.purge.masked
sed -i 's/ .*//'      VI_19_004.nuc.purge.masked
```

## Step 2 — Scaffold against the reference (nucmer)

Auto-generate a RagTag command per isolate, then run 4 in parallel:

```bash
#!/bin/bash
set -euo pipefail

REF="VI_19_031.chr.ref.nuc.fasta.masked"
CMD_FILE="ragtag.cmds.txt"

> "$CMD_FILE"   # start fresh

for f in *.masked; do
    [[ "$f" == "$REF" ]] && continue
    isolate=$(basename "$f" .nuc.purge.masked)
    echo "ragtag.py scaffold \
  --aligner nucmer \
  -f 5000 \
  -u \
  -o ${isolate}.ragtag \
  ${REF} \
  ${f}" >> "$CMD_FILE"
done

echo "Commands written to $CMD_FILE"
```

Run in parallel:

```bash
source /programs/miniconda3/bin/activate ragtag
parallel -j 4 < ragtag.cmds.txt
```

Prefix every output file with its isolate name and combine for easier review:

```bash
#!/bin/bash
set -euo pipefail

mkdir -p combined_prefixed

for dir in *.ragtag; do
    [[ ! -d "$dir" ]] && continue
    isolate=$(basename "$dir" .ragtag)
    echo "Processing $isolate ..."
    for f in "$dir"/*; do
        [[ -f "$f" ]] || continue
        base=$(basename "$f")
        cp "$f" "combined_prefixed/${isolate}_${base}"
    done
done

echo "All prefixed RagTag files copied to combined_prefixed/"
```

Pull the `.delta` alignment files aside for inspection:

```bash
mkdir delta
cp *.delta delta
```

## Step 3 — Special case: VI_18_037 (misassembly)


```bash
ragtag.py scaffold \
  --aligner nucmer \
  -f 5000 -u \
  -o VI_18_037.2.ragtag \
  VI_19_031.chr.ref.nuc.fasta.masked \
  VI_18_037.nuc.purge.2.masked

for f in *; do mv "$f" "VI_18_037_${f}"; done
```

## Step 4 — Special case: VI_EUNL24 (misassembly)

Run `ragtag correct` first (with the isolate's own reads) to identify misassembly breakpoints:

```bash
grep '>' VI_EUNL24.nuc.purge.masked > skip.VI_EUNL24.txt

ragtag.py correct \
  --aligner nucmer \
  -f 5000 \
  -j skip.VI_EUNL24.txt \
  -o VI_EUNL24.ragtag.correct \
  -u \
  -R VI_EUNL24.fastq.gz \
  -T corr \
  VI_19_031.chr.ref.nuc.fasta.masked \
  VI_EUNL24.nuc.purge.masked
```

Inspect the alignment to find the split point on `ptg000001l`:

```bash
show-coords -rclT ragtag.scaffold.asm.delta | grep "ptg000001l" > ptg000001l.coords
sort -k3,3n ptg000001l.coords > ptg000001l.coords.sorted
# Split point identified at 5,014,854
```

Manually split the contig, rename the two halves, and rebuild the assembly:

```bash
seqkit grep -n -p ptg000001l VI_EUNL24.nuc.purge.masked > ptg000001l.fa
seqkit subseq -r 1:5014853         ptg000001l.fa > ptg000001l_part1.fa
seqkit subseq -r 5014854:8998627   ptg000001l.fa > ptg000001l_part2.fa
sed -i 's/>ptg000001l/>ptg000001l_a/' ptg000001l_part1.fa
sed -i 's/>ptg000001l/>ptg000001l_b/' ptg000001l_part2.fa

seqkit grep -v -n -p ptg000001l VI_EUNL24.nuc.purge.masked > VI_EUNL24_wo_ptg000001l.fa
cat VI_EUNL24_wo_ptg000001l.fa ptg000001l_part1.fa ptg000001l_part2.fa \
  > VI_EUNL24.nuc.purge.2.masked
```

Re-scaffold the corrected assembly:

```bash
ragtag.py scaffold \
  --aligner nucmer \
  -f 5000 -u \
  -o VI_EUNL24.2.ragtag \
  VI_19_031.chr.ref.nuc.fasta.masked \
  VI_EUNL24.nuc.purge.2.masked

for f in *; do mv "$f" "VI_EUNL24_${f}"; done
```

## Step 5 — Special case: VI_19_011 (misassembly)

Same approach — `ragtag correct` first, identify split point on `ptg000002l`, manually split and re-scaffold:

```bash
grep '>' VI_19_011.nuc.purge.masked > skip.VI_19_011.txt

ragtag.py correct \
  --aligner nucmer \
  -f 5000 \
  -j skip.VI_19_011.txt \
  -o VI_19_011.ragtag.correct \
  -u \
  -R VI_19_011.fastq.gz \
  -T corr \
  VI_19_031.chr.ref.nuc.fasta.masked \
  VI_19_011.nuc.purge.masked

show-coords -rclT ragtag.scaffold.asm.delta | grep "ptg000002l" > ptg000002l.coords
sort -k3,3n ptg000002l.coords > ptg000002l.coords.sorted

seqkit grep -n -p ptg000002l VI_19_011.nuc.purge.masked > ptg000002l.fa
seqkit subseq -r 1:462668           ptg000002l.fa > ptg000002l_part1.fa
seqkit subseq -r 2282836:4496502    ptg000002l.fa > ptg000002l_part2.fa
sed -i 's/>ptg000002l/>ptg000002l_part1/' ptg000002l_part1.fa
sed -i 's/>ptg000002l/>ptg000002l_part2/' ptg000002l_part2.fa

seqkit grep -v -n -p ptg000002l VI_19_011.nuc.purge.masked > VI_19_011_wo_ptg000002l.fa
cat VI_19_011_wo_ptg000002l.fa ptg000002l_part1.fa ptg000002l_part2.fa \
  > VI_19_011.nuc.purge.2.masked

ragtag.py scaffold \
  --aligner nucmer \
  -f 5000 -u \
  -o VI_19_011.2.ragtag \
  VI_19_031.chr.ref.nuc.fasta.masked \
  VI_19_011.nuc.purge.2.masked

# Also produced a run using the un-corrected input, for comparison:
ragtag.py scaffold \
  --aligner nucmer \
  -f 5000 -u \
  -o VI_19_011.1_22_26.ragtag \
  VI_19_031.chr.ref.nuc.fasta.masked \
  VI_19_011.nuc.purge.masked

for f in *; do mv "$f" "VI_19_011_${f}"; done
```

## Step 6 — Rename scaffolded FASTAs to PanSeq format

Only run this **after** checking the `.delta` files. Replaces the VI_19_031 prefix inherited from the reference, strips `_RagTag`, and adds the PanSeq isolate tag:

```bash
#!/bin/bash

REF_PREFIX="VI_19_031#"

for f in *_ragtag.scaffold.fasta; do
    isolate=$(echo "$f" | sed -E 's/_ragtag\.scaffold\.fasta//')
    echo "Processing $isolate ..."

    # 1. Replace reference prefix with isolate prefix
    sed "s/^>${REF_PREFIX}/>${isolate}#/" "$f" > "${isolate}_renamed.fasta"

    # 2. Remove RagTag suffix in headers
    sed -i 's/_RagTag//g' "${isolate}_renamed.fasta"

    # 3. Add isolate tag before contig IDs for pangenome naming
    sed "s/^>ptg/>${isolate}#1#ptg/" "${isolate}_renamed.fasta" > "${isolate}_panseq.fasta"

    rm -f "${isolate}_renamed.fasta"
done

echo "All renamed and pangenome-formatted FASTAs created."
```

## Step 7 — QC

### Chromosome sizes

```bash
for f in *_panseq.fasta; do
    samtools faidx "$f"
done

for f in *_panseq.fasta; do
    iso=$(basename "$f" _panseq.fasta)
    awk -v iso="$iso" '{print iso"\t"$1"\t"$2}' "${f}.fai"
done > all_chrom_sizes.tsv
```

### seqkit stats (before and after scaffolding)

```bash
export PATH=/programs/seqkit-0.15.0:$PATH
seqkit stats -a *_panseq.fasta > all_seqstats.txt
seqkit stats -a *.masked       > all_seqstats.before.ragtag.txt
```


### BUSCO

```bash
source /programs/miniconda3/bin/activate busco-5.5.0

touch buscos.cmds

for f in *_panseq.fasta; do
    iso=$(basename "$f" _panseq.fasta)
    echo "busco -f -i $f -l ascomycota_odb10 -o ${iso}.busco.g -m geno -c 8" >> buscos.cmds
done

parallel -j 4 < buscos.cmds

# Spot-check on a single scaffolded assembly:
busco -f -i ragtag.scaffold.fasta -l ascomycota_odb10 -o check -m geno -c 8
```

### Merqury (k-mer completeness + QV)

Build meryl databases in parallel:

```bash
touch meryl.cmds
export MERQURY=/merqury

for f in /workdir/hf332/merqury/*.fastq; do
    base=$(basename "$f" .fastq)
    echo "singularity run \
      -B /workdir/hf332/merqury:/data \
      --pwd /data \
      /programs/merqury-1.3/merqury.sif \
      meryl count k=21 threads=16 output ${base}.meryl /data/${base}.fastq" >> meryl.cmds
done

parallel -j 2 < meryl.cmds
```

Run merqury.sh per isolate:

```bash
touch merqury.cmds
for asm in *_panseq.fasta; do
    prefix=${asm%_panseq.fasta}
    mkdir -p merqury_out/$prefix

    echo "cd merqury_out/$prefix && \
singularity run \
  --env MERQURY=/merqury \
  --bind /workdir/$USER:/data \
  --pwd /data/merqury_out/$prefix \
  /programs/merqury-1.3/merqury.sif \
  merqury.sh /data/${prefix}.meryl /data/${prefix}_panseq.fasta ${prefix}.merqury" \
>> merqury.cmds
done

parallel -j 2 < merqury.cmds
```

Summarise QV and completeness across all isolates:

```bash
for f in *merqury.qv; do
    printf "%s\t" "$f"
    cat "$f"
done > merqury_qv_summary.tsv

for f in *completeness.stats; do
    printf "%s\t" "$f"
    cat "$f"
done > merqury_kmer_summary.tsv
```

## Next step

Each isolate now has a chromosome-scaffolded assembly in PanSeq format (`{isolate}_panseq.fasta`). These feed into the pangenome graph build.

Repeat masking was re-run on the PanSeq FASTAs (to keep GFF3 coordinates consistent with the graph-input FASTAs) — see `repeat_masking/`.

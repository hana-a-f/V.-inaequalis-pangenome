# Effector Family Assignment

This folder documents how *V. inaequalis* secreted proteins (from `effector_prediction/`) are assigned to **expanded effector families defined by Rocafort et al. (2023)**. Each family is turned into an HMM profile and used to search the *V. inaequalis* secretome; BLASTP is run in parallel against the same family sequences. Hits from both searches are then combined to give each secreted protein its family assignment, which is finally mapped back to orthogroups.


## Pipeline overview

1. **Combine secretomes** into one BLAST database across all 18 isolates
2. **Assemble family sequences** from Rocafort et al. (2023) — one FASTA per family, only families with >3 members
3. **Signal peptide processing** on family sequences (SignalP → trim to mature peptide)
4. **Align** each family with MAFFT (`--auto`)
5. **Trim** alignments with trimAl (`-gt 0.3 -cons 50`)
6. **Build HMMs** and **search** with HMMER (`hmmbuild` + `hmmsearch --domtblout`)
7. **BLASTP** each family against the combined secretome (E-value 1e-5, ≥50% query coverage)
8. **AvrRvi6 special case** — separate BLASTP against the secretome
9. **Filter + combine** HMMER (i-evalue < 1e-5) and BLASTP hits → per-protein family assignment
10. **Map assignments to orthogroups**

## Inputs

- `<isolate>.locus.filt.secretome.fa` (18 files) — from `effector_prediction/`
- `families/*.fa` — per-family protein sequences from Rocafort et al. (2023) supplementary data (only families with >3 members)
- `avr6.fa` — AvrRvi6 protein sequence (from ref genome or published)
- Known MAX effectors (for family assembly reference):
  - PDB 2MYV (https://www.rcsb.org/structure/2MYV)
  - UniProt L7J9X3 (https://www.uniprot.org/uniprotkb/L7J9X3/entry)
  - PDB 6R5J (https://www.rcsb.org/structure/6R5J)
  - PDB 2LW6 (https://www.rcsb.org/structure/2LW6)

## Tools used

| Tool | Version / source |
| --- | --- |
| makeblastdb / BLASTP | 2.16.0 (system) — Altschul et al., 1997 |
| SignalP | 5.0 (`/programs/signalp-5.0/bin/signalp`) |
| MAFFT | 7.520 (`/programs/mafft/bin`) — Katoh & Standley, 2013 |
| trimAl | 1.4 (`/programs/trimal-1.4/source`) — Capella-Gutiérrez et al., 2009 |
| HMMER | 3.4 (`hmmbuild`, `hmmsearch`) — Finn et al., 2011 |
| GNU parallel | (system) |

## Files in this folder

| File | Purpose |
| --- | --- |
| `families/` | Per-family FASTA made from (from Rocafort et al. 2023) |
| `avr6.fa` | AvrRvi6 protein query for Step 8 (Sannier et al. 2025)|
| `Venturia_all_secretome.fa` | Concatenated secretome from all 18 isolates |

---

## Step 1 — Combine secretomes into one BLAST database

```bash
cat *.locus.filt.secretome.fa > Venturia_all_secretome.fa
makeblastdb -in Venturia_all_secretome.fa -dbtype prot -out Venturia_secretome
```

## Step 2 — Assemble family sequences from Rocafort et al. (2023)

Download the family sequences from Rocafort et al. (2023) — see paper: https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1011294 — keeping only families with **more than three members**. One FASTA per family in `families/`.

## Step 3a — SignalP on family sequences

```bash
mkdir -p signalp
: > signalp.cmds

for f in families/*.fa; do
    fam=$(basename "$f" .fa)
    echo "/programs/signalp-5.0/bin/signalp \
      -org euk \
      -batch 20000 \
      -fasta $f \
      -gff3 \
      -prefix signalp/${fam}" >> signalp.cmds
done

parallel -j 5 < signalp.cmds
```

## Step 3b — Trim the signal peptide to give mature-peptide FASTAs

Reads each SignalP GFF3, identifies the `signal_peptide` feature end coordinate per protein, and strips those N-terminal residues from the FASTA:

```bash
mkdir -p families_mature

for fa in families/*.fa; do
    fam=$(basename "$fa" .fa)
    gff="signalp/${fam}.gff3"
    out="families_mature/${fam}.mature.fa"

    echo "Trimming signal peptides for $fam"

    awk '
    BEGIN { FS="\t" }

    # First pass: read SignalP GFF
    FNR==NR && $3=="signal_peptide" {
        sp_end[$1] = $5
        next
    }

    # Second pass: read FASTA, strip signal peptide
    /^>/ {
        if (seq_id != "") {
            end = (seq_id in sp_end ? sp_end[seq_id] : 0)
            trimmed = substr(seq, end + 1)
            if (length(trimmed) > 0)
                print ">" seq_id "\n" trimmed
        }
        seq_id = substr($0, 2)
        seq = ""
        next
    }
    { seq = seq $0 }

    END {
        if (seq_id != "") {
            end = (seq_id in sp_end ? sp_end[seq_id] : 0)
            trimmed = substr(seq, end + 1)
            if (length(trimmed) > 0)
                print ">" seq_id "\n" trimmed
        }
    }
    ' "$gff" "$fa" > "$out"
done

echo "=== SIGNAL PEPTIDE TRIMMING COMPLETE ==="
```

**Sanity check** — original vs. trimmed sequence counts per family:

```bash
for f in families/*.fa; do
    fam=$(basename "$f" .fa)
    echo -n "$fam  "
    grep -c "^>" "$f"
    grep -c "^>" "families_mature/${fam}.mature.fa"
done
```

## Step 4 — Align each family with MAFFT

```bash
mkdir -p msa
export PATH=/programs/mafft/bin:$PATH

for f in families_mature/*.fa; do
    fam=$(basename "$f" .mature.fa)
    mafft --auto "$f" > "msa/${fam}.aln.fa"
done
```

## Step 5 — Trim alignments with trimAl

`-gt 0.3` = drop columns present in <30% of sequences; `-cons 50` = keep at least 50% of columns as a conservation threshold:

```bash
export PATH=/programs/trimal-1.4/source:$PATH
mkdir -p msa_trimmed

for aln in msa/*.aln.fa; do
    fam=$(basename "$aln" .aln.fa)
    echo "Trimming alignment for $fam"
    trimal \
      -in "$aln" \
      -out "msa_trimmed/${fam}.aln.trim.fa" \
      -gt 0.3 \
      -cons 50
done

echo "=== ALIGNMENT TRIMMING COMPLETE ==="
```

## Step 6 — Build HMMs and search the secretome

```bash
mkdir -p hmmer_build hmmer_results

for aln in msa_trimmed/*.aln.trim.fa; do
    fam=$(basename "$aln" .aln.trim.fa)

    echo "Building HMM for $fam"
    hmmbuild "hmmer_build/${fam}.hmm" "$aln"

    echo "Searching secretome with $fam"
    hmmsearch \
      --domtblout "hmmer_results/${fam}.domtbl" \
      -o "hmmer_results/${fam}.out" \
      --cpu 8 \
      "hmmer_build/${fam}.hmm" \
      Venturia_all_secretome.fa
done

echo "=== ALL HMMER RUNS COMPLETE ==="
```

## Step 7 — BLASTP each family against the secretome

```bash
mkdir -p blastp

for f in families_mature/*.fa; do
    fam=$(basename "$f" .mature.fa)
    blastp \
      -query "$f" \
      -db Venturia_secretome \
      -evalue 1e-5 \
      -qcov_hsp_perc 50 \
      -seg no \
      -max_target_seqs 10 \
      -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
      -num_threads 8 \
      -out "blastp/${fam}.tsv"
done
```

## Step 8 — AvrRvi6 special case (targeted BLASTP)

```bash
mkdir -p avr6
blastp \
  -query avr6.fa \
  -db Venturia_secretome \
  -evalue 1e-5 \
  -qcov_hsp_perc 50 \
  -seg no \
  -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
  -out avr6/Avr6_vs_secretome.blastp.tsv
```

## Step 9 — Combine BLASTP + HMMER hits (add script)

*done in R

## Step 10 — Map family assignments to orthogroups (add script)

Cross-references your `family_assignments.tsv` (from Step 9) with `orthology/Orthogroups.tsv` (from OrthoFinder). Each family-assigned protein has an orthogroup ID; roll up to per-orthogroup family assignments so a downstream analysis can say "orthogroup OG0001234 is family X."

Suggested output: `family_by_orthogroup.tsv` with columns `Orthogroup`, `Family`, `N_proteins_in_family`.

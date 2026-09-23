#!/usr/bin/env python3

from Bio import SeqIO
import sys
import os

# -----------------------------
# arguments
# -----------------------------
# python collapse_longest_isoform_by_locus.py \
#   VI_18_019.merged.gff3 \
#   VI_18_019.longest.prot.fa \
#   VI_18_019

gff_file = sys.argv[1]
prot_fa = sys.argv[2]
prefix = sys.argv[3]

out_fa = f"{prefix}.locus_rep.prot.fa"
out_bed = f"{prefix}.locus_rep.bed"

# -----------------------------
# load proteins
# -----------------------------
proteins = SeqIO.to_dict(SeqIO.parse(prot_fa, "fasta"))
protein_ids = set(proteins.keys())

print(f"[INFO] Loaded {len(protein_ids)} longest-isoform proteins")

# -----------------------------
# parse GFF: locus → transcripts + coords
# -----------------------------
locus_data = []
tx_to_locus = {}

with open(gff_file) as gff:
    for line in gff:
        if line.startswith("#"):
            continue

        parts = line.rstrip().split()
        if len(parts) < 9:
            continue

        if parts[2] != "locus":
            continue

        chrom = parts[0]
        start = int(parts[3]) - 1   # BED is 0-based
        end = parts[4]
        info = parts[8]

        if "transcripts=" not in info:
            continue

        transcripts = (
            info.split("transcripts=")[1]
            .split(";")[0]
            .split(",")
        )

        locus_data.append((chrom, start, end, transcripts))

# -----------------------------
# choose representative protein per locus
# -----------------------------
kept = 0

with open(out_fa, "w") as fa_out, open(out_bed, "w") as bed_out:
    for chrom, start, end, transcripts in locus_data:

        # transcripts that survived AGAT
        candidates = [t for t in transcripts if t in protein_ids]
        if not candidates:
            continue

        # priority 1: SEC (miniprot / secreted)
        chosen = None
        for t in candidates:
            if t.startswith("SEC"):
                chosen = t
                break

        # priority 2: longest protein
        if chosen is None:
            chosen = max(candidates, key=lambda x: len(proteins[x].seq))

        # write FASTA
        SeqIO.write(proteins[chosen], fa_out, "fasta")

        # write BED
        bed_out.write(f"{chrom}\t{start}\t{end}\t{chosen}\n")

        kept += 1

print(f"[DONE] Wrote {kept} locus-representative proteins")
print(f"[DONE] FASTA: {out_fa}")
print(f"[DONE] BED:   {out_bed}")


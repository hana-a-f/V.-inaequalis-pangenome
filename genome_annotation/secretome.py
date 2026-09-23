#!/usr/bin/env python3

import sys
from Bio import SeqIO

def main(signalp_fasta, tm_fasta, output_fasta):
    # Read transmembrane protein IDs into a set
    tm_ids = set(rec.id for rec in SeqIO.parse(tm_fasta, "fasta"))

    # Filter signal peptide proteins that are NOT in transmembrane set
    signalp_records = []
    for rec in SeqIO.parse(signalp_fasta, "fasta"):
        if rec.id not in tm_ids:
            signalp_records.append(rec)

    # Write the secretome proteins to the output FASTA
    SeqIO.write(signalp_records, output_fasta, "fasta")

    print(f"Secretome predicted: {len(signalp_records)} proteins written to {output_fasta}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python secretome.py signalp_proteins.fa transmembrane_proteins.fa secretome_proteins.fa")
        sys.exit(1)

    signalp_fasta = sys.argv[1]
    tm_fasta = sys.argv[2]
    output_fasta = sys.argv[3]

    main(signalp_fasta, tm_fasta, output_fasta)



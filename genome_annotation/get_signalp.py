import sys
from Bio import SeqIO

#  signalp.gff, protein.fasta, new fasta name


def filter_fasta(gff_file, fasta_file, output_fasta):
    # Step 1: Extract IDs from the GFF file
    gene_ids = set()
    with open(gff_file, "r") as gff:
        for line in gff:
            if line.strip():  # Ignore empty lines
                gene_id = line.split("\t")[0]  # First column
                gene_ids.add(gene_id)

    # Step 2: Filter the FASTA file
    with open(output_fasta, "w") as output:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in gene_ids:
                SeqIO.write(record, output, "fasta")

    print(f"Filtered sequences saved to {output_fasta}")

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 4:
        print("Usage: python get_signalp.py input.gff proteins.fasta filtered_proteins.fasta")
        sys.exit(1)

    # Assign arguments to variables
    gff_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_fasta = sys.argv[3]

    # Call the filtering function
    filter_fasta(gff_file, fasta_file, output_fasta)


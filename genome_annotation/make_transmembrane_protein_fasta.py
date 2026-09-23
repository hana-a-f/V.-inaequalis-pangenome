import sys
from Bio import SeqIO

###this is a script which will match the proteins with a transmembrane domain with the gene names in a fasta file##

tmhmm_file = sys.argv[1]
protein_file = sys.argv[2]
output_fasta = sys.argv[3]

#collect the gene ids from the tmhmm file

gene_ids = set()
with open(tmhmm_file, "r") as trans:
	for line in trans:
		gene_id = line.strip().split()[0] #the first column is the gene id
		gene_ids.add(gene_id) #add in the gene id to the set defined above

#filter the fasta protein file with the gene_ids

with open(output_fasta, "w") as outfile:
	for record in SeqIO.parse(protein_file, "fasta"): #seqIO parse creates an iterator that foes through each sequence
		if record.id in gene_ids: #record.id is the sequence id 
			SeqIO.write(record, outfile, "fasta") #this file writes the current record, a sequence that passeed the filtering condition into the outfile in fasta format

print(f"Filtered FASTA file with transmembrane domain proteins written to {output_fasta}")
	

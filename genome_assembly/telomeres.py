from Bio import SeqIO
from Bio.Seq import Seq
import sys
import re

# === CONFIG ===
fasta_path = sys.argv[1]
repeat_motif = "TTAGGG"
rc_motif = str(Seq(repeat_motif).reverse_complement())
window_size = 1000
min_hits = 4

def find_repeat_positions(seq, motif):
    return [m.start() for m in re.finditer(f'(?={motif})', seq)]

with open("telomeres_precise.bed", "w") as outbed:
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        seq_len = len(seq)

        # === 5' END ===
        head = seq[:window_size]
        hits_fwd = find_repeat_positions(head, repeat_motif)
        hits_rev = find_repeat_positions(head, rc_motif)
        all_hits = sorted(hits_fwd + hits_rev)

        if len(all_hits) >= min_hits:
            start = all_hits[0]
            end = all_hits[-1] + len(repeat_motif)
            outbed.write(f"{record.id}\t{start}\t{end}\t{repeat_motif}_telomere\t{len(all_hits)}\t5prime\n")

        # === 3' END ===
        tail = seq[-window_size:]
        hits_fwd = find_repeat_positions(tail, repeat_motif)
        hits_rev = find_repeat_positions(tail, rc_motif)
        all_hits = sorted(hits_fwd + hits_rev)

        if len(all_hits) >= min_hits:
            # Adjust coordinates to whole sequence
            start = seq_len - window_size + all_hits[0]
            end = seq_len - window_size + all_hits[-1] + len(repeat_motif)
            outbed.write(f"{record.id}\t{start}\t{end}\t{repeat_motif}_telomere\t{len(all_hits)}\t3prime\n")

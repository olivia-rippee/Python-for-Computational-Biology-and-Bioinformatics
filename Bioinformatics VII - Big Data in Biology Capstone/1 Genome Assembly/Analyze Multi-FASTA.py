import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics/Bioinformatics VII - Big Data in Biology Capstone/1 Genome Assembly")

from Bio import SeqIO


# Read sequences
# -------------------------
fasta_sequences = list(SeqIO.parse("Data/SPAdes_contigs.fasta", "fasta"))


# Number of sequences
# -------------------------
len(fasta_sequences) # Output: 3


# Length of sequences
# -------------------------
lengths = {}
for sequence in fasta_sequences:
    lengths[sequence.id] = len(sequence.seq)

# ordered = sorted(lengths, key = lengths.get)

for seq_id in lengths:
    print(f"{seq_id}: {lengths[seq_id]} bp")


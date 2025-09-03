import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics/Bioinformatics VII - Big Data in Biology Capstone/1 Genome Assembly/Data")

from Bio import SeqIO

records = SeqIO.parse("8595-1.txt", "fastq")
count = SeqIO.write(records, "Phred_score_8595-1.txt", "qual")
print("Converted %i records" % count)
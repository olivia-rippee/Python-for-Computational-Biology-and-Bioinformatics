import os
import glob
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import pandas as pd

os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics/Bioinformatics VII - Big Data in Biology Capstone/2 Variant Detection and Classification")

SPIKE_START = 21562  # 21563 in 1-based
SPIKE_END = 25384
SPIKE_LENGTH = SPIKE_END - SPIKE_START + 1


def ExtractSpike(seq):
    '''Extract spike region from full genome sequence.'''
    spike_seq = seq[SPIKE_START:SPIKE_END + 1]
    return spike_seq if len(spike_seq) == SPIKE_LENGTH else None

def LabelVariantByAA(spike_seq):
    '''Label variant based on amino acid mutations in spike sequence.'''
    if spike_seq is None:
        return 'Other'

    # Remove gaps and Ns
    spike_clean = spike_seq.replace('-', '').replace('N', '').upper()

    # Translate to AAs
    try:
        spike_protein = str(Seq(spike_clean).translate())
    except Exception:
        return 'Other'

    if len(spike_protein) < 716:
        return 'Other'

    # Omicron: A67V
    if spike_protein[66] == 'V':
        return 'Omicron'

    # Delta: P681R
    if spike_protein[680] == 'R':
        return 'Delta'

    # Alpha: check A570D, P681H, and T716I
    alpha_mutations = [spike_protein[569] == 'D',  # A570D
                       spike_protein[680] == 'H',  # P681H
                       spike_protein[715] == 'I']  # T716I
    if sum(alpha_mutations) >= 1:
        return 'Alpha'

    return 'Other'


variant_data = defaultdict(Counter)

for file_path in glob.glob('Data/UK-Genomes/*/*_A.fasta'):
    date = os.path.basename(os.path.dirname(file_path))  # extract date from folder name

    for record in SeqIO.parse(file_path, 'fasta'):
        spike = ExtractSpike(str(record.seq))
        label = LabelVariantByAA(spike)
        variant_data[date][label] += 1


df = pd.DataFrame.from_dict(variant_data, orient='index').fillna(0).astype(int)
df = df.sort_index()
df_percent = df.div(df.sum(axis=1), axis=0) * 100

for col in ['Omicron', 'Alpha', 'Delta', 'Other']:
    if col not in df_percent.columns:
        df_percent[col] = 0

df_percent = df_percent[['Alpha', 'Delta', 'Omicron', 'Other']]


df_percent.plot(kind='bar', stacked=True, colormap='tab10', figsize=(14, 7))
plt.title('Percentage of SARS-CoV-2 Variants Over Time')
plt.ylabel('Percentage (%)')
plt.xlabel('Date')
plt.xticks(rotation=45)
plt.legend(title='Variant')
plt.tight_layout()
plt.savefig('Data/Variant_Percentage_Over_Time.jpg', format='jpg', dpi=300)
plt.show()

df_percent.to_csv("Data/variant_percentages.csv")

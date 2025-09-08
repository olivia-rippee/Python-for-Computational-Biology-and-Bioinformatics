import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics/Bioinformatics VII - Big Data in Biology Capstone/3 Protein Structure Analysis")

# Calculate the RMSD between a predicted protein model using SWISS 
# and an experimentally determined model from the Protein Data Bank 
# ------------------------------------------------------------------

from prody import *

def printMatch(match):
    print('Chain 1 : {}'.format(match[0]))
    print('Chain 2 : {}'.format(match[1]))
    print('Length : {}'.format(len(match[0])))
    print('Seq identity: {}'.format(match[2]))
    print('Seq overlap : {}'.format(match[3]))
    print('RMSD : {}\n'.format(calcRMSD(match[0], match[1])))

struct1 = parsePDB('Data/6CRX.pdb')
struct2 = parsePDB('Data/6vxx.pdb')

matches = matchChains(struct1,struct2,seqid=75,overlap=80)
for match in matches:
    printMatch(match)


# Calculating RMSD of two chains
# ------------------------------------------------------------------
first_ca = matches[0][0] 
second_ca = matches[0][1]
calcTransformation(first_ca, second_ca).apply(first_ca);
calcRMSD(first_ca,second_ca)



# Merging multiple chains to compute RMSD of an overall structure
# ------------------------------------------------------------------
first_ca = matches[0][0] + matches[4][0] + matches[8][0]
second_ca = matches [0][1] + matches[4][1] + matches[8][1]
calcTransformation(first_ca, second_ca).apply(first_ca);
calcRMSD(first_ca, second_ca)


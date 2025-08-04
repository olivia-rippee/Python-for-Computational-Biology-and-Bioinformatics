import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics I - Hidden Messages in DNA/Data")



# -----------------------------------------------------
# Get Randomized Motifs
# -----------------------------------------------------

import random


def RandomMotifs(Dna, k, t):
    '''Get randomized motifs.
    Input:  A list of strings Dna, and integers k and t.
    Output: RandomMotifs(Dna, k, t).
    '''

    random_motifs = []
    t = len(Dna)
    l = len(Dna[0])
    for i in range(t):
        ran_num = random.randint(0, l-k)
        random_motifs.append(Dna[i][ran_num:ran_num+k])
    return random_motifs


# Example
# --------------
Dna = ("TTACCTTAAC", "GATGTCTGTC", "ACGGCGTTAG", "CCCTAACGAG", "CGTCAGAGGT")
print(RandomMotifs(Dna, k=3, t=5))
# Output: ['ACC', 'GAT', 'TAG', 'TAA', 'AGA']



# -----------------------------------------------------
#  Randomized Motif Search
# -----------------------------------------------------

import random

def RandomizedMotifSearch(Dna, k, t):
    '''Input:  Positive integers k and t, followed by a list of strings Dna
    Output: RandomizedMotifSearch(Dna, k, t)'''
    
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs
       

def RandomMotifs(Dna, k, t):
    random_motifs = []
    t = len(Dna)
    l = len(Dna[0])
    for i in range(t):
        ran_num = random.randint(0, l-k)
        random_motifs.append(Dna[i][ran_num:ran_num+k])
    return random_motifs


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = CountWithPseudocounts(Motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j] / (t+4)) # each nt adds 1
    return profile    


def CountWithPseudocounts(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def Motifs(Profile, Dna):
    k = len(Profile["A"])
    most_probable_kmers = []
    for text in Dna:
        most_probable = ProfileMostProbablePattern(text, k, Profile)
        most_probable_kmers.append(most_probable)
    return most_probable_kmers


def ProfileMostProbablePattern(text, k, profile):
    max_probability = -0.1 #initialize below 0 in case all kmers have probability 0
    most_probable_kmer = ""
    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
        probability = Pr(pattern, profile)  
        if probability > max_probability:
            max_probability = probability
            most_probable_pattern = pattern
    return most_probable_pattern


def Pr(Text, Profile):
    p=1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return  p


def Consensus(Motifs):
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def Score(Motifs):
    consensus=Consensus(Motifs)
    t=len(Motifs)
    k=len(Motifs[0])
    score=0
    for i in range(t):
        for j in range(k):
            if Motifs[i][j]!=consensus[j]:
                score+=1
    return score



# Example
# -------------
Dna = ("CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA")
print(RandomizedMotifSearch(Dna, k=8, t=5))
# Output (random each time but something like): ['TCTCGGGG', 'CCAAGGTG', 'TACAGGCG', 'TTCAGGTG', 'TCCACGTG']



# Run RandomizedMotifSearch 1000 times and select the best motifs
with open("dataset_30307_5.txt", "r") as file:
    k, t = file.readline().split()
    k, t = int(k), int(t)
    DNA = file.readline().split()

all_best_motifs = []
for _ in range(1000):
    best_motifs = RandomizedMotifSearch(DNA, k, t)
    all_best_motifs.append(best_motifs)
best_motifs_overall = min(all_best_motifs, key=Score)

print(" ".join(best_motifs_overall))
# Output: AGATTGGTCATCACC ATACTCGGCACCAAA ATATCGGCTTCCAAA 
#         ATATCGGGCTAAAAA TTATCGGGCACCATC ATATCGGGCATACAA 
#         ATATCGGGCACGCGA ATAAATGGCACCAAA ATATCGCATACCAAA 
#         TGCTCGGGCACCAAA ATCCGGGGCACCAAA ATATCGGGTGGCAAA 
#         ATATCGGGCACCTTC TAATCGGGCACCAAC AATGCGGGCACCAAA 
#         ATATCGTCTACCAAA ATATGCAGCACCAAA ATATCACCCACCAAA 
#         ATATCATACACCAAA ATATTTCGCACCAAA




# -----------------------------------------------------
# Compute the probability that ten randomly selected 15-mers 
# from the ten 600-nucleotide long strings in the Subtle Motif 
# Problem capture at least one implanted 15-mer.
# -----------------------------------------------------

kmer = "ACGATCGATCGATAA" #random 15-mer

# 1. Compute p1: probability of not capturing the implanted k-mer (15-mer) in one string.
p1 = (600-len(kmer)) / (600-len(kmer)+1)

# 2. Notice for the entire problem we have to deal with ten similar 
# cases, i.e. you have to multiply p1 * p2... *p10, where 
# p1 = p2 = ... = p10. So you just compute p1 to the 10th power:
pC = pow(p1,10)

# 3. Then you just compute the 'opposite' probability, i.e. 
# the probability that from ten 600-length nucleotide string, 
# we capture at least one implanted 15-mer! 
Answer1 = 1 - pC
# Output: 0.01693439692911025



# -----------------------------------------------------
# Compute the probability that ten randomly selected 15-mers 
# from ten 600-nucleotide long strings (as in the Subtle
# Motif Problem) capture at least two implanted 15-mers. 
# -----------------------------------------------------

p_first_string = 1/586
p_not_first_string = 585/586
p2 = p_first_string * (p_not_first_string ** 9)
Answer2 = Answer1 - (p2 * 10)
print(Answer2)
# Output: 0.0001298567056762373



# -----------------------------------------------------
#  Normalize Probabilities
# -----------------------------------------------------
# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities

def Normalize(Probabilities):
    total_probability = sum(Probabilities.values())
    normalized_probabilities = {key: value / total_probability for key, value in Probabilities.items()}
    return normalized_probabilities




# Examples
# --------------
Probabilities = {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1}
print(Normalize(Probabilities))
# Output: {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}


Probabilities = {'A': 0.15, 'B': 0.6, 'C':0.225, 'D':0.225, 'E':0.3}
print(Normalize(Probabilities))
# Output {'A': 0.09999999999999999, 'B': 0.39999999999999997, 'C': 0.15, 'D': 0.15, 'E': 0.19999999999999998}


Probabilities = {'A': 0.45, 'B': 0.63, 'C':0.09, 'D': 0.27, 'E': 0.36 }
print(Normalize(Probabilities))             
# Output: {'A': 0.24999999999999997, 'B': 0.35, 'C': 0.04999999999999999, 'D': 0.15, 'E': 0.19999999999999996}                 
# Rounded: 0.25 0.35 0.05 0.15 0.20                 
          
           
          
            
           
# -----------------------------------------------------
#  Weighted Dice Roll
# -----------------------------------------------------

def WeightedDie(Probabilities):
    p = random.uniform(0, 1)
    
    # Iterate through k-mers and their probabilities
    for kmer, probability in Probabilities.items():
        p -= probability
        if p <= 0:
            return kmer


# Example
# ---------------
Probabilities = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
print(WeightedDie(Probabilities))
# Output: one of the four bases, with equal chance




# -----------------------------------------------------
# Randomly generate a k-mer from Text whose probabilities are generated from profile
# -----------------------------------------------------

import random

def Pr(Text, Profile):
    p=1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return  p

def Normalize(Probabilities):
    total_probability = sum(Probabilities.values())
    normalized_probabilities = {key: value / total_probability for key, value in Probabilities.items()}
    return normalized_probabilities

def WeightedDie(Probabilities):
    p = random.uniform(0, 1)
    for kmer, probability in Probabilities.items():
        p -= probability
        if p <= 0:
            return kmer
        
# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)

def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {} 
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)


# Example
# -------------
Text = "AAACCCAAACCC"
profile = {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}
k = 2
random_kmer = ProfileGeneratedString(Text, profile, k)
print(random_kmer)
# Output: one of AA, AC, CA, or CC




# -----------------------------------------------------
# Gibbs Sampling
# -----------------------------------------------------

import random

# Input:  Integers k, t, and N, followed by a collection of strings Dna
# Output: GibbsSampler(Dna, k, t, N)

def GibbsSampler(Dna, k, t, N):
    Motifs=RandomMotifs(Dna, k, t)
    BestMotifs=Motifs
    for j in range(N):
        i=random.randint(1,t)
        text=Motifs.pop(i-1)
        profile=ProfileWithPseudocounts(Motifs)
        kmer=ProfileMostProbableKmer(text, k, profile)
        Motifs.insert(i-1,kmer)
        if Score(Motifs)<Score(BestMotifs):
            BestMotifs=Motifs
    return BestMotifs


def RandomMotifs(Dna, k, t):
    random_motifs = []
    t = len(Dna)
    l = len(Dna[0])
    for i in range(t):
        ran_num = random.randint(0, l-k)
        random_motifs.append(Dna[i][ran_num:ran_num+k])
    return random_motifs


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = CountWithPseudocounts(Motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j] / (t+4)) # each nt adds 1
    return profile    


def CountWithPseudocounts(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def Motifs(Profile, Dna):
    k = len(Profile["A"])
    most_probable_kmers = []
    for text in Dna:
        most_probable = ProfileMostProbablePattern(text, k, Profile)
        most_probable_kmers.append(most_probable)
    return most_probable_kmers


def ProfileMostProbableKmer(text, k, profile):
    max_probability = -0.1 #initialize below 0 in case all kmers have probability 0
    most_probable_kmer = ""
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        probability = Pr(kmer, profile)  
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer
    return most_probable_kmer


def Pr(Text, Profile):
    p=1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return  p


def Consensus(Motifs):
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def Score(Motifs):
    consensus=Consensus(Motifs)
    t=len(Motifs)
    k=len(Motifs[0])
    score=0
    for i in range(t):
        for j in range(k):
            if Motifs[i][j]!=consensus[j]:
                score+=1
    return score



# Example
# ------------
Dna = ("CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA")
print(GibbsSampler(Dna, k=8, t=5, N=100))
# Output: something like ['TCTCGGGG', 'CCAAGGTG', 'TACAGGCG', 'TTCAGGTG', 'TCCACGTG']







# Gibbs Sampler Simulation
# -------------------------

import random

def GibbsSampler(Dna, k, t, N):
    """Implements the Gibbs Sampler algorithm."""
    Motifs = []
    for dna_string in Dna:
        start_index = random.randint(0, len(dna_string) - k)
        Motifs.append(dna_string[start_index:start_index + k])
    BestMotifs = Motifs
    for _ in range(N):
        i = random.randint(0, t - 1)
        Profile = ProfileWithPseudocounts(Motifs[:i] + Motifs[i+1:])  # Profile without Motifi
        Motifi = ProfileRandomlyGeneratedKmer(Dna[i], k, Profile)
        Motifs = Motifs[:i] + [Motifi] + Motifs[i+1:]
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def ProfileWithPseudocounts(Motifs):
    """Creates a profile matrix with pseudocounts."""
    k = len(Motifs[0])
    t = len(Motifs)
    profile = {'A': [1] * k, 'C': [1] * k, 'G': [1] * k, 'T': [1] * k}
    for motif in Motifs:
        for i, nucleotide in enumerate(motif):
            profile[nucleotide][i] += 1
    for nucleotide in profile:
        for i in range(k):
            profile[nucleotide][i] /= (t + 4)
    return profile

def ProfileRandomlyGeneratedKmer(Text, k, Profile):
    """Generates a k-mer from Text based on probabilities in Profile."""
    probabilities = []
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i + k]
        probability = 1
        for j, nucleotide in enumerate(kmer):
            probability *= Profile[nucleotide][j]
        probabilities.append(probability)
    # Normalize probabilities
    total_probability = sum(probabilities)
    probabilities = [p / total_probability for p in probabilities]
    # Choose a k-mer based on probabilities
    random_index = random.choices(range(len(probabilities)), probabilities)[0]
    return Text[random_index:random_index + k]

def Score(Motifs):
    """Calculates the score of a set of motifs."""
    k = len(Motifs[0])
    t = len(Motifs)
    consensus = Consensus(Motifs)
    score = 0
    for motif in Motifs:
        for i in range(k):
            if motif[i] != consensus[i]:
                score += 1
    return score

def Consensus(Motifs):
    """Generates the consensus string from a set of motifs."""
    k = len(Motifs[0])
    t = len(Motifs)
    consensus = ''
    for i in range(k):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in Motifs:
            counts[motif[i]] += 1
        max_count = 0
        max_nucleotide = ''
        for nucleotide, count in counts.items():
            if count > max_count:
                max_count = count
                max_nucleotide = nucleotide
        consensus += max_nucleotide
    return consensus

with open("dataset_30309_11.txt", "r") as file:
    k, t, N = file.readline().split()
    k, t, N = int(k), int(t), int(N)
    Dna = file.readline().split()
    
all_best_motifs = []
for _ in range(20):
    best_motifs = GibbsSampler(Dna, k, t, N)
    all_best_motifs.append(best_motifs)
best_motifs_overall = min(all_best_motifs, key=Score)

print(" ".join(best_motifs_overall))
# Output: TTCCCCACGTAGTCC TTCTAGCCCTAATTC TCGGTCCCCTAATTC 
#         TTCCTAGTCTAATTC AGCCTCCCCTAATTT TTCCTCCCCTGCATC 
#         TTCGCACCCTAATTC TTCCTCCCCTAAGAT TTCCCGTCCTAATTC 
#         TTCCGAACCTAATTC TTCCTCCCCTAGGAC TTTTACCCCTAATTC 
#         TTCCTCTGATAATTC TTCCTCCCCATTTTC TTCCTCCCTGCATTC 
#         GTCCTCCCCTAATGG TTCCTCTGTTAATTC TTCCTCCAACAATTC 
#         CGTCTCCCCTAATTC TTCCTGGACTAATTC




# -----------------------------------------------------
# Quiz question: What are the 3-mers after one iteration of RandomizedMotifSearch?  
# In other words, what are the 3-mers Motifs(Profile(Motifs), Dna)?
# -----------------------------------------------------

import random

def CountWithPseudocounts(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = [1] * k  # pseudocount = 1
    for motif in Motifs:
        for j, symbol in enumerate(motif):
            count[symbol][j] += 1
    return count

def ProfileWithPseudocounts(Motifs):
    count = CountWithPseudocounts(Motifs)
    t = len(Motifs) + 4  # include pseudocounts
    profile = {}
    for symbol in "ACGT":
        profile[symbol] = [count[symbol][j] / t for j in range(len(Motifs[0]))]
    return profile

def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return p

def ProfileMostProbablePattern(text, k, profile):
    max_probability = -1
    most_probable_pattern = text[0:k]
    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
        probability = Pr(pattern, profile)
        if probability > max_probability:
            max_probability = probability
            most_probable_pattern = pattern
    return most_probable_pattern

def Motifs(Profile, Dna):
    k = len(Profile["A"])
    return [ProfileMostProbablePattern(text, k, Profile) for text in Dna]

def RandomizedMotifSearch_One_Iteration(Dna, initial_motifs):
    profile = ProfileWithPseudocounts(initial_motifs)
    new_motifs = Motifs(profile, Dna)
    return new_motifs




# Examples
# ------------------
Dna = ["AAGCCAAA", "AATCCTGG", "GCTACTTG", "ATGTTTTG"]
initial_motifs = ["CCA", "CCT", "CTT", "TTG"]
result = RandomizedMotifSearch_One_Iteration(Dna, initial_motifs)
print(" ".join(result))
# Output: CCA CCT CTT TTT



Dna = ["ATGAGGTC", "GCCCTAGA", "AAATAGAT", "TTGTGCTA"]
initial_motifs = ["GTC", "CCC", "ATA", "GCT"]
result = RandomizedMotifSearch_One_Iteration(Dna, initial_motifs)
print(" ".join(result))
# Output: GTC GCC ATA GCT




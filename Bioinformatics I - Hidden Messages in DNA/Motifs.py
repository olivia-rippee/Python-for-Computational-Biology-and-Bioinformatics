import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics I - Hidden Messages in DNA/Data")



# What is the expected number of occurrences of a 9-mer in 500 random DNA strings, 
# each of length 1000? Assume that the sequences are formed by selecting each 
# nucleotide (A, C, G, T) with the same probability (0.25).

num_strings = 500
string_length = 1000
kmer_length = 9
probability_per_nucleotide = 0.25

probability_kmer = probability_per_nucleotide ** kmer_length
num_positions_per_string = string_length - kmer_length + 1
expected_occurrences = num_positions_per_string * num_strings * probability_kmer
print(expected_occurrences)
# Output: 1.89208984375




# -----------------------------------------------------
# Brute Force Method
# -----------------------------------------------------
# Enumerate all (k, d)-motifs in a collection of strings.

def HammingDist(p, q):
    """
    Calculate the Hamming distance between two sequences p and q.
    The Hamming distance between two strings is the number of positions at which the corresponding symbols are different.
    
    Parameters:
    p (sequence): The first sequence.
    q (sequence): The second sequence.
    
    Returns:
    int: The Hamming distance between p and q.
    
    Note:
    The lengths of the two sequences do not need to be the same. The function compares elements at corresponding positions until the end of the shorter sequence is reached.
    """
    # Use the zip function to pair up elements from both sequences and generate a list of booleans indicating whether the elements are unequal
    # Sum the booleans (True values count as 1) to calculate the total number of differences, which represents the Hamming distance
    return sum(p != q for p, q in zip(p, q))

def ImmediateNeighbor(Pattern):
    """
    Finds all immediate neighbor patterns for a given pattern.

    Parameters:
    Pattern (str): The input pattern string.

    Returns:
    list: A list containing all immediate neighbor patterns.
    """
    # Initialize the neighborhood list with the input pattern itself
    Neighborhood = [Pattern]
    # Iterate over each position in the pattern
    for i in range(len(Pattern)):
        # Iterate over all nucleotides
        for j in ["A", "C", "G", "T"]:
            # If the current nucleotide is not equal to j, replace it and add to the neighborhood list
            if Pattern[i] != j:
                Neighborhood.append(Pattern[:i] + j + Pattern[i+1:])
    return Neighborhood

def Neighbors(Pattern, d):
    """
    Finds all neighbors of a given pattern within a specified Hamming distance.

    Parameters:
    Pattern (str): The input pattern string.
    d (int): Maximum Hamming distance.

    Returns:
    list: A list containing all patterns that are within the maximum Hamming distance d of the input pattern.
    """
    # Initialize the neighborhood with the input pattern itself
    Neighborhood = [Pattern]
    # For each distance from 0 to d, find the immediate neighbors of the current patterns
    for j in range(d):
        # For each current neighbor, find its immediate neighbors and add them to the neighborhood
        for Neighbor in Neighborhood:
            Neighborhood = Neighborhood + ImmediateNeighbor(Neighbor)
        # Remove duplicate neighbors to ensure uniqueness in the neighborhood
        Neighborhood = list(set(Neighborhood))
    # Remove any remaining duplicates before returning the final neighborhood
    Neighborhood = list(set(Neighborhood))
    return Neighborhood

def MultiNeighbors(k, d, Dna):
    """
    Calculate the d-neighborhood set of k-mers for each DNA sequence.

    Parameters:
    k -- int, the value of k in k-mers
    d -- int, the maximum distance for the neighborhood
    Dna -- list, a list of input DNA sequences

    Returns:
    Neighborhood -- dict, the d-neighborhood set of k-mers for each DNA sequence
    """
    # Initialize an empty dictionary to store the neighborhood sets for each DNA sequence
    Neighborhood = {}
    # Iterate over each DNA sequence in the input
    for String in Dna:
        # Initialize an empty list for the current DNA sequence's neighborhood set
        Neighborhood[String] = []
        # Iterate over all possible k-mers in the current DNA sequence
        for i in range(len(String) - k + 1):
            # Extract the k-mer at the current position
            Pattern = String[i:i+k]
            # Add the d-neighborhood of the current k-mer to the neighborhood set
            Neighborhood[String] = Neighborhood[String] + Neighbors(Pattern, d)
        # Remove duplicates from the neighborhood set
        Neighborhood[String] = list(set(Neighborhood[String]))
    # Return the neighborhood sets for each DNA sequence
    return Neighborhood


def MotifEnumeration(Dna, k, d):
    """
    Main function: enumerates all common motifs in DNA sequences given certain conditions.

    Parameters:
    Dna - A list containing multiple DNA sequences.
    k - The length of the motif.
    d - The allowed Hamming distance.

    Returns:
    A list containing all motifs that satisfy the conditions.
    """
    # Initialize the motifs list
    Motifs = []
    # Build a dictionary of k-mers and their neighborhoods for each DNA sequence
    Neighborhood = MultiNeighbors(k, d, Dna)
    
    # Iterate through each DNA sequence
    for String in Dna:
        # Iterate through all possible k-mers and their neighborhoods for the current DNA sequence
        for Pattern in Neighborhood[String]:
            # Ensure the current pattern has not been added to the motifs list
            if Pattern not in Motifs:
                # Initialize a counter to count the number of matching DNA sequences
                Count = 0
                # Iterate through each DNA sequence
                for String2 in Dna:
                    # Iterate through all possible k-mers for the current DNA sequence
                    for i in range(len(String2) - k + 1):
                        # If the Hamming distance between the current pattern and the k-mer is no more than d, increment the counter
                        if HammingDist(Pattern, String2[i:i+k]) <= d:
                            Count += 1
                            break
                # If all DNA sequences match, add the current pattern to the motifs list
                if Count == len(Dna):
                    Motifs.append(Pattern)
    
    # Remove duplicate motifs
    Motifs = list(set(Motifs))
    # Return all motifs that satisfy the conditions
    return(Motifs)




# Examples
# -----------
with open("dataset_30302_8.txt", "r") as file:
    k, d = file.readline().split()
    k, d = int(k), int(d)
    Dna = file.readline().split()
print(*MotifEnumeration(Dna, k, d))
# Output: TCCGA GGAAC GGTAC GGGAC GGCAC






# -----------------------------------------------------
# Count occurrences of each nucleotide in each column of the motif matrix
# -----------------------------------------------------
# Count(Motifs) stores the number of times that nucleotide i appears in column j of Motifs
# Input:  A set of kmers Motifs
# Output: Count(Motifs)

def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count



# -----------------------------------------------------
# Profile matrix
# -----------------------------------------------------
# element (i, j) = frequency of the i-th nucleotide in the j-th column of the motif matrix
# number of occurrences of the i-th nucleotide divided by t the number of nucleotides in the column

def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.

def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = Count(Motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j] / t)
    return profile



# Example
# --------------
Motifs = ("AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG")
print(Count(Motifs))

# Output: {'A': [0.2, 0.4, 0.2, 0.0, 0.0, 0.4],
#          'C': [0.4, 0.2, 0.8, 0.4, 0.0, 0.0],
#          'G': [0.2, 0.2, 0.0, 0.4, 0.2, 0.2],
#          'T': [0.2, 0.2, 0.0, 0.2, 0.8, 0.4]}



# -----------------------------------------------------
# Consensus string (most popular nucleotides in each column of the motif matrix)
# -----------------------------------------------------

def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
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


# Example
# --------------
Motifs = ("AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG")
print(Consensus(Motifs))

# Output: CACCTA







# -----------------------------------------------------
# Compute the score
# -----------------------------------------------------
# Sum the number of symbols in the j-th column of Motifs that
# do not match the symbol in position j of the consensus string


def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
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

def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
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


Motifs = ("AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG")
print(Score(Motifs))

# Output: 14




# -----------------------------------------------------
# Scoring via "entropy"
# -----------------------------------------------------

# -------------------
# From profile matrix
# -------------------
import numpy as np

def ProfileToMotifEntropy(Profile):
    entropy = 0
    for i in Profile:
        entropy += -i*(np.log2(i))
    return float(entropy)

    
# Examples
# -----------
Profile = [0.2,0.2,0.6,0.2,0.1,0.7,0.9,0.1,0.9,0.1,0.9,0.1,0.1,0.4,0.5,0.1,0.1,0.8,0.1,0.2,0.7,0.3,0.4,0.3,0.6,0.4]
# to simplify the calculation, ignore all 0s and 1s (they don't impact the result)
ProfileToMotifEntropy(Profile)



ProfileA = [0.5, 0.5]         # [0.5, 0, 0, 0.5]
ProfileB = [0.25, 0.25, 0.25, 0.25]
ProfileC = [1]                      # [0, 0, 0, 1]
ProfileD = [0.25, 0.5, 0.25]        # [0.25, 0, 0.5, 0.25]

ProfileToMotifEntropy(ProfileA) # Output: 1.0
ProfileToMotifEntropy(ProfileB) # Output: 2.0
ProfileToMotifEntropy(ProfileC) # Output: 0.0
ProfileToMotifEntropy(ProfileD) # Output: 1.5
# C < A < D < B




# -------------------
# From list of motifs
# -------------------

def CountMotifPercent(Motifs):
    count = {}
    columns = []
    for i in range(len(Motifs[0])):
        columns.append([motif[i] for motif in Motifs])
    for i in range(len(columns)):
        count[i] = {'A': columns[i].count('A')/len(columns[i]), 
                    'C': columns[i].count('C')/len(columns[i]), 
                    'G': columns[i].count('G')/len(columns[i]), 
                    'T': columns[i].count('T')/len(columns[i])}
    return count


import math

def MotifEntropy(Motifs):
    entropy = 0
    percents = CountMotifPercent(Motifs)
    for i in range(len(percents)):
        for nucleotide in percents[i]:
            if percents[i][nucleotide] != 0:
                entropy += percents[i][nucleotide] * math.log2(percents[i][nucleotide])
    return -entropy



# Example
# -----------
Motifs = ("TCGGGGGTTTTT", "CCGGTGACTTAC", "ACGGGGATTTTC", "TTGGGGACTTTT", 
          "AAGGGGACTTCC", "TTGGGGACTTCC", "TCGGGGATTCAT", "TCGGGGATTCCT",
          "TAGGGGAACTAC", "TCGGGTATAACC")

print(MotifEntropy(Motifs))
# Output: 9.916290005356972




# -----------------------------------------------------
# Calculate Distance Between Pattern And Strings
# -----------------------------------------------------

def HammingDistance(pattern1, pattern2):
    return sum(c1 != c2 for c1, c2 in zip(pattern1, pattern2))

def DistanceBetweenPatternAndStrings(pattern, dna):
    k = len(pattern)
    distance = 0
    for text in dna:
        hamming_distance = float('inf')  # Initialize with infinity
        for i in range(len(text) - k + 1):
            pattern_prime = text[i:i + k]
            current_distance = HammingDistance(pattern, pattern_prime)
            hamming_distance = min(hamming_distance, current_distance)
        distance += hamming_distance
    return distance


# Examples
# ----------
pattern = "AAA"
DNA = ['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT']
print(DistanceBetweenPatternAndStrings(pattern, DNA))
# Output: 5


with open("dataset_30312_1.txt", "r") as file:
    lines = file.readlines()
    pattern = lines[0][:-1]
    DNA = lines[1].split(' ')
    DNA[-1] = DNA[-1][:-1]
print(DistanceBetweenPatternAndStrings(pattern, DNA))
# Output: 65




# -----------------------------------------------------
# Median String
# -----------------------------------------------------
# find a k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern


def HammingDistance(pattern1, pattern2):
    return sum(c1 != c2 for c1, c2 in zip(pattern1, pattern2))

def MedianString(k, dna):
    distance = float('inf')
    pattern = ""

    for i in range(len(dna[0]) - k + 1):
        current_pattern = dna[0][i:i+k]
        current_distance = sum(min(HammingDistance(current_pattern, seq[i:i+k]) for seq in dna) for i in range(len(dna[0]) - k + 1))
        if current_distance < distance:
            distance = current_distance
            pattern = current_pattern
    return pattern


# Examples
# ----------
k = 3
DNA = ("AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTTCGGGACAG")
print(MedianString(k, DNA))
# Output: GAC


with open("dataset_30304_9.txt", "r") as file:
    k = int(file.readline())
    DNA = file.readline().split()
print(MedianString(k, DNA))
# Output: TCCTTT





# -----------------------------------------------------
# Find all Median Strings of length k
# -----------------------------------------------------

from itertools import product

def HammingDistance(p, q):
    return sum(pc != qc for pc, qc in zip(p, q))

def DistanceBetweenPatternAndStrings(pattern, dna):
    k = len(pattern)
    return sum(
        min(HammingDistance(pattern, dna_seq[i:i+k]) for i in range(len(dna_seq) - k + 1))
        for dna_seq in dna
    )

def MedianStrings(k, dna):
    all_kmers = [''.join(p) for p in product('ACGT', repeat=k)]
    min_distance = float('inf')
    median_kmers = []

    for kmer in all_kmers:
        dist = DistanceBetweenPatternAndStrings(kmer, dna)
        if dist < min_distance:
            min_distance = dist
            median_kmers = [kmer]
        elif dist == min_distance:
            median_kmers.append(kmer)
    return median_kmers




# Example: Which of these 7-mers are median strings for this motif matrix?
# ---------------------
DNA = ["CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC",
    "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC",
    "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"]
k=7
candidates = ["GATGAGT", "TCTGAAG", "AACGCTG", "GTCAGCG", "AATCCTA", "GTAGGAA"]
medians = MedianStrings(k, DNA)

print("Median strings found:")
print(medians)
print("\nCandidates that are median strings:")
for c in candidates:
    if c in medians:
        print(c)
# Output: AATCCTA GTAGGAA




# -----------------------------------------------------
# Probability of bases
# -----------------------------------------------------

# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    # insert your code here
    p=1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return  p



# Examples
# ----------------
Profile = {'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0], 
           'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6], 
           'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0], 
           'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]}

Text = "TCGTGGATTTCC"
Pr(Text, Profile)
# Output: 0.0


Text = "ACGGGGATTACC"
Pr(Text, Profile)
# Ouput: 0.0008398080000000002



Profile =  {'A': [0.4, 0.3, 0.0, 0.1, 0.0, 0.9], 
            'C': [0.2, 0.3, 0.0, 0.4, 0.0, 0.1], 
            'G': [0.1, 0.3, 1.0, 0.1, 0.5, 0.0], 
            'T': [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]}
print(Pr("GAGCTA", Profile))
# Output: 0.0054

print(Pr("CAGTGA", Profile))
# Output: 0.0108

print(Pr("TCGGTA", Profile))
# Output: 0.00405




# -----------------------------------------------------
# Compute  probability of every k-mer in a string Text 
# and find a profile-most probable k-mer
# -----------------------------------------------------

def Pr(Text, Profile):
    # insert your code here
    p=1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return  p

# The profile matrix assumes that the first row corresponds to A, the second corresponds to C,
# the third corresponds to G, and the fourth corresponds to T.
# You should represent the profile matrix as a dictionary whose keys are 'A', 'C', 'G', and 'T' and whose values are lists of floats
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




# Examples
# ----------
Text = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
k = 5
Profile =  {'A': [0.2, 0.2, 0.3, 0.2, 0.3], 
            'C': [0.4, 0.3, 0.1, 0.5, 0.1], 
            'G': [0.3, 0.3, 0.5, 0.2, 0.4], 
            'T': [0.1, 0.2, 0.1, 0.1, 0.2]}
print(ProfileMostProbableKmer(Text, k, Profile))
# Output: CCGAG



with open("dataset_30305_3.txt", 'r') as file:
    lines = file.read().strip().split('\n')
    Text = lines[0].strip()
    k = int(lines[1].strip())
    
    profile_lines = lines[2:6]
    Profile = {}
    bases = ['A', 'C', 'G', 'T']
    for base, line in zip(bases, profile_lines):
        Profile[base] = list(map(float, line.strip().split()))

print(ProfileMostProbableKmer(Text, k, Profile))
# Output: AGCTCGCAGCAATT




# -----------------------------------------------------
# Greedy Motif Search
# -----------------------------------------------------
# heavily influenced by zeros

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

def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = Count(Motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j] / t)
    return profile

def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
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
    # insert your code here
    p=1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return  p

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs



# Examples
# -----------
DNA = ("GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG")
print(*GreedyMotifSearch(DNA, k=3, t=5))
# Output: CAG CAG CAA CAA CAA


with open ("dataset_30305_5.txt", "r") as file:
    k, t = file.readline().split()
    k, t = int(k), int(t)
    DNA = file.readline().split()
print(*GreedyMotifSearch(DNA, k, t))
# Output: CGCTCATGGTGG CAGATCATGTTT GTTCATGTGCGC etc




# Example
# -----------
Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG", "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC", "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC", "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG", "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA", "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA", "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG", "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG", "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]
t = len(Dna)
k = 15
Motifs = GreedyMotifSearch(Dna, k, t)
print(*Motifs)
print(Score(Motifs))

# Output: GTTAGGGCCGGAAGT CCGATCGGCATCACT ACCGTCGATGTGCCC 
#         GGGTCAGGTATATTT GTGACCGACGTCCCC CTGTTCGCCGGCAGC 
#         CTGTTCGATATCACC GTACATGTCCAGAGC GCGATAGGTGAGATT 
#         CTCATCGCTGTCATC
# 64




# -----------------------------------------------------
# Use pseudocounts to remove zeros (all kmers possible)
# -----------------------------------------------------
# Add one to all cells based on Laplace’s Rule of Succession
# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)

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


# Example
# ---------------
Motifs = ("AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG")
print(CountWithPseudocounts)
# Output: {'A': [2, 3, 2, 1, 1, 3], 
#          'C': [3, 2, 5, 3, 1, 1], 
#          'G': [2, 2, 1, 3, 2, 2], 
#          'T': [2, 2, 1, 2, 5, 3]}



# -----------------------------------------------------
# Profile matrix using pseudocounts
# -----------------------------------------------------
# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)

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


# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
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



# Example
# --------------
Motifs = ("AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG")
print(ProfileWithPseudocounts(Motifs))
# Output: 
#{'A': [0.2222222222222222, 0.3333333333333333, 0.2222222222222222, 0.1111111111111111, 0.1111111111111111, 0.3333333333333333], 
#'C': [0.3333333333333333, 0.2222222222222222, 0.5555555555555556, 0.3333333333333333, 0.1111111111111111, 0.1111111111111111], 
#'G': [0.2222222222222222, 0.2222222222222222, 0.1111111111111111, 0.3333333333333333, 0.2222222222222222, 0.2222222222222222], 
#'T': [0.2222222222222222, 0.2222222222222222, 0.1111111111111111, 0.2222222222222222, 0.5555555555555556, 0.3333333333333333]}







# -----------------------------------------------------
# Greedy motif search using pseudocounts
# -----------------------------------------------------

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


# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
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


# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)

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
    # insert your code here
    p=1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return  p



# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearchWithPseudocounts(Dna, k, t)

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs




# Examples
# --------------
Dna = ("GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG")
print(*GreedyMotifSearchWithPseudocounts(Dna, k=3, t=5))
# Output: TTC ATC TTC ATC TTC



Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",
       "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",
       "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",
       "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",
       "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",
       "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",
       "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",
       "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",
       "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",
       "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

Motifs = GreedyMotifSearchWithPseudocounts(Dna, k=15, t=len(Dna))
print(*Motifs)
print(Score(Motifs))

# Output: GGACTTCAGGCCCTA GGTCAAACGACCCTA GGACGTAAGTCCCTA GGATTACCGACCGCA 
#         GGCCGAACGACCCTA GGACCTTCGGCCCCA GGACTTCTGTCCCTA GGACTTTCGGCCCTG 
#         GGACTAACGGCCCTC GGACCGAAGTCCCCG
# 35



with open("dataset_30306_9.txt", "r") as file:
    k, t = file.readline().split()
    k, t = int(k), int(t)
    DNA = file.readline().split()
print(*GreedyMotifSearchWithPseudocounts(DNA, k, t))
# Output: GACCCATCGTAT AAGCCATCGTTT GACCCATCGTCT etc






# -----------------------------------------------------
# Find profile-most probable k-mers in each string of a list
# -----------------------------------------------------
# Input:  A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna)

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




# Example 
# --------------
Profile = {'A': [0.8, 0.0, 0.0, 0.2], 'C': [0.0, 0.6, 0.2, 0.0], 
           'G': [0.2, 0.2, 0.8, 0.0], 'T': [0.0, 0.2, 0.0, 0.8]}
Dna = ("TTACCTTAAC", "GATGTCTGTC", "ACGGCGTTAG", "CCCTAACGAG", "CGTCAGAGGT")
print(Motifs(Profile, Dna))
# Output: ['ACCT', 'ATGT', 'GCGT', 'ACGA', 'AGGT']





# -----------------------------------------------------
# Get Randomized Motifs
# -----------------------------------------------------

import random

# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)

def RandomMotifs(Dna, k, t):
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

# Input:  Positive integers k and t, followed by a list of strings Dna
# Output: RandomizedMotifSearch(Dna, k, t)

def RandomizedMotifSearch(Dna, k, t):
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




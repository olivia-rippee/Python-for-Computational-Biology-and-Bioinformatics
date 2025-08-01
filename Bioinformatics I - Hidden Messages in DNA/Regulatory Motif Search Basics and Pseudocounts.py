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
    '''Compute the number of mismatches (Hamming distance) between two strings.
    
    Input:  Two strings p and q.
    Output: An integer value representing the Hamming Distance between p and q.
    Note: The lengths of the two sequences do not need to be the same. The function 
    compares elements at corresponding positions until the end of the shorter 
    sequence is reached.'''

    # Use the zip function to pair up elements from both sequences and generate a list of booleans indicating whether the elements are unequal
    # Sum the booleans (True values count as 1) to calculate the total number of differences, which represents the Hamming distance
    return sum(p != q for p, q in zip(p, q))

def ImmediateNeighbor(Pattern):
    """Finds all immediate neighbor patterns for a given pattern.
    
    Input: A pattern string.
    Output: A list containing all immediate neighbor patterns."""
    
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
    '''Find the set of all k-mers whose Hamming distance from Pattern does not exceed d.
    
    Input: A pattern string and Maximum Hamming distance d.
    Output: A list containing all patterns that are within the maximum 
    Hamming distance d of the input pattern.
    '''
    
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

    Input: The length of the motifs (k), the maximum distance for the neighborhood (d), and
    a list of input DNA sequences (DNA).

    Output:
    Neighborhood -- dict, the d-neighborhood set of k-mers for each DNA sequence
    """
    # Initialize an empty dictionary to store the neighborhood sets for each DNA sequence
    Neighborhood = {}
    # Iterate over each DNA sequence in the input
    for String in DNA:
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
    """Enumerates all common motifs in DNA sequences given certain conditions.

    Input: The length of the motifs (k), the maximum distance for the neighborhood (d), and
    a list of input DNA sequences (DNA)..
    Output: A list containing all motifs that satisfy the conditions.
    """

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


def Count(Motifs):
    '''Calculate the number of times that nucleotide i appears in column j of Motifs.
    
    Input:  A set of kmers Motifs
    Output: Count(Motifs)'''
    
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
    '''Calculates an array as the number of occurrences of the i-th nucleotide 
    divided by t the number of nucleotides in the column.
    Element (i, j) is frequency of the i-th nucleotide in the j-th column of 
    the motif matrix.
    
    Input:  A list of kmers Motifs
    Output: the profile matrix of Motifs, as a dictionary of lists.'''
    
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


def Consensus(Motifs):
    '''Input:  A set of kmers Motifs
    Output: A consensus string of Motifs.'''
    
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

def Score(Motifs):
    '''Sum the number of symbols in the j-th column of Motifs that
    do not match the symbol in position j of the consensus string.
    
    Input:  A set of k-mers Motifs
    Output: The score of these k-mers.'''
    
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
# Probability of sequence
# -----------------------------------------------------

def Pr(Text, Profile):
    '''Calculate the probability of the sequence.
    
    Input:  String Text and profile matrix Profile.
    Output: Pr(Text, Profile).'''

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
# Compute probability of every k-mer in a string Text 
# and find a profile-most probable k-mer
# -----------------------------------------------------

def Pr(Text, Profile):
    p=1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return  p


def ProfileMostProbableKmer(text, k, profile):
    '''Compute the probability of every k-mer in a string Text 
    and find a profile-most probable k-mer.
    
    The profile matrix assumes that the first row corresponds to A, 
    the second corresponds to C,the third corresponds to G, and the 
    fourth corresponds to T. You should represent the profile matrix 
    as a dictionary whose keys are 'A', 'C', 'G', and 'T' and whose 
    values are lists of floats.
    '''
    
    max_probability = -0.1       #initialize below 0 in case all kmers have probability 0
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


def GreedyMotifSearch(Dna, k, t):
    '''Input:  A list of kmers Dna, and integers k and t 
    (where t is the number of kmers in Dna).
    Output: GreedyMotifSearch(Dna, k, t).
    Heavily influenced by zeros.
    '''
    
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

def CountWithPseudocounts(Motifs):
    '''Use pseudocounts to remove the heavy influence of zeros on
    the probability of the kmer, since all kmers possible even if 
    the probability is small.
    
    Adds one to all cells based on Laplace’s Rule of Succession.
    
    # Input: A set of kmers Motifs
    # Output: CountWithPseudocounts(Motifs)
    '''
    
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

def ProfileWithPseudocounts(Motifs):
    '''Generate the Profile matrix using pseudocounts.
    
    Input:  A set of kmers Motifs.
    Output: ProfileWithPseudocounts(Motifs).
    '''
    
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
    max_probability = -0.1       #initialize below 0 in case all kmers have probability 0
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


def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    '''Conduct a greedy motif search using pseudocounts.
    Input:  A list of kmers Dna, and integers k and t 
    (where t is the number of kmers in DNA).
    Output: GreedyMotifSearchWithPseudocounts(DNA, k, t).
    '''
    
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(DNA[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(DNA[j], k, P))
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

def Motifs(Profile, Dna):
    '''Find profile-most probable k-mers in each string of a list.
    
    Input:  A profile matrix Profile and a list of strings Dna
    Output: Motifs(Profile, Dna)
    '''
    
    k = len(Profile["A"])
    most_probable_kmers = []
    for text in Dna:
        most_probable = ProfileMostProbablePattern(text, k, Profile)
        most_probable_kmers.append(most_probable)
    return most_probable_kmers


def ProfileMostProbablePattern(text, k, profile):
    max_probability = -0.1     #initialize below 0 in case all kmers have probability 0
    most_probable_pattern = ""
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


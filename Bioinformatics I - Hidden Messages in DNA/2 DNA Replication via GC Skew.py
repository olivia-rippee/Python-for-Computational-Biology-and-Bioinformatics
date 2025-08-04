import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics I - Hidden Messages in DNA/Data")

# -----------------------------------------------------    
# Difference (skew) array for difference in bases
# -----------------------------------------------------

def SkewArray(Genome):
    
    '''Generate the difference (skew) array for difference in bases.
    The difference between G and C is increasing for the forward half-strand 
    and decreasing for the reverse half-strand.
    
    Input:  A String Genome
    Output: The skew array of Genome as a list.'''
    
    skew = [0]
    count_G = 0
    count_C = 0
    
    for i in range(len(Genome)):
        if Genome[i] == 'G':
            count_G += 1
        elif Genome[i] == 'C':
            count_C += 1
        skew.append(count_G - count_C)
    return skew


# Examples
# ------------------
Genome = "CATGGGCATCGGCCATACGCC"
print(*SkewArray(Genome))
# Output: 0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2



# -----------------------------------------------------
# Find position in genome where skew diagram attains a minimum (locate ori)
# -----------------------------------------------------


def MinimumSkew(Genome):
    
    '''Find position in genome where skew diagram attains a minimum (locate ori).
    Skew decreases along the reverse half-strand and increases along the forward half-strand. 
    Skew should achieve a minimum where the reverse half-strand ends and the forward half-strand begins (ori).
    
    Input:  A DNA string Genome
    Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
    '''
    
    positions = [] 
    skew = SkewArray(Genome)
    min_skew = min(skew)
    for i in range(len(skew)):
        if skew[i] == min_skew:
            positions.append(i)
    return positions


def SkewArray(Genome):
    
    '''Generate the difference (skew) array for difference in bases.
    The difference between G and C is increasing for the forward half-strand 
    and decreasing for the reverse half-strand.
    
    Input:  A String Genome
    Output: The skew array of Genome as a list.'''
    
    skew = {}
    skew = [0]
    count_G = 0
    count_C = 0
    for i in range(len(Genome)):
        if Genome[i] == 'G':
            count_G += 1
        elif Genome[i] == 'C':
            count_C += 1
        skew.append(count_G - count_C)
    return skew



# Examples
# -----------------
print(MinimumSkew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"))
# Output: 11 24


print(MinimumSkew(Genome="GATACACTTCCCGAGTAGGTACTG"))
# Output: 12


with open("dataset_30277_10.txt", "r") as file:
    Genome = file.readline().strip()
print(*MinimumSkew(Genome))
# Output: 55673 55898




# -----------------------------------------------------
# Report i for which the skew array of the genome attains a maximum
# -----------------------------------------------------

def MaximumSkew(Genome):
    
    '''Find position in genome where skew diagram attains a maximum (locate ter).
    Skew decreases along the reverse half-strand and increases along the forward half-strand. 
    Skew should achieve a maximum where the forward half-strand ends and the reverse 
    half-strand begins (ter).
    
    Input:  A DNA string Genome
    Output: A list containing all integers i maximizing Skew(Prefix_i(Text)) 
    over all values of i (from 0 to length(Genome)).
    '''
    
    positions = [] 
    skew = SkewArray(Genome)
    max_skew = max(skew)
    for i in range(len(skew)):
        if skew[i] == max_skew:
            positions.append(i)
    return positions


def SkewArray(Genome):
    
    '''Generate the difference (skew) array for difference in bases.
    The difference between G and C is increasing for the forward half-strand 
    and decreasing for the reverse half-strand.
    
    Input:  A String Genome
    Output: The skew array of Genome as a list.'''
    
    skew = {}
    skew = [0]
    count_G = 0
    count_C = 0
    for i in range(len(Genome)):
        if Genome[i] == 'G':
            count_G += 1
        elif Genome[i] == 'C':
            count_C += 1
        skew.append(count_G - count_C)
    return skew


# Example
# ------------
Genome = "GATACACTTCCCGAGTAGGTACTG"
print(MaximumSkew(Genome))
# Ouput: [1, 2, 3, 4]




# -----------------------------------------------------
# Compute the number of mismatches (Hamming distance) between two strings
# -----------------------------------------------------

def HammingDistance(p, q):
    '''Compute the number of mismatches (Hamming distance) between two strings.
    
    Input:  Two strings p and q.
    Output: An integer value representing the Hamming Distance between p and q.'''
    
    if len(p) != len(q):
        raise ValueError("Input strings must have the same length")
    distance=0

    # Iterate through the characters of both strings and compare them
    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance


# Examples
# -----------------
HammingDistance(p="GGGCCGTTGGT", q="GGACCGTTGAC")
# Output: 3


p ="TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC"
q = "GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA"
HammingDistance(p, q)
# Output: 50


p ="CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT"
q = "CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG"
HammingDistance(p, q)
# Output: 36




with open("dataset_30278_3.txt", "r") as file:
    p = file.readline().strip()
    q = file.readline().strip()
HammingDistance(p, q)
# Output: 895




# -----------------------------------------------------
# Find all approximate occurrences of a pattern in a string
# -----------------------------------------------------



def ApproximatePatternMatching(Text, Pattern, d):
    '''Find all approximate occurrences of a pattern in a string.
    
    Input:  Strings Pattern and Text along with an integer d
    Output: A list containing all starting positions where Pattern appears
    as a substring of Text with at most d mismatches.'''
    
    positions = []
    n = len(Text)
    k = len(Pattern)
    for i in range(n-k+1):
        if HammingDistance(Text[i:i+k], Pattern) <= d:
            positions.append(i)
    return positions
    
    
def HammingDistance(p, q):

    if len(p) != len(q):
        raise ValueError("Input strings must have the same length")
    distance=0

    # Iterate through the characters of both strings and compare them
    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance




# Examples
# ---------------
Pattern = "ATTCTGGA"
Text = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
d=2
print(*ApproximatePatternMatching(Text, Pattern, d))
# Ouput: 6 7 26 27


with open("dataset_30278_4.txt", "r") as file:
    Pattern = file.readline().strip()
    Text = file.readline().strip()
    d = int(file.readline().strip())
print(*ApproximatePatternMatching(Text, Pattern, d))
# Output: 3 22 46 49 74 82 90 92 98 127 149 164 236 etc





# -----------------------------------------------------
# Compute number of occurrences of Pattern in Text with at most d mismatches
# -----------------------------------------------------


def ApproximatePatternCount(Pattern, Text, d):
    '''Compute the number of occurrences of Pattern in Text with at most d mismatches.
    
    Input:  Strings Pattern and Text, and an integer d
    Output: The number of times Pattern appears in Text with at most d mismatches
    '''

    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Pattern, Text[i:i+len(Pattern)]) <= d:
            count+=1
    return count

def HammingDistance(p, q):
    if len(p) != len(q):
        raise ValueError("Input strings must have the same length")
    distance=0

    # Iterate through the characters of both strings and compare them
    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance



# Examples
# -------------
ApproximatePatternCount(Text="TTTAGAGCCTTCAGAGG", Pattern="GAGG", d=2)
# Output: 4

ApproximatePatternCount(Text="AACAAGCTGATAAACATTTAAAGAG", Pattern="AAAAA", d=2)
# Output: 11

# Compute Count2(CATGCCATTCGCATTGTCCCAGTGA, CCC).
ApproximatePatternCount(Text="CATGCCATTCGCATTGTCCCAGTGA", Pattern="CCC", d=2)
# Output: 15


with open("dataset_30278_6.txt", "r") as file:
    Pattern = file.readline().strip()
    Text = file.readline().strip()
    d = int(file.readline().strip())
print(len(ApproximatePatternMatching(Text, Pattern, d)))
# Output: 60




# -----------------------------------------------------
# Immediate Neighbors
# -----------------------------------------------------

def Neighbors(Pattern, d):
    '''Find the set of all k-mers whose Hamming distance from Pattern does not exceed d.
    '''
    
    if d == 0:
        return {Pattern}
    if len(Pattern) == 0:
        return {""}
    nucleotides = ['A', 'C', 'G', 'T']
    neighborhood = set()

    # Recursive call for suffix
    suffix_neighbors = Neighbors(Pattern[1:], d)

    for neighbor in suffix_neighbors:
        # Case 1: No change to the first symbol
        if HammingDistance(Pattern[1:], neighbor) < d:
            for nucleotide in nucleotides:
                neighborhood.add(nucleotide + neighbor)
                
        # Case 2: Change the first symbol
        else:
            neighborhood.add(Pattern[0] + neighbor)
            
    return neighborhood

 
def HammingDistance(p, q):
    if len(p) != len(q):
        raise ValueError("Input strings must have the same length")
    distance=0

    # Iterate through the characters of both strings and compare them
    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance




# Examples
# -------------
print(*Neighbors(Pattern = "GTAGTTTGG", d = 2))
# Output: GCAGTTTAG GTAGTGTCG GTAGTCTGC GTAGGCTGG etc

# How many 5-mers are in the 2-neighborhood of Pattern = TGCAT?
print(len(Neighbors(Pattern = "TGCAT", d = 2)))
# Output: 106


with open("dataset_30282_4.txt", "r") as file:
    Pattern = file.readline().strip()
    d = int(file.readline().strip())
print(*Neighbors(Pattern, d))
# Output: GCATGAGCGA GGCTGAGCAA TGATGAACGA AAATGAGCGA etc




# -----------------------------------------------------
# Frequent Words with Mismatches
# -----------------------------------------------------

def HammingDistance(p, q):
    if len(p) != len(q):
        raise ValueError("Input strings must have the same length")
    distance=0

    # Iterate through the characters of both strings and compare them
    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance


def Neighbors(Pattern, d):
    if d == 0:
        return {Pattern}
    if len(Pattern) == 0:
        return {""}
    nucleotides = ['A', 'C', 'G', 'T']
    neighborhood = set()

    # Recursive call for suffix
    suffix_neighbors = Neighbors(Pattern[1:], d)

    for neighbor in suffix_neighbors:
        # Case 1: No change to the first symbol
        if HammingDistance(Pattern[1:], neighbor) < d:
            for nucleotide in nucleotides:
                neighborhood.add(nucleotide + neighbor)
                
        # Case 2: Change the first symbol
        else:
            neighborhood.add(Pattern[0] + neighbor)
            
    return neighborhood

def FrequentWordsMismatches(Text, k, d):
    '''Find frequent words in Text of length k with number of mismatches at most d.'''
    patterns = {}
    for i in range(len(Text) - k + 1):
        pattern = Text[i:i+k]
        neighborhood = Neighbors(pattern, d)
        for neighbor in neighborhood:
            patterns[neighbor] = patterns.get(neighbor, 0) + 1
    max_count = max(patterns.values())
    return [pattern for pattern, count in patterns.items() if count == max_count]



# Example
# --------------
with open("dataset_30278_9.txt", "r") as file:
    Text = file.readline().strip()
    line2 = file.readline().strip()
    k, d = map(int, line2.split())
print(*FrequentWordsMismatches(Text, k, d))
# Output: AAAAC




# -----------------------------------------------------
# Frequent Words with Mismatches and Reverse Complements
# -----------------------------------------------------

def HammingDistance(p, q):
    if len(p) != len(q):
        raise ValueError("Input strings must have the same length")
    distance=0

    # Iterate through the characters of both strings and compare them
    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance

def Neighbors(Pattern, d):
    if d == 0:
        return {Pattern}
    if len(Pattern) == 0:
        return {""}
    nucleotides = ['A', 'C', 'G', 'T']
    neighborhood = set()

    # Recursive call for suffix
    suffix_neighbors = Neighbors(Pattern[1:], d)

    for neighbor in suffix_neighbors:
        # Case 1: No change to the first symbol
        if HammingDistance(Pattern[1:], neighbor) < d:
            for nucleotide in nucleotides:
                neighborhood.add(nucleotide + neighbor)
                
        # Case 2: Change the first symbol
        else:
            neighborhood.add(Pattern[0] + neighbor)
            
    return neighborhood

def ReverseComplement(Pattern):
    Pattern = Reverse(Pattern) # reverse all letters in a string
    Pattern = Complement(Pattern) # complement each letter in a string
    return Pattern

def Reverse(Pattern):
    # your code here
    z=""
    for char in Pattern:
        z = char + z
    return z

def Complement(Pattern):
    string=''
    dict={'A':'T','T':'A','G':'C','C':'G'}
    for char in Pattern:
        string = string + dict[char]
    return string
    

def FrequentWordsMismatchesReverseComplements(Text, k, d):
    '''Find frequent words in Text of length k with number of mismatches at most d.
    Include reverse complements in the search.'''
    
    freq = defaultdict(int)
    for i in range(len(Text) - k + 1):
        pattern = Text[i:i+k]
        neighborhood = neighbors(pattern, d)
        for neighbor in neighborhood:
            freq[neighbor] += 1
            freq[ReverseComplement(neighbor)] += 1  # Count reverse complement too
    max_freq = max(freq.values())
    return [pattern for pattern, count in freq.items() if count == max_freq]


# Example
# ---------------
with open("dataset_30278_10.txt", "r") as file:
    Text = file.readline().strip()
    line2 = file.readline().strip()
    k, d = map(int, line2.split())
print(*FrequentWordsMismatchesReverseComplements(Text, k, d))
# Output: CCCCCCC GGGGGGG




# -----------------------------------------------------
# Find a DnaA box in S. enterica
# -----------------------------------------------------

def MinimumSkew(Genome):
    positions = [] 
    skew = SkewArray(Genome)
    min_skew = min(skew)
    for i in range(len(skew)):
        if skew[i] == min_skew:
            positions.append(i)
    return positions

# Input:  A String Genome
# Output: SkewArray(Genome)
def SkewArray(Genome):
    skew = {}
    skew = [0]
    count_G = 0
    count_C = 0
    for i in range(len(Genome)):
        if Genome[i] == 'G':
            count_G += 1
        elif Genome[i] == 'C':
            count_C += 1
        skew.append(count_G - count_C)
    return skew

def HammingDistance(p, q):
    if len(p) != len(q):
        raise ValueError("Input strings must have the same length")
    distance=0

    # Iterate through the characters of both strings and compare them
    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance

def Neighbors(Pattern, d):
    if d == 0:
        return {Pattern}
    if len(Pattern) == 0:
        return {""}
    nucleotides = ['A', 'C', 'G', 'T']
    neighborhood = set()

    # Recursive call for suffix
    suffix_neighbors = Neighbors(Pattern[1:], d)

    for neighbor in suffix_neighbors:
        # Case 1: No change to the first symbol
        if HammingDistance(Pattern[1:], neighbor) < d:
            for nucleotide in nucleotides:
                neighborhood.add(nucleotide + neighbor)
                
        # Case 2: Change the first symbol
        else:
            neighborhood.add(Pattern[0] + neighbor)
            
    return neighborhood

def ReverseComplement(Pattern):
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return Pattern

def Reverse(Pattern):
    z=""
    for char in Pattern:
        z = char + z
    return z

def Complement(Pattern):
    string=''
    dict={'A':'T','T':'A','G':'C','C':'G'}
    for char in Pattern:
        string = string + dict[char]
    return string

from collections import defaultdict

def FrequentWordsMismatchesReverseComplements(Text, k, d):
    freq = defaultdict(int)
    for i in range(len(Text) - k + 1):
        pattern = Text[i:i+k]
        neighborhood = Neighbors(pattern, d)
        for neighbor in neighborhood:
            freq[neighbor] += 1
            freq[ReverseComplement(neighbor)] += 1  # Count reverse complement too
    max_freq = max(freq.values())
    return [pattern for pattern, count in freq.items() if count == max_freq]




with open("Salmonella_enterica.txt", "r") as file:
    lines = file.readlines()
    Info = lines[0].strip()        
    Genome = ''.join(line.strip() for line in lines[1:])
    
print(MinimumSkew(Genome))
# 3764856 3764858

Text = Genome[3764857-250: 3764857+250]

print(*FrequentWordsMismatchesReverseComplements(Text, k=9, d=1))
# Output: CCGGAAGCT AGCTTCCGG

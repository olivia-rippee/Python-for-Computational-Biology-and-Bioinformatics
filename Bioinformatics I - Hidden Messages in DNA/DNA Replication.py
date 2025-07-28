import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics I - Hidden Messages in DNA/Data")


# -----------------------------------------------------
# Transcribe DNA to RNA
# -----------------------------------------------------

#file = open('sample_dna.txt', 'r')
#DNA = file.read()
#print("DNA String: ", DNA)


from Bio.Seq import Seq

DNA_sequence =  Seq("ACGTACGTCCCAGCC")
RNA_sequence = DNA_sequence.transcribe()
print(RNA_sequence)





# -----------------------------------------------------
# Transcribe + Translate DNA to Protein
# -----------------------------------------------------

from Bio.Seq import Seq

DNA_sequence =  Seq("ACGTACGTCCCAGCC")
protein_sequence = DNA_sequence.translate()
print(protein_sequence)





# -----------------------------------------------------
# Translate RNA to Protein
# -----------------------------------------------------

from Bio.Seq import Seq
rna_sequence = Seq("CCAAGUACAGAGAUUAAC")
rna_sequence = Seq("CCACGUACUGAAAUUAAC")
rna_sequence = Seq("CCCAGGACUGAGAUCAAU")
rna_sequence = Seq("CCCAGUACCGAGAUGAAU")
protein_sequence = rna_sequence.translate()
print(protein_sequence)




# -----------------------------------------------------
# Count occurrences of a pattern in a string
# -----------------------------------------------------

def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count


# Example: oriC of Vibrio cholerae
# ---------------------------------
Text = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"
Pattern = "TGATCA"
print(PatternCount(Text, Pattern))
# Output: 8


# Example: oriC region from T petrophila
# ---------------------------------------
Text = "AACTCTATACCTCCTTTTTGTCGAATTTGTGTGATTTATAGAGAAAATCTTATTAACTGAAACTAAAATGGTAGGTTTGGTGGTAGGTTTTGTGTACATTTTGTAGTATCTGATTTTTAATTACATACCGTATATTGTATTAAATTGACGAACAATTGCATGGAATTGAATATATGCAAAACAAACCTACCACCAAACTCTGTATTGACCATTTTAGGACAACTTCAGGGTGGTAGGTTTCTGAAGCTCTCATCAATAGACTATTTTAGTCTTTACAAACAATATTACCGTTCAGATTCAAGATTCTACAACGCTGTTTTAATGGGCGTTGCAGAAAACTTACCACCTAAAATCCAGTATCCAAGCCGATTTCAGAGAAACCTACCACTTACCTACCACTTACCTACCACCCGGGTGGTAAGTTGCAGACATTATTAAAAACCTCATCAGAAGCTTGTTCAAAAATTTCAATACTCGAAACCTACCACCTGCGTCCCCTATTATTTACTACTACTAATAATAGCAGTATAATTGATCTGA"
count_1 = PatternCount("ATGATCAAG", Text)
count_2 = PatternCount("CTTGATCAT", Text)
print(count_1 + count_2)
# Output: 0


# Examples
# ---------------------------------------
print(PatternCount(Text="ACTGTACGATGATGTGTGTCAAAG", Pattern="TGT"))
# Output: 4

print(PatternCount(Text="GCGCG", Pattern="GCG"))
# Output: 2

with open("dataset_30272_6.txt", "r") as file:
    Text = file.readline().strip()
    Pattern = file.readline().strip()
print(Text)
print(Pattern)
print(PatternCount(Text, Pattern))
# Output: 35






# -----------------------------------------------------
# Find Frequent Patterns (of length k) Function
# -----------------------------------------------------

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n - k + 1):
        kmer = Text[i:i + k]
        if kmer in freq:
            freq[kmer] += 1
        else:
            freq[kmer] = 1
    return freq

def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words



# Example: oriC of Vibrio cholerae
# -----------------------------------
Text = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"
k = 10
print(FrequentWords(Text, k))


# Examples
# --------------
print(*FrequentWords(Text="TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT", k=3))
# Output: GTG


with open("dataset_30272_13.txt", "r") as file:
    Text = file.readline().strip()
    k = int(file.readline().strip())
print(Text)
print(k)
print(FrequentWords(Text, k))
# Output: ['AGCCCAACCGG']





# -----------------------------------------------------
# Find the reverse complement of a sequence
# -----------------------------------------------------
## Input:  A DNA string Pattern
## Output: The reverse complement of Pattern
    
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
    


# Examples
# -------------
print(ReverseComplement("GCTAGCT"))
# Output: AGCTAGC


print(ReverseComplement("TTGTGTC"))
# Output: GACACAA


with open("dataset_30273_2.txt", "r") as file:
    Pattern = file.readline().strip()
print(Pattern)
print(len(Pattern))
print(ReverseComplement(Pattern))
# Output: 9385 base DNA sequence





    
# -----------------------------------------------------    
# Locate all occurrences of a pattern in a string
# -----------------------------------------------------

def PatternMatching(Pattern, Genome):
    positions = []  # Output variable
    Pattern_length = len(Pattern)
    Genome_length = len(Genome)
    for i in range(Genome_length - Pattern_length + 1):
        if Genome[i:i + Pattern_length] == Pattern:
            positions.append(i)
    return positions



# Examples
# ---------------------
Pattern = 'ATAT'
Genome = 'GATATATGCATATACTT'
print(PatternMatching(Pattern, Genome))
# Output: [1, 3, 9]


Pattern = "CTTGATCAT"
Genome = "CTTGATCATACGTCCTTGATCATAGTCAGTCGCCTTGATCATAAAGGTTC"
positions = PatternMatching(Pattern, Genome)
print(positions)
# Output: [0, 14, 33]


Text = "GACGATATACGACGATA"
Pattern = "ATA"
print(*PatternMatching(Pattern, Text))
# Output: 4 6 14



with open("dataset_30273_5.txt", "r") as file:
    Pattern = file.readline().strip()
    Genome = file.readline().strip()
print(*PatternMatching(Pattern, Genome))
# Output: 1 16 36 139 171 202 422 456 473 497 522 529 536 554 561 604 611 etc.
# * prints as list with no [], (), or commas



with open("Vibrio_cholerae.txt", "r") as file:
    Genome = file.readline().strip()
    
Pattern = "CTTGATCAT"
print(*PatternMatching(Pattern, Genome))
# Output: 60039 98409 129189 152283 152354 152411 163207 197028 200160 357976 376771 392723 532935 600085 622755 1065555

Pattern = "ATGATCAAG"
print(*PatternMatching(Pattern, Genome))
# Output: 116556 149355 151913 152013 152394 186189 194276 200076 224527 307692 479770 610980 653338 679985 768828 878903 985368






# -----------------------------------------------------    
# Find every k-mer that forms a clump
# -----------------------------------------------------

# A k-mer forms an (L, t)-clump inside a (longer) string Genome if there is a 
# substring of Genome of length L in which this k-mer appears at least t times.
# One window of length L contains t instances of the sequence


def PatternCount(Text,Pattern):
    count=0
    len_pattern=len(Pattern)
    for i in range(len(Text)-len(Pattern)+1):
        if Pattern == Text[i:i+len_pattern]:
            count+=1
    return count


# Example: "TGCA" forms a (25, 3)-clump in the Genome:
# "gatcagcataagggtccCTGCAATGCATGACAAGCCTGCAGTtgttttac"
# One window of 25 consecutive bases contains 3 instances of the sequence TGCA
# ---------------------------------
sequence="tgca"
genome="gatcagcataagggtccctgcaatgcatgacaagcctgcagttgttttac"
window_length=25

for i in range(len(genome)-window_length+1):
    window=genome[i:i+window_length]
    count_seq_window=PatternCount(window, sequence)
    print(f"The sequence {sequence} is observed {count_seq_window} time(s) in window:{window}")
    

# Example 
# --------------------
PatternCount("CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC", "CGCG")
# Output: 5




# -----------------------------------------------------    
# Find every k-mer that forms a clump
# -----------------------------------------------------

def FindClumps(Text, k, L, t):
     Patterns = []
     n = len(Text)
     for i in range(n - L + 1):
          window = Text[i:i + L]
          freqMap = FrequencyMap(window, K)
          for x in freqMap:
               if freqMap[x] >= t and x not in Patterns:
                    Patterns.append(x)
     return Patterns

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n - k + 1):
        kmer = Text[i:i + k]
        if kmer in freq:
            freq[kmer] += 1
        else:
            freq[kmer] = 1
    return freq



# Examples
# --------------
with open("dataset_30274_5.txt", "r") as file:
    Genome = file.readline().strip()
    line2 = file.readline().strip()
    k, L, t = map(int, line2.split())
print(*FindClumps(Genome, k, L, t))
# Output: GTCTATTTCC TTTTCCGTTG CGTTAAAGTC TCGTTAAAGT TTTGGCGCCC TCGTACGACT ACTCGTACGA CTCGTACGAC TCCTGCCATT CCTGCCATTC GCGCGCTTAT CGCGCTTATG





# -----------------------------------------------------    
# Count the number of occurrences of a base in a window of ExtendedGenome
# -----------------------------------------------------
# (first x bases added to the end because bacterial DNA is circular)

# Input:  Strings Genome and symbol
# Output: SymbolArray(Genome, symbol)

def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

# Input:  Strings Text and Pattern
# Output: The number of times Pattern appears in Text
def PatternCount(symbol, ExtendedGenome):
    count = 0
    for i in range(len(ExtendedGenome)):
        if ExtendedGenome[i:i+(len(symbol))] == symbol:
            count += 1
    return count


# Example
# -------------------
ExtendedGenome = 'AAAAGGGG'
symbol = 'A'
print(PatternCount(symbol, ExtendedGenome))




# -----------------------------------------------------    
# Difference (skew) array for difference in bases
# -----------------------------------------------------
# Diff btwn G and C:
# increasing = forward half-strand; decreasing = reverse half-strand

# Input:  A String Genome
# Output: The skew array of Genome as a list.

def SkewArray(Genome):
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
# Skew decreases along the reverse half-strand and increases along the forward half-strand. 
# Skew should achieve a minimum where the reverse half-strand ends and the forward half-strand begins (ori)

# Input:  A DNA string Genome
# Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)

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
    positions = [] 
    skew = SkewArray(Genome)
    max_skew = max(skew)
    for i in range(len(skew)):
        if skew[i] == max_skew:
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


# Example
# ------------
Genome = "GATACACTTCCCGAGTAGGTACTG"
print(MaximumSkew(Genome))
# Ouput: [1, 2, 3, 4]




# -----------------------------------------------------
# Compute the number of mismatches (Hamming distance) between two strings
# -----------------------------------------------------

# Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q.

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

# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches

def ApproximatePatternMatching(Text, Pattern, d):
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

# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches

def ApproximatePatternCount(Pattern, Text, d):
    count = 0 # initialize count variable
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
# Find the set of all k-mers whose Hamming distance from Pattern does not exceed d


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
    patterns = {}
    for i in range(len(Text) - k + 1):
        pattern = Text[i:i+k]
        neighborhood = neighbors(pattern, d)
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
# counts the same k-mers twice if they are reverse palindromes


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





# -----------------------------------------------------
# Quiz Questions
# -----------------------------------------------------

8.7 // 4
# Ouput: 2.0



x=0
for y in range(0,5):
    x+=y
print(x)
# Ouput: 10 (doesn't include 5)



a=list(range(5))
b=a
a[2]=12
print(b)
# Output: b = [0,1,12,3,4]


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
    '''Count occurrences of a pattern in a string.'''
    
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
    '''Calculate the frequency of each k-mer in Text.
    
    Input: Text and an integer k.
    Output: A dictionary with unique k-mers as keys and their frequencies in counts.'''
    
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
    '''Find the most frequent patterns (of length k) in Text.
    
    Input: Text and an integer k.
    Output: A list of most frequent k-mers in Text.'''
    
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
print(*FrequentWords(Text, k))
# Output: CTCTTGATCA TCTTGATCAT


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
    
def ReverseComplement(Pattern):
    '''Generate the reverse complemement of a sequence.
    
    Input:  A DNA string Pattern
    Output: The reverse complement of Pattern'''
    
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return Pattern

def Reverse(Pattern):
    '''Reverse all letters in a sequence.'''
    
    z=""
    for char in Pattern:
        z = char + z
    return z

def Complement(Pattern):
    '''Complement each letter in a string.'''
    
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
# Output: a 9385 base DNA sequence





    
# -----------------------------------------------------    
# Locate all occurrences of a pattern in a string
# -----------------------------------------------------

def PatternMatching(Pattern, Genome):
    '''Locate all occurrences of a pattern in a string.'''
    
    positions = []
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

def PatternCount(Text,Pattern):
    '''Find every k-mer that forms a clump.
    
    A k-mer forms an (L, t)-clump inside a (longer) string Genome if there is a 
    substring of Genome of length L in which this k-mer appears at least t times.
    One window of length L contains t instances of the sequence.
    
    Input:  Strings Text and Pattern
    Output: The number of times Pattern appears in Text
    '''
    
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
    '''Find every k-mer that forms a clump.'''
    
    Patterns = []
    n = len(Text)
    for i in range(n - L + 1):
         window = Text[i:i + L]
         freqMap = FrequencyMap(window, k)
         for x in freqMap:
             if freqMap[x] >= t and x not in Patterns:
                 Patterns.append(x)
    return Patterns


def FrequencyMap(Text, k):
    '''Calculate the frequency of each k-mer in Text.
    
    Input: Text and an integer k.
    Output: A dictionary with unique k-mers as keys and their frequencies in counts.'''
    
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

def SymbolArray(Genome, symbol):
    '''Count the number of occurrences of a base in a window of the genome.
    The first x bases are added to the end because bacterial DNA is circular.
    
    Input:  Strings Genome and symbol.
    Output: SymbolArray(Genome, symbol).'''
    
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


def PatternCount(symbol, ExtendedGenome):
    '''Find every k-mer that forms a clump.
    
    A k-mer forms an (L, t)-clump inside a (longer) string Genome if there is a 
    substring of Genome of length L in which this k-mer appears at least t times.
    One window of length L contains t instances of the sequence.
    
    Input:  Strings Text and Pattern
    Output: The number of times Pattern appears in Text
    '''
    
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
# Basic Python Quiz Questions
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


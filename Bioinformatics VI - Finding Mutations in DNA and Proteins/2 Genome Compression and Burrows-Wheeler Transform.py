import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics VI - Finding Mutations in DNA and Proteins/Data")

# -----------------------------------------------
# Suffix Arrays
# -----------------------------------------------

def SuffixArrayConstruction(text):
    '''Construct the suffix array of a string (memory-efficient alternative to suffix trees).

    Input: A string Text.
    Output: SuffixArray(Text), as a space-separated collection of integers.'''
    
    # Create a list of (suffix, index) pairs
    suffixes = [(text[i:], i) for i in range(len(text))]
    
    # Sort the suffixes lexicographically based on the suffix
    sorted_suffixes = sorted(suffixes)
    
    # Extract and return just the indices from the sorted suffixes
    return [index for (suffix, index) in sorted_suffixes]



# Example 1
# -----------
text = "AACGATAGCGGTAGA$"
suffix_array = SuffixArrayConstruction(text)
print(" ".join(map(str, suffix_array)))
# Output: 15 14 0 1 12 6 4 2 8 13 3 7 9 10 11 5


# Example 2
# -----------
with open("dataset_30231_2.txt", "r") as file:
    text = file.read().strip()
    
suffix_array = SuffixArrayConstruction(text)
with open("output.txt", "w") as file:
    file.write(" ".join(map(str, suffix_array)))
# Output: 966 214 650 190 885 215 309 914 844 651 191 431 568 ... 144 326 143



# -----------------------------------------------
# Construct the Burrows-Wheeler Transform
# -----------------------------------------------
    
def ConstructBWT(Text):
    '''Construct the Burrows-Wheeler transform of a string.

    BWT: Form all possible cyclic rotations of Text, then order them lexicographically 
    to form a |Text| × |Text| matrix of symbols denoted M(Text). The last column of 
    M(Text) is called the Burrows-Wheeler transform of Text. This is the only column 
    of M(Text) that can reconstruct Text from BWT(Text).
    
    Input: A string Text.
    Output: BWT(Text).'''
    
    # Generate all cyclic rotations of the text
    rotations = [Text[i:] + Text[:i] for i in range(len(Text))]
    
    # Sort the rotations lexicographically
    sorted_rotations = sorted(rotations)
    
    # Construct the BWT by taking the last character of each sorted rotation
    bwt = ''.join(rotation[-1] for rotation in sorted_rotations)
    
    return bwt



# Example 1
# -----------
text = "GCGTGCCTGGTCA$"
bwt = ConstructBWT(text)
print(bwt) # Output: ACTGGCT$TGCGGC



# Example 2
# -----------
with open("dataset_30223_5.txt", "r") as file:
    text =  file.read().strip()
bwt = ConstructBWT(text)
print(bwt)
with open("output.txt", "w") as file:
    file.write(bwt)
# Output: GGGGCAAACAT...GCCTTTTCAATC$AGAGCGA...AGGCTGTGA



# -----------------------------------------------
# Inverse Burrows-Wheeler
# -----------------------------------------------

def InverseBWT(Transform):
    '''Reconstruct a string from its Burrows-Wheeler transform.
    
    Uses the First-Last property: The k-th occurrence of a symbol in 
    FirstColumn and the k-th occurrence of this symbol in LastColumn 
    correspond to the same position of this symbol in Text.

    Input: A string Transform (with a single "$$" symbol).
    Output: The string Text such that BWT(Text) = Transform.'''
    
    # Last column is the BWT transform
    LastColumn = list(Transform)
    
    # First column is the sorted BWT transform
    FirstColumn = sorted(LastColumn)
    
    # Map character occurrences to their positions in Last and First columns
    def CountOccurrences(column):
        counts = {}
        result = []
        for char in column:
            if char not in counts:
                counts[char] = 0
            result.append((char, counts[char]))
            counts[char] += 1
        return result

    LastColumnIndexed = CountOccurrences(LastColumn)
    FirstColumnIndexed = CountOccurrences(FirstColumn)

    # Build a mapping from LastColumn to FirstColumn
    LastToFirst = {}
    for i, pair in enumerate(LastColumnIndexed):
        LastToFirst[i] = FirstColumnIndexed.index(pair)

    # Reconstruct the text using the LastToFirst mapping
    index = LastColumn.index('$')
    original = []

    for _ in range(len(Transform)):
        char = LastColumn[index]
        original.append(char)
        index = LastToFirst[index]

    return ''.join(reversed(original))


# Example 1
# -----------
bwt = "TTCCTAACG$A"
text = InverseBWT(bwt)
print(text) # Output: TACATCACGT$


# Example 3
# -----------
with open("dataset_30225_10.txt", "r") as file:
    bwt = file.read().strip()

text = InverseBWT(bwt)
with open("output.txt", "w") as file:
    file.write(text)
# Output: TTTTCCAGTCTCTGGGGTATGCCAAT...AACAAAAGTGAA$



# -----------------------------------------------
# Pattern Matching with Burrows-Wheeler
# -----------------------------------------------

def BWMatching(LastColumn, Pattern, LastToFirst):
    '''Counts the total number of matches of Pattern in Text.
    
    Input: A string BWT(Text), followed by a space-separated collection of Patterns.
    Output: A space-separated list of integers, where the i-th integer corresponds to the 
    number of substring matches of the i-th member of Patterns in Text.'''
    
    top = 0
    bottom = len(LastColumn) - 1

    while top <= bottom:
        if Pattern:
            symbol = Pattern[-1]
            Pattern = Pattern[:-1]

            # Get range in LastColumn from top to bottom
            subColumn = LastColumn[top:bottom+1]
            if symbol in subColumn:
                # Find first and last index of symbol in subColumn
                topIndex = subColumn.index(symbol) + top
                bottomIndex = bottom - subColumn[::-1].index(symbol)
                top = LastToFirst[topIndex]
                bottom = LastToFirst[bottomIndex]
            else:
                return 0
        else:
            return bottom - top + 1
    return 0

def BuildLastToFirst(LastColumn):
    # Count ranks of characters in LastColumn and FirstColumn
    from collections import defaultdict

    def CountRanks(column):
        counts = defaultdict(int)
        ranks = []
        for char in column:
            ranks.append((char, counts[char]))
            counts[char] += 1
        return ranks

    FirstColumn = sorted(LastColumn)
    lastRanks = CountRanks(LastColumn)
    firstRanks = CountRanks(FirstColumn)

    # Map each (char, rank) in LastColumn to its index in FirstColumn
    lastToFirst = {}
    for index, pair in enumerate(firstRanks):
        lastToFirst[pair] = index

    return [lastToFirst[rank] for rank in lastRanks]

def BWMatchingDriver(BWTAndPatterns):
    parts = BWTAndPatterns.strip().split()
    BWT = parts[0]
    Patterns = parts[1:]

    LastColumn = list(BWT)
    LastToFirst = BuildLastToFirst(LastColumn)

    results = []
    for pattern in Patterns:
        count = BWMatching(LastColumn, pattern, LastToFirst)
        results.append(str(count))
    return ' '.join(results)


# Example 1
# -----------
input_data = "TCCTCTATGAGATCCTATTCTATGAAACCTTCA$GACCAAAATTCTCCGGC CCT CAC GAG CAG ATC"
print(BWMatchingDriver(input_data))
# Output: 2 1 1 0 1


# Example 2
# -----------
with open("dataset_30226_8.txt", "r") as file:
    input_data = file.read().strip()

matches = BWMatchingDriver(input_data)
with open("output.txt", "w") as file:
    file.write(matches + '\n')

# Output: 1 1 1 0 0 0 1 0 1 1 1 0 1 1 0 0 1 0 1 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 ...


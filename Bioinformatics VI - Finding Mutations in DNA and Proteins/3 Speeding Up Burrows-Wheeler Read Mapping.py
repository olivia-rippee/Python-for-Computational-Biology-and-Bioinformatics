import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics VI - Finding Mutations in DNA and Proteins/Data")

# -----------------------------------------------
# Faster Burrows-Wheeler Matching
# -----------------------------------------------

from collections import defaultdict

def BetterBWMatching(bwt, patterns):
    '''Counts the total number of matches of Pattern in Text (but does not say where). This version 
    is faster than the previous version as it eliminates the need for both FirstColumn and LastToFirst.
    
    Input: A string BWT(Text) followed by a space-separated collection of strings Patterns.
    Output: A space-separated list of integers, where the i-th integer corresponds to the 
    number of substring matches of the i-th member of Patterns in Text.'''

    first_occurrence, counts = PreprocessBWT(bwt)
    result = []
    for pattern in patterns:
        count = CountOccurrences(pattern, bwt, first_occurrence, counts)
        result.append(str(count))
    return ' '.join(result)

def PreprocessBWT(bwt):
    # Build FirstOccurrence
    first_col = ''.join(sorted(bwt))
    first_occurrence = {}
    for i, char in enumerate(first_col):
        if char not in first_occurrence:
            first_occurrence[char] = i

    # Build Count dictionary
    counts = defaultdict(lambda: [0] * (len(bwt) + 1))
    for i in range(1, len(bwt) + 1):
        char = bwt[i - 1]
        for c in counts:
            counts[c][i] = counts[c][i - 1]
        counts[char][i] += 1

    return first_occurrence, counts

def CountOccurrences(pattern, bwt, first_occurrence, counts):
    top = 0
    bottom = len(bwt) - 1
    while top <= bottom:
        if pattern:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            if counts[symbol][bottom + 1] - counts[symbol][top] > 0:
                top = first_occurrence[symbol] + counts[symbol][top]
                bottom = first_occurrence[symbol] + counts[symbol][bottom + 1] - 1
            else:
                return 0
        else:
            return bottom - top + 1
    return 0



# Example 1
# -----------
bwt = "GGCGCCGC$TAGTCACACACGCCGTA"
patterns_str = "ACC CCG CAG"
patterns = patterns_str.split()
matches = BetterBWMatching(bwt, patterns)
print(matches) # Output: 1 2 1


# Example 2
# -----------
with open("dataset_30227_7.txt", "r") as file:
    lines = file.read().strip().split('\n')
    bwt = lines[0]
    patterns = lines[1].split()
    
matches = BetterBWMatching(bwt, patterns)
with open("output.txt", "w") as file:
    file.write(matches)
    
# Output: 0 1 0 1 0 1 0 0 1 1 0 1 0 0 0 0 0 0 0 1 0 0 1 1 1 1 0 0 0 1 0 1 0 ...


# -----------------------------------------------
# Partial Suffix Arrays
# -----------------------------------------------

def SuffixArrayK(Text, k):
    '''Construct a partial suffix array to save memory by retaining only the 
    elements of full array that are divisible by k, along with their indices i.

    Input: A string Text and a positive integer k.
    Output: SuffixArrayk(Text), in the form of a list of ordered pairs 
    (i, SuffixArray(i)) for all nonempty entries in the partial suffix array.'''

    
    # Step 1: Build full suffix array
    suffixes = [(Text[i:], i) for i in range(len(Text))]
    suffixes.sort()
    full_suffix_array = [index for (suffix, index) in suffixes]
    
    # Step 2: Create partial suffix array for indices divisible by k
    partial_suffix_array = []
    for i, suffix_index in enumerate(full_suffix_array):
        if suffix_index % k == 0:
            partial_suffix_array.append((i, suffix_index))
    
    return partial_suffix_array


# Example 1
# -----------
Text = "PANAMABANANAS$"
k = 5
result = SuffixArrayK(Text, k)
for i, val in result:
    print(i, val)

# Output:
    # 1 5
    # 11 10
    # 12 0
    

# Example 2
# -----------
with open("dataset_30234_2.txt", "r") as file:
    Text = file.readline().strip()
    k = int(file.readline().strip())
    
partial_suffix_array = SuffixArrayK(Text, k)
with open("output.txt", "w") as file:
    for i, val in partial_suffix_array:
            file.write(f"{i} {val}\n")

# Output:
    # 0 20000   1 55   13 5360   14 12695   17 2465 ... 19991 10220    19996 8765


# -----------------------------------------------
# Multiple Pattern Matching
# -----------------------------------------------

def MultiplePatternMatching(text, patterns):
    '''Identify all occurrences of multiple pattern strings within a given text. It is designed to 
    efficiently search for many short patterns (e.g., DNA reads) within a large string (e.g., genome), 
    using a memory-efficient version of the BetterBWMatching algorithm with checkpoint arrays.
    
    Input: A string Text followed by a space-separated collection of strings Patterns.
    Output: For each string Pattern in Patterns, the string Pattern followed by a colon, followed 
    by a space-separated list of all starting positions in Text where Pattern appears as a substring.'''
    
    bwt, suffix_array = BuildBWT(text)
    C = 100  # checkpoint interval
    first_occurrence = BuildFirstOccurrence(bwt)
    checkpoints, _ = BuildCheckpointArrays(bwt, C)

    results = {}
    for pattern in patterns:
        match_positions = BetterBWMatching(bwt, pattern, first_occurrence, C, checkpoints)
        results[pattern] = sorted([suffix_array[i] for i in match_positions])
    return results

def BetterBWMatching(bwt, pattern, first_occurrence, C, checkpoints):
    top = 0
    bottom = len(bwt) - 1

    while top <= bottom:
        if pattern:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            if symbol in bwt[top:bottom + 1]:
                top = first_occurrence[symbol] + CountSymbol(bwt, symbol, top, C, checkpoints)
                bottom = first_occurrence[symbol] + CountSymbol(bwt, symbol, bottom + 1, C, checkpoints) - 1
            else:
                return []
        else:
            return list(range(top, bottom + 1))
    return []

def BuildBWT(text):
    text += "$"
    suffixes = [(text[i:], i) for i in range(len(text))]
    suffixes.sort()
    bwt = ''.join(text[i-1] if i > 0 else '$' for (_, i) in suffixes)
    suffix_array = [i for (_, i) in suffixes]
    return bwt, suffix_array

def BuildFirstOccurrence(bwt):
    first_col = sorted(bwt)
    first_occurrence = {}
    for i, char in enumerate(first_col):
        if char not in first_occurrence:
            first_occurrence[char] = i
    return first_occurrence

def BuildCheckpointArrays(bwt, C):
    counts = {}
    total = {}
    for char in set(bwt):
        total[char] = 0
        counts[char] = []

    for i in range(len(bwt)):
        if i % C == 0:
            for char in total:
                counts[char].append(total[char])
        total[bwt[i]] += 1

    return counts, total

def CountSymbol(bwt, symbol, pos, C, checkpoints):
    count = 0
    checkpoint_index = pos // C
    count = checkpoints[symbol][checkpoint_index]
    start = checkpoint_index * C
    for i in range(start, pos):
        if bwt[i] == symbol:
            count += 1
    return count


# Example 1
# ----------
text = "AATCGGGTTCAATCGGGGT"
patterns = ["ATCG", "GGGT"]

matches = MultiplePatternMatching(text, patterns)
for pattern in patterns:
    positions = ' '.join(map(str, matches[pattern]))
    print(f"{pattern}: {positions}")

# Output:
    # ATCG: 1 11
    # GGGT: 4 15


# Example 2
# ----------
with open("dataset_30229_4.txt", "r") as file:
    lines = file.read().strip().split()
    text = lines[0]
    patterns = lines[1:]

matches = MultiplePatternMatching(text, patterns)

with open("output.txt", "w") as file:
    for pattern in patterns:
        positions = ' '.join(map(str, matches[pattern]))
        file.write(f"{pattern}: {positions}\n")

# Output: 
    # GTACTGGGT: 2358
    # ACCCGAGAC: 
    # CAGGACGCA: 
    # AAGCGTAAA: 
    # TATTCTTTA: 5856
    # CACCCCCCA: 
    # ACTAACCAC: 
    # TACAGCTTA: 
    # CCGAACGCC: 5585 ...


# -----------------------------------------------
# Multiple Approximate Pattern Matching
# -----------------------------------------------

def ApproximatePatternMatchingBWT(Text, Patterns, d):
    '''Find all approximate (at most d mismatches) occurrences of a collection of patterns in a text.

    Input: A string Text, followed by a collection of space-separated strings Patterns, and an integer d.
    Output: For each string Pattern in Patterns, the string Pattern followed by a colon, followed by a 
    space-separated collection of all positions where Pattern appears as a substring of Text with at most 
    d mismatches.'''

    BWT, SuffixArray = BuildBWT(Text)
    FirstOccurrence = BuildFirstOccurrence(BWT)
    C = 100
    Checkpoints, _ = BuildCheckpointArrays(BWT, C)

    def BWExactMatchLocations(Seed):
        top, bottom = 0, len(BWT) - 1
        for char in reversed(Seed):
            if char not in FirstOccurrence:
                return []
            top = FirstOccurrence[char] + CountSymbol(BWT, char, top, C, Checkpoints)
            bottom = FirstOccurrence[char] + CountSymbol(BWT, char, bottom + 1, C, Checkpoints) - 1
            if top > bottom:
                return []
        return list(range(top, bottom + 1))

    def HammingDistance(s1, s2):
        return sum(a != b for a, b in zip(s1, s2))

    Results = {}
    for Pattern in Patterns:
        n = len(Pattern)
        k = n // (d + 1)
        Matches = set()
        for i in range(d + 1):
            Seed = Pattern[i * k: (i + 1) * k]
            SeedMatches = BWExactMatchLocations(Seed)
            for loc in SeedMatches:
                Start = SuffixArray[loc] - i * k
                if Start < 0 or Start + n > len(Text):
                    continue
                Window = Text[Start:Start + n]
                if HammingDistance(Window, Pattern) <= d:
                    Matches.add(Start)
        Results[Pattern] = sorted(Matches)

    return Results

def BuildBWT(Text):
    Text += "$"
    Suffixes = [(Text[i:], i) for i in range(len(Text))]
    Suffixes.sort()
    BWT = ''.join(Text[i - 1] if i > 0 else '$' for (_, i) in Suffixes)
    SuffixArray = [i for (_, i) in Suffixes]
    return BWT, SuffixArray

def BuildFirstOccurrence(BWT):
    FirstCol = sorted(BWT)
    FirstOccurrence = {}
    for i, char in enumerate(FirstCol):
        if char not in FirstOccurrence:
            FirstOccurrence[char] = i
    return FirstOccurrence

def BuildCheckpointArrays(BWT, C):
    Counts = {}
    Total = {}
    for char in set(BWT):
        Total[char] = 0
        Counts[char] = []
    for i in range(len(BWT)):
        if i % C == 0:
            for char in Total:
                Counts[char].append(Total[char])
        Total[BWT[i]] += 1
    return Counts, Total

def CountSymbol(BWT, Symbol, Pos, C, Checkpoints):
    CheckpointIndex = Pos // C
    Count = Checkpoints[Symbol][CheckpointIndex]
    Start = CheckpointIndex * C
    for i in range(Start, Pos):
        if BWT[i] == Symbol:
            Count += 1
    return Count



# Example 1
# -----------
text = "ACATGCTACTTT"
patterns = ["ATT", "GCC", "GCTA", "TATT"]
d = 1
matches = ApproximatePatternMatchingBWT(text, patterns, d)
for pattern in patterns:
    positions = ' '.join(map(str, matches[pattern]))
    print(f"{pattern}: {positions}")

# Output:
    # ATT: 2 7 8 9
    # GCC: 4
    # GCTA: 4
    # TATT: 6



# Example 2
# -----------
with open("dataset_30230_10.txt", "r") as file:
    lines = file.read().strip().split()
    text = lines[0]
    patterns = lines[1:-1]
    d = int(lines[-1])
    
matches = ApproximatePatternMatchingBWT(text, patterns, d)

with open("output.txt", "w") as file:
    for pattern in patterns:
        positions = ' '.join(map(str, matches[pattern]))
        file.write(f"{pattern}: {positions}\n")


# -----------------------------------------------
# Largest Guaranteed Shared Kmer
# -----------------------------------------------

def LargestSharedKmer(n, d):
    '''Find the largest value of k such that any two strings of length n 
    with at most d mismatches are guaranteed to share a k-mer.
    
    Input: Length of the stringsn and maximum number of mismatches d.
    Output: Largest k satisfying the condition.'''
    
    for k in range(n, 0, -1):
        if n // k > d:
            return k
    return 0  # In case no such k is found


# Say that you know that two strings of length 363 match with at most 5 
# mismatches, but you don’t know what the strings are.  What is the largest 
# value of k such that we can guarantee that the two strings share a k-mer?
# ---------------------------------------------------------------------------
print(LargestSharedKmer(n=363, d=5)) # Output: 60

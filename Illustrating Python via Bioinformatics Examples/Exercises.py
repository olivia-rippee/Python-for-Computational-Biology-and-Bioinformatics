# -----------------------------------------------------------
# Exercise 1: Find pairs of characters
# -----------------------------------------------------------

# Write a function count_pairs(DNA, pair) that returns the number of 
# occurrences of a pair of characters (pair) in a DNA string (DNA). 


def count_pairs(DNA, pair):
    count = 0
    for i in range(len(DNA) - 1):
        if DNA[i:i+2] == pair:
            count += 1
    return count


# Example:
# ---------
DNA = 'ACTGCTATCCATT'
pair = 'AT'
print(count_pairs(DNA, pair))  
# Output: 2

# -----------------------------------------------------------
# Exercise 2: Count substrings
# -----------------------------------------------------------

# This is an extension of Exercise 1: Find pairs of characters: 
# count how many times a certain string appears in another string. 

# Hint. For each match of the first character of the substring 
# in the main string, check if the next n characters in the main 
# string matches the substring, where n is the length of the substring. 
# Use slices like s[3:9] to pick out a substring of s.


def count_substring(DNA, substring):
    count = 0
    sub_len = len(substring)
    for i in range(len(DNA) - sub_len + 1):
        if DNA[i:i+sub_len] == substring:
            count += 1
    return count

# Example
# -------------
DNA = 'ACGTTACGGAACG'
substring = 'ACG'
print(count_substring(DNA, substring))  
# Output: 3



# -----------------------------------------------------------
# Exercise 3: Allow different types for a function argument
# -----------------------------------------------------------

# Consider the family of find_consensus_v* functions from the section 
# Analyzing the Frequency Matrix. The different versions work on different 
# representations of the frequency matrix. Make a unified find_consensus 
# function that accepts different data structures for the frequency_matrix. 
# Test on the type of data structure and perform the necessary actions. 

def find_consensus(frequency_matrix, DNA_length=None):
    """
    Determines the consensus string from a frequency matrix.
    Supports multiple internal representations:
      - list of lists
      - dict of dicts
      - defaultdict of defaultdicts (requires DNA_length)
    """

    consensus = ''
    
    # List of lists (v1, v2)
    if isinstance(frequency_matrix, list) and \
       isinstance(frequency_matrix[0], list):

        base2index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        DNA_length = len(frequency_matrix[0])

        for i in range(DNA_length):
            max_freq = -1
            max_freq_base = None

            for base in 'ACGT':
                freq = frequency_matrix[base2index[base]][i]
                if freq > max_freq:
                    max_freq = freq
                    max_freq_base = base
                elif freq == max_freq:
                    max_freq_base = '-'

            consensus += max_freq_base

    # Dict of dicts (v3)
    elif isinstance(frequency_matrix, dict) and \
         isinstance(frequency_matrix.get('A'), dict):

        DNA_length = len(frequency_matrix['A'])

        for i in range(DNA_length):
            max_freq = -1
            max_freq_base = None

            for base in 'ACGT':
                freq = frequency_matrix[base][i]
                if freq > max_freq:
                    max_freq = freq
                    max_freq_base = base
                elif freq == max_freq:
                    max_freq_base = '-'

            consensus += max_freq_base

    # Defaultdicts (v4) – need explicit DNA_length
    elif hasattr(frequency_matrix, '__getitem__') and \
         hasattr(frequency_matrix.get('A', {}), '__getitem__'):

        if DNA_length is None:
            raise ValueError('DNA_length must be provided for defaultdict input')

        for i in range(DNA_length):
            max_freq = -1
            max_freq_base = None

            for base in 'ACGT':
                freq = frequency_matrix[base][i]
                if freq > max_freq:
                    max_freq = freq
                    max_freq_base = base
                elif freq == max_freq:
                    max_freq_base = '-'

            consensus += max_freq_base

    else:
        raise TypeError("Unsupported frequency_matrix structure")

    return consensus

# List of lists (v1, v2)
# ------------------------
matrix_v1 = [
    [5, 1, 0],  # A
    [0, 3, 0],  # C
    [0, 0, 6],  # G
    [1, 2, 0]   # T
]
print(find_consensus(matrix_v1))  
# Output: ACG



# Dict of dicts (v3)
# ------------------------
matrix_v3 = {
    'A': {0: 5, 1: 1, 2: 0},
    'C': {0: 0, 1: 3, 2: 0},
    'G': {0: 0, 1: 0, 2: 6},
    'T': {0: 1, 1: 2, 2: 0}
}
print(find_consensus(matrix_v3))  
# Output: ACG



# Defaultdict of defaultdicts (v4)
# ------------------------
from collections import defaultdict

matrix_v4 = defaultdict(lambda: defaultdict(int))
matrix_v4['A'][0] = 5
matrix_v4['A'][1] = 1
matrix_v4['A'][2] = 0
matrix_v4['C'][1] = 3
matrix_v4['G'][2] = 6
matrix_v4['T'][0] = 1
matrix_v4['T'][1] = 2

print(find_consensus(matrix_v4, DNA_length=3))  
# Output: ACG

# -----------------------------------------------------------
# Exercise 4: Make a function more robust
# -----------------------------------------------------------

# Consider the function get_base_counts(DNA), which counts how many times 
# A, C, G, and T appears in the string DNA. Unfortunately, this function 
# crashes if other letters appear in DNA. Write an enhanced function 
# get_base_counts2 which solves this problem. 


def get_base_counts2(DNA):
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for base in DNA:
        if base in counts:
            counts[base] += 1
        else:
            pass  # Ignore unknown characters
    return counts


# Example
# -----------
DNA = 'ADLSTTLLD'
result = get_base_counts2(DNA)
print(result)


# -----------------------------------------------------------
# Exercise 5: Find proportion of bases inside/outside exons
# -----------------------------------------------------------

# Consider the lactase gene. What is the proportion of base A
#  inside and outside exons of the lactase gene?

# Hint. Write a function get_exons, which returns all the substrings 
# of the exon regions concatenated. Also write a function get_introns, 
# which returns all the substrings between the exon regions concatenated. 
# The function get_base_frequencies from the section Finding Base Frequencies 
# can then be used to analyze the frequencies of bases A, C, G, and T in 
# the two strings.

def get_exons(DNA, exons):
    """Return concatenated exon regions from DNA."""
    return ''.join(DNA[start:end] for start, end in exons)

def get_introns(DNA, exons):
    """Return concatenated intron regions from DNA."""
    introns = []
    prev_end = 0
    for start, end in sorted(exons):
        introns.append(DNA[prev_end:start])
        prev_end = end
    introns.append(DNA[prev_end:])
    return ''.join(introns)

def get_base_frequencies_v2(DNA):
    """Return proportion of each base in a DNA string."""
    return {base: DNA.count(base) / float(len(DNA)) for base in 'ATGC'}


# Dummy example
# -----------
lactase_DNA = "ATGCGTACGTAGCTAGCTAGCTGATCGTAGCTAGCTAGC"
exon_regions = [(0, 10), (20, 30)]  

exon_seq = get_exons(lactase_DNA, exon_regions)
intron_seq = get_introns(lactase_DNA, exon_regions)

exon_freqs = get_base_frequencies_v2(exon_seq)
intron_freqs = get_base_frequencies_v2(intron_seq)

print(f"Proportion of A in exons: {exon_freqs['A']:.3f}")
print(f"Proportion of A in introns: {intron_freqs['A']:.3f}")



# -----------------------------------------------------------
# Exercise 6: Speed up Markov chain mutation
# -----------------------------------------------------------

# The functions transition and mutate_via_markov_chain were made 
# for being easy to read and understand. Upon closer inspection, 
# we realize that the transition function constructs the interval_limits 
# every time a random transition is to be computed, and we want to run a 
# large number of transitions. By merging the two functions, pre-computing 
# interval limits for each from_base, and adding a loop over N mutations, 
# one can reduce the computation of interval limits to a minimum. Perform 
# such an efficiency enhancement. Measure the CPU time of this new function 
# versus the mutate_via_markov_chain function for 1 million mutations.


# Efficient Function
# ----------------------
import random
import time

def create_markov_chain():
    markov_chain = {}
    for from_base in 'ATGC':
        slice_points = sorted([0] + [random.random() for _ in range(3)] + [1])
        transition_probabilities = [
            slice_points[i+1] - slice_points[i] for i in range(4)
        ]
        markov_chain[from_base] = {
            base: p for base, p in zip('ATGC', transition_probabilities)
        }
    return markov_chain

def check_transition_probabilities(markov_chain):
    for from_base in 'ATGC':
        s = sum(markov_chain[from_base][to_base] for to_base in 'ATGC')
        if abs(s - 1) > 1E-15:
            raise ValueError(f"Wrong sum: {s} for '{from_base}'")

def precompute_cumulative_probs(markov_chain):
    """
    Create a mapping from each base to a list of (limit, base) tuples for fast drawing.
    """
    cumulative = {}
    for base in 'ATGC':
        probs = markov_chain[base]
        limits = []
        total = 0.0
        for to_base in 'ATGC':
            total += probs[to_base]
            limits.append((total, to_base))
        cumulative[base] = limits
    return cumulative

def mutate_many(DNA, cumulative_probs, N):
    DNA_list = list(DNA)
    length = len(DNA_list)

    for _ in range(N):
        i = random.randint(0, length - 1)
        from_base = DNA_list[i]
        r = random.random()
        for limit, to_base in cumulative_probs[from_base]:
            if r < limit:
                DNA_list[i] = to_base
                break
    return ''.join(DNA_list)


# Original Function
# ----------------------
def draw(discrete_probdist):
    limit = 0
    r = random.random()
    for value in discrete_probdist:
        limit += discrete_probdist[value]
        if r < limit:
            return value

def mutate_via_markov_chain(DNA, markov_chain):
    DNA_list = list(DNA)
    mutation_site = random.randint(0, len(DNA_list) - 1)
    from_base = DNA_list[mutation_site]
    to_base = draw(markov_chain[from_base])
    DNA_list[mutation_site] = to_base
    return ''.join(DNA_list)

def mutate_many_original(DNA, markov_chain, N):
    for _ in range(N):
        DNA = mutate_via_markov_chain(DNA, markov_chain)
    return DNA


# Runtime Comparison
# ----------------------
if __name__ == '__main__':
    # Setup
    DNA = ''.join(random.choice('ATGC') for _ in range(1000))
    markov_chain = create_markov_chain()
    check_transition_probabilities(markov_chain)

    # Optimized version
    cumulative_probs = precompute_cumulative_probs(markov_chain)
    N = 1000000
    start_time = time.time()
    mutate_many(DNA, cumulative_probs, N)
    elapsed_optimized = time.time() - start_time
    print(f'Optimized version: {elapsed_optimized:.2f} seconds')
    # Output: 0.88 seconds

    # Original version
    start_time = time.time()
    mutate_many_original(DNA, markov_chain, N)
    elapsed_original = time.time() - start_time
    print(f'Original version: {elapsed_original:.2f} seconds')
    # Output: 12.03 seconds


# -----------------------------------------------------------
# Exercise 7: Extend the constructor in class Gene
# -----------------------------------------------------------

# Modify the constructor in class Gene such that giving no arguments 
# to the constructor makes the class call up the generate_string method 
# (from the DNA_functions module) which generates a random DNA sequence.

from DNA_functions import generate_string, download, read_DNAfile, read_exon_regions, get_base_frequencies, format_frequencies
from region import Region  # Assuming Region is defined in a separate module

class Gene(object):
    def __init__(self, DNA=None, exon_regions=None):
        """
        DNA: string, (urlbase, filename) tuple, or None (generates random DNA)
        exon_regions: list of (start,end) tuples, (urlbase,filename), or None
        """
        if DNA is None:
            DNA = generate_string()  # Generate random DNA
        elif isinstance(DNA, (list, tuple)) and len(DNA) == 2 and all(isinstance(e, str) for e in DNA):
            download(urlbase=DNA[0], filename=DNA[1])
            DNA = read_DNAfile(DNA[1])
        elif isinstance(DNA, str):
            pass
        else:
            raise TypeError(
                f'DNA={DNA} is not a valid string or (urlbase, filename) tuple'
            )

        self._DNA = DNA

        if exon_regions is None:
            self._exon_regions = []
            self._exons = []
            self._introns = [Region(DNA, 0, len(DNA))]  # all is intron
        else:
            if isinstance(exon_regions, (list, tuple)) and \
               len(exon_regions) == 2 and all(isinstance(e, str) for e in exon_regions):
                download(urlbase=exon_regions[0], filename=exon_regions[1])
                exon_regions = read_exon_regions(exon_regions[1])
            elif isinstance(exon_regions, (list, tuple)) and \
                 isinstance(exon_regions[0], (list, tuple)) and \
                 isinstance(exon_regions[0][0], int) and \
                 isinstance(exon_regions[0][1], int):
                pass  # OK
            else:
                raise TypeError(
                    f'exon_regions={exon_regions} is not valid'
                )

            self._exon_regions = exon_regions
            self._exons = [Region(DNA, start, end) for start, end in exon_regions]

            # Compute introns
            self._introns = []
            prev_end = 0
            for start, end in exon_regions:
                if start > prev_end:
                    self._introns.append(Region(DNA, prev_end, start))
                prev_end = end
            if prev_end < len(DNA):
                self._introns.append(Region(DNA, prev_end, len(DNA)))

    def write(self, filename, chars_per_line=70):
        """Write DNA sequence to file."""
        tofile_with_line_sep(self._DNA, filename, chars_per_line)

    def count(self, base):
        """Count occurrences of a base."""
        return self._DNA.count(base)

    def get_base_frequencies(self):
        """Return base frequencies."""
        return get_base_frequencies(self._DNA)

    def format_base_frequencies(self):
        """Return formatted base frequencies."""
        return format_frequencies(self.get_base_frequencies())



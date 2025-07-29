# -----------------------------------------------------------
# Counting Letters in DNA Strings
# -----------------------------------------------------------

# List Iteration
# ------------------------------------------------

def count_v1(dna, base):
    dna = list(dna)  # convert string to list of letters
    i = 0            # counter
    for c in dna:
        if c == base:
            i += 1
    return i



# String Iteration
# ------------------------------------------------

def count_v2(dna, base):
    i = 0 # counter
    for c in dna:
        if c == base:
            i += 1
    return i

dna = 'ATGCGGACCTAT'
base = 'C'
n = count_v2(dna, base)

# printf-style formatting
print '%s appears %d times in %s' % (base, n, dna)

# (new) format string syntax
print '{base} appears {n} times in {dna}'.format(
    base=base, n=n, dna=dna)



# Program Flow
# ------------------------------------------------

def count_v2_demo(dna, base):
    print 'dna:', dna
    print 'base:', base
    i = 0 # counter
    for c in dna:
        print 'c:', c
        if c == base:
            print 'True if test'
            i += 1
    return i

n = count_v2_demo('ATGCGGACCTAT', 'C')



#ipdb> s
#> /some/disk/user/bioinf/src/count_v2_demo.py(2)count_v2_demo()
#1     1 def count_v1_demo(dna, base):
#----> 2     print 'dna:', dna
#      3     print 'base:', base

#ipdb> s
#dna: ATGCGGACCTAT
#> /some/disk/user/bioinf/src/count_v2_demo.py(3)count_v2_demo()
#      2     print 'dna:', dna
#----> 3     print 'base:', base
#      4     i = 0 # counter

#ipdb> print base
#C




# Index Iteration
# ------------------------------------------------

def count_v3(dna, base):
    i = 0 # counter
    for j in range(len(dna)):
        if dna[j] == base:
            i += 1
    return i




# While Loops
# ------------------------------------------------

def count_v4(dna, base):
    i = 0 # counter
    j = 0 # string index
    while j < len(dna):
        if dna[j] == base:
            i += 1
        j += 1
    return i



# Summing a Boolean List
# ------------------------------------------------

def count_v5(dna, base):
    m = []   # matches for base in dna: m[i]=True if dna[i]==base
    for c in dna:
        if c == base:
            m.append(True)
        else:
            m.append(False)
    return sum(m)



# Inline If Test
# ------------------------------------------------

def count_v6(dna, base):
    m = []   # matches for base in dna: m[i]=True if dna[i]==base
    for c in dna:
        m.append(True if c == base else False)
    return sum(m)



# Using Boolean Values Directly
# ------------------------------------------------

def count_v7(dna, base):
    m = []   # matches for base in dna: m[i]=True if dna[i]==base
    for c in dna:
        m.append(c == base)
    return sum(m)



# List Comprehensions
# ------------------------------------------------

def count_v8(dna, base):
    m = [c == base for c in dna]
    return sum(m)

# As a single line:

def count_v9(dna, base):
    return sum([c == base for c in dna])



# Using a Sum Iterator
# ------------------------------------------------

def count_v10(dna, base):
    return sum(c == base for c in dna)



# Extracting Indices
# ------------------------------------------------
# Instead of making a boolean list with elements expressing whether 
# a letter matches the given base or not, we may collect all the 
# indices of the matches. This can be done by adding an if test to 
# the list comprehension:

def count_v11(dna, base):
    return len([i for i in range(len(dna)) if dna[i] == base])

# In an interactive Python shell:
# >>> dna = 'AATGCTTA'
# >>> base = 'A'
# >>> indices = [i for i in range(len(dna)) if dna[i] == base]
# >>> indices
# [0, 1, 7]
# >>> print dna[0], dna[1], dna[7]  # check
# A A A



# Using Python's Library
# ------------------------------------------------

def count_v12(dna, base):
    return dna.count(base)

def compare_efficiency():

    
# -----------------------------------------------------------
# Efficiency Assessment
# -----------------------------------------------------------

# Generating Random DNA Strings
# ------------------------------------------------

# Generating a long string of repeated character
N = 1000000
dna = 'A'*N
# Output: string 'AAA...A, of length N


# Generating a DNA string with a random composition of letters
# make a list of random letters and then join all those letters to a string
import random
alphabet = list('ATGC')
dna = [random.choice(alphabet) for i in range(N)]
dna = ''.join(dna)  # join the list elements to a string


import random
def generate_string(N, alphabet='ACGT'):
    return ''.join([random.choice(alphabet) for i in range(N)])
dna = generate_string(50)



# Measuring CPU Time
# ------------------------------------------------

# Measuring the time spent in a program can be done by the time module

import time
t0 = time.clock()
# do stuff
t1 = time.clock()
cpu_time = t1 - t0



# Running through all our functions made so far and recording timings
import time
functions = [count_v1, count_v2, count_v3, count_v4,
             count_v5, count_v6, count_v7, count_v8,
             count_v9, count_v10, count_v11, count_v12]
timings = []  # timings[i] holds CPU time for functions[i]

for function in functions:
    t0 = time.clock()
    function(dna, 'A')
    t1 = time.clock()
    cpu_time = t1 - t0
    timings.append(cpu_time)


# Iterate over timings and functions simultaneously via zip
for cpu_time, function in zip(timings, functions):
    print '{f:<9s}: {cpu:.2f} s'.format(
        f=function.func_name, cpu=cpu_time)




# -----------------------------------------------------------
# Verifying the Implementations
# -----------------------------------------------------------

# New function that first computes a certainly correct answer to a 
# counting problem and then calls all the count_* functions, stored 
# in the list functions, to check that each call has the correct result
# ------------------------------------------------

def test_count_all():
    dna = 'ATTTGCGGTCCAAA'
    exact = dna.count('A')
    for f in functions:
        if f(dna, 'A') != exact:
            print f.__name__, 'failed'

    # where the correct answer is dna.count('A')




# Revised test function to follow convention
# ------------------------------------------------

def test_count_all():
    dna = 'ATTTGCGGTCCAAA'
    expected = dna.count('A')

    functions = [count_v1, count_v2, count_v3, count_v4,
                 count_v5, count_v6, count_v7, count_v8,
                 count_v9, count_v10, count_v11, count_v12]
    for f in functions:
        success = f(dna, 'A') == expected
        msg = '%s failed' % f.__name__
        assert success, msg



# -----------------------------------------------------------
# Computing Frequencies
# -----------------------------------------------------------


# Separate Frequency Lists
# ------------------------------------------------

def freq_lists(dna_list):
    n = len(dna_list[0])
    A = [0]*n
    T = [0]*n
    G = [0]*n
    C = [0]*n
    for dna in dna_list:
        for index, base in enumerate(dna):
            if base == 'A':
                A[index] += 1
            elif base == 'C':
                C[index] += 1
            elif base == 'G':
                G[index] += 1
            elif base == 'T':
                T[index] += 1
    return A, C, G, T

# >>> for index, base in enumerate(['t', 'e', 's', 't']):
# ...   print index, base
# ...
# 0 t
# 1 e
# 2 s
# 3 t


# Example
dna_list = ['GGTAG', 'GGTAC', 'GGTGC']
A, C, G, T = freq_lists(dna_list)
print A
print C
print G
print T

# Output:
# [0, 0, 0, 2, 0]
# [0, 0, 0, 0, 2]
# [3, 3, 0, 1, 1]
# [0, 0, 3, 0, 0]



# Nested List
# ------------------------------------------------

# >>> frequency_matrix = [A, C, G, T]
# >>> frequency_matrix[2][3]
# 2
# >>> G[3]  # same element
# 2


def freq_list_of_lists_v1(dna_list):
    # Create empty frequency_matrix[i][j] = 0
    # i=0,1,2,3 corresponds to A,T,G,C
    # j=0,...,length of dna_list[0]
    frequency_matrix = [[0 for v in dna_list[0]] for x in 'ACGT']

    for dna in dna_list:
      for index, base in enumerate(dna):
          if base == 'A':
              frequency_matrix[0][index] +=1
          elif base == 'C':
              frequency_matrix[1][index] +=1
          elif base == 'G':
              frequency_matrix[2][index] +=1
          elif base == 'T':
              frequency_matrix[3][index] +=1

    return frequency_matrix

# Example
dna_list = ['GGTAG', 'GGTAC', 'GGTGC']
frequency_matrix = freq_list_of_lists_v1(dna_list)
print frequency_matrix
# Output: [[0, 0, 0, 2, 0], [0, 0, 0, 0, 2], [3, 3, 0, 1, 1], [0, 0, 3, 0, 0]]




# Dictionary for More Convenient Indexing
# ------------------------------------------------

# >>> base2index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
# >>> base2index['G']
# 2

def freq_list_of_lists_v2(dna_list):
    frequency_matrix = [[0 for v in dna_list[0]] for x in 'ACGT']
    base2index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for dna in dna_list:
        for index, base in enumerate(dna):
            frequency_matrix[base2index[base]][index] += 1

    return frequency_matrix



# Numerical Python Array
# ------------------------------------------------

import numpy as np

def freq_numpy(dna_list):
    frequency_matrix = np.zeros((4, len(dna_list[0])), dtype=np.int)
    base2index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for dna in dna_list:
        for index, base in enumerate(dna):
            frequency_matrix[base2index[base]][index] += 1

    return frequency_matrix


# The resulting frequency_matrix object can be indexed as [b][i] or [b,i], 
# with integers b and i. Typically, b will be something line base2index['C'].



# Dictionary of Lists
# ------------------------------------------------

def freq_dict_of_lists_v1(dna_list):
    n = max([len(dna) for dna in dna_list])
    frequency_matrix = {
        'A': [0]*n,
        'C': [0]*n,
        'G': [0]*n,
        'T': [0]*n,
        }
    for dna in dna_list:
        for index, base in enumerate(dna):
            frequency_matrix[base][index] += 1
    return frequency_matrix


# Example
frequency_matrix = freq_dict_of_lists_v1(dna_list)

import pprint   # for nice printout of nested data structures
pprint.pprint(frequency_matrix)

# Output:
#{'A': [0, 0, 0, 2, 0],
# 'C': [0, 0, 0, 0, 2],
# 'G': [3, 3, 0, 1, 1],
# 'T': [0, 0, 3, 0, 0]}


# Compact initilization:
# dict = {key: value for key in some_sequence}

# Example
frequency_matrix = {base: [0]*n for base in 'ACGT'}

# More compact function:
    
def freq_dict_of_lists_v2(dna_list):
    n = max([len(dna) for dna in dna_list])
        # It is also possible to write:
        # n = max(len(dna) for dna in dna_list)
        # n = max(dna_list, key=len)
    frequency_matrix = {base: [0]*n for base in 'ACGT'}
    for dna in dna_list:
        for index, base in enumerate(dna):
            frequency_matrix[base][index] += 1
    return frequency_matrix




# Dictionary of Dictionaries
# ------------------------------------------------

def freq_dict_of_dicts_v1(dna_list):
    n = max([len(dna) for dna in dna_list])
    frequency_matrix = {base: {index: 0 for index in range(n)}
                        for base in 'ACGT'}
    for dna in dna_list:
        for index, base in enumerate(dna):
            frequency_matrix[base][index] += 1
    return frequency_matrix




# Using Dictionaries with Default Values
# ------------------------------------------------

frequency_matrix = {base: {index: 0 for index in range(n)}
                    for base in 'ACGT'}


from collections import defaultdict

def freq_dict_of_dicts_v2(dna_list):
    n = max([len(dna) for dna in dna_list])
    frequency_matrix = {base: defaultdict(lambda: 0)
                        for base in 'ACGT'}
    for dna in dna_list:
        for index, base in enumerate(dna):
            frequency_matrix[base][index] += 1
    return frequency_matrix



# Using Arrays and Vectorization
# ------------------------------------------------

def freq_dict_of_arrays_v1(dna_list):
    n = max([len(dna) for dna in dna_list])
    frequency_matrix = {base: np.zeros(n, dtype=np.int)
                        for base in 'ACGT'}
    for dna in dna_list:
        for index, base in enumerate(dna):
            frequency_matrix[base][index] += 1
    return frequency_matrix


# >>> dna = 'ACAT'
# >>> dna = np.array(dna, dtype='c')
# >>> dna
# array(['A', 'C', 'A', 'T'], dtype='|S1')

# >>> b = dna == 'A'
# >>> b
# array([ True, False,  True, False], dtype=bool)

# >>> i = np.asarray(b, dtype=np.int)
# >>> i
# array([1, 0, 1, 0])
# >>> frequency_matrix['A'] = frequency_matrix['A'] + i

for dna in dna_list:
    dna = np.array(dna, dtype='c')
    for base in 'ACGT':
        b = dna == base
        i = np.asarray(b, dtype=np.int)
        frequency_matrix[base] = frequency_matrix[base] + i


# Final vectorized function:
    
def freq_dict_of_arrays_v2(dna_list):
    n = max([len(dna) for dna in dna_list])
    frequency_matrix = {base: np.zeros(n, dtype=np.int)
                        for base in 'ACGT'}
    for dna in dna_list:
        dna = np.array(dna, dtype='c')
        for base in 'ACCT':
            frequency_matrix[base] += dna == base
    return frequency_matrix


# -----------------------------------------------------------
# Analyzing the Frequency Matrix
# -----------------------------------------------------------


# List of Lists Frequency Matrix
# ------------------------------------------------

def find_consensus_v1(frequency_matrix):
    if isinstance(frequency_matrix, list) and \
       isinstance(frequency_matrix[0], list):
        pass # right type
    else:
        raise TypeError('frequency_matrix must be list of lists')
        
    base2index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    consensus = ''
    dna_length = len(frequency_matrix[0])

    for i in range(dna_length):  # loop over positions in string
        max_freq = -1            # holds the max freq. for this i
        max_freq_base = None     # holds the corresponding base
        for base in 'ATGC':
            if frequency_matrix[base2index[base]][i] > max_freq:
                max_freq = frequency_matrix[base2index[base]][i]
                max_freq_base = base
            elif frequency_matrix[base2index[base]][i] == max_freq:
                max_freq_base = '-' # more than one base as max
        consensus += max_freq_base  # add new base with max freq
    return consensus



       
# Dict of Dicts Frequency Matrix
# ------------------------------------------------

# 1. The base2index dict is no longer needed.
# 2. Access of sublist, frequency_matrix[0], to test for type and length 
#    of the strings, must be replaced by frequency_matrix['A'].


def find_consensus_v3(frequency_matrix):
    if isinstance(frequency_matrix, dict) and \
       isinstance(frequency_matrix['A'], dict):
        pass # right type
    else:
        raise TypeError('frequency_matrix must be dict of dicts')

    consensus = ''
    dna_length = len(frequency_matrix['A'])

    for i in range(dna_length):  # loop over positions in string
        max_freq = -1            # holds the max freq. for this i
        max_freq_base = None     # holds the corresponding base

        for base in 'ACGT':
            if frequency_matrix[base][i] > max_freq:
                max_freq = frequency_matrix[base][i]
                max_freq_base = base
            elif frequency_matrix[base][i] == max_freq:
                max_freq_base = '-' # more than one base as max

        consensus += max_freq_base  # add new base with max freq
    return consensus


# Example
frequency_matrix = freq_dict_of_dicts_v1(dna_list)
pprint.pprint(frequency_matrix)
print find_consensus_v3(frequency_matrix)

# Output:
#{'A': {0: 0, 1: 0, 2: 0, 3: 2, 4: 0},
# 'C': {0: 0, 1: 0, 2: 0, 3: 0, 4: 2},
# 'G': {0: 3, 1: 3, 2: 0, 3: 1, 4: 1},
# 'T': {0: 0, 1: 0, 2: 3, 3: 0, 4: 0}}
# Consensus string: GGTAC


def find_consensus_v4(frequency_matrix, dna_length):
    if isinstance(frequency_matrix, dict) and \
       isinstance(frequency_matrix['A'], dict):
        pass # right type
    else:
        raise TypeError('frequency_matrix must be dict of dicts')

    consensus = ''
    dna_length = len(frequency_matrix['A'])

    for i in range(dna_length):  # loop over positions in string
        max_freq = -1            # holds the max freq. for this i
        max_freq_base = None     # holds the corresponding base

        for base in 'ACGT':
            if frequency_matrix[base][i] > max_freq:
                max_freq = frequency_matrix[base][i]
                max_freq_base = base
            elif frequency_matrix[base][i] == max_freq:
                max_freq_base = '-' # more than one base as max

        consensus += max_freq_base  # add new base with max freq
    return consensus
    
# See Exercise 3
 

# -----------------------------------------------------------
# Dot Plots from Pair of DNA Sequences
# -----------------------------------------------------------

# Using Lists of Lists
# ------------------------------------------------

def dotplot_list_of_lists(dna_x, dna_y):
    dotplot_matrix = [['0' for x in dna_x] for y in dna_y]
    for x_index, x_value in enumerate(dna_x):
        for y_index, y_value in enumerate(dna_y):
            if x_value == y_value:
                dotplot_matrix[y_index][x_index] = '1'
    return dotplot_matrix


# Example
dna_x = 'TAATGCCTGAAT'
dna_y = 'CTCTATGCC'

M = dotplot_list_of_lists(dna_x, dna_x)
for row in M:
    for column in row:
        print column,
    print
    
    
plot = '\n'.join([' '.join(row) for row in dotplot_matrix])

def make_string_expanded(dotplot_matrix):
    rows = []
    for row in dotplot_matrix:
        row_string = ' '.join(row)
        rows.append(row_string)
    plot = '\n'.join(rows)
    return plot

M2 = [['1', '1', '0', '1'],
     ['1', '1', '1', '1'],
     ['0', '0', '1', '0'],
     ['0', '0', '0', '1']]

s = make_string_expanded(M2)




# Using Numerical Python Arrays
# ------------------------------------------------

def dotplot_numpy(dna_x, dna_y):
    dotplot_matrix = np.zeros((len(dna_y), len(dna_x)), np.int)
    for x_index, x_value in enumerate(dna_x):
        for y_index, y_value in enumerate(dna_y):
            if x_value == y_value:
                dotplot_matrix[y_index,x_index] = 1
    return dotplot_matrix

print dotplot_numpy(dna_x, dna_y)



def dotplot_list_of_lists(dna_x, dna_y):
    dotplot_matrix = [['0' for x in dna_x] for y in dna_y]
    for x_index, x_value in enumerate(dna_x):
        for y_index, y_value in enumerate(dna_y):
            if x_value == y_value:
                dotplot_matrix[y_index][x_index] = '1'
    return dotplot_matrix

dna_x = 'TAATGCCTGAAT'
dna_y = 'CTCTATGCC'

M = dotplot_list_of_lists(dna_x, dna_x)
for row in M:
    for column in row:
        print column,
    print

def make_string(dotplot_matrix):
    return '\n'.join([' '.join(row) for row in dotplot_matrix])


def make_string_expanded(dotplot_matrix):
    rows = []
    for row in dotplot_matrix:
        row_string = ' '.join(row)
        rows.append(row_string)
    plot = '\n'.join(rows)
    return plot

M2 = [['1', '1', '0', '1'],
     ['1', '1', '1', '1'],
     ['0', '0', '1', '0'],
     ['0', '0', '0', '1']]

s = make_string_expanded(M2)

# end of testing join operations in detail
print '-------------'
M = dotplot_list_of_lists(dna_x, dna_x)
print make_string(M)
print repr(make_string(M))

import numpy as np

def dotplot_numpy(dna_x, dna_y):
    dotplot_matrix = np.zeros((len(dna_y), len(dna_x)), np.int)
    for x_index, x_value in enumerate(dna_x):
        for y_index, y_value in enumerate(dna_y):
            if x_value == y_value:
                dotplot_matrix[y_index,x_index] = 1
    return dotplot_matrix

print dotplot_numpy(dna_x, dna_y)

dna_x = 'ATTGCAGCTTAAGGAATCGTGCAGATTAAAGGCACCACGAATTAAGACCAGGGACATAA'
dna_y = 'ATTGCAGCGGAAGGAATACTGCAGGTTAAAGTCACCAGGAATTCAGACCAGTTTCATAA'
dna_y = 'ATTGCAGCCTAAGGAATCGTGCAGATTAAACTCACCACGAATTAAGACCAGGGACATAT'
p = dotplot_numpy(dna_x, dna_y)
#from matplotlib.pyplot import plot, show
#p = np.eye(60) # test
#x = [i for i in xrange(p.shape[1]) for j in xrange(p.shape[0]) if p[i,j] == 1]
#y = [j for i in xrange(p.shape[1]) for j in xrange(p.shape[0]) if p[i,j] == 1]
#plot(x, y, 'ro')
#show()

# Test
M = dotplot_list_of_lists(dna_x, dna_y)
for i in xrange(p.shape[0]):
    for j in xrange(p.shape[1]):
        assert p[i,j] == int(M[i][j]), '%s vs %s' % (p[i,j], int(M[i][j]))


        
# -----------------------------------------------------------
# Finding Base Frequencies
# -----------------------------------------------------------

def get_base_counts(dna):
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for base in dna:
        counts[base] += 1
    return counts

def get_base_frequencies_v1(dna):
    counts = get_base_counts(dna)
    return {base: count*1.0/len(dna)
            for base, count in counts.items()}

def get_base_frequencies_v2(dna):
        return {base: dna.count(base)/float(len(dna))
                for base in 'ATGC'}
    
# Example
dna = 'ACCAGAGT'
frequencies = get_base_frequencies_v2(dna)

def format_frequencies(frequencies):
    return ', '.join(['%s: %.2f' % (base, frequencies[base])
                      for base in frequencies])

print "Base frequencies of sequence '%s':\n%s" % \
      (dna, format_frequencies(frequencies))
    
# Output: Base frequencies of sequence 'ACCAGAGT':
#         A: 0.38, C: 0.25, T: 0.12, G: 0.25


# Example
import urllib, os
urlbase = 'http://hplgit.github.com/bioinf-py/data/'
yeast_file = 'yeast_chr1.txt'
if not os.path.isfile(yeast_file):
    url = urlbase + yeast_file
    urllib.urlretrieve(url, filename=yeast_file)


def read_dnafile_v1(filename):
    lines = open(filename, 'r').readlines()
    # Remove newlines in each line (line.strip()) and join
    dna = ''.join([line.strip() for line in lines])
    return dna


def read_dnafile_v2(filename):
    dna = ''
    for line in open(filename, 'r'):
        dna += line.strip()
    return dna

dna = read_dnafile_v2(yeast_file)
yeast_freq = get_base_frequencies_v2(dna)
print "Base frequencies of yeast DNA (length %d):\n%s" % \
      (len(dna), format_frequencies(yeast_freq))

# Output: Base frequencies of yeast DNA (length 230208):
# A: 0.30, C: 0.19, T: 0.30, G: 0.20




# ----------------------------------------------------------- 
# Translating Genes into Proteins
# ----------------------------------------------------------- 

# Download the genetic code
def download(urlbase, filename):
    if not os.path.isfile(filename):
        url = urlbase + filename
        try:
            urllib.urlretrieve(url, filename=filename)
        except IOError as e:
            raise IOError('No Internet connection')
        # Check if downloaded file is an HTML file, which
        # is what github.com returns if the URL is not existing
        f = open(filename, 'r')
        if 'DOCTYPE html' in f.readline():
            raise IOError('URL %s does not exist' % url)
           
            
# Make a dictionary of this file that maps the code (first column) 
# on to the 1-letter name (second column):
def read_genetic_code_v1(filename):
    infile = open(filename, 'r')
    genetic_code = {}
    for line in infile:
        columns = line.split()
        genetic_code[columns[0]] = columns[1]
    return genetic_code 
            

# Download the file, read it, and makE the dictionary
urlbase = 'http://hplgit.github.com/bioinf-py/data/'
genetic_code_file = 'genetic_code.tsv'
download(urlbase, genetic_code_file)
code = read_genetic_code_v1(genetic_code_file)

    
def read_genetic_code_v2(filename):
    return dict([line.split()[0:2] for line in open(filename, 'r')])

def read_genetic_code_v3(filename):
    genetic_code = {}
    for line in open(filename, 'r'):
        columns = line.split()
        genetic_code[columns[0]] = {}
        genetic_code[columns[0]]['1-letter']   = columns[1]
        genetic_code[columns[0]]['3-letter']   = columns[2]
        genetic_code[columns[0]]['amino acid'] = columns[3]
    return genetic_code

def read_genetic_code_v4(filename):
    genetic_code = {}
    for line in open(filename, 'r'):
        c = line.split()
        genetic_code[c[0]] = {
            '1-letter': c[1], '3-letter': c[2], 'amino acid': c[3]}
    return genetic_code


def read_exon_regions_v1(filename):
    positions = []
    infile = open(filename, 'r')
    for line in infile:
        start, end = line.split()
        start, end = int(start), int(end)
        positions.append((start, end))
    infile.close()
    return positions


def read_exon_regions_v2(filename):
    return [tuple(int(x) for x in line.split())
            for line in open(filename, 'r')]

lactase_exon_regions = read_exon_regions_v2(lactase_exon_file)


def create_mRNA(gene, exon_regions):
    mrna = ''
    for start, end in exon_regions:
        mrna += gene[start:end].replace('T','U')
    return mrna

mrna = create_mRNA(lactase_gene, lactase_exon_regions)


def tofile_with_line_sep_v1(text, filename, chars_per_line=70):
    outfile = open(filename, 'w')
    for i in xrange(0, len(text), chars_per_line):
        start = i
        end = start + chars_per_line
        outfile.write(text[start:end] + '\n')
    outfile.close()

output_folder = 'output'
if not os.path.isdir(output_folder):
    os.mkdir(output_folder)
filename = os.path.join(output_folder, 'lactase_mrna.txt')
tofile_with_line_sep_v1(mrna, filename)


def tofile_with_line_sep_v2(text, foldername, filename,
                            chars_per_line=70):
    if not os.path.isdir(foldername):
        os.makedirs(foldername)
    filename = os.path.join(foldername, filename)
    outfile = open(filename, 'w')

    if chars_per_line == 'inf':
        outfile.write(text)
    else:
        for i in xrange(0, len(text), chars_per_line):
            start = i
            end = start + chars_per_line
            outfile.write(text[start:end] + '\n')
    outfile.close()


def create_protein_fixed(mrna, genetic_code):
    protein_fixed = ''
    trans_start_pos = mrna.find('AUG')
    for i in range(len(mrna[trans_start_pos:])/3):
        start = trans_start_pos + i*3
        end = start + 3
        amino = genetic_code[mrna[start:end]]
        if amino == 'X':
            break
        protein_fixed += amino
    return protein_fixed

protein = create_protein_fixed(mrna, genetic_code)
tofile_with_line_sep_v2(protein, 'output',
                        'lactase_protein_fixed.txt', 70)

print '10 last amino acids of the correct lactase protein: ', \
      protein[-10:]
print 'Length of the correct protein: ', len(protein)

# Output: 10 last amino acids of the correct lactase protein:  QQELSPVSSF
#         Length of the correct protein:  1927



# -----------------------------------------------------------
# Some Humans Can Drink Milk, While Others Cannot
# -----------------------------------------------------------

# Congenital lactase deficiency is caused by mutation of T --> A 
# in position 30049 in lactase gene

def congential_lactase_deficiency(
    lactase_gene,
    genetic_code,
    lactase_exon_regions,
    output_folder=os.curdir,
    mrna_file=None,
    protein_file=None):

    pos = 30049
    mutated_gene = lactase_gene[:pos] + 'A' + lactase_gene[pos+1:]
    mutated_mrna = create_mRNA(mutated_gene, lactase_exon_regions)

    if mrna_file is not None:
        tofile_with_line_sep_v2(
            mutated_mrna, output_folder, mrna_file)

    mutated_protein = create_protein_fixed(
        mutated_mrna, genetic_code)

    if protein_file:
        tofile_with_line_sep_v2(
            mutated_protein, output_folder, protein_file)
    return mutated_protein

mutated_protein = congential_lactase_deficiency(
    lactase_gene, genetic_code, lactase_exon_regions,
    output_folder='output',
    mrna_file='mutated_lactase_mrna.txt',
    protein_file='mutated_lactase_protein.txt')

print '10 last amino acids of the mutated lactase protein:', \
      mutated_protein[-10:]
print 'Length of the mutated lactase protein:', \
      len(mutated_protein)


# Output: 10 last amino acids of the mutated lactase protein: GFIWSAASAA
#         Length of the mutated lactase protein: 1389

# See genes2proteins.py



# -----------------------------------------------------------
# Random Mutations of Genes
# -----------------------------------------------------------


# A Simple Mutation Model
# ------------------------------------------------

import random

def mutate_v1(dna):
    dna_list = list(dna)
    mutation_site = random.randint(0, len(dna_list) - 1)
    dna_list[mutation_site] = random.choice(list('ATCG'))
    return ''.join(dna_list)


# Example
dna = 'ACGGAGATTTCGGTATGCAT'
print 'Starting DNA:', dna
print format_frequencies(get_base_frequencies_v2(dna))

nmutations = 10000
for i in range(nmutations):
    dna = mutate_v1(dna)

print 'DNA after %d mutations:' % nmutations, dna
print format_frequencies(get_base_frequencies_v2(dna))

# Output: 
# Starting DNA: ACGGAGATTTCGGTATGCAT
#A: 0.25, C: 0.15, T: 0.30, G: 0.30
# DNA after 10000 mutations: AACCAATCCGACGAGGAGTG
# A: 0.35, C: 0.25, T: 0.10, G: 0.30


# Vectorized Version
# ------------------------------------------------

import numpy as np
# Use integers in random numpy arrays and map these to characters according to
i2c = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

def mutate_v2(dna, N):
    dna = np.array(dna, dtype='c')  # array of characters
    mutation_sites = np.random.random_integers(
        0, len(dna) - 1, size=N)
    # Must draw bases as integers
    new_bases_i = np.random.random_integers(0, 3, size=N)
    # Translate integers to characters
    new_bases_c = np.zeros(N, dtype='c')
    for i in i2c:
        new_bases_c[new_bases_i == i] = i2c[i]
    dna[mutation_sites] = new_bases_c
    return ''.join(dna.tolist())


def generate_string_v1(N, alphabet='ACGT'):
    return ''.join([random.choice(alphabet) for i in range(N)])


def generate_string_v2(N, alphabet='ACGT'):
    # Draw random integers 0,1,2,3 to represent bases
    dna_i = np.random.random_integers(0, 3, N)
    # Translate integers to characters
    dna = np.zeros(N, dtype='c')
    for i in i2c:
        dna[dna_i == i] = i2c[i]
    return ''.join(dna.tolist())



# A Markov Chain Mutation Model
# ------------------------------------------------
        
slice_points = sorted(
    [0] + [random.random() for i in range(3)] + [1])
transition_probabilities = [slice_points[i+1] - slice_points[i]
                            for i in range(4)]

# markov_chain['A'] = {'A': ..., 'C': ..., 'G': ..., 'T': ...}
markov_chain['A'] = {base: p for base, p in
                     zip('ACGT', transition_probabilities)}


def draw(discrete_probdist):
    """
    Draw random value from discrete probability distribution
    represented as a dict: P(x=value) = discrete_probdist[value].
    """
    # Method:
    # http://en.wikipedia.org/wiki/Pseudo-random_number_sampling
    limit = 0
    r = random.random()
    for value in discrete_probdist:
        limit += discrete_probdist[value]
        if r < limit:
            return value
        
        
def create_markov_chain():
    markov_chain = {}
    for from_base in 'ATGC':
        # Generate random transition probabilities by dividing
        # [0,1] into four intervals of random length
       slice_points = sorted(
           [0] + [random.random()for i in range(3)] + [1])
       transition_probabilities = \
           [slice_points[i+1] - slice_points[i] for i in range(4)]
       markov_chain[from_base] = {base: p for base, p
                         in zip('ATGC', transition_probabilities)}
    return markov_chain

mc = create_markov_chain()
print mc
print mc['A']['T'] # probability of transition from A to T
        
        
def check_transition_probabilities(markov_chain):
    for from_base in 'ATGC':
        s = sum(markov_chain[from_base][to_base]
                for to_base in 'ATGC')
        if abs(s - 1) > 1E-15:
            raise ValueError('Wrong sum: %s for "%s"' % \
                             (s, from_base))
        
        
def check_draw_approx(discrete_probdist, N=1000000):
    """
    See if draw results in frequencies approx equal to
    the probability distribution.
    """
    frequencies = {value: 0 for value in discrete_probdist}
    for i in range(N):
        value = draw(discrete_probdist)
        frequencies[value] += 1
    for value in frequencies:
        frequencies[value] /= float(N)
    print ', '.join(['%s: %.4f (exact %.4f)' % \
                     (v, frequencies[v], discrete_probdist[v])
                     for v in frequencies])


def mutate_via_markov_chain(dna, markov_chain):
    dna_list = list(dna)
    mutation_site = random.randint(0, len(dna_list) - 1)
    from_base = dna[mutation_site]
    to_base = draw(markov_chain[from_base])
    dna_list[mutation_site] = to_base
    return ''.join(dna_list)
        
        
# See Exercise 6

dna = 'TTACGGAGATTTCGGTATGCAT'
print 'Starting DNA:', dna
print format_frequencies(get_base_frequencies_v2(dna))

mc = create_markov_chain()
import pprint
print 'Transition probabilities:\n', pprint.pformat(mc)
nmutations = 10000
for i in range(nmutations):
    dna = mutate_via_markov_chain(dna, mc)

print 'DNA after %d mutations (Markov chain):' % nmutations, dna
print format_frequencies(get_base_frequencies_v2(dna))

# Example Output:
# Starting DNA: TTACGGAGATTTCGGTATGCAT
# A: 0.23, C: 0.14, T: 0.36, G: 0.27
# Transition probabilities:
# {'A': {'A': 0.4288890546751146,
#        'C': 0.4219086988655296,
#        'G': 0.00668870644455688,
#        'T': 0.14251354001479888},
#  'C': {'A': 0.24999667668640035,
#        'C': 0.04718309085408834,
#        'G': 0.6250440975238185,
#        'T': 0.0777761349356928},
#  'G': {'A': 0.16022955651881965,
#        'C': 0.34652746609882423,
#        'G': 0.1328031742612512,
#        'T': 0.3604398031211049},
#  'T': {'A': 0.20609823213950174,
#        'C': 0.17641112746655452,
#        'G': 0.010267621176125452,
#        'T': 0.6072230192178183}}
# DNA after 10000 mutations (Markov chain): GGTTTAAGTCAGCTATGATTCT
# A: 0.23, C: 0.14, T: 0.41, G: 0.23


def transition_into_bases(markov_chain):
    return {to_base: sum(markov_chain[from_base][to_base]
                         for from_base in 'ATGC')/4.0
            for to_base in 'ATGC'}

print transition_into_bases(mc)
# Output: {'A': 0.26, 'C': 0.25, 'T': 0.30, 'G': 0.19}

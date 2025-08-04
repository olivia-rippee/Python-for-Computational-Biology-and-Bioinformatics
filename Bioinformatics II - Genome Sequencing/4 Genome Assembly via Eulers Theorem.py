import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics II - Genome Sequencing/Data")

# Hamiltonian cycle: a cycle that visits every node in a graph once
# Eulerian cycle: a cycle that visits every edge in a graph once

# -----------------------------------------------
# Eulerian Cycle
# -----------------------------------------------

from collections import defaultdict

def ParseInputEulerianCycle(data_source):
    '''Parses an adjacency list from either a multiline string or a file path.
    Returns a graph as a defaultdict(list).
    '''
    graph = defaultdict(list)

    # Determine if input is a file path or string
    if os.path.isfile(data_source):
        with open(data_source, 'r') as f:
            lines = f.readlines()
    else:
        lines = data_source.strip().split("\n")

    for line in lines:
        if ':' in line:
            node, neighbors = line.split(":")
            node = int(node.strip())
            neighbors = list(map(int, neighbors.strip().split()))
            graph[node].extend(neighbors)

    return graph


def EulerianCycle(graph):
    '''Generate the Eulerian cycle (if it exists) for a graph.
    An Eulerian cycle visits every edge in a graph once.
    
    Input: graph expressed as a matrix
    Output: Eulerian cycle
    '''
    
    graph_copy = {node: neighbors[:] for node, neighbors in graph.items()}
    start_node = next(iter(graph))  # Start from any node
    stack = [start_node]
    cycle = []

    while stack:
        curr = stack[-1]
        if graph_copy[curr]:
            next_node = graph_copy[curr].pop()
            stack.append(next_node)
        else:
            cycle.append(stack.pop())
    
    return cycle[::-1]  # reverse to get correct order



# Examples
# -------------
input_str = '''
0: 3
1: 0
2: 1 6
3: 2
4: 2
5: 4
6: 5 8
7: 9
8: 7
9: 6
'''
graph = ParseInputEulerianCycle(input_str)
cycle = EulerianCycle(graph)
print(*cycle)
# Output: 0 3 2 6 8 7 9 6 5 4 2 1 0


graph = ParseInputEulerianCycle("dataset_30187_2.txt")
cycle = EulerianCycle(graph)
with open("output.txt", 'w') as out_file:
    out_file.write(" ".join(map(str, cycle)) + "\n")
# Output: 0 29 916 918 917 29 145 198 969 968 967 198 etc.



# -----------------------------------------------
# Eulerian Path
# -----------------------------------------------

from collections import defaultdict, deque

def ParseInputEulerianPath(input_data):
    '''Parses a multiline string or a file path into an adjacency list.'''
    adj_list = defaultdict(list)

    # If it's a valid file path, read from file
    if os.path.isfile(input_data):
        with open(input_data, 'r') as f:
            lines = f.readlines()
    else:
        lines = input_data.strip().split('\n')

    for line in lines:
        if ':' in line:
            node, edges = line.split(':')
            adj_list[int(node.strip())] = list(map(int, edges.strip().split()))

    return adj_list


def EulerianPath(adj_list):
    '''
    Input: adjacency list
    Output: Eulerian path
    '''
    
    # Build the graph
    graph = defaultdict(deque)
    in_deg = defaultdict(int)
    out_deg = defaultdict(int)
    
    for node in adj_list:
        for neighbor in adj_list[node]:
            graph[node].append(neighbor)
            out_deg[node] += 1
            in_deg[neighbor] += 1
            
    # Find start node
    start = None
    for node in set(list(in_deg.keys()) + list(out_deg.keys())):
        if out_deg[node] - in_deg[node] == 1:
            start = node
            break
    
    if not start:
        start = next(iter(adj_list))  # arbitrary start if Eulerian circuit

    path = []
    stack = [start]
    
    while stack:
        while graph[stack[-1]]:
            next_node = graph[stack[-1]].popleft()
            stack.append(next_node)
        path.append(stack.pop())
    
    return path[::-1]  # reverse to get the correct path



# Examples
# -------------
input_str = '''
0: 2
1: 3
2: 1
3: 0 4
6: 3 7
7: 8
8: 9
9: 6
'''

adj_list = ParseInputEulerianPath(input_str)
path = EulerianPath(adj_list)
print(*path)
# Output: 6 7 8 9 6 3 0 2 1 3 4

    
adj_list = ParseInputEulerianPath("dataset_30187_6.txt")
path = EulerianPath(adj_list)
with open("output.txt", "w") as out_file:
    out_file.write(" ".join(map(str, path)))
# Output: 2 1 4 8 7 9 4 6 5 1 0


# -----------------------------------------------
# String Reconstruction (Genome Assembly)
# -----------------------------------------------

from collections import defaultdict, deque

def ParseInputFileReconstruction(input_data):
    with open(input_data, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    k = int(lines[0])
    patterns = lines[1].split()
    return k, patterns

def DeBruijnGraphFromKmers(k, kmers):
    adj_list = defaultdict(list)
    for pattern in kmers:
        prefix = pattern[:-1]
        suffix = pattern[1:]
        adj_list[prefix].append(suffix)
    return adj_list

def EulerianPath(adj_list):
    graph = defaultdict(deque)
    in_deg = defaultdict(int)
    out_deg = defaultdict(int)
    
    for node in adj_list:
        for neighbor in adj_list[node]:
            graph[node].append(neighbor)
            out_deg[node] += 1
            in_deg[neighbor] += 1
    
    # Find the start node
    start = None
    for node in set(list(in_deg.keys()) + list(out_deg.keys())):
        if out_deg[node] - in_deg[node] == 1:
            start = node
            break
    if not start:
        start = next(iter(adj_list))

    path = []
    stack = [start]

    while stack:
        while graph[stack[-1]]:
            stack.append(graph[stack[-1]].popleft())
        path.append(stack.pop())

    return path[::-1]

def PathToGenome(path):
    text = path[0]
    for i in range(1, len(path)):
        text += path[i][-1]
    return text

def StringReconstruction(k, patterns):
    '''The String Reconstruction Problem reduces to finding an 
    Eulerian path in the de Bruijn graph generated from reads/
    
    Input: An integer k followed by a list of k-mers patterns.
    Output: A string Text with k-mer composition equal to patterns.
    '''
    
    dB = DeBruijnGraphFromKmers(k, patterns)
    path = EulerianPath(dB)
    text = PathToGenome(path)
    return text



# Examples
# -----------
k = 4
patterns = ("CTTA", "ACCA", "TACC", "GGCT", "GCTT", "TTAC")
print(StringReconstruction(k, patterns))
# Output: GGCTTACCA


k, patterns = ParseInputFileReconstruction("dataset_30187_7.txt")
result = StringReconstruction(k, patterns)
with open("output.txt", "w") as out_file:
    out_file.write(result.replace(" ", ""))
# Output: CGCGTATGGGGAAGTTGCCTTTCAAGGCGAAATTGATAGAGAATAATG...



# -----------------------------------------------
# Count k-Universal Circular Strings for k
# -----------------------------------------------

def CountkUniversalStrings(k, n=2):     
    '''Returns the number of k-universal circular strings over an alphabet of size n.
    Uses formula: ((n!)^(n^(k-1))) / n^k
    For n=2 (binary), simplified to 2^(2^(k-1) - k).
    '''
    from math import factorial

    if n == 2:
        # Use simplified formula for binary alphabet
        return 2**(2**(k-1) - k)
    else:
        numerator = factorial(n)**(n**(k-1))
        denominator = n**k
        return numerator // denominator


# Example
# ---------
print(CountkUniversalStrings(k=3))
# Output: 2



# -----------------------------------------------
# Find a k-Universal Circular String
# -----------------------------------------------
    
from collections import defaultdict, deque
from itertools import product

def DeBruijnGraphAllKmers(k):
    '''Generate all binary k-mers.'''
    
    kmers = [''.join(p) for p in product('01', repeat=k)]
    
    adj_list = defaultdict(list)
    for pattern in kmers:
        prefix = pattern[:-1]
        suffix = pattern[1:]
        adj_list[prefix].append(suffix)
    return adj_list

def EulerianCycle(adj_list):
    graph = defaultdict(deque)
    for node in adj_list:
        for neighbor in adj_list[node]:
            graph[node].append(neighbor)

    start = next(iter(adj_list))  # arbitrary start node
    path = []
    stack = [start]

    while stack:
        if graph[stack[-1]]:
            stack.append(graph[stack[-1]].popleft())
        else:
            path.append(stack.pop())

    return path[::-1]

def CycleToUniversalString(cycle):
    '''The universal string is constructed by concatenating 
    first node + last chars of all subsequent nodes.
    Cycle is list of (k-1)-mers representing nodes of the De Bruijn graph.
    '''
    
    text = cycle[0]
    for node in cycle[1:]:
        text += node[-1]
    # Since it's circular, trim last k-1 characters to keep length 2^k
    k = len(cycle[0]) + 1
    return text[:2**k]

def kUniversalCircularString(k):
    '''Generate a k-universal circular string (De Bruijn sequence) over binary alphabet {0,1}.
    
    Input: An integer k.
    Output: A k-universal circular string.
    '''
    graph = DeBruijnGraphAllKmers(k)
    cycle = EulerianCycle(graph)
    universal_string = CycleToUniversalString(cycle)
    return universal_string



# Example
# ----------
print(kUniversalCircularString(k=3))
# Output: 00010111


with open("dataset_30187_11.txt", "r") as file:
    k = int(file.readline())
print(kUniversalCircularString(k))
# Output: 000000000100000001100000010100000011100000100100000101
# 10000011010000011110000100010000100110000101010000101110000110
# 01000011011000011101000011111000100011000100101000100111000101
# 00100010101100010110100010111100011001100011010100011011100011
# 10010001110110001111010001111110010010010110010011010010011110
# 01010011001010101001010111001011011001011101001011111001100111
# 00110101100110110100110111100111010100111011100111101100111110
# 10011111110101010110101011110101101110101110110101111110110110
# 111110111011110111111111
    
    
    
# -----------------------------------------------
# (k,d)-mer Composition
# -----------------------------------------------
    
def kdmerComposition(text, k, d):
    '''Generate the (k,d)-mer composition of a string, including repeats, in lexicographic order.
    
    Input: A string of DNA text, length of substrings k, and distance bewteen k-mers d.
    Output: (k,d)-mer composition as a string of tuples separated by spaces.
    '''
    pairs = []
    n = len(text)
    
    # Loop over starting indices for first k-mer
    for i in range(n - (k + d + k) + 1):
        first_kmer = text[i:i + k]
        second_kmer = text[i + k + d:i + k + d + k]
        pair = f"({first_kmer}|{second_kmer})"
        pairs.append(pair)

    # Sort lexicographically
    pairs.sort()
    return ' '.join(pairs)


# Example
# ----------
print(kdmerComposition(text="TAATGCCATGGGATGTT", k=3, d=2))
# Output: (AAT|CAT) (ATG|ATG) (ATG|ATG) (CAT|GAT) (CCA|GGA) (GCC|GGG) (GGG|GTT) (TAA|CCA) (TGC|TGG) (TGG|TGT)



# -----------------------------------------------
# String Spelled by a Gapped Genome Path
# -----------------------------------------------

def StringSpelledByPatterns(patterns, k):
    '''Concatenates the patterns into a string by overlapping k-1 characters.'''
    result = patterns[0]
    for pattern in patterns[1:]:
        result += pattern[-1]
    return result

def StringSpelledByGappedPatterns(gapped_patterns, k, d):
    '''Returns a string spelled by gapped patterns if it exists.
    
    Input: Integers k and d followed by a sequence of (k, d)-mers (a1|b1), … , (an|bn) such that Suffix(ai|bi) = Prefix(ai+1|bi+1) for 1 ≤ i ≤ n-1.
    Output: A string Text of length k + d + k + n - 1 such that the i-th (k, d)-mer in Text is equal to (ai|bi)  for 1 ≤ i ≤ n (if such a string exists).
    '''
    
    # Separate the patterns into two lists: FirstPatterns and SecondPatterns
    first_patterns = [pattern.split('|')[0] for pattern in gapped_patterns]
    second_patterns = [pattern.split('|')[1] for pattern in gapped_patterns]

    # Construct the prefix and suffix strings
    prefix_string = StringSpelledByPatterns(first_patterns, k)
    suffix_string = StringSpelledByPatterns(second_patterns, k)

    # Check consistency of the overlap
    for i in range(k + d, len(prefix_string)):
        if prefix_string[i] != suffix_string[i - k - d]:
            return "there is no string spelled by the gapped patterns"

    # Concatenate the valid parts to form the final string
    return prefix_string + suffix_string[-(k + d):]



# Examples
# -----------
input_patterns = ["GACC|GCGC", "ACCG|CGCC", "CCGA|GCCG","CGAG|CCGG", "GAGC|CGGA"]
print(StringSpelledByGappedPatterns(input_patterns, k=4, d=2))
# Output: GACCGAGCGCCGGA



# Reconstruct the genome spelled by the following Eulerian path of (2,1)-mers: 
# (AG|AG)→(GC|GC)→(CA|CT)→(AG|TG)→(GC|GC)→(CT|CT)→(TG|TG)→(GC|GC)→(CT|CA).
# --------------------------------------------------------------------
input_patterns = ["AG|AG", "GC|GC", "CA|CT", "AG|TG", "GC|GC", "CT|CT", "TG|TG", "GC|GC", "CT|CA"]
print(StringSpelledByGappedPatterns(input_patterns, k=2, d=1))
# Output: AGCAGCTGCTGCA



# The graph has another Eulerian path: (AG|AG)→(GC|GC)→(CT|CT)→(TG|TG)→(GC|GC)→(CA|CT)→(AG|TG)→(GC|GC)→(CT|CA). 
# Can you reconstruct a genome spelled by this path?
# --------------------------------------------------------------------
input_patterns = ["AG|AG", "GC|GC", "CT|CT", "TG|TG", "GC|GC", "CA|CT", "AG|TG", "GC|GC", "CT|CA"]
print(StringSpelledByGappedPatterns(input_patterns, k=2, d=1))
# Output: there is no string spelled by the gapped patterns

# Not every Eulerian path in the paired de Bruijn graph constructed from (k,d)-mer 
# composition spells out a solution of the String Reconstruction from Read-Pairs Problem.



with open("dataset_30208_4.txt") as file:
    k, d = file.readline().split()
    k, d = int(k), int(d)
    PairedReads = file.readline().split()
result = StringSpelledByGappedPatterns(PairedReads, k, d)
    
with open("output.txt", "w") as out_file:
    out_file.write(result)
# Output: CACAGAGAACTCAGCGTCCGCCTCCCTACGTCATTCTTCGCCTTATTTCAGGAATCTCAT....



# -----------------------------------------------
# String Reconstruction from Read-Pairs
# -----------------------------------------------

from collections import defaultdict, deque

def BuildPairedDeBruijnGraph(PairedReads):
    '''Build the paired de Bruijn graph from paired k-mers.'''
    
    graph = defaultdict(list)
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)

    for pair in PairedReads:
        a, b = pair.split('|')
        prefix = (a[:-1], b[:-1])
        suffix = (a[1:], b[1:])
        graph[prefix].append(suffix)
        out_degree[prefix] += 1
        in_degree[suffix] += 1
        
        # Ensure suffix node is in the graph dictionary (even if no outgoing edges)
        if suffix not in graph:
            graph[suffix] = []

        # Also make sure all nodes appear in degree dicts
        if prefix not in in_degree:
            in_degree[prefix] = in_degree[prefix]
        if suffix not in out_degree:
            out_degree[suffix] = out_degree[suffix]

    return graph, in_degree, out_degree


def FindStartNode(graph, in_degree, out_degree):
    '''Find the start node for Eulerian path: 
    node with out_degree = in_degree + 1 (if exists), otherwise any node with edges.
    '''
    
    start = None
    for node in graph:
        outdeg = out_degree[node]
        indeg = in_degree[node]
        if outdeg - indeg == 1:
            return node  # Start node for Eulerian path
    # If no node has out_degree - in_degree = 1, just return any node with outgoing edges
    return next(iter(graph))


def EulerianPath(graph, in_degree, out_degree):
    '''Find an Eulerian path in the graph using Hierholzer’s algorithm.'''
    
    graph_copy = {node: edges[:] for node, edges in graph.items()}
    start_node = FindStartNode(graph, in_degree, out_degree)
    path = []
    stack = [start_node]

    while stack:
        curr = stack[-1]
        if graph_copy[curr]:
            next_node = graph_copy[curr].pop()
            stack.append(next_node)
        else:
            path.append(stack.pop())

    return path[::-1]  # Reverse path to get correct order


def StringSpelledByPatterns(patterns, k):
    '''Concatenates the patterns into a string by overlapping k-1 characters.'''
    
    result = patterns[0]
    for pattern in patterns[1:]:
        result += pattern[-1]
    return result


def StringSpelledByGappedPatterns(gapped_patterns, k, d):
    '''Return string spelled by gapped patterns.'''
    
    first_patterns = [pattern.split('|')[0] for pattern in gapped_patterns]
    second_patterns = [pattern.split('|')[1] for pattern in gapped_patterns]

    prefix_string = StringSpelledByPatterns(first_patterns, k)
    suffix_string = StringSpelledByPatterns(second_patterns, k)

    for i in range(k + d, len(prefix_string)):
        if prefix_string[i] != suffix_string[i - k - d]:
            return "there is no string spelled by the gapped patterns"

    return prefix_string + suffix_string[-(k + d):]


def StringReconstructionFromReadPairs(k, d, PairedReads):
    '''Reconstruct a string from its paired composition.
    Input: Integers k and d followed by a collection of paired k-mers PairedReads.
    Output: A string Text with (k, d)-mer composition equal to PairedReads.
    '''
    graph, in_degree, out_degree = BuildPairedDeBruijnGraph(PairedReads)
    path = EulerianPath(graph, in_degree, out_degree)
    
    # The path is a list of tuples (prefixes of pairs), convert it back to (k,d)-mers
    gapped_patterns = []
    for i in range(len(path) - 1):
        a = path[i][0] + path[i+1][0][-1]
        b = path[i][1] + path[i+1][1][-1]
        gapped_patterns.append(a + "|" + b)
    
    return StringSpelledByGappedPatterns(gapped_patterns, k, d)
    
    
    
# Examples
# ------------
k=4
d=2
PairedReads = ["GAGA|TTGA", "TCGT|GATG", "CGTG|ATGT", "TGGT|TGAG", "GTGA|TGTT", 
               "GTGG|GTGA", "TGAG|GTTG", "GGTC|GAGA", "GTCG|AGAT"]
print(StringReconstructionFromReadPairs(k, d, PairedReads))
# Output: GTGGTCGTGAGATGTTGA


with open("dataset_30188_16.txt") as file:
    k, d = file.readline().split()
    k, d = int(k), int(d)
    PairedReads = file.readline().split()
result = StringReconstructionFromReadPairs(k, d, PairedReads)
    
with open("output.txt", "w") as out_file:
    out_file.write(result)
# Output: CCTCTCTGCTTCCTGTTCTACGCTCACCTGCCGATTAGGCACGAGCGAGAAGAT...


# -----------------------------------------------
# Generate Contigs
# -----------------------------------------------

from collections import defaultdict

def ParseAdjacencyList(input_data):
    '''
    Parses an adjacency list from a multiline string or file.
    Returns a dictionary {node: [neighbors]} with integer keys/values.
    '''
    graph = defaultdict(list)

    # Check if input_data is a filename by testing if it exists as a file
    if isinstance(input_data, str) and os.path.isfile(input_data):
        with open(input_data, 'r') as f:
            lines = f.readlines()
    elif isinstance(input_data, str):
        # Treat input_data as multiline string
        lines = input_data.strip().split('\n')
    else:
        # Assume input_data is an iterable of lines
        lines = [line.strip() for line in input_data if line.strip()]

    for line in lines:
        line = line.strip()
        if not line or ':' not in line:
            continue
        src, dsts = line.split(':', 1)
        src = int(src.strip())
        neighbors = list(map(int, dsts.strip().split()))
        graph[src].extend(neighbors)

    return graph


def InOutDegrees(graph):
    in_deg = defaultdict(int)
    out_deg = defaultdict(int)
    
    for node in graph:
        out_deg[node] = len(graph[node])
        for neighbor in graph[node]:
            in_deg[neighbor] += 1

    # Ensure all nodes are in both dictionaries
    all_nodes = set(graph.keys()) | set(in_deg.keys())
    for node in all_nodes:
        in_deg[node] = in_deg.get(node, 0)
        out_deg[node] = out_deg.get(node, 0)
    
    return in_deg, out_deg


def MaximalNonBranchingPaths(graph):
    '''Generate all non-branching paths in a graph, iterate through all nodes of the 
    graph that are not 1-in-1-out nodes and generate all non-branching paths starting 
    at each such node, and find all isolated cycles in the graph.
    
    Input: The adjacency list of a graph whose nodes are integers.
    Output: The collection of all maximal nonbranching paths in this graph.
    '''
    
    in_deg, out_deg = InOutDegrees(graph)
    paths = []
    visited_edges = defaultdict(set)

    # Step 1: Non-1-in-1-out node paths
    for node in graph:
        if not (in_deg[node] == 1 and out_deg[node] == 1):
            for neighbor in graph[node]:
                if neighbor not in visited_edges[node]:
                    path = [node, neighbor]
                    visited_edges[node].add(neighbor)
                    current = neighbor

                    while (
                        current in graph and
                        in_deg[current] == 1 and
                        out_deg[current] == 1 and
                        graph[current]  # Ensure it has an outgoing edge
                    ):
                        next_node = graph[current][0]
                        path.append(next_node)
                        visited_edges[current].add(next_node)
                        current = next_node

                    paths.append(path)

    # Step 2: Isolated cycles (1-in-1-out nodes not already visited)
    visited_nodes = set()
    for node in graph:
        if in_deg[node] == 1 and out_deg[node] == 1 and node not in visited_nodes:
            cycle = [node]
            visited_nodes.add(node)
            current = graph[node][0]

            while current != node and current not in visited_nodes:
                cycle.append(current)
                visited_nodes.add(current)
                if current not in graph or not graph[current]:
                    break  # No outgoing edge
                current = graph[current][0]

            if current == node:
                cycle.append(current)  # Close the cycle
                paths.append(cycle)

    return paths


# Example
# ------------
input_str = '''
1: 2
2: 3
3: 4 5
6: 7
7: 6
'''
graph = ParseAdjacencyList(input_str)
paths = MaximalNonBranchingPaths(graph)
for path in paths:
    print(" ".join(map(str, path)))
# Output:
# 1 2 3
# 3 4
# 3 5
# 6 7 6


graph = ParseAdjacencyList("input_6.txt")
paths = MaximalNonBranchingPaths(graph)

with open('output.txt', 'w') as out_file:
    for path in paths:
        line = " ".join(map(str, path))
        out_file.write(line + "\n")
# Output:
# 1 2
# 2 5 10 2
# 7 3 4 8 9 7
# 5 10 2 5
# 16 111 16



# -----------------------------------------------
# Generate Contigs
# -----------------------------------------------

from collections import defaultdict

def BuildDeBruijnGraph(kmers):
    '''Build deBruijn graph from k-mers.
    Returns adjacency list: {prefix: [suffixes]} where prefix and suffix are (k-1)-mers.
    '''
    graph = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph[prefix].append(suffix)
    return graph

def InOutDegrees(graph):
    in_deg = defaultdict(int)
    out_deg = defaultdict(int)
    
    for node in graph:
        out_deg[node] = len(graph[node])
        for neighbor in graph[node]:
            in_deg[neighbor] += 1

    # Ensure all nodes are in both dictionaries
    all_nodes = set(graph.keys()) | set(in_deg.keys())
    for node in all_nodes:
        in_deg[node] = in_deg.get(node, 0)
        out_deg[node] = out_deg.get(node, 0)
    
    return in_deg, out_deg


def MaximalNonBranchingPaths(graph):
    '''Generate all non-branching paths in a graph, iterate through all nodes of the graph that 
    are not 1-in-1-out nodes and generate all non-branching paths starting at each such node, 
    and find all isolated cycles in the graph.
    
    Input: The adjacency list of a graph whose nodes are integers.
    Output: The collection of all maximal nonbranching paths in this graph.
    '''
    in_deg, out_deg = InOutDegrees(graph)
    paths = []

    # Count edges
    edge_count = defaultdict(lambda: defaultdict(int))
    for node in graph:
        for neighbor in graph[node]:
            edge_count[node][neighbor] += 1

    visited_edges = defaultdict(lambda: defaultdict(int))

    for node in graph:
        if not (in_deg[node] == 1 and out_deg[node] == 1):
            for neighbor in graph[node]:
                # Visit edges while not all duplicates visited
                while visited_edges[node][neighbor] < edge_count[node][neighbor]:
                    path = [node, neighbor]
                    visited_edges[node][neighbor] += 1
                    current = neighbor

                    while (
                        current in graph and
                        in_deg[current] == 1 and
                        out_deg[current] == 1 and
                        graph[current]
                    ):
                        next_node = graph[current][0]
                        # Check if next edge visited less than count
                        if visited_edges[current][next_node] < edge_count[current][next_node]:
                            path.append(next_node)
                            visited_edges[current][next_node] += 1
                            current = next_node
                        else:
                            break

                    paths.append(path)

    return paths


def PathToGenome(path):
    '''Reconstruct the genome from its genome path.
    
    Input: genome path
    Output: genome
    '''
    
    Text = path[0]
    for i in range(1, len(path)):
        Text += path[i][-1]
    return Text


def GenerateContigs(kmers):
    '''Generate the contigs from a collection of reads (with imperfect coverage).
    
    Input: A collection of k-mers Patterns.
    Output: All contigs in DeBruijn(Patterns).
    '''
    
    graph = BuildDeBruijnGraph(kmers)
    paths = MaximalNonBranchingPaths(graph)
    contigs = [PathToGenome(path) for path in paths]
    return contigs
    
    
    
# Examples
# ---------------
kmers = ["ATG", "ATG", "TGT", "TGG", "CAT", "GGA", "GAT", "AGA"]
contigs = GenerateContigs(kmers)
contigs.sort()
print(" ".join(contigs))
# Output: AGA ATG ATG CAT GAT TGGA TGT



with open('dataset_30189_5.txt', 'r') as file:
    kmers = []
    for line in file:
        kmers.extend(line.strip().split())
contigs = GenerateContigs(kmers)

with open('output.txt', 'w') as out_file:
    out_file.write(" ".join(contigs))
    
# Output: ATATTTCTACATAATAAACAAACAAGATGTTCCATCCGAGCTGTGGAGAAGGCTCTTGTCTAGACACC 
# ATATTTCTACATAATAAACAAACAAGATGTTCCATCCGAGCTGTGGAGAAGGCTCTTGTCTAGACACC 
# CTACAAAGTTTAGATACCCTTGCCGCACTTAACAGTCTCTAGAGTTAACAGGCGCAAGTCACCCCGAA 
# CTACAAAGTTTAGATACCCTTGCCGCACTTAACAGTCTCTAGAGTTAACAGGCGCAAGTCACCCCGAA etc.


kmers = ["AAAT", "AATG", "ACCC", "ACGC", "ATAC", "ATCA", "ATGC", "CAAA", "CACC", "CATA", "CATC", "CCAG", "CCCA", "CGCT", "CTCA", "GCAT", "GCTC", "TACG", "TCAC", "TCAT", "TGCA"]
contigs = GenerateContigs(kmers)
print(contigs[0])




# -----------------------------------------------
# Minimum number of edges to add to make each node balanced
# -----------------------------------------------

from collections import defaultdict

def MinEdgesToBalance(graph):
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)

    # Calculate out-degree and in-degree
    for node, neighbors in graph.items():
        out_degree[node] += len(neighbors)
        for neighbor in neighbors:
            in_degree[neighbor] += 1

    # Get all nodes (some nodes may only appear as destinations)
    all_nodes = set(in_degree.keys()) | set(out_degree.keys())

    # Calculate total positive imbalance
    imbalance_sum = 0
    for node in all_nodes:
        imbalance = out_degree[node] - in_degree[node]
        if imbalance > 0:
            imbalance_sum += imbalance

    return imbalance_sum

# Example
# -----------
graph = {1: [2, 3, 5],
         2: [1, 4],
         3: [2, 5],
         4: [1, 2, 5],
         5: [3, 4]}
print(MinEdgesToBalance(graph))
# Output: 2

graph = {1: [2, 3, 5],
         2: [4],
         3: [2, 5],
         4: [1, 2, 5],
         5: [3]}

print(MinEdgesToBalance(graph))
# Output: 4



# There is a single (linear) string with the following (3,1)-mer composition.  Find it.
# -----------------------------------------------
from collections import defaultdict, deque, Counter

def BuildGraph(paired_kmers):
    graph = defaultdict(list)
    in_degree = Counter()
    out_degree = Counter()
    for left, right in paired_kmers:
        prefix = (left[:-1], right[:-1])
        suffix = (left[1:], right[1:])
        graph[prefix].append(suffix)
        out_degree[prefix] += 1
        in_degree[suffix] += 1
    return graph, in_degree, out_degree

def FindStartNode(graph, in_degree, out_degree):
    start = None
    for node in graph:
        outdeg = out_degree[node]
        indeg = in_degree[node]
        if outdeg - indeg == 1:
            return node  # start node for Eulerian path
        if outdeg > 0:
            start = node
    return start

def EulerianPath(graph, in_degree, out_degree):
    start = FindStartNode(graph, in_degree, out_degree)
    stack = [start]
    path = []
    local_graph = {u: deque(v) for u,v in graph.items()}
    
    while stack:
        u = stack[-1]
        if u in local_graph and local_graph[u]:
            v = local_graph[u].popleft()
            stack.append(v)
        else:
            path.append(stack.pop())
    path.reverse()
    return path

def ReconstructStringFromPath(path, k, d):
    left_string = path[0][0]
    right_string = path[0][1]
    for i in range(1, len(path)):
        left_string += path[i][0][-1]
        right_string += path[i][1][-1]
    # Now merge left and right strings with gap d
    overlap = k + d
    for i in range(overlap, len(left_string)):
        if left_string[i] != right_string[i - overlap]:
            return None  # inconsistent
    return left_string + right_string[-overlap:]

def ParsePairedKmers(raw):
    paired_kmers = []
    # Format: (XXX|XXX)
    items = raw.strip().split()
    for item in items:
        item = item.strip("()")
        left, right = item.split("|")
        paired_kmers.append((left, right))
    return paired_kmers


raw_input = "(ACC|ATA) (ACT|ATT) (ATA|TGA) (ATT|TGA) (CAC|GAT) (CCG|TAC) (CGA|ACT) (CTG|AGC) (CTG|TTC) (GAA|CTT) (GAT|CTG) (GAT|CTG) (TAC|GAT) (TCT|AAG) (TGA|GCT) (TGA|TCT) (TTC|GAA)"


paired_kmers = ParsePairedKmers(raw_input)
k = 3
d = 1
graph, in_degree, out_degree = BuildGraph(paired_kmers)
path = EulerianPath(graph, in_degree, out_degree)
result = ReconstructStringFromPath(path, k, d)
if result:
    print("Reconstructed string:")
    print(result)
else:
    print("No consistent string found.")

# Output: CACCGATACTGATTCTGAAGCTT




# Give a linear string having the following 4-mer composition.
# -----------------------------------------------
from collections import defaultdict, deque

def BuildDeBruijnGraph(kmers):
    graph = defaultdict(list)
    indegree = defaultdict(int)
    outdegree = defaultdict(int)

    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph[prefix].append(suffix)
        outdegree[prefix] += 1
        indegree[suffix] += 1

    return graph, indegree, outdegree

def FindStartNode(graph, indegree, outdegree):
    start = None
    for node in graph:
        # start node: outdegree = indegree + 1 (if exists)
        if outdegree[node] - indegree[node] == 1:
            return node
    # if no such node, return any node with edges
    return next(iter(graph))

def EulerianPath(graph, start):
    stack = [start]
    path = []
    local_graph = {node: deque(neigh) for node, neigh in graph.items()}

    while stack:
        current = stack[-1]
        if current in local_graph and local_graph[current]:
            next_node = local_graph[current].popleft()
            stack.append(next_node)
        else:
            path.append(stack.pop())
    return path[::-1]

def ReconstructString(path):
    # path is a list of 3-mers
    # first node full, then add last char of each subsequent node
    string = path[0]
    for node in path[1:]:
        string += node[-1]
    return string



# Give a linear string having the following 4-mer composition.
# ---------
kmers = ["AAAT", "AATG", "ACCC", "ACGC", "ATAC", "ATCA", "ATGC",
    "CAAA", "CACC", "CATA", "CATC", "CCAG", "CCCA", "CGCT",
    "CTCA", "GCAT", "GCTC", "TACG", "TCAC", "TCAT", "TGCA"]

graph, indegree, outdegree = BuildDeBruijnGraph(kmers)
start_node = FindStartNode(graph, indegree, outdegree)
path = EulerianPath(graph, start_node)
print(ReconstructString(path))
# Output: s

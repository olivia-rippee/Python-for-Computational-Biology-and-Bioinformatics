import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics III - Comparing Genes, Proteins, and Genomes/Data")


# -----------------------------------------------
# Number of permutations of length n
# -----------------------------------------------
import math

def CountPermutations(n):
    permutations = math.factorial(n) * 2**n
    return permutations


# Example
# ----------
CountPermutations(7) # Output: 645120


# -----------------------------------------------
# Number of different reversals that can be applied to a permutation of length n
# -----------------------------------------------

def CountReversals(n):
    n = int(n*(n+1) / 2)
    return n


# Example
# -------
CountReversals(100) # Output: 5050



# Random Breakage model follows an exponential distribution.
# Identity permutation: [+1 +2 +3 ... +n]


# -----------------------------------------------
# Sorting by Reversals (Greedy Algorithm)
# -----------------------------------------------

def GreedySorting(permutation):
    '''Sort by reversals left to right down the sequence. Greedy heuristic provides 
    a poor approximation for the reversal distance.
    
    Input: A permutation (list of integers).
    Output: The sequence of permutations corresponding to applying GreedySorting to P,
    ending with the identity permutation.
    
    Each step of the sequence reflects a reversal applied during the sorting process.
    '''
    
    result = []
    permutation = permutation[:]  # Make a copy to avoid modifying the original list

    for k in range(len(permutation)):
        if abs(permutation[k]) != k + 1:
            # Find index where k+1 or -(k+1) is located
            target = k + 1
            for i in range(k, len(permutation)):
                if abs(permutation[i]) == target:
                    # Reverse and negate the segment from k to i
                    permutation[k:i+1] = [-x for x in reversed(permutation[k:i+1])]
                    result.append(permutation[:])
                    break

        # After reversing, check if the element is -k-1
        if permutation[k] == -(k + 1):
            permutation[k] = -permutation[k]
            result.append(permutation[:])

    return result


# Examples
# ----------
permutation = [-3, +4, +1, +5, -2]
steps = GreedySorting(permutation)

for step in steps:
    print(' '.join(f'{x:+d}' for x in step))

# Output:
    # -1 -4 +3 +5 -2
    # +1 -4 +3 +5 -2
    # +1 +2 -5 -3 +4
    # +1 +2 +3 +5 +4
    # +1 +2 +3 -4 -5
    # +1 +2 +3 +4 -5
    # +1 +2 +3 +4 +5
    
    
with open("dataset_30161_4.txt", "r") as file:
    line = file.readline()
    permutation = list(map(int, line.strip().replace('(', '').replace(')', '').split()))
    
steps = GreedySorting(permutation)

with open("output.txt", "w") as out_file:
    for step in steps:
        out_file.write(' '.join(f'{x:+d}' for x in step) + '\n')

# Output:
    # +1 +109 -82 -73 -122 +77 +54 -110 +10 +11 -79 +45 -27 +94 +99 +85 -37 ...
    # +1 +2 +62 -118 +20 +76 -14 +115 -23 -64 -133 +104 +130 -72 +86 -47 -44 ...
    # +1 +2 +3 +108 -81 +51 -140 -60 -117 -69 -135 -95 +111 +109 -82 -73 -122 ...
    
    
# Compute the number of steps required by GreedySorting to sort the following 
# permutation (i.e., to transform it into the identity permutation.)
# -----------------------------------------------
permutation = [+2, +6, -8, -17, +7, -14, +18, +3, -5, -16, -11, -19, -4, +10, +13, -20, -1, +9, -12, +15]
steps = GreedySorting(permutation)
len(steps)
# Output: 33


permutation = [-16, -20, +12, +18, -14, -17, -15, -6, -8, -19, -11, +13, -10, +4, -5, -2, +7, -3, +1, -9]
steps = GreedySorting(permutation)
len(steps)
# Output: 29





# What is the largest number of reversals that GreedySorting 
# could ever require to sort a permutation of length 100?
# ----------------------------------------------- 
# Swap each k into right position, then swap sign. Then last element is sorted.
n = 100
2*n - 1 # Output: 199



# Adjacency: p[i+1] − p[i] = 1       (e.g. +11 +12, -10 -9, but not +10 +9)
# Breakpoint: p[i+1] − p[i] ==/ 1


# How many permutations on n elements have n + 1 adjacencies?
    # 1; only the identity permutation



# How many permutations of length 200 have exactly 199 adjacencies?

# From the identity perumation (n+1 adjacencies), insert two breakpoints (single reversals).
import math
math.comb(201, 2)    # Output: 20100




# -----------------------------------------------
# Number of Breakpoints in a Permutation
# -----------------------------------------------

def CountBreakpoints(permutation):
    n = len(permutation)
    
    # Add 0 at the beginning and (n + 1) at the end
    extended_permutation = [0] + permutation + [n + 1]

    breakpoints = 0
    for i in range(len(extended_permutation) - 1):
        if extended_permutation[i + 1] - extended_permutation[i] != 1:
            breakpoints += 1

    return breakpoints



# Example
# ---------
permutation_str = "+3 +4 +5 -12 -8 -7 -6 +1 +2 +10 +9 -11 +13 +14"
permutation_list = permutation_str.split()
permutation = [int(x) for x in permutation_list]
print(CountBreakpoints(permutation))  # Output: 8



with open("dataset_30162_6.txt", "r") as file:
    line = file.readline().strip()
    permutation_list = line.split()

permutation = [int(x) for x in permutation_list]
print(CountBreakpoints(permutation))  # Output: 173



permutation_str = "-16 -20 +11 +12 -14 -13 -15 -6 -8 -19 -18 -17 -10 +4 -5 -2 +7 -3 +1 -9"
permutation_list = permutation_str.split()
permutation = [int(x) for x in permutation_list]
print(CountBreakpoints(permutation))  # Output: 17


permutation_str = "+6 -12 -9 +17 +18 -4 +5 -3 +11 +19 +20 +10 +8 +15 -14 -13 +2 +7 -16 -1"
permutation_list = permutation_str.split()
permutation = [int(x) for x in permutation_list]
print(CountBreakpoints(permutation))  # Output: 18



# -----------------------------------------------
# Number of Breakpoints Between Two Permutations
# -----------------------------------------------

def CountBreakpointsBetweenPermutations(perm1, perm2):
    n = len(perm1)
    
    # Step 1: Map each element in perm2 to its position
    pos_in_perm2 = [0] * (n + 1)          # 1-based permutation elements
    for idx, val in enumerate(perm2):
        pos_in_perm2[val] = idx
    
    # Step 2: Count breakpoints in P relative to Q
    breakpoints = 0
    for i in range(n - 1):
        # Check adjacency in perm2
        if abs(pos_in_perm2[perm1[i + 1]] - pos_in_perm2[perm1[i]]) != 1:
            breakpoints += 1
    
    return breakpoints



# How many permutations of length 10 have the property that no 
# reversal applied to P decreases Breakpoints(P)?
math.factorial(10) # Output: 3628800



# -----------------------------------------------
# Transform a circular chromosome into a cycle (sequence of nodes)
# -----------------------------------------------

def ChromosomeToCycle(Chromosome):
    '''
    Input: A chromosome Chromosome containing n synteny blocks.
    Output: The sequence Nodes of integers between 1 and 2n resulting from applying ChromosomeToCycle to Chromosome.
    '''

    Nodes = [0] * (2 * len(Chromosome))  # initialize list of size 2n
    for j in range(1, len(Chromosome) + 1):
        i = Chromosome[j - 1]  # Adjust index for 0-based Python lists
        if i > 0:
            Nodes[2 * j - 2] = 2 * i - 1
            Nodes[2 * j - 1] = 2 * i
        else:
            Nodes[2 * j - 2] = -2 * i
            Nodes[2 * j - 1] = -2 * i - 1
    return Nodes



# Examples
# ----------
chromosome = (1, -2, -3, 4)
print(*ChromosomeToCycle(chromosome))
# Output: 1 2 4 3 6 5 7 8



with open("dataset_30165_4.txt", "r") as file:
    line = file.readline().strip()
    line = line.replace('(', '').replace(')', '')  # Remove parentheses
    chromosome = tuple(int(x) for x in line.split())
    
nodes = []
for i in chromosome:
    if i > 0:
        nodes.append(2 * i - 1)
        nodes.append(2 * i)
    else:
        nodes.append(-2 * i)
        nodes.append(-2 * i - 1)

output = ' '.join(map(str, nodes))

with open("output.txt", "w") as out_file:
    out_file.write(output + '\n')

# Output: 2 1 3 4 6 5 8 7 9 10 12 11 13 14 15 16 18 17 19 20 22 21 24 23 
# 25 26 28 27 30 29 32 31 34 33 36 35 37 38 39 40 41 42 43 44 46 45 48 47 
# 49 50 51 52 53 54 55 56 58 57 59 60 61 62 63 64 66 65 67 68 70 69 72 71 
# 74 73 75 76 78 77 80 79 82 81 84 83 85 86 88 87 89 90 91 92 94 93 96 95 
# 97 98 100 99 102 101 104 103 105 106 108 107 109 110 111 112 114 113 116 
# 115 117 118 119 120 122 121 124 123 125 126 127 128 129 130


# -----------------------------------------------
# Transform a cycle (sequence of nodes) into a circular chromosome
# -----------------------------------------------

def CycleToChromosome(Nodes):
    '''
    Input: A sequence Nodes of integers between 1 and 2n.
    Output: The chromosome Chromosome containing n synteny blocks resulting from applying CycleToChromosome to Nodes.
    '''

    n = len(Nodes) // 2
    Chromosome = [0] * n
    for j in range(1, n + 1):
        if Nodes[2 * j - 2] < Nodes[2 * j - 1]:
            Chromosome[j - 1] = Nodes[2 * j - 1] // 2
        else:
            Chromosome[j - 1] = - (Nodes[2 * j - 2] // 2)
    return Chromosome


# Examples
# ----------
input_str = "(1 2 4 3 6 5 7 8)"
input_str = input_str.strip().replace('(', '').replace(')', '')
Nodes = [int(x) for x in input_str.split()]

chromosome = CycleToChromosome(Nodes)
print('(' + ' '.join(f"+{x}" if x > 0 else str(x) for x in chromosome) + ')')
# Output: (+1 -2 -3 +4)



with open("dataset_30165_5.txt", "r") as file:
    line = file.readline().strip()
    line = line.replace('(', '').replace(')', '')
    Nodes = [int(x) for x in line.split()]

chromosome = CycleToChromosome(Nodes)
output = '(' + ' '.join(f"+{x}" if x > 0 else str(x) for x in chromosome) + ')'

with open("output.txt", "w") as out_file:
    out_file.write(output + '\n')
    
# Output: (-1 +2 +3 +4 +5 +6 +7 -8 -9 -10 -11 -12 +13 +14 +15 +16 -17 -18 +19 
#         +20 -21 -22 -23 -24 +25 -26 -27 -28 +29 +30 +31 -32 -33 -34 +35 +36 
#         -37 +38 +39 -40 -41 +42 -43 -44 -45 +46 -47 -48 -49 -50 +51 -52 +53 
#         -54 +55 +56 -57 +58 +59 +60 +61 -62 -63)

    
    
# -----------------------------------------------
# Construct Colored Edges of Genome Graph
# -----------------------------------------------

def ChromosomeToCycle(Chromosome):
    Nodes = [0] * (2 * len(Chromosome))
    for j in range(1, len(Chromosome) + 1):
        i = Chromosome[j - 1]
        if i > 0:
            Nodes[2 * j - 2] = 2 * i - 1
            Nodes[2 * j - 1] = 2 * i
        else:
            Nodes[2 * j - 2] = -2 * i
            Nodes[2 * j - 1] = -2 * i - 1
    return Nodes


def ColoredEdges(P):
    '''
    Input: A genome P.
    Output: The collection of colored edges in the genome graph of P in the form (x, y).
    '''
    Edges = []
    for Chromosome in P:
        Nodes = ChromosomeToCycle(Chromosome)
        n = len(Chromosome)
        for j in range(n):
            start = Nodes[2 * j + 1]
            end = Nodes[(2 * j + 2) % (2 * n)]
            Edges.append((start, end))
    return Edges


def ParseGenome(genome_str):
    genome_str = genome_str.strip()
    genome_str = genome_str.replace(')(', ')|(')  # split chromosomes
    chromosomes_str = genome_str.split('|')
    chromosomes = []
    for c_str in chromosomes_str:
        c_str = c_str.strip().replace('(', '').replace(')', '')
        chromosome = [int(x) for x in c_str.split()]
        chromosomes.append(chromosome)
    return chromosomes



# Examples
# ----------
sample_input = "(+1 -2 -3)(+4 +5 -6)"
P = ParseGenome(sample_input)
edges = ColoredEdges(P)
print(', '.join(f"({x}, {y})" for x, y in edges))
# Output: (2, 4), (3, 6), (5, 1), (8, 9), (10, 12), (11, 7)



with open("dataset_30165_7.txt", "r") as file:
    genome_str = file.readline().strip()
P = ParseGenome(genome_str)
edges = ColoredEdges(P)

output = ', '.join(f"({x}, {y})" for x, y in edges)
with open("output.txt", "w") as out_file:
    out_file.write(output + '\n')
    
# Output: (2, 3), (4, 6), (5, 8), (7, 9), (10, 11), (12, 14), (13, 16), 
# (15, 18), (17, 20), (19, 21), (22, 24), (23, 26), (25, 27), (28, 30), 
# (29, 31), (32, 34), (33, 35), (36, 37), (38, 39), (40, 41), (42, 44), 
# (43, 46), (45, 48), (47, 1), (50, 52), (51, 54), (53, 55), (56, 58), ...
    
    
    
# -----------------------------------------------
# Graph to Genome
# -----------------------------------------------

from collections import defaultdict

def CycleToChromosome(Nodes):
    n = len(Nodes) // 2
    Chromosome = [0] * n
    for j in range(1, n + 1):
        if Nodes[2 * j - 2] < Nodes[2 * j - 1]:
            Chromosome[j - 1] = Nodes[2 * j - 1] // 2
        else:
            Chromosome[j - 1] = - (Nodes[2 * j - 2] // 2)
    return Chromosome


def GraphToGenome(GenomeGraph):
    '''Convert a genome graph back into a genome.
    
    Input: The colored edges ColoredEdges of a genome graph.
    Output: The genome P corresponding to this genome graph.
    '''
    adj = defaultdict(list)
    for a, b in GenomeGraph:
        adj[a].append(b)
        adj[b].append(a)

    visited = set()
    chromosomes = []

    for node in adj:
        if node in visited:
            continue

        current = node
        cycle = []

        while True:
            visited.add(current)
            cycle.append(current)

            # Get neighbors
            neighbors = adj[current]
            # Choose next node that's unvisited
            for neighbor in neighbors:
                if neighbor not in visited:
                    current = neighbor
                    break
            else:
                break  # No unvisited neighbors — cycle is complete

        if len(cycle) % 2 == 0:
            chromosome = CycleToChromosome(cycle)
            chromosomes.append(chromosome)

    return chromosomes


def FormatGenome(P):
    return ''.join(
        '(' + ' '.join(f"+{x}" if x > 0 else f"{x}" for x in chromosome) + ')'
        for chromosome in P)


# Example
# ---------
GenomeGraph = [(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)]
result = GraphToGenome(GenomeGraph)
print(FormatGenome(result))
    

    
    
# -----------------------------------------------
# Calculate Two-Break Distance
# -----------------------------------------------
    
import re

def ParseGenome(genome_str):
    # Parse genome from string like "(+1 +2 +3 +4 +5 +6)(+7 +8)"
    chromosomes = []
    for chromosome in re.findall(r'\((.*?)\)', genome_str):
        blocks = list(map(lambda x: int(x), chromosome.split()))
        chromosomes.append(blocks)
    return chromosomes


def ReadGenomesFromFile(filename):
    with open(filename, "r") as file:
        P_str = file.readline().strip()
        Q_str = file.readline().strip()
    return P_str, Q_str


def GetEdges(chromosomes):
    edges = []
    for chromosome in chromosomes:
        nodes = []
        for block in chromosome:
            head = 2 * abs(block) - 1
            tail = 2 * abs(block)
            if block > 0:
                nodes.append((head, tail))
            else:
                nodes.append((tail, head))
        
        for i in range(len(nodes)):
            # Connect tail of current block to head of next block (circular)
            current = nodes[i]
            next_block = nodes[(i + 1) % len(nodes)]
            edges.append((current[1], next_block[0]))
    return edges


def BuildBreakpointGraph(P, Q):
    '''Nodes are integers representing block ends.'''
    
    graph = {}
    
    def AddEdge(u, v):
        graph.setdefault(u, []).append(v)
        graph.setdefault(v, []).append(u)
    
    P_edges = GetEdges(P)
    Q_edges = GetEdges(Q)
    
    for u, v in P_edges:
        AddEdge(u, v)
    for u, v in Q_edges:
        AddEdge(u, v)
    return graph

def CountCycles(graph):
    visited = set()
    cycles = 0
    
    def Dfs(node):
        stack = [node]
        while stack:
            current = stack.pop()
            if current not in visited:
                visited.add(current)
                for neighbor in graph[current]:
                    if neighbor not in visited:
                        stack.append(neighbor)
    
    for node in graph:
        if node not in visited:
            Dfs(node)
            cycles += 1
    return cycles


def CountBlocks(genomes):
    # Count total blocks (genes) from genome P (or Q, same)
    total = 0
    for chromosome in genomes:
        total += len(chromosome)
    return total


def TwoBreakDistance(P_str, Q_str):
    P = ParseGenome(P_str)
    Q = ParseGenome(Q_str)
    
    graph = BuildBreakpointGraph(P, Q)
    cycles = CountCycles(graph)
    blocks = CountBlocks(P)
    
    return blocks - cycles



# Example
# ---------
P = "(+1 +2 +3 +4 +5 +6)"
Q = "(+1 -3 -6 -5)(+2 -4)"
print(TwoBreakDistance(P, Q))  # Output: 3


P, Q = ReadGenomesFromFile("dataset_30163_4.txt")
print(TwoBreakDistance(P, Q))  # Output: 9578



# -----------------------------------------------
# Two-Break Genome Graph
# -----------------------------------------------

import re

def TwoBreakOnGenomeGraph(GenomeGraph, i1, i2, i3, i4):
    """
    Input: The colored edges of a genome graph GenomeGraph, followed by indices i1, i2, i3, and i4.
    Output: The colored edges of the genome graph resulting from applying the 2-break operation 
            2-BreakOnGenomeGraph(GenomeGraph, i1, i2, i3, i4).
    """
    # Remove edges (i1, i2) and (i3, i4) - order doesn't matter (undirected edges)
    newEdges = []
    for edge in GenomeGraph:
        if (edge == (i1, i2) or edge == (i2, i1)) or (edge == (i3, i4) or edge == (i4, i3)):
            continue
        newEdges.append(edge)
    
    # Add new edges (i1, i3) and (i2, i4)
    newEdges.append((i1, i3))
    newEdges.append((i2, i4))
    
    return newEdges


def ReadInput(filename):
    with open(filename, "r") as file:
        edges_line = file.readline().strip()
        indices_line = file.readline().strip()
    
    edges = re.findall(r'\((\d+), (\d+)\)', edges_line)
    GenomeGraph = [(int(a), int(b)) for a, b in edges]
    
    i1, i2, i3, i4 = map(int, indices_line.split(', '))
    return GenomeGraph, i1, i2, i3, i4



# Examples
# ---------
input_edges_str = "(2, 4), (3, 8), (7, 5), (6, 1)"
indices_str = "1, 6, 3, 8"

edges = re.findall(r'\((\d+), (\d+)\)', input_edges_str)
GenomeGraph = [(int(a), int(b)) for a,b in edges]
i1, i2, i3, i4 = map(int, indices_str.split(', '))

result = TwoBreakOnGenomeGraph(GenomeGraph, i1, i2, i3, i4)
print(', '.join(f'({a}, {b})' for a,b in result))

# Output: (2, 4), (3, 1), (7, 5), (6, 8)



GenomeGraph, i1, i2, i3, i4 = ReadInput("dataset_30166_2.txt")
edges = TwoBreakOnGenomeGraph(GenomeGraph, i1, i2, i3, i4)
with open("output.txt", "w") as out_file:
    out_file.write(', '.join(f'({a}, {b})' for a,b in edges) + '\n')

# Output: (1, 3), (4, 5), (6, 7), (8, 10), (9, 11), (12, 14), (13, 15), 
# (16, 17), (18, 20), (19, 22), (21, 24), (23, 26), (25, 28), (27, 29), 
# (30, 32), (31, 34), (33, 36), (35, 38), (37, 40), (39, 41), (42, 44), ...



# -----------------------------------------------
# Two-Break On Genome
# -----------------------------------------------

import re

def ChromosomeToCycle(chromosome):
    nodes = []
    for block in chromosome:
        if block > 0:
            nodes.extend([2*block - 1, 2*block])
        else:
            block = -block
            nodes.extend([2*block, 2*block - 1])
    return nodes


def CycleToChromosome(cycle):
    chromosome = []
    for i in range(0, len(cycle), 2):
        head, tail = cycle[i], cycle[i+1]
        if head < tail:
            chromosome.append(tail//2)
        else:
            chromosome.append(-head//2)
    return chromosome


def BlackEdges(P):
    edges = []
    for chromosome in P:
        nodes = ChromosomeToCycle(chromosome)
        length = len(nodes)
        for j in range(0, length, 2):
            edge = (nodes[j+1], nodes[(j+2) % length])
            edges.append(edge)
    return edges


def ColoredEdges(P):
    '''Colored edges are the black edges of genome P.'''
    return BlackEdges(P)


def GraphToGenome(GenomeGraph):
    # Build adjacency dict from edges
    adj = {}
    for (a, b) in GenomeGraph:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)

    chromosomes = []
    visited = set()

    for node in adj:
        if node in visited:
            continue
        cycle = []
        current = node
        while True:
            visited.add(current)
            cycle.append(current)
            # Find the next node in adj[current] not visited yet, or the one which completes the cycle
            neighbors = adj[current]
            next_node = None
            for nbr in neighbors:
                if nbr not in visited:
                    next_node = nbr
                    break
            if next_node is None:
                # no unvisited neighbors left; cycle complete
                break
            current = next_node
        chromosome = CycleToChromosome(cycle)
        chromosomes.append(chromosome)
    return chromosomes


def TwoBreakOnGenomeGraph(GenomeGraph, i1, i2, i3, i4):
    newEdges = []
    for edge in GenomeGraph:
        if (edge == (i1, i2) or edge == (i2, i1)) or (edge == (i3, i4) or edge == (i4, i3)):
            continue
        newEdges.append(edge)
    newEdges.append((i1, i3))
    newEdges.append((i2, i4))
    return newEdges


def TwoBreakOnGenome(P, i1, i2, i3, i4):
    GenomeGraph = BlackEdges(P) + ColoredEdges(P)  # typically, the graph is black + colored edges
    GenomeGraph = TwoBreakOnGenomeGraph(GenomeGraph, i1, i2, i3, i4)
    P_new = GraphToGenome(GenomeGraph)
    return P_new


def ParseGenome(genome_str):
    '''Parses genome strings like "(+1 -2 -4 +3)(+5 +6)" into list of chromosomes.'''
    chromosomes = []
    chrom_strs = re.findall(r'\(([^)]+)\)', genome_str)
    for c in chrom_strs:
        blocks = list(map(int, c.strip().split()))
        chromosomes.append(blocks)
    return chromosomes


def FormatGenome(P):
    return ''.join('(' + ' '.join(f'{("+" if x > 0 else "")}{x}' for x in chrom) + ')' for chrom in P)


def ReadInput(filename):
    with open(filename, "r") as file:
        genome_line = file.readline().strip()
        indices_line = file.readline().strip()
    P = ParseGenome(genome_line)
    i1, i2, i3, i4 = map(int, indices_line.split(', '))
    return P, i1, i2, i3, i4


def WriteOutput(filename, P):
    with open(filename, "w") as file:
        file.write(FormatGenome(P) + "\n")
        
        

# Examples
# -------------
P = [[+1, -2, -4, +3]]
i1, i2, i3, i4 = 1, 6, 3, 8
P_prime = TwoBreakOnGenome(P, i1, i2, i3, i4)
print(FormatGenome(P_prime))
# Output: (+1 -2)(-3 +4)



P_file, i1_file, i2_file, i3_file, i4_file = ReadInput(".txt")
P_prime_file = TwoBreakOnGenome(P_file, i1_file, i2_file, i3_file, i4_file)
WriteOutput("output.txt", P_prime_file)



# -----------------------------------------------
# Two-Break Sorting
# -----------------------------------------------
import sys
import numpy as np
import copy

class TwoBreakSorting:
    def __init__(self):
        genomes = self.InputGenomes()
        result = self.ShortestRearrangement(genomes[0], genomes[1])
        self.SaveResult(result)

    def InputGenomes(self):
        '''Read genomes from stdin, parse, and return them.'''
        data = sys.stdin.read().strip().split('\n')
        return self.ParseGenomes(data)

    def ReadGenomesFromFile(self, filename='input.txt'):
        '''Read genomes from a file, parse, and return them.'''
        with open(filename, 'r') as f:
            data = [line.strip() for line in f]
        return self.ParseGenomes(data)

    def ParseGenomes(self, data):
        genomes = []
        for g in data:
            g = g.split(')(')
            genome = []
            for d in g:
                d = d.split()
                genome.append(
                    [int(d[0][1:] if d[0].startswith('(') else d[0])] + 
                    [int(e) for e in d[1:-1]] + 
                    [int(d[-1][:-1] if d[-1].endswith(')') else d[-1])]
                )
            genomes.append(genome)
        return genomes
    
    def SaveResult(self, result, filename='result.txt'):
        with open(filename, 'w') as f:
            for r in result:
                f.write(self.PrintGenome(r) + '\n')

    def ChromosomeToCycle(self, chromosome):
        l = len(chromosome)
        nodes = [0] * (2 * l)
        for j in range(l):
            i = chromosome[j]
            if i > 0:
                nodes[2*j] = 2*i - 1
                nodes[2*j + 1] = 2*i
            else:
                nodes[2*j] = -2*i
                nodes[2*j + 1] = -2*i - 1
        return nodes
        
    def CycleToChromosome(self, nodes):
        l = len(nodes) // 2
        chromosome = [0] * l
        for j in range(l):
            if nodes[2*j] < nodes[2*j + 1]:
                chromosome[j] = nodes[2*j + 1] // 2
            else:
                chromosome[j] = -nodes[2*j] // 2
        return chromosome
    
    def PrintChromosome(self, chromosome):
        print('(' + ' '.join(['+' + str(e) if e > 0 else str(e) for e in chromosome]) + ')')

    def ColoredEdges(self, genome):
        edges = set()
        for chromosome in genome:
            nodes = self.ChromosomeToCycle(chromosome)
            nodes.append(nodes[0])
            for j in range(len(chromosome)):
                edges.add((nodes[2*j + 1], nodes[2*j + 2]))
        return edges
        
    def PrintGenome(self, genome):
        result = ''.join(
            '(' + ' '.join(['+' + str(e) if e > 0 else str(e) for e in chromosome]) + ')'
            for chromosome in genome
        )
        print(result)
        return result

    def TwoBreakOnGraph(self, edges, i0, i1, j0, j1):
        edges.discard((i0, i1))
        edges.discard((i1, i0))
        edges.discard((j0, j1))
        edges.discard((j1, j0))
        edges.add((i0, j0))
        edges.add((i1, j1))
        return edges

    def GroupNodes(self, edges):
        parent = {}
        rank = {}
        for e in edges:
            parent.setdefault(e[0], e[0])
            parent.setdefault(e[1], e[1])
            rank.setdefault(e[0], 0)
            rank.setdefault(e[1], 0)

        def FindParent(i):
            if i != parent[i]:
                parent[i] = FindParent(parent[i])
            return parent[i]
        
        def Union(i, j):
            i_id = FindParent(i)
            j_id = FindParent(j)
            if i_id == j_id:
                return
            if rank[i_id] > rank[j_id]:
                parent[j_id] = i_id
            else:
                parent[i_id] = j_id
                if rank[i_id] == rank[j_id]:
                    rank[j_id] += 1
        
        def UnionEdges(edge):
            Union(edge[0], edge[1])
            Union(edge[0], edge[0] + 1 if edge[0] % 2 == 1 else edge[0] - 1)
            Union(edge[1], edge[1] + 1 if edge[1] % 2 == 1 else edge[1] - 1)

        for e in edges:
            UnionEdges(e)

        nodesID = {}
        nodesSets = set()

        for e in edges:
            id = FindParent(e[0])
            nodesID[e[0]] = id
            nodesID[e[1]] = id
            nodesSets.add(id)
        
        return nodesSets, nodesID
    
    def BuildEdgeDict(self, edges, nodesSet, nodesID):
        edgeDict = {}
        for e in edges:
            id = nodesID[e[0]]
            if id not in edgeDict:
                edgeDict[id] = {}
            edgeDict[id][e[0]] = e[1]
            edgeDict[id][e[1]] = e[0]
        return edgeDict
            
    def TwoBreakOnGenome(self, genome, i0, i1, j0, j1):
        edges = self.TwoBreakOnGraph(self.ColoredEdges(genome), i0, i1, j0, j1)
        nodesSet, nodesID = self.GroupNodes(edges)
        edgeDict = self.BuildEdgeDict(edges, nodesSet, nodesID)
        nodesDict = {}
        for id, eDict in edgeDict.items():
            nodesDict[id] = []
            currNode0 = list(eDict)[0]
            while eDict:
                nodesDict[id].append(currNode0)
                currNode1 = currNode0 + 1 if currNode0 % 2 == 1 else currNode0 - 1
                nodesDict[id].append(currNode1)
                newNode = eDict[currNode1]
                del eDict[currNode0]
                del eDict[currNode1]
                currNode0 = newNode
        newGenome = {id: self.CycleToChromosome(nodes) for id, nodes in nodesDict.items()}
        return sorted(newGenome.values(), key=lambda x: abs(x[0]))
    
    def EdgeFromNontrivialCycle(self, edges, redEdges, blueEdges, blocks):
        parent = {}
        rank = {}
        for e in edges:
            parent.setdefault(e[0], e[0])
            parent.setdefault(e[1], e[1])
            rank.setdefault(e[0], 0)
            rank.setdefault(e[1], 0)

        def FindParent(i):
            if i != parent[i]:
                parent[i] = FindParent(parent[i])
            return parent[i]
        
        def Union(i, j):
            i_id = FindParent(i)
            j_id = FindParent(j)
            if i_id == j_id:
                return
            if rank[i_id] > rank[j_id]:
                parent[j_id] = i_id
            else:
                parent[i_id] = j_id
                if rank[i_id] == rank[j_id]:
                    rank[j_id] += 1

        for e in edges:
            Union(e[0], e[1])

        nodesID = {}
        nodesSets = set()

        for e in edges:
            id = FindParent(e[0])
            nodesID[e[0]] = id
            nodesID[e[1]] = id
            nodesSets.add(id)
        
        cycles = len(nodesSets)
        hasNontrivialCycle = False
        edge = None
        removedRedEdges = []
        if cycles != blocks:
            hasNontrivialCycle = True
            edgeDict = {}
            redEdgeDict = {}
            for e in edges:
                id = nodesID[e[0]]
                if id not in edgeDict:
                    edgeDict[id] = {}
                edgeDict[id][e[0]] = e[1]
                edgeDict[id][e[1]] = e[0]
                if edge is None and len(edgeDict[id]) > 2 and e in blueEdges:
                    edge = (e[0], e[1])
                    edgeID = id
                if e in redEdges:
                    if id not in redEdgeDict:
                        redEdgeDict[id] = {}
                    redEdgeDict[id][e[0]] = e[1]
                    redEdgeDict[id][e[1]] = e[0]
            removedRedEdges.append((edge[0], redEdgeDict[edgeID][edge[0]]))
            removedRedEdges.append((edge[1], redEdgeDict[edgeID][edge[1]]))
        return hasNontrivialCycle, removedRedEdges        

    def ShortestRearrangement(self, P, Q):
        blocks = sum(len(a) for a in P)
        result = [P]
        redEdges = self.ColoredEdges(P)
        blueEdges = self.ColoredEdges(Q)
        breakpointGraph = redEdges.union(blueEdges)
        hasNontrivialCycle, removedRedEdges = self.EdgeFromNontrivialCycle(breakpointGraph, redEdges, blueEdges, blocks)
        while hasNontrivialCycle:
            redEdges = self.TwoBreakOnGraph(redEdges, removedRedEdges[0][0], removedRedEdges[0][1], removedRedEdges[1][0], removedRedEdges[1][1])
            breakpointGraph = redEdges.union(blueEdges)
            P = self.TwoBreakOnGenome(P, removedRedEdges[0][0], removedRedEdges[0][1], removedRedEdges[1][0], removedRedEdges[1][1])
            hasNontrivialCycle, removedRedEdges = self.EdgeFromNontrivialCycle(breakpointGraph, redEdges, blueEdges, blocks)
            result.append(P)
        return result

            
# Examples
# -------------
P = [[1, -2, -3, 4]]
Q = [[1, 2, -4, -3]]

twoBreak = TwoBreakSorting.__new__(TwoBreakSorting)  
results = twoBreak.ShortestRearrangement(P, Q)
for genome in results:
    twoBreak.PrintGenome(genome)
# Output:
    # (+1 -2 -3 +4)
    # (-1 +3 +2)(+4)
    # (-1 +3 -2)(-4)
    # (+4 -2 -1 +3)


twoBreak = TwoBreakSorting.__new__(TwoBreakSorting)  
genomes = twoBreak.ReadGenomesFromFile("dataset_30163_5.txt")
results = twoBreak.ShortestRearrangement(genomes[0], genomes[1])
twoBreak.SaveResult(results, "output.txt")



# -----------------------------------------------
# Identify and Count Shared kmers
# -----------------------------------------------

from collections import Counter

def ReverseComplement(kmer):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(kmer))


def GetKmerCounts(seq, k):
    counts = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        counts[kmer] += 1
    return counts


def SharedKmers(seq1, seq2, k):
    '''Find all shared kmers between two strings, including reverse complements.'''
    
    seq1_kmers = GetKmerCounts(seq1, k)
    seq2_kmers = GetKmerCounts(seq2, k)

    shared_kmers = []

    for kmer1, count1 in seq1_kmers.items():
        rc1 = ReverseComplement(kmer1)
        count_kmer2 = seq2_kmers.get(kmer1, 0)
        count_rc2 = seq2_kmers.get(rc1, 0)

        # Avoid double counting if kmer1 == rc1 (palindrome)
        if kmer1 == rc1:
            total_matches = count_kmer2  # don't add count_rc2 twice
        else:
            total_matches = count_kmer2 + count_rc2

        # Append kmer1 repeated for every match instance
        shared_kmers.extend([kmer1] * (count1 * total_matches))

    return shared_kmers



# Example 
# -----------
seq1, seq2 = "AAACTCATC", "TTTCAAATC"
kmers = SharedKmers(seq1, seq2, k=2)
print(kmers) 
# Output: ['AA', 'AA', 'AA', 'AA', 'AA', 'AA', 'AA', 
#          'AA', 'TC', 'TC', 'TC', 'TC', 'CA', 'AT']
print(len(kmers)) # Output: 14



seq1, seq2 = "TGGCCTGCACGGTAG", "GGACCTACAAATGGC"
kmers = SharedKmers(seq1, seq2, k=3)
print(len(kmers)) # Output: 7


seq1, seq2 = "TGCCCCGGTGGTGAG", "AAGGTCGCACCTCGT"
kmers = SharedKmers(seq1, seq2, k=3)
print(len(kmers)) # Output: 8



# -----------------------------------------------
# Position of Shared kmers
# -----------------------------------------------

def ReverseComplement(kmer):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(kmer))


def SharedKmersPositions(seq1, seq2, k):
    '''Given two strings, find all their shared k-mers or their reverse complements.

    Input: An integer k and two strings

    Output:  A list of ordered pairs (x, y), where x is the starting position of 
    the k-mer in seq1 and y is the starting position of the matching k-mer or 
    reverse complement in seq2.
    '''
    
    kmer_positions_seq1 = {}
    kmer_positions_seq2 = {}

    # Store all k-mers positions in seq1
    for i in range(len(seq1) - k + 1):
        kmer = seq1[i:i+k]
        if kmer not in kmer_positions_seq1:
            kmer_positions_seq1[kmer] = []
        kmer_positions_seq1[kmer].append(i)

    # Store all k-mers positions in seq2
    for j in range(len(seq2) - k + 1):
        kmer = seq2[j:j+k]
        if kmer not in kmer_positions_seq2:
            kmer_positions_seq2[kmer] = []
        kmer_positions_seq2[kmer].append(j)

    shared_positions = []

    # For each kmer in seq1, check for kmer or its reverse complement in seq2
    for kmer1, positions1 in kmer_positions_seq1.items():
        rc_kmer1 = ReverseComplement(kmer1)
        
        # Check if kmer or reverse complement is in seq2
        seq2_positions = []
        if kmer1 in kmer_positions_seq2:
            seq2_positions.extend(kmer_positions_seq2[kmer1])
        if rc_kmer1 in kmer_positions_seq2 and rc_kmer1 != kmer1:
            seq2_positions.extend(kmer_positions_seq2[rc_kmer1])
        
        # Add all position pairs
        for i in positions1:
            for j in seq2_positions:
                shared_positions.append((i, j))

    return shared_positions



# Examples
# -----------
seq1, seq2 = "AAACTCATC", "TTTCAAATC"
shared_positions = SharedKmersPositions(seq1, seq2, k=3)
for pair in shared_positions:
    print(pair)
# Output: (0, 4) (0, 0) (4, 2) (6, 6)



with open("dataset_30164_5.txt", "r") as file:
    k = int(file.readline().strip())
    seq1 = file.readline().strip()
    seq2 = file.readline().strip()

shared_positions = SharedKmersPositions(seq1, seq2, k)

with open("output.txt", "w") as out_file:
    for pair in shared_positions:
        out_file.write(f"{pair}\n")
# Output: (182, 79) (182, 2111) (7851, 79) (7851, 2111) (10020, 79) ...

 


# How many shared 30-mers do the E. coli and S. enterica genomes have?
# ---------------------------------------------------------------------
with open("E_coli.txt", "r") as file:
     E_coli = file.readline().strip()

with open("Salmonella_enterica.txt", "r") as file:
    S_enterica = file.readline().strip()

shared_positions = SharedKmersPositions(E_coli, S_enterica, k=30)
print(len(shared_positions)) # Output: 268101


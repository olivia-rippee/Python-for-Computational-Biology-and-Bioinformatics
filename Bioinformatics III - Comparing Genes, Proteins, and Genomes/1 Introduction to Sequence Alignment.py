import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics III - Comparing Genes, Proteins, and Genomes/Data")


# -----------------------------------------------
# Manhattan Tourist Problem
# -----------------------------------------------


# How many different paths are there from source to 
# sink in a 16 × 12 rectangular grid?
# -----------------------------------------------
import math
print(math.comb(28, 16))  # or math.comb(28, 12)
# Output: 30421755



# Imagine a hypothetical world in which there are two amino acids, 
# X and Z, having respective masses 2 and 3. How many linear peptides 
# can be formed from these amino acids having mass equal to 25? 
# (Remember that the order of amino acids matters.)
# -----------------------------------------------

def CountPeptides(Mass):
    '''{X: 2, Z: 3}'''
    if Mass == 0:
        return 1
    if Mass < 0:
        return 0
    return CountPeptides(Mass - 2) + CountPeptides(Mass - 3)

total_mass = 25
print(CountPeptides(total_mass))
# Output: 465



# -----------------------------------------------
# Changing Money via Recursion
# -----------------------------------------------

def RecursiveChange(Money, Coins):
    '''Completely impractical because it recalculates the optimal 
    coin combination for a given value of money over and over again.'''
    
    if Money == 0:
        return 0
    MinNumCoins = float('inf')
    for coin in Coins:
        if Money >= coin:
            NumCoins = RecursiveChange(Money - coin, Coins)
            if NumCoins + 1 < MinNumCoins:
                MinNumCoins = NumCoins + 1
    return MinNumCoins


Money = 76
Coins = (5, 4, 1)
# print(RecursiveChange(Money, Coins))



# -----------------------------------------------
# Use dynamic programming to fill in the next ten 
# values of MinNumCoins(m) (i.e., for 13 ≤ m ≤ 22).
# -----------------------------------------------

def DynamicChangeList(Money, Coins):
    '''
    Input: An integer money and an array Coins = (coin1, ..., coind).
    Output: The minimum number of coins with denominations Coins that changes money.
    '''
  
    MinNumCoins = [0] + [float('inf')] * Money
    for m in range(1, Money + 1):
        for coin in Coins:
            if m >= coin:
                MinNumCoins[m] = min(MinNumCoins[m], MinNumCoins[m - coin] + 1)
    return MinNumCoins


coins = [1, 4, 5]
min_num_coins_1to13 = [0, 1, 2, 3, 1, 1, 2, 3, 2, 2, 2, 3, 
                 3, 4, 4, 3, 4, 3, 3, 3, 4, 4, 4]

min_coins = DynamicChangeList(22, coins)
print(*min_coins[1:13]) # matches above
print(*min_coins[13:23])
# Output: 3 3 3 4 4 4 4 4 5 5



# -----------------------------------------------
# Changing Money via Dynamic Programming
# -----------------------------------------------

def DynamicChange(Money, Coins):
    '''
    Input: An integer money and an array Coins = (coin1, ..., coind).
    Output: The minimum number of coins with denominations Coins that changes money,
    where key is the denomination and value is the number of coins.
    '''
    MinNumCoins = [0]  

    for m in range(1, Money + 1): 
        MinNumCoins.append(float('inf'))

        for i in range(len(Coins)): 
            if m >= Coins[i]:
                if MinNumCoins[m - Coins[i]] + 1 < MinNumCoins[m]:  
                    MinNumCoins[m] = MinNumCoins[m - Coins[i]] + 1 

    return MinNumCoins[Money]



# Examples
# --------------
Money = 40
Coins = [50, 25, 20, 10, 5, 1]
print(DynamicChange(Money, Coins))
# Output: 2


with open("dataset_30195_10.txt", "r") as file:
    Money = int(file.readline())
    Coins = list(map(int, file.readline().strip().split()))
print(DynamicChange(Money, Coins))
# Output: 824



# -----------------------------------------------
# Optimized DynamicChange
# -----------------------------------------------

from collections import Counter

def DynamicChange(Money, Coins):
    '''Optimized such that array size required does not exceed the 
    value of the largest coin denomination and that it not only computes 
    the minimum number of coins but also returns these coins.
    
    Input: Integer Money and list Coins.
    Output: (MinNumCoins, dictionary of coins used to form Money)
    '''
    
    MinNumCoins = {0: 0}
    LastCoinUsed = {}

    for m in range(1, Money + 1):
        MinNumCoins[m] = float('inf')
        for coin in Coins:
            if m >= coin:
                if MinNumCoins[m - coin] + 1 < MinNumCoins[m]:
                    MinNumCoins[m] = MinNumCoins[m - coin] + 1
                    LastCoinUsed[m] = coin

    coin_counts = Counter()
    current = Money
    while current > 0:
        coin = LastCoinUsed.get(current)
        if coin is None:
            return float('inf'), {}  # No solution
        coin_counts[coin] += 1
        current -= coin

    return MinNumCoins[Money], dict(coin_counts)


# Example
# ----------
with open("dataset_30195_10.txt", "r") as file:
    Money = int(file.readline())
    Coins = list(map(int, file.readline().strip().split()))
print(*DynamicChange(Money, Coins))
# Output: 824 {22: 818, 21: 5, 12: 1}



# -----------------------------------------------
# Length of Longest Path (Manhattan Tourist)
# -----------------------------------------------

def ManhattanTourist(n, m, Down, Right):
    '''Find the length of a longest path in the Manhattan Tourist Problem.

    Input: Integers n and m, followed by an n × (m + 1) matrix Down and an 
    (n + 1) × m matrix Right. The two matrices are separated by the "-" symbol.
    Output: The length of a longest path from source (0, 0) to sink (n, m) in 
    the rectangular grid whose edges are defined by the matrices Down and Right.
    '''
    
    # Initialize the score matrix s
    s = [[0] * (m + 1) for _ in range(n + 1)]

    # Fill in the first column
    for i in range(1, n + 1):
        s[i][0] = s[i - 1][0] + Down[i - 1][0]

    # Fill in the first row
    for j in range(1, m + 1):
        s[0][j] = s[0][j - 1] + Right[0][j - 1]

    # Fill in the rest of the grid
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i][j] = max(
                s[i - 1][j] + Down[i - 1][j],
                s[i][j - 1] + Right[i][j - 1]
            )

    return s[n][m]


with open("dataset_30205_10.txt", "r") as file:
    lines = [line.strip() for line in file.readlines() if line.strip()]
    n, m = map(int, lines[0].split())
    sep_index = lines.index('-')
    Down = [list(map(int, lines[i].split())) for i in range(1, sep_index)]
    Right = [list(map(int, lines[i].split())) for i in range(sep_index + 1, 
                                                             sep_index + 1 + n + 1)]
length_longest_path = ManhattanTourist(n, m, Down, Right)
print(length_longest_path)
# Output: 82



# -----------------------------------------------
# Length of Longest Path (Manhattan Tourist) Including Diagonals
# -----------------------------------------------

def GetGridDimensions(Down, Right):
    '''Calculate the dimnesions n and m given weight 
    rectangular matrices Down and Right.
    '''
    # number of rows in Down (vertical edges)
    n = len(Down)
    
    # number of columns in Right (horizontal edges)
    m = len(Right[0]) if Right else 0

    return n, m


def ManhattanTourist(n, m, Down, Right, Diagonal):
    '''Find the length of a longest path in the Manhattan Tourist 
    Problem with diagonal weights.

    Input:
    - n, m: grid dimensions
    - Down: n × (m+1) matrix of downward edge weights
    - Right: (n+1) × m matrix of rightward edge weights
    - Diagonal: n × m matrix of diagonal edge weights

    Output:
    - length of longest path from (0,0) to (n,m) considering down, 
    right, and diagonal moves
    '''

    # Initialize the score matrix s
    s = [[0] * (m + 1) for _ in range(n + 1)]

    # Fill in the first column (only downward moves possible)
    for i in range(1, n + 1):
        s[i][0] = s[i - 1][0] + Down[i - 1][0]

    # Fill in the first row (only rightward moves possible)
    for j in range(1, m + 1):
        s[0][j] = s[0][j - 1] + Right[0][j - 1]

    # Fill in the rest of the grid considering diagonal moves
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i][j] = max(
                s[i - 1][j] + Down[i - 1][j],            # from above
                s[i][j - 1] + Right[i][j - 1],           # from left
                s[i - 1][j - 1] + Diagonal[i - 1][j - 1] # from diagonal
            )

    return s[n][m]


# Example
# ------------
Down = [[1, 0, 2, 4, 3],
        [4, 6, 5, 2, 1],
        [4, 4, 5, 2, 1],
        [5, 6, 8, 5, 3]]

Right = [[3, 2, 4, 0],
         [3, 2, 4, 2],
         [0, 7, 3, 4],
         [3, 3, 0, 2],
         [1, 3, 2, 2]]

Diagonal = [[5, 0, 2, 1],
            [8, 4, 3, 0],
            [10, 8, 9, 5],
            [5, 6, 4, 7]]

n, m = GetGridDimensions(Down, Right)
length_longest_path = ManhattanTourist(n, m, Down, Right, Diagonal)
print(length_longest_path)
# Output: 35



# -----------------------------------------------
# Count Possible Topological Orderings
# -----------------------------------------------

from collections import defaultdict

class DirectedGraph:
    def __init__(self, adjacency_list):
        self.graph = defaultdict(list)

        for vertex, neighbors in adjacency_list.items():
            self.graph[vertex] = neighbors

    def WithoutVertex(self, vertex):
        new_graph = self.graph.copy()
        del new_graph[vertex]

        for neighbors in new_graph.values():
            if vertex in neighbors:
                neighbors.remove(vertex)

        return DirectedGraph(new_graph)


    def Sources(self):
        all_vertices = set(self.graph.keys())
        all_neighbors = set(neighbor for neighbors in self.graph.values() for neighbor in neighbors)
        
        return all_vertices - all_neighbors

 
    def NumberOfTopologicalSortings(self):
        if not self.graph:
            return 1

        else:
            return sum(self.WithoutVertex(source).NumberOfTopologicalSortings()
                       for source in self.Sources())
    
    
    def AllTopologicalSortings(self):
        if not self.graph:
            return [[]]  

        all_sortings = []
        for source in self.Sources():
            smaller_graph = self.WithoutVertex(source)
            for order in smaller_graph.AllTopologicalSortings():
                all_sortings.append([source] + order)

        return all_sortings



# Batman Example 
# ----------------
graph = {'tights': ['leotard', 'boots'],
         'leotard': ['shorts','cape','gloves'],
         'shorts': ['boots','belt'],
         'boots': [],
         'cape': ['hood'],
         'gloves': [],
         'belt':[],
         'hood':[]}

dg = DirectedGraph(graph)
dg.NumberOfTopologicalSortings()
# Output: 120



graph = {'a': ['b', 'c', 'e', 'f'],
         'b': ['c', 'f'],
         'c': ['d'],
         'd': [],
         'e': ['d', 'f'],
         'f': []}

dg = DirectedGraph(graph)
all_sortings = dg.AllTopologicalSortings()

for order in all_sortings:
    print(", ".join(order))
    
# Output:
# a, b, c, e, d, f
# a, b, c, e, f, d
# a, b, e, c, d, f
# a, b, e, c, f, d
# a, b, e, f, c, d
# a, e, b, c, d, f
# a, e, b, c, f, d
# a, e, b, f, c, d
    


# -----------------------------------------------
# Longest Common Subsequence
# -----------------------------------------------

def LCSBackTrack(v, w):
    '''Backtrack instead of going forward to avoid getting stuck in nodes.
    '''
    n, m = len(v), len(w)
    s = [[0] * (m + 1) for _ in range(n + 1)]
    backtrack = [[''] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        s[i][0] = 0
    for j in range(m + 1):
        s[0][j] = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = 1 if v[i - 1] == w[j - 1] else 0
            choices = [
                s[i - 1][j],                # down
                s[i][j - 1],                # right
                s[i - 1][j - 1] + match     # diagonal
            ]
            s[i][j] = max(choices)

            if s[i][j] == s[i - 1][j]:
                backtrack[i][j] = "↓"
            elif s[i][j] == s[i][j - 1]:
                backtrack[i][j] = "→"
            elif s[i][j] == s[i - 1][j - 1] + match:
                backtrack[i][j] = "↘"

    return backtrack


def OutputLCS(backtrack, v, i, j):
    '''
    Input: Two sequence strings.
    Output: A longest common subsequence. 
    
    Note: more than one solution may exist.
    '''
    
    if i == 0 or j == 0:
        return ""
    if backtrack[i][j] == "↓":
        return OutputLCS(backtrack, v, i - 1, j)
    elif backtrack[i][j] == "→":
        return OutputLCS(backtrack, v, i, j - 1)
    else: # diagonal ↘
        return OutputLCS(backtrack, v, i - 1, j - 1) + v[i - 1]


# Examples
# ----------
v = "AACCTTGG"
w = "ACACTGTGA"
backtrack = LCSBackTrack(v, w)
lcs = OutputLCS(backtrack, v, len(v), len(w))
print(lcs)
# Output: AACTGG


v = "CTCGAT"
w = "TACGTC"
backtrack = LCSBackTrack(v, w)
lcs = OutputLCS(backtrack, v, len(v), len(w))
print(*lcs)
# Output: TCGT


with open("dataset_30197_5.txt", "r") as file:
    v = file.readline()
    w = file.readline()

backtrack = LCSBackTrack(v, w)
lcs = OutputLCS(backtrack, v, len(v), len(w))
print(lcs)
# Output: GAGTAACTTAAACGTATTTATGAGTATTAAAATGGTTTGGGCCAATCACAGGTCAATAGTGGTCCAAC
# CTGTATTACAGGTCCGGCGGCTTCGGTTGCTAGCCAATCATGTTAGGATCCAGGAATGCTCGGGCTTGGATTAGAG
# GCCTCTAAGACACATGGGTCCCACTGTTTACTCTAACGGCGTCATCTGTTCAGAGCATCCGGCCACTCCTCTTGGG
# GTGTAGTAAGCTCTGGTGCCAGCGTGAATTCTCTATGGAAAGGTCCTTCGCCAAAATAGCATATCGGCGTAATACG
# CTTGAGAGGTTCGTTGTTATACTCGTTTCATGTCCTTGCGGCAACAAGGCACAATGTTGTAAGTTGATCGCACCGC
# AGAGAACAAGCTACTCGTTACTCTAGGAAGCTCTGTTTACCGCATGCAATGAAATACTCGCTATATACTCTCGTAG
# CCATTGGACCCCCTTTTTTTTCCGTTCCGTTTATTTTTAGGCGAGGCGTGAAGCACAATGCGATTGGACGTAGCAG
# GCTGCCCCGCCGGCAGCTATGCTGGATTAATCATGGTAA




# -----------------------------------------------
# Longest Path in an Arbitrary Directed Acyclic Graph
# -----------------------------------------------

from collections import defaultdict

def LongestPathInDAG(start, end, edges):
    '''Find the longest path in a Directed Acyclic Graph.

    Input: An integer representing the starting node to consider in a graph, followed by an integer representing the ending node to consider, followed by a list of edges in the graph. The edge notation "0 1 7" indicates that an edge connects node 0 to node 1 with weight 7.  You may assume a given topological order corresponding to nodes in increasing order.
    Output: The length of a longest path in the graph, followed by a longest path as a sequence of space-separated node labels. (If multiple longest paths exist, you may return any one.)
    '''
    
    adj = defaultdict(list)
    nodes = set()
    for u, v, w in edges:
        adj[u].append((v, w))
        nodes.update([u, v])

    maxNode = max(nodes)
    dist = [-float('inf')] * (maxNode + 1)
    prev = [None] * (maxNode + 1)
    dist[start] = 0

    for u in range(maxNode + 1):  # Given topological order
        for v, w in adj[u]:
            if dist[u] + w > dist[v]:
                dist[v] = dist[u] + w
                prev[v] = u

    if dist[end] == -float('inf'):
        print("No path found")
        return

    # Reconstruct path
    path = []
    curr = end
    while curr is not None:
        path.append(curr)
        curr = prev[curr]
    path.reverse()

    return dist[end], path


# Examples
# -----------
inputData = '''0 4
0 1 7
0 2 4
2 3 2
1 4 1
3 4 3'''

lines = inputData.strip().split('\n')
start, end = map(int, lines[0].split())
edges = [tuple(map(int, line.split())) for line in lines[1:]]

length, path = LongestPathInDAG(start, end, edges)
print(length)
print(' '.join(map(str, path)))
# Output:
# 9
# 0 2 3 4



with open("dataset_30197_7.txt", "r") as file:
    lines = file.read().strip().split('\n')
    start, end = map(int, lines[0].split())
    edges = [tuple(map(int, line.split())) for line in lines[1:]]

length, path = LongestPathInDAG(start, end, edges)

with open("output.txt", "w") as out_file:
    out_file.write(str(length) + "\n")
    out_file.write(' '.join(map(str, path)) + "\n")
    
# Output: 
# 194
# 0 1 3 4 6 9 10 12 14 24 33 38 41 42 47 48 49



# ----------------------------------------------
# Longest Path Given Graph
# ----------------------------------------------

def LongestPathDAG(nodes, edges, start):
    dist = {node: float('-inf') for node in nodes}
    pred = {node: None for node in nodes}
    dist[start] = 0

    for u in nodes:
        for v, w in edges.get(u, []):
            if dist[u] + w > dist[v]:
                dist[v] = dist[u] + w
                pred[v] = u

    # Find the node with the maximum distance
    max_node = max(dist, key=dist.get)

    # Reconstruct the longest path
    path = []
    current = max_node
    while current is not None:
        path.append(current)
        current = pred[current]

    path.reverse()
    return dist[max_node], path




# Example
# ------------
nodes = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
start = nodes[0]
edges = {'a': [('b', 5), ('c', 6), ('d', 5)],
         'b': [('c', 2), ('f', 4)],
         'c': [('e', 4), ('f', 3), ('g', 5)],
         'd': [('e', 6), ('f', 8)],
         'e': [('g', 2)],
         'f': [('g', 1)],
         'g': []}

length, path = LongestPathDAG(nodes, edges, start)
print(' '.join(map(str, path)))
# Output: a d f g




# -----------------------------------------------
# Every LCS of Two Strings
# -----------------------------------------------

def LCSBackTrack(v, w):
    '''Find all longest common substrings between two strings.'''
    
    n, m = len(v), len(w)
    s = [[0] * (m + 1) for _ in range(n + 1)]
    backtrack = [[[] for _ in range(m + 1)] for _ in range(n + 1)]  # list of directions

    for i in range(n + 1):
        s[i][0] = 0
    for j in range(m + 1):
        s[0][j] = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = 1 if v[i - 1] == w[j - 1] else 0
            choices = [
                s[i - 1][j],                # down
                s[i][j - 1],                # right
                s[i - 1][j - 1] + match     # diagonal
            ]
            max_score = max(choices)
            s[i][j] = max_score

            # Store all directions that yield max_score
            if max_score == s[i - 1][j]:
                backtrack[i][j].append("↓")
            if max_score == s[i][j - 1]:
                backtrack[i][j].append("→")
            if max_score == s[i - 1][j - 1] + match:
                backtrack[i][j].append("↘")

    return backtrack


def OutputLCSSet(backtrack, v, i, j):
    '''Return the set all longest common substrings between two strings.'''
    
    if i == 0 or j == 0:
        return {""}  # return a set with empty string to allow combination

    lcs_set = set()

    for direction in backtrack[i][j]:
        if direction == "↓":
            lcs_set.update(OutputLCS(backtrack, v, i - 1, j))
        elif direction == "→":
            lcs_set.update(OutputLCS(backtrack, v, i, j - 1))
        else:  # diagonal ↘
            prev_set = OutputLCS(backtrack, v, i - 1, j - 1)
            for seq in prev_set:
                lcs_set.add(seq + v[i - 1])

    return lcs_set


# Example
# --------
v = "AGGGGCTC"
w = "ACTGGTCA"

backtrack = LCSBackTrack(v, w)
print(*OutputLCSSet(backtrack, v, len(v), len(w)))
# Output: AGGGGTC AGGTC AGGGTC

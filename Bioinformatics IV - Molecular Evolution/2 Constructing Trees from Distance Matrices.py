import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics IV - Molecular Evolution/Data")

# -----------------------------------------------
# Calculate Discrepancy (SSE) 
# -----------------------------------------------

def Discrepancy(T, D):
    '''Calculate Discrepancy (SSE) between non-additive matrix (tree_dist) D and 
    approximation T (matrix_dist).
    Discrepancy(T, D) = sum_{1<=i<=j<=n} (d_T[i][j] - D[i][j])**2
    
    Input:
        T (list of list of int): Distance matrix from the tree (n x n)
        D (list of list of int): Given distance matrix (n x n)

    Output: The squared discrepancy.
    '''
    
    n = len(D)
    discrepancy = 0
    for i in range(n):
        for j in range(i + 1, n):
            diff = D[i][j] - T[i][j]
            discrepancy += diff * diff
    return discrepancy


# Example 
# ---------
T = [[0, 7, 9, 10],
     [7, 0, 10, 11],
     [9, 10, 0, 3],
     [10, 11, 3, 0]]

D = [[0, 3, 4, 3],
     [3, 0, 4, 5],
     [4, 4, 0, 2],
     [3, 5, 2, 0]]

result = Discrepancy(T, D)
# Output: 163



# -----------------------------------------------
# UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
# -----------------------------------------------
import sys

def UPGMA(D):
    '''Find the adjacency list for a non-additive distance matrix using an
    Unweighted Pair Group Method with Arithmetic Mean (UPGMA) algorithm.
    Note: Generates incorrect trees from additive matrices.
    
    Input: Initial distance matrix D (n x n).
    Output: Adjacency list with edge lengths.
    '''
    n = len(D)
    # Each cluster: (node_id, members set, height)
    clusters = {i: {'members': {i}, 'height': 0.0} for i in range(n)}
    # Distance between clusters (keys: frozenset of cluster ids)
    dist = {}
    for i in range(n):
        for j in range(i+1, n):
            dist[frozenset({i,j})] = D[i][j]

    next_node = n
    # adjacency list: node -> list of (neighbor, weight)
    adj = {i: [] for i in range(2*n-1)}

    while len(clusters) > 1:
        pair, dmin = min(dist.items(), key=lambda x: x[1])
        c1, c2 = tuple(pair)
    
        members = clusters[c1]['members'].union(clusters[c2]['members'])
        height = dmin / 2.0
    
        clusters[next_node] = {'members': members, 'height': height}
    
        h1 = clusters[c1]['height']
        h2 = clusters[c2]['height']
        w1 = height - h1
        w2 = height - h2
    
        adj[next_node].append((c1, w1))
        adj[next_node].append((c2, w2))
        adj[c1].append((next_node, w1))
        adj[c2].append((next_node, w2))
    
        # Calculate distances BEFORE removing old distances
        for c_other in clusters:
            if c_other != c1 and c_other != c2 and c_other != next_node:
                size1 = len(clusters[c1]['members'])
                size2 = len(clusters[c2]['members'])
                d1 = dist[frozenset({c1, c_other})]
                d2 = dist[frozenset({c2, c_other})]
                new_d = (size1 * d1 + size2 * d2) / (size1 + size2)
                dist[frozenset({next_node, c_other})] = new_d
    
        # Now remove old distances involving c1, c2
        keys_to_remove = [k for k in dist if c1 in k or c2 in k]
        for k in keys_to_remove:
            dist.pop(k)
    
        # Remove merged clusters
        clusters.pop(c1)
        clusters.pop(c2)
    
        next_node += 1

    # Output adjacency list
    # We only output edges from each node to neighbors with weight rounded to 3 decimals
    # Sort output by node number, then neighbor number for consistent output
    for node in range(next_node):
        adj[node].sort(key=lambda x: x[0])

    lines = []
    for node in range(next_node):
        for (nbr, w) in adj[node]:
            # print edge once in each direction, so just print all
            # format weight to 3 decimals as in sample
            lines.append(f"{node}->{nbr}:{w:.3f}")

    return lines


def ReadDistanceMatrix(filename):
    with open(filename, "r") as f:
        n = int(f.readline().strip())
        D = []
        for _ in range(n):
            row = list(map(float, f.readline().strip().split()))
            D.append(row)
    return D

def ReadDistanceMatrixNamed(filename):
    with open(filename, "r") as f:
        labels = f.readline().strip().split('\t')[1:]  # skip first empty cell
        D = []
        for line in f:
            parts = line.strip().split('\t')
            # Skip the first column (row label)
            row_values = list(map(float, parts[1:]))
            D.append(row_values)
    return labels, D


# Example
# ---------
D = [[0, 20, 17, 11],
     [20, 0, 20, 13],
     [17, 20, 0, 10],
     [11, 13, 10, 0]]
tree = UPGMA(D)
with open("output.txt", "w") as file:
    for line in tree:
        file.write(line + "\n")

# Output: 0->5:7.000 1->6:8.833 2->4:5.000 3->4:5.000 4->2:5.000 4->3:5.000 
# 4->5:2.000 5->0:7.000 5->4:2.000 5->6:1.833 6->1:8.833 6->5:1.833



D = ReadDistanceMatrix("dataset_30288_8.txt")
tree = UPGMA(D)
with open("output.txt", "w") as file:
    for line in tree:
        file.write(line + "\n")
        
# Output: 0->30:324.500 1->28:322.500 2->33:345.500 3->26:314.500 
# 4->28:322.500 5->48:509.521 6->39:393.250 7->29:323.500 8->25:313.500



labels, D = ReadDistanceMatrixNamed("coronavirus_distance_matrix_nonadditive.txt")
tree = UPGMA(D)
with open("output.txt", "w") as file:
    for line in tree:
        file.write(line + "\n")
        
# Output: 0->10:147.500 1->10:147.500 2->11:153.500 3->12:247.167 4->13:409.000
# 5->13:409.000 6->15:496.500 7->9:8.000 8->9:8.000 9->7:8.000 9->8:8.000
# 9->14:470.625 10->0:147.500 10->1:147.500 10->11:6.000 11->2:153.500
# 11->10:6.000 11->12:93.667 12->3:247.167 12->11:93.667 12->14:231.458
# 13->4:409.000 13->5:409.000 13->16:129.857 14->9:470.625 14->12:231.458
# 14->15:17.875 15->6:496.500 15->14:17.875 15->16:42.357 16->13:129.857 16->15:42.357



# -----------------------------------------------
# Neighbor-Joining Algorithm
# -----------------------------------------------
import numpy as np

def TotalDistance(D, idx):
    """Sum of distances from node idx to all other nodes."""
    return sum(D[idx])

def NeighborJoiningMatrix(D):
    """Construct the neighbor-joining matrix D* from distance matrix D."""
    n = len(D)
    total_dist = [TotalDistance(D, i) for i in range(n)]
    D_star = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                D_star[i, j] = (n - 2) * D[i][j] - total_dist[i] - total_dist[j]
            else:
                D_star[i, j] = np.inf  # ignore diagonal
    return D_star

def NeighborJoining(D, labels=None, next_node=None):
    '''
    Input: An integer n, followed by an n x n distance matrix.
    Output: An adjacency list for the tree resulting from applying the neighbor-joining algorithm. Edge-weights should be accurate to two decimal places (they are provided to three decimal places in the sample output below).
    '''

    n = len(D)
    if labels is None:
        labels = list(range(n))  # initial leaf labels 0..n-1
    if next_node is None:
        next_node = n  # next internal node label

    if n == 2:
        # Base case: tree with single edge
        # Return edges in both directions for output format
        edges = {}
        i, j = labels[0], labels[1]
        length = D[0][1]
        edges[(i, j)] = length
        edges[(j, i)] = length
        return edges, next_node

    D_star = NeighborJoiningMatrix(D)
    
    # Find pair (i,j) with minimal D_star[i,j]
    i, j = np.unravel_index(np.argmin(D_star), D_star.shape)
    
    total_dist = [TotalDistance(D, k) for k in range(n)]
    delta = (total_dist[i] - total_dist[j]) / (n - 2)
    
    limb_length_i = 0.5 * (D[i][j] + delta)
    limb_length_j = 0.5 * (D[i][j] - delta)
    
    # Create new node m
    m = next_node
    next_node += 1
    
    # Compute distances from other nodes k to m
    D_m = []
    nodes = [x for x in range(n) if x != i and x != j]
    for k in nodes:
        dist = 0.5 * (D[k][i] + D[k][j] - D[i][j])
        D_m.append(dist)
    
    # Build new distance matrix D_new (size n-1)
    D_new = np.zeros((n-1, n-1))
    new_labels = [labels[k] for k in nodes] + [m]
    
    # Fill D_new
    for a in range(n-2):
        for b in range(n-2):
            D_new[a][b] = D[nodes[a]][nodes[b]]
    for a in range(n-2):
        D_new[a][n-2] = D_m[a]
        D_new[n-2][a] = D_m[a]
    D_new[n-2][n-2] = 0.0
    
    # Recursive call
    T, next_node = NeighborJoining(D_new.tolist(), new_labels, next_node)
    
    # Add limbs connecting i and j to new node m
    leaf_i, leaf_j = labels[i], labels[j]
    # Add edges both ways for output
    T[(m, leaf_i)] = limb_length_i
    T[(leaf_i, m)] = limb_length_i
    T[(m, leaf_j)] = limb_length_j
    T[(leaf_j, m)] = limb_length_j
    
    return T, next_node



# Examples
# ---------
n = 4
D = [[0, 23, 27, 20],
     [23, 0, 30, 28],
     [27, 30, 0, 30],
     [20, 28, 30, 0]]
T, _ = NeighborJoining(D)
for (a, b), length in sorted(T.items()):
    print(f"{a}->{b}:{length:.3f}".rstrip('0').rstrip('.'))

# Output: 0->4:8 1->5:13.5 2->5:16.5 3->4:12 4->0:8 
# 4->3:12 4->5:2 5->1:13.5 5->2:16.5 5->4:2




with open("dataset_30289_7.txt", "r") as file:
    n = int(file.readline())
    D = []
    for _ in range(n):
        row = file.readline().strip().split()
        D.append([float(x) for x in row])

T, _ = NeighborJoining(D)

with open("output.txt", "w") as file:
    for (a, b), length in sorted(T.items()):
        line = f"{a}->{b}:{length:.3f}".rstrip('0').rstrip('.')
        file.write(line + '\n')
        
# Output: 0->38:537.115 1->37:547.94 2->39:528.261 3->37:501.06 4->43:557.132
# 5->42:562.775 6->50:484.953 7->41:537.155 8->46:595.625 9->40:478.097 ...




with open("coronavirus_distance_matrix_nonadditive.txt", "r") as file:
    header = file.readline().strip().split('\t')
    labels = header[1:]  # skip the empty first column header
    D = []
    for _ in range(len(labels)):
        row = file.readline().strip().split('\t')[1:]  # skip first column label
        D.append([float(x) for x in row])
    
T, _ = NeighborJoining(D, labels=list(range(len(labels))))
    
def LabelName(x):
    if x < len(labels):
        return labels[x]
    else:
        return str(x)
    
with open("output.txt", "w") as file:
    for (a, b), length in sorted(T.items()):
        line = f"{LabelName(a)}->{LabelName(b)}:{length:.3f}".rstrip('0').rstrip('.')
        file.write(line + '\n')

# Output: Pig->12:146.219 Horse->12:148.781 Mouse->13:149.156 Dog->11:252.979 Cat->8:401.083
# Turkey->8:416.917 Civet->9:488.8 Human->10:463.719 8->Cat:401.083 8->Turkey:416.917
# 8->9:163.7 9->Civet:488.8 9->8:163.7 9->10:22.531 10->Human:463.719 10->9:22.531
# 10->11:249.771 11->Dog:252.979 11->10:249.771 11->13:86.406 12->Pig:146.219 
# 12->Horse:148.781 12->13:10.344 13->Mouse:149.156 13->11:86.406 13->12:10.344

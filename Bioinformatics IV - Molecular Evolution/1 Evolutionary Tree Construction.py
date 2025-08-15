import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics IV - Molecular Evolution/Data")

# -----------------------------------------------
# Distances Between Leaves
# -----------------------------------------------

from collections import defaultdict
import heapq

def ParseInputToAdjacencyList(inputLines):
    adjacency = defaultdict(list)
    for line in inputLines:
        src, rest = line.split("->")
        dst, weight = rest.split(":")
        src, dst, weight = int(src), int(dst), int(weight)
        adjacency[src].append((dst, weight))
    return adjacency

def Dijkstra(startNode, adjacency, n):
    distances = {}
    heap = [(0, startNode)]
    while heap:
        currDist, node = heapq.heappop(heap)
        if node in distances:
            continue
        distances[node] = currDist
        for neighbor, weight in adjacency[node]:
            if neighbor not in distances:
                heapq.heappush(heap, (currDist + weight, neighbor))
    return [distances[i] for i in range(n)]

def ComputeLeafDistances(n, inputLines):
    '''
    Input:  An integer n followed by the adjacency list of a weighted tree with n leaves.
    Output: An n x n matrix (di,j), where di,j is the length of the path between leaves i and j.'''
    
    adjacency = ParseInputToAdjacencyList(inputLines)
    resultMatrix = []
    for i in range(n):
        distances = Dijkstra(i, adjacency, n)
        resultMatrix.append(distances)
    return resultMatrix

def PrintDistanceMatrix(matrix):
    for row in matrix:
        print("\t".join(map(str, row)))

def ReadInputFromFile(filePath):
    with open(filePath, 'r') as file:
        lines = [line.strip() for line in file.readlines()]
    n = int(lines[0])
    edges = lines[1:]
    return n, edges

def WriteOutputToFile(filePath, matrix):
    with open(filePath, 'w') as file:
        for row in matrix:
            file.write("\t".join(map(str, row)) + "\n")



# Example 1
# -----------
inputData = ["0->4:11", "1->4:2", "2->5:6", "3->5:7", "4->0:11", 
             "4->1:2", "4->5:4", "5->4:4", "5->3:7", "5->2:6"]

n = 4
matrix = ComputeLeafDistances(n, inputData)
PrintDistanceMatrix(matrix)

# Output:    
    # 0	    13	21	22
    # 13	0	12	13
    # 21	12	0	13
    # 22	13	13	0


# Example 2
# -----------
n, inputLines = ReadInputFromFile("dataset_30284_12.txt")
matrix = ComputeLeafDistances(n, inputLines)
WriteOutputToFile("output.txt", matrix)

# Output: 
    # 0	31 52 100 29 208 89	67	37	54	124	103	66	76	104	169	178	85	214	150	
        # 95 129 31	68	203	108	187	187	50	51	139	145
    # 31 0 45 93 32	211	92	70	40	57	127	106	59	79	97	172	181	78	217 153 
        # 88 132 24	61	206	111	190	190	53	54	142	148
    # ...


# -----------------------------------------------
# Limb Length
# -----------------------------------------------

def LimbLength(n, j, D):
    ''''Compute the length of a limb in a tree defined by an additive distance matrix 
    (assuming i and j are neighbors).
    Input: An integer n (n = len(D)), an integer j between 0 and n - 1, and a space-separated 
    additive distance matrix D (whose elements are integers).
    Output: The limb length of the leaf in Tree(D) corresponding to row j of this distance 
    matrix (use 0-based indexing).
    '''
    
    minLimbLength = float('inf')
    for i in range(n):
        if i == j:
            continue
        for k in range(n):
            if k == j or i == k:
                continue
            limbLength = (D[i][j] + D[j][k] - D[i][k]) // 2
            if limbLength < minLimbLength:
                minLimbLength = limbLength
    return minLimbLength

def ReadFromFile(filePath):
    with open(filePath, "r") as file:
        lines = [line.strip() for line in file if line.strip()]
    n = int(lines[0])
    j = int(lines[1])
    matrix = [list(map(int, line.split())) for line in lines[2:]]
    return n, j, matrix

def WriteToFile(filePath, value):
    with open(filePath, "w") as file:
        file.write(str(value) + '\n')
        
        
# Example 1
# -----------
n, j = 4, 1
D = [[0, 13, 21, 22],
     [13, 0, 12, 13],
     [21, 12, 0, 13],
     [22, 13, 13, 0]]
print(LimbLength(n, j, D))
# Output: 2



# Example 2
# -----------
n, j, D = ReadFromFile("dataset_30285_11.txt")
print(LimbLength(n, j, D))
# Output: 259


# Example 3
# -----------
D = [[0, 20, 9, 11],
     [20, 0, 17, 11],
     [9, 17, 0, 8],
     [11, 11, 8, 0]]
n = len(D)
k = 2  # Matrix indices: i  j  k  l
print(LimbLength(n, k, D))
# Output: 3


# Example 4
# -----------
D = [[0, 14, 17, 17],  # i
     [14, 0, 7, 13],   # j
     [17, 7, 0, 16],   # k
     [17, 13, 16, 0]]   # l
n = len(D)
k = 2 
print(LimbLength(n, k, D))
# Output: 5


# Example 5
# -----------
D = [[0, 13, 16, 10],  # i
     [13, 0, 21, 15],   # j
     [16, 21, 0, 18],   # k
     [10, 15, 18, 0]]   # l
n = len(D)
i = 0
print(LimbLength(n, i, D))
# Output: 4



# -----------------------------------------------
# Additive Phylogeny (most matrices are not additive)
# -----------------------------------------------

def LimbLength(n, j, D):
    minLimbLength = float('inf')
    for i in range(n):
        if i == j:
            continue
        for k in range(n):
            if k == j or i == k:
                continue
            limbLength = (D[i][j] + D[j][k] - D[i][k]) // 2
            if limbLength < minLimbLength:
                minLimbLength = limbLength
    return minLimbLength


def AdditivePhylogeny(D, n, nextNode):
    '''Returns adjacency list and current nextNode index.
    Inputs:
      - D: distance matrix
      - n: current number of leaves (working size)
      - nextNode: next internal node id to assign (starts at n)
    Outputs:
      - tree: adjacency list {node: [(neighbor,weight), ...]}
      - nextNode: updated next internal node id
    '''

    if n == 2:
        # Base case: only two leaves, connect directly
        tree = {0: [(1, D[0][1])], 1: [(0, D[0][1])]}
        return tree, nextNode

    # Find limb length for leaf n-1
    limbLen = LimbLength(n, n - 1, D)

    # Adjust the distance matrix by subtracting limbLen from last leaf distances
    for j in range(n - 1):
        D[j][n - 1] -= limbLen
        D[n - 1][j] = D[j][n - 1]

    # Find two nodes i,k s.t. D[i][k] = D[i][n-1] + D[n-1][k]
    # i and k in range 0 to n-2 (excluding leaf n-1)
    x = None
    i = k = None
    for a in range(n - 1):
        for b in range(n - 1):
            if a == b:
                continue
            if D[a][b] == D[a][n - 1] + D[n - 1][b]:
                i, k = a, b
                x = D[a][n - 1]
                break
        if x is not None:
            break

    # Remove leaf n-1 from D to form D'
    Dprime = [row[:n-1] for row in D[:n-1]]

    # Recursively build tree for n-1 leaves
    tree, nextNode = AdditivePhylogeny(Dprime, n - 1, nextNode)

    # Add back leaf n-1 by inserting it at distance x along path i-k

    # Find path from i to k in current tree
    path = FindPath(tree, i, k)

    # Traverse path to find where to insert new node at distance x from i
    distFromI = 0
    for idx in range(len(path) - 1):
        u = path[idx]
        v = path[idx + 1]
        w = GetEdgeWeight(tree, u, v)
        if distFromI + w == x:
            # Insert leaf directly connected to v
            AttachLeaf(tree, nextNode, n - 1, v, limbLen)
            nextNode += 1
            break
        elif distFromI + w > x:
            # Split edge (u,v) to insert internal node at distance x
            distToSplit = x - distFromI
            # Insert new internal node
            newNode = nextNode
            nextNode += 1
            # Modify edges:
            # Remove edge u-v
            RemoveEdge(tree, u, v)
            # Add edges u-newNode, newNode-v
            AddEdge(tree, u, newNode, distToSplit)
            AddEdge(tree, newNode, v, w - distToSplit)
            # Attach leaf n-1 to newNode
            AttachLeaf(tree, newNode, n - 1, limbLen)
            break
        distFromI += w

    return tree, nextNode


def FindPath(tree, start, end, path=None, visited=None):
    '''DFS to find any path from start to end'''
    if path is None:
        path = []
    if visited is None:
        visited = set()
    path.append(start)
    visited.add(start)
    if start == end:
        return path[:]
    for neighbor, _ in tree.get(start, []):
        if neighbor not in visited:
            res = FindPath(tree, neighbor, end, path, visited)
            if res is not None:
                return res
    path.pop()
    return None


def GetEdgeWeight(tree, u, v):
    '''Return weight of edge u-v'''
    for neigh, w in tree[u]:
        if neigh == v:
            return w
    return None


def RemoveEdge(tree, u, v):
    '''Remove edge u-v and v-u from tree'''
    tree[u] = [(n, w) for n, w in tree[u] if n != v]
    tree[v] = [(n, w) for n, w in tree[v] if n != u]


def AddEdge(tree, u, v, w):
    '''Add edge u-v with weight w bidirectionally'''
    if u not in tree:
        tree[u] = []
    if v not in tree:
        tree[v] = []
    tree[u].append((v, w))
    tree[v].append((u, w))


def AttachLeaf(tree, attachNode, leaf, limbLen):
    '''Attach leaf to attachNode with edge weight limbLen'''
    AddEdge(tree, attachNode, leaf, limbLen)


def ParseInput():
    n = int(input())
    D = []
    for _ in range(n):
        row = list(map(int, input().split()))
        D.append(row)
    return n, D


def ReadNamedMatrix(filename):
    with open(filename, "r") as file:
        lines = file.readlines()

    # Extract column headers (skip first empty cell)
    headers = lines[0].strip().split('\t')[1:]

    matrix = {}
    for line in lines[1:]:
        parts = line.strip().split('\t')
        row_name = parts[0]
        values = list(map(int, parts[1:]))

        # Map column header to value
        matrix[row_name] = dict(zip(headers, values))

    return matrix

def PrintTree(tree):
    '''Print adjacency list as per problem format'''
    for u in sorted(tree.keys()):
        for v, w in tree[u]:
            print(f"{u}->{v}:{w}")

def ReadInputFromFile(filename):
    with open(filename, "r") as file:
        n = int(file.readline().strip())
        D = []
        for _ in range(n):
            row = list(map(int, file.readline().strip().split()))
            D.append(row)
    return n, D

def WriteOutputToFile(filename, tree):
    with open(filename, "w") as file:
        for u in sorted(tree.keys()):
            for v, length in sorted(tree[u], key=lambda x: x[0]):
                file.write(f"{u}->{v}:{length}\n")


# Example 1
# -----------
n = 4
D = [[0, 13, 21, 22],
    [13, 0, 12, 13],
    [21, 12, 0, 13],
    [22, 13, 13, 0]]

tree, _ = AdditivePhylogeny(D, n, nextNode=n)

for u in sorted(tree.keys()):
    for (v, length) in sorted(tree[u], key=lambda x: x[0]):
        print(f"{u}->{v}:{length}")

# Output:
    # 0->4:11  1->4:2    2->5:6    3->5:0    4->0:11    
    # 4->1:2   4->5:4    5->2:6    5->3:0    5->4:4


# Example 2
# -----------
n, D = ReadInputFromFile("dataset_30286_6.txt")
tree, _ = AdditivePhylogeny(D, n, nextNode=n)
WriteOutputToFile("output.txt", tree)

# Output:
    # 0->29:465 1->51:66 2->49:154 3->30:539 4->31:987 5->32:476 6->33:229
    # 7->34:134 8->35:427 9->36:342 10->37:722 11->38:275 12->39:423 ...


# Example 3
# -----------
D = ReadNamedMatrix("coronavirus_distance_matrix_additive.txt")
n = len(D)
tree, _ = AdditivePhylogeny(D, n, nextNode=n)
WriteOutputToFile("output.txt", tree)

# Output:
    # 0->29:465 1->51:66 2->49:154 3->30:539 4->31:987 ...


# ----------------------------------------------- 
# Check whether matrix is additive
# -----------------------------------------------

def IsAdditiveMatrix(matrix):
    if not matrix or not matrix[0]:
        return True  # Empty matrix is trivially additive
    
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Initialize row vector r and column vector c
    r = [0] * rows
    c = [0] * cols
    
    # Fix r[0] = 0 (arbitrary, to solve system uniquely)
    # Then c[j] = M[0][j]
    for j in range(cols):
        c[j] = matrix[0][j]
    
    # Calculate r[i] = M[i][0] - c[0]
    for i in range(rows):
        r[i] = matrix[i][0] - c[0]
    
    # Check if M[i][j] == r[i] + c[j] for all i,j
    for i in range(rows):
        for j in range(cols):
            if matrix[i][j] != r[i] + c[j]:
                return False
    
    return True


# Example 1
# -----------
matrix1 = [[1, 2, 3],
           [2, 3, 4],
           [3, 4, 5]]
print(IsAdditiveMatrix(matrix1))  # True


# Example 2
# -----------
matrix2 = [[1, 2],
           [3, 5]]
print(IsAdditiveMatrix(matrix2))  # False


# ----------------------------------------------- 
# Determine Matrices Fit By a Tree
# -----------------------------------------------

from collections import defaultdict
import heapq

def BuildTree(edges):
    '''Builds an adjacency list from edge list.'''
    tree = defaultdict(list)
    for u, v, w in edges:
        tree[u].append((v, w))
        tree[v].append((u, w))
    return tree

def Dijkstra(start, tree):
    '''Computes shortest paths from a start node using Dijkstra's algorithm.'''
    dist = {node: float('inf') for node in tree}
    dist[start] = 0
    heap = [(0, start)]
    while heap:
        d, u = heapq.heappop(heap)
        for v, w in tree[u]:
            if dist[v] > d + w:
                dist[v] = d + w
                heapq.heappush(heap, (d + w, v))
    return dist

def ComputeDistanceMatrix(leaves, tree):
    '''Computes the distance matrix between all pairs of leaves.'''
    matrix = {}
    for leaf in leaves:
        dists = Dijkstra(leaf, tree)
        matrix[leaf] = {other: dists[other] for other in leaves}
    return matrix

def MatricesMatch(mat1, mat2, leaves):
    '''Checks if two distance matrices match exactly for a given leaf set.'''
    for a in leaves:
        for b in leaves:
            if mat1[a][b] != mat2[a][b]:
                return False
    return True

def FindMatchingMatrixForTree(edges, leaves, matrices):
    '''Builds the tree, computes the distance matrix, and checks for a match.'''
    tree = BuildTree(edges)
    computedMatrix = ComputeDistanceMatrix(leaves, tree)
    
    for candidate in matrices:
        if MatricesMatch(candidate['matrix'], computedMatrix, leaves):
            print(f"The tree fits {candidate['name']}")
            return candidate['name']
    
    print("No matching matrix found.")
    return None


# Example 1
# -----------
edges = [('i', 'center1', 2),
         ('j', 'center1', 4),
         ('center1', 'center2', 6),
         ('center2', 'k', 1),
         ('center2', 'l', 6)]

leaves = ['i', 'j', 'k', 'l']

matrices = [
    {'name': 'Matrix 1',
     'matrix': {
         'i': {'i': 0, 'j': 6, 'k':10, 'l':14},
         'j': {'i': 6, 'j': 0, 'k':12, 'l':16},
         'k': {'i':10, 'j':12, 'k': 0, 'l': 6},
         'l': {'i':14, 'j':16, 'k': 6, 'l': 0},}},
    {'name': 'Matrix 2',
     'matrix': {
         'i': {'i': 0, 'j': 7, 'k':10, 'l':14},
         'j': {'i': 7, 'j': 0, 'k':11, 'l':15},
         'k': {'i':10, 'j':11, 'k': 0, 'l': 6},
         'l': {'i':14, 'j':15, 'k': 6, 'l': 0},}},
    {'name': 'Matrix 3',
     'matrix': {
         'i': {'i': 0, 'j': 6, 'k': 9, 'l':13},
         'j': {'i': 6, 'j': 0, 'k':11, 'l':15},
         'k': {'i': 9, 'j':11, 'k': 0, 'l': 6},
         'l': {'i':13, 'j':15, 'k': 6, 'l': 0},}},
    {'name': 'Matrix 4',
     'matrix': {
         'i': {'i': 0, 'j': 6, 'k': 9, 'l':14},
         'j': {'i': 6, 'j': 0, 'k':11, 'l':16},
         'k': {'i': 9, 'j':11, 'k': 0, 'l': 7},
         'l': {'i':14, 'j':16, 'k': 7, 'l': 0},}}]

FindMatchingMatrixForTree(edges, leaves, matrices)
# Output: Matrix 4


# Example 2
# -----------
edges = [('i', 'center1', 2),
        ('j', 'center1', 4),
        ('center1', 'center2', 6),
        ('center2', 'l', 5),
        ('center2', 'k', 1)]

leaves = ['i', 'j', 'k', 'l']

matrices = [
    {'name': 'Matrix 1',
     'matrix': {
         'i': {'i': 0, 'j': 6, 'k':10, 'l':13},
         'j': {'i': 6, 'j': 0, 'k':12, 'l':15},
         'k': {'i':10, 'j':12, 'k': 0, 'l': 7},
         'l': {'i':13, 'j':15, 'k': 7, 'l': 0},}},
    
    {'name': 'Matrix 2',
     'matrix': {
         'i': {'i': 0, 'j': 6, 'k': 9, 'l':14},
         'j': {'i': 6, 'j': 0, 'k':11, 'l':16},
         'k': {'i': 9, 'j':11, 'k': 0, 'l': 7},
         'l': {'i':14, 'j':16, 'k': 7, 'l': 0},}},
    
    {'name': 'Matrix 3',
     'matrix': {
         'i': {'i': 0, 'j': 6, 'k':10, 'l':14},
         'j': {'i': 6, 'j': 0, 'k':12, 'l':16},
         'k': {'i':10, 'j':12, 'k': 0, 'l': 6},
         'l': {'i':14, 'j':16, 'k': 6, 'l': 0},}},
    
    {'name': 'Matrix 4',
     'matrix': {
         'i': {'i': 0, 'j': 6, 'k': 9, 'l':13},
         'j': {'i': 6, 'j': 0, 'k':11, 'l':15},
         'k': {'i': 9, 'j':11, 'k': 0, 'l': 6},
         'l': {'i':13, 'j':15, 'k': 6, 'l': 0},}}]

FindMatchingMatrixForTree(edges, leaves, matrices)
# Output: Matrix 4
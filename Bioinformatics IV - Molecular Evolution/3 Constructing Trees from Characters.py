import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics IV - Molecular Evolution/Data")


# -----------------------------------------------
# Small Parsimony for Rooted Trees
# -----------------------------------------------

from collections import defaultdict

def SmallParsimonyRooted(n, edges):
    '''Find the most parsimonious labeling of the internal nodes in a rooted tree.
    
    Input: An integer n followed by an adjacency list for a rooted binary tree with n 
    leaves labeled by DNA strings.
    Output: The minimum parsimony score of this tree, followed by the adjacency list 
    of a tree corresponding to labeling internal nodes by DNA strings in order to 
    minimize the parsimony score of the tree.
    '''
    
    children = defaultdict(list)
    parent = {}
    label = {}

    next_leaf_node = 1000000  # large start for unique leaf IDs

    for line in edges:
        u, v = line.split("->")
        u = int(u)
        # Check if v is leaf label or internal node
        if all(ch in "ACGT" for ch in v):  # leaf label
            leaf_node = next_leaf_node
            next_leaf_node += 1
            label[leaf_node] = v
            children[u].append(leaf_node)
            parent[leaf_node] = u
        else:
            v = int(v)
            children[u].append(v)
            parent[v] = u

    all_nodes = set(children.keys()) | set(parent.keys()) | set(label.keys())

    # Find root (node with no parent)
    root = None
    for node in all_nodes:
        if node not in parent:
            root = node
            break

    leaves = [node for node in label]
    length = len(label[leaves[0]])
    Bases = ['A', 'C', 'G', 'T']

    final_labels = {node: [''] * length for node in all_nodes}
    total_score = 0

    for pos in range(length):
        S = {node: [float('inf')] * 4 for node in all_nodes}
        Tag = {node: 0 for node in all_nodes}

        # Initialize leaves
        for leaf in leaves:
            Tag[leaf] = 1
            for i, b in enumerate(Bases):
                S[leaf][i] = 0 if label[leaf][pos] == b else float('inf')

        # Ripe node loop
        while True:
            ripe_nodes = []
            for node in all_nodes:
                if Tag[node] == 0:
                    if all(Tag[child] == 1 for child in children[node]):
                        ripe_nodes.append(node)
            if not ripe_nodes:
                break

            for v in ripe_nodes:
                for k in range(4):
                    cost_sum = 0
                    for c in children[v]:
                        costs = [S[c][i] + (0 if i == k else 1) for i in range(4)]
                        cost_sum += min(costs)
                    S[v][k] = cost_sum
                Tag[v] = 1

        min_score = min(S[root])
        total_score += min_score

        assigned = {}

        def AssignChar(node, chosen_char_idx):
            assigned[node] = chosen_char_idx
            final_labels[node][pos] = Bases[chosen_char_idx]
            for c in children[node]:
                costs = [S[c][i] + (0 if i == chosen_char_idx else 1) for i in range(4)]
                min_i = costs.index(min(costs))
                AssignChar(c, min_i)

        root_char_idx = S[root].index(min_score)
        AssignChar(root, root_char_idx)

    # Convert labels
    for node in final_labels:
        final_labels[node] = ''.join(final_labels[node])

    def HammingDistance(s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    edges_out = []
    for u in children:
        for v in children[u]:
            dist = HammingDistance(final_labels[u], final_labels[v])
            edges_out.append((final_labels[u], final_labels[v], dist))

    # Avoid duplicates and output symmetric edges
    seen = set()
    unique_edges = []
    for u_label, v_label, dist in edges_out:
        edge_key = (u_label, v_label)
        if edge_key not in seen:
            seen.add(edge_key)
            unique_edges.append(f"{u_label}->{v_label}:{dist}")
            unique_edges.append(f"{v_label}->{u_label}:{dist}")

    unique_edges.sort()

    return total_score, unique_edges


# Example 1
# -----------
n = 4
edges = ["4->CAAATCCC", "4->ATTGCGAC", "5->CTGCGCTG", "5->ATGGACGA", "6->4", "6->5"]

score, out_edges = SmallParsimonyRooted(n, edges)
print(score)
for e in out_edges:
    print(e)
    
# Output: 16
    # ATAGACAA->ATAGACAC:1  ATAGACAA->ATGGACAA:1A  TAGACAC->ATAGACAA:1  ATAGACAC->ATTGCGAC:3
    # ATAGACAC->CAAATCCC:5  ATGGACAA->ATAGACAA:1  ATGGACAA->ATGGACGA:1  ATGGACAA->CTGCGCTG:5
    # ATGGACGA->ATGGACAA:1  ATTGCGAC->ATAGACAC:3  CAAATCCC->ATAGACAC:5  CTGCGCTG->ATGGACAA:5
        

# Example 2
# -----------
with open("dataset_30291_9.txt", "r") as file:
    lines = [line.strip() for line in file if line.strip()]
    n = int(lines[0])
    edges = lines[1:]

score, result_edges = SmallParsimonyRooted(n, edges)

with open("output.txt", "w") as file:
    file.write(str(score) + '\n')
    for edge in result_edges:
        file.write(edge + '\n')

# Output: 11831
    # AAAACCAACCAAGCAGTAAGATTCCACTCGCAACAATGCGGGGGTATCAGTCGGACTATGGGGAAAATAACAAACAT
    # ACAATACTCCCTGGGATTTATGAAAAAGACGGCCGTGAAAAAAACGATGAGATGAATCGCTCGGAAGCAGCTCTAAC
    # TGGATGCTACACCAA->AAGACCAGACCAGCAGTAAGATTCCACTCGCAACAGTGCAGAGCCATCAGTCGGACTATG
    # GGGATAATAAAATACACACAATACACCCTGGTGTTTATGAAAAAGACGGCCGAGAAAAAAAGGCTAAGATGAATCGC
    # TCGGAAGCAGCTCTAACTGGATGCTACAGCAA:21
    # AAAACCAACCAAGCAGTAAGATTCCACTCGCAACAATGCGGGGGTATCAGTCGGACTATGGGGAAAATAACAAACAT
    # ACAATACTCCCTGGGATTTATGAAAAAGACGGCCGTGAAAAAAACGATGAGATGAATCGCTCGGAAGCAGCTCTAAC
    # TGGATGCTACACCAA->CCCGGTCACAAAGCATGAACGTTCTACTCGCCAGGATGCGGGGGTTCACTTGTCAGTAAT
    # TGGCCCCGAACTAATTTCCAATAGTCTGCGGGATAGTTTAGGCTTGCGAGCCTGTAAACAACTTTGAGCCGCGTCAC
    # TCGCAAGCGGCTCCAAGCCGATGCTAGGCCGA:72 ...



# -----------------------------------------------
# Small Parsimony for Unrooted Trees
# -----------------------------------------------
import sys
import queue
import numpy as np

class SmallParsimony: # unrooted
    '''Find the most parsimonious labeling of the internal nodes in a rooted tree.
    
    Input: An integer n followed by an adjacency list for an unrooted binary tree with n leaves labeled by DNA strings.
    Output: The minimum parsimony score of this tree, followed by the adjacency list of the tree corresponding to labeling 
    internal nodes by DNA strings in order to minimize the parsimony score of the tree.'''

    def __init__(self):
        n, adj, nodes, lastEdge = self._input()
        s = self.runSmallParsimony(n, adj, nodes, lastEdge)
        self.printResults(s, adj, nodes)
        

    def _input(self):
        with open("dataset_30291_11.txt", "r") as file:
            data = [line.strip() for line in file if line.strip()]
        
        n = int(data[0])
        adj = [dict() for _ in range(n)]
        nodes = ['' for _ in range(n)]
        currNode = 0
        for d in data[1:]:
            d = d.split('->')
            try:
                p = int(d[0])
            except:
                p = currNode
                nodes[p] = d[0]
                currNode += 1
            try:
                c = int(d[1])
            except:
                continue
            if p > len(adj)-1 or c > len(adj)-1:
                adj.extend([dict() for _ in range(max([p,c])-len(adj)+1)])
            adj[p][c] = 0
            adj[c][p] = 0
        nodes.extend(['' for _ in range(len(adj)-n+1)])
        lastEdge = [int(i) for i in data[-1].split('->')]
        return n, adj, nodes, lastEdge
    
    def printResults(self, s, adj, nodes):
        print(s)
        for i, d in enumerate(adj):
            for j, w in d.items():
                print(nodes[i]+'->'+nodes[j]+':'+str(w))
    
    def charIndConversion(self):
        char2ind = {'A':0, 'C':1, 'G':2, 'T':3}
        ind2char = {0:'A', 1:'C', 2:'G', 3:'T'}
        return char2ind, ind2char
    
    def singleSmallParsimony(self, n, adjC, adjP, adj, nodes, char2ind, ind2char, charInd):
        s = [[np.inf]*4 for _ in range(len(adjC))]
        backtrack = [[(-1, -1) for _ in range(4)] for __ in range(len(adjC))]
        processed = [0 for _ in range(len(adjC))]
        ripe = set()
        for i in range(n):
            s[i][char2ind[nodes[i][charInd]]] = 0
            processed[i] = 1
            if len(adjP[i]) > 0:
                ripe.add(adjP[i][0])
        
        while len(ripe) > 0:
            v = ripe.pop()
            for k in range(4):
                l = [s[adjC[v][0]][i] + (0 if k == i else 1) for i in range(4)]
                r = [s[adjC[v][1]][i] + (0 if k == i else 1) for i in range(4)]
                largmin = np.argmin(l)
                rargmin = np.argmin(r)
                backtrack[v][k] = (largmin, rargmin)
                s[v][k] = l[largmin] + r[rargmin]
            processed[v] = 1
            if len(adjP[v]) > 0 and all([processed[u] for u in adjC[adjP[v][0]]]):
                ripe.add(adjP[v][0])
        
        ind = np.argmin(s[v])
        nodes[v] += ind2char[ind]
        smin = s[v][ind]

        q = queue.Queue()
        q.put((v, ind))
        while not q.empty():
            v, k = q.get()
            if len(adjC[v]) > 0:
                u, w = adjC[v]
                l, r = backtrack[v][k]
                
                if k != l:
                    adj[v][u] += 1
                    adj[u][v] += 1
                if k != r:
                    adj[v][w] += 1
                    adj[w][v] += 1
                if len(adjC[u]) > 0:
                    nodes[u] += ind2char[l]
                    nodes[w] += ind2char[r]
                    q.put((u, l))
                    q.put((w, r))        
        return smin
    
    def runSmallParsimony(self, n, adj, nodes, lastEdge):
        def dist(v, w):
            d = 0
            l = len(v)
            for i in range(l):
                if v[i] != w[i]:
                    d += 1
            return d

        char2ind, ind2char = self.charIndConversion()
        root = len(adj)
        del adj[lastEdge[0]][lastEdge[1]]
        del adj[lastEdge[1]][lastEdge[0]]
        adj.append(dict())
        adj[root][lastEdge[0]] = 0
        adj[lastEdge[0]][root] = 0
        adj[root][lastEdge[1]] = 0
        adj[lastEdge[1]][root] = 0
        adjC = [[] for _ in range(len(adj))]
        adjP = [[] for _ in range(len(adj))]
        for p in range(n, len(adj)):
            c = sorted(list(adj[p].keys()))
            adjC[p].append(c[0])
            adjC[p].append(c[1])
            adjP[c[0]].append(p)
            adjP[c[1]].append(p)
        s = 0
        for i in range(len(nodes[0])):
            s += self.singleSmallParsimony(n, adjC, adjP, adj, nodes, char2ind, ind2char, i)
        d = dist(nodes[lastEdge[0]], nodes[lastEdge[1]])
        del adj[root]
        del adj[lastEdge[0]][root]
        del adj[lastEdge[1]][root]
        adj[lastEdge[0]][lastEdge[1]] = d
        adj[lastEdge[1]][lastEdge[0]] = d
        return s



# Example 1
# -----------
data = ["4", "TCGGCCAA->4", "4->TCGGCCAA", "CCTGGCTG->4", "4->CCTGGCTG",
    "CACAGGAT->5", "5->CACAGGAT", "TGAGTACC->5", "5->TGAGTACC", "4->5", "5->4"]
SmallParsimony() 
# Output: 17
    # TCGGCCAA->CCAGGCAA:3 CCTGGCTG->CCAGGCAA:3 CACAGGAT->CAAGGAAA:4
    # TGAGTACC->CAAGGAAA:5 CCAGGCAA->TCGGCCAA:3 CCAGGCAA->CCTGGCTG:3
    # CCAGGCAA->CAAGGAAA:2 CAAGGAAA->CACAGGAT:4 CAAGGAAA->TGAGTACC:5
    # CAAGGAAA->CCAGGCAA:2


# Example 2
# -----------
with open("dataset_30291_11.txt", "r") as file:
    data = [line.strip() for line in file if line.strip()]
SmallParsimony() 
# Output: 515
    # ATCCCTCACATGCAGGCAAGCGACTGCACG->ATCCCCCCCATGCACACAACCCACTCCAAT:9
    # ACTCTCCCATGGTCCAGTTCTCAGGCTTAT->ATCCCCCCCATGCACACAACCCACTCCAAT:16
    # TTGTATTTGGTTCTTGATGGCCAGCGTACC->TTTCCATAGGTTCACTAAAGCCACTCCACA:15
    # TCTCCACAACTCGAGTGCTGCTTCTCCCCA->TTTCCATAGGTTCACTAAAGCCACTCCACA:13
    # CTAACGTCTGAATGACAAGTTCCTTCAATT->ATTCCCTCGGTTCGCTAAACACACTCAAAT:16



# -----------------------------------------------
# Find Two Nearest Neighbors of Tree
# -----------------------------------------------

def ParseInput(manual_lines=None):
    if manual_lines is not None:
        # First line: two nodes a and b separated by space
        first_line = manual_lines[0].strip()
        a, b = map(int, first_line.split())
        lines = manual_lines[1:]
    else:
        # Read from stdin
        a, b = map(int, input().split())
        import sys
        lines = [line.strip() for line in sys.stdin if line.strip()]
    
    adj = {}
    for line in lines:
        node_str, neighbor_str = line.split('->')
        node = int(node_str)
        neighbor = int(neighbor_str)
        if node not in adj:
            adj[node] = []
        adj[node].append(neighbor)
    return a, b, adj

def NearestNeighborInterchangesNeighbors(a, b, adj):
    a_neighbors = [x for x in adj[a] if x != b]
    b_neighbors = [x for x in adj[b] if x != a]

    if len(a_neighbors) != 2 or len(b_neighbors) != 2:
        raise ValueError("Invalid input: nodes a and b must have exactly 3 neighbors each.")

    a1, a2 = a_neighbors
    b1, b2 = b_neighbors

    def SwapSubtrees(adj, a, b, swap_a, swap_b):
        from copy import deepcopy
        new_adj = deepcopy(adj)

        new_adj[a].remove(swap_a)
        new_adj[swap_a].remove(a)

        new_adj[b].remove(swap_b)
        new_adj[swap_b].remove(b)

        new_adj[a].append(swap_b)
        new_adj[swap_b].append(a)

        new_adj[b].append(swap_a)
        new_adj[swap_a].append(b)

        return new_adj

    neighbor1 = SwapSubtrees(adj, a, b, a2, b1)
    neighbor2 = SwapSubtrees(adj, a, b, a2, b2)

    return neighbor1, neighbor2

def PrintAdj(adj):
    nodes = sorted(adj.keys())
    lines = []
    for node in nodes:
        for neighbor in sorted(adj[node]):
            lines.append(f"{node}->{neighbor}")
    print('\n'.join(lines))

def ParseInputFromFile(filename):
    with open(filename, "r") as file:
        lines = [line.strip() for line in file if line.strip()]
    a, b = map(int, lines[0].split())
    edges = lines[1:]
    adj = {}
    for line in edges:
        node_str, neighbor_str = line.split('->')
        node = int(node_str)
        neighbor = int(neighbor_str)
        if node not in adj:
            adj[node] = []
        adj[node].append(neighbor)
    return a, b, adj

def WriteOutputToFile(filename, neighbor1, neighbor2):
    with open(filename, "w") as file:
        def WriteAdj(adj):
            for node in sorted(adj.keys()):
                for neighbor in sorted(adj[node]):
                    file.write(f"{node}->{neighbor}\n")
        WriteAdj(neighbor1)
        file.write('\n')
        WriteAdj(neighbor2)


# Example 1
# -----------
manual_input = ["5 4", "0->4", "4->0", "1->4", "4->1", "2->5", 
                "5->2", "3->5", "5->3", "4->5", "5->4"]

a, b, adj = ParseInput(manual_input)
neighbor1, neighbor2 = NearestNeighborInterchangesNeighbors(a, b, adj)
PrintAdj(neighbor1)
print()
PrintAdj(neighbor2)

# Output: 0->5 1->4 2->5 3->4 4->1 4->3 4->5 5->0 5->2 5->4
        # 0->4 1->5 2->5 3->4 4->0 4->3 4->5 5->1 5->2 5->4



# Example 2
# -----------
a, b, adj = ParseInputFromFile("dataset_30292_6.txt")
neighbor1, neighbor2 = NearestNeighborInterchangesNeighbors(a, b, adj)
WriteOutputToFile("output.txt", neighbor1, neighbor2)

# Output: 0->32 1->32 2->33 3->33 4->34 5->34 6->35 7->35 8->36 9->36 10->37
    # 11->37 12->38 13->38 14->39 15->39 16->40 17->40 18->41 19->41 20->42 21->42
    # ... 61->60
    # 
    # 0->321->32 2->33 3->33 4->34 5->34 6->35 7->35 8->36 9->36 10->37 11->37
    # 12->38 13->38 14->39 15->39 16->40 17->40 18->41 19->4 20->42 ... 61->60


# -----------------------------------------------
# Nearest Neighbors Interchange for Large Parsimony
# -----------------------------------------------
import sys
import queue
import numpy as np
from copy import deepcopy

class LargeParsimony:
    '''Implement the nearest neighbor interchange heuristic for the Large Parsimony Problem.

    Input: An integer n, followed by an adjacency list for an unrooted binary tree whose n leaves are labeled by DNA strings and
    whose internal nodes are labeled by integers.
    Output: The parsimony score and unrooted labeled tree obtained after every step of the nearest neighbor interchange heuristic.
    Each step should be separated by a blank line.
    '''

    def __init__(self, input_data=None, filename=None):
        if input_data:
            n, adj, nodes, lastEdge = self.ParseInput(input_data)
        elif filename:
            n, adj, nodes, lastEdge = self.ReadFromFile(filename)
        else:
            raise ValueError("Provide either input_data (list) or filename (string)")
    
        trees = self.RunNearestNeighborInterchange(n, adj, nodes, lastEdge)
        if filename:
            self.SaveTrees(trees)  # save to file if input from file
        else:
            self.PrintTrees(trees)  # print if manual input

    def _Input(self):
        data = sys.stdin.read().strip().split('\n')
        n = int(data[0])
        adj = [dict() for _ in range(n)]
        nodes = ['' for _ in range(n)]
        currNode = 0
        for d in data[1:]:
            d = d.split('->')
            try:
                p = int(d[0])
            except:
                p = currNode
                nodes[p] = d[0]
                currNode += 1
            try:
                c = int(d[1])
            except:
                continue
            if p > len(adj)-1 or c > len(adj)-1:
                adj.extend([dict() for _ in range(max([p, c])-len(adj)+1)])
            adj[p][c] = 0
            adj[c][p] = 0
        nodes.extend(['' for _ in range(len(adj)-n+1)])
        lastEdge = [int(i) for i in data[-1].split('->')]
        return n, adj, nodes, lastEdge

    def ParseInput(self, data):
        data = [line.strip() for line in data if line.strip()]
        n = int(data[0])
        adj = [dict() for _ in range(n)]
        nodes = ['' for _ in range(n)]
        currNode = 0
        for d in data[1:]:
            d = d.split('->')
            try:
                p = int(d[0])
            except:
                p = currNode
                nodes[p] = d[0]
                currNode += 1
            try:
                c = int(d[1])
            except:
                continue
            if p > len(adj) - 1 or c > len(adj) - 1:
                adj.extend([dict() for _ in range(max([p, c]) - len(adj) + 1)])
            adj[p][c] = 0
            adj[c][p] = 0
        nodes.extend(['' for _ in range(len(adj) - n + 1)])
        lastEdge = [int(i) for i in data[-1].split('->')]
        return n, adj, nodes, lastEdge

    def ReadFromFile(self, filename):
        with open(filename, "r") as file:
            data = [line.strip() for line in file if line.strip()]
        return self.ParseInput(data)

    def PrintResults(self, s, adj, nodes):
        print(s)
        for i, d in enumerate(adj):
            for j, w in d.items():
                print(nodes[i] + '->' + nodes[j] + ':' + str(w))

    def PrintTrees(self, trees):
        for s, adj, nodes in trees:
            print(s)
            for i, d in enumerate(adj):
                for j, w in d.items():
                    print(nodes[i] + '->' + nodes[j] + ':' + str(w))
            print('')

    def SaveTrees(self, trees):
        with open("output.txt", "w") as file:
            for s, adj, nodes in trees:
                file.write(str(s) + '\n')
                for i, d in enumerate(adj):
                    for j, w in d.items():
                        file.write(nodes[i] + '->' + nodes[j] + ':' + str(w) + '\n')
                file.write('\n')

    def CharIndConversion(self):
        char2ind = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        ind2char = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
        return char2ind, ind2char

    def SingleSmallParsimony(self, n, adjC, adjP, adj, nodes, char2ind, ind2char, charInd):
        s = [[np.inf] * 4 for _ in range(len(adjC))]
        backtrack = [[(-1, -1) for _ in range(4)] for __ in range(len(adjC))]
        processed = [0 for _ in range(len(adjC))]
        ripe = set()
        for i in range(n):
            s[i][char2ind[nodes[i][charInd]]] = 0
            processed[i] = 1
            if len(adjP[i]) > 0:
                ripe.add(adjP[i][0])

        while len(ripe) > 0:
            v = ripe.pop()
            for k in range(4):
                l = [s[adjC[v][0]][i] + (0 if k == i else 1) for i in range(4)]
                r = [s[adjC[v][1]][i] + (0 if k == i else 1) for i in range(4)]
                largmin = np.argmin(l)
                rargmin = np.argmin(r)
                backtrack[v][k] = (largmin, rargmin)
                s[v][k] = l[largmin] + r[rargmin]
            processed[v] = 1
            if len(adjP[v]) > 0 and all([processed[u] for u in adjC[adjP[v][0]]]):
                ripe.add(adjP[v][0])

        ind = np.argmin(s[v])
        nodes[v] += ind2char[ind]
        smin = s[v][ind]

        q = queue.Queue()
        q.put((v, ind))
        while not q.empty():
            v, k = q.get()
            if len(adjC[v]) > 0:
                u, w = adjC[v]
                l, r = backtrack[v][k]

                if k != l:
                    adj[v][u] += 1
                    adj[u][v] += 1
                if k != r:
                    adj[v][w] += 1
                    adj[w][v] += 1

                if len(adjC[u]) > 0:
                    nodes[u] += ind2char[l]
                    q.put((u, l))
                if len(adjC[w]) > 0:
                    nodes[w] += ind2char[r]
                    q.put((w, r))

        return smin

    def RunSmallParsimony(self, n, adj, nodes, lastEdge):
        def Dist(v, w):
            return sum(1 for i in range(len(v)) if v[i] != w[i])

        char2ind, ind2char = self.CharIndConversion()
        root = len(adj)
        del adj[lastEdge[0]][lastEdge[1]]
        del adj[lastEdge[1]][lastEdge[0]]
        adj.append(dict())
        adj[root][lastEdge[0]] = 0
        adj[lastEdge[0]][root] = 0
        adj[root][lastEdge[1]] = 0
        adj[lastEdge[1]][root] = 0
        adjC = [[] for _ in range(len(adj))]
        adjP = [[] for _ in range(len(adj))]
        q = queue.Queue()
        q.put(root)
        visited = [False for _ in range(len(adj))]
        visited[root] = True
        while not q.empty():
            curr = q.get()
            for v in adj[curr].keys():
                if not visited[v]:
                    adjP[v].append(curr)
                    visited[v] = True
                    q.put(v)
        for u, d in enumerate(adjP):
            for v in d:
                adjC[v].append(u)
        s = 0
        for i in range(len(nodes[0])):
            s += self.SingleSmallParsimony(n, adjC, adjP, adj, nodes, char2ind, ind2char, i)
        d = Dist(nodes[lastEdge[0]], nodes[lastEdge[1]])
        del adj[root]
        del adj[lastEdge[0]][root]
        del adj[lastEdge[1]][root]
        adj[lastEdge[0]][lastEdge[1]] = d
        adj[lastEdge[1]][lastEdge[0]] = d
        return s, adj, nodes

    def FindNearestNeighbors(self, edge, adj):
        adj1 = deepcopy(adj)
        adj2 = deepcopy(adj)

        del adj1[edge[0]][edge[1]]
        del adj1[edge[1]][edge[0]]
        e0 = list(adj1[edge[0]].keys())
        e1 = list(adj1[edge[1]].keys())
        adj1[edge[0]][e1[0]] = 0
        adj1[edge[1]][e0[0]] = 0
        adj1[e1[0]][edge[0]] = 0
        adj1[e0[0]][edge[1]] = 0
        del adj1[e1[0]][edge[1]]
        del adj1[e0[0]][edge[0]]
        del adj1[edge[0]][e0[0]]
        del adj1[edge[1]][e1[0]]
        adj1[edge[0]][edge[1]] = 0
        adj1[edge[1]][edge[0]] = 0

        adj2[edge[0]][e1[1]] = 0
        adj2[edge[1]][e0[0]] = 0
        adj2[e1[1]][edge[0]] = 0
        adj2[e0[0]][edge[1]] = 0
        del adj2[e1[1]][edge[1]]
        del adj2[e0[0]][edge[0]]
        del adj2[edge[0]][e0[0]]
        del adj2[edge[1]][e1[1]]
        return adj1, adj2

    def RunNearestNeighborInterchange(self, n, adj, nodes, lastEdge):
        trees = []
        score = np.inf
        newScore, newAdj, newNodes = self.RunSmallParsimony(n, adj, deepcopy(nodes), lastEdge)
        while newScore < score:
            score = newScore
            adj = newAdj
            visited = set()
            for v in range(n, len(adj)):
                for u in adj[v].keys():
                    if u >= n and not (v, u) in visited:
                        adj1, adj2 = self.FindNearestNeighbors([v, u], adj)
                        for i, a in enumerate(adj1):
                            adj1[i] = dict.fromkeys(a, 0)
                        for i, a in enumerate(adj2):
                            adj2[i] = dict.fromkeys(a, 0)
                        neighborScore, neighborAdj, neighborNodes = self.RunSmallParsimony(n, adj1, deepcopy(nodes), [v, u])
                        if neighborScore < newScore:
                            newScore = neighborScore
                            newAdj = neighborAdj
                            newNodes = neighborNodes
                        neighborScore, neighborAdj, neighborNodes = self.RunSmallParsimony(n, adj2, deepcopy(nodes), [v, u])
                        if neighborScore < newScore:
                            newScore = neighborScore
                            newAdj = neighborAdj
                            newNodes = neighborNodes
                        visited.add((v, u))
                        visited.add((u, v))
            if newScore < score:
                trees.append((newScore, newAdj, newNodes))
        return trees


# Example 1
# ----------
input_data = ["5", "GCAGGGTA->5", "TTTACGCG->5", "CGACCTGA->6", "GATTCCAC->6", "5->TTTACGCG",
    "5->GCAGGGTA", "5->7", "TCCGTAGT->7", "7->5", "7->6", "7->TCCGTAGT",
    "6->GATTCCAC", "6->CGACCTGA", "6->7"]
LargeParsimony(input_data)

# Output: 22
    # GCAGGGTA->GCAGCGGA:2 TTTACGCG->TCAGCGGA:5 CGACCTGA->GAACCCGA:3 GATTCCAC->GAACCCGA:4
    # TCCGTAGT->TCAGCGGA:4 TCAGCGGA->TTTACGCG:5 TCAGCGGA->TCCGTAGT:4 TCAGCGGA->GCAGCGGA:1
    # GAACCCGA->CGACCTGA:3 GAACCCGA->GATTCCAC:4 GAACCCGA->GCAGCGGA:3 GCAGCGGA->GAACCCGA:3
    # GCAGCGGA->GCAGGGTA:2 GCAGCGGA->TCAGCGGA:1

    # 21
    # GCAGGGTA->GCAGCGGA:2 TTTACGCG->TCTGCGGA:4 CGACCTGA->GCAGCGGA:4 GATTCCAC->GCTGCGGA:5
    # TCCGTAGT->TCTGCGGA:4 TCTGCGGA->TTTACGCG:4 TCTGCGGA->TCCGTAGT:4 TCTGCGGA->GCTGCGGA:1
    # GCTGCGGA->GATTCCAC:5 GCTGCGGA->TCTGCGGA:1 GCTGCGGA->GCAGCGGA:1 GCAGCGGA->GCAGGGTA:2
    # GCAGCGGA->CGACCTGA:4 GCAGCGGA->GCTGCGGA:1


# Example 2
# ----------
lp = LargeParsimony(filename = "dataset_30292_8.txt")

# Output: 157
    # GCTTAGCTACTCACAAAGATAGGCTTAGAGAAGCTCACCG->ACACCCAGGCTCTCAAAGGTAGGCTTGGACATGACGTCCG:17
    # ATCCCTGGGGGCTAGACTGCCCTCTAGTCCCTAACGTCCG->ACACCTAGGAGCTAAACGGCAGTCTAGTCCATAACGTCCG:9
    # AAATTCTTGAACTTACTTCTGATTTCTCTCAGCATGGATG->ACACGTATCAACTAACCACGAGTTTATTCCATCAAGGCAG:19 ...

    # 156
    # GCTTAGCTACTCACAAAGATAGGCTTAGAGAAGCTCACCG->ACACCCAAACTAACAAAGGTAGGCTTAGACAAGACGTCCG:14
    # ATCCCTGGGGGCTAGACTGCCCTCTAGTCCCTAACGTCCG->ACACCTAGGAGATAAACCGCAGACTAGTCCATCACGTCCG:12 ...


    # 155
    # GCTTAGCTACTCACAAAGATAGGCTTAGAGAAGCTCACCG->ACACCCATACTAACAAAGATAGGCGTAGACAAGACGACCG:12
    # ATCCCTGGGGGCTAGACTGCCCTCTAGTCCCTAACGTCCG->ACACCTATGAAATAAACTCGCGTCTAGTCCATCACGTCCG:13

    # 153
    # GCTTAGCTACTCACAAAGATAGGCTTAGAGAAGCTCACCG->ACACCCATTCTAACAAAGATAGGCGTAGAGACGACGACCG:13
    # ATCCCTGGGGGCTAGACTGCCCTCTAGTCCCTAACGTCCG->ACACCTATCAGATAGACTCGCGTCTAGTCCATCACGTCCG:12

    # 150
    # GCTTAGCTACTCACAAAGATAGGCTTAGAGAAGCTCACCG->ACACCCCTTCTAAAAAAGATAGGCGTAGAGACGACGACCG:13
    # ATCCCTGGGGGCTAGACTGCCCTCTAGTCCCTAACGTCCG->ACACCTAGGAGATAGACTCTCGTCTAGTCCATCACGTCCG:10

    # 149
    # GCTTAGCTACTCACAAAGATAGGCTTAGAGAAGCTCACCG->ACTAAACTTCTATAGAAGTTAGACGTGGAGACGATGACCG:15
    # ATCCCTGGGGGCTAGACTGCCCTCTAGTCCCTAACGTCCG->ACACCTAAGAGATAGACTCTCGTCTAGTCCATCACGTCCG:11


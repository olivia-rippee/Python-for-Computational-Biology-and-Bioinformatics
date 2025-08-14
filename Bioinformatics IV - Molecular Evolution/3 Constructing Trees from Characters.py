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


# Examples
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

from collections import defaultdict

def SmallParsimonyUnrooted(n, edges):
    adj = defaultdict(list)
    labels = {}
    node_id = 0
    node_map = {}

    # Step 1: Convert all nodes and labels to integer nodes
    for line in edges:
        u, v = line.split("->")

        # Map u
        try:
            u_int = int(u)
        except ValueError:
            if u not in node_map:
                node_map[u] = node_id
                labels[node_id] = u
                node_id += 1
            u_int = node_map[u]

        # Map v
        try:
            v_int = int(v)
        except ValueError:
            if v not in node_map:
                node_map[v] = node_id
                labels[node_id] = v
                node_id += 1
            v_int = node_map[v]

        adj[u_int].append(v_int)
        adj[v_int].append(u_int)

    # Step 2: Pick arbitrary edge to root the tree
    found = False
    for node in adj:
        for neighbor in adj[node]:
            u, v = node, neighbor
            found = True
            break
        if found:
            break

    # Step 3: Insert a new root node
    new_root = max(adj.keys()) + 1
    adj[new_root] = [u, v]
    adj[u].remove(v)
    adj[v].remove(u)
    adj[u].append(new_root)
    adj[v].append(new_root)

    # Step 4: Build edge list for rooted tree
    def BuildEdges(start_node):
        result = []
        visited = set()
        stack = [(start_node, None)]

        while stack:
            node, parent = stack.pop()
            visited.add(node)

            for neighbor in adj[node]:
                if neighbor == parent:
                    continue
                if neighbor in labels:
                    result.append(f"{node}->{labels[neighbor]}")
                else:
                    result.append(f"{node}->{neighbor}")
                stack.append((neighbor, node))
        return result

    rooted_edges = BuildEdges(new_root)

    # Step 5: Call rooted version
    score, full_edges = SmallParsimonyRooted(new_root + 1, rooted_edges)

    # Step 6: Remove artificial root from output
    filtered_edges = []
    for line in full_edges:
        u, rest = line.split("->")
        v, dist = rest.split(":")
        if u == str(new_root) or v == str(new_root):
            continue
        filtered_edges.append(line)

    return score, filtered_edges


def SmallParsimonyRooted(n, edges):
    '''Find the most parsimonious labeling of the internal nodes in a rooted tree.
    
    Input: An integer n followed by an adjacency list for an unrooted binary tree with n leaves labeled by DNA strings.
    Output: The minimum parsimony score of this tree, followed by the adjacency list of the tree corresponding to labeling 
    internal nodes by DNA strings in order to minimize the parsimony score of the tree.'''

    children = defaultdict(list)
    parent = {}
    label = {}

    next_leaf_node = 1000000

    for line in edges:
        u, v = line.split("->")
        u = int(u)
        if all(ch in "ACGT" for ch in v):
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

        for leaf in leaves:
            Tag[leaf] = 1
            for i, b in enumerate(Bases):
                S[leaf][i] = 0 if label[leaf][pos] == b else float('inf')

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

    for node in final_labels:
        final_labels[node] = ''.join(final_labels[node])

    def HammingDistance(s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    edges_out = []
    for u in children:
        for v in children[u]:
            dist = HammingDistance(final_labels[u], final_labels[v])
            edges_out.append((final_labels[u], final_labels[v], dist))

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



# Examples
# ---------
n = 4
edges = ["TCGGCCAA->4", "4->TCGGCCAA", "CCTGGCTG->4", "4->CCTGGCTG",
        "CACAGGAT->5", "5->CACAGGAT", "TGAGTACC->5", "5->TGAGTACC", "4->5", "5->4"]

score, edges = SmallParsimonyUnrooted(n, edges)
print(score)
for line in edges:
    print(line)







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


# Examples
# -----------
manual_input = ["5 4", "0->4", "4->0", "1->4", "4->1", "2->5", 
                "5->2", "3->5", "5->3", "4->5", "5->4"]

a, b, adj = ParseInput(manual_input)
neighbor1, neighbor2 = NearestNeighborInterchangesNeighbors(a, b, adj)
PrintAdj(neighbor1)
print()
PrintAdj(neighbor2)

# Output: 0->5 1->4 2->5 3->4 4->1 4->3 4->5 5->0 5->2 5->4
#         0->4 1->5 2->5 3->4 4->0 4->3 4->5 5->1 5->2 5->4




input_filename = "dataset_30292_6.txt"
output_filename = "output.txt"
a, b, adj = ParseInputFromFile(input_filename)
neighbor1, neighbor2 = NearestNeighborInterchangesNeighbors(a, b, adj)
WriteOutputToFile(output_filename, neighbor1, neighbor2)

# Output: 0->32 1->32 2->33 3->33 4->34 5->34 6->35 7->35 8->36 9->36 10->37
# 11->37 12->38 13->38 14->39 15->39 16->40 17->40 18->41 19->41 20->42 21->42
# ... 61->60
# 
# 0->321->32 2->33 3->33 4->34 5->34 6->35 7->35 8->36 9->36 10->37 11->37
# 12->38 13->38 14->39 15->39 16->40 17->40 18->41 19->4 20->42 ... 61->60




# -----------------------------------------------
# Nearest Neighbors Interchange for Large Parsimony
# -----------------------------------------------
from collections import defaultdict, deque
import copy

def HammingDistance(s1, s2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def ReadInputFromFile(filename):
    with open(filename, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    n = int(lines[0])
    adj_lines = lines[1:]
    return n, adj_lines

def WriteOutputToFile(filename, output_lines):
    with open(filename, "w") as f:
        f.write('\n'.join(output_lines))

def ParseInput(n, lines):
    adj = defaultdict(list)
    labels = {}
    for line in lines:
        left, right = line.split('->')
        if left.isdigit():
            u = int(left)
        else:
            u = left
        if right.isdigit():
            v = int(right)
        else:
            v = right
        adj[u].append(v)
        adj[v].append(u)
        if not str(u).isdigit():
            labels[u] = u
        if not str(v).isdigit():
            labels[v] = v
    return adj, labels

def SmallParsimony(adj, labels):
    nodes = list(adj.keys())
    n = len(next(iter(labels.values())))
    S = {v: [set() for _ in range(n)] for v in nodes}

    # Find a root: pick an internal node (integer) if any, else any node
    root = None
    for v in nodes:
        if isinstance(v, int):
            root = v
            break
    if root is None:
        root = nodes[0]

    parent = {root: None}
    order = []
    queue = deque([root])
    visited = set([root])  # Track visited nodes to prevent cycles

    while queue:
        u = queue.popleft()
        order.append(u)
        for w in adj[u]:
            if w not in visited:
                parent[w] = u
                visited.add(w)
                queue.append(w)

    # Initialize S for leaves and internal nodes
    for v in nodes:
        if v in labels:  # Leaf labeled by DNA string
            for i, ch in enumerate(labels[v]):
                S[v][i] = {ch}
        else:  # Internal node: initially all possible nucleotides at each position
            for i in range(n):
                S[v][i] = {'A', 'C', 'G', 'T'}

    # Bottom-up phase: postorder traversal
    for v in reversed(order):
        if v in labels:
            continue  # Leaf already assigned
        children = [w for w in adj[v] if w != parent[v]]
        for i in range(n):
            intersect = S[children[0]][i].intersection(S[children[1]][i])
            if intersect:
                S[v][i] = intersect
            else:
                S[v][i] = S[children[0]][i].union(S[children[1]][i])

    # Top-down phase: assign sequences starting from root
    seq = {}
    seq[root] = ''.join(sorted(S[root][i])[0] for i in range(n))
    queue = deque([root])
    while queue:
        v = queue.popleft()
        for w in adj[v]:
            if w == parent.get(v):
                continue
            chars = []
            for i in range(n):
                if seq[v][i] in S[w][i]:
                    chars.append(seq[v][i])
                else:
                    chars.append(sorted(S[w][i])[0])
            seq[w] = ''.join(chars)
            queue.append(w)

    # Compute total parsimony score
    total_score = 0
    counted = set()
    for v in adj:
        for w in adj[v]:
            if (v, w) not in counted and (w, v) not in counted:
                total_score += sum(ch1 != ch2 for ch1, ch2 in zip(seq[v], seq[w]))
                counted.add((v, w))

    # Create adjacency with sequences as labels for output
    new_adj = {}
    for v in adj:
        v_label = seq[v] if isinstance(v, int) else v
        new_adj[v_label] = []
        for w in adj[v]:
            w_label = seq[w] if isinstance(w, int) else w
            new_adj[v_label].append(w_label)

    return total_score, new_adj, seq

def FindInternalEdges(adj):
    internal_edges = []
    for v in adj:
        for w in adj[v]:
            # Check both are integers before comparing
            if isinstance(v, int) and isinstance(w, int):
                if v < w:
                    internal_edges.append((v, w))
    return internal_edges

def SwapSubtrees(adj, a, b, swap_a, swap_b):
    new_adj = copy.deepcopy(adj)
    new_adj[a].remove(swap_a)
    new_adj[swap_a].remove(a)
    new_adj[b].remove(swap_b)
    new_adj[swap_b].remove(b)
    new_adj[a].append(swap_b)
    new_adj[swap_b].append(a)
    new_adj[b].append(swap_a)
    new_adj[swap_a].append(b)
    return new_adj

def GenerateNNINeighbors(adj, edge):
    a, b = edge
    a_neighbors = [x for x in adj[a] if x != b]
    b_neighbors = [x for x in adj[b] if x != a]
    neighbors = []
    neighbors.append(SwapSubtrees(adj, a, b, a_neighbors[0], b_neighbors[0]))
    neighbors.append(SwapSubtrees(adj, a, b, a_neighbors[0], b_neighbors[1]))
    neighbors.append(SwapSubtrees(adj, a, b, a_neighbors[1], b_neighbors[0]))
    neighbors.append(SwapSubtrees(adj, a, b, a_neighbors[1], b_neighbors[1]))
    return neighbors

def FormatSolution(score, adj):
    lines = [str(score)]
    edges_printed = set()
    for v in sorted(adj.keys()):
        for w in sorted(adj[v]):
            if (v,w) not in edges_printed and (w,v) not in edges_printed:
                dist = HammingDistance(v, w)
                lines.append(f"{v}->{w}:{dist}")
                edges_printed.add((v,w))
    lines.append("")  # blank line to separate steps
    return lines

def RunNNI(adj, labels):
    output_lines = []
    improved = True
    current_adj = adj
    current_labels = labels
    score, labeled_adj, seq = SmallParsimony(current_adj, current_labels)
    output_lines.extend(FormatSolution(score, labeled_adj))
    while improved:
        improved = False
        internal_edges = FindInternalEdges(current_adj)
        best_score = score
        best_adj = current_adj
        for edge in internal_edges:
            neighbors = GenerateNNINeighbors(current_adj, edge)
            for neighbor_adj in neighbors:
                neighbor_score, labeled_neighbor_adj, _ = SmallParsimony(neighbor_adj, current_labels)
                if neighbor_score < best_score:
                    best_score = neighbor_score
                    best_adj = neighbor_adj
                    improved = True
        if improved:
            current_adj = best_adj
            score, labeled_adj, seq = SmallParsimony(current_adj, current_labels)
            output_lines.extend(FormatSolution(score, labeled_adj))
    return output_lines


# Examples
# ----------
n = 5
adj = ["GCAGGGTA->5", "TTTACGCG->5", "CGACCTGA->6", "GATTCCAC->6", "5->TTTACGCG",
    "5->GCAGGGTA", "5->7", "TCCGTAGT->7", "7->5", "7->6", "7->TCCGTAGT",
    "6->GATTCCAC", "6->CGACCTGA", "6->7"]
RunNNI(n, adj)




n, adj_lines = ReadInputFromFile("dataset_30292_8.txt")
adj, labels = ParseInput(n, adj_lines)
output_lines = RunNNI(adj, labels)
WriteOutputToFile("output.txt", output_lines)



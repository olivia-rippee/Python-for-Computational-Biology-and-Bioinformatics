import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics VI - Finding Mutations in DNA and Proteins/Data")

# -----------------------------------------------
# Trie Construction
# -----------------------------------------------

def TrieConstruction(Patterns):
    '''Construct a trie (directed tree of patterns).

    Input: A space-separated collection of strings Patterns.
    Output: The adjacency list corresponding to Trie(Patterns), in the following format. 
    If Trie(Patterns) has n nodes, first label the root with 0 and then label the remaining 
    nodes with the integers 1 through n - 1 in any order you like. Each edge of the adjacency 
    list of Trie(Patterns) will be encoded by a triple: the first two members of the triple 
    must be the integers labeling the initial and terminal nodes of the edge, respectively; 
    the third member of the triple must be the symbol labeling the edge.'''

    trie = {}
    nodeId = 0

    for pattern in Patterns:
        currentNode = 0
        for symbol in pattern:
            found = False
            if currentNode in trie:
                for childNode, edgeSymbol in trie[currentNode]:
                    if edgeSymbol == symbol:
                        currentNode = childNode
                        found = True
                        break
            else:
                trie[currentNode] = []

            if not found:
                nodeId += 1
                if currentNode not in trie:
                    trie[currentNode] = []
                trie[currentNode].append((nodeId, symbol))
                currentNode = nodeId

    return trie

def PrintTrie(Trie, OutputFile=None):
    if OutputFile:
        with open(OutputFile, "w") as file:
            for node in Trie:
                for child, symbol in Trie[node]:
                    file.write(f"{node} {child} {symbol}\n")
    else:
        for node in Trie:
            for child, symbol in Trie[node]:
                print(f"{node} {child} {symbol}")



# Example 1
# -----------
input_patterns = "ATAGA ATC GAT".split()
trie = TrieConstruction(input_patterns)
PrintTrie(trie)
# Output:
    # 0 1 A
    # 0 7 G
    # 1 2 T
    # 2 3 A
    # 2 6 C
    # 3 4 G
    # 4 5 A
    # 7 8 A
    # 8 9 T


# Example 2
# -----------
with open("dataset_30220_4.txt", "r") as file:
    patterns = file.read().strip().split()
trie = TrieConstruction(patterns)
PrintTrie(trie, "output.txt")
# Output: 
    # 0 1 C
    # 0 93 G
    # 0 552 T
    # 0 635 A
    # 1 2 C
    # 1 273 G
    # 1 721 A
    # 1 1673 T
    
    
# -----------------------------------------------
# Trie Matching (Multiple Pattern Matching)
# -----------------------------------------------

def TrieMatching(Text, Patterns):
    '''Find whether any strings in Patterns match a substring of Text. 
    
    Input: A string Text and a space-separated collection of strings Patterns.
    Output: All starting positions in Text where a string from Patterns appears as a substring.'''
    
    trie = TrieConstruction(Patterns)
    result = {pattern: [] for pattern in Patterns}
    for i in range(len(Text)):
        matchIndex = PrefixTrieMatching(Text[i:], trie, i)
        if matchIndex is not None:
            # Figure out which pattern it was by matching it against patterns
            for pattern in Patterns:
                if Text[i:i+len(pattern)] == pattern:
                    result[pattern].append(i)
    return result

def PrefixTrieMatching(Text, Trie, StartIndex):
    '''Check whether any string from Patterns matches a prefix of Text.'''
    
    index = 0
    symbol = Text[index] if index < len(Text) else ''
    v = 0
    while True:
        if IsLeaf(v, Trie):
            return StartIndex
        elif v in Trie:
            matched = False
            for child, edgeSymbol in Trie[v]:
                if symbol == edgeSymbol:
                    index += 1
                    symbol = Text[index] if index < len(Text) else ''
                    v = child
                    matched = True
                    break
            if not matched:
                return None
        else:
            return None

def TrieConstruction(Patterns):
    trie = {}
    nodeId = 0

    for pattern in Patterns:
        currentNode = 0
        for symbol in pattern:
            found = False
            if currentNode in trie:
                for childNode, edgeSymbol in trie[currentNode]:
                    if edgeSymbol == symbol:
                        currentNode = childNode
                        found = True
                        break
            else:
                trie[currentNode] = []

            if not found:
                nodeId += 1
                if currentNode not in trie:
                    trie[currentNode] = []
                trie[currentNode].append((nodeId, symbol))
                currentNode = nodeId

    return trie

def IsLeaf(Node, Trie):
    return Node not in Trie or len(Trie[Node]) == 0

def PrintMatches(Matches, OutputFile=None):
    if OutputFile:
        with open(OutputFile, 'w') as f:
            for pattern in Matches:
                positions = ' '.join(map(str, Matches[pattern]))
                f.write(f"{pattern}: {positions}\n")
    else:
        for pattern in Matches:
            positions = ' '.join(map(str, Matches[pattern]))
            print(f"{pattern}: {positions}")


# Example 1
# ------------
Text = "AATCGGGTTCAATCGGGGT"
Patterns = ["ATCG", "GGGT"]

matches = TrieMatching(Text, Patterns)
PrintMatches(matches)
# Output:
    # ATCG: 1 11
    # GGGT: 4 15


# Example 2
# ------------
with open("dataset_30220_8.txt", "r") as f:
    lines = f.read().strip().split('\n')
    text = lines[0]
    patterns = lines[1].split()

matches = TrieMatching(text, patterns)
PrintMatches(matches, "output.txt")
# Output:
    # CGCGAGGCG: 2472 3715 4806 4813 6859 6866 9196
    # ACTGCATAC: 189 5048 7049 7056 7773
    # GGCTCACGG: 173 180 1510 1517 2255 5114 5121
    # GTGACGCGT: 2075 2082 2842 4468 9107
    # GCAAAATGC: 1094 1396 1403 2965 2972 3886 3893 5912 5919
    # TCGCCTCTC: 2205 2751 2758 4777 4784 7606 8905 8912
    # TATTATCTA: 837 1889 1896 3659 3666 6474 6481



# -----------------------------------------------
# Maximal Non-Branching Paths
# -----------------------------------------------

from collections import defaultdict

def MaximalNonBranchingPaths(adj_list):
    '''Generate all non-branching paths in a graph. 
    Iterates through all nodes of the graph that are not 1-in-1-out nodes,  
    generates all non-branching paths starting at each such node, and finds 
    all isolated cycles in the graph.

    Input: The adjacency list of a graph whose nodes are integers.
    Output: The collection of all maximal nonbranching paths in this graph.
    '''
    
    # Calculate indegree and outdegree for each node
    indegree = defaultdict(int)
    outdegree = defaultdict(int)
    
    for v in adj_list:
        outdegree[v] = len(adj_list[v])
        for w in adj_list[v]:
            indegree[w] += 1
     
    def Is1In1Out(node):
        '''Check if a node is 1-in-1-out'''
        return indegree[node] == 1 and outdegree[node] == 1
    
    paths = []
    visited_edges = set()  # Keep track of edges visited
    
    # 1. Find non-branching paths starting at nodes that are NOT 1-in-1-out
    for v in adj_list:
        if not Is1In1Out(v):
            if outdegree[v] > 0:
                for w in adj_list[v]:
                    edge = (v, w)
                    if edge not in visited_edges:
                        path = [v, w]
                        visited_edges.add(edge)
                        current = w
                        while Is1In1Out(current):
                            next_node = adj_list[current][0]  # only one outgoing edge
                            edge = (current, next_node)
                            if edge in visited_edges:
                                break  # Avoid infinite loops
                            path.append(next_node)
                            visited_edges.add(edge)
                            current = next_node
                        paths.append(path)
    
    # 2. Find isolated cycles composed entirely of 1-in-1-out nodes
    for v in adj_list:
        if Is1In1Out(v):
            for w in adj_list[v]:
                if (v, w) not in visited_edges:
                    cycle = [v, w]
                    visited_edges.add((v, w))
                    current = w
                    while current != v:
                        next_node = adj_list[current][0]
                        if (current, next_node) in visited_edges:
                            break
                        cycle.append(next_node)
                        visited_edges.add((current, next_node))
                        current = next_node
                    if current == v:
                        paths.append(cycle)
    
    return paths

def ReadGraphFromFile(filename):
    graph = defaultdict(list)
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            node_part, neighbors_part = line.split(':')
            node = int(node_part.strip())
            neighbors = list(map(int, neighbors_part.strip().split()))
            graph[node].extend(neighbors)
    return graph


# Example 1
# ----------
adj_list = {1: [2], 2: [3], 3: [4, 5], 6: [7], 7: [6]}

nonbranching_paths = MaximalNonBranchingPaths(adj_list)
for path in nonbranching_paths:
    print(' '.join(map(str, path)))
    
# Output:
    # 1 2 3
    # 3 4
    # 3 5
    # 6 7 6


# -----------------------------------------------
# Suffix Tree Construction
# -----------------------------------------------

class Node:
    def __init__(self):
        self.children = {}  # symbol -> (Node, position, length)
        self.leaf_label = None
        
def ModifiedSuffixTrieConstruction(Text):
    '''Construct the suffix trie of a string.
    Input: A string Text.
    Output: The root of the trie.'''
    
    root = Node()
    n = len(Text)
    for i in range(n):
        currentNode = root
        for j in range(i, n):
            currentSymbol = Text[j]
            if currentSymbol in currentNode.children:
                child, pos, length = currentNode.children[currentSymbol]
                currentNode = child
            else:
                newNode = Node()
                currentNode.children[currentSymbol] = (newNode, j, 1)
                currentNode = newNode
        # label the leaf with starting position of suffix
        currentNode.leaf_label = i
    return root

def ModifiedSuffixTreeConstruction(root):
    '''Construct the suffix tree of a string.
    Input: The root of the trie.
    Output: A space-separated list of the edge labels of SuffixTree(Text).'''
    
    # Compress non-branching paths
    def Compress(node):
        for symbol, (child, pos, length) in list(node.children.items()):
            # Compress child if it has only one child and is not a leaf
            while len(child.children) == 1 and child.leaf_label is None:
                next_symbol, (grandchild, gpos, glen) = next(iter(child.children.items()))
                # Combine edge labels by extending length
                length += glen
                child = grandchild
            node.children[symbol] = (child, pos, length)
            Compress(child)
    Compress(root)
    return root

def CollectEdges(node, Text, edges):
    for symbol, (child, pos, length) in node.children.items():
        edges.append(Text[pos:pos+length])
        CollectEdges(child, Text, edges)


# Example 1
# -----------
Text = "ATAAATG$"
root = ModifiedSuffixTrieConstruction(Text)
root = ModifiedSuffixTreeConstruction(root)
edges = []

CollectEdges(root, Text, edges)
print(' '.join(edges))
# Output: A T AAATG$ G$ A ATG$ TG$ T AAATG$ G$ G$ $


# Example 2
# -----------
with open("dataset_30222_4.txt", "r") as file:
    Text = file.read().strip()

root = ModifiedSuffixTrieConstruction(Text)
root = ModifiedSuffixTreeConstruction(root)
edges = []

CollectEdges(root, Text, edges)
with open("output.txt", "w") as file:
    file.write(' '.join(edges) + '\n')
# Output: 
    # C G T C A CTCTAGAATCCATGTGAC.... ... $


# How many leaves will SuffixTree("TCTGAGCCCTACTGTCGAGAAATATGTATCTCGCCCCCGCAGCTT$") have?
# ----------------------------------------------------------------------------------------
Text = "TCTGAGCCCTACTGTCGAGAAATATGTATCTCGCCCCCGCAGCTT$"
len(Text) # Output: 46
# number of leaves = length of genome



# -----------------------------------------------
# Longest Repeat
# -----------------------------------------------

class Node:
    total = 0
    def __init__(self, startIdx=0, depth=0):
        Node.total += 1
        self.id = self.total
        self.startIdx = startIdx
        self.depth = depth

class Edge:
    def __init__(self, startIdx, endIdx, text, startNode, endNode=None, leafLabel=None):
        self.startIdx = startIdx
        self.endIdx = endIdx
        self.text = text
        self.startNode = startNode
        self.endNode = endNode
        self.leafLabel = leafLabel

    def Length(self):
        return self.endIdx - self.startIdx + 1
    
    def Str(self):
        return self.text[self.startIdx:self.endIdx+1]
    
    def StartChar(self):
        return self.text[self.startIdx]    
    
    def __str__(self):
        return self.StartChar()

class SuffixTree:  # Naive algorithm
    def __init__(self, text):
        self.root = Node()
        self.text = text
        self.tree = dict()
        self.Build(self.root, self.text)
        lr = self.LongestRepeat(self.tree, self.text)
        print(lr)
        with open("output.txt", "w") as file:
            file.write(lr)

    def Match(self, i, root, text):
        l = len(text)
        currNode = root
        atNode = True
        for j in range(i, l):
            if atNode:
                currPos = 0
                if not text[j] in self.tree[currNode]:
                    return (currNode, None, j, -1)
                else:
                    currEdge = self.tree[currNode][text[j]]
                    currString = currEdge.Str()
                    lrString = len(currString) - 1
                    if lrString == 0:
                        currNode = currEdge.endNode
                        continue
                    else:
                        atNode = False                    
            else:
                currPos += 1
                if text[j] != currString[currPos]:
                    return (currNode, currEdge, j, currEdge.startIdx + currPos)
                else:
                    lrString -= 1
                    if lrString == 0:
                        currNode = currEdge.endNode
                        atNode = True

    def AddEdge(self, node, startIdx, endIdx, leafLabel):
        newEdge = Edge(startIdx, endIdx, self.text, node, None, leafLabel)
        self.tree[node][newEdge.StartChar()] = newEdge

    def SplitEdge(self, edge, startIdx, endIdx, cutIdx, leafLabel):
        newNode = Node(leafLabel, startIdx - leafLabel)
        newEdge = Edge(startIdx, endIdx, self.text, newNode, None, leafLabel)
        self.tree[newNode] = dict()
        self.tree[newNode][newEdge.StartChar()] = newEdge
        edge2 = Edge(cutIdx, edge.endIdx, self.text, newNode, edge.endNode)
        self.tree[newNode][edge2.StartChar()] = edge2
        self.tree[edge.startNode][edge.StartChar()].endIdx = cutIdx - 1
        self.tree[edge.startNode][edge.StartChar()].endNode = newNode

    def Build(self, root, text):
        l = len(text)
        edge1 = Edge(0, l - 1, text, root)
        self.tree[root] = dict()
        self.tree[root][edge1.StartChar()] = edge1
        for i in range(1, l):
            currNode, currEdge, j, cutIdx = self.Match(i, root, text)
            if not currEdge:
                self.AddEdge(currNode, j, l - 1, i)
            else:
                self.SplitEdge(currEdge, j, l - 1, cutIdx, i)
    
    def ExploreEdges(self, tree):
        results = []
        for node in tree.keys():
            for edge in tree[node].values():
                results.append(edge.Str())
        return results
    
    def PrintTree(self, tree):
        for node in tree.keys():
            print(node.id)
            for edge in tree[node].values():
                if edge.endNode is None:
                    e = None
                else:
                    e = edge.endNode.id
                print(edge.startIdx, edge.endIdx, edge.startNode.id, e, edge.Str())
        print('')
    
    def SaveEdges(self, tree):
        f = open("output.txt", "w")
        f.write('\n'.join(self.ExploreEdges(tree)))
        f.close()
    
    def LongestRepeat(self, tree, text):
        best = (0, 0)
        for node in tree.keys():
            if node.depth > best[1]:
                best = (node.startIdx, node.depth)
        return text[best[0]:best[0] + best[1]]



# Example 1
# ----------
Text = "ATATCGTTTTATCGTT"
SuffixTree(Text + '$') 
# Output: TATCGTT


# Example 2
# ----------
file = open("dataset_30222_5.txt", "r")
with open("dataset_30222_5.txt", "r") as file:
    for line in file:
        data = line.strip()
        
SuffixTree(data + '$')
# Output: CGGACGAATACTGTCGTCGAAGGTAGGGCAGTTCCGGCCCGACCAAACTACTCGGAGGGGCATGTCGTGTATC


# -----------------------------------------------
# Tree Coloring
# -----------------------------------------------

def TreeColoring(adjacency_list, leaf_colors):
    '''Color the internal nodes of a tree given the colors of its leaves.
    
    We color a leaf in this suffix tree blue if it is labeled by the starting 
    position of a suffix starting in Text1; we color a leaf red if it is labeled 
    by the starting position of a suffix starting in Text2.
    
    Input: An adjacency list, followed by color labels for leaf nodes.
    Output: Color labels for all nodes, in any order.'''

    colors = {}
    
    # Assign leaf colors first
    for node, color in leaf_colors.items():
        colors[node] = color
    
    # Initially, all non-leaf nodes are gray
    all_nodes = set(adjacency_list.keys())
    leaf_nodes = set(leaf_colors.keys())
    internal_nodes = all_nodes - leaf_nodes
    for node in internal_nodes:
        colors[node] = 'gray'
    
    def IsRipe(node):
        # A node is ripe if it is gray and has no gray children
        if colors[node] != 'gray':
            return False
        children = adjacency_list[node]
        return all(colors[child] != 'gray' for child in children)
    
    # Repeat while there are ripe nodes
    while True:
        ripe_nodes = [node for node in internal_nodes if IsRipe(node)]
        if not ripe_nodes:
            break
        
        for node in ripe_nodes:
            children_colors = {colors[child] for child in adjacency_list[node]}
            if len(children_colors) > 1:
                colors[node] = 'purple'
            else:
                colors[node] = children_colors.pop()
    
    return colors

def ReadInputFromFile(filename):
    adjacency_list = {}
    leaf_colors = {}

    with open(filename, "r") as file:
        # Read adjacency list until the dash '-'
        for line in file:
            line = line.strip()
            if line == '-':
                break
            if not line:
                continue
            # Parse lines like "2: 0 1"
            node_part, children_part = line.split(':')
            node = int(node_part.strip())
            if children_part.strip() == '':
                children = []
            else:
                children = list(map(int, children_part.strip().split()))
            adjacency_list[node] = children
        
        # Read leaf colors
        for line in file:
            line = line.strip()
            if not line:
                continue
            node_str, color = line.split()
            node = int(node_str)
            leaf_colors[node] = color
    
    return adjacency_list, leaf_colors



# Example 1
# -----------
adjacency_list = {0: [], 1: [], 2: [0, 1], 3: [], 4: [], 5: [2, 3], 6: [], 7: [4, 5, 6],}
leaf_colors = {0: 'red', 1: 'red', 3: 'blue', 4: 'blue', 6: 'red',}

colored_tree = TreeColoring(adjacency_list, leaf_colors)
for node in sorted(colored_tree):
    print(f"{node} {colored_tree[node]}")

# Output:
    # 0 red   1 red   2 red   3 blue   4 blue   5 purple   6 red   7 purple


# Example 2
# -----------
adjacency_list, leaf_colors = ReadInputFromFile("dataset_30233_6.txt")
colored_tree = TreeColoring(adjacency_list, leaf_colors)
with open("output.txt", "w") as file:
    for node in sorted(colored_tree):
        file.write(f"{node} {colored_tree[node]}\n")


# -----------------------------------------------
# Longest Shared Substring
# -----------------------------------------------

class Node:
    def __init__(self):
        self.children = {}  # edge_label_start_char -> child Node
        self.suffix_link = None
        self.start = -1
        self.end = -1
        self.index = -1  # leaf index if leaf
        self.color = 'gray'
        self.parent = None
        
def LongestSharedSubstring(text1, text2):
    '''Find the longest substring shared by two strings.

    Input: Strings Text1 and Text2.
    Output: The longest substring that occurs in both Text1 and Text2.
    Multiple solutions may exist, in which case the function returns one.'''
    
    root, text, len1, len2 = BuildSuffixTrie(text1, text2)
    AssignLeafColors(root, text, len1, len2)

    max_length_info = [0, None]  # [max_length, node]
    FindLongestPurpleNode(root, text, max_length_info)

    if max_length_info[1] is None:
        return ""

    return GetSubstring(max_length_info[1], text)

def BuildSuffixTrie(text1, text2):
    '''Append unique terminators for suffix tree.'''
    
    text = text1 + '#' + text2 + '$'
    root = Node()
    root.start = -1
    root.end = -1

    for i in range(len(text)):
        current = root
        for j in range(i, len(text)):
            c = text[j]
            if c not in current.children:
                new_node = Node()
                new_node.start = j
                new_node.end = j  # only one char per edge
                new_node.parent = current
                current.children[c] = new_node
            current = current.children[c]
        current.index = i  # mark leaf

    return root, text, len(text1), len(text2)

def AssignLeafColors(node, text, len1, len2):
    '''Recursively assign colors to leaves:
       red = suffix from text1 (index < len1)
       blue = suffix from text2 (index > len1)'''
    
    if node.index != -1:
        if node.index < len1:
            node.color = 'red'
        elif node.index > len1:
            node.color = 'blue'
        else:
            # index == len1, terminator, ignore or assign any color
            node.color = 'gray'
        return node.color

    children_colors = set()
    for child in node.children.values():
        child_color = AssignLeafColors(child, text, len1, len2)
        children_colors.add(child_color)

    if len(children_colors) == 1:
        node.color = children_colors.pop()
    else:
        # If children have different colors, mark purple
        if 'red' in children_colors and 'blue' in children_colors:
            node.color = 'purple'
        else:
            node.color = children_colors.pop()

    return node.color

def GetSubstring(node, text):
    '''Reconstruct substring from root to this node.'''
    
    substr = []
    while node.parent is not None:
        substr.append(text[node.start:node.end+1])
        node = node.parent
    return ''.join(reversed(substr))

def FindLongestPurpleNode(node, text, max_length_info, cur_depth=0):
    '''Find deepest purple node (longest shared substring).'''
    
    if node.color == 'purple':
        length = cur_depth
        if length > max_length_info[0]:
            max_length_info[0] = length
            max_length_info[1] = node

    for child in node.children.values():
        edge_length = child.end - child.start + 1
        FindLongestPurpleNode(child, text, max_length_info, cur_depth + edge_length)



# Example 1
# -----------
text1 = "TCGGTAGATTGCGCCCACTC"
text2 = "AGGGGCTCGCAGTGTAAGAA"

longest_substring = LongestSharedSubstring(text1, text2)
print(longest_substring) # Output: TCG


# Example 2
# -----------
with open("dataset_30222_6.txt", "r") as file:
    text1 = file.readline().strip()
    text2 = file.readline().strip()

longest_substring = LongestSharedSubstring(text1, text2)
print(longest_substring) # Output: CGAGTAAACC



# -----------------------------------------------
# Shortest Non-Shared Substring
# -----------------------------------------------

from collections import deque

class Node:
    def __init__(self):
        self.children = {}
        self.parent = None
        self.char_from_parent = ''
        self.depth = 0  # length of substring from root to this node

def ShortestNonSharedSubstring(text1, text2):
    '''Find the shortest substring of one string that does not appear in another string.

    Input: Strings Text1 and Text2.
    Output: The shortest substring of Text1 that does not appear in Text2.
    Multiple solutions may exist, in which case the function returns one.'''
    
    root = BuildSuffixTrie(text1)
    queue = deque([root])

    while queue:
        node = queue.popleft()
        # root corresponds to empty substring, skip it
        if node != root:
            substr = GetSubstringFromNode(node)
            if not SubstringExistsInText(text2, substr):
                return substr
        for child in node.children.values():
            queue.append(child)

    return ""  # fallback if all substrings appear in text2

def BuildSuffixTrie(text):
    root = Node()
    for i in range(len(text)):
        current = root
        for j in range(i, len(text)):
            c = text[j]
            if c not in current.children:
                new_node = Node()
                new_node.parent = current
                new_node.char_from_parent = c
                new_node.depth = current.depth + 1
                current.children[c] = new_node
            current = current.children[c]
    return root

def SubstringExistsInText(text, substring):
    # Check if substring occurs in text (simple search)
    return substring in text

def GetSubstringFromNode(node):
    chars = []
    while node.parent is not None:
        chars.append(node.char_from_parent)
        node = node.parent
    return ''.join(reversed(chars))



# Example 1
# -----------
text1 = "CCAAGCTGCTAGAGG"
text2 = "CATGCTGGGCTGGCT"

result = ShortestNonSharedSubstring(text1, text2)
print(result) # Output: CC


# Example 2
# -----------
with open("dataset_30222_7.txt", "r") as file:
    text1 = file.readline().strip()
    text2 = file.readline().strip()

shortest_substring = ShortestNonSharedSubstring(text1, text2)
print(shortest_substring) # Output: AAATT

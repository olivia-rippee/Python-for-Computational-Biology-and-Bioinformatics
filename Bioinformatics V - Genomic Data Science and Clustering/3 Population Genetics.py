# -----------------------------------------------
# Mitochondrial Survival Over Generations
# -----------------------------------------------

def ExpectedMitochondrialSurvival(generations):
    '''Returns the expected fraction of original mitochondrial genomes
    that remain after a given number of generations, assuming every woman 
    in this population has exactly two children.
    
    If a woman has no daughters, her mtDNA is lost. Otherwise, it survives.
    This is a 1/4 probability of loss, or 3/4 probability of survival.

    Input: Number of generations.
    Output: Expected fraction (rounded to 3 decimal places).
    '''
    
    survival_fraction = (3 / 4) ** generations
    return round(survival_fraction, 3)


# Example
# ---------
ExpectedMitochondrialSurvival(generations = 6) # 0.178



# -----------------------------------------------
# Difference Between SNPs
# -----------------------------------------------

def Diff(s, t):
    '''Compute the Diff(s, t) value between two single-nucleotide
    polymorphism vectors s and t.

    We want to determine how well a SNP s explains another SNP t. We average 
    this behavior over all pairs of individuals.

    Diff(s, t) is defined as the ratio:
        (# of pairs (i, j) such that s[i] ≠ s[j] and t[i] ≠ t[j]) /
        (# of pairs (i, j) such that t[i] ≠ t[j])

    The larger the value of Diff(s,t), the more that any variance exhibited in t can be
    explained by s. If Diff(s,t) = 1, then we say that s is fully informative for t.         

    Input: A binary vector representing SNP values and a binary vector representing 
    target or group labels. Must be the same length.

    Output: A decimal value between 0 and 1 (rounded to 3 decimal places) representing
    how often SNP differences correspond with group differences.'''
    
    if len(s) != len(t):
        raise ValueError("s and t must be the same length.")
    
    n = len(s)
    numerator = 0
    denominator = 0

    for i in range(n):
        for j in range(i + 1, n):
            if t[i] != t[j]:
                denominator += 1
                if s[i] != s[j]:
                    numerator += 1

    if denominator == 0:
        return 0.0  # Avoid division by zero

    return round(numerator / denominator, 3)


# Example 1
# ------------
s = (1, 0, 0, 0, 0)
t = (0, 1, 0, 1, 1)
Diff(s,t) # Output: 0.5


# Example 2
# ------------
s = (0, 0, 1, 1, 0, 0, 1, 0) 
t = (1, 1, 0, 0, 1, 1, 1, 1)
Diff(s,t) # Output: 0.833



# -----------------------------------------------
# Quantify how well a collection of SNPs explains a single SNP
# -----------------------------------------------

from itertools import combinations

def ExtendedDiff(S_prime, t):
    '''Computes Diff(S', t) where S' is a list of SNP vectors, and t is a target SNP.

    Input: SNPs in S' and target SNP
    Output: Diff(S', t), rounded to 3 decimal places.'''
    
    n = len(t)
    numerator = 0
    denominator = 0

    for i, j in combinations(range(n), 2):
        if t[i] != t[j]:
            denominator += 1
            if any(s[i] != s[j] for s in S_prime):
                numerator += 1

    if denominator == 0:
        return 0.0

    return round(numerator / denominator, 3)


# Example
# ----------
S_prime = [(1, 0, 0, 0, 0), (1, 1, 1, 0, 0)]
t = (0, 1, 0, 1, 1)
print(ExtendedDiff(S_prime, t)) # Output: 0.833



# -----------------------------------------------
# k-Most Informative SNPs Heuristic
# -----------------------------------------------

import random
from itertools import combinations

def Diff(S_prime, t):
    '''Compute Diff(S', t) as per the definition:
    Numerator: #pairs (i,j) with t_i != t_j and for some s in S', s_i != s_j
    Denominator: #pairs (i,j) with t_i != t_j'''
    
    n = len(t)
    numerator = 0
    denominator = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            if t[i] != t[j]:
                denominator += 1
                if any(s[i] != s[j] for s in S_prime):
                    numerator += 1
    
    if denominator == 0:
        return 0.0
    
    return numerator / denominator

def RandomizedHaplotypeSearch(S, T, k, max_iter=1000):
    '''Identify the best subset of k SNPs explaining T.
    
    Input: SNP matrices corresponding to two collections of SNPs S and T, along with an integer k.
        max_iter (int): Maximum iterations to avoid infinite loops.
    Output: Best subset of k SNPs from S minimizing Diff(S', T).'''

    # Start with random k SNPs
    bestSNPs = random.sample(S, k)
    
    for _ in range(max_iter):
        currentSNPs = bestSNPs.copy()
        improved = False
        
        # SNPs not in current subset
        not_in_current = [s for s in S if s not in currentSNPs]
        
        for i, s in enumerate(currentSNPs):
            for s_prime in not_in_current:
                new_subset = currentSNPs.copy()
                new_subset[i] = s_prime
                
                if Diff(new_subset, T) < Diff(bestSNPs, T):
                    bestSNPs = new_subset
                    improved = True
                    break  # break inner loop to restart improvement search
            if improved:
                break  # break outer loop if improved
        
        if not improved:
            # No improvement, return result
            return bestSNPs
    
    # If max_iter reached without convergence
    return bestSNPs


# Example
# ---------
S = [(1, 0, 0, 0, 0),
     (1, 1, 1, 0, 0),
     (0, 1, 0, 1, 0),
     (1, 1, 0, 0, 1),
     (0, 0, 1, 1, 1)]
T = (0, 1, 0, 1, 1)
k = 2

best_subset = RandomizedHaplotypeSearch(S, T, k)
print("Best subset of SNPs:", best_subset)
print("Diff value:", Diff(best_subset, T))

# Output:
    # Best subset of SNPs: [(1, 0, 0, 0, 0), (0, 0, 1, 1, 1)]
    # Diff value: 0.6666666666666666


# -----------------------------------------------
# Compatible Columns
# -----------------------------------------------

from itertools import product

def CountCompatibleColumns(v):
    '''Counts how many binary columns of the same length as v are compatible 
    with v. Two binary columns are compatible if their 1s are either subset, 
    superset, or disjoint with respect to each other.

    Input: A binary tuple representing the reference column (e.g. (1, 0, 1, 0, 1)).
    Output: Number of compatible binary columns.'''
    
    n = len(v)
    Ov = {i for i, val in enumerate(v) if val == 1}
    compatible_count = 0

    for w in product([0, 1], repeat=n):
        Ow = {i for i, val in enumerate(w) if val == 1}

        if Ow.issubset(Ov) or Ov.issubset(Ow) or Ov.isdisjoint(Ow):
            compatible_count += 1

    return compatible_count


# Example
# ---------
v = (1, 0, 1, 0, 1)
print(CountCompatibleColumns(v))  # Output: 14



# -----------------------------------------------
# Perfect Phylogeny
# -----------------------------------------------

from collections import defaultdict
from typing import List, Dict

def PerfectPhylogeny(matrix):
    '''Reconstruct a perfect phylogeny from a binary character matrix.

    Perfect phylogeny:
        1. The SNP vector at the root consists of only zeros (as it corresponds to an ancestor
        possessing none of the characters corresponding to the SNPs).
        2. Each row of A is the SNP vector of exactly one leaf.
        3. For any column j of A, there is a single node v such that every node in Tv , the
        subtree of T rooted at v, contains a 1 at the j-th position, and every other node in
        T contains a 0 at the j-th position.
    
    If a tree satisfying the perfect phylogeny assumption exists for a given matrix A, then we say 
    that the tree perfectly fits A, and we call A phylogenetic.
    
    Input: A binary character matrix.
    Output: A tree that fits the matrix perfectly, if such a tree exists.'''

    if not AllColumnsCompatible(matrix):
        raise ValueError("Input matrix is not phylogenetic: incompatible column pair found.")

    return BuildCharacterTree(matrix)

def BuildCharacterTree(matrix):
    '''Construct a tree from the binary character matrix.
    Returns a nested dictionary representing the tree.'''
    
    n = len(matrix)
    m = len(matrix[0])
    
    # Build a character-to-indices map
    char_to_indices = defaultdict(set)
    for j in range(m):
        for i in range(n):
            if matrix[i][j] == 1:
                char_to_indices[j].add(i)

    # Sort characters by inclusion (subset) relationship
    characters = list(range(m))
    characters.sort(key=lambda c: len(char_to_indices[c]))

    # Build tree: root starts with all zeros
    tree = {"Id": "Root", "Vector": [0]*m, "Children": []}
    node_map = {"": tree}  # path as key

    for c in characters:
        added = False
        for path, node in list(node_map.items()):
            node_indices = set()
            if path == "":
                node_indices = set(range(n))
            else:
                path_chars = list(map(int, path.strip().split("-")))
                node_indices = set.intersection(*(char_to_indices[ch] for ch in path_chars))

            if char_to_indices[c].issubset(node_indices):
                new_path = path + f"-{c}" if path else f"{c}"
                new_node = {
                    "Id": f"Char_{c}",
                    "Vector": [1 if i in list(map(int, new_path.split("-"))) else 0 for i in range(m)],
                    "Children": []
                }
                node["Children"].append(new_node)
                node_map[new_path] = new_node
                added = True
                break
        if not added:
            raise ValueError("Could not insert character into tree — matrix might not be compatible.")

    return tree

def IsCompatible(col1, col2):
    seen = set()
    for a, b in zip(col1, col2):
        seen.add((a, b))
        if len(seen) == 4 or ((1, 0) in seen and (0, 1) in seen and (1, 1) in seen):
            return False
    return True

def AllColumnsCompatible(matrix):
    num_cols = len(matrix[0])
    for i in range(num_cols):
        for j in range(i + 1, num_cols):
            col_i = [row[i] for row in matrix]
            col_j = [row[j] for row in matrix]
            if not IsCompatible(col_i, col_j):
                return False
    return True


# Example
# --------
matrix = [[1, 0, 0],
          [1, 1, 0],
          [1, 1, 1]]
tree = PerfectPhylogeny(matrix)
import json
print(json.dumps(tree, indent=2))

# Output: {
    #  "Id": "Root",
    #  "Vector": [
    #    0,
    #    0,
    #    0
    #  ],
    #  "Children": [
    #    {
    #      "Id": "Char_2",
    #      "Vector": [
    #        0,
    #        0,
    #        1
    #      ],
    #      "Children": []
    #    },
    #    {
    #      "Id": "Char_1",
    #      "Vector": [
    #        0,
    #        1,
    #        0
    #      ],
    #      "Children": []
    #    },
    #    {
    #      "Id": "Char_0",
    #      "Vector": [
    #        1,
    #        0,
    #        0
    #      ],
    #      "Children": []
    #    }
    #  ]
    # }



# -----------------------------------------------
# Augmented Perfect Phylogeny
# -----------------------------------------------
# Augment tree with the data from a newly genotyped individual

class TreeNode:
    def __init__(self, mutation=None):
        self.mutation = mutation  # mutation index for the edge leading here
        self.children = []
        self.parent = None
        self.label = set()  # mutations along the path from root to this node
        self.is_leaf = False
        self.individual_ids = []  # store individuals assigned to this leaf

def BuildLabels(root, parentLabel=set()):
    root.label = set(parentLabel)
    if root.mutation is not None:
        root.label.add(root.mutation)
    for child in root.children:
        BuildLabels(child, root.label)

def PrintTree(node, indent=0):
    prefix = "  " * indent
    if node.mutation is None:
        label = "root"
    else:
        label = f"Mutation {node.mutation}"
    if node.is_leaf:
        label += " (leaf: " + ", ".join(node.individual_ids) + ")"
    print(prefix + label)
    for child in node.children:
        PrintTree(child, indent + 1)

def FourGameteTest(A):
    # Checks if matrix A admits a perfect phylogeny (no recombination)
    n, m = len(A), len(A[0])
    for i in range(m):
        for j in range(i+1, m):
            patterns = set()
            for k in range(n):
                patterns.add((A[k][i], A[k][j]))
            if len(patterns) == 4:
                return False  # four gametes present → no perfect phylogeny
    return True

def BuildPerfectPhylogenyTree(A, individual_ids):
    '''
    Build a perfect phylogeny tree T from SNP matrix A and list of individual_ids.
    Returns root node of the tree.
    Assumes A admits a perfect phylogeny (no conflicts).
    '''

    if not FourGameteTest(A):
        raise ValueError("Matrix does not admit a perfect phylogeny")

    root = TreeNode()
    BuildLabels(root, set())

    # Insert each individual as a path of mutations
    for idx, row in enumerate(A):
        current = root
        # mutations that this individual has
        muts = [i for i, val in enumerate(row) if val == 1]

        # For each mutation, try to find or create a child with that mutation
        for mut in muts:
            # Look for child with this mutation
            child = next((c for c in current.children if c.mutation == mut), None)
            if child is None:
                # Create new node
                child = TreeNode(mut)
                child.parent = current
                current.children.append(child)
            current = child
        # At the end of path, mark leaf and add individual id
        current.is_leaf = True
        current.individual_ids.append(individual_ids[idx])

    BuildLabels(root)
    return root

def FindInsertionPoint(root, c):
    stack = [root]
    insertionNode = root
    while stack:
        node = stack.pop()
        if node.label.issubset(c):
            insertionNode = node
            for child in node.children:
                stack.append(child)
    return insertionNode

def ViolatesUniqueness(root, mutations, insertionNode):
    used = set()

    def Dfs(node):
        if node != insertionNode and node.mutation is not None:
            used.add(node.mutation)
        for child in node.children:
            Dfs(child)

    Dfs(root)
    return any(m in used for m in mutations)

def InsertIndividual(root, c, individualId):
    BuildLabels(root)
    cMutations = set(i for i, v in enumerate(c) if v == 1)
    insertionNode = FindInsertionPoint(root, cMutations)
    newMutations = cMutations - insertionNode.label

    if ViolatesUniqueness(root, newMutations, insertionNode):
        print(f"Cannot insert individual {individualId} — violates perfect phylogeny.")
        return False

    current = insertionNode
    for mutation in newMutations:
        newNode = TreeNode(mutation)
        newNode.parent = current
        current.children.append(newNode)
        current = newNode

    leaf = TreeNode()
    leaf.is_leaf = True
    leaf.individual_ids.append(individualId)
    leaf.parent = current
    current.children.append(leaf)

    print(f"Inserted individual {individualId} as a leaf.")
    return True



# Example
# ----------
A = [[0,0,0],  # individual A
     [1,0,0],  # individual B
     [1,1,0],] # individual C
individual_ids = ["A", "B", "C"]

print("Building perfect phylogeny tree from matrix A...")
root = BuildPerfectPhylogenyTree(A, individual_ids)
print("\nInitial tree:")
PrintTree(root)

# Output:
    # Initial tree:
    # root (leaf: A)
      # Mutation 0 (leaf: B)
        # Mutation 1 (leaf: C)



# Add new individual 
c = [1,1,0]  # mutations 0 and 1
print("\nInserting new individual D with SNP vector:", c)
InsertIndividual(root, c, individualId="D")

print("\nAugmented tree:")
PrintTree(root)

# Output: 
    # Inserting new individual D with SNP vector: [1, 1, 0]
    # Inserted individual D as a leaf.

    # Augmented tree:
    # root (leaf: A)
      # Mutation 0 (leaf: B)
        # Mutation 1 (leaf: C)
         # root (leaf: D)
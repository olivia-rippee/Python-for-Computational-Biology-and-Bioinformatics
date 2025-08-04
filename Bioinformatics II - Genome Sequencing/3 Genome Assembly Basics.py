import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics II - Genome Sequencing/Data")



# -----------------------------------------------
# String Composition
# -----------------------------------------------

def Composition(k, Text):
    '''Generate the k-mer composition of a string
    
    Input: An integer k and a string Text.
    Output: Composition(Text).
    '''
    
    kmers = []
    n = len(Text)
    i = 0
    while i <= (n - k):
        kmers.append(Text[i:i+k])
        i += 1
    return kmers



# Examples
# ------------
k=5
Text="CAATCCAAC"
print(*Composition(k, Text))
# Output: CAATC AATCC ATCCA TCCAA CCAAC


with open("dataset_30153_3.txt", "r") as file:
    k = int(file.readline())
    Text = file.readline().strip()
    
with open(r"output.txt", 'w') as file:
       result = Composition(k, Text)
       for item in result:
            file.write(f"{item} ")



            
            
# -----------------------------------------------
# String Spelled by a Genome Path
# -----------------------------------------------     

def PathToGenome(path):
    '''Reconstruct the genome from its genome path.
    
    Input: Genome path
    Output: Genome
    '''
    
    Text = path[0]
    for i in range(1, len(path)):
        Text += path[i][-1]
    return Text



# Examples
# ------------------
path = ['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT']
print(PathToGenome(path))
# Output: ACCGAAGCT


with open("dataset_30182_3.txt", "r") as file:
    content = file.read()
    path = content.strip().split()
result = PathToGenome(path)
with open(r"output.txt", 'w') as file:
       file.write(result)
# Output: CCTACTACACCCCTTACACATAGCCATCTGAGAGTTGGCCGTCC...



# -----------------------------------------------
# Overlap Graph Problem
# -----------------------------------------------     

def OverlapGraph(patterns):
    '''Construct the overlap graph of a collection of k-mers.
    
    Input: A collection patterns of k-mers.
    Output: The overlap graph Overlap(patterns).
    '''

    patterns = patterns.split()
    adjacency_list = {}    # Create an entry for each pattern in the adjacency list
    for pattern in patterns:
        adjacency_list[pattern] = []
    for pattern in patterns:
        suffix = pattern[1:]  # Suffix of the current pattern
        for other_pattern in patterns:
            if pattern != other_pattern and other_pattern.startswith(suffix):
                adjacency_list[pattern].append(other_pattern)
    return adjacency_list
    
    
    
# Example
# ------------
with open("dataset_30182_10.txt", "r") as file:
    patterns = file.read().strip()
    
overlap_graph = OverlapGraph(patterns)

with open("output.txt", "w") as out_file:
    for node, edges in overlap_graph.items():
        if edges:  # Only write nodes with edges
            for edge in edges:
                out_file.write(f"{node}: {edge}\n")
                
# Output: CACTGGGGTACTGATCGGGGTGCTA: ACTGGGGTACTGATCGGGGTGCTAG
# AACCGCCCACATTTTTTGGCACATA: ACCGCCCACATTTTTTGGCACATAC
# CGGCATCGAGTCCGTGGAGATAATT: GGCATCGAGTCCGTGGAGATAATTG
# ATTAGCCATGAAGGAGCATTCCTCC: TTAGCCATGAAGGAGCATTCCTCCT etc.



# -----------------------------------------------
# De Bruijn Graph from a String
# -----------------------------------------------

from collections import defaultdict

def Composition(k, text):
    result = []
    for i in range(0, len(text) - k + 1):
        chunk = text[i:i + k]
        result.append(chunk)
    return result

def DeBruijnGraph(k, text):
    '''Construct de Bruijn graph by gluing identically labeled nodes in PathGraph.
    
    Input: An integer k and a string Text.
    Output: DeBruijn_k(Text).
    '''
    
    composition = Composition(k, text)
    dic = defaultdict(list)
    for i in composition:
        prefix = i[:k - 1]
        sufix = i[1:]
        dic[prefix].append(sufix)
    return dic


# Examples
# ------------
k=4
Text = "AAGATTCTCTAAGA"
deBruijn_graph = DeBruijnGraph(k, Text)
print(deBruijn_graph)

with open("output.txt", "w") as out_file:
    for node, edges in deBruijn_graph.items():
        if edges:  # Only write nodes with edges
            if len(edges) > 1:
                out_file.write(f"{node}: {' '.join(edges)}\n")
            else:
                out_file.write(f"{node}: {edges[0]}\n")
# Output: 
# AAG: AGA AGA
# AGA: GAT
# GAT: ATT
# ATT: TTC
# TTC: TCT
# TCT: CTC CTA
# CTC: TCT
# CTA: TAA
# TAA: AAG



with open("dataset_30183_6.txt", "r") as file:
    k = int(file.readline())
    Text = file.read().strip()
    
deBruijn_graph = DeBruijnGraph(k, Text)

with open("output.txt", "w") as out_file:
    for node, edges in deBruijn_graph.items():
        if edges:  # Only write nodes with edges
            if len(edges) > 1:
                out_file.write(f"{node}: {' '.join(edges)}\n")
            else:
                out_file.write(f"{node}: {edges[0]}\n")
                
# Output: 
# TGAGGCTTCTT: GAGGCTTCTTC
# GAGGCTTCTTC: AGGCTTCTTCA
# AGGCTTCTTCA: GGCTTCTTCAT
# GGCTTCTTCAT: GCTTCTTCATC
# GCTTCTTCATC: CTTCTTCATCG
# CTTCTTCATCG: TTCTTCATCGG etc.



# -----------------------------------------------
# DeBruijn Graph from k-mers
# -----------------------------------------------

def DeBruijnGraphFromKmers(kmers):
    '''Generate the deBruijn graph of a list of k-mers.
    Input: A collection of k-mers patterns.
    Output: The adjacency list of the de Bruijn graph DeBruijn(patterns).
    '''
    
    adj_list = {}
    for pattern in kmers:
        prefix = pattern[:-1]
        suffix = pattern[1:]
        if prefix in adj_list:
            adj_list[prefix].append(suffix)
        else:
            adj_list[prefix] = [suffix]
    return adj_list


# Examples
# -----------
kmers = ("GAGG", "CAGG", "GGGG", "GGGA", "CAGG", "AGGG", "GGAG")
deBruijn_graph = DeBruijnGraphFromKmers(kmers)

with open("output.txt", "w") as out_file:
    for node, edges in deBruijn_graph.items():
            if len(edges) > 1:
                out_file.write(f"{node}: {' '.join(edges)}\n")
            else:
                out_file.write(f"{node}: {edges[0]}\n")

# Output: 
# AGG: GGG
# CAG: AGG AGG
# GAG: AGG
# GGA: GAG
# GGG: GGA GGG




with open("dataset_30184_8.txt", "r") as file:
    kmers = file.readline().strip().split()
deBruijn_graph = DeBruijnGraphFromKmers(kmers)

with open("output.txt", "w") as out_file:
    for node, edges in deBruijn_graph.items():
        if len(edges) > 1:
            out_file.write(f"{node}: {' '.join(edges)}\n")
        else:
            out_file.write(f"{node}: {edges[0]}\n")

# Output:
# TCGATCAGCCCCTCGGTGT: CGATCAGCCCCTCGGTGTC
# GACAAAGCCGTCCCTCACA: ACAAAGCCGTCCCTCACAG
# GAGCACCACTTGTATGGGA: AGCACCACTTGTATGGGAT
# CTACACGACAAAGCCGTCC: TACACGACAAAGCCGTCCC TACACGACAAAGCCGTCCC TACACGACAAAGCCGTCCC
# AGCAGACCATCTGACCCAG: GCAGACCATCTGACCCAGT etc.

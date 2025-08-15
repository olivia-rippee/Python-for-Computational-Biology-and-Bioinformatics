import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics IV - Molecular Evolution/Data")


# -----------------------------------------------
# Graph of Spectrum
# -----------------------------------------------

def ConstructSpectrumGraph(spectrum):
    '''Construct the graph of a spectrum.

    Input: A space-delimited list of integers Spectrum.
    Output: Graph(Spectrum).'''
    
    spectrum = sorted([0] + spectrum)  # Include starting mass 0
    edges = []
    
    AA_mass_table = {}
    with open("amino_acid_integer_mass_table.txt", "r") as file:
        for line in file:
            aa, mass = line.strip().split()
            AA_mass_table[int(mass)] = aa


    for i in range(len(spectrum)):
        for j in range(i + 1, len(spectrum)):
            diff = spectrum[j] - spectrum[i]
            if diff in AA_mass_table:
                aa = AA_mass_table[diff]
                edges.append(f"{spectrum[i]}->{spectrum[j]}:{aa}")

    return edges


def ReadSpectrumFromFile(filename):
    with open(filename, "r") as file:
        line = file.readline().strip()
        return list(map(int, line.split()))

def WriteGraphToFile(edges, filename):
    with open(filename, "w") as file:
        for edge in edges:
            file.write(edge + '\n')



# Examples
# ----------
spectrumInput = "57 71 154 185 301 332 415 429 486"
spectrum = list(map(int, spectrumInput.split()))
edges = ConstructSpectrumGraph(spectrum)

for edge in edges:
    print(edge)

# Output: 0->57:G  0->71:A  57->154:P  57->185:Q  71->185:N  154->301:F 185->332:F  
    # 301->415:N  301->429:Q  332->429:P  415->486:A  429->486:G
        
        

spectrum = ReadSpectrumFromFile("dataset_30262_5.txt")
edges = ConstructSpectrumGraph(spectrum)
WriteGraphToFile(edges, "output.txt")
        
# Output: 0->101:T  0->114:N  101->214:L  114->227:L  214->285:A  227->390:Y
    # 285->384:V  384->513:E  390->546:R  513->610:P  546->674:Q  610->707:P ...



# -----------------------------------------------
# Decode Ideal Spectrum
# -----------------------------------------------

massTable = {}
with open("amino_acid_integer_mass_table.txt", "r") as file:
    for line in file:
        aa, mass = line.strip().split()
        massTable[int(mass)] = aa

def ConstructSpectrumGraph(spectrum):
    spectrum = sorted([0] + spectrum)
    graph = {}  # dict: node -> list of (neighbor, aa)
    
    for i in range(len(spectrum)):
        for j in range(i + 1, len(spectrum)):
            diff = spectrum[j] - spectrum[i]
            if diff in massTable:
                aa = massTable[diff]
                u = spectrum[i]
                v = spectrum[j]
                if u not in graph:
                    graph[u] = []
                graph[u].append((v, aa))
    return graph

def IdealSpectrum(peptide, massTable):
    prefixMass = [0]
    for aa in peptide:
        prefixMass.append(prefixMass[-1] + massTable[aa])
    spectrum = set(prefixMass)
    peptideLength = len(peptide)
    for i in range(peptideLength):
        for j in range(i + 1, peptideLength + 1):
            subpeptideMass = prefixMass[j] - prefixMass[i]
            spectrum.add(subpeptideMass)
    return sorted(spectrum)

def DFS(graph, current, target, path, peptide, allPaths):
    if current == target:
        allPaths.append(''.join(peptide))
        return
    if current not in graph:
        return
    for neighbor, aa in graph[current]:
        DFS(graph, neighbor, target, path + [neighbor], peptide + [aa], allPaths)

def DecodeIdealSpectrum(spectrum):
    '''
    Input: A space-delimited list of integers Spectrum.
    Output: An amino acid string that explains Spectrum.'''
    
    start = 0
    end = max(spectrum)
    graph = ConstructSpectrumGraph(spectrum)
    allPaths = []
    DFS(graph, start, end, [], [], allPaths)
    for peptide in allPaths:
        return peptide

def ReadSpectrumFromFile(filename):
    with open(filename, "r") as file:
        line = file.readline().strip()
        return sorted([int(x) for x in line.split()])



# Example 1
# ----------
spectrum_str = "57 71 154 185 301 332 415 429 486"
spectrum = list(map(int, spectrum_str.split()))
print(DecodeIdealSpectrum(spectrum))
# Output: GPFNA


# Example 2
# ----------
spectrum = ReadSpectrumFromFile("dataset_30262_8.txt")
peptide = DecodeIdealSpectrum(spectrum)
print(peptide)
# Output: LGGDLTQSPTFRYYSTHMYFHQHL



# -----------------------------------------------
# Convert Peptide to Peptide Vector
# -----------------------------------------------

def ConvertPeptideToVector(Peptide, massTable):
    ''' Convert a peptide into a peptide vector.

    Input: An amino acid string P.
    Output: The peptide vector of P (in the form of space-separated integers).'''
    
    # Calculate prefix masses
    prefix_masses = []
    current_mass = 0
    for amino_acid in Peptide:
        current_mass += massTable[amino_acid]
        prefix_masses.append(current_mass)
    
    total_mass = prefix_masses[-1]
    
    # Initialize the peptide vector
    peptide_vector = [0] * total_mass
    
    # Set positions corresponding to prefix masses to 1
    for mass in prefix_masses:
        peptide_vector[mass - 1] = 1  # -1 because list is 0-indexed
    
    return ' '.join(map(str, peptide_vector))


# Example 1
# -----------
massTable = {"X":4, "Z":5}

peptide = "XZZXX"
vector = ConvertPeptideToVector(peptide, massTable)
print(vector)  # Output: 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0 0 0 1



# Example 2
# -----------
massTable = {"G":57, "A":71, "S":87, "P":97, "V":99, "T":101, "C":103, 
             "I":113, "L":113, "N":114, "D":115, "K":128, "Q":128, "E":129, 
             "M":131, "H":137, "F":147, "R":156, "Y":163, "W":186}

with open("dataset_30264_5.txt", "r") as file:
    peptide = file.readline().strip()
vector = ConvertPeptideToVector(peptide, massTable)

with open("output.txt", "w") as file:
    file.write(vector)

# Output: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    # 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    # 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    # 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    # 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    


# -----------------------------------------------
# Convert Peptide Vector to Peptide
# -----------------------------------------------

def ConvertVectorToPeptide(Vector, massTable):
    '''Convert a peptide vector into a peptide.
    
    Input: A space-delimited binary vector P.
    Output: An amino acid string whose binary peptide vector matches P.'''
    
    # Step 1: Parse binary vector string into list of ints
    vector = list(map(int, Vector.strip().split()))
    
    # Step 2: Find prefix masses (1-based index)
    prefix_masses = [i + 1 for i, val in enumerate(vector) if val == 1]
    
    # Step 3: Get differences between prefix masses
    peptide_masses = [prefix_masses[0]]  # First mass is from 0 to first prefix
    for i in range(1, len(prefix_masses)):
        peptide_masses.append(prefix_masses[i] - prefix_masses[i - 1])
    
    # Step 4: Convert masses to amino acids
    peptide = ""
    for mass in peptide_masses:
        if mass not in massTable:
            raise ValueError(f"No amino acid with mass {mass}")
        peptide += massTable[mass][0]  # Pick the first valid amino acid
    
    return peptide


# Example 1
# ----------
massTable = {4:"X", 5:"Z"}
Vector = "0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0 0 0 1"
peptide = ConvertVectorToPeptide(Vector, massTable)
print(peptide)  # Output: XZZXX



# Example 2
# ----------
massTable = {57:"G", 71:"A", 87:"S", 97:"P", 99:"V", 101:"T", 103:"C", 
             113:"L", 114:"N", 115:"D", 128:"Q", 129:"E", 131:"M", 
             137:"H", 147:"F", 156:"R", 163:"Y", 186:"W"}

with open("dataset_30264_6.txt", "r") as file:
    Vector = file.readline().strip()
    
peptide = ConvertVectorToPeptide(Vector, massTable)
print(peptide)  #Output: HQPCRWNDRSVSLLEAWRYVGMPQDNALDFVVC




# -----------------------------------------------
# Peptide Sequencing
# -----------------------------------------------

def PeptideSequencing(SpectralVector, massTable):
    '''Given a spectral vector Spectrum' = s1 ... sm, the goal is to find a peptide whose binary peptide vector 
       maximizes the dot product with Spectrum'. The peptide vector must have the same length as Spectrum', and 
       the peptide's mass must match the length of the spectrum. We model this as a weighted DAG path problem 
       from node 0 to node m, where each node represents a prefix mass, and edges (i → j) exist if (j - i) is 
       a valid amino acid mass. The weight of a node is the spectral score at that position.
       
       Input: A space-delimited spectral vector Spectrum'.
       Output: An amino acid string with maximum score against Spectrum'. For masses with more than one amino acid, any choice may be used.
       '''

    # Parse spectral vector
    spectrum = list(map(int, SpectralVector.strip().split()))
    spectrum = [0] + spectrum  # s0 = 0
    
    n = len(spectrum)
    dp = [float('-inf')] * n      # dp[i] = best score to reach node i
    backtrack = [None] * n        # backtrack[i] = (previous_node, amino_acid)
    
    dp[0] = 0  # start at node 0

    # Dynamic programming over all nodes
    for i in range(1, n):
        for mass in massTable:
            j = i - mass
            if j >= 0:
                score = dp[j] + spectrum[i]
                if score > dp[i]:
                    dp[i] = score
                    backtrack[i] = (j, massTable[mass][0])  # pick any valid amino acid

    # Reconstruct peptide from backtrack
    peptide = []
    i = n - 1
    while i != 0:
        j, aa = backtrack[i]
        peptide.append(aa)
        i = j

    return ''.join(reversed(peptide))


# Example 1
# ----------
massTable = {4:"X", 5:"Z"}
SpectralVector = "0 0 0 4 -2 -3 -1 -7 6 5 3 2 1 9 3 -8 0 3 1 2 1 8"
peptide = PeptideSequencing(SpectralVector, massTable)
print(peptide)  # Output: XZZXX



# Example 2
# ----------
massTable = {57:"G", 71:"A", 87:"S", 97:"P", 99:"V", 101:"T", 103:"C", 
             113:"L", 114:"N", 115:"D", 128:"Q", 129:"E", 131:"M", 
             137:"H", 147:"F", 156:"R", 163:"Y", 186:"W"}

with open("dataset_30264_13.txt", "r") as file:
    SpectralVector = file.readline().strip()

peptide = PeptideSequencing(SpectralVector, massTable)
print(peptide)  # Output: SGSVGGSGDAP

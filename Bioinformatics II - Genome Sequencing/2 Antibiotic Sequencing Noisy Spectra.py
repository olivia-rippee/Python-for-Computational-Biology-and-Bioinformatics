import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics II - Genome Sequencing/Data")


# -----------------------------------------------
# Cyclopeptide Scoring
# -----------------------------------------------

from collections import Counter

# integer mass table
AminoAcidMass = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97,
    'V': 99, 'T': 101, 'C': 103, 'I': 113,
    'L': 113, 'N': 114, 'D': 115, 'K': 128,
    'Q': 128, 'E': 129, 'M': 131, 'H': 137,
    'F': 147, 'R': 156, 'Y': 163, 'W': 186
}

def CyclicSpectrum(Peptide):
    '''Generate the cyclic spectrum of a peptide.'''
    
    PrefixMass = [0]
    for AminoAcid in Peptide:
        PrefixMass.append(PrefixMass[-1] + AminoAcidMass[AminoAcid])
    
    PeptideMass = PrefixMass[-1]
    Spectrum = [0]

    for Start in range(len(Peptide)):
        for End in range(Start + 1, len(Peptide) + 1):
            Mass = PrefixMass[End] - PrefixMass[Start]
            Spectrum.append(Mass)
            if Start > 0 and End < len(Peptide):
                CyclicMass = PeptideMass - Mass
                Spectrum.append(CyclicMass)
    
    return sorted(Spectrum)

def Score(Peptide, Spectrum):
    '''Compute the score of a cyclic peptide against a spectrum.
    Accounts for false or missing masses.'''
    
    TheoreticalSpectrum = Counter(CyclicSpectrum(Peptide))
    ExperimentalSpectrum = Counter(Spectrum)
    Score = 0

    for Mass in TheoreticalSpectrum:
        Score += min(TheoreticalSpectrum[Mass], ExperimentalSpectrum.get(Mass, 0))
    
    return Score



# Examples 
# ---------
Peptide = "NQEL"
Spectrum = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
Score(Peptide, Spectrum)
# Output: 11


with open("dataset_30244_3.txt", "r") as file:
    lines = file.read().strip().split('\n')
    Peptide = lines[0].strip()
    Spectrum = list(map(int, lines[1].strip().split()))
Score(Peptide, Spectrum)
# Output: 621


Peptide = "MAMA" 
Spectrum = [0, 71, 98, 99, 131, 202, 202, 202, 202, 202, 299, 333, 333, 333, 503]
Score(Peptide, Spectrum)
# Output: 9


# -----------------------------------------------
# Trim candidate peptides to top N
# -----------------------------------------------

from collections import Counter

def LinearSpectrum(Peptide):
    '''Generate the linear spectrum of a peptide (no wrap-around).'''
    
    PrefixMass = [0]
    for AminoAcid in Peptide:
        PrefixMass.append(PrefixMass[-1] + AminoAcidMass[AminoAcid])
    
    Spectrum = [0]
    for Start in range(len(Peptide)):
        for End in range(Start + 1, len(Peptide) + 1):
            Mass = PrefixMass[End] - PrefixMass[Start]
            Spectrum.append(Mass)
    
    return sorted(Spectrum)

def LinearScore(Peptide, Spectrum):
    '''Compute the linear score of a peptide with respect to a spectrum.'''
    
    TheoreticalSpectrum = Counter(LinearSpectrum(Peptide))
    ExperimentalSpectrum = Counter(Spectrum)
    Score = 0

    for Mass in TheoreticalSpectrum:
        Score += min(TheoreticalSpectrum[Mass], ExperimentalSpectrum.get(Mass, 0))
    
    return Score


# Examples
# --------------------
Peptide = "NQEL"
Spectrum = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
LinearScore(Peptide, Spectrum)
# Output: 8


with open("dataset_30249_1.txt", "r") as file:
    lines = file.read().strip().split('\n')
    Peptide = lines[0].strip()
    Spectrum = list(map(int, lines[1].strip().split()))
LinearScore(Peptide, Spectrum)
# Output: 294


Peptide = "PEEP" 
Spectrum = [0, 97, 97, 129, 194, 196, 226, 226, 244, 258, 323, 323, 452]
LinearScore(Peptide, Spectrum)
# Output: 8


# -----------------------------------------------
# Cyclopeptide Sequencing for Spectra with Errors
# -----------------------------------------------

from collections import Counter

def LinearSpectrum(Peptide):
    PrefixMass = [0]
    for AminoAcid in Peptide:
        PrefixMass.append(PrefixMass[-1] + AminoAcidMass[AminoAcid])
    
    Spectrum = [0]
    for Start in range(len(Peptide)):
        for End in range(Start + 1, len(Peptide) + 1):
            Mass = PrefixMass[End] - PrefixMass[Start]
            Spectrum.append(Mass)
    
    return sorted(Spectrum)

def LinearScore(Peptide, Spectrum):
    TheoreticalSpectrum = Counter(LinearSpectrum(Peptide))
    ExperimentalSpectrum = Counter(Spectrum)
    Score = 0

    for Mass in TheoreticalSpectrum:
        Score += min(TheoreticalSpectrum[Mass], ExperimentalSpectrum.get(Mass, 0))
    
    return Score


def Trim(Leaderboard, Spectrum, N):
    '''Compute the score of a linear peptide with respect to a spectrum.
    Keep the top N peptides in Leaderboard based on LinearScore (with ties).

    Input: A collection of peptides Leaderboard, a collection of integers Spectrum, and an integer N.
    Output: The N highest-scoring linear peptides on Leaderboard with respect to Spectrum.
    '''
    
    AminoAcidMass = {}
    with open("amino_acid_integer_mass_table.txt", "r") as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                key, value = parts[0], int(parts[1])
                AminoAcidMass[key] = value
    
    Alphabet = list(AminoAcidMass.keys())
    
    if not Leaderboard:
        return []
    LinearScores = []

    # For j ← 1 to |Leaderboard|:
    for Peptide in Leaderboard:
        Score = LinearScore(Peptide, Spectrum)
        LinearScores.append(Score)

    # Sort Leaderboard and LinearScores by descending score
    Combined = list(zip(Leaderboard, LinearScores))
    Combined.sort(key=lambda X: X[1], reverse=True)
    
    SortedLeaderboard = [Pair[0] for Pair in Combined]
    SortedScores = [Pair[1] for Pair in Combined]

    # For j ← N + 1 to |Leaderboard|:
    for Index in range(N, len(SortedLeaderboard)):
        if SortedScores[Index] < SortedScores[N - 1]:
            # Remove all peptides starting from the j-th peptide
            return SortedLeaderboard[:Index]

    return SortedLeaderboard



# Example
# -------------
Leaderboard = ["LAST", "ALST", "TLLT", "TQAS"]
Spectrum = [0, 71, 87, 101, 113, 158, 184, 188, 259, 271, 372]
N = 2
print(*Trim(Leaderboard, Spectrum, N))
# Output: LAST ALST



with open("dataset_30249_3.txt", "r") as file:
        lines = file.read().strip().split('\n')
        Leaderboard = lines[0].strip().split()
        Spectrum = list(map(int, lines[1].strip().split()))
        N = int(lines[2].strip())

final_leaderboard = Trim(Leaderboard, Spectrum, N)

with open("output.txt", "w") as out_file:
    out_file.write(' '.join(final_leaderboard) + '\n')
    
# Output: RGYDCSLVINCQHCGYMTMEGQSLGTPCPAQMQLTLGDCTEPLP 
# DHNWYMQCFYNESQLRIGTACTMGNHFSLGGWMTDPEFYDTPWY 
# VWLYLCVVAHTNYGDEHYLVTGCMMMWPHSNHCPRNWPLATRAI 
# AQYWPMKMRNQDLARLFSIFQFHGVAPDCYKDHHFEQGTMGLVC 
# EPFQSDRWTWFQPTIFNQNSHLISRQGFWDNTYGWFGWIYYFVD 
# QVDKVHDSEVHWNYDVLVQNMCRCLQMKEQCNFVDHLRYMHFPS




# -----------------------------------------------
# Leaderboard Cyclopeptide Sequencing
# -----------------------------------------------

from collections import Counter

def LoadAminoAcidMassTable(FileName="amino_acid_integer_mass_table.txt"):
    AminoAcidMass = {}
    with open(FileName, "r") as file:
        for line in file:
            Parts = line.strip().split()
            if len(Parts) == 2:
                Key, Value = Parts[0], int(Parts[1])
                AminoAcidMass[Key] = Value
    return AminoAcidMass

AminoAcidMass = LoadAminoAcidMassTable()
Alphabet = list(AminoAcidMass.keys())

def LinearSpectrum(Peptide):
    PrefixMass = [0]
    for AminoAcid in Peptide:
        PrefixMass.append(PrefixMass[-1] + AminoAcidMass[AminoAcid])
    Spectrum = [0]
    for Start in range(len(Peptide)):
        for End in range(Start + 1, len(Peptide) + 1):
            Spectrum.append(PrefixMass[End] - PrefixMass[Start])
    return sorted(Spectrum)

def CyclicSpectrum(Peptide):
    PrefixMass = [0]
    for AminoAcid in Peptide:
        PrefixMass.append(PrefixMass[-1] + AminoAcidMass[AminoAcid])
    PeptideMass = PrefixMass[-1]
    Spectrum = [0]
    for Start in range(len(Peptide)):
        for End in range(Start + 1, len(Peptide) + 1):
            Spectrum.append(PrefixMass[End] - PrefixMass[Start])
            if Start > 0 and End < len(Peptide):
                Spectrum.append(PeptideMass - (PrefixMass[End] - PrefixMass[Start]))
    return sorted(Spectrum)

def LinearScore(Peptide, Spectrum):
    Theoretical = Counter(LinearSpectrum(Peptide))
    Experimental = Counter(Spectrum)
    return sum(min(Theoretical[M], Experimental[M]) for M in Theoretical)

def Score(Peptide, Spectrum):
    Theoretical = Counter(CyclicSpectrum(Peptide))
    Experimental = Counter(Spectrum)
    return sum(min(Theoretical[M], Experimental[M]) for M in Theoretical)

def Mass(Peptide):
    return sum(AminoAcidMass[A] for A in Peptide)

def Expand(Leaderboard):
    Expanded = []
    for Peptide in Leaderboard:
        for AminoAcid in Alphabet:
            Expanded.append(Peptide + AminoAcid)
    return Expanded

def Trim(Leaderboard, Spectrum, N):
    if not Leaderboard:
        return []

    LinearScores = []
    for Peptide in Leaderboard:
        ScoreValue = LinearScore(Peptide, Spectrum)
        LinearScores.append(ScoreValue)

    Combined = list(zip(Leaderboard, LinearScores))
    Combined.sort(key=lambda X: X[1], reverse=True)

    SortedLeaderboard = [P[0] for P in Combined]
    SortedScores = [P[1] for P in Combined]

    for Index in range(N, len(SortedLeaderboard)):
        if SortedScores[Index] < SortedScores[N - 1]:
            return SortedLeaderboard[:Index]
    return SortedLeaderboard


def CyclicPermutations(Peptide):
    return [peptide[i:] + peptide[:i] for i in range(len(Peptide))]


def CanonicalCyclicForm(Peptide):
    mass_sequence = [AminoAcidMass[a] for a in peptide]
    rotations = [mass_sequence[i:] + mass_sequence[:i] for i in range(len(mass_sequence))]
    return tuple(min(rotations))

def LeaderboardCyclopeptideSequencing(N, Spectrum, return_all=False):
    '''
    Find cyclic peptide(s) having maximum score against an experimental spectrum.
    Cyclopeptide sequencing for spectra with errors.
    Retain N highest-scoring linear peptides including ties.

    Input: An integer N and a collection of integers Spectrum.
    Output: LeaderPeptide.
    '''
    
    Leaderboard = [""]
    LeaderPeptides = []
    MaxScore = 0
    ParentMass = max(Spectrum)

    while Leaderboard:
        Leaderboard = Expand(Leaderboard)
        NewLeaderboard = []
        for Peptide in Leaderboard:
            PepMass = Mass(Peptide)
            if PepMass == ParentMass:
                PeptideScore = Score(Peptide, Spectrum)
                if PeptideScore > MaxScore:
                    MaxScore = PeptideScore
                    LeaderPeptides = [Peptide]
                elif PeptideScore == MaxScore:
                    LeaderPeptides.append(Peptide)
            elif PepMass < ParentMass:
                NewLeaderboard.append(Peptide)
        Leaderboard = Trim(NewLeaderboard, Spectrum, N)

    if return_all:
        return LeaderPeptides

    else:
        return LeaderPeptides[0] if LeaderPeptides else ""



# Example
# -----------
N = 10
Spectrum = [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460]

Peptide = LeaderboardCyclopeptideSequencing(N, Spectrum)
print('-'.join(str(AminoAcidMass[A]) for A in Peptide))

# Output: 71-147-113-129



with open("dataset_30244_8.txt", "r") as file:
    N = int(file.readline().strip())
    Spectrum = list(map(int, file.readline().strip().split()))

Peptide = LeaderboardCyclopeptideSequencing(N, Spectrum)
print('-'.join(str(AminoAcidMass[a]) for a in Peptide))

# Output: 71-147-99-156-113-71-101-99-115-156-103-
# 163-115-147-114-71-128-115-97-128-114-103-113-113





# Run LeaderboardCyclopeptideSequencing on Spectrum25 with N = 1000. 
# You should find 38 linear peptides of maximum score 83 (corresponding 
# to 15 different cyclic peptides). What are they? 
# ----------------------------------------------------------------------
from collections import Counter

def LoadAminoAcidMassTable(FileName="amino_acid_integer_mass_table.txt"):
    AminoAcidMass = {}
    with open(FileName, "r") as file:
        for line in file:
            Parts = line.strip().split()
            if len(Parts) == 2:
                Key, Value = Parts[0], int(Parts[1])
                AminoAcidMass[Key] = Value
    return AminoAcidMass

AminoAcidMass = LoadAminoAcidMassTable()
 
def LinearSpectrum(Peptide):
    prefix_mass = [0]
    
    for aa in peptide:
        prefix_mass.append(prefix_mass[-1] + aa)
    spectrum = [0]

    for i in range(len(Peptide)):
        for j in range(i + 1, len(Peptide) + 1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])

    return sorted(spectrum)


def CyclicSpectrum(Peptide):
    spectrum = [0]
    n = len(Peptide)
    extended_peptide = peptide + peptide
    peptide_mass = sum(Peptide)

    for i in range(n):
        for j in range(i + 1, i + n):
            sub_peptide = extended_peptide[i:j]
            spectrum.append(sum(sub_peptide))

    spectrum.append(peptide_mass)

    return sorted(spectrum)


def Score(Peptide, Spectrum, cyclic=True):
    theo = CyclicSpectrum(Peptide) if cyclic else LinearSpectrum(Peptide)
    theo_count = Counter(theo)
    spec_count = Counter(Spectrum)
    return sum(min(theo_count[x], spec_count[x]) for x in theo_count)


def Expand(Peptides):
    return [pep + [mass] for pep in Peptides for mass in AminoAcidMasses]

 
def Trim(Leaderboard, Spectrum, N):
    scored = [(pep, Score(pep, Spectrum, cyclic=False)) for pep in Leaderboard]
    scored.sort(key=lambda x: x[1], reverse=True)
    
    if len(scored) <= N:
        return [pep for pep, _ in scored]
    
    cutoff = scored[N - 1][1]
    
    return [pep for pep, s in scored if s >= cutoff]


def LeaderboardCyclopeptideSequencing(Spectrum, N):
    leaderboard = [[]]
    leader_peptides = []
    leader_score = 0
    parent_mass = max(spectrum)

    while leaderboard:
        leaderboard = Expand(leaderboard)
        new_leaderboard = []

        for peptide in leaderboard:
            mass = sum(Peptide)
            
            if mass == parent_mass:
                pep_score = Score(Peptide, Spectrum, cyclic=True)
                
                if pep_score > leader_score:
                    leader_peptides = [peptide]
                    leader_score = pep_score

                elif pep_score == leader_score:
                    leader_peptides.append(peptide)
                
                new_leaderboard.append(peptide)

            elif mass < parent_mass:
                new_leaderboard.append(peptide)

        leaderboard = Trim(new_leaderboard, Spectrum, N)

    return ['-'.join(map(str, pep)) for pep in leader_peptides]

 

# Example
# -------------
with open("dataset_30245_2.txt", "r") as file:
    N = int(file.readline())
    Spectrum = list(map(int, file.readline().strip().split()))

Peptides = LeaderboardCyclopeptideSequencing(Spectrum, N)
print(len(Peptides)) # Output: 34

with open("output.txt", "w") as file:
    output_lines = []
    for peptide_str in Peptides:
        # Split peptide string into masses (strings)
        masses = peptide_str.split('-')  
        
        # Convert each to int and then to str of the mass itself (or get from AminoAcidMass if needed)
        # The peptides already store masses, so just join them:
        output_lines.append('-'.join(masses))
    
    file.write(' '.join(output_lines) + "\n")



with open("Tyrocidine_B1_Spectrum_25.txt", "r") as file:
    Spectrum = list(map(int, file.readline().strip().split()))
N = 1000

Peptides = LeaderboardCyclopeptideSequencing(Spectrum, N)
print(len(Peptides)) # Output: 38

with open("output.txt", "w") as file:
    file.write(' '.join('-'.join(str(AminoAcidMass[a]) for a in peptide) 
                        for peptide in Peptides) + "\n")


# The same thing
# -------------------
with open("dataset_30244_10.txt", "r") as file:
    N = int(file.readline())
    Spectrum = list(map(int, file.readline().strip().split()))

Peptides = LeaderboardCyclopeptideSequencing(N, Spectrum, return_all=True)
print(len(Peptides)) # Output: 38

with open("output.txt", "w") as file:
    file.write(' '.join('-'.join(str(AminoAcidMass[a]) for a in peptide) 
                        for peptide in Peptides) + "\n")


# -----------------------------------------------
# Leaderboard Cyclopeptide Sequencing on the extended amino acid alphabet 
# -----------------------------------------------

from collections import Counter

AminoAcidMasses = list(range(57, 201))
 
def LinearSpectrum(Peptide):
    prefix_mass = [0]
    
    for aa in peptide:
        prefix_mass.append(prefix_mass[-1] + aa)
    spectrum = [0]

    for i in range(len(Peptide)):
        for j in range(i + 1, len(Peptide) + 1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])

    return sorted(spectrum)


def CyclicSpectrum(Peptide):
    spectrum = [0]
    n = len(Peptide)
    extended_peptide = peptide + peptide
    peptide_mass = sum(Peptide)

    for i in range(n):
        for j in range(i + 1, i + n):
            sub_peptide = extended_peptide[i:j]
            spectrum.append(sum(sub_peptide))

    spectrum.append(peptide_mass)

    return sorted(spectrum)


def Score(Peptide, spectrum, cyclic=True):
    theo = CyclicSpectrum(Peptide) if cyclic else LinearSpectrum(Peptide)
    theo_count = Counter(theo)
    spec_count = Counter(spectrum)
    return sum(min(theo_count[x], spec_count[x]) for x in theo_count)


def Expand(peptides):
    return [pep + [mass] for pep in peptides for mass in AminoAcidMasses]

 
def Trim(leaderboard, spectrum, n):
    scored = [(pep, Score(pep, spectrum, cyclic=False)) for pep in leaderboard]
    scored.sort(key=lambda x: x[1], reverse=True)
    
    if len(scored) <= n:
        return [pep for pep, _ in scored]
    
    cutoff = scored[n - 1][1]
    
    return [pep for pep, s in scored if s >= cutoff]


def LeaderboardCyclopeptideSequencing(Spectrum, N):
    leaderboard = [[]]
    leader_peptides = []
    leader_score = 0
    parent_mass = max(Spectrum)

    while leaderboard:
        leaderboard = Expand(leaderboard)
        new_leaderboard = []

        for peptide in leaderboard:
            mass = sum(Peptide)
            
            if mass == parent_mass:
                pep_score = Score(Peptide, Spectrum, cyclic=True)
                
                if pep_score > leader_score:
                    leader_peptides = [peptide]
                    leader_score = pep_score

                elif pep_score == leader_score:
                    leader_peptides.append(Peptide)
                
                new_leaderboard.append(Peptide)

            elif mass < parent_mass:
                new_leaderboard.append(Peptide)

        leaderboard = Trim(new_leaderboard, Spectrum, n)

    return ['-'.join(map(str, pep)) for pep in leader_peptides]

 

# Example
# -------------
with open("dataset_30245_2.txt", "r") as file:
    N = int(file.readline())
    Spectrum = list(map(int, file.readline().strip().split()))

Peptides = LeaderboardCyclopeptideSequencing(Spectrum, N)
print(len(Peptides)) # Output: 34

with open("output.txt", "w") as file:
    output_lines = []
    for peptide_str in Peptides:
        # Split peptide string into masses (strings)
        masses = peptide_str.split('-')  
        
        # Convert each to int and then to str of the mass itself (or get from AminoAcidMass if needed)
        # The peptides already store masses, so just join them:
        output_lines.append('-'.join(masses))
    
    file.write(' '.join(output_lines) + "\n")

# Output: 113-147-97-186-147-114-128-98-65-99-128 
# 147-97-186-147-114-128-98-65-99-128-113 
# 128-113-147-97-186-147-114-128-98-65-99 
# 114-147-186-97-147-113-128-99-65-98-128  ...



# -----------------------------------------------
# Spectral Convolution 
# -----------------------------------------------

def SpectralConvolution(Spectrum):
    '''Compute the convolution of a spectrum.

    Input: A collection of integers Spectrum in increasing order.
    Output: The list of elements in the convolution of Spectrum. 
    If an element has multiplicity k, it should appear exactly k times.
    '''
    
    convolution = []

    for i in range(len(Spectrum)):
        for j in range(len(Spectrum)):
            diff = Spectrum[j] - Spectrum[i]
            if diff > 0:
                convolution.append(diff)

    return convolution



# Examples
# ------------
Spectrum = [0, 137, 186, 323]
result = SpectralConvolution(Spectrum)
print(' '.join(map(str, result)))
# Output: 137 186 323 49 186 137



with open("dataset_30246_4.txt", "r") as file:
    line = file.readline().strip()
    Spectrum = list(map(int, line.split()))
    
convolution = SpectralConvolution(Spectrum)

with open("output.txt", "w") as out_file:
    out_file.write(' '.join(map(str, convolution)))

# Output: 71 71 97 99 101 128 128 137 137 142 196 199 
# 208 225 229 229 236 270 274 279 300 324 326 333 345 ...




# What element has largest multiplicity in the spectral 
# convolution of Spectrum = {0 57 118 179 236 240 301?
​# -------------------------------------------------------

Spectrum = [0, 57, 118, 179, 236, 240, 301]
result = SpectralConvolution(Spectrum)

def CountMultiplicity(lst):
    counts = {}
    for item in lst:
        counts[item] = counts.get(item, 0) + 1
    return counts

CountMultiplicity(result)

# Output:
# {57: 2, 118: 2, 179: 2, 236: 1, 240: 1, 301: 1, 
# 61: 4, 122: 3, 183: 2, 244: 1, 4: 1, 65: 1}

# Answer: 61



# -----------------------------------------------
# Convolution Cyclopeptide Sequencing
# -----------------------------------------------

from collections import Counter
from itertools import combinations

 
def SpectralConvolution(Spectrum):
    convolution = []

    for i in range(len(Spectrum)):

        for j in range(i):
            diff = Spectrum[i] - Spectrum[j]

            if 57 <= diff <= 200:
                convolution.append(diff)

    return convolution


def LinearSpectrum(Peptide):
    prefix_mass = [0]

    for mass in peptide:
        prefix_mass.append(prefix_mass[-1] + mass)

    linear_spec = [0]

    for i in range(len(prefix_mass)):
        for j in range(i + 1, len(prefix_mass)):
            linear_spec.append(prefix_mass[j] - prefix_mass[i])

    return sorted(linear_spec)


def CyclicSpectrum(Peptide):
    prefix_mass = [0]

    for mass in peptide:
        prefix_mass.append(prefix_mass[-1] + mass)

    peptide_mass = prefix_mass[-1]
    cyclic_spec = [0]

    for i in range(len(Peptide)):
        
        for j in range(i + 1, len(Peptide) + 1):
            cyclic_spec.append(prefix_mass[j] - prefix_mass[i])

            if i > 0 and j < len(Peptide):
                cyclic_spec.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))

    return sorted(cyclic_spec)


def Score(Peptide, Spectrum, cyclic=True):
    if cyclic:
        pep_spec = CyclicSpectrum(Peptide)

    else:
        pep_spec = LinearSpectrum(Peptide)

    spectrum_count = Counter(Spectrum)
    pep_spec_count = Counter(pep_spec)

    return sum(min(pep_spec_count[m], Spectrum_count[m]) for m in pep_spec_count)


def Trim(leaderboard, Spectrum, N):
    scored = [(pep, Score(pep, Spectrum, cyclic=False)) for pep in leaderboard]
    scored.sort(key=lambda x: x[1], reverse=True)
    
    if len(scored) <= N:
        return [x[0] for x in scored]

    threshold = scored[N-1][1]
    
    return [pep for pep, s in scored if s >= threshold]


def LeaderboardCyclopeptideSequencing(Spectrum, N, AminoAcidMasses):
    leaderboard = [[]]
    leader_peptide = []
    parent_mass = max(spectrum)

    while leaderboard:
        new_leaderboard = []

        for peptide in leaderboard:

            for mass in AminoAcidMasses:
                new_peptide = peptide + [mass]

                if sum(new_peptide) == parent_mass:

                    if Score(new_peptide, Spectrum) > Score(leader_peptide, Spectrum):
                        leader_peptide = new_peptide

                    new_leaderboard.append(new_peptide)

                elif sum(new_peptide) < parent_mass:
                    new_leaderboard.append(new_peptide)

        leaderboard = Trim(new_leaderboard, Spectrum, N)

    return leader_peptide


def ConvolutionCyclopeptideSequencing(Spectrum, M, N):
    '''Compute the convolution of an experimental spectrum. Then select the M most 
    frequent elements between 57 and 200 in the convolution to form an extended 
    alphabet of candidate amino acid masses. In order to be fair, we should include 
    the top M elements of the convolution "with ties". Finally, run the algorithm 
    LeaderboardCyclopeptideSequencing, where the amino acid masses are restricted 
    to this alphabet.
    
    Input: An integer M, an integer N, and a collection of (possibly repeated) 
    integers Spectrum.
    Output: A cyclic peptide LeaderPeptide with amino acids taken only from 
    the top M elements (and ties) of the convolution of Spectrum that fall 
    between 57 and 200, and where the size of Leaderboard is restricted to 
    the top N (and ties).
    '''
    
    Spectrum.sort()
    conv = SpectralConvolution(Spectrum)
    freq = Counter(conv)
    filtered = [(mass, count) for mass, count in freq.items() if 57 <= mass <= 200]
    sorted_filtered = sorted(filtered, key=lambda x: x[1], reverse=True)

    result_masses = []

    if sorted_filtered:
        threshold = sorted_filtered[min(M, len(sorted_filtered))-1][1]

        for mass, count in sorted_filtered:

            if count >= threshold:
                result_masses.append(mass)


    return LeaderboardCyclopeptideSequencing(Spectrum, N, result_masses)

 

# Examples
# -----------
M = 20
N = 60
Spectrum = [57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285,
            299, 307, 323, 356, 364, 394, 422, 493]
peptide = ConvolutionCyclopeptideSequencing(Spectrum, M, N)

print('-'.join(map(str, peptide)))
# Output: 99-71-137-57-57-72


with open("dataset_30246_7.txt", "r") as file:
    M = int(file.readline().strip())
    N = int(file.readline().strip())
    Spectrum = list(map(int, file.readline().strip().split()))
    
peptide = ConvolutionCyclopeptideSequencing(Spectrum, M, N)

with open("output.txt", "w") as file:
    file.write('-'.join(map(str, peptide)) + '\n')
# Output: 147-113-87-128-156-115-97-113-87-103-156-147-163-129



# Run ConvolutionCyclopeptideSequencing on Spectrum25 with N = 1000 
# and M = 20. Identify the 86 highest-scoring linear peptides.
# ------------------------------------------------------------------

from collections import Counter

def ConvultionOfSpectra(Spectrum):
    """(list,list) -> list"""
    Convultion = []
    n = len(Spectrum)
    for i in range(1, n):
        for j in range(i):
            diff = Spectrum[i] - Spectrum[j]
            if diff > 0:  # Exclude negative values and zeros
                Convultion.append(diff)
    return sorted(Convultion)

def LinearSpectrum(Peptide):
    PrefixMass = [0]
    for Mass in Peptide:
        PrefixMass.append(PrefixMass[-1] + Mass)
    LinearSpec = [0]
    for i in range(len(PrefixMass)):
        for j in range(i + 1, len(PrefixMass)):
            LinearSpec.append(PrefixMass[j] - PrefixMass[i])
    return sorted(LinearSpec)

def CyclicSpectrum(Peptide):
    PrefixMass = [0]
    for Mass in Peptide:
        PrefixMass.append(PrefixMass[-1] + Mass)
    PeptideMass = PrefixMass[-1]
    CyclicSpec = [0]
    for i in range(len(Peptide)):
        for j in range(i + 1, len(Peptide) + 1):
            CyclicSpec.append(PrefixMass[j] - PrefixMass[i])
            if i > 0 and j < len(Peptide):
                CyclicSpec.append(PeptideMass - (PrefixMass[j] - PrefixMass[i]))
    return sorted(CyclicSpec)

def LinearpeptideScoring(Peptide, Spectrum):
    """(list, list) -> int"""
    Theo = Counter(LinearSpectrum(Peptide))
    Exp = Counter(Spectrum)
    Score = sum(min(Theo[m], Exp[m]) for m in Theo)
    return Score

def CyclopeptideScoring(Peptide, Spectrum):
    Theo = Counter(CyclicSpectrum(Peptide))
    Exp = Counter(Spectrum)
    Score = sum(min(Theo[m], Exp[m]) for m in Theo)
    return Score

def Trim(LeaderBord, Spectrum, N):
    """(list of int, list of int, int) -> list of str"""
    k = len(LeaderBord)
    LinearScores = [0] * k
    for j in range(k):
        Peptide = LeaderBord[j]
        LinearScores[j] = LinearpeptideScoring(Peptide, Spectrum)
    sorted_LeaderBord = [x for _, x in sorted(list(zip(LinearScores, LeaderBord)), reverse=True)]
    sorted_LinearScores = sorted(LinearScores, reverse=True)

    if N >= k:
        return sorted_LeaderBord
    
    for j in range(N, k):
        if sorted_LinearScores[j] < sorted_LinearScores[N - 1]:
            return sorted_LeaderBord[:j]
    return sorted_LeaderBord


def ConvolutionCyclopeptideSequencingMax(Spectrum, M, N):
    conv = ConvultionOfSpectra(Spectrum)
    freq = Counter(conv)
    filtered = [(m, c) for m, c in freq.items() if 57 <= m <= 200]
    filtered.sort(key=lambda x: x[1], reverse=True)
    
    if len(filtered) == 0:
        AminoAcidMasses = []
    else:
        cutoff = filtered[M - 1][1] if M <= len(filtered) else filtered[-1][1]
        AminoAcidMasses = [m for m, c in filtered if c >= cutoff]

    ParentMass = max(Spectrum)
    Leaderboard = [[]]
    BestScore = 0
    BestPeptides = set()

    while Leaderboard:
        NextBoard = []
        for pep in Leaderboard:
            for mass in AminoAcidMasses:
                NewPep = pep + [mass]
                s = sum(NewPep)
                if s == ParentMass:
                    sc = CyclopeptideScoring(NewPep, Spectrum)
                    if sc > BestScore:
                        BestScore = sc
                        BestPeptides = {tuple(NewPep)}
                    elif sc == BestScore:
                        BestPeptides.add(tuple(NewPep))
                    NextBoard.append(NewPep)
                elif s < ParentMass:
                    NextBoard.append(NewPep)
        Leaderboard = Trim(NextBoard, Spectrum, N)

    return BestScore, [list(p) for p in BestPeptides]


# Example
# -------------

with open("dataset_30246_8.txt", 'r') as file:
    M = int(file.readline().strip())
    N = int(file.readline().strip())
    Spectrum_line = file.readline().strip()
    Spectrum = list(map(int, Spectrum_line.split()))

BestScore, BestPeptides = result

with open("output.txt", "w") as out_file:
    for pep in BestPeptides:
        out_file.write("-".join(map(str, pep)) + "\n")

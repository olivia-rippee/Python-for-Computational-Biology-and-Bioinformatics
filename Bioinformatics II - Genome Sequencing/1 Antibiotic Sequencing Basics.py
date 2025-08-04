import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics II - Genome Sequencing/Data")


# -----------------------------------------------
# Protein Translation
# -----------------------------------------------

with open("RNA_codon_table_1.txt", "r") as file:
    GeneticCode = {parts[0]: parts[1] if len(parts) > 1 else '' 
                   for parts in (line.strip().split() for line in file)}


def Translate(RNA, GeneticCode):
    '''Translate an RNA string into an amino acid string.
    Input: An RNA string Pattern and the array (dictionary) GeneticCode.
    Output: The translation of Pattern into an amino acid string Peptide.
    '''
    protein = ''
    for i in range(0, len(RNA), 3):
        codon = RNA[i:i+3]
        amino_acid = GeneticCode.get(codon, '')
        if amino_acid == '':
            break  # Stop translation at stop codon
        protein += amino_acid
    return protein



# Example
# ---------
RNA = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
print(Translate(RNA, GeneticCode))
# Output: MAMAPRTEINSTRING


with open("dataset_30213_4.txt") as file:
    RNA = file.readline()
    
protein = Translate(RNA, GeneticCode)

with open("output.txt", "w") as out_file:
    out_file.write(protein)
    
# Output: MRNFCVLAGSTRGCNRSEGSSVACPQKKVPAEGYRCVRCRALRIASMPYEA...


RNA1 = "CCAAGUACAGAGAUUAAC"
print(Translate(RNA1, GeneticCode)) # PSTEIN

RNA2 = "CCAAGAACAGAUAUCAAU"
print(Translate(RNA2, GeneticCode)) # PRTDIN

RNA3 = "CCACGUACUGAAAUUAAC"
print(Translate(RNA3, GeneticCode)) # PRTEIN

RNA4 = "CCGAGGACCGAAAUCAAC"
print(Translate(RNA4, GeneticCode)) # PRTEIN




# -----------------------------------------------
# Number of DNA strings (with length k) transcribe and translate into a given sequence
# -----------------------------------------------

# Protein key to codon(s)
GeneticCode = {
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],               # Alanine
    'C': ['UGU', 'UGC'],                             # Cysteine
    'D': ['GAU', 'GAC'],                             # Aspartic Acid
    'E': ['GAA', 'GAG'],                             # Glutamic Acid
    'F': ['UUU', 'UUC'],                             # Phenylalanine
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],               # Glycine
    'H': ['CAU', 'CAC'],                             # Histidine
    'I': ['AUU', 'AUC', 'AUA'],                      # Isoleucine
    'K': ['AAA', 'AAG'],                             # Lysine
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'], # Leucine
    'M': ['AUG'],                                    # Methionine (Start)
    'N': ['AAU', 'AAC'],                             # Asparagine
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],               # Proline
    'Q': ['CAA', 'CAG'],                             # Glutamine
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], # Arginine
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], # Serine
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],               # Threonine
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],               # Valine
    'W': ['UGG'],                                    # Tryptophan
    'Y': ['UAU', 'UAC']                              # Tyrosine
}

RNAtoDNA = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}


from itertools import product

def CountDNASequences(AminoAcidSequence, k):
    if k % 3 != 0:
        raise ValueError("k must be a multiple of 3.")
    if len(AminoAcidSequence) != k // 3:
        raise ValueError("Length of amino acid sequence must match k / 3.")
    
    Total = 1
    for AminoAcid in AminoAcidSequence:
        Codons = GeneticCode.get(AminoAcid)
        if not Codons:
            raise ValueError(f"Invalid or unsupported amino acid: {AminoAcid}")
        Total *= len(Codons)
    return Total


def RNAtoDNASequence(RNACodon):
    return ''.join(RNAtoDNA[Base] for Base in RNACodon)


def GenerateDNASequences(AminoAcidSequence, k, MaxSequences=5):
    if k% 3 != 0 or len(AminoAcidSequence) != k // 3:
        raise ValueError("Mismatch between K and amino acid sequence length.")
    
    CodonLists = [GeneticCode[AminoAcid] for AminoAcid in AminoAcidSequence]
    Sequences = []
    
    for CodonCombo in product(*CodonLists):
        DNASeq = ''.join(RNAtoDNASequence(Codon) for Codon in CodonCombo)
        Sequences.append(DNASeq)
        if len(Sequences) >= MaxSequences:
            break
    return Sequences




# Example
# ------------
TyrocidineB1 = "VKLFPWFNQY"  # Val-Lys-Leu-Phe-Pro-Trp-Phe-Asn-Gln-Tyr
print(CountDNASequences(TyrocidineB1, k=30))
print(*GenerateDNASequences(TyrocidineB1, k=30))

# Output: 6144      CAAUUUAAUAAAGGAACCAAAUUAGUUAUA 
# CAAUUUAAUAAAGGAACCAAAUUAGUUAUG CAAUUUAAUAAAGGAACCAAAUUAGUCAUA 
# CAAUUUAAUAAAGGAACCAAAUUAGUCAUG CAAUUUAAUAAAGGAACCAAAUUGGUUAUA

# Product of different ways to get each AA: 4×2×6×2×4×1×2×2×2×2=6144


# How many DNA strings transcribe and translate into the amino acid string CYCLIC?
AminoAcidSequence="CYCLIC"
k = 3*len(AminoAcidSequence)
print(CountDNASequences(AminoAcidSequence, k))



# -----------------------------------------------
# Section of Genome Encoding a Peptide
# -----------------------------------------------

# Codon key to amino acid
GeneticCode = {
    "AUG": "M", "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y",
    "UAC": "Y", "UGU": "C", "UGC": "C", "UGG": "W", "CUU": "L",
    "CUC": "L", "CUA": "L", "CUG": "L", "CCU": "P", "CCC": "P",
    "CCA": "P", "CCG": "P", "CAU": "H", "CAC": "H", "CAA": "Q",
    "CAG": "Q", "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AUU": "I", "AUC": "I", "AUA": "I", "ACU": "T", "ACC": "T",
    "ACA": "T", "ACG": "T", "AAU": "N", "AAC": "N", "AAA": "K",
    "AAG": "K", "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V", "GCU": "A",
    "GCC": "A", "GCA": "A", "GCG": "A", "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E", "GGU": "G", "GGC": "G", "GGA": "G",
    "GGG": "G"
}

def TranscribeDNAToRNA(DNA):
    return DNA.replace('T', 'U')


def ReverseComplement(DNA):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(DNA))


def TranslateRNAToPeptide(RNA, GeneticCode):
    peptide = ''
    for i in range(0, len(RNA), 3):
        codon = RNA[i:i+3]
        if len(codon) < 3 or codon not in GeneticCode:
            return None
        peptide += GeneticCode[codon]
    return peptide


def PeptideEncoding(DNA, Peptide, GeneticCode):
    '''
    Find substrings of a genome encoding a given amino acid sequence.
    Input: A DNA string Text, an amino acid string Peptide, and the array GeneticCode.
    Output: All substrings of Text encoding Peptide (if any such substrings exist).
    '''
    DNA = DNA.upper()
    peptide_length = len(Peptide) * 3
    substrings = []

    for i in range(len(DNA) - peptide_length + 1):
        fragment = DNA[i:i+peptide_length]

        # Transcribe and translate original strand
        rna = TranscribeDNAToRNA(fragment)
        if TranslateRNAToPeptide(rna, GeneticCode) == Peptide:
            substrings.append(fragment)

        # Transcribe and translate reverse complement strand
        rev_comp = ReverseComplement(fragment)
        rev_rna = TranscribeDNAToRNA(rev_comp)
        if TranslateRNAToPeptide(rev_rna, GeneticCode) == Peptide:
            substrings.append(fragment)

    return substrings


# Example
# -------------
DNA = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
Peptide = "MA"
print(*PeptideEncoding(DNA, Peptide, GeneticCode))
# Output: ATGGCC GGCCAT ATGGCC



with open("dataset_30213_7.txt") as file:
    DNA = file.readline().strip()
    Peptide = file.readline().strip()
    
results = PeptideEncoding(DNA, Peptide, GeneticCode)

with open("output.txt", "w") as out_file:
    out_file.write(' '.join(results))   
    
# Output: ATGTGTAATGCTACTTATCAACCATGCATG ATGTGTAATGCCACCTACCAGCCGTGCATG 
# CATGCAAGGTTGATAAGTGGCGTTACACAT ATGTGCAACGCAACTTATCAGCCTTGCATG 
# CATGCACGGTTGGTACGTAGCATTGCACAT CATACATGGTTGGTAAGTAGCGTTGCACAT 
# ATGTGCAACGCGACGTATCAACCGTGTATG CATGCAGGGCTGATACGTGGCATTACACAT ...



# How many starting positions in Bacillus brevis encode Tyrocidine B1?
# --------------------------------------------------------------------
TyrocidineB1 = "VKLFPWFNQY"  # Val-Lys-Leu-Phe-Pro-Trp-Phe-Asn-Gln-Tyr

with open("Bacillus_brevis.txt", "r") as file:
    lines = file.readlines()
    Genome = ''.join(line.strip() for line in lines)
    
results = PeptideEncoding(Genome, TyrocidineB1, GeneticCode)
print(len(results))
# Output: 0
# Nonribosomal peptides (NRPs) do not adhere to the Central Dogma, so we cannot infer them from the genome




# -----------------------------------------------
# Number of Subpeptides for Sequencing - Cyclic
# -----------------------------------------------

def CountSubpeptidesCyclic(n):
    '''Calculate how many subpeptides a cyclic peptide of length n has.
    Input: An integer n.
    Output: The number of subpeptides of a cyclic peptide of length n.
    '''
    return n*(n-1)


# Example
# ----------
with open("dataset_30215_3.txt", "r") as file:
    n = int(file.readline())
print(CountSubpeptidesCyclic(n)) 
# Output: 670732302





# -----------------------------------------------
# Linear Spectrum of a Peptide
# -----------------------------------------------

def LinearSpectrum(Peptide):
    '''Generate the linear spectrum of peptide fragments in increasing order.
    Input: An amino acid string Peptide.
    Output: The linear spectrum of Peptide.
    '''
    AminoAcidMass = {}
    with open("amino_acid_integer_mass_table.txt", "r") as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                key, value = parts[0], int(parts[1])
                AminoAcidMass[key] = value
    
    Alphabet = list(AminoAcidMass.keys())
    
    PrefixMass = [0]
    for i in range(len(Peptide)):
        for symbol in Alphabet:
            if symbol == Peptide[i]:
                PrefixMass.append(PrefixMass[i] + AminoAcidMass[symbol])
                break
    
    LinearSpectrum = [0]
    for i in range(len(Peptide)):
        for j in range(i + 1, len(Peptide) + 1):
            LinearSpectrum.append(PrefixMass[j] - PrefixMass[i])
    
    return sorted(LinearSpectrum)



# Examples
# -----------
Peptide = "NQEL"
print(*LinearSpectrum(Peptide))
# Output: 0 113 114 128 129 242 242 257 370 371 484



with open("dataset_30248_2.txt", "r") as file:
    Peptide = file.readline().strip()

spectrum = LinearSpectrum(Peptide)

with open("output.txt", "w") as out_file:
    out_file.write(' '.join(str(fragment) for fragment in spectrum))





# -----------------------------------------------
# Generating Theoretical Spectrum of Cyclic Peptide
# -----------------------------------------------

def GenerateTheoreticalSpectrum(Peptide):
    '''Generate the spectrum of a cyclic peptide fragments in increasing order.
    Input: An amino acid string Peptide.
    Output: The cyclic spectrum of Peptide.
    '''
    AminoAcidMass = {}
    with open("amino_acid_integer_mass_table.txt", "r") as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                key, value = parts[0], int(parts[1])
                AminoAcidMass[key] = value

    Alphabet = list(AminoAcidMass.keys())

    PrefixMass = [0]
    for i in range(len(Peptide)):
        for symbol in Alphabet:
            if symbol == Peptide[i]:
                PrefixMass.append(PrefixMass[i] + AminoAcidMass[symbol])
                break

    peptideMass = PrefixMass[len(Peptide)]
    
    CyclicSpectrum = [0]
    for i in range(len(Peptide)):
        for j in range(i + 1, len(Peptide) + 1):
            mass = PrefixMass[j] - PrefixMass[i]
            CyclicSpectrum.append(mass)
            # Add cyclic subpeptide mass
            if i > 0 and j < len(Peptide):
                CyclicSpectrum.append(peptideMass - mass)

    return sorted(CyclicSpectrum)

    
# Examples
# ----------
Peptide = "LEQN"  
cyclospectrum = GenerateTheoreticalSpectrum(Peptide)
print(*cyclospectrum)
# Output: 0 113 114 128 129 227 242 242 257 355 356 370 371 484



with open("dataset_30215_4.txt", "r") as file:
    Peptide = file.readline().strip()

cyclospectrum = GenerateTheoreticalSpectrum(Peptide)

with open("output.txt", "w") as out_file:
    out_file.write(' '.join(str(fragment) for fragment in cyclospectrum))

# Output: 0 71 71 99 101 103 103 113 113 113 115 115 128 128 129 129 142 ...



# -----------------------------------------------
# Counting Peptides with Given Mass Problem: 
# -----------------------------------------------

def CountPeptidesWithMass(AminoAcidMasses, TargetMass):
    '''Compute the number of peptides of given mass.
    Input: An integer m.
    Output: The number of linear peptides having integer mass m.
    '''
    dp = [0] * (TargetMass + 1)
    dp[0] = 1  # base case: one way to make mass 0 (empty peptide)
    for mass in range(1, TargetMass + 1):
        for a in AminoAcidMasses:
            if mass - a >= 0:
                dp[mass] += dp[mass - a]
    return dp[TargetMass]


# Example
# ---------------
AminoAcidMasses = [57, 71, 87, 97, 99, 101, 103, 113, 114,
    115, 128, 129, 131, 137, 147, 156, 163, 186]

print(CountPeptidesWithMass(AminoAcidMasses, TargetMass=1253))
# Output: 7972200649745

# BruteForceCyclopeptideSequencing is completely impractical
# Number of peptides ~= k · C**m




# -----------------------------------------------
# Number of Subpeptides for Sequencing - Linear
# -----------------------------------------------

def CountSubpeptidesLinear(n):
    '''Calculate how many subpeptides a linear peptide of length n has.
    Input: An integer n.
    Output: The number of subpeptides of a cyclic peptide of length n.
    '''
    return (n * (n + 1))//2 + 1


# Examples
# -------------
print(CountSubpeptidesLinear(n=4))
# Output: 11

print(CountSubpeptidesLinear(n=16801))
# Output: 141145202



# -----------------------------------------------
# Branch-and-Bound Algorithm for Cyclopeptide Sequencing
# -----------------------------------------------

from collections import Counter

def Cyclospectrum(Peptide):
    '''Generate the cyclospectrum of a cyclic peptide.
    Input: A list of integer masses representing the peptide.
    Output: The sorted cyclospectrum (list of masses) of the peptide.
    '''
    N = len(Peptide)
    Spectrum = [0, sum(Peptide)]

    for Length in range(1, N):
        for Start in range(N):
            Subpeptide = Peptide[Start:Start + Length]
            if Start + Length > N:
                Subpeptide += Peptide[:(Start + Length) % N]
            Spectrum.append(sum(Subpeptide))

    return sorted(Spectrum)


def LinearSpectrum(Peptide):
    '''Generate the linear spectrum of a peptide (used for consistency checking).
    Input: A list of integer masses representing the peptide.
    Output: The sorted linear spectrum of the peptide.
    '''
    PrefixMass = [0]
    for MassValue in Peptide:
        PrefixMass.append(PrefixMass[-1] + MassValue)
    
    Spectrum = [0]
    for I in range(len(Peptide)):
        for J in range(I + 1, len(Peptide) + 1):
            Spectrum.append(PrefixMass[J] - PrefixMass[I])
    
    return sorted(Spectrum)


def IsConsistent(Peptide, Spectrum):
    '''
    Check if the peptide's linear spectrum is consistent with the given spectrum.
    Input: A candidate peptide and the experimental spectrum.
    Output: Boolean (True if consistent, False otherwise).
    '''
    PeptideSpectrum = Counter(LinearSpectrum(Peptide))
    SpectrumCounter = Counter(Spectrum)
    
    for Mass in PeptideSpectrum:
        if PeptideSpectrum[Mass] > SpectrumCounter[Mass]:
            return False
    return True


def CyclopeptideSequencing(Spectrum):
    '''Perform cyclopeptide sequencing using a branch-and-bound algorithm.
    Input: The experimental mass spectrum.
    Output: All possible peptides that match the given spectrum.
    '''
    CandidatePeptides = [[]]
    FinalPeptides = []
    ParentMass = max(Spectrum)
    AminoAcidMasses = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115,
        128, 129, 131, 137, 147, 156, 163, 186]

    while CandidatePeptides:
        # Expand logic inlined here
        ExpandedPeptides = []
        for Peptide in CandidatePeptides:
            for Mass in AminoAcidMasses:
                ExpandedPeptides.append(Peptide + [Mass])
        CandidatePeptides = []

        for Peptide in ExpandedPeptides:
            PeptideMass = sum(Peptide)  # Mass inlined
            if PeptideMass == ParentMass:
                if Cyclospectrum(Peptide) == sorted(Spectrum) and Peptide not in FinalPeptides:
                    FinalPeptides.append(Peptide)
            elif IsConsistent(Peptide, Spectrum):
                CandidatePeptides.append(Peptide)

    return FinalPeptides



# Examples
# -------------
Spectrum = [0, 113, 128, 186, 241, 299, 314, 427]
result = CyclopeptideSequencing(Spectrum)
print(' '.join(['-'.join(map(str, peptide)) for peptide in result]))
# Output: 113-128-186 113-186-128 128-113-186 128-186-113 186-113-128 186-128-113



with open("dataset_30217_6.txt", "r") as file:
    line = file.readline().strip()
    Spectrum = list(map(int, line.split()))

result = CyclopeptideSequencing(Spectrum)
formatted = ['-'.join(map(str, peptide)) for peptide in result]

with open("output.txt", "w") as out_file:
    out_file.write(' '.join(formatted))
    
# Output: 71-103-131-163-103-113-131-113-147-101-114 
# 71-114-101-147-113-131-113-103-163-131-103 
# 101-114-71-103-131-163-103-113-131-113-147 ...


    
    
# Which of the following cyclic peptides could have generated the theoretical spectrum 
# 0 71 101 113 131 184 202 214 232 285 303 315 345 416?
# -------------------------------------------------------------------------------------  

MassTable = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
    'T': 101, 'C': 103, 'I': 113, 'L': 113,
    'N': 114, 'D': 115, 'K': 128, 'Q': 128,
    'E': 129, 'M': 131, 'H': 137, 'F': 147,
    'R': 156, 'Y': 163, 'W': 186
}

def PeptideMass(Peptide):
    return sum(MassTable[AA] for AA in Peptide)


def CyclicSpectrum(Peptide):
    PrefixMass = [0]
    for AA in Peptide:
        PrefixMass.append(PrefixMass[-1] + MassTable[AA])
    PeptideMassValue = PrefixMass[-1]
    Spectrum = [0]
    
    Length = len(Peptide)
    for I in range(Length):
        for J in range(I + 1, Length + 1):
            Spectrum.append(PrefixMass[J] - PrefixMass[I])
            if I > 0 and J < Length:
                Spectrum.append(PeptideMassValue - (PrefixMass[J] - PrefixMass[I]))
    
    return sorted(Spectrum)


def MatchSpectrum(Peptide, TargetSpectrum):
    ComputedSpectrum = CyclicSpectrum(Peptide)
    return sorted(ComputedSpectrum) == sorted(TargetSpectrum)


def FindMatchingPeptides(Peptides, TargetSpectrum):
    Matches = []
    for Peptide in Peptides:
        if MatchSpectrum(Peptide, TargetSpectrum):
            Matches.append(Peptide)
    return Matches


CandidatePeptides = ["TMIA", "TAIM", "IAMT", "MTAL", "MAIT", "MLAT"]
TargetSpectrum = [0, 71, 101, 113, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416]
print(*FindMatchingPeptides(CandidatePeptides, TargetSpectrum))
# Output: IAMT MAIT
    
    
    

# Which of the following linear peptides is consistent with 
# Spectrum = {0 71 99 101 103 128 129 199 200 204 227 230 231 298 303 328 330 332 333}?
# -------------------------------------------------------------------------------------

AminoAcidMasses = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115,
                   128, 129, 131, 137, 147, 156, 163, 186]

MassTable = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
    'T': 101, 'C': 103, 'I': 113, 'L': 113,
    'N': 114, 'D': 115, 'K': 128, 'Q': 128,
    'E': 129, 'M': 131, 'H': 137, 'F': 147,
    'R': 156, 'Y': 163, 'W': 186
}

def LinearSpectrum(Peptide):
    PrefixMass = [0]
    for AA in Peptide:
        PrefixMass.append(PrefixMass[-1] + MassTable[AA])
    Spectrum = [0]
    for I in range(len(Peptide)):
        for J in range(I + 1, len(Peptide) + 1):
            Spectrum.append(PrefixMass[J] - PrefixMass[I])
    return sorted(Spectrum)

def IsConsistent(Peptide, Spectrum):
    PeptideSpectrum = LinearSpectrum(Peptide)
    SpectrumCopy = Spectrum.copy()
    for Mass in PeptideSpectrum:
        if Mass in SpectrumCopy:
            SpectrumCopy.remove(Mass)
        else:
            return False
    return True

def FindConsistentPeptides(Peptides, Spectrum):
    Matches = []
    for Peptide in Peptides:
        if IsConsistent(Peptide, Spectrum):
            Matches.append(Peptide)
    return Matches



CandidatePeptides = ["ETC", "QCV", "TCQ", "TVQ", "TCE", "AVQ"]
TargetSpectrum = [0, 71, 99, 101, 103, 128, 129, 199, 200, 204, 227, 230, 231, 298, 303, 328, 330, 332, 333]
print(*FindConsistentPeptides(CandidatePeptides, TargetSpectrum))
# Output: ETC TCQ TVQ

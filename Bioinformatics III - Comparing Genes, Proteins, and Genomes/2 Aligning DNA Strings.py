import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics III - Comparing Genes, Proteins, and Genomes/Data")


# -----------------------------------------------
# Global Alignment
# -----------------------------------------------

def GlobalAlignment(s1, s2, match, mismatch, indel):
    '''
    Input: A match reward, a mismatch penalty, an indel penalty, and two nucleotide strings.
    Output: The maximum alignment score of these strings followed by an alignment achieving this maximum score.
    '''

    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    backtrack = [[None] * (n + 1) for _ in range(m + 1)]

    # Initialize first row and column with indel penalties
    for i in range(1, m + 1):
        dp[i][0] = -i * indel
        backtrack[i][0] = 'U'  # Up
    for j in range(1, n + 1):
        dp[0][j] = -j * indel
        backtrack[0][j] = 'L'  # Left

    # Fill in the DP table and the backtrack table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s1[i - 1] == s2[j - 1]:
                score_diag = dp[i - 1][j - 1] + match
            else:
                score_diag = dp[i - 1][j - 1] - mismatch
            score_up = dp[i - 1][j] - indel
            score_left = dp[i][j - 1] - indel

            dp[i][j] = max(score_diag, score_up, score_left)

            if dp[i][j] == score_diag:
                backtrack[i][j] = 'D'  # Diagonal
            elif dp[i][j] == score_up:
                backtrack[i][j] = 'U'
            else:
                backtrack[i][j] = 'L'

    # Traceback to get aligned strings
    aligned_s1 = []
    aligned_s2 = []
    i, j = m, n
    while i > 0 or j > 0:
        if backtrack[i][j] == 'D':
            aligned_s1.append(s1[i - 1])
            aligned_s2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif backtrack[i][j] == 'U':
            aligned_s1.append(s1[i - 1])
            aligned_s2.append('-')
            i -= 1
        else:  # 'L'
            aligned_s1.append('-')
            aligned_s2.append(s2[j - 1])
            j -= 1

    # Reverse the aligned strings
    aligned_s1.reverse()
    aligned_s2.reverse()

    # Return score and alignments
    return dp[m][n], ''.join(aligned_s1), ''.join(aligned_s2)



# Examples
# -----------
score, a1, a2 = GlobalAlignment("GAGA", "GAT", 1, 1, 2)
print(score)
print(a1)
print(a2)




with open("dataset_30199_3.txt", "r") as file:
    match, mismatch, indel = map(int, file.readline().split())
    s1 = file.readline().strip()
    s2 = file.readline().strip()
    
score, a1, a2 = GlobalAlignment(s1, s2, match, mismatch, indel)
with open("output.txt", "w") as out_file:
    out_file.write(f"{score}\n")
    out_file.write(f"{a1}\n")
    out_file.write(f"{a2}\n")

# Output: 13
# AAGTGATTTGACGATTTACATAAATGGCAATATATCGTTCCACTTTTATCGCCCCACGGAGTGCAGGATTCCAATCC...
# ---T-A-GTG-CGATTTACATAAATGGCAATATATCGTCCCACTTTTGT-G---CA-GGATTCCA--A-T-C-GCGG...




with open("dataset_30202_14.txt", "r") as file:
    match, mismatch, gap = map(int, file.readline().split())
    s1 = file.readline().strip()
    s2 = file.readline().strip()

# Shorter = v, longer = w
if len(s1) >= len(s2):
    v, w = s2, s1
else:
    v, w = s1, s2

score, aligned_v, aligned_w  = GlobalAlignment(v, w, match, mismatch, gap)

with open("output.txt", "w") as out_file:
    out_file.write(f"{score}\n")
    out_file.write(f"{aligned_v}\n")
    out_file.write(f"{aligned_w}\n")



# -----------------------------------------------
# Local Alignment
# -----------------------------------------------

def ReadScoringMatrix(filename):
    with open(filename) as f:
        lines = f.readlines()
        header = lines[0].split()
        matrix = {}
        for line in lines[1:]:
            parts = line.strip().split()
            row_char = parts[0]
            scores = list(map(int, parts[1:]))
            for col_char, score in zip(header, scores):
                matrix[(row_char, col_char)] = score
        return matrix
    

def LocalAlignment(s1, s2, scoreMatrix, indel):
    '''
    Input: A scoring matrix (dict), an indel penalty (int), and two protein strings.
    Output: The maximum local alignment score and one optimal alignment.
    '''
    
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    backtrack = [[None] * (n + 1) for _ in range(m + 1)]

    maxScore = 0
    maxPos = (0, 0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            char1 = s1[i - 1]
            char2 = s2[j - 1]
            score_diag = dp[i - 1][j - 1] + scoreMatrix[(char1, char2)]
            score_up = dp[i - 1][j] - indel
            score_left = dp[i][j - 1] - indel
            dp[i][j] = max(0, score_diag, score_up, score_left)

            if dp[i][j] == 0:
                backtrack[i][j] = None
            elif dp[i][j] == score_diag:
                backtrack[i][j] = 'D'
            elif dp[i][j] == score_up:
                backtrack[i][j] = 'U'
            else:
                backtrack[i][j] = 'L'

            if dp[i][j] > maxScore:
                maxScore = dp[i][j]
                maxPos = (i, j)

    # Traceback from max score position
    aligned_s1 = []
    aligned_s2 = []
    i, j = maxPos
    while dp[i][j] != 0:
        if backtrack[i][j] == 'D':
            aligned_s1.append(s1[i - 1])
            aligned_s2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif backtrack[i][j] == 'U':
            aligned_s1.append(s1[i - 1])
            aligned_s2.append('-')
            i -= 1
        else:  # 'L'
            aligned_s1.append('-')
            aligned_s2.append(s2[j - 1])
            j -= 1

    aligned_s1.reverse()
    aligned_s2.reverse()

    return maxScore, ''.join(aligned_s1), ''.join(aligned_s2)


# Examples
# -----------------
scoreMatrix = ReadScoringMatrix("PAM250.txt")
indel = 5

s1 = "MEANLY"
s2 = "PENALTY"
score, a1, a2 = LocalAlignment(s1, s2, scoreMatrix, indel)
print(score, a1, s2)
# Output: 15 
# EANL-Y 
# PENALTY




# -----------------------------------------------
# Edit Distance
# -----------------------------------------------

def EditDistance(s1, s2):
    m, n = len(s1), len(s2)

    # Create DP table with (m+1) x (n+1)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize base cases
    for i in range(m + 1):
        dp[i][0] = i  # cost of deleting all characters of s1
    for j in range(n + 1):
        dp[0][j] = j  # cost of inserting all characters of s2

    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s1[i - 1] == s2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]  # no operation needed
            else:
                dp[i][j] = min(
                    dp[i - 1][j] + 1,    # deletion
                    dp[i][j - 1] + 1,    # insertion
                    dp[i - 1][j - 1] + 1 # substitution
                )
    return dp[m][n]



# Examples
# -----------
s1 = "GAGA"
s2 = "GAT"
EditDistance(s1, s2)
# Output: 2


with open("dataset_30200_3.txt", "r") as file:
    s1 = file.readline()
    s2 = file.readline()
EditDistance(s1, s2)
# Output: 317



# -----------------------------------------------
# Fitting Alignment
# -----------------------------------------------

def ReadScoringMatrix(filename):
    with open(filename) as f:
        lines = f.readlines()
        header = lines[0].split()
        matrix = {}
        for line in lines[1:]:
            parts = line.strip().split()
            row_char = parts[0]
            scores = list(map(int, parts[1:]))
            for col_char, score in zip(header, scores):
                matrix[(row_char, col_char)] = score
        return matrix


def FittingAlignment(v, w, scoreMatrix, indel):
    m, n = len(v), len(w)
    dp = [[float('-inf')] * (n + 1) for _ in range(m + 1)]
    backtrack = [[None] * (n + 1) for _ in range(m + 1)]

    # Initialization: free gaps at start of v, so first column is 0 for all i
    for i in range(m + 1):
        dp[i][0] = 0

    # Fill DP matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i-1][j-1] + scoreMatrix[(v[i-1], w[j-1])]
            delete = dp[i-1][j] - indel
            insert = dp[i][j-1] - indel
            max_score = max(match, delete, insert)
            dp[i][j] = max_score

            if max_score == match:
                backtrack[i][j] = 'D'  # diagonal
            elif max_score == delete:
                backtrack[i][j] = 'U'  # up
            else:
                backtrack[i][j] = 'L'  # left

    # Find max score in last column j=n (w fully aligned)
    max_score = float('-inf')
    max_i = None
    for i in range(1, m + 1):
        if dp[i][n] > max_score:
            max_score = dp[i][n]
            max_i = i

    # Traceback from (max_i, n)
    v_aligned = []
    w_aligned = []
    i, j = max_i, n
    while j > 0:
        if backtrack[i][j] == 'D':
            v_aligned.append(v[i-1])
            w_aligned.append(w[j-1])
            i -= 1
            j -= 1
        elif backtrack[i][j] == 'U':
            v_aligned.append(v[i-1])
            w_aligned.append('-')
            i -= 1
        else:  # 'L'
            v_aligned.append('-')
            w_aligned.append(w[j-1])
            j -= 1

    # Reverse aligned strings
    v_aligned = v_aligned[::-1]
    w_aligned = w_aligned[::-1]

    return max_score, ''.join(v_aligned), ''.join(w_aligned)


# Example
# ------------
scoreMatrix = ReadScoringMatrix("BLOSUM62.txt")
indel = 1
v = "DISCREPANTLY"
w = "PATENT"

score, v_aligned, w_aligned = FittingAlignment(v, w, scoreMatrix, indel)
print(score)
print(v_aligned)
print(w_aligned)
# Output: 20
# PA--NT
# PATENT



# -----------------------------------------------
# Overlap Alignment
# -----------------------------------------------

def OverlapAlignment(v, w, match, mismatch, indel):
    n, m = len(v), len(w)

    # DP and backtrack matrices
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    backtrack = [[None] * (m + 1) for _ in range(n + 1)]

    # Initialize
    for j in range(1, m + 1):
        dp[0][j] = -indel * j  # must align prefix of w
    for i in range(1, n + 1):
        dp[i][0] = 0  # allow alignment to start at any suffix of v

    # Fill DP matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if v[i - 1] == w[j - 1]:
                score = match
            else:
                score = -mismatch

            options = [
                (dp[i - 1][j] - indel, 'up'),
                (dp[i][j - 1] - indel, 'left'),
                (dp[i - 1][j - 1] + score, 'diag')
            ]

            dp[i][j], backtrack[i][j] = max(options)

    # Find the best alignment score in the last column (w is fully aligned)
    max_score = float('-inf')
    max_i = -1
    for i in range(1, n + 1):
        if dp[i][m] > max_score:
            max_score = dp[i][m]
            max_i = i

    # Backtrack from (max_i, m)
    i, j = max_i, m
    v_aligned = []
    w_aligned = []

    while j > 0:
        move = backtrack[i][j]
        if move == 'diag':
            v_aligned.append(v[i - 1])
            w_aligned.append(w[j - 1])
            i -= 1
            j -= 1
        elif move == 'up':
            v_aligned.append(v[i - 1])
            w_aligned.append('-')
            i -= 1
        elif move == 'left':
            v_aligned.append('-')
            w_aligned.append(w[j - 1])
            j -= 1

    v_aligned = ''.join(reversed(v_aligned))
    w_aligned = ''.join(reversed(w_aligned))

    return max_score, v_aligned, w_aligned




# Examples
# ---------
match, mismatch, indel = 1, 1, 2
v, w = "GAGA", "GAT"

score, v_aligned, w_aligned = OverlapAlignment(v, w, match, mismatch, indel)
print(score)
print(v_aligned)
print(w_aligned)
# Output: 1
# GAG
# GAT



with open("dataset_30200_7.txt", "r") as file:
    match, mismatch, indel = map(int, file.readline().split()) 
    v = file.readline().strip()
    w = file.readline().strip()
    
score, v_aligned, w_aligned = OverlapAlignment(v, w, match, mismatch, indel)

with open("output.txt", "w") as out_file:
    out_file.write(f"{score}\n")
    out_file.write(f"{v_aligned}\n")
    out_file.write(f"{w_aligned}\n")

# Output: 149
# TAG-TCTGAAGGT--CCTAACTGAGGTGAGTAGC-GGAAGG-GCGA-T-GCCCAACCACCCT----ATCAGGAAGC...
# GAGGTCTGAAGTTAGCCTAAAAGAGGTGAGTAGAAGGGCGATGCCAATAGAAAAACCACCGTGGTTATCAGGAAGC...


s1, s2 = "GATACACT", "ACGACCACAGATACCGCTATTCACTATATCGTT"
match, mismatch, indel = 1, 1, 1

score, v_aligned, w_aligned = OverlapAlignment(v, w, match, mismatch, indel)
print(score)
# Output: 1


# -----------------------------------------------
# Score an Alignment
# -----------------------------------------------

def AlignmentScore(s1, s2, match, mismatch, indel):
    score = 0
    for a, b in zip(s1, s2):
        if a == '-' and b == '-':
            # both gaps, no score
            continue
        elif a == b:
            # match reward
            score += match
        elif a == '-' or b == '-':
            score -= indel
        else:
            score -= mismatch
    return score


# Examples
# ----------
match, mismatch, indel = 1, 1, 2
s1, s2 = "TCGAC--ATT", "CC---GAA-T"
print(AlignmentScore(s1, s2, match, mismatch, indel))
# Output: -10


# TACTATTTACAGTAGACACGT
# AACAGAC-ATAC-AGATACCT
# What is the score of the bold portion of this alignment as a local alignment 
# if the match score is 1, the mismatch penalty is 3, and the indel penalty is 1?
# ----------------------------------------------------------------
match, mismatch, indel = 1, 3, 1
s1, s2 = "ACAGTAGACAC", "ATAC-AGATAC"
print(AlignmentScore(s1, s2, match, mismatch, indel))
# Output: -3


# Say that the match score is 1, the mismatch penalty is 0, and the indel penalty is 2. 
# Score the following overlap alignment.
# AGTACATCAGAGGAGTT-ACATACTAACG
#              AGTTCACAGGCTA-CGTACAGATATTACGACAGGCAGA
 # ----------------------------------------------------------------
match, mismatch, indel = 1, 0, 2
s1, s2 = "AGTT-ACATACTAACG", "AGTTCACAGGCTA-CG"
print(AlignmentScore(s1, s2, match, mismatch, indel))
# Output: 8





# -----------------------------------------------
# Score for Two Unaligned Strings
# -----------------------------------------------

def AlignmentScore(s1, s2, match, mismatch, indel):
    score = 0
    for a, b in zip(s1, s2):
        if a == '-' and b == '-':
            # both gaps, no score
            continue
        elif a == b:
            # match
            score += match
        elif a == '-' or b == '-':
            score -= indel
        else:
            score -= mismatch
    return score


def FittingAlignment(s1, s2, match, mismatch, indel):
    """
    match, mismatch, indel: positive integers (penalty magnitudes)
    match is a positive reward,
    mismatch and indel are penalties, so they are subtracted internally.
    """
    m, n = len(s1), len(s2)
    dp = [[float('-inf')] * (n + 1) for _ in range(m + 1)]
    
    # Initialization
    for j in range(n + 1):
        dp[0][j] = 0  # free to start anywhere in s2
    
    for i in range(1, m + 1):
        dp[i][0] = dp[i-1][0] - indel  # gaps in s2 to align s1 prefix
    
    # Fill dp table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Calculate score for match or mismatch
            if s1[i-1] == s2[j-1]:
                score_sub = match
            else:
                score_sub = -mismatch  # subtract mismatch penalty
            
            dp[i][j] = max(
                dp[i-1][j-1] + score_sub,  # match/mismatch
                dp[i-1][j] - indel,        # deletion (gap in s2)
                dp[i][j-1] - indel         # insertion (gap in s1)
            )
    
    # The best fitting alignment score is max in last row (whole s1 aligned)
    max_score = max(dp[m])
    return max_score


# Example 
# ----------
s1, s2 = "GATACACT", "ACGACCACAGATACCGCTATTCACTATATCGTT"
match, mismatch, indel = 1, 1, 1

print(FittingAlignment(s1, s2, match, mismatch, indel))
# Output: 5

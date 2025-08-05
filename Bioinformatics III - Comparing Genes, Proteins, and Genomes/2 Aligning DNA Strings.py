import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics III - Comparing Genes, Proteins, and Genomes/Data")


# -----------------------------------------------
# Global Alignment
# -----------------------------------------------

def GlobalAlignment(match, mismatch, indel, s1, s2):
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
score, a1, a2 = GlobalAlignment(1, 1, 2, "GAGA", "GAT")
print(score)
print(a1)
print(a2)



with open("dataset_30199_3.txt", "r") as file:
    match, mismatch, indel = map(int, file.readline().split())
    s1 = file.readline().strip()
    s2 = file.readline().strip()
    
score, a1, a2 = GlobalAlignment(match, mismatch, indel, s1, s2)
with open("output.txt", "w") as out_file:
    out_file.write(f"{score}\n")
    out_file.write(f"{a1}\n")
    out_file.write(f"{a2}\n")

# Output: 13
# AAGTGATTTGACGATTTACATAAATGGCAATATATCGTTCCACTTTTATCGCCCCACGGAGTGCAGGATTCCAATCC...
# ---T-A-GTG-CGATTTACATAAATGGCAATATATCGTCCCACTTTTGT-G---CA-GGATTCCA--A-T-C-GCGG...




    
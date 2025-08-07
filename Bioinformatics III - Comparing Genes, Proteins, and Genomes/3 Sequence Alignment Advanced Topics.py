import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics III - Comparing Genes, Proteins, and Genomes/Data")


# -----------------------------------------------
# Global Alignment with Affine Gap Penalties
# -----------------------------------------------

def AffineAlignmentGlobal(match, mismatch, gap_open, gap_extend, v, w):
    '''Gap: contiguous sequence of spaces in a row of an alignment
    
    Input: A match reward, a mismatch penalty, a gap opening penalty, a gap extension penalty, and two nucleotide strings.
    Output: The maximum alignment score between v and w, followed by an alignment of v and w achieving this maximum score.
    '''
    
    n, m_len = len(v), len(w)

    # Initialize matrices
    middle = [[float('-inf')] * (m_len + 1) for _ in range(n + 1)]  # Match/mismatch
    lower = [[float('-inf')] * (m_len + 1) for _ in range(n + 1)]  # Gap in w (vertical)
    upper = [[float('-inf')] * (m_len + 1) for _ in range(n + 1)]  # Gap in v (horizontal)

    # Initialization
    middle[0][0] = 0
    for i in range(1, n + 1):
        lower[i][0] = -gap_open - (i - 1) * gap_extend
        middle[i][0] = lower[i][0]
    for j in range(1, m_len + 1):
        upper[0][j] = -gap_open - (j - 1) * gap_extend
        middle[0][j] = upper[0][j]

    # Fill in matrices
    for i in range(1, n + 1):
        for j in range(1, m_len + 1):
            s = match if v[i - 1] == w[j - 1] else -mismatch

            lower[i][j] = max(lower[i - 1][j] - gap_extend, middle[i - 1][j] - gap_open - gap_extend)
            upper[i][j] = max(upper[i][j - 1] - gap_extend, middle[i][j - 1] - gap_open - gap_extend)
            middle[i][j] = max(middle[i - 1][j - 1] + s, lower[i][j], upper[i][j])


    final_score = max(middle[n][m_len], lower[n][m_len], upper[n][m_len])
    if final_score == lower[n][m_len]:
        current_matrix = 'lower'
    elif final_score == upper[n][m_len]:
        current_matrix = 'upper'
    else:
        current_matrix = 'middle'

    # Backtrack
    aligned_v = []
    aligned_w = []
    i, j = n, m_len

    while i > 0 or j > 0:
        if current_matrix == 'middle':
            s = match if v[i - 1] == w[j - 1] else -mismatch
            if i > 0 and j > 0 and middle[i][j] == middle[i - 1][j - 1] + s:
                aligned_v.append(v[i - 1])
                aligned_w.append(w[j - 1])
                i -= 1
                j -= 1
                current_matrix = 'middle'
            elif i > 0 and middle[i][j] == lower[i][j]:
                current_matrix = 'lower'
            elif j > 0 and middle[i][j] == upper[i][j]:
                current_matrix = 'upper'
        elif current_matrix == 'lower':
            aligned_v.append(v[i - 1])
            aligned_w.append('-')
            if lower[i][j] == lower[i - 1][j] - gap_extend:
                current_matrix = 'lower'
            else:
                current_matrix = 'middle'
            i -= 1
        elif current_matrix == 'upper':
            aligned_v.append('-')
            aligned_w.append(w[j - 1])
            if upper[i][j] == upper[i][j - 1] - gap_extend:
                current_matrix = 'upper'
            else:
                current_matrix = 'middle'
            j -= 1

    return final_score, ''.join(reversed(aligned_v)), ''.join(reversed(aligned_w))


# Example
# -----------
match, mismatch, gap_open, gap_extend = 1, 3, 2, 1
v, w = "GA", "GTTA"

score, aligned_v, aligned_w = AffineAlignmentGlobal(match, mismatch, gap_open, gap_extend, v, w)
print(score)
print(aligned_v)
print(aligned_w)

# Output: -2
# G--A
# GTTA



# -----------------------------------------------
# Score Global Alignment with Affine Gap Penalties
# -----------------------------------------------

def AffineAlignmentScore(aligned_v, aligned_w, match, mismatch, gapOpen, gapExtend):
    '''
    aligned_v, aligned_w: aligned sequences of equal length (strings)
    match: positive score for match (e.g., 1)
    mismatch: penalty for mismatch (e.g., 1 means -1 penalty)
    gapOpen: penalty for gap opening (e.g., 4 means -4 penalty)
    gapExtend: penalty for gap extension (e.g., 1 means -1 penalty)
    '''

    score = 0
    in_gap_v = False
    in_gap_w = False

    for i in range(len(aligned_v)):
        a = aligned_v[i]
        b = aligned_w[i]

        if a != '-' and b != '-':
            # Both nucleotides
            if a == b:
                score += match
            else:
                score -= mismatch
            in_gap_v = False
            in_gap_w = False

        elif a == '-' and b != '-':
            # Gap in aligned_v
            if not in_gap_v:
                score -= gapOpen
                in_gap_v = True
            else:
                score -= gapExtend
            in_gap_w = False

        elif b == '-' and a != '-':
            # Gap in aligned_w
            if not in_gap_w:
                score -= gapOpen
                in_gap_w = True
            else:
                score -= gapExtend
            in_gap_v = False

        else:
            # Both gaps - no penalty or reward
            pass

    return score




# Examples
# ---------
aligned_v, aligned_w = "TCGAC--ATT", "CC---GAA-T"
match, mismatch, gapOpen, gapExtend = 1, 1, 4, 1
print(AffineAlignmentScore(aligned_v, aligned_w, match, mismatch, gapOpen, gapExtend))
# Output: -13


aligned_v, aligned_w = "A-C--GTTAC", "ATGCAG---T"
match, mismatch, gapOpen, gapExtend = 1, 1, 4, 1
print(AffineAlignmentScore(aligned_v, aligned_w, match, mismatch, gapOpen, gapExtend))
# Output: -15




# -----------------------------------------------
# Local with Affine Gap Penalties
# -----------------------------------------------

def AffineAlignmentLocal(match, mismatch, gapOpen, gapExtend, v, w):
    '''
    Compute optimal local alignment with affine gap penalties.
    
    Input: A match reward, a mismatch penalty, a gap opening penalty, a gap extension penalty, and two nucleotide strings.
    Output: The maximum local alignment score between v and w, followed by an alignment achieving this score.
    Aligned v and w are the same string for local alignment.
    '''

    n, m_len = len(v), len(w)

    # Initialize matrices
    middle = [[0] * (m_len + 1) for _ in range(n + 1)]  # Match/mismatch
    lower = [[0] * (m_len + 1) for _ in range(n + 1)]  # Gap in w (vertical)
    upper = [[0] * (m_len + 1) for _ in range(n + 1)]  # Gap in v (horizontal)

    max_score = 0
    max_pos = (0, 0)
    max_matrix = 'middle'

    for i in range(1, n + 1):
        for j in range(1, m_len + 1):
            s = match if v[i - 1] == w[j - 1] else -mismatch

            # Lower matrix (gap in w / vertical gap)
            lower[i][j] = max(
                lower[i - 1][j] - gapExtend,
                middle[i - 1][j] - gapOpen - gapExtend,
                0
            )

            # Upper matrix (gap in v / horizontal gap)
            upper[i][j] = max(
                upper[i][j - 1] - gapExtend,
                middle[i][j - 1] - gapOpen - gapExtend,
                0
            )

            # Middle matrix (match/mismatch)
            middle[i][j] = max(
                middle[i - 1][j - 1] + s,
                lower[i][j],
                upper[i][j],
                0
            )

            # Keep track of max score
            if middle[i][j] >= max_score:
                max_score = middle[i][j]
                max_pos = (i, j)
                if middle[i][j] == lower[i][j]:
                    max_matrix = 'lower'
                elif middle[i][j] == upper[i][j]:
                    max_matrix = 'upper'
                else:
                    max_matrix = 'middle'

    # Backtracking
    aligned_v = []
    aligned_w = []
    i, j = max_pos
    current_matrix = max_matrix

    while i > 0 or j > 0:
        if current_matrix == 'middle':
            if middle[i][j] == 0:
                break
            s = match if v[i - 1] == w[j - 1] else -mismatch
            if i > 0 and j > 0 and middle[i][j] == middle[i - 1][j - 1] + s:
                aligned_v.append(v[i - 1])
                aligned_w.append(w[j - 1])
                i -= 1
                j -= 1
                current_matrix = 'middle'
            elif i > 0 and middle[i][j] == lower[i][j]:
                current_matrix = 'lower'
            elif j > 0 and middle[i][j] == upper[i][j]:
                current_matrix = 'upper'
        elif current_matrix == 'lower':
            if lower[i][j] == 0:
                break
            aligned_v.append(v[i - 1])
            aligned_w.append('-')
            if lower[i][j] == lower[i - 1][j] - gapExtend:
                current_matrix = 'lower'
            else:
                current_matrix = 'middle'
            i -= 1
        elif current_matrix == 'upper':
            if upper[i][j] == 0:
                break
            aligned_v.append('-')
            aligned_w.append(w[j - 1])
            if upper[i][j] == upper[i][j - 1] - gapExtend:
                current_matrix = 'upper'
            else:
                current_matrix = 'middle'
            j -= 1

    return max_score, ''.join(reversed(aligned_v)), ''.join(reversed(aligned_w))


# Example
# ---------
match, mismatch, gapOpen, gapExtend = 1, 3, 2, 1
v, w = "ACGTGCAACTGAC", "ACTGACAGT"

score, aligned_v, aligned_w = AffineAlignmentLocal(match, mismatch, gapOpen, gapExtend, v, w)
print(score)
print(aligned_v)
print(aligned_w)

# Output: 6
# ACTGAC
# ACTGAC



# -----------------------------------------------
# Middle Node(s) of Longest Path(s)
# -----------------------------------------------

def FindMiddleNode(v, w):
    n = len(v)
    m = len(w)
    middle = m // 2

    def ForwardDP():
        previous = [0] * (n + 1)
        for j in range(1, middle + 1):
            current = [0] * (n + 1)
            for i in range(1, n + 1):
                if v[i - 1] == w[j - 1]:
                    current[i] = previous[i - 1] + 1
                else:
                    current[i] = max(previous[i], current[i - 1])
            previous = current
        return previous

    def ReverseDP():
        reversed_v = v[::-1]
        reversed_w = w[middle:][::-1]
        prev = [0] * (n + 1)
        for j in range(1, len(reversed_w) + 1):
            curr = [0] * (n + 1)
            for i in range(1, n + 1):
                if reversed_v[i - 1] == reversed_w[j - 1]:
                    curr[i] = prev[i - 1] + 1
                else:
                    curr[i] = max(prev[i], curr[i - 1])
            prev = curr
        return prev[::-1]

    forward = ForwardDP()
    reverse = ReverseDP()

    max_len = -1
    middle_nodes = []

    for i in range(len(forward)):
        total_len = forward[i] + reverse[i]
        if total_len > max_len:
            max_len = total_len
            middle_nodes = [(i, middle)]
        elif total_len == max_len:
            middle_nodes.append((i, middle))

    return middle_nodes



# Example
# --------
v, w = "ACCA", "CAAC"
print(*FindMiddleNode(v, w))
# Output: (0, 2) (1, 2) (2, 2) (3, 2) (4, 2)




# -----------------------------------------------
# Find Middle Edge (in Linear Space)
# -----------------------------------------------

def FindMiddleEdge(v, w, matchScore, mismatchPenalty, gapPenalty):
    ''' 
    Input: Two nucleotide strings, A match score, a mismatch penalty, and an indel penalty.
    Output: A middle edge in the alignment graph.
    
    Note: Enter the shorter string first (v), then the longer string (w).
    '''

    def Score(a, b):
        return matchScore if a == b else -mismatchPenalty

    def ComputeMiddleColumnScores(s1, s2):
        previous = [i * -gapPenalty for i in range(len(s1) + 1)]
        for j in range(1, len(s2) + 1):
            current = [j * -gapPenalty]
            for i in range(1, len(s1) + 1):
                match = previous[i - 1] + Score(s1[i - 1], s2[j - 1])
                delete = previous[i] - gapPenalty
                insert = current[i - 1] - gapPenalty
                current.append(max(match, delete, insert))
            previous = current
        return previous

    def ReverseString(s):
        return s[::-1]

    lenV, lenW = len(v), len(w)
    mid = lenW // 2

    scoreLeft = ComputeMiddleColumnScores(v, w[:mid])
    scoreRight = ComputeMiddleColumnScores(ReverseString(v), ReverseString(w[mid:]))

    totalScores = [l + r for l, r in zip(scoreLeft, reversed(scoreRight))]
    maxIndex = totalScores.index(max(totalScores))

    startNode = (maxIndex, mid)

    diagonalScore = (
        scoreLeft[maxIndex]
        + (Score(v[maxIndex], w[mid]) if maxIndex < lenV and mid < lenW else float('-inf'))
        + (scoreRight[lenV - (maxIndex + 1)] if maxIndex < lenV and mid < lenW else float('-inf'))
    )

    downScore = (
        (scoreLeft[maxIndex + 1] - gapPenalty if maxIndex + 1 <= lenV else float('-inf'))
        + (scoreRight[lenV - (maxIndex + 1)] if maxIndex + 1 <= lenV else float('-inf'))
    )

    rightScore = (
        (scoreLeft[maxIndex] - gapPenalty if mid + 1 <= lenW else float('-inf'))
        + (scoreRight[lenV - maxIndex] if mid + 1 <= lenW else float('-inf'))
    )

    bestScore = max(diagonalScore, downScore, rightScore)
    if bestScore == diagonalScore:
        endNode = (maxIndex + 1, mid + 1)
    elif bestScore == downScore:
        endNode = (maxIndex + 1, mid)
    else:
        endNode = (maxIndex, mid + 1)

    return startNode, endNode




# Examples
# ----------
match, mismatch, gap = 1, 1, 2
v, w = "GAT", "GAGA"

start, end = FindMiddleEdge(v, w, match, mismatch, gap)
print(f"{start[0]} {start[1]}")
print(f"{end[0]} {end[1]}")
# Output:
# 2 2
# 2 3


with open("dataset_30202_12.txt", "r") as file:
    match, mismatch, gap = map(int, file.readline().split()) 
    v = file.readline().strip()
    w = file.readline().strip()

middle_edge = FindMiddleEdge(w, v, match, mismatch, gap)
    
print(f"{middle_edge[0][0]} {middle_edge[0][1]}")
print(f"{middle_edge[1][0]} {middle_edge[1][1]}")




# -----------------------------------------------
# Find Global Alignment (in Linear Space)
# -----------------------------------------------

def GlobalAlignmentLinearSpace(v, w, match, mismatch, gap):
    """
    Implement LinearSpaceAlignment to solve the Global Alignment Problem for a large dataset.

    Input: A match reward, a mismatch penalty, an indel penalty, and two long (2000 nucleotide) DNA strings.
    Output: The maximum alignment score of these strings, followed by an alignment achieving this maximum score.
    """

    # Scoring helper function
    def Score(a, b):
        if a == b:
            return match
        else:
            return -mismatch

    # Compute score vectors for linear space DP from top to bottom for substring v[top:bottom], w[left:right]
    def ScoreVector(v, w, top, bottom, left, right):
        # DP array: previous and current rows
        prev = [i * -gap for i in range(bottom - top + 1)]
        for j in range(left + 1, right + 1):
            current = [0] * (bottom - top + 1)
            current[0] = (j - left) * -gap
            for i in range(1, bottom - top + 1):
                match_mismatch = prev[i - 1] + Score(v[top + i - 1], w[j - 1])
                delete = prev[i] - gap
                insert = current[i - 1] - gap
                current[i] = max(match_mismatch, delete, insert)
            prev = current
        return prev

    # Find the middle node i in the middle column (middle = (left + right)//2)
    def MiddleNode(v, w, top, bottom, left, right):
        middle = (left + right) // 2
        scoreL = ScoreVector(v, w, top, bottom, left, middle)
        scoreR = ScoreVector(v[::-1], w[::-1], len(v) - bottom, len(v) - top, len(w) - right, len(w) - middle)
        scoreR = scoreR[::-1]
        max_score = float('-inf')
        max_i = top
        for i in range(bottom - top + 1):
            s = scoreL[i] + scoreR[i]
            if s > max_score:
                max_score = s
                max_i = top + i
        return max_i

    # Find the middle edge symbol →, ↓, or ↘ and the vertical coordinate of its initial node
    def MiddleEdge(v, w, top, bottom, left, right):
        middle = (left + right) // 2
        # Scores from left half
        scoreL = ScoreVector(v, w, top, bottom, left, middle)
        # Scores from right half (reverse)
        scoreR = ScoreVector(v[::-1], w[::-1], len(v) - bottom, len(v) - top, len(w) - right, len(w) - middle)
        scoreR = scoreR[::-1]

        max_score = float('-inf')
        max_i = top

        for i in range(bottom - top + 1):
            s = scoreL[i] + scoreR[i]
            if s > max_score:
                max_score = s
                max_i = top + i

        # Now determine which edge from node (max_i, middle) gives max score
        if max_i < bottom and middle < right:
            diag = scoreL[max_i - top] + scoreR[max_i - top + 1] + Score(v[max_i], w[middle])
        else:
            diag = float('-inf')
        if max_i < bottom:
            down = scoreL[max_i - top + 1] + scoreR[max_i - top] - gap
        else:
            down = float('-inf')
        if middle < right:
            right_score = scoreL[max_i - top] + scoreR[max_i - top] - gap
        else:
            right_score = float('-inf')

        max_edge_score = max(diag, down, right_score)
        if max_edge_score == diag:
            return '↘', max_i
        elif max_edge_score == down:
            return '↓', max_i
        else:
            return '→', max_i

    def LinearSpaceAlignment(v, w, top, bottom, left, right):
        if left == right:
            return 'V' * (bottom - top)
        if top == bottom:
            return 'H' * (right - left)

        middle = (left + right) // 2
        midEdge, midNode = MiddleEdge(v, w, top, bottom, left, right)
        left_path = LinearSpaceAlignment(v, w, top, midNode, left, middle)

        if midEdge == '→':
            mid_symbol = 'H'
        elif midEdge == '↓':
            mid_symbol = 'V'
        elif midEdge == '↘':
            mid_symbol = 'D'
        else:
            raise ValueError(f"Unknown middle edge symbol: {midEdge}")

        if midEdge in ('→', '↘'):
            middle += 1
        if midEdge in ('↓', '↘'):
            midNode += 1

        right_path = LinearSpaceAlignment(v, w, midNode, bottom, middle, right)
        return left_path + mid_symbol + right_path

    # Reconstruct alignment strings from path string
    def ReconstructAlignment(v, w, path):
        aligned_v = []
        aligned_w = []
        i = 0
        j = 0
        for p in path:
            if p == 'D':  # diagonal
                aligned_v.append(v[i])
                aligned_w.append(w[j])
                i += 1
                j += 1
            elif p == 'V':  # vertical (gap in w)
                aligned_v.append(v[i])
                aligned_w.append('-')
                i += 1
            elif p == 'H':  # horizontal (gap in v)
                aligned_v.append('-')
                aligned_w.append(w[j])
                j += 1
        return ''.join(aligned_v), ''.join(aligned_w)

    # Calculate final alignment score from aligned sequences
    def CalculateScore(aligned_v, aligned_w):
        score_sum = 0
        for a, b in zip(aligned_v, aligned_w):
            if a == '-' or b == '-':
                score_sum -= gap
            elif a == b:
                score_sum += match
            else:
                score_sum -= mismatch
        return score_sum

    # Run LinearSpaceAlignment on full strings
    path = LinearSpaceAlignment(v, w, 0, len(v), 0, len(w))

    aligned_v, aligned_w = ReconstructAlignment(v, w, path)
    max_score = CalculateScore(aligned_v, aligned_w)

    return max_score, aligned_v, aligned_w




# Examples
# ----------
match, mismatch, gap = 1, 1, 2
v, w = "GAGA", "GAT"

score, aligned_v, aligned_w  = GlobalAlignmentLinearSpace(v, w, match, mismatch, gap)
print(score)
print(aligned_v)
print(aligned_w)

# Output: -1
# GAGA
# GA-T




with open("dataset_30202_14.txt", "r") as file:
    match, mismatch, gap = map(int, file.readline().split())
    s1 = file.readline().strip()
    s2 = file.readline().strip()

# Shorter = v, longer = w
if len(s1) >= len(s2):
    v, w = s2, s1
else:
    v, w = s1, s2

score, aligned_v, aligned_w  = GlobalAlignmentLinearSpace(v, w, match, mismatch, gap)

with open("output.txt", "w") as out_file:
    out_file.write(f"{score}\n")
    out_file.write(f"{aligned_v}\n")
    out_file.write(f"{aligned_w}\n")




# -----------------------------------------------
# Multiple Longest Common Subsequence (Global Alignment)
# -----------------------------------------------

def MultipleLCSGlobal(s1, s2, s3):
    n1, n2, n3 = len(s1), len(s2), len(s3)
    
    # DP table and backtrack table
    dp = [[[0]*(n3+1) for _ in range(n2+1)] for _ in range(n1+1)]
    back = [[[None]*(n3+1) for _ in range(n2+1)] for _ in range(n1+1)]

    # Directions: 7 options (match, or gaps in any of the 3)
    directions = [
        (-1, -1, -1),  # match or mismatch
        (-1,  0,  0),
        ( 0, -1,  0),
        ( 0,  0, -1),
        (-1, -1,  0),
        (-1,  0, -1),
        ( 0, -1, -1)
    ]

    # Fill DP table
    for i in range(n1+1):
        for j in range(n2+1):
            for k in range(n3+1):
                best = -1
                move = None

                for di, dj, dk in directions:
                    ni, nj, nk = i + di, j + dj, k + dk
                    if ni < 0 or nj < 0 or nk < 0:
                        continue

                    score = dp[ni][nj][nk]
                    if di == dj == dk == -1 and s1[ni] == s2[nj] == s3[nk]:
                        score += 1

                    if score > best:
                        best = score
                        move = (di, dj, dk)

                if move:
                    dp[i][j][k] = best
                    back[i][j][k] = move

    # Backtrack to construct alignment
    i, j, k = n1, n2, n3
    a1, a2, a3 = [], [], []

    while i > 0 or j > 0 or k > 0:
        di, dj, dk = back[i][j][k]
        ni, nj, nk = i + di, j + dj, k + dk

        a1.append(s1[ni] if di == -1 else '-')
        a2.append(s2[nj] if dj == -1 else '-')
        a3.append(s3[nk] if dk == -1 else '-')

        i, j, k = ni, nj, nk

    # Reverse to get correct order
    aligned1 = ''.join(reversed(a1))
    aligned2 = ''.join(reversed(a2))
    aligned3 = ''.join(reversed(a3))

    return dp[n1][n2][n3], aligned1, aligned2, aligned3


# Examples
# ----------
s1 = "ATATCCG"
s2 = "TCCGA"
s3 = "ATGTACTG"

score, a1, a2, a3 = MultipleLCSGlobal(s1, s2, s3)

print(score)
print(a1)
print(a2)
print(a3)




with open("dataset_30203_5.txt", "r") as file:
    s1 = file.readline().strip()
    s2 = file.readline().strip()
    s3 = file.readline().strip()

score, a1, a2, a3 = MultipleLCSGlobal(s1, s2, s3)

print(score)
print(a1)
print(a2)
print(a3)

# Output: 3
# ---GT---AGTGGC-
# -GCGT---A---CCA
# --TTACAA---GC-



# -----------------------------------------------
# Multiple Longest Common Subsequence (Local Alignment)
# -----------------------------------------------

def MultipleLCSLocal(S1, S2, S3):
    n1, n2, n3 = len(S1), len(S2), len(S3)
    
    # DP and backtrack tables
    dp = [[[0] * (n3 + 1) for _ in range(n2 + 1)] for _ in range(n1 + 1)]
    back = [[[None] * (n3 + 1) for _ in range(n2 + 1)] for _ in range(n1 + 1)]

    # 7 directions: match/mismatch and possible gaps
    directions = [
        (-1, -1, -1),
        (-1,  0,  0),
        ( 0, -1,  0),
        ( 0,  0, -1),
        (-1, -1,  0),
        (-1,  0, -1),
        ( 0, -1, -1)
    ]

    max_score = 0
    max_pos = (0, 0, 0)

    # Fill DP table
    for i in range(1, n1 + 1):
        for j in range(1, n2 + 1):
            for k in range(1, n3 + 1):
                best = 0
                move = None

                for di, dj, dk in directions:
                    ni, nj, nk = i + di, j + dj, k + dk
                    if ni < 0 or nj < 0 or nk < 0:
                        continue

                    score = dp[ni][nj][nk]

                    # Add +1 for match
                    if di == dj == dk == -1 and S1[ni] == S2[nj] == S3[nk]:
                        score += 1

                    if score > best:
                        best = score
                        move = (di, dj, dk)

                dp[i][j][k] = best
                back[i][j][k] = move

                # Update max score and position
                if best > max_score:
                    max_score = best
                    max_pos = (i, j, k)

    # Backtrack from max position
    i, j, k = max_pos
    a1, a2, a3 = [], [], []

    while back[i][j][k] and dp[i][j][k] > 0:
        di, dj, dk = back[i][j][k]
        ni, nj, nk = i + di, j + dj, k + dk

        a1.append(S1[ni] if di == -1 else '-')
        a2.append(S2[nj] if dj == -1 else '-')
        a3.append(S3[nk] if dk == -1 else '-')

        i, j, k = ni, nj, nk

    # Reverse to get correct alignment
    aligned1 = ''.join(reversed(a1))
    aligned2 = ''.join(reversed(a2))
    aligned3 = ''.join(reversed(a3))

    return max_score, aligned1, aligned2, aligned3


def LCSBackTrack(v, w):
    '''Backtrack instead of going forward to avoid getting stuck in nodes.
    '''
    n, m = len(v), len(w)
    s = [[0] * (m + 1) for _ in range(n + 1)]
    backtrack = [[''] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        s[i][0] = 0
    for j in range(m + 1):
        s[0][j] = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = 1 if v[i - 1] == w[j - 1] else 0
            choices = [
                s[i - 1][j],                # down
                s[i][j - 1],                # right
                s[i - 1][j - 1] + match     # diagonal
            ]
            s[i][j] = max(choices)

            if s[i][j] == s[i - 1][j]:
                backtrack[i][j] = "↓"
            elif s[i][j] == s[i][j - 1]:
                backtrack[i][j] = "→"
            elif s[i][j] == s[i - 1][j - 1] + match:
                backtrack[i][j] = "↘"

    return backtrack


def OutputLCS(backtrack, v, i, j):
    '''
    Input: Two sequence strings.
    Output: A longest common subsequence. 
    
    Note: more than one solution may exist.
    '''
    
    if i == 0 or j == 0:
        return ""
    if backtrack[i][j] == "↓":
        return OutputLCS(backtrack, v, i - 1, j)
    elif backtrack[i][j] == "→":
        return OutputLCS(backtrack, v, i, j - 1)
    else: # diagonal ↘
        return OutputLCS(backtrack, v, i - 1, j - 1) + v[i - 1]



# Example
# ----------
s1, s2, s3 = "CCAATACGAC", "GCCTTACGCT", "CCCTAGCGGC"
score, a1, a2, a3 = MultipleLCSLocal(s1, s2, s3)
backtrack = LCSBackTrack(a1, a2)
lcs = ''.join(c1 for c1, c2, c3 in zip(a1, a2, a3) if c1 == c2 == c3 and c1 != '-')
print(lcs)
# Output: CCTACGC



s1, s2, s3 = "TCTAGCGAAC", "ATTACCGATC", "TTCACTGACG"

score, a1, a2, a3 = MultipleLCSLocal(s1, s2, s3)
backtrack = LCSBackTrack(a1, a2)
lcs = ''.join(c1 for c1, c2, c3 in zip(a1, a2, a3) if c1 == c2 == c3 and c1 != '-')
print(lcs)
# Output: TTACGAC



s1, s2, s3 = "CGGAACTGGT", "TGAGACGGTA", "TGCGACGGCT"
score, a1, a2, a3 = MultipleLCSLocal(s1, s2, s3)
backtrack = LCSBackTrack(a1, a2)
lcs = ''.join(c1 for c1, c2, c3 in zip(a1, a2, a3) if c1 == c2 == c3 and c1 != '-')
print(lcs)
# Output: GGACGGT



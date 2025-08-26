import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics VI - Finding Mutations in DNA and Proteins/Data")

# -----------------------------------------------
# Profile Hidden Markov Model (HMM)
# -----------------------------------------------

import numpy as np

def ConstructProfileHMM(theta, alphabet, alignment):
    '''Construct a profile HMM from a multiple alignment.
    
    Input: A threshold θ, followed by an alphabet Σ, followed by a multiple alignment 
    whose strings are formed from Σ.
    Output: The transition matrix followed by the emission matrix of HMM(Alignment, θ).'''

    n, l = alignment.shape
    threshold = theta * n
    kept = [True] * l
    for i in range(l):
        if sum(alignment[:, i] == '-') >= threshold:
            kept[i] = False

    levels = [0] * l
    for i in range(l):
        levels[i] = levels[i - 1] if i > 0 else 0
        if kept[i]:
            levels[i] += 1

    alphabetDict = ({alphabet[i]: i for i in range(len(alphabet))}, {i: alphabet[i] for i in range(len(alphabet))})

    def GetIndex(level, state):
        # states: 0 = I, 1 = M, 2 = D/E; level = profile column number (0-based)
        if level == 0:
            return 1 if state == 0 else 0
        if state != 0:
            return 3 * level - 2 + state
        else:
            return 3 * level + 1

    size = sum(kept) * 3 + 3
    transition = np.zeros((size, 3), dtype=float)
    emission = np.zeros((size, len(alphabet)), dtype=float)

    for i in range(n):
        lastLevel = 0
        lastState = -1
        lastInd = GetIndex(lastLevel, lastState)
        for j in range(l):
            currLevel = levels[j]
            if kept[j]:
                currState = 2 if alignment[i, j] == '-' else 1
                currInd = GetIndex(currLevel, currState)
                transition[lastInd, currState] += 1
                if currState == 1:
                    emission[currInd, alphabetDict[0][alignment[i, j]]] += 1
                lastInd = currInd
            else:
                if alignment[i, j] != '-':
                    currState = 0
                    currInd = GetIndex(currLevel, currState)
                    transition[lastInd, currState] += 1
                    emission[currInd, alphabetDict[0][alignment[i, j]]] += 1
                    lastInd = currInd
        transition[lastInd, 2] += 1

    # Normalize transitions and emissions
    for i in range(transition.shape[0]):
        s = sum(transition[i, :])
        if s != 0:
            transition[i, :] /= s
        s = sum(emission[i, :])
        if s != 0:
            emission[i, :] /= s

    # Build full transition matrix from the compact form
    fullTransition = np.zeros((size, size), dtype=float)
    fullTransition[0, 1:4] = transition[0, :]

    for i in range(1, size - 4):
        mod = i % 3
        if mod == 1:
            fullTransition[i, i:i + 3] = transition[i, :]
        elif mod == 2:
            fullTransition[i, i + 2:i + 5] = transition[i, :]
        else:
            fullTransition[i, i + 1:i + 4] = transition[i, :]

    for i in range(-4, -1):
        fullTransition[i, -2] = transition[i, 0]
        fullTransition[i, -1] = transition[i, -1]

    # Construct states list
    states = [''] * size
    states[0] = 'S'
    states[-1] = 'E'
    states[1] = 'I0'
    s = ('M', 'D', 'I')
    for i in range(2, size - 1):
        states[i] = s[(i + 1) % 3] + str((i + 1) // 3)

    return states, np.round(fullTransition, 3), np.round(emission, 3)

def ReadFromFile(filename):
    with open(filename, "r") as file:
        data = file.read().split()
    theta = float(data[0])
    ind = [i for i in range(len(data)) if data[i] == '--------']
    alphabet = data[ind[0] + 1:ind[1]]
    alignment = np.array([[*s] for s in data[ind[1] + 1:]])
    return theta, alphabet, alignment

def PrintTransitionAndEmission(alphabet, states, fullTransition, emission):
    print('\t' + '\t'.join(states))
    for i in range(fullTransition.shape[0]):
        print('\t'.join([states[i]] + list(map(str, fullTransition[i, :]))))
    print('--------')
    print('\t' + '\t'.join(alphabet))
    for i in range(emission.shape[0]):
        print('\t'.join([states[i]] + list(map(str, emission[i, :]))))

def WriteTransitionAndEmission(alphabet, states, fullTransition, emission):
    with open("output.txt", "w") as file:
        file.write('\t' + '\t'.join(states) + '\n')
        for i in range(fullTransition.shape[0]):
            file.write('\t'.join([states[i]] + list(map(str, fullTransition[i, :]))) + '\n')
        file.write('--------\n')
        file.write('\t' + '\t'.join(alphabet) + '\n')
        for i in range(emission.shape[0]):
            file.write('\t'.join([states[i]] + list(map(str, emission[i, :]))) + '\n')



# Example 1
# ----------
theta = 0.289
alphabet = ['A', 'B', 'C', 'D', 'E']
input_alignment = ["EBA", "E-D", "EB-", "EED", "EBD", "EBE", "E-D", "E-D"]
alignment = np.array([list(row) for row in input_alignment])

states, transition_matrix, emission_matrix = ConstructProfileHMM(theta, alphabet, alignment)
PrintTransitionAndEmission(alphabet, states, transition_matrix, emission_matrix)

# Output:
    # 	    S	I0	M1	D1	I1	M2	D2	I2	E
    # S	    0.0	0.0	1.0	0.0	0.0	0.0	0.0	0.0	0.0
    # I0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
    # M1	0.0	0.0	0.0	0.0	0.625	0.375	0.0	0.0	0.0
    # D1	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
    # I1	0.0	0.0	0.0	0.0	0.0	0.8	0.2	0.0	0.0
    # M2	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	1.0
    # D2	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	1.0
    # I2	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
    # E	    0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
    # --------
    # 	    A	B	C	D	E
    # S	    0.0	0.0	0.0	0.0	0.0
    # I0	0.0	0.0	0.0	0.0	0.0
    # M1	0.0	0.0	0.0	0.0	1.0
    # D1	0.0	0.0	0.0	0.0	0.0
    # I1	0.0	0.8	0.0	0.0	0.2
    # M2	0.143	0.0	0.0	0.714	0.143
    # D2	0.0	0.0	0.0	0.0	0.0
    # I2	0.0	0.0	0.0	0.0	0.0
    # E	    0.0	0.0	0.0	0.0	0.0
    
    
# Example 2
# ----------
theta, alphabet, alignment = ReadFromFile("dataset_30331_15.txt")
states, transition_matrix, emission_matrix = ConstructProfileHMM(theta, alphabet, alignment)
WriteTransitionAndEmission(alphabet, states, transition_matrix, emission_matrix)

# Output: 
    # 	S	I0	M1	D1	I1	M2	D2	I2	M3	D3	I3	M4	D4	I4	M5	D5	I5	M6	D6	I6	M7	D7	I7	M8	D8	I8	E
# S 	0.0	0.0	0.857	0.143	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
# I0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
# M1	0.0	0.0	0.0	0.0	0.0	1.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
# ...
# -------
# 	    A	B	C	D	E
# S	    0.0	0.0	0.0	0.0	0.0
# I0	0.0	0.0	0.0	0.0	0.0
# M1	0.167	0.0	0.0	0.0	0.833
# ...

        
# Example 3
# ----------
theta = 0.35
alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
alignment_strings = ["VKKLGEQFR-NKTIFNQPSGGDLEIVMHSFNCGGEFFYCNTTQLFN----------NSTES---DTITL",
                     "VKKLGEQFR-NKTIFNQPSGGDLEIVMHSFNCGGEFFYCNTQLFN----------NSTDNG---DTITL",
                     "VKKLGEQFR-NKTIFNQPSGGDLEIVMHSFNCGGEFFYCNTQLF----------NSTESNN--DTITL",
                     "VDKLREQFGKNKTIFNQPSGGDLEIVMHTFNCGGEFFYCNTTQLFNSTWNS---TGNGTESYNGQENGT",
                     "VDKLREQFGKNKTIFNQPSGGDLEIVMHTFNCGGEFFYCNTTQLFNSTWNG---TNTT--GLDG--NDTITL",
                     "VKLREQFGKNKTIFNQPSGGDLEIVTHFNCAGGEFFYCNTTQLFNSNWTC---NSTE--GLHG--DDTITL",
                     "VKLREQFG-GNKTIFNQPSGGGLEIVMHSFNCGGEFFYCNTTQLFNN--TR------NSTESNGQGNDTITL",
                     "VKLREQFGKNKTIFKQSGGDLEIVTHFNCAGEFFYCNTTQLFNSNWTE------NSITGLDG--NDTITL",
                     "VGKLREQFGK-NKTIFNQPSGGDLEIVMHSFNCQGEFFYCNTRLENSTWDNSTWNSTGKDKEGN--NDTITL"]
max_len = max(len(seq) for seq in alignment_strings)
padded_alignment = [seq.ljust(max_len, '-') for seq in alignment_strings]
alignment = np.array([list(seq) for seq in padded_alignment])

states, transition_matrix, emission_matrix = ConstructProfileHMM(theta, alphabet, alignment)
WriteTransitionAndEmission(alphabet, states, transition_matrix, emission_matrix)
# Output is saved as HIVProfileHMM.txt


# -----------------------------------------------
# Profile HMM with Pseudocounts
# -----------------------------------------------

import numpy as np

def ConstructProfileHMMWithPseudocounts(theta, sigma, alphabet, alignment):
    '''Construct a profile HMM with pseudocounts from a multiple alignment.

    Input: A threshold θ and a pseudocount σ, followed by an alphabet Σ, 
    followed by a multiple alignment Alignment whose strings are formed from Σ.
    Output: The transition and emission matrices of HMM(Alignment, θ, σ).'''

    alignment = np.array([list(row) for row in alignment])
    n, l = alignment.shape
    threshold = theta * n
    
    # Determine which columns to keep as match states
    kept = [True] * l
    for i in range(l):
        if sum(alignment[:, i] == '-') >= threshold:
            kept[i] = False
    
    # Level mapping (profile column index)
    levels = [0] * l
    for i in range(l):
        levels[i] = levels[i - 1] if i > 0 else 0
        if kept[i]:
            levels[i] += 1
    
    alphabetDict = {ch: idx for idx, ch in enumerate(alphabet)}
    
    def GetIndex(level, state):
        # States: 0 = I, 1 = M, 2 = D/E
        if level == 0:
            return 1 if state == 0 else 0
        if state != 0:
            return 3 * level - 2 + state
        else:
            return 3 * level + 1
    
    size = sum(kept) * 3 + 3
    transition = np.zeros((size, 3), dtype=float)  # I, M, D/E
    emission = np.zeros((size, len(alphabet)), dtype=float)
    
    for i in range(n):
        lastLevel = 0
        lastState = -1
        lastInd = GetIndex(lastLevel, lastState)
        for j in range(l):
            currLevel = levels[j]
            if kept[j]:
                currState = 2 if alignment[i, j] == '-' else 1
                currInd = GetIndex(currLevel, currState)
                transition[lastInd, currState] += 1
                if currState == 1:
                    emission[currInd, alphabetDict[alignment[i, j]]] += 1
                lastInd = currInd
            else:
                if alignment[i, j] != '-':
                    currState = 0
                    currInd = GetIndex(currLevel, currState)
                    transition[lastInd, currState] += 1
                    emission[currInd, alphabetDict[alignment[i, j]]] += 1
                    lastInd = currInd
        transition[lastInd, 2] += 1  # End with Delete/End
    
    # Normalize transition and emission counts to probabilities
    for i in range(transition.shape[0]):
        s = sum(transition[i, :])
        if s != 0:
            transition[i, :] /= s
        s = sum(emission[i, :])
        if s != 0:
            emission[i, :] /= s
    
    # Add pseudocounts and renormalize transitions
    for i in range(transition.shape[0] - 4):
        transition[i, :] += sigma
        transition[i, :] /= sum(transition[i, :])
    for i in range(-4, -1):
        transition[i, 0] += sigma
        transition[i, 2] += sigma
        transition[i, :] /= sum(transition[i, :])
    
    # Add pseudocounts and renormalize emissions for I and M states only
    for i in range(1, emission.shape[0] - 1):
        if i % 3 != 0:  # M or I states (skip D states)
            emission[i, :] += sigma
            emission[i, :] /= sum(emission[i, :])
    
    # Build full transition matrix from compact transition
    fullTransition = np.zeros((size, size), dtype=float)
    fullTransition[0, 1:4] = transition[0, :]
    
    for i in range(1, size - 4):
        mod = i % 3
        if mod == 1:  # M state
            fullTransition[i, i:i + 3] = transition[i, :]
        elif mod == 2:  # D state
            fullTransition[i, i + 2:i + 5] = transition[i, :]
        else:  # I state
            fullTransition[i, i + 1:i + 4] = transition[i, :]
    
    for i in range(-4, -1):
        fullTransition[i, -2] = transition[i, 0]
        fullTransition[i, -1] = transition[i, 2]
    
    # Build list of states names
    states = [''] * size
    states[0] = 'S'
    states[-1] = 'E'
    states[1] = 'I0'
    stateTypes = ('M', 'D', 'I')
    for i in range(2, size - 1):
        states[i] = stateTypes[(i + 1) % 3] + str((i + 1) // 3)
    
    return states, np.round(fullTransition, 3), np.round(emission, 3)

def ReadFromFilePseudocounts(filename):
    with open(filename, "r") as file:
        data = file.read().split()
    theta = float(data[0])
    sigma = float(data[1])
    ind = [i for i in range(len(data)) if data[i] == '--------']
    alphabet = data[ind[0] + 1:ind[1]]
    alignment = np.array([[*s] for s in data[ind[1] + 1:]])
    return theta, sigma, alphabet, alignment

def PrintTransitionAndEmission(alphabet, states, fullTransition, emission):
    print('\t' + '\t'.join(states))
    for i in range(fullTransition.shape[0]):
        print('\t'.join([states[i]] + list(map(str, fullTransition[i, :]))))
    print('--------')
    print('\t' + '\t'.join(alphabet))
    for i in range(emission.shape[0]):
        print('\t'.join([states[i]] + list(map(str, emission[i, :]))))

def WriteTransitionAndEmission(alphabet, states, fullTransition, emission):
    with open("output.txt", "w") as file:
        file.write('\t' + '\t'.join(states) + '\n')
        for i in range(fullTransition.shape[0]):
            file.write('\t'.join([states[i]] + list(map(str, fullTransition[i, :]))) + '\n')
        file.write('--------\n')
        file.write('\t' + '\t'.join(alphabet) + '\n')
        for i in range(emission.shape[0]):
            file.write('\t'.join([states[i]] + list(map(str, emission[i, :]))) + '\n')

# Example 1
# ----------
theta, sigma = 0.358, 0.01
alphabet = ["A", "B", "C", "D", "E"]
input_alignment = ["A-A", "ADA", "ACA", "A-C", "-EA", "D-A"]
alignment = np.array([list(row) for row in input_alignment])

states, fullTransition, emission = ConstructProfileHMMWithPseudocounts(theta, sigma, alphabet, alignment)
PrintTransitionAndEmission(alphabet, states, fullTransition, emission)

# Output:
    # 	    S	I0	M1	D1	I1	M2	D2	I2	E
    # S	    0.0	0.01	0.819	0.172	0.0	0.0	0.0	0.0	0.0
    # I0	0.0	0.333	0.333	0.333	0.0	0.0	0.0	0.0	0.0
    # M1	0.0	0.0	0.0	0.0	0.398	0.592	0.01	0.0	0.0
    # D1	0.0	0.0	0.0	0.0	0.981	0.01	0.01	0.0	0.0
    # I1	0.0	0.0	0.0	0.0	0.01	0.981	0.01	0.0	0.0
    # M2	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.01	0.99
    # D2	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.5	0.5
    # I2	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.5	0.5
    # E	    0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
    # --------
    # 	    A	B	C	D	E
    # S	    0.0	0.0	0.0	0.0	0.0
    # I0	0.2	0.2	0.2	0.2	0.2
    # M1	0.771	0.01	0.01	0.2	0.01
    # D1	0.0	0.0	0.0	0.0	0.0
    # I1	0.01	0.01	0.327	0.327	0.327
    # M2	0.803	0.01	0.168	0.01	0.01
    # D2	0.0	0.0	0.0	0.0	0.0
    # I2	0.2	0.2	0.2	0.2	0.2
    # E  	0.0	0.0	0.0	0.0	0.0


# Example 2
# ----------
theta, sigma, alphabet, alignment = ReadFromFilePseudocounts("dataset_30332_5.txt")
states, fullTransition, emission = ConstructProfileHMMWithPseudocounts(theta, sigma, alphabet, alignment)
WriteTransitionAndEmission(alphabet, states, fullTransition, emission)

# Output:
    # 	    S	I0	M1	D1	I1	M2	D2	I2	E
    # S	    0.0	0.01	0.786	0.204	0.0	0.0	0.0	0.0	0.0
    # I0	0.0	0.333	0.333	0.333	0.0	0.0	0.0	0.0	0.0
    # M1	0.0	0.0	0.0	0.0	0.495	0.252	0.252	0.0	0.0
    # D1	0.0	0.0	0.0	0.0	0.01	0.981	0.01	0.0	0.0
    # I1	0.0	0.0	0.0	0.0	0.01	0.981	0.01	0.0	0.0
    # M2	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.01	0.99
    # D2	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.01	0.99
    # I2	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.5	0.5
    # E	    0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
    # --------
    # 	    A	B	C	D	E
    # S	    0.0	0.0	0.0	0.0	0.0
    # I0	0.2	0.2	0.2	0.2	0.2
    # M1	0.01	0.01	0.01	0.962	0.01
    # D1	0.0	0.0	0.0	0.0	0.0
    # I1	0.01	0.01	0.962	0.01	0.01
    # M2	0.01	0.248	0.486	0.248	0.01
    # D2	0.0	0.0	0.0	0.0	0.0
    # I2	0.2	0.2	0.2	0.2	0.2
    # E	    0.0	0.0	0.0	0.0	0.0


# -----------------------------------------------
# Sequence Alignment with Profile HMM
# -----------------------------------------------

import numpy as np

def SequenceAlignmentProfileHMM(states, fullT, emit, alphabet, text):
    '''Align a new sequence to a family of sequences using a profile HMM via the Viterbi algorithm.
    
    The choice of alignment path is based on varying transition and emission probabilities. 
    These individual scoring parameters for each column allow us to detect subtle similarities 
    that traditional sequence alignment often misses. 
    
    Input: A string x followed by a threshold θ and a pseudocount σ, followed by an 
    alphabet Σ, followed by a multiple alignment whose strings are formed from Σ.
    Output: An optimal hidden path emitting x in HMM(Alignment, θ, σ).'''
    
    N = len(states)
    L = len(text)
    alpha_idx = {ch: idx for idx, ch in enumerate(alphabet)}
    
    # Use log-probs to avoid underflow
    logT = np.log(fullT + 1e-300)
    logE = np.log(emit + 1e-300)
    
    V = np.full((N, L + 1), -np.inf)
    B = np.zeros((N, L + 1), dtype=int)
    
    V[0, 0] = 0  # start
    
    for i in range(L + 1):
        for s in range(1, N):  # skip Start
            maxProb = -np.inf
            argMax = -1
            if emit[s].sum() == 0:  # D state → no emission
                for prev in range(N):
                    if fullT[prev, s] > 0 and V[prev, i] > -np.inf:
                        prob = V[prev, i] + logT[prev, s]
                        if prob > maxProb:
                            maxProb = prob
                            argMax = prev
                V[s, i] = maxProb
                B[s, i] = argMax
            elif i > 0:  # M or I state → consumes a character
                symbol = text[i - 1]
                for prev in range(N):
                    if fullT[prev, s] > 0 and V[prev, i - 1] > -np.inf:
                        prob = V[prev, i - 1] + logT[prev, s] + logE[s, alpha_idx[symbol]]
                        if prob > maxProb:
                            maxProb = prob
                            argMax = prev
                V[s, i] = maxProb
                B[s, i] = argMax

    # Termination
    lastProbs = V[:, L] + logT[:, -1]
    bestEnd = np.argmax(lastProbs)
    cur = bestEnd
    i = L
    path = []
    while cur != 0:
        path.append(states[cur])
        prev = B[cur, i]
        if emit[cur].sum() > 0:
            i -= 1  # Only move left when emitting
        cur = prev
    path.reverse()

    return path

def ConstructProfileHMMWithPseudocounts(theta, sigma, alphabet, alignment):
    alignment = np.array([list(row) for row in alignment])
    n, l = alignment.shape
    threshold = theta * n
    
    # Determine which columns to keep as match states
    kept = [True] * l
    for i in range(l):
        if sum(alignment[:, i] == '-') >= threshold:
            kept[i] = False
    
    # Level mapping (profile column index)
    levels = [0] * l
    for i in range(l):
        levels[i] = levels[i - 1] if i > 0 else 0
        if kept[i]:
            levels[i] += 1
    
    alphabetDict = {ch: idx for idx, ch in enumerate(alphabet)}
    
    def GetIndex(level, state):
        # States: 0 = I, 1 = M, 2 = D/E
        if level == 0:
            return 1 if state == 0 else 0
        if state != 0:
            return 3 * level - 2 + state
        else:
            return 3 * level + 1
    
    size = sum(kept) * 3 + 3
    transition = np.zeros((size, 3), dtype=float)  # I, M, D/E
    emission = np.zeros((size, len(alphabet)), dtype=float)
    
    for i in range(n):
        lastLevel = 0
        lastState = -1
        lastInd = GetIndex(lastLevel, lastState)
        for j in range(l):
            currLevel = levels[j]
            if kept[j]:
                currState = 2 if alignment[i, j] == '-' else 1
                currInd = GetIndex(currLevel, currState)
                transition[lastInd, currState] += 1
                if currState == 1:
                    emission[currInd, alphabetDict[alignment[i, j]]] += 1
                lastInd = currInd
            else:
                if alignment[i, j] != '-':
                    currState = 0
                    currInd = GetIndex(currLevel, currState)
                    transition[lastInd, currState] += 1
                    emission[currInd, alphabetDict[alignment[i, j]]] += 1
                    lastInd = currInd
        transition[lastInd, 2] += 1  # End with Delete/End
    
    # Normalize transition and emission counts to probabilities
    for i in range(transition.shape[0]):
        s = sum(transition[i, :])
        if s != 0:
            transition[i, :] /= s
        s = sum(emission[i, :])
        if s != 0:
            emission[i, :] /= s
    
    # Add pseudocounts and renormalize transitions
    for i in range(transition.shape[0] - 4):
        transition[i, :] += sigma
        transition[i, :] /= sum(transition[i, :])
    for i in range(-4, -1):
        transition[i, 0] += sigma
        transition[i, 2] += sigma
        transition[i, :] /= sum(transition[i, :])
    
    # Add pseudocounts and renormalize emissions for I and M states only
    for i in range(1, emission.shape[0] - 1):
        if i % 3 != 0:  # M or I states (skip D states)
            emission[i, :] += sigma
            emission[i, :] /= sum(emission[i, :])
    
    # Build full transition matrix from compact transition
    fullTransition = np.zeros((size, size), dtype=float)
    fullTransition[0, 1:4] = transition[0, :]
    
    for i in range(1, size - 4):
        mod = i % 3
        if mod == 1:  # M state
            fullTransition[i, i:i + 3] = transition[i, :]
        elif mod == 2:  # D state
            fullTransition[i, i + 2:i + 5] = transition[i, :]
        else:  # I state
            fullTransition[i, i + 1:i + 4] = transition[i, :]
    
    for i in range(-4, -1):
        fullTransition[i, -2] = transition[i, 0]
        fullTransition[i, -1] = transition[i, 2]
    
    # Build list of states names
    states = [''] * size
    states[0] = 'S'
    states[-1] = 'E'
    states[1] = 'I0'
    stateTypes = ('M', 'D', 'I')
    for i in range(2, size - 1):
        states[i] = stateTypes[(i + 1) % 3] + str((i + 1) // 3)
    
    return states, np.round(fullTransition, 3), np.round(emission, 3)

def ReadFromFileProfileHMMAlignment(filename):
    with open(filename, "r") as file:
        data = file.read().splitlines()
    data = [line.strip() for line in data if line.strip()]

    text = data[0]
    theta, sigma = map(float, data[2].split())
    alphabet = data[4].split()
    alignment = np.array([list(row.strip()) for row in data[6:]])
    
    return text, theta, sigma, alphabet, alignment



# Example 1
# ----------
text = "AEFDFDC"
theta, sigma = 0.4, 0.01
alphabet = ["A", "B", "C", "D", "E", "F"]
alignment = ["ACDEFACADF", "AFDA---CCF", "A--EFD-FDC", "ACAEF--A-C", "ADDEFAAADF"]

states, fullT, emit = ConstructProfileHMMWithPseudocounts(theta, sigma, alphabet, alignment)
path = SequenceAlignmentProfileHMM(states, fullT, emit, alphabet, text)
print(*path) # Output: M1 D2 D3 M4 M5 I5 M6 M7 M8


# Example 2
# ---------
text, theta, sigma, alphabet, alignment = ReadFromFileProfileHMMAlignment("dataset_30332_14.txt")
states, fullT, emit = ConstructProfileHMMWithPseudocounts(theta, sigma, alphabet, alignment)
path = SequenceAlignmentProfileHMM(states, fullT, emit, alphabet, text)
print(*path) # Output: M1 I1 M2 M3 M4 M5 D6 I6 I6 I6 I6 I6 M7 M8 M9 M10 M11 
    # I11 M12 M13 I13 I13 I13 I13 I13 M14 I14 M15 M16 I16 M17 I17 M18 M19 
    # M20 M21 M22 M23 I23 M24 M25 I25 M26 M27 M28 I28 M29 M30 I30 I30


# -----------------------------------------------
# Probability of Emission
# -----------------------------------------------

def ProbabilityEmission(alignment, match_position, amino_acid):
    '''Computes the emission probability of a symbol from a given match state in an alignment.

    Input: List of aligned sequences (strings of equal length), 1-based index of the match 
    state (i.e., column in the alignment), and the one-letter amino acid whose emission
    probability we want.
    Output: The emission probability.'''
    
    if not alignment or not alignment[0]:
        raise ValueError("Alignment must not be empty.")
    
    if match_position < 1 or match_position > len(alignment[0]):
        raise ValueError("match_position out of range.")

    column_index = match_position - 1
    residues = [seq[column_index] for seq in alignment if seq[column_index] != '-']
    
    if not residues:
        return 0.0

    amino_acid_count = residues.count(amino_acid)
    total_count = len(residues)
    return round(amino_acid_count / total_count, 3)


# What is the probability that serine (S) will be emitted from M(8), the 
# 8th match state of the profile HMM corresponding to this alignment? 
# -----------------------------------------------------------------------
alignment = ["M--QKCASHLE-AR",
             "MSNL-C-APD-LER",
             "MSAPNCARKYDI-R",
             "MS-SSCADED-IIR",
             "M--TKC-SKLEIDR"]

result = ProbabilityEmission(alignment, match_position=8, amino_acid='S')
print(result)  # Output: 0.4






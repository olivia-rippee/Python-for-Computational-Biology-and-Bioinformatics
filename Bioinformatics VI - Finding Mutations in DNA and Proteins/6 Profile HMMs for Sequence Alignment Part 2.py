import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics VI - Finding Mutations in DNA and Proteins/Data")

# -----------------------------------------------
# HMM Parameter Estimation
# -----------------------------------------------

import numpy as np
from collections import defaultdict

def HMMParameterEstimation(x, alphabet, path, states):
    '''Find optimal parameters explaining the emitted string and hidden path of an HMM.
    
    Input: A string x of symbols emitted from an HMM, followed by the HMM's alphabet Σ, 
    followed by a path π, followed by the collection of states of the HMM.
    Output: A transition matrix Transition followed by an emission matrix Emission that
    maximize Pr(x, π) over all possible transition and emission matrices.'''
    
    # Initialize count dictionaries
    transition_counts = defaultdict(lambda: defaultdict(int))
    emission_counts = defaultdict(lambda: defaultdict(int))

    # Count transitions
    for i in range(len(path) - 1):
        state_from = path[i]
        state_to = path[i + 1]
        transition_counts[state_from][state_to] += 1

    # Count emissions
    for i in range(len(x)):
        state = path[i]
        symbol = x[i]
        emission_counts[state][symbol] += 1

    # Convert counts to probabilities
    # Transition matrix
    transition_matrix = {}
    for state_from in states:
        total = sum(transition_counts[state_from][s] for s in states)
        transition_matrix[state_from] = {}
        for state_to in states:
            prob = (transition_counts[state_from][state_to] / total) if total > 0 else 1 / len(states)
            transition_matrix[state_from][state_to] = round(prob, 3)

    # Emission matrix
    emission_matrix = {}
    for state in states:
        total = sum(emission_counts[state][symbol] for symbol in alphabet)
        emission_matrix[state] = {}
        for symbol in alphabet:
            prob = (emission_counts[state][symbol] / total) if total > 0 else 1 / len(alphabet)
            emission_matrix[state][symbol] = round(prob, 3)

    return transition_matrix, emission_matrix

def ReadHMMEstimationInput(filename):
    with open(filename, "r") as file:
        lines = [line.strip() for line in file if line.strip()]
    
    x = lines[0]
    alphabet = lines[2].split()
    path = lines[4]
    states = lines[6].split()
    
    return x, alphabet, path, states

def PrintHMMParameters(transition_matrix, emission_matrix, states, alphabet):
    # Print transition matrix
    print('\t' + '\t'.join(states))
    for state_from in states:
        row = [state_from] + [str(transition_matrix[state_from][state_to]) for state_to in states]
        print('\t'.join(row))

    print('--------')
    
    # Print emission matrix
    print('\t' + '\t'.join(alphabet))
    for state in states:
        row = [state] + [str(emission_matrix[state][symbol]) for symbol in alphabet]
        print('\t'.join(row))

def WriteHMMParameters(transition_matrix, emission_matrix, states, alphabet):
    with open("output.txt", "w") as file:
        # Write transition matrix header
        file.write('\t' + '\t'.join(states) + '\n')
        for from_state in states:
            row = [from_state] + [f"{transition_matrix[from_state][to_state]:.3f}" for to_state in states]
            file.write('\t'.join(row) + '\n')

        file.write('--------\n')

        # Write emission matrix header
        file.write('\t' + '\t'.join(alphabet) + '\n')
        for state in states:
            row = [state] + [f"{emission_matrix[state][symbol]:.3f}" for symbol in alphabet]
            file.write('\t'.join(row) + '\n')


# Example 1
# ----------
x = "yzzzyxzxxx"
alphabet = ["x", "y", "z"]
path = "BBABABABAB"
states = ["A", "B", "C"]

transition, emission = HMMParameterEstimation(x, alphabet, path, states)
PrintHMMParameters(transition, emission, states, alphabet)
# Output: 
    #  	A	B	C
    # A	0.0	1.0	0.0
    # B	0.8	0.2	0.0
    # C	0.333	0.333	0.333
    # --------
    # 	x	y	z
    # A	0.25	0.25	0.5
    # B	0.5	0.167	0.333
    # C	0.333	0.333	0.333


# Example 2
# ----------
x, alphabet, path, states = ReadHMMEstimationInput("dataset_30333_4.txt")
transition, emission = HMMParameterEstimation(x, alphabet, path, states)
WriteHMMParameters(transition, emission, states, alphabet)
# Output:
    # 	A	B	C	D
    # A	0.333	0.190	0.190	0.286
    # B	0.115	0.346	0.269	0.269
    # C	0.208	0.167	0.292	0.333
    # D	0.179	0.321	0.214	0.286
    # --------
    # 	x	y	z
    # A	0.381	0.238	0.381
    # B	0.462	0.192	0.346
    # C	0.333	0.417	0.250
    # D	0.448	0.379	0.172


# -----------------------------------------------
# Viterbi Learning for HMM Parameter Estimation
# -----------------------------------------------

import numpy as np
from collections import defaultdict

def ViterbiLearning(j, x, alphabet, states, transition, emission):
    '''Estimate the parameters of an HMM explaining an emitted string using Viterbi learning.
    
    Viterbi learning is dependent on the initial guess for Parameters, so it may become 
    stuck in a local optimum. Like other heuristics, it is often run many times, retaining the 
    best choice of Parameters.
    
    Input: A number of iterations j, followed by a string x of symbols emitted by an HMM, 
    followed by the HMM's alphabet Σ, followed by the HMM's states, followed by initial 
    transition and emission matrices for the HMM.
    Output: Emission and transition matrices resulting from applying Viterbi learning for j iterations.'''
    
    for _ in range(j):
        path = Viterbi(x, alphabet, states, transition, emission)
        transition, emission = HMMParameterEstimation(x, path, states, alphabet)  # fixed argument order
    return transition, emission

def Viterbi(x, alphabet, states, transition, emission):
    n = len(x)
    k = len(states)
    dp = np.full((k, n), -np.inf)
    backtrack = np.zeros((k, n), dtype=int)
    index_state = {i: s for i, s in enumerate(states)}

    # Initialization
    for i, s in enumerate(states):
        dp[i][0] = np.log(1e-10 + emission[s].get(x[0], 0))  # assume uniform initial prob

    # Dynamic programming
    for j in range(1, n):
        for curr in range(k):
            curr_state = states[curr]
            for prev in range(k):
                prev_state = states[prev]
                trans_prob = transition[prev_state].get(curr_state, 0)
                emis_prob = emission[curr_state].get(x[j], 0)
                prob = dp[prev][j - 1] + np.log(1e-10 + trans_prob) + np.log(1e-10 + emis_prob)
                if prob > dp[curr][j]:
                    dp[curr][j] = prob
                    backtrack[curr][j] = prev

    # Termination
    last = np.argmax(dp[:, -1])
    path = [last]
    for j in range(n - 1, 0, -1):
        last = backtrack[last][j]
        path.append(last)
    path = path[::-1]

    # Debugging length check
    # print(f"Viterbi path length: {len(path)}, emitted string length: {len(x)}") 

    return [index_state[i] for i in path]

def HMMParameterEstimation(x, path, states, alphabet):
    # Initialize count dictionaries
    transition_counts = defaultdict(lambda: defaultdict(int))
    emission_counts = defaultdict(lambda: defaultdict(int))

    # Count transitions
    for i in range(len(path) - 1):
        state_from = path[i]
        state_to = path[i + 1]
        transition_counts[state_from][state_to] += 1

    # Count emissions
    for i in range(len(x)):
        state = path[i]
        symbol = x[i]
        emission_counts[state][symbol] += 1

    # Convert counts to probabilities
    # Transition matrix
    transition_matrix = {}
    for state_from in states:
        total = sum(transition_counts[state_from][s] for s in states)
        transition_matrix[state_from] = {}
        for state_to in states:
            prob = (transition_counts[state_from][state_to] / total) if total > 0 else 1 / len(states)
            transition_matrix[state_from][state_to] = round(prob, 3)

    # Emission matrix
    emission_matrix = {}
    for state in states:
        total = sum(emission_counts[state][symbol] for symbol in alphabet)
        emission_matrix[state] = {}
        for symbol in alphabet:
            prob = (emission_counts[state][symbol] / total) if total > 0 else 1 / len(alphabet)
            emission_matrix[state][symbol] = round(prob, 3)

    return transition_matrix, emission_matrix

def ReadViterbiLearningInput(filename):
    with open(filename, "r") as file:
        lines = [line.strip() for line in file if line.strip()]

    # Find all divider indices
    divider_indices = [i for i, line in enumerate(lines) if line == '--------']

    j = int(lines[0])
    x = lines[divider_indices[0] + 1]         # line after first '--------' is x
    alphabet = lines[divider_indices[1] + 1].split()
    states = lines[divider_indices[2] + 1].split()

    transition_lines = lines[divider_indices[3] + 1: divider_indices[4]]
    emission_lines = lines[divider_indices[4] + 1:]

    # Parse transition matrix
    headers = transition_lines[0].split()
    transition = {}
    for line in transition_lines[1:]:
        parts = line.split()
        row_state = parts[0]
        transition[row_state] = {}
        for i, val in enumerate(parts[1:]):
            transition[row_state][headers[i]] = float(val)

    # Parse emission matrix
    headers = emission_lines[0].split()
    emission = {}
    for line in emission_lines[1:]:
        parts = line.split()
        row_state = parts[0]
        emission[row_state] = {}
        for i, val in enumerate(parts[1:]):
            emission[row_state][headers[i]] = float(val)

    return j, x, alphabet, states, transition, emission

def PrintHMMParameters(transition_matrix, emission_matrix, states, alphabet):
    # Print transition matrix
    print('\t' + '\t'.join(states))
    for state_from in states:
        row = [state_from] + [str(transition_matrix[state_from][state_to]) for state_to in states]
        print('\t'.join(row))

    print('--------')
    
    # Print emission matrix
    print('\t' + '\t'.join(alphabet))
    for state in states:
        row = [state] + [str(emission_matrix[state][symbol]) for symbol in alphabet]
        print('\t'.join(row))

def WriteHMMParameters(transition_matrix, emission_matrix, states, alphabet):
    with open("output.txt", "w") as file:
        # Write transition matrix header
        file.write('\t' + '\t'.join(states) + '\n')
        for from_state in states:
            row = [from_state] + [f"{transition_matrix[from_state][to_state]:.3f}" for to_state in states]
            file.write('\t'.join(row) + '\n')

        file.write('--------\n')

        # Write emission matrix header
        file.write('\t' + '\t'.join(alphabet) + '\n')
        for state in states:
            row = [state] + [f"{emission_matrix[state][symbol]:.3f}" for symbol in alphabet]
            file.write('\t'.join(row) + '\n')


# Example 1
# ----------
j = 100
x = "zyzxzxxxzz"
alphabet = ['x', 'y', 'z']
states = ['A', 'B']
transition = {'A': {'A': 0.599, 'B': 0.401},
              'B': {'A': 0.294, 'B': 0.706}}
emission = {'A': {'x': 0.424, 'y': 0.367, 'z': 0.209},
            'B': {'x': 0.262, 'y': 0.449, 'z': 0.289}}

transition_estimate, emission_estimate = ViterbiLearning(j, x, alphabet, states, transition, emission)
PrintHMMParameters(transition_estimate, emission_estimate, states, alphabet)

# Output:
    # 	A	B
    # A	0.5	0.5
    # B	0.0	1.0
    # --------
    # 	x	y	z
    # A	0.333	0.333	0.333
    # B	0.4	0.1	0.5


# Example 2
# ----------
j, x, alphabet, states, transition, emission = ReadViterbiLearningInput("dataset_30333_8.txt")
transition_estimate, emission_estimate = ViterbiLearning(j, x, alphabet, states, transition, emission)
WriteHMMParameters(transition_estimate, emission_estimate, states, alphabet)

# Output:
    # 	A	B	C	D
    # A	0.364	0.636	0.000	0.000
    # B	0.000	0.694	0.000	0.306
    # C	0.467	0.533	0.000	0.000
    # D	0.000	0.000	0.625	0.375
    # --------
    # 	x	y	z
    # A	1.000	0.000	0.000
    # B	0.000	0.440	0.560
    # C	0.000	0.533	0.467
    # D	1.000	0.000	0.000


# -----------------------------------------------
# Soft Decoding
# -----------------------------------------------

def SoftDecoding(x, alphabet, states, transition, emission):
    '''Find the probability that an HMM was in a particular state at a particular 
    moment given its emitted string.
    
    Input: A string x, followed by the alphabet Σ from which x was constructed, followed by 
    the states States, transition matrix Transition, and emission matrix Emission of an HMM 
    (Σ, States, Transition, Emission).
    Output: An |x| x |States| matrix whose (i, k)-th element holds the conditional probability 
    Pr(πi = k|x).'''
    
    n = len(x)
    
    # Forward matrix
    forward = [{state: 0.0 for state in states} for _ in range(n)]
    for state in states:
        forward[0][state] = (1 / len(states)) * emission[state][x[0]]
    for i in range(1, n):
        for curr_state in states:
            forward[i][curr_state] = sum(
                forward[i - 1][prev_state] * transition[prev_state][curr_state]
                for prev_state in states
            ) * emission[curr_state][x[i]]

    # Backward matrix
    backward = [{state: 0.0 for state in states} for _ in range(n)]
    for state in states:
        backward[-1][state] = 1.0
    for i in range(n - 2, -1, -1):
        for curr_state in states:
            backward[i][curr_state] = sum(
                transition[curr_state][next_state] * emission[next_state][x[i + 1]] * backward[i + 1][next_state]
                for next_state in states
            )

    # Total probability
    total_prob = sum(forward[-1][state] for state in states)

    # Soft decoding result
    result = []
    for i in range(n):
        row = [forward[i][state] * backward[i][state] / total_prob for state in states]
        result.append(row)

    return result

def ParseSoftDecodingInputFromFile(filename):
    with open(filename) as file:
        lines = [line.strip() for line in file if line.strip()]

    # Parse sections based on "--------"
    sections = []
    current = []
    for line in lines:
        if line.startswith("--------"):
            if current:
                sections.append(current)
                current = []
        else:
            current.append(line)
    if current:
        sections.append(current)

    # Unpack sections
    x = sections[0][0]
    alphabet = sections[1][0].split()
    states = sections[2][0].split()

    # Parse transition matrix
    transition = {}
    header_states = sections[3][0].split()
    for row in sections[3][1:]:
        parts = row.split()
        state = parts[0]
        transition[state] = {
            header_states[i]: float(parts[i + 1]) for i in range(len(header_states))}

    # Parse emission matrix
    emission = {}
    header_symbols = sections[4][0].split()
    for row in sections[4][1:]:
        parts = row.split()
        state = parts[0]
        emission[state] = {
            header_symbols[i]: float(parts[i + 1]) for i in range(len(header_symbols))}

    return x, alphabet, states, transition, emission


# Example 1
# ----------
x = "zyxxxxyxzz"
alphabet = ["x", "y", "z"]
states = ["A", "B"]
transition = {"A": {"A": 0.911, "B": 0.089},
              "B": {"A": 0.228, "B": 0.772}}
emission = {"A": {"x": 0.356, "y": 0.191, "z": 0.453},
            "B": {"x": 0.040, "y": 0.467, "z": 0.493}}

result = SoftDecoding(x, alphabet, states, transition, emission)
print("A\tB")
for row in result:
    print('\t'.join(f"{val:.4f}" for val in row))

# Output:
    # A	        B
    # 0.5438	0.4562
    # 0.6492	0.3508
    # 0.9647	0.0353
    # 0.9936	0.0064
    # 0.9957	0.0043
    # 0.9891	0.0109
    # 0.9154	0.0846
    # 0.9640	0.0360
    # 0.8737	0.1263
    # 0.8167	0.1833



# Example 2
# ----------
x, alphabet, states, transition, emission = ParseSoftDecodingInputFromFile("dataset_30334_5.txt")
result = SoftDecoding(x, alphabet, states, transition, emission)

with open("output.txt", "w") as file:
    file.write('\t'.join(states) + '\n')
    for row in result:
        file.write('\t'.join(f"{val:.4f}" for val in row) + '\n')

# Output:
    # A	        B
    # 0.1352	0.8648
    # 0.0430	0.9570
    # 0.3292	0.6708
    # 0.3072	0.6928
    # 0.2611	0.7389
    # 0.3339	0.6661
    # 0.0405	0.9595
    # 0.2890	0.7110
    # 0.2677	0.7323
    # 0.2729	0.7271


# -----------------------------------------------
# Baum-Welch Learning
# -----------------------------------------------

def BaumWelchLearning(x, alphabet, states, transition, emission, iterations):
    '''Performs the Baum-Welch algorithm (Expectation-Maximization) to estimate the 
    parameters of a Hidden Markov Model (HMM).

    This algorithm alternates between:
    - E-step: Computing the expected probabilities of being in each state at each time (γ),
              and the expected probabilities of transitioning between states at each step (ξ),
              using the forward-backward algorithm.
    - M-step: Updating the transition and emission probabilities based on these expectations.
    
    Input: A sequence of emitted symbols x = x1 . . . xn in an alphabet A, generated 
    by a k-state HMM with unknown transition and emission probabilities, initial 
    Transition and Emission matrices, and a number of iterations j.
    Output: A matrix of transition probabilities Transition and a matrix of emission 
    probabilities Emission that maximizes Pr(x, π) over all possible transition and 
    emission matrices and over all hidden paths π.'''
    
    n = len(x)
    A = alphabet
    S = states

    for _ in range(iterations):
        # Forward
        fwd = [{s: 0.0 for s in S} for _ in range(n)]
        for s in S:
            fwd[0][s] = 1 / len(S) * emission[s][x[0]]
        for i in range(1, n):
            for s in S:
                fwd[i][s] = sum(fwd[i - 1][s_prev] * transition[s_prev][s] for s_prev in S) * emission[s][x[i]]

        # Backward
        bwd = [{s: 0.0 for s in S} for _ in range(n)]
        for s in S:
            bwd[-1][s] = 1.0
        for i in range(n - 2, -1, -1):
            for s in S:
                bwd[i][s] = sum(transition[s][s_next] * emission[s_next][x[i + 1]] * bwd[i + 1][s_next] for s_next in S)

        # Total probability of the string
        prob_x = sum(fwd[n - 1][s] for s in S)

        # γ[i][s] = Pr(π_i = s | x)
        gamma = [{s: fwd[i][s] * bwd[i][s] / prob_x for s in S} for i in range(n)]

        # ξ[i][s1][s2] = Pr(π_i = s1 and π_{i+1} = s2 | x)
        xi = [{s1: {s2: (fwd[i][s1] * transition[s1][s2] * emission[s2][x[i + 1]] * bwd[i + 1][s2]) / prob_x
                    for s2 in S}
                for s1 in S}
            for i in range(n - 1)]

        # Update transition matrix
        new_transition = {s1: {} for s1 in S}
        for s1 in S:
            denom = sum(gamma[i][s1] for i in range(n - 1))
            for s2 in S:
                numer = sum(xi[i][s1][s2] for i in range(n - 1))
                new_transition[s1][s2] = numer / denom if denom > 0 else 0.0

        # Update emission matrix
        new_emission = {s: {a: 0.0 for a in A} for s in S}
        for s in S:
            denom = sum(gamma[i][s] for i in range(n))
            for a in A:
                numer = sum(gamma[i][s] for i in range(n) if x[i] == a)
                new_emission[s][a] = numer / denom if denom > 0 else 0.0

        # Assign updated parameters
        transition = new_transition
        emission = new_emission

    return transition, emission

def ParseBaumWelchInput(filename):
    with open(filename) as file:
        lines = [line.strip() for line in file if line.strip()]

    iterations = int(lines[0])
    sections = []
    current = []
    for line in lines[1:]:
        if line.startswith("--------"):
            if current:
                sections.append(current)
                current = []
        else:
            current.append(line)
    if current:
        sections.append(current)

    x = sections[0][0]
    alphabet = sections[1][0].split()
    states = sections[2][0].split()

    # Transition matrix
    transition = {}
    header_states = sections[3][0].split()
    for row in sections[3][1:]:
        parts = row.split()
        transition[parts[0]] = {header_states[i]: float(parts[i + 1]) for i in range(len(header_states))}

    # Emission matrix
    emission = {}
    header_symbols = sections[4][0].split()
    for row in sections[4][1:]:
        parts = row.split()
        emission[parts[0]] = {header_symbols[i]: float(parts[i + 1]) for i in range(len(header_symbols))}

    return iterations, x, alphabet, states, transition, emission

def PrintBaumWelch(transition, emission, states, alphabet):
    print('\t' + '\t'.join(states))
    for s1 in states:
        print(s1 + '\t' + '\t'.join(f"{transition[s1][s2]:.3f}" for s2 in states))
    print('--------')
    print('\t' + '\t'.join(alphabet))
    for s in states:
        print(s + '\t' + '\t'.join(f"{emission[s][a]:.3f}" for a in alphabet))

def WriteBaumWelch(transition, emission, states, alphabet):
    with open("output.txt", "w") as file:
        # Write transition matrix
        file.write('\t' + '\t'.join(states) + '\n')
        for s1 in states:
            file.write(s1 + '\t' + '\t'.join(f"{transition[s1][s2]:.3f}" for s2 in states) + '\n')
        file.write('--------\n')

        # Write emission matrix
        file.write('\t' + '\t'.join(alphabet) + '\n')
        for s in states:
            file.write(s + '\t' + '\t'.join(f"{emission[s][a]:.3f}" for a in alphabet) + '\n')


# Example 1
# ----------
iterations = 10
x = "xzyyzyzyxy"
alphabet = ["x", "y", "z"]
states = ["A", "B"]
transition = {"A": {"A": 0.019, "B": 0.981},
              "B": {"A": 0.668, "B": 0.332}}
emission = {"A": {"x": 0.175, "y": 0.003, "z": 0.821},
            "B": {"x": 0.196, "y": 0.512, "z": 0.293}}

transition_estimate, emission_estimate = BaumWelchLearning(x, alphabet, states, transition, emission, iterations)
PrintBaumWelch(transition_estimate, emission_estimate, states, alphabet)

# Output:
    # 	A	B
    # A	0.000	1.000
    # B	0.786	0.214
    # --------
    # 	x	y	z
    # A	0.242	0.000	0.758
    # B	0.172	0.828	0.000


# Example 2
# ----------
iterations, x, alphabet, states, transition, emission = ParseBaumWelchInput("dataset_30335_5.txt")
transition_estimate, emission_estimate = BaumWelchLearning(x, alphabet, states, transition, emission, iterations)
WriteBaumWelch(transition_estimate, emission_estimate, states, alphabet)

# Output:
    # 	A	B	C	D
    # A	0.000	0.089	0.910	0.001
    # B	0.536	0.000	0.000	0.464
    # C	0.000	0.916	0.084	0.000
    # D	0.233	0.415	0.352	0.000
    # --------
    # 	x	y	z
    # A	0.580	0.239	0.182
    # B	0.596	0.366	0.037
    # C	0.000	0.118	0.882
    # D	0.113	0.887	0.000

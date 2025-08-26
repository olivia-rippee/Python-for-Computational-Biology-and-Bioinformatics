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

def SaveHMMParameters(transition_matrix, emission_matrix, states, alphabet):
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
SaveHMMParameters(transition, emission, states, alphabet)
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

def SaveHMMParameters(transition_matrix, emission_matrix, states, alphabet):
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
SaveHMMParameters(transition_estimate, emission_estimate, states, alphabet)

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

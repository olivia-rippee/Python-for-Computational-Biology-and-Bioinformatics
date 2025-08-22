import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics VI - Finding Mutations in DNA and Proteins/Data")

# -----------------------------------------------
# Probability of a Hidden Path
# -----------------------------------------------

def ProbabilityHiddenPath(HiddenPath, States, TransitionMatrix):
    '''Compute the probability of a hidden path.
    
    Input: A hidden path π followed by the states States and transition matrix 
    Transition of an HMM (Σ, States, Transition, Emission).
    Output: The probability of this path, Pr(π).'''
    
    num_states = len(States)
    initial_prob = 1 / num_states  # Equal probability to start in any state
    probability = initial_prob

    for i in range(1, len(HiddenPath)):
        prev_state = HiddenPath[i - 1]
        curr_state = HiddenPath[i]
        probability *= TransitionMatrix[prev_state][curr_state]

    return probability

def ParseInputHiddenPath(source, from_file=False):
    if from_file:
        with open(source, "r") as file:
            lines = file.read().strip().splitlines()
    else:
        lines = source.strip().splitlines()

    divider_indices = [i for i, line in enumerate(lines) if '--------' in line]

    hidden_path = lines[0].strip()
    states = lines[divider_indices[0] + 1].strip().split()

    header = lines[divider_indices[1] + 1].strip().split()
    matrix_lines = lines[divider_indices[1] + 2:]

    transition_matrix = {}
    for line in matrix_lines:
        parts = line.strip().split()
        row_state = parts[0]
        values = list(map(float, parts[1:]))
        transition_matrix[row_state] = dict(zip(header, values))

    return hidden_path, states, transition_matrix



# Example 1
# -----------
input_string = """ABABBBAAAA
--------
A B
--------
\tA\tB
A\t0.377\t0.623
B\t0.26\t0.74"""

path, states, transition = ParseInputHiddenPath(input_string)
probability = ProbabilityHiddenPath(path, states, transition)
print(probability) # Output: 0.0003849286917546758


# Example 2
# -----------
path, states, transition = ParseInputHiddenPath("dataset_30326_8.txt", from_file=True)
probability = ProbabilityHiddenPath(path, states, transition)
print(probability) # Output: 7.324774418398683e-18


# -----------------------------------------------
# Probability of an Outcome Given a Hidden Path
# -----------------------------------------------

def ProbabilityOutcomeGivenPath(x, path, emission_matrix):
    '''Compute the probability that an HMM will emit a string given its hidden path.
    
    Input: A string x, followed by the alphabet from which x was constructed, followed 
    by a hidden path π, followed by the states States and emission matrix Emission of 
    an HMM (Σ, States, Transition, Emission).
    Output: The conditional probability Pr(x|π) that x will be emitted given that the 
    HMM follows the hidden path π.'''
    
    prob = 1.0
    for symbol, state in zip(x, path):
        prob *= emission_matrix[state][symbol]
    return prob

def ParseInputEmission(source, from_file=False):
    if from_file:
        with open(source, "r") as file:
            lines = file.read().strip().splitlines()
    else:
        lines = source.strip().splitlines()

    divider_indices = [i for i, line in enumerate(lines) if '--------' in line]

    x = lines[0].strip()
    # alphabet = lines[divider_indices[0] + 1].strip().split()
    hidden_path = lines[divider_indices[1] + 1].strip()
    # states = lines[divider_indices[2] + 1].strip().split()

    header = lines[divider_indices[3] + 1].strip().split()
    matrix_lines = lines[divider_indices[3] + 2:]

    emission_matrix = {}
    for line in matrix_lines:
        parts = line.strip().split()
        state = parts[0]
        probs = list(map(float, parts[1:]))
        emission_matrix[state] = dict(zip(header, probs))

    # No need to return alphabet and states; the information is contained in x and emission_matrix
    return x, hidden_path, emission_matrix 


# Example 1
# -----------
input_str = """zzzyxyyzzx
--------
x y z
--------
BAAAAAAAAA
--------
A B
--------
\tx\ty\tz
A\t0.176\t0.596\t0.228
B\t0.225\t0.572\t0.203"""

x, path, emission_matrix = ParseInputEmission(input_str)
probability = ProbabilityOutcomeGivenPath(x, path, emission_matrix)
print(probability) # Output: 3.5974895474624624e-06


# Example 2
# -----------
x, path, emission_matrix = ParseInputEmission("dataset_30326_10.txt", from_file=True)
probability = ProbabilityOutcomeGivenPath(x, path, emission_matrix)
print(probability) # Output: 3.665653574261883e-30



# -----------------------------------------------
# Find the Most Likely Hidden Path (Viterbi)
# -----------------------------------------------

def FindMostLikelyHiddenPath(observations, states, initial_prob, transition_prob, emission_prob):
    '''Find the maximum product weight path.
    
    Input: Observations x, list of states, and dictionaries of initial probabilities,
    transition probabilities, and emission probabilities.
    Output: String of state symbols with length equal to the legnth of x.'''
    
    n = len(observations)
    
    # Initialize Viterbi table and backpointer
    viterbi = {state: [0.0] * n for state in states}
    backpointer = {state: [None] * n for state in states}

    # Initialization (t = 0)
    for state in states:
        viterbi[state][0] = initial_prob[state] * emission_prob[state][observations[0]]

    # Recursion
    for t in range(1, n):
        for curr_state in states:
            max_prob = 0.0
            best_prev_state = None
            for prev_state in states:
                prob = (
                    viterbi[prev_state][t - 1]
                    * transition_prob[prev_state][curr_state]
                    * emission_prob[curr_state][observations[t]]
                )
                if prob > max_prob:
                    max_prob = prob
                    best_prev_state = prev_state
            viterbi[curr_state][t] = max_prob
            backpointer[curr_state][t] = best_prev_state

    # Termination
    final_state = max(states, key=lambda s: viterbi[s][n - 1])
    path = [final_state]

    # Backtracking
    for t in range(n - 1, 0, -1):
        final_state = backpointer[final_state][t]
        path.insert(0, final_state)

    return ''.join(path)


# Example
# ---------
observations = "HHTT"
states = ['F', 'B']
initial_prob = {'F': 0.5, 'B': 0.5}
transition_prob = {'F': {'F': 0.9, 'B': 0.1},
                   'B': {'F': 0.1, 'B': 0.9}}
emission_prob = {'F': {'H': 0.5, 'T': 0.5},
                 'B': {'H': 0.75, 'T': 0.25}}

path = FindMostLikelyHiddenPath(observations, states, initial_prob, transition_prob, emission_prob)
print(path)  # Output: FFFF



# -----------------------------------------------
# Decoding: Finding an Optimal Hidden Path
# -----------------------------------------------

import math

def ViterbiAlgorithm(x, states, transition, emission):
    '''Find an optimal hidden path in an HMM given a string of its emitted symbols using
    the Viterbi algorithm.
    
    Input: A string x, followed by the states States, transition matrix Transition, and 
    emission matrix Emission of an HMM (Σ, States, Transition, Emission).
    Output: A path that maximizes the (unconditional) probability Pr(x, π) over all possible paths π.'''


    n = len(x)
    # Initialize DP tables
    # Use log probabilities to avoid underflow
    V = {state: [-math.inf] * n for state in states}  # log prob table
    backpointer = {state: [None] * n for state in states}
    
    # Initial probabilities - uniform as per problem statement
    initial_prob = math.log(1 / len(states))
    
    # Initialization
    for state in states:
        emit_prob = emission[state].get(x[0], 1e-10)
        V[state][0] = initial_prob + math.log(emit_prob)
        backpointer[state][0] = None  # no previous at first position

    
    # Recursion
    for i in range(1, n):
        for curr_state in states:
            max_prob = -math.inf
            max_state = None
            emit_prob = emission[curr_state].get(x[i], 1e-10)
            log_emit_prob = math.log(emit_prob)
    
            for prev_state in states:
                trans_prob = transition[prev_state].get(curr_state, 0)
                if trans_prob == 0:
                    continue
                log_trans_prob = math.log(trans_prob)
    
                prob = V[prev_state][i - 1] + log_trans_prob + log_emit_prob
    
                if prob > max_prob:
                    max_prob = prob
                    max_state = prev_state
    
            V[curr_state][i] = max_prob
            backpointer[curr_state][i] = max_state
    
    # Termination - find max probability at last position
    max_final_prob = -math.inf
    max_final_state = None
    for state in states:
        if V[state][n - 1] > max_final_prob:
            max_final_prob = V[state][n - 1]
            max_final_state = state
    
    # Backtrack path
    path = [max_final_state]
    for i in range(n - 1, 0, -1):
        path.insert(0, backpointer[path[0]][i])
    
    return ''.join(path)

def ParseInputViterbi(source, from_file=False):
    if from_file:
        with open(source, "r") as file:
            lines = file.read().strip().splitlines()
    else:
        lines = source.strip().splitlines()

    divider_indices = [i for i, line in enumerate(lines) if '--------' in line]
    x = lines[0].strip()
    # alphabet = lines[divider_indices[0] + 1].strip().split()
    states = lines[divider_indices[1] + 1].strip().split()

    transition_header = lines[divider_indices[2] + 1].strip().split()
    transition_lines = lines[divider_indices[2] + 2: divider_indices[3]]
    transition_matrix = {}
    for line in transition_lines:
        parts = line.strip().split()
        state_from = parts[0]
        probs = list(map(float, parts[1:]))
        transition_matrix[state_from] = dict(zip(transition_header, probs))

    emission_header = lines[divider_indices[3] + 1].strip().split()
    emission_lines = lines[divider_indices[3] + 2:]
    emission_matrix = {}
    for line in emission_lines:
        parts = line.strip().split()
        state = parts[0]
        probs = list(map(float, parts[1:]))
        emission_matrix[state] = dict(zip(emission_header, probs))

    return x, states, transition_matrix, emission_matrix


# Example 1
# ----------
input_data = """xyxzzxyxyy
--------
x y z
--------
A B
--------
    A   B
A   0.641   0.359
B   0.729   0.271
--------
    x   y   z
A   0.117   0.691   0.192
B   0.097   0.42    0.483"""

x, states, transition, emission = ParseInputViterbi(input_data)
path = ViterbiAlgorithm(x, states, transition, emission)
print(path) # Output: AAABBAAAAA


# Example 2
# ----------
x, states, transition, emission = ParseInputViterbi("dataset_30327_7.txt", from_file=True)
path = ViterbiAlgorithm(x, states, transition, emission)
print(path) # Output: ABABAABAAABABBBAAABABBAAAABBABBBBABABABABBAAAABBABBABAABAABAAABAAAAAAAAAAAABAAAAAAABAAAABAAAAAAAABBA



# Example 3 - Apply to CG Islands on X Chromosome
# --------------------------------------------------
def ReadFASTA(filename, max_length=1_000_000):
    seq = []
    with open(filename, "r") as file:
        for line in file:
            if line.startswith('>'):
                continue
            seq.append(line.strip().upper())
            if sum(len(s) for s in seq) >= max_length:
                break
    return ''.join(seq)[:max_length]

def CountCGIslands(path):
    count = 0
    in_island = False
    for state in path:
        if state == 'CG' and not in_island:
            count += 1
            in_island = True
        elif state == 'NCG':
            in_island = False
    return count

# Input
# -----
sequence = ReadFASTA("chrX.txt")
states = ['CG', 'NCG']
alphabet = ['A', 'C', 'G', 'T']
transition = {'CG': {'CG': 0.999, 'NCG': 0.001},
              'NCG': {'CG': 0.0001, 'NCG': 0.9999}}
emission = {'CG': {'A': 0.15, 'C': 0.35, 'G': 0.35, 'T': 0.15},
            'NCG': {'A': 0.30, 'C': 0.20, 'G': 0.20, 'T': 0.30}}

# path = ViterbiAlgorithm(sequence, states, transition, emission) # long run
# num_islands = CountCGIslands(path)
# print(num_islands) # Expected Output: 62


# -----------------------------------------------
# Outcome Likelihood
# -----------------------------------------------
    
def LikelihoodOutcome(x, states, transition, emission):
    '''Find the probability that an HMM emits a given string.
    
    Input: A string x, followed by the states States, transition matrix Transition, 
    and emission matrix Emission of an HMM (Σ, States, Transition, Emission).
    Output: The probability Pr(x) that the HMM emits x.'''

    n = len(x)
    # Forward probabilities table: forward[state][i] = Pr(prefix of length i emitted and ends in state)
    forward = {state: [0.0] * n for state in states}

    # Initial probabilities: assume uniform start
    initial_prob = 1 / len(states)

    # Initialization: probability of starting in state and emitting x[0]
    for state in states:
        emit_prob = emission[state].get(x[0], 1e-10)
        forward[state][0] = initial_prob * emit_prob

    # Recursion
    for i in range(1, n):
        for curr_state in states:
            total_prob = 0.0
            emit_prob = emission[curr_state].get(x[i], 1e-10)
            if emit_prob == 0:
                # No emission, probability is 0 for this path
                forward[curr_state][i] = 0.0
                continue
            for prev_state in states:
                trans_prob = transition[prev_state].get(curr_state, 1e-10)
                total_prob += forward[prev_state][i - 1] * trans_prob
            forward[curr_state][i] = total_prob * emit_prob

    # Termination: sum probabilities over all states at last position
    prob_x = sum(forward[state][n - 1] for state in states)
    return prob_x

def ParseInputLikelihood(source, from_file=False):
    if from_file:
        with open(source, "r") as file:
            lines = file.read().strip().splitlines()
    else:
        lines = source.strip().splitlines()

    divider_indices = [i for i, line in enumerate(lines) if '--------' in line]
    x = lines[0].strip()
    # alphabet = lines[divider_indices[0] + 1].strip().split()
    states = lines[divider_indices[1] + 1].strip().split()

    transition_header = lines[divider_indices[2] + 1].strip().split()
    transition_lines = lines[divider_indices[2] + 2 : divider_indices[3]]
    transition = {}
    for line in transition_lines:
        parts = line.strip().split()
        row_state = parts[0]
        probs = list(map(float, parts[1:]))
        transition[row_state] = dict(zip(transition_header, probs))

    emission_header = lines[divider_indices[3] + 1].strip().split()
    emission_lines = lines[divider_indices[3] + 2 :]
    emission = {}
    for line in emission_lines:
        parts = line.strip().split()
        row_state = parts[0]
        probs = list(map(float, parts[1:]))
        emission[row_state] = dict(zip(emission_header, probs))

    return x, states, transition, emission


# Example 1
# ----------
input_data = """xzyyzzyzyy
    --------
    x y z
    --------
    A B
    --------
    	A	B
    A	0.303	0.697 
    B	0.831	0.169 
    --------
    	x	y	z
    A	0.533	0.065	0.402 
    B	0.342	0.334	0.324"""

x, states, transition, emission = ParseInputLikelihood(input_data)
probability = LikelihoodOutcome(x, states, transition, emission)
print(probability) # Output: 1.1005510319694847e-06


# Example 2
# ----------
x, states, transition, emission = ParseInputLikelihood("dataset_30328_4.txt", from_file=True)
probability = LikelihoodOutcome(x, states, transition, emission)
print(probability) # Output: 1.933558587497655e-49



# -----------------------------------------------
# Most Likely Outcome
# -----------------------------------------------

import math

def MostLikelyOutcome(n, alphabet, states, transition, emission):
    '''Find a most likely string emitted by an HMM.

    Input: An HMM (Σ, States, Transition, Emission) and an integer n.
    Output: A most likely string x = x1 . . . xn emitted by this HMM, 
    a string maximizing the probability Pr(x) that the HMM will emit x.'''

    dp = {state: [-math.inf] * n for state in states}
    backpointer = {state: [None] * n for state in states}
    emit_choice = {state: [None] * n for state in states}

    initial_log_prob = math.log(1 / len(states))  # Uniform initial probability

    # Initialization (i = 0)
    for state in states:
        max_prob = -math.inf
        best_symbol = None
        for symbol in alphabet:
            emit_prob = emission[state].get(symbol, 1e-10)
            prob = initial_log_prob + math.log(emit_prob)
            if prob > max_prob:
                max_prob = prob
                best_symbol = symbol
        dp[state][0] = max_prob
        emit_choice[state][0] = best_symbol

    # Recursion
    for i in range(1, n):
        for curr_state in states:
            max_prob = -math.inf
            best_prev = None
            best_symbol = None
            for prev_state in states:
                trans_prob = transition[prev_state].get(curr_state, 1e-10)
                log_trans = math.log(trans_prob)
                for symbol in alphabet:
                    emit_prob = emission[curr_state].get(symbol, 1e-10)
                    log_emit = math.log(emit_prob)
                    prob = dp[prev_state][i - 1] + log_trans + log_emit
                    if prob > max_prob:
                        max_prob = prob
                        best_prev = prev_state
                        best_symbol = symbol
            dp[curr_state][i] = max_prob
            backpointer[curr_state][i] = best_prev
            emit_choice[curr_state][i] = best_symbol

    # Termination
    max_final_prob = -math.inf
    last_state = None
    for state in states:
        if dp[state][n - 1] > max_final_prob:
            max_final_prob = dp[state][n - 1]
            last_state = state

    # Backtrack
    emitted_string = [None] * n
    curr_state = last_state
    for i in range(n - 1, -1, -1):
        emitted_string[i] = emit_choice[curr_state][i]
        curr_state = backpointer[curr_state][i]

    return ''.join(emitted_string)


# Example
# --------
alphabet = ['H', 'T']
states = ['F', 'B']
transition = {'F': {'F': 0.6, 'B': 0.4},
              'B': {'F': 0.5, 'B': 0.5}}
emission = {'F': {'H': 0.5, 'T': 0.5},
            'B': {'H': 0.75, 'T': 0.25}}

outcome = MostLikelyOutcome(4, alphabet, states, transition, emission)
print(outcome) # Output: HHHH

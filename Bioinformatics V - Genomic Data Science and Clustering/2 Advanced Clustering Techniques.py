import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics V - Genomic Data Science and Clustering/Data")

# -----------------------------------------------
# Probability of Weighted Coin
# -----------------------------------------------

def ComputeProbabilitySequence(sequence, theta):
    prob = 1.0
    for flip in sequence:
        if flip == 'H':
            prob *= theta
        elif flip == 'T':
            prob *= (1 - theta)
        else:
            raise ValueError("Invalid character in sequence. Use only 'H' and 'T'.")
    return prob

# Example 1
# ----------
sequence = "HTHHT"
theta = 0.7
print(ComputeProbabilitySequence(sequence, theta))  # Output: 0.031


# Example 2
# ----------
sequence = "HHTHHH"
theta = 0.8
print(ComputeProbabilitySequence(sequence, theta))  # Output: 0.066


# -----------------------------------------------
# Hidden Matrix
# -----------------------------------------------

import numpy as np

def ComputeResponsibilityMatrix(data, centers, beta):
    '''Compute the responsibility matrix for all data points and centers.

    Input: data point coordinates as np.array of shape (2,), centers coordinates 
    as np.array of shape (num_centers, 2), and stiffness parameter
    Output: Responsibility matrix.'''

    num_data = data.shape[0]
    num_centers = centers.shape[0]
    
    responsibility_matrix = np.zeros((num_data, num_centers))
    
    for i in range(num_data):
        # Compute distances from data point i to all centers
        dists = np.linalg.norm(centers - data[i], axis=1)
        dists = np.maximum(dists, 1e-10)  # avoid division by zero
        
        # Compute inverse powered distances
        inv_power_dists = 1 / (dists ** (2 * beta))
        
        # Normalize to get responsibilities
        responsibility_matrix[i] = inv_power_dists / np.sum(inv_power_dists)
    
    return responsibility_matrix



# Given the following Data and Centers, compute HiddenMatrix_1,2 (i.e., 
# the responsibility of the first center for the second datapoint) using 
# the partition function with stiffness parameter equal to 1.
# ----------------------------------------------------------------------
data = np.array([[2,8], [2,5], [6,9], [7,5], [5,2]])
centers = np.array([[3,5], [5,4]])
beta = 1
responsibility_matrix = ComputeResponsibilityMatrix(data, centers, beta)

data_index = 1   # second data point
center_index = 0 # first center
result = responsibility_matrix[data_index, center_index]
print(f"Responsibility of center {center_index+1} for data point {data_index+1}: {result:.3f}")
# Output: 0.909





# -----------------------------------------------
# Weighted Center of Gravity
# -----------------------------------------------

def WeightedCenterofGravity(data, weights):
    total_weight = sum(weights)
    weighted_x = sum(w * x for (x, y), w in zip(data, weights))
    weighted_y = sum(w * y for (x, y), w in zip(data, weights))
    center_x = weighted_x / total_weight
    center_y = weighted_y / total_weight
    return center_x, center_y


# Compute the weighted center of gravity corresponding to the first row of HiddenMatrix.  
# ----------------------------------------------------------------------------------------
data = [(2, 6), (4, 9), (5, 7), (6, 5), (8, 3)]
weights = [0.6, 0.1, 0.8, 0.5, 0.7]

center = WeightedCenterofGravity(data, weights)
print(f"{round(center[0], 3)} {round(center[1], 3)}")
# Output: 5.259 5.444



# Compute the weighted center of gravity corresponding to the second row of HiddenMatrix.
# ----------------------------------------------------------------------------------------
data = [(2,8), (2,5), (6,9), (7,5), (5,2)]
weights = [0.5, 0.7, 0.2, 0.6, 0.1]

center = WeightedCenterofGravity(data, weights)
print(f"{round(center[0], 3)} {round(center[1], 3)}")
# Output: 3.952 5.952



# -----------------------------------------------
# Soft k-Means Clustering
# -----------------------------------------------

import math

def SoftkMeans(k, m, beta, data, steps=100):
    '''Soft clustering is needed when data points don't clearly belong to a single 
    cluster, or when there's uncertainty in cluster assignment.
    
    Input: Integers k and m, followed by a stiffness parameter β, followed by a set 
    of points Data in m-dimensional space.
    Output: A set Centers consisting of k points (centers) resulting from applying 
    the soft k-means clustering algorithm. Select the first k points from Data as 
    the first centers for the algorithm and run the algorithm for 100 steps.
    '''
    
    n = len(data)
    centers = [point[:] for point in data[:k]]  # Initialize centers as first k points

    for step in range(steps):
        # E-step: compute responsibilities
        responsibilities = []
        for i in range(n):
            distances = [math.exp(-beta * EuclideanDistance(data[i], centers[j])) for j in range(k)]
            total = sum(distances)
            responsibilities.append([d / total for d in distances])

        # M-step: update centers
        for j in range(k):
            numerator = [0.0] * m
            denominator = 0.0
            for i in range(n):
                for dim in range(m):
                    numerator[dim] += responsibilities[i][j] * data[i][dim]
                denominator += responsibilities[i][j]
            centers[j] = [x / denominator for x in numerator]

    return centers

def EuclideanDistance(p1, p2):
    return sum((a - b) ** 2 for a, b in zip(p1, p2))


# Example
# ---------
raw_input = """2 2
               2.7
               1.3 1.1
               1.3 0.2
               0.6 2.8
               3.0 3.2
               1.2 0.7
               1.4 1.6
               1.2 1.0
               1.2 1.1
               0.6 1.5
               1.8 2.6
               1.2 1.3
               1.2 1.0
               0.0 1.9"""

lines = raw_input.strip().split('\n')
k, m = map(int, lines[0].split())
beta = float(lines[1])
data = [list(map(float, line.split())) for line in lines[2:]]

centers = SoftkMeans(k, m, beta, data)
for center in centers:
    print(' '.join(f'{coord:.3f}' for coord in center))

# Output:
    # 1.803 2.857
    # 1.059 1.142


# -----------------------------------------------
# Distance Between Clusters
# -----------------------------------------------

def AverageDistance(D, C1, C2, labels):
    # Map cluster labels to their indices in D matrix
    idx_C1 = [labels.index(x) for x in C1]
    idx_C2 = [labels.index(x) for x in C2]

    total = 0
    count = 0
    for i in idx_C1:
        for j in idx_C2:
            total += D[i][j]
            count += 1
    return total / count if count > 0 else 0


def MinimumDistance(D, C1, C2, labels):
    idx_C1 = [labels.index(x) for x in C1]
    idx_C2 = [labels.index(x) for x in C2]

    distances = []
    for i in idx_C1:
        for j in idx_C2:
            distances.append(D[i][j])
    return min(distances) if distances else 0


# If C1 = {i, l} and C2 = {j, k}, compute Davg(C1, C2).
# ------------------------------------------------------
D = [[0, 20, 9, 11],
     [20, 0, 17, 11],
     [9, 17, 0, 8],
     [11, 11, 8, 0]]

labels = ['i', 'j', 'k', 'l']
C1 = ['i', 'l']
C2 = ['j', 'k']

davg = AverageDistance(D, C1, C2, labels)
print(round(davg, 3))  # Output: 12

dmin = MinimumDistance(D, C1, C2, labels)
print(round(dmin, 3))  # Output: 8


# -----------------------------------------------
# Hierarchical Clustering
# -----------------------------------------------

import math

def HierarchicalClustering(n, D, linkage='average'):
    '''Implements Hierarchical Clustering using average linkage (Davg), accounting for
    subclusters within clusters.
    
    Generates n different partitions of the underlying data into clusters, all represented 
    by a tree in which each node is labeled by a cluster of genes. The first partition has 
    n single-element clusters represented by the leaves of the tree, with each element forming 
    its own cluster. The second partition merges the two “closest” clusters into a single cluster 
    consisting of two elements. In general, the i-th partition merges the two closest clusters 
    from the (i - 1)-th partition and has n - i + 1 clusters.

    Input: An integer n, followed by an n x n distance matrix.
    Output: The result of applying HierarchicalClustering to this distance matrix (using Davg),
    with each newly created cluster listed on each line.'''
    

    clusters = [[i] for i in range(n)]
    result = []

    def AverageDistance(c1, c2):
        total = 0
        count = 0
        for i in c1:
            for j in c2:
                if i != j:
                    total += D[i][j]
                    count += 1
        return total / count if count > 0 else 0

    def MinimumDistance(c1, c2):
        return min(D[i][j] for i in c1 for j in c2 if i != j)

    # Choose the appropriate distance function
    if linkage == 'average':
        Distance = AverageDistance
    elif linkage == 'minimum':
        Distance = MinimumDistance
    else:
        raise ValueError("Unsupported linkage type. Use 'average' or 'minimum'.")

    while len(clusters) > 1:
        min_dist = math.inf
        to_merge = (0, 1)

        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                d = Distance(clusters[i], clusters[j])
                if d < min_dist:
                    min_dist = d
                    to_merge = (i, j)

        i, j = to_merge
        new_cluster = sorted(clusters[i] + clusters[j])
        result.append(' '.join(str(x + 1) for x in new_cluster))

        # Merge clusters
        clusters.pop(max(i, j))
        clusters.pop(min(i, j))
        clusters.append(new_cluster)

    return result


# Example 1
# -----------
n = 7
D = [[0.00, 0.74, 0.85, 0.54, 0.83, 0.92, 0.89],
     [0.74, 0.00, 1.59, 1.35, 1.20, 1.48, 1.55],
     [0.85, 1.59, 0.00, 0.63, 1.13, 0.69, 0.73],
     [0.54, 1.35, 0.63, 0.00, 0.66, 0.43, 0.88],
     [0.83, 1.20, 1.13, 0.66, 0.00, 0.72, 0.55],
     [0.92, 1.48, 0.69, 0.43, 0.72, 0.00, 0.80],
     [0.89, 1.55, 0.73, 0.88, 0.55, 0.80, 0.00]]

clusters = HierarchicalClustering(n, D, linkage='average')
print('\n'.join(clusters) + '\n')
# Output: 
    # 4 6
    # 5 7
    # 3 4 6
    # 1 2
    # 3 4 5 6 7
    # 1 2 3 4 5 6 7
    
    
    
# Example 2
# ----------
with open("dataset_30177_7.txt", "r") as file:
    lines = file.read().strip().split('\n')
    n = int(lines[0])
    D = [list(map(float, line.strip().split())) for line in lines[1:n+1]]
    
clusters = HierarchicalClustering(n, D, linkage='average')

with open("output.txt", "w") as file:
    file.write('\n'.join(clusters) + '\n')

# Output: 
    # 3 4
    # 8 16
    # 1 11
    # 8 16 19
    # 8 13 16 19
    # 12 20
    # 3 4 9
    # 8 13 14 16 19
    # 15 17
    # 10 18
    # 1 10 11 18
    # 2 6
    # 5 15 17
    # 7 8 13 14 16 19
    # 2 3 4 6 9
    # 2 3 4 6 9 12 20
    # 5 7 8 13 14 15 16 17 19
    # 1 5 7 8 10 11 13 14 15 16 17 18 19
    # 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20

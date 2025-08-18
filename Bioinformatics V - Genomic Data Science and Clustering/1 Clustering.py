import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics V - Genomic Data Science and Clustering/Data")

# -----------------------------------------------
# k-Centers Clustering: Largest Distance
# -----------------------------------------------

import math

def MaxDistance(Data, Centers):
    '''Compute the maximum distance from any point in Data to its closest center in Centers.
    This function calculates the maximum of the minimum Euclidean distances 
    between each point in the Data set and the nearest point in the Centers set.

    Input: A list of m-dimensional data points and a list of m-dimensional center points.
    Output: The maximum distance from a data point to its closest center.'''
    
    def EuclideanDistance(point1, point2):
        '''Compute the Euclidean distance between two points in m-dimensional space.

        Input: Two points represented as a tuple or list of coordinates.
        Output: The Euclidean distance between Point1 and Point2.'''
        
        return math.sqrt(sum((v - w) ** 2 for v, w in zip(point1, point2)))

    def DistanceToCenters(DataPoint, Centers):
        '''Find the minimum distance from a given data point to any of the center points.

        Input: A single point from the Data set and a list of center points.
        Output: The minimum Euclidean distance from DataPoint to any point in Centers.'''
        
        return min(EuclideanDistance(DataPoint, center) for center in Centers)

    return max(DistanceToCenters(point, Centers) for point in Data)


# Example 1
# -----------
Data = [(1, 6), (1, 3), (3, 4), (5, 6), (5, 2), (7, 1), (8, 7), (10, 3)]
Centers =  [(2, 4), (6, 7), (7, 3)]
print(MaxDistance(Data, Centers)) # Output: 3


# Example 2
# -----------
Data = [(2, 8), (2, 5), (6, 9), (7, 5), (5, 2)]
Centers = [(3, 5), (5, 4)]
print(MaxDistance(Data, Centers)) # Output: 5



# -----------------------------------------------
# k-Centers Clustering
# -----------------------------------------------

import itertools
import math

def kCenterClustering(Data, k):
    '''Solves the k-Center Clustering Problem using brute-force search over all possible center combinations.
    Only use for small datasets; otherwise use a greedy heuristic.

    Input: A list of data points in m-dimensional space and an integer specifying the number of centers to select.
    Output: A list of k centers from Data that minimizes MaxDistance(Data, Centers).
    '''
    BestCenters = None
    MinMaxDist = float('inf')

    for CentersCandidate in itertools.combinations(Data, k):
        CurrentMaxDist = MaxDistance(Data, CentersCandidate)
        if CurrentMaxDist < MinMaxDist:
            MinMaxDist = CurrentMaxDist
            BestCenters = CentersCandidate

    return list(BestCenters)


def MaxDistance(Data, Centers):
    '''Compute the maximum distance from any point in Data to its closest center in Centers.
    This function calculates the maximum of the minimum Euclidean distances 
    between each point in the Data set and the nearest point in the Centers set.

    Input: A list of m-dimensional data points and a list of m-dimensional center points.
    Output: The maximum distance from a data point to its closest center.'''
    
    def EuclideanDistance(point1, point2):
        '''Compute the Euclidean distance between two points in m-dimensional space.

        Input: Two points represented as a tuple or list of coordinates.
        Output: The Euclidean distance between Point1 and Point2.'''
        
        return math.sqrt(sum((v - w) ** 2 for v, w in zip(point1, point2)))

    def DistanceToCenters(DataPoint, Centers):
        '''Find the minimum distance from a given data point to any of the center points.

        Input: A single point from the Data set and a list of center points.
        Output: The minimum Euclidean distance from DataPoint to any point in Centers.'''
        
        return min(EuclideanDistance(DataPoint, center) for center in Centers)

    return max(DistanceToCenters(point, Centers) for point in Data)


# Example
# ----------
Data = [(1, 6), (1, 3), (3, 4), (5, 6), (5, 2), (7, 1), (8, 7), (10, 3)]
k = 3

OptimalCenters = kCenterClustering(Data, k)
print('Optimal Centers:', OptimalCenters)
print('Max Distance:', MaxDistance(Data, OptimalCenters))


# -----------------------------------------------
# Find Center of Gravity
# -----------------------------------------------

def CenterOfGravity(Data):
    n = len(Data)
    sum_x = sum(point[0] for point in Data)
    sum_y = sum(point[1] for point in Data)
    sum_z = sum(point[2] for point in Data)
    
    centroid = round(sum_x / n), round(sum_y / n), round(sum_z / n)
    return f"({centroid[0]}, {centroid[1]}, {centroid[2]})"


# Example 
# ---------
Data = [(17, 0, -4), (3, 14, 23), (9, 7, 16), (7, 3, 5)]
print(CenterOfGravity(Data)) # Output: (9, 6, 10)



# -----------------------------------------------
# Farthest Travel Heuristic (for k-Center Clustering)
# -----------------------------------------------

import random

def FarthestFirstTraversal(Data, k, m):
    '''Applies the Farthest-First Traversal heuristic to select k centers from the dataset.
    The first point in Data is used to initialize the Centers.
    
    k-Center Clustering is rarely used for gene expression analysis, since biologists are usually 
    interested in analyzing typical rather than maximum deviations (the latter may correspond 
    to outliers representing experimental errors).

    Input: A list of points in m-dimensional space (each point is a list or tuple of length m)
    and an integer specifying the number of centers to select.

    Output: A list of k points selected from Data as centers using the Farthest-First Traversal algorithm.
    '''
    
    # Initialize Centers with the first point in Data
    Centers = [Data[0]]

    # Continue selecting farthest points until we have K centers
    while len(Centers) < k:
        FarthestPoint = max(Data,
            key=lambda point: min(
                sum((v - w) ** 2 for v, w in zip(point, center)) ** 0.5
                for center in Centers))
        Centers.append(FarthestPoint)

    return Centers

def PrintPoints(points):
    for point in points:
        print(' '.join(str(coord) for coord in point))



# Example 1
# ------------
k, m = 3, 2
Data = [(0.0, 0.0), (5.0, 5.0), (0.0, 5.0), (1.0, 1.0), (2.0, 2.0), (3.0, 3.0), (1.0, 2.0)]
Centers = FarthestFirstTraversal(Data, k, m)
PrintPoints(Centers) # Output: 0.0 0.0    5.0 5.0    0.0 5.0


# Example 2
# ------------
with open("dataset_30181_2.txt", "r") as file:
        first_line = file.readline().strip()
        k, m = map(int, first_line.split())

        Data = []
        for line in file:
            coords = list(map(float, line.strip().split()))
            if len(coords) != m:
                raise ValueError(f'Expected {m} coordinates, got {len(coords)} in line: {line}')
            Data.append(tuple(coords))
Centers = FarthestFirstTraversal(Data, k, m)
PrintPoints(Centers) 
# Output: 29.9 21.2 1.3 10.7
    # 1.6 5.2 31.8 4.3
    # 0.1 2.6 1.3 22.7
    # 0.2 26.1 7.0 1.9
    # 27.0 0.2 18.2 7.0
    # 4.7 17.1 24.1 26.1
    # 22.8 28.1 8.9 32.4
    


# -----------------------------------------------
# Calculate Distortion
# -----------------------------------------------

def Distortion(Data, Centers):
    def SquaredDistance(P1, P2):
        return sum((a - b) ** 2 for a, b in zip(P1, P2))
    
    n = len(Data)
    total = 0
    for point in Data:
        distances = [SquaredDistance(point, center) for center in Centers]
        total += min(distances)
    return total / n


# Example
# ---------
Data = [(2, 8), (2, 5), (6, 9), (7, 5), (5, 2)]
Centers = [(3, 5), (5, 4)]
print(Distortion(Data, Centers)) # Output: 9



# -----------------------------------------------
# Calculate Squared Error Distortion
# -----------------------------------------------

def SquaredErrorDistortion(Data, Centers):
    '''Compute the mean squared distance from each data point to its nearest center
    for k-Centers clustering.

    Input: List of n data points and ist of k center points (tuples/lists).
    Output: The squared error distortion.'''
    
    n = len(Data)
    total_squared_distance = 0.0

    for point in Data:
        # Compute squared distances to all centers
        squared_distances = [
            sum((v - w) ** 2 for v, w in zip(point, center)) for center in Centers]
        # Add the minimum squared distance for this point
        total_squared_distance += min(squared_distances)

    return total_squared_distance / n

def ReadSquaredErrorInput(lines):
    '''Parse input lines for the Squared Error Distortion problem.

    Input: List of strings representing lines of input.
    Output: Number of centers, dimension of points, list of k center points, and lsit of data points.
    '''
    
    k, m = map(int, lines[0].strip().split())

    Centers = []
    for i in range(1, k+1):
        coords = tuple(map(float, lines[i].strip().split()))
        if len(coords) != m:
            raise ValueError(f'Expected {m} coordinates for center, got {len(coords)}')
        Centers.append(coords)

    Data = []
    for line in lines[k+1:]:
        coords = tuple(map(float, line.strip().split()))
        if len(coords) != m:
            raise ValueError(f'Expected {m} coordinates for data point, got {len(coords)}')
        Data.append(coords)

    return k, m, Centers, Data

def ReadSquaredErrorInputFromFile(filename):
    '''Reads k, m, Centers, and Data from a file for the squared error distortion problem.
    Lines with '--------' are ignored.

    Input: Path to input file.
    Output: Number of centers, dimension of points, list of k center points, and lsit of data points.
    '''
    
    with open(filename, "r") as file:
        lines = [line.strip() for line in file if line.strip() and line.strip() != '--------']

    k, m = map(int, lines[0].split())

    Centers = []
    for i in range(1, k + 1):
        coords = tuple(map(float, lines[i].split()))
        if len(coords) != m:
            raise ValueError(f'Expected {m} coordinates for center, got {len(coords)}')
        Centers.append(coords)

    Data = []
    for line in lines[k + 1:]:
        coords = tuple(map(float, line.split()))
        if len(coords) != m:
            raise ValueError(f'Expected {m} coordinates for data point, got {len(coords)}')
        Data.append(coords)

    return k, m, Centers, Data


# Example 1
# ----------
input_lines = ["2 2", "2.31 4.55", "5.96 9.08", "--------", "3.42 6.03",
    "6.23 8.25", "4.76 1.64", "4.47 4.33", "3.95 7.61", "8.93 2.97",
    "9.74 4.03", "1.73 1.28", "9.72 5.01", "7.27 3.77"]

filtered_lines = [line for line in input_lines if line.strip() != '--------']
k, m, Centers, Data = ReadSquaredErrorInput(filtered_lines)

distortion = SquaredErrorDistortion(Data, Centers)
print(distortion) # Output: 18.245559999999994


# Example 2
# ----------
k, m, Centers, Data = ReadSquaredErrorInputFromFile("dataset_30170_3.txt")
distortion = SquaredErrorDistortion(Data, Centers)
print(distortion) # Output: 40.84936570086142



# -----------------------------------------------
# Clustering in One-Dimensional Space for Any k
# -----------------------------------------------

def KMeans1D(points, k):
    '''Solve the 1D k-means clustering problem o.
    
    Input: List of 1D points and number of clusters.
    Output: minimum squared error distortion and list of tuples 
    (start_index, end_index) for each cluster.'''
    
    points = sorted(points)
    n = len(points)
    cost = PrecomputeCosts(points)
    
    DP = [[float('inf')] * (k + 1) for _ in range(n + 1)]
    backtrack = [[-1] * (k + 1) for _ in range(n + 1)]
    
    DP[0][0] = 0.0
    
    for i in range(1, n + 1):
        for j in range(1, min(k, i) + 1):
            for m in range(j - 1, i):
                current_cost = DP[m][j - 1] + cost[m][i - 1]
                if current_cost < DP[i][j]:
                    DP[i][j] = current_cost
                    backtrack[i][j] = m
    
    # Reconstruct clusters
    clusters = []
    idx = n
    for j in range(k, 0, -1):
        start = backtrack[idx][j]
        clusters.append((start, idx - 1))
        idx = start
    clusters.reverse()
    
    return DP[n][k], clusters, points

def PrecomputeCosts(points):
    '''Precompute cost(l, r): squared error of clustering points[l:r+1] into one cluster.
    Uses prefix sums and prefix squared sums for O(1) interval cost calculation.
    '''
    n = len(points)
    prefix_sum = [0.0] * (n + 1)
    prefix_sq_sum = [0.0] * (n + 1)
    
    for i in range(1, n + 1):
        prefix_sum[i] = prefix_sum[i-1] + points[i-1]
        prefix_sq_sum[i] = prefix_sq_sum[i-1] + points[i-1] ** 2
    
    cost = [[0.0] * n for _ in range(n)]
    
    for l in range(n):
        for r in range(l, n):
            length = r - l + 1
            s = prefix_sum[r+1] - prefix_sum[l]
            sq_s = prefix_sq_sum[r+1] - prefix_sq_sum[l]
            mean = s / length
            # sum of squared distances = sum(x_i^2) - 2*mean*sum(x_i) + length * mean^2
            cost[l][r] = sq_s - 2 * mean * s + length * mean * mean
    
    return cost


# Example
# ------------
points_1d = [0.0, 5.0, 0.0, 1.0, 2.0, 3.0, 1.0, 2.0]
k = 3
distortion, clusters, sorted_points = KMeans1D(points_1d, k)
print(f'Minimum squared error distortion: {distortion:.4f}')
print('Clusters (index ranges in sorted data):')
for start, end in clusters:
    cluster_points = sorted_points[start:end+1]
    print(f'Indices {start} to {end}: points {cluster_points}')

# Output: Minimum squared error distortion: 1.6667
# Clusters (index ranges in sorted data):
    # Indices 0 to 3: points [0.0, 0.0, 1.0, 1.0]
    # Indices 4 to 6: points [2.0, 2.0, 3.0]
    # Indices 7 to 7: points [5.0]



# -----------------------------------------------
# K-Means Clustering via Lloyd's Algorithm
# -----------------------------------------------

def LloydAlgorithm(k, m, Data):
    ''' Run Lloyd's algorithm for k-means clustering.

    Input: Number of clusters, dimension of points, and List of data points.
    Output: Final list of k centers after convergence.'''
    
    Centers = Data[:k]

    while True:
        clusters = AssignPointsToClusters(Data, Centers)
        new_centers = ComputeNewCenters(clusters, m)

        # Check convergence
        if all(
            all(abs(a - b) < 1e-6 for a, b in zip(c1, c2))
            for c1, c2 in zip(Centers, new_centers)):
            break
        Centers = new_centers

    return Centers

def ParseInputLines(lines):
    k, m = map(int, lines[0].strip().split())
    Data = [tuple(map(float, line.strip().split())) for line in lines[1:]]
    return k, m, Data

def ReadInputFromFile(filename):
    with open(filename, "r") as file:
        lines = file.readlines()

    data_lines = lines[1:]
    Data = []
    for line in data_lines:
        parts = line.strip().split()
        try:
            ratios = list(map(float, parts[3:10]))
            if len(ratios) == 7:  # only keep valid 7D points
                Data.append(tuple(ratios))
        except ValueError:
            continue  # skip malformed lines
    return Data

def AssignPointsToClusters(Data, Centers):
    '''Assigns each data point to the closest center.

    Input: List of data points and list of current center points.
    Output: A list of k lists, where each list holds the points assigned to one cluster.
    '''
    
    clusters = [[] for _ in Centers]
    for point in Data:
        distances = [sum((x - c) ** 2 for x, c in zip(point, center)) for center in Centers]
        closest_index = distances.index(min(distances))
        clusters[closest_index].append(point)
    return clusters

def ComputeNewCenters(clusters, m):
    '''Computes new centers as the mean of each cluster.

    Input: List of clusters (each a list of points) and dimension of the points.
    Output: List of new center points.'''
    
    new_centers = []
    for cluster in clusters:
        if cluster:
            mean = tuple(
                sum(point[i] for point in cluster) / len(cluster)
                for i in range(m))
            new_centers.append(mean)
        else:
            # Empty cluster, leave unchanged (could also reinitialize randomly)
            new_centers.append((0.0,) * m)
    return new_centers

def PrintCenters(Centers):
    for center in Centers:
        print(' '.join(f'{coord:.3f}' for coord in center))


# Example 1
# -----------
input_lines = ["2 2", "1.3 1.1", "1.3 0.2", "0.6 2.8", "3.0 3.2", "1.2 0.7", "1.4 1.6", 
               "1.2 1.0", "1.2 1.1", "0.6 1.5", "1.8 2.6", "1.2 1.3", "1.2 1.0", "0.0 1.9"]
k, m, Data = ParseInputLines(input_lines)
Centers = LloydAlgorithm(k, m, Data)
PrintCenters(Centers)
# Output:
    # 1.800 2.867
    # 1.060 1.140


# Example 2
# -----------
k, m, Data = ReadInputFromFile("dataset_30171_3.txt")
Centers = LloydAlgorithm(k, m, Data)

with open("output.txt", "w") as file:
    for center in Centers:
        line = ' '.join(f'{coord:.3f}' for coord in center)
        file.write(f'{line}\n')
# Output:
    # 18.947 11.999 10.771
    # 6.337 20.298 5.807
    # 4.159 3.215 4.690
    # 15.406 4.758 3.927
    # 4.628 10.615 5.053
    # 10.993 3.966 13.155
    # 3.875 7.569 16.789



# -----------------------------------------------
# Probability of An Empty Cluster
# -----------------------------------------------

import math

def ProbabilityAtLeastOneEmpty(numClusters, numPoints):
    '''Estimate the probability that at least one cluster has no points when placing
    `numPoints` balls into `numClusters` bins uniformly at random with replacement.
    This is equivalent to computing the probability that not all bins receive at least one ball.

    The probability that all bins are hit is:
        P(all bins hit) = (n! * S(k, n)) / n^k
        
        where:
            - n is the number of bins
            - k is the number of balls
            - S(k, n) is the Stirling number of the second kind

    Finally, we return:
        P(at least one bin is empty) = 1 - P(all bins hit)

    Input: Number of clusters (e.g., bins) and number of points (e.g., balls).

    Output: Probability that at least one bin is empty.
    '''
    
    def StirlingSecondKind(k, n):
        '''Compute the Stirling number of the second kind S(k, n).'''
        sum_ = 0
        for j in range(n + 1):
            sum_ += (-1) ** (n - j) * math.comb(n, j) * j ** k
        return sum_ // math.factorial(n)

    if numPoints < numClusters:
        # It's impossible for all bins to get a ball if there are fewer balls
        return 1.0

    totalAssignments = numClusters ** numPoints
    ontoFunctions = math.factorial(numClusters) * StirlingSecondKind(numPoints, numClusters)
    probAllBinsHit = ontoFunctions / totalAssignments
    return 1 - probAllBinsHit


# Example
# --------
prob = ProbabilityAtLeastOneEmpty(5, 5)
print(f"Probability at least one clump gets no center: {prob:.4f}")
# Output: 0.9616




# -----------------------------------------------
# Ways to Cluster x Points in 
# -----------------------------------------------

from math import comb

def StirlingSecondKind(n, k):
    '''Compute the number of ways are there to cluster 5 points into 2 clusters'''
    
    def Factorial(X):
        if X == 0 or X == 1:
            return 1
        else:
            return X * Factorial(X - 1)
    
    total = 0
    for j in range(k + 1):
        total += ((-1) ** (k - j)) * comb(k, j) * (j ** n)
    return total // Factorial(k)


# Example
# ---------
n, k = 5, 2
result = StirlingSecondKind(n, k)
print(f"The number of ways to cluster {n} points into {k} clusters is: {result}")
# Output: 15



# -----------------------------------------------
# k-Means++Initializer
# -----------------------------------------------

import numpy as np

def kMeansPlusPlusInitializer(data, k):
    n_samples = data.shape[0]
    centers = []

    # Step 1: Choose one center randomly from data
    first_center_idx = np.random.choice(n_samples)
    centers.append(data[first_center_idx])

    while len(centers) < k:
        # Compute squared distances from each point to the nearest center
        dist_sq = np.array([
            min(np.sum((x - center)**2) for center in centers)
            for x in data
        ])

        # Compute probabilities proportional to squared distances
        prob = dist_sq / dist_sq.sum()

        # Choose the next center index based on the computed probabilities
        next_center_idx = np.random.choice(n_samples, p=prob)
        centers.append(data[next_center_idx])

    return np.array(centers)


# Example
# ---------
data = np.array([[1, 2], [3, 4], [5, 6], [7, 8], [9, 10]])
k = 3
centers = kMeansPlusPlusInitializer(data, k)
print("Centers:\n", *centers) # Output: [1 2] [7 8] [3 4]

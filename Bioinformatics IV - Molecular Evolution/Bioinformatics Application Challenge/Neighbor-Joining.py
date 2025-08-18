import numpy as np
import pandas as pd

# Original distance matrix
# -------------------------
labels = ['W', 'X', 'Y', 'Z']
D = pd.DataFrame([[0, 11, 2, 16],
                  [11, 0, 13, 15],
                  [2, 13, 0, 9],
                  [16, 15, 9, 0]], 
                  index=labels, columns=labels)

print("Original Distance Matrix D:")
print(D, "\n")


# Construct the neighbor-joining matrix D* from the distance matrix D given above.
# ---------------------------------------------------------------------------------

def RowSums(df):
    return df.sum(axis=1)

def ConstructNJMatrix(D):
    n = len(D)
    r = RowSums(D)
    DStar = pd.DataFrame(index=D.index, columns=D.columns, dtype=float)
    for i in D.index:
        for j in D.columns:
            if i != j:
                DStar.loc[i, j] = (n - 2) * D.loc[i, j] - r[i] - r[j]
            else:
                DStar.loc[i, j] = np.nan
    return DStar


DStar = ConstructNJMatrix(D)
print("Neighbor-Joining Matrix D*:")
print(DStar, "\n")



# Construct a 3x3 distance matrix D2 using the neighbor-joining matrix D*. 
# “A” refers to joining W and Y as neighbors of the tree. You will need to 
# fill in the two branches that were not joined in addition to all distances 
# in the resulting 3x3 distance matrix.
# ---------------------------------------------------------------------------

def JoinPair(D, i, j, newLabel):
    r = RowSums(D)
    n = len(D)

    # Branch lengths from i and j to new node
    LimbI = 0.5 * D.loc[i, j] + (r[i] - r[j]) / (2 * (n - 2))
    LimbJ = D.loc[i, j] - LimbI

    # Distances from new node to others
    NewDistances = {}
    for k in D.index:
        if k != i and k != j:
            d = (D.loc[i, k] + D.loc[j, k] - D.loc[i, j]) / 2
            NewDistances[k] = d

    # New labels: new node and all other non-joined nodes
    NewLabels = [newLabel] + [k for k in D.index if k != i and k != j]
    DNew = pd.DataFrame(index=NewLabels, columns=NewLabels, dtype=float)

    for k in DNew.index:
        for l in DNew.columns:
            if k == l:
                DNew.loc[k, l] = 0
            elif k == newLabel:
                DNew.loc[k, l] = NewDistances[l]
            elif l == newLabel:
                DNew.loc[k, l] = NewDistances[k]
            else:
                DNew.loc[k, l] = D.loc[k, l]

    return DNew, LimbI, LimbJ


D2, LimbW, LimbY = JoinPair(D, 'W', 'Y', 'A')
print("New 3x3 Distance Matrix D2 (after joining W and Y into A):")
print(D2, "\n")
print(f"Branch lengths:\n  W — A: {LimbW}\n  Y — A: {LimbY}\n")



# Construct the neighbor-joining matrix D2* from the distance matrix D2.
# ----------------------------------------------------------------------
D2Star = ConstructNJMatrix(D2)
print("Neighbor-Joining Matrix D2*:")
print(D2Star, "\n")



# Construct a 2x2 distance matrix D3 using the neighbor-joining matrix D2*. 
# “B” refers to the new leaf formed by joining A and X in the neighbor-joining 
# algorithm. You will need to fill in the branch that was not joined in addition 
# to the distance in the resulting 2x2 distance matrix.
# ----------------------------------------------------------------------------
D3, LimbA, LimbX = JoinPair(D2, 'A', 'X', 'B')
print("New 2x2 Distance Matrix D3 (after joining A and X into B):")
print(D3, "\n")
print(f"Branch lengths:\n  A — B: {LimbA}\n  X — B: {LimbX}\n")




# Tree Diagram
# -------------------------------------------------
import networkx as nx
import matplotlib.pyplot as plt

G = nx.Graph()

# Add edges with branch lengths
G.add_edge('W', 'A', weight=LimbW)
G.add_edge('Y', 'A', weight=LimbY)
G.add_edge('A', 'B', weight=LimbA)
G.add_edge('X', 'B', weight=LimbX)
G.add_edge('B', 'Z', weight=D3.loc['B', 'Z'])

# Position the nodes for a tree layout
pos = nx.spring_layout(G, seed=42)

# Draw nodes
nx.draw_networkx_nodes(G, pos, node_size=1000, node_color='lightblue')

# Draw edges
nx.draw_networkx_edges(G, pos)

# Draw labels
nx.draw_networkx_labels(G, pos, font_size=12, font_weight='bold')

# Draw edge labels (branch lengths)
edge_labels = nx.get_edge_attributes(G, 'weight')
formatted_labels = {edge: f"{w:.2f}" for edge, w in edge_labels.items()}
nx.draw_networkx_edge_labels(G, pos, edge_labels=formatted_labels, font_size=10)

# Show the tree
plt.title("Neighbor-Joining Tree")
plt.axis('off')
plt.tight_layout()
plt.show()





import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

#Part 1: Defining the graph of protein folding states

#Example graph of 6 states (S1-S6)
#ΔG (kcal/mol) = energy change from one folding step to the next
edges = [
    ('U', 'A', -2.0), #'U' = Unfolded protein state / Starting point
    ('A', 'B', -1.5),
    ('B', 'C', -3.0),
    ('C', 'F', -2.0),
    ('B', 'D', 1.0),
    ('D', 'E', -4.0),
    ('E', 'F', 2.0),
    ('C', 'A', -1.5), #possible unstable back-loop
]

#Creating a directed weighted graph
G = nx.DiGraph()
G.add_weighted_edges_from(edges)

#Printing graph contents
print("\nNodes:", list(G.nodes))
print("Edges (with ΔG values):")
for u, v, w in G.edges(data='weight'):
    print(f" {u} -> {v}: {w} kcal/mol")

#Part 2: Visualizing the graph

plt.figure(figsize=(8, 6))
pos = nx.spring_layout(G, seed=42)
nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=1200, arrowsize=20)
edge_labels = nx.get_edge_attributes(G, 'weight')
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='red')
plt.title("Protein Folding State Graph (ΔG in kcal/mol)")
plt.show()

#Part 3: Implementing the Bellman-Ford Algorithm

def bellman_ford(graph, source):
    #Initialize distances and predecessors
    distance = {node: float('inf') for node in graph.nodes}
    predecessor = {node: None for node in graph.nodes}
    distance[source] = 0

    #Get the number of nodes
    N = len(graph.nodes)

    #Main relaxation loop
    for _ in range(N-1):
        for u, v, data in graph.edges(data=True):
            weight = data['weight']
            if distance[u] + weight < distance[v]:
                distance[v] = distance[u] + weight
                predecessor[v] = u
    
    #Check for negative weight cycles
    negative_cycle = []
    for u, v, data in graph.edges(data=True):
        weight = data['weight']
        if distance[u] + weight < distance[v]:
            negative_cycle.append((u, v))
    
    return distance, predecessor, negative_cycle

#Part 4: Running the algorithm

source = "U" #unfolded state
distance, predecessor, negative_cycle = bellman_ford(G, source)

#Output results
print("Minimum energy distances (ΔG):")
for node, dist in distance.items():
    print(f"{source} -> {node}: {dist:.2f}")

if negative_cycle:
    print("\n Negative Energy Cycle Detected:")
    for edge in negative_cycle:
        print(f" Cycle edge: {edge}")
else:
    print("\nNo negative cycles detected - stable folding pathway")

#Part 5: Reconstructing and visualizing the lowest-energy path

def reconstruct_path(predecessor, start, end):
    path = []
    while end is not None:
        path.insert(0, end)
        end = predecessor[end]
    if path[0] == start:
        return path
    else:
        return None

target = 'F' #Folded / Final state
path = reconstruct_path(predecessor, source, target)

print(f"\nLowest energy folding pathway from {source} to {target}: {path}")
print(f"Total ΔG = {distance[target]:.2f} kcal/mol")

#Highlighting the path
if path:
    edge_path = list(zip(path[:-1], path[1:]))
    plt.figure(figsize=(8, 6))
    nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=1200, arrowsize=20)
    nx.draw_networkx_edges(G, pos, edgelist=edge_path, edge_color='green', width=3)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='red')
    plt.title(f"Lowest-Energy Folding Pathway: {path}\nTotal ΔG = {distance[target]:.2f} kcal/mol")
    plt.show()

#Part 6: Complexity and Analysis
num_nodes = len(G.nodes)
num_edges = len(G.edges)
complexity = f"O(V*E) = O({num_nodes}*{num_edges}) = O({num_nodes*num_edges})"
print(f"\nAlgorithm complexity: {complexity}")
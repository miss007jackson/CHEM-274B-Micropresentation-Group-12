import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio.PDB import PDBParser


# Step 1: Define the ProteinGraph class
class ProteinGraph:
    def __init__(self, edges=None, pdb_file=None, name="ProteinGraph"):
        self.G = nx.DiGraph()
        self.name = name
        self.pos = None

        # Initialize graph from edges or PDB
        if edges:
            self.G.add_weighted_edges_from(edges)
        elif pdb_file:
            self.build_from_pdb(pdb_file)
        else:
            raise ValueError("Must provide either edges or pdb_file")

        # Initialize edge labels for visualization
        self.edge_labels = {(u, v): f"{w:+.1f} kcal/mol"
                            for u, v, w in self.G.edges(data="weight")}

    # Step 2: Build a graph from PDB file
    def build_from_pdb(self, pdb_file):
        structure = PDBParser(QUIET=True).get_structure("protein", pdb_file)
        residues = [res for res in structure.get_residues() if res.has_id('CA')]

        self.G.add_weighted_edges_from(
            [(residues[i].get_id()[1], residues[i+1].get_id()[1], -0.5)
             for i in range(len(residues) - 1)]
        )


    # Step 3: Bellman-Ford with correct negative cycle detection
    def bellman_ford(self, source):
        nodes = list(self.G.nodes)
        distance = {node: float('inf') for node in nodes}
        predecessor = {node: None for node in nodes}
        distance[source] = 0

        # Standard relaxation (|V|-1 iterations)
        for _ in range(len(nodes) - 1):
            for u, v, w in self.G.edges(data="weight"):
                if distance[u] + w < distance[v]:
                    distance[v] = distance[u] + w
                    predecessor[v] = u

        # Detect any edge that can still be relaxed → negative cycle exists
        cycle_vertex = None
        for u, v, w in self.G.edges(data="weight"):
            if distance[u] + w < distance[v]:
                cycle_vertex = v
                break

        # No negative cycle detected
        if cycle_vertex is None:
            return distance, predecessor, []

        # Move cycle_vertex backward |V| steps to ensure we are inside the cycle
        for _ in range(len(nodes)):
            cycle_vertex = predecessor[cycle_vertex]

        # Reconstruct the actual cycle by following predecessors
        cycle_edges = []
        cur = cycle_vertex
        while True:
            prev = predecessor[cur]
            cycle_edges.append((prev, cur))
            cur = prev
            if cur == cycle_vertex:
                break

        return distance, predecessor, cycle_edges

    # Step 4: Visualize the graph
    def visualize(self, filename="graph.png"):
        if not self.pos:
            self.pos = nx.spring_layout(self.G, seed=42, k=1.2)

        # Default colors based on sign of weight
        edge_colors = ["red" if w < 0 else "blue"
                       for _, _, w in self.G.edges(data="weight")]

        # Highlight the actual negative cycle edges in black
        _, _, neg_cycle = self.bellman_ford(next(iter(self.G.nodes)))
        neg_cycle_set = set(neg_cycle)

        edge_colors = [
            "black" if (u, v) in neg_cycle_set else c
            for (u, v), c in zip(self.G.edges(), edge_colors)
        ]

        # Draw the graph
        fig, ax = plt.subplots(figsize=(12, 8))
        nx.draw(self.G, pos=self.pos, node_color='#A7D5F2', node_size=1800,
                edgecolors='black', linewidths=1.3, ax=ax, with_labels=True,
                font_weight="bold", font_size=14)
        nx.draw_networkx_edges(self.G, self.pos, edge_color=edge_colors,
                               width=2, ax=ax, arrowstyle='-|>', arrowsize=22)
        nx.draw_networkx_edge_labels(self.G, self.pos,
                                     edge_labels=self.edge_labels,
                                     font_color="darkred", ax=ax)

        legend = [mpatches.Patch(color=c, label=l) for c, l in
                  [("red", "Negative Energy"),
                   ("blue", "Positive Energy"),
                   ("black", "Negative Cycle")]]
        ax.legend(handles=legend, loc="best", fontsize=12, frameon=False)
        ax.set_title(f"{self.name} (ΔG in kcal/mol)", fontsize=18)

        plt.savefig(filename, dpi=300)
        plt.close()

    # Step 5: Compute and visualize lowest-energy path
    def lowest_energy_path(self, source, target, plot_file="lowest_energy_path.png"):
        distance, predecessor, neg_cycle = self.bellman_ford(source)

        if neg_cycle:
            print("Negative cycle detected! Cannot compute guaranteed lowest-energy path.")
            return

        # Reconstruct path
        path = self.reconstruct_path(predecessor, source, target)
        if not path:
            print(f"No path exists from {source} to {target}")
            return

        print(f"Lowest-energy path from {source} to {target}: {path}, ΔG={distance[target]:.2f}")

        edge_path = list(zip(path[:-1], path[1:]))

        fig, ax = plt.subplots(figsize=(8, 6))
        edge_colors = ["red" if w < 0 else "blue"
                       for _, _, w in self.G.edges(data="weight")]

        nx.draw(self.G, pos=self.pos, node_color='lightblue', node_size=1200,
                ax=ax, with_labels=True)
        nx.draw_networkx_edges(self.G, self.pos, edge_color=edge_colors,
                               width=2, ax=ax)
        nx.draw_networkx_edges(self.G, self.pos, edgelist=edge_path,
                               edge_color='green', width=3, ax=ax)
        nx.draw_networkx_edge_labels(self.G, self.pos,
                                     edge_labels=self.edge_labels,
                                     font_color='darkred', ax=ax)

        legend = [mpatches.Patch(color=c, label=l) for c, l in
                  [("green", "Lowest-Energy Path"),
                   ("red", "Negative ΔG"),
                   ("blue", "Positive ΔG")]]
        ax.legend(handles=legend, loc="best", fontsize=12, frameon=False)
        ax.set_title(f"{self.name} - Lowest-Energy Pathway\nΔG = {distance[target]:.2f}",
                     fontsize=16)

        plt.savefig(plot_file, dpi=300)
        plt.close()

    # Helper: reconstruct path
    def reconstruct_path(self, predecessor, start, end):
        path = []
        while end is not None:
            path.insert(0, end)
            end = predecessor[end]
        return path if path and path[0] == start else None
    
# Step 6: Example usage
if __name__ == "__main__":
    # Example 1: Graph with a negative cycle
    neg_edges_example = [
        ('U', 'A', -2.0), ('A', 'B', -1.5), ('B', 'C', -3.0), ('C', 'F', -2.0),
        ('B', 'D', 1.0), ('D', 'E', -4.0), ('E', 'F', 2.0),
        ('C', 'A', -2.5)  # negative cycle edge
    ]
    example_n_graph = ProteinGraph(edges=neg_edges_example, name="Example 6-State Negative Cycle Graph")
    example_n_graph.visualize("example_neg_cycle_graph.png")
    example_n_graph.lowest_energy_path("U", "F", "example_neg_lowest_energy_path.png")

    # Example 2: Graph without a negative cycle
    edges_example = [
        ('U', 'A', -2.0), ('A', 'B', -1.5), ('B', 'C', -3.0), ('C', 'F', -2.0),
        ('B', 'D', 1.0), ('D', 'E', -4.0), ('E', 'F', 2.0)
    ]
    example_graph = ProteinGraph(edges=edges_example, name="Example 6-State Graph")
    example_graph.visualize("example_graph.png")
    example_graph.lowest_energy_path("U", "F", "example_lowest_energy_path.png")

    # Example 3: Protein from PDB (Ubiquitin)
    ubiquitin_graph = ProteinGraph(pdb_file="1ubq.pdb", name="Ubiquitin")
    ubiquitin_graph.visualize("ubiquitin_graph.png")
    residues = sorted(ubiquitin_graph.G.nodes)
    ubiquitin_graph.lowest_energy_path(residues[0], residues[-1], "ubiquitin_lowest_energy_path.png")
# Protein Folding Graph Analysis

This project models protein folding pathways as **directed, weighted graphs**, where nodes represent folding states (or residues) and edges represent Gibbs free energy changes (ΔG in kcal/mol).

It supports two main input types:

- **Custom graphs:** defined by lists of edges with ΔG weights
- **Protein graphs:** automatically generated from PDB files using only Cα atoms

The goal is to find **lowest-energy folding pathways** and detect **negative-energy cycles** using the Bellman–Ford algorithm.

---

## Features

- Parses PDB files for residue sequence graphs
- Custom input of folding graphs with arbitrary states
- Uses Bellman–Ford algorithm to find minimum-energy paths
- Detects and highlights negative cycles, if present
- Visualizes graphs with color-coded edges and paths

---

## Graph Representation

- Nodes represent folding states or residue positions
- Edges represent ΔG energy changes along transitions
- Protein graphs assign a uniform ΔG = –0.5 kcal/mol between sequential residues

---

## Algorithm

The Bellman–Ford algorithm is used to:

- Compute the lowest-energy paths from a source node
- Detect negative-energy cycles

If a negative cycle is detected, the program:

- Prints a warning:  
  `Negative cycle detected! Cannot compute guaranteed lowest-energy path.`
- Identifies and highlights the negative cycle edges in black in the visualization
- Still outputs the best path found before the cycle detection

---

## Visualization

Edges are color-coded as follows:

- **Red:** Negative ΔG (favorable transitions)
- **Blue:** Positive ΔG (unfavorable transitions)
- **Black:** Edges part of a detected negative cycle
- **Green:** The lowest-energy path (if no negative cycle exists)

Nodes are shown as light-blue circles.

---

## Example Output

### Example 1: 6-State Graph With a Negative Cycle
When running the program with arbitrary states containing a negative cycle, we get the following:
```
Negative cycle detected! Cannot compute guaranteed lowest-energy path.
```
**Lowest-energy folding pathway:**
```
U -> A -> B -> C -> F, ΔG=-8.50
```
This means a negative cycle exists, preventing a guaranteed minimal path calculation. However, we have the algorithm still report the best path found for visualization purposes. 

*Graph visualization highlights the negative cycle edges in black.*

---

### Example 2: 6-State Graph without Negative Cycle
When running the program with arbitrary states that *do not* include a negative cycle, we get the following:

**Minimum energy distances (ΔG):**
```
U -> U: 0.00 kcal/mol
U -> A: -2.00 kcal/mol
U -> B: -3.50 kcal/mol
U -> C: -6.50 kcal/mol
U -> D: -2.50 kcal/mol
U -> E: -6.50 kcal/mol
U -> F: -4.50 kcal/mol
```

**Lowest-energy folding pathway:**
```
U -> A -> B -> C -> F, ΔG=-8.50
```

*Graph visualization highlights the lowest-energy path in green. Negative ΔG edges appear red; positive ΔG edges appear blue.*

### Example 3: Ubiquitin (PDB ID: 1UBQ)
When running the program with data using a real-life protein (ubiquitin), we get the following:

**Minimum energy distances (ΔG):**  
*(ΔG values are illustrative; actual per-residue ΔG depends on Cα sequential edges)*  
```
1 -> 1: 0.00 kcal/mol
1 -> 2: -0.50 kcal/mol
1 -> 3: -1.00 kcal/mol
...
1 -> 76: -37.50 kcal/mol
```

**Lowest-energy folding pathway:**  

```
1 -> 2 -> 3 -> ... -> 76
```

*Graph visualization highlights the lowest-energy path in green. Negative ΔG edges appear red; positive ΔG edges appear blue.*

---

## Approach & Methodology

1. **Protein Structure Data**  
   - Extracted from the Protein Data Bank (PDB), e.g., ubiquitin (PDB ID: 1UBQ)  
   - Only alpha-carbon (Cα) atoms are used for a simplified backbone graph  

2. **Graph Modeling**  
   - Nodes = residue or folding state  
   - Edges = ΔG energy change between states  
   - Negative ΔG = favorable folding, Positive ΔG = unfavorable  

3. **Energy Pathway Analysis**  
   - Bellman-Ford algorithm identifies lowest-energy paths  
   - Flags negative-energy cycles (indicative of unphysical folding)

---

## Tools & Libraries

- **Python**
  - [NetworkX](https://networkx.org/) – Graph structure
  - [Matplotlib](https://matplotlib.org/) – Visualization
- **BioPython** – PDBParser to extract residue information
- **Graphviz** – Optional hierarchical graph layout

---

## References

- Berman et al., 2000, *Nucleic Acids Research*, 28, 235–242  
- Cock et al., 2009, *Bioinformatics*, 25, 1422–1423  
- Cormen et al., 2009, *Introduction to Algorithms*, MIT Press

---

## Usage

```bash
# Create environment, install dependencies, and run script
make

# Clean generated files and virtual environment
make clean
```
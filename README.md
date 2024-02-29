# PerturboNetKit Tutorial

## Version: 1.0
## Date: 2024-02-13
## Author: Stephan Buchner

This tutorial provides a step-by-step guide on using the Reusablenetworkanalysis package for network analysis. The package is used to analyze biological networks, and it includes methods for calculating various network metrics and visualizing the results.
Moreover, the package provides code for the analyses of perturbomes and calculation of interactions between between nodes.
### Installation

Before getting started, make sure you have the Reusablenetworkanalysis package installed. You can install it using pip:
not available yet - need to publish on Pypi first
```bash
pip install PerturboNetKit
```

Other possible installation procedure:
```bash
git clone https://github.com/Bunibal/Softwaredevelopmentinternship
cd PerturboNetKit
python setup.py install
```

### Tutorial

```python
# Import necessary libraries
import PerturboNetKit as RNA
import networkx as nx
import pandas as pd
import numpy as np

# Step 1: Load the Human interactome network
PPI = nx.read_gml("path/to/Human_Interactome.gml")

# Step 2: Subsample the network for shorter runtime
PPI = PPI.subgraph(list(PPI.nodes)[0:5000])

# Step 3: Initialize PerturboNetKit analysis
analysis = RNA.Networkpipeline(PPI,
                               targets_node_file="path/to/CLOUD_All_Targets.csv",
                               verbose=True)

# Step 4: Analyze network properties

# Uncomment the following lines to get degree-related information
# degrees, degree_distribution, cumulative_degree_distribution = analysis.get_degree_distribution()
# analysis.plot_degree_distributions()

# Obtain a subnetwork of only the specified targets
subnetwork = analysis.get_subnetwork_nodes(overwrite_network=True)
# Potentially plot the subnetwork
analysis.plot_network(node_size=10, only_LCC=True)

# Calculate the closeness-centrality of the nodes
analysis.calculate_centrality_node()
# Potentially plot the centrality
analysis.plot_centralities()

# Compare the size of the largest connected component with random nodes
analysis.compare_lcc_size_against_random()
# Potentially plot the results
analysis.plot_lcc_size_results()
# Or view the results in the terminal
for key in analysis.node_lcc_results:
    print(key, analysis.node_lcc_results[key])

# Calculate the shortest path between the targets
analysis.check_shortest_path_between_targets_against_random()
# Potentially plot the results
analysis.plot_shortest_path_between_targets()

# Calculate the overlap between the nodes
analysis.calculate_distances_between_clouds()
# View the results in the terminal
for key in analysis.overlap_between_nodes:
    print(key, analysis.overlap_between_nodes[key])

# Step 5: Interaction between nodes

# Import necessary classes
from PerturboNetKit import Perturbome, CalculateInteractions

# Read perturbation data from CSV
perts = pd.read_csv("path/to/test_data_combi_seq_paper.csv",
                    sep=";", header=0, index_col=0)

# Create a dictionary of perturbations and interactions
perturbations = {}
interactions = {}

# Create the Perturbome object
perturbome = Perturbome(perturbations=perturbations, interactions=interactions)

# Extract no_treatment_samples as pd.DataFrame
no_treatment_samples = pd.DataFrame(perts.loc[:, "no_treatment"])

# Calculate interactions between the nodes
ints = CalculateInteractions(perturbome, no_treatment_samples)
# Get interaction values for all perturbations
ints.get_interaction_values()
# Categorize the interactions
ints.categorize_interactions()
# Potentially plot the number of different interactions
ints.plot_interactions_histogram()
```

For further information on the methods and their parameters, please refer to the docstrings of the respective methods or the api.html file in the docs folder.


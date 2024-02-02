# Reusablenetworkanalysis Tutorial

## Version: 0.1
## Date: 2023-12-29
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
import PerturboNetKit as PNK
import networkx as nx

# Load the Human interactome network
PPI = nx.read_gml("/path/to/Human_Interactome.gml")

# Possibly subsample the network for shorter runtime
PPI = PPI.subgraph(list(PPI.nodes)[0:5000])

# Initialize the Networkpipeline with the Human interactome network
analysis = PNK.Networkpipeline(
    PPI,
    targets_node_file="/path/to/CLOUD_All_Targets.csv",
    verbose=True
)

# Calculate the network metrics available in the Analysisclass. 
# The code was reframed and updated from the original code by Caldera et al. (2019).
degrees, degree_distribution, cumulative_degree_distribution = analysis.get_degree_distribution()
analysis.plot_degree_distributions()
subnetwork = analysis.get_subnetwork_nodes(overwrite_network=True)
analysis.calculate_centrality_node()
analysis.compare_lcc_size_against_random()
analysis.check_shortest_path_between_targets_against_random()
analysis.calculate_overlap_between_nodes()

# Results visualization methods are available in the Plottingclass
analysis.plot_network(node_size=10, only_LCC=True)
analysis.plot_degree_distribution()
analysis.plot_centralities()
analysis.plot_lcc_size_results()
analysis.plot_shortest_path_between_targets()

# Calculate distances between nodes. 
# Code was updated from the original code by Guney et al. (2016).

# TODO: Add example for distance calculation

# Calculate the perturbome of the network
from PerturboNetKit import Perturbome, CalculateInteractions

# Initialize the Perturbome class / load your perturbome
perturbome = Perturbome(
)
# TODO: Add example for perturbome/interaction calculation
```

For further information on the methods and their parameters, please refer to the docstrings of the respective methods or the api.html file in the build folder.

### References

Caldera, M. et al. (2019) ‘Mapping the perturbome network of cellular perturbations’, Nature Communications, 10, p. 5140. Available at: https://doi.org/10.1038/s41467-019-13058-9.

Guney, E. et al. (2016) ‘Network-based in silico drug efficacy screening’, Nature Communications, 7(1), p. 10331. Available at: https://doi.org/10.1038/ncomms10331.

Menche, J. et al. (2015) ‘Uncovering disease-disease relationships through the incomplete interactome’, Science, 347(6224), p. 1257601. Available at: https://doi.org/10.1126/science.1257601.

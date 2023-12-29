# Version: 1.0
# Date: 2023-12-29
# Author: Stephan Buchner

import Reusablenetworkanalysis as RNA
import networkx as nx

# load the Human interactome network

PPI = nx.read_gml("/home/bunibal/PycharmProjects/Softwaredevelopmentinternship/data/Human_Interactome.gml")

# subsample the network for shorter runtime
PPI = PPI.subgraph(list(PPI.nodes)[0:5000])

analysis = RNA.Networkpipeline(PPI,
                               targets_node_file="/home/bunibal/PycharmProjects/Softwaredevelopmentinternship/data/CLOUD_All_Targets.csv",
                               verbose=True)

# All of the available Analysis methods are listed below. Uncomment the ones you do not want to use.
degrees, degree_distribution, cumulative_degree_distribution = analysis.get_degree_distribution()
subnetwork = analysis.get_subnetwork_nodes(overwrite_network=True)
analysis.calculate_centrality_node()
analysis.compare_lcc_size_against_random()
analysis.check_shortest_path_between_targets_against_random()
analysis.calculate_overlap_between_nodes()

# Plottings of the results are available in the Plottingclass. Uncomment the ones you do not want to use.

analysis.plot_network(node_size=10, only_LCC=True)
analysis.plot_degree_distribution()

# Version: 1.0
# Date: 2023-12-29
# Author: Stephan Buchner

import PerturboNetKit as RNA
import networkx as nx

# load the Human interactome network

PPI = nx.read_gml("/home/bunibal/PycharmProjects/Softwaredevelopmentinternship/data/Human_Interactome.gml")

# subsample the network for shorter runtime
PPI = PPI.subgraph(list(PPI.nodes)[0:5000])

analysis = RNA.Networkpipeline(PPI,
                               targets_node_file="/home/bunibal/PycharmProjects/Softwaredevelopmentinternship/data/CLOUD_All_Targets.csv",
                               verbose=True)

# All the available Analysis methods are listed below. Uncomment the ones you do not want to use.
degrees, degree_distribution, cumulative_degree_distribution = analysis.get_degree_distribution()
# if you do not provide networks as input you should ue the following plotting method before getting the subnetwork
analysis.plot_degree_distributions()

# obtain a subnetwork of only the targets provieded in node_to_target dictionary
subnetwork = analysis.get_subnetwork_nodes(overwrite_network=True)
# potentially plot the subnetwork
analysis.plot_network(node_size=10, only_LCC=True)

# calculate the closeness-centrality of the nodes
analysis.calculate_centrality_node()
# potentially plot the centrality
analysis.plot_centralities()

# get the largest connected component and compare it to the lcc of randomly drawn nodes
analysis.compare_lcc_size_against_random()
# potentially plot the results
analysis.plot_lcc_size_results()
# or view the results in the terminal
for key in analysis.node_lcc_results:
    print(key, analysis.node_lcc_results[key])

# calculate the shortest path between the targets
analysis.check_shortest_path_between_targets_against_random()
# potentially plot the results
analysis.plot_shortest_path_between_targets()

# calculate the overlap between the nodes
analysis.calculate_overlap_between_nodes()
# view the results in the terminal
for key in analysis.overlap_between_nodes:
    print(key, analysis.overlap_between_nodes[key])

# Interaction between nodes
from PerturboNetKit import Perturbome, CalculateInteractions
import pandas as pd
import numpy as np

perts = pd.read_csv("/home/bunibal/PycharmProjects/Softwaredevelopmentinternship/data/test_data_combi_seq_paper.csv",
                    sep=";", header=0, index_col=0)

#replace the commas with dots
perts = perts.apply(lambda x: x.str.replace(",", "."))

# get the perturbations into correct format for the perturbome object
ind = lambda: zip([drug.split("_")[3] for drug in perts.index],
          [drug.split("_")[4].removesuffix("ReadsPerGene") for drug in perts.index])

perturbations = {}
interactions = {}
for i, (a, b) in enumerate(ind()):
    if a == "DMSO" and b != "DMSO":
        perturbations[b] = perts.iloc[i].to_numpy(dtype=np.float64)
    elif b == "DMSO" and a != "DMSO":
        perturbations[a] = perts.iloc[i].to_numpy(dtype=np.float64)
    elif a != "DMSO" and b != "DMSO" and a != b:
        interactions[(a, b)] = perts.iloc[i].to_numpy(dtype=np.float64)

# create the perturbome object
perturbome = Perturbome(perturbations=perturbations, interactions=interactions)

# no_treatment_samples

no_treatment_samples = perts[[a == "DMSO" and b == "DMSO" for a, b in ind()]].apply(pd.to_numeric, errors="coerce")


# calculate the interactions between the nodes
ints = CalculateInteractions(perturbome, no_treatment_samples)
# get the interaction values for all the perturbations
ints.get_interaction_values()
# categorize the interactions
ints.categorize_interactions()
# potentially plot the # number of different interactions
ints.plot_interactions_histogram()

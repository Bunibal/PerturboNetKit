import Reusablenetworkanalysis as RNA
import networkx as nx

# load the Human interactome network

PPI = nx.read_gml("/home/bunibal/PycharmProjects/Softwaredevelopmentinternship/data/Human_Interactome.gml")

# subsample the network for shorter runtime
PPI = PPI.subgraph(list(PPI.nodes)[0:5000])
analysis = RNA.Analysis(PPI,
                        targets_node_file="/home/bunibal/PycharmProjects/Softwaredevelopmentinternship/data/CLOUD_All_Targets.csv",
                        verbose=True)

degrees, degree_distribution, cumulative_degree_distribution = analysis.get_degree_distribution()
subnetwork = analysis.get_subnetwork_nodes(overwrite_network=True)
# analysis.calculate_centrality_node()
analysis.compare_lcc_size_against_random()


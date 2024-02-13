import pandas as pd
import os
import numpy as np
import networkx as nx
import random as rd
from scipy.stats import mannwhitneyu
from .Distancesclass import NodeDistances

flatten = lambda *n: (e for a in n
                      for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))
"""
Flattens a list and all sublists recursively.

:param *n: unflattened list
:type *n: list

:return: flattened list
:rtype: list
"""


class Analysis(NodeDistances):
    """
    This class is used to perform network analysis on a given network.
    It comprises the following steps of analysis:

    - Degree distribution
    - Centrality (closeness-centrality)
    - LCC size
    - shortest path

    Attributes:
    -----------
    network: networkx.Graph
        The network which should be analyzed

    node_to_target: dict
        dictionary mapping targets to their node (disease,drug, ....)

    disease_similarties_bin_means : None
        Placeholder for disease similarity bin means.

    similarties_bin_means : None
        Placeholder for similarity bin means.

    overlap_between_nodes : None
        Placeholder for node overlap results.

    cloud_ShortestPaths_results : None
        Placeholder for cloud shortest paths results.

    node_lcc_results : None
        Placeholder for node largest connected component results.

    node_centrality_results : None
        Placeholder for node centrality results.

    targets_node_file : str or None
        A CSV file containing the targets and their corresponding nodes.

    verbose : bool
        A boolean indicating if the class should print information.
    """

    def __init__(self, network, targets_node_file=None, node_to_target=None, verbose=False):  # TODO: specify how csv file has to look like with example
        """
        Class constructor initializing attributes and calling the superclass constructor.

        :param network: The network to be analyzed in NetworkX format.
        :type network: networkx.Graph

        :param targets_node_file: A CSV file containing the targets and their corresponding nodes.
        :type targets_node_file: str or None

        :param node_to_target: Mapping of nodes to their corresponding targets.
        :type node_to_target: dict or None

        :param verbose: A boolean indicating if the class should print information.
        :type verbose: bool
        """
        self.disease_similarties_bin_means = None
        self.similarties_bin_means = None
        self.overlap_between_nodes = None
        self.cloud_ShortestPaths_results = None
        self.node_lcc_results = None
        self.node_centrality_results = None
        self.targets_node_file = targets_node_file  # A csv file containing the targets and their corresponding nodes
        self.verbose = verbose  # A boolean indicating if the class should print information
        super().__init__(network, node_to_target)
        self.check_input()  # Check if the input is correct

    # function to check the input
    def check_input(self):
        """
        This function checks if the network is a networkx graph or not.
        It also checks if files exist that you have given the class and loads the targets_node_file if given.
        """
        if not isinstance(self.network, nx.Graph):
            raise TypeError('The network is not a networkx graph.')
        if self.targets_node_file is not None:
            if not os.path.exists(self.targets_node_file):
                raise FileNotFoundError('The provided targets node file does not exist.')
            self.get_targets_node()
        else:
            print("Please provide a targets_node_file or add it from another source.")

    def get_targets_node(self, targets_node_file=None, columnname="EntrezIDs"):
        """
        This function reads a csv file containing the targets and their corresponding nodes and returns a dictionary.
        Alternatively you can save a dictionary of the map into the attribute node_to_target of the Analysisobject or Networkpipeline.

        :param targets_node_file: You can provide a csv file to read the map from with the given columnname (default: EntrezIDs). Should match the nodenames in the network
        :type targets_node_file: file.csv

        :return: a map connecting the node (key) to the targets (values)
        :rtype: dict
        """
        if targets_node_file is None and self.targets_node_file is None:
            print("You need to provide a filename to read from!")
        elif targets_node_file is not None:
            self.targets_node_file = targets_node_file
        self.node_to_target = {}
        df = pd.read_csv(self.targets_node_file, sep=',', header=0, index_col=0)
        for index, row in df.iterrows():
            if type(row[columnname]) == str:
                self.node_to_target[index] = row[columnname].split(';')
            else:
                print(f"{row[columnname]} with the type {type(row[columnname])} cannot be read.")
        if self.verbose:
            print('Number of different targets: %d' % len(self.node_to_target))
        return self.node_to_target

    def get_degree_distribution(self, network=None):
        """
        Return the degree distribution of the network.

        :param network: The network to be analyzed in NetworkX format.
        :type network: networkx.Graph

        :return: A tuple containing three lists:
                 - degrees (list): List of degrees for each node in the network.
                 - degree_distribution (list): List of frequencies of each degree in the network.
                 - cumulative_degree_distribution (list): List of cumulative frequencies of degrees.
        :rtype: tuple
        """
        if network is None:
            network = self.network
        degrees_network = [x[1] for x in network.degree()]
        degrees_network_unique = list(set(degrees_network))
        degrees_network_unique.sort()
        degrees = []
        degree_distribution = []
        cumulative_degree_distribution = []
        for degree in degrees_network_unique:
            degrees.append(degree)
            degree_distribution.append(degrees_network.count(degree) / float(len(degrees_network)))
            cumulative_degree_distribution.append(
                len([x for x in degrees_network if x >= degree]) / float(len(degrees_network)))
        if self.verbose:
            print(f"Mean network degree: {round(np.mean(degrees_network))}")
        return degrees, degree_distribution, cumulative_degree_distribution

    def get_subnetwork_nodes(self, overwrite_network=False):
        """
        Return a subnetwork containing only the targets.

        :param overwrite_network: Indicates whether you want to replace the network in the analysis with the subnetwork
        :type overwrite_network: networkx.Graph

        :return: A NetworkX graph representing the subnetwork containing only targets.
        :rtype: networkx.Graph
        """
        if self.node_to_target is None:
            raise ValueError('No nodes to target dictionary provided.')
        targets = list(set(flatten(list(self.node_to_target.values()))))
        subnetwork_node = self.network.subgraph(targets)
        if overwrite_network:
            self.network = subnetwork_node
        return subnetwork_node

    def calculate_centrality(self, targets, network=None):
        """
        Calculate closeness centralities for nodes in the network.

        :param targets: List of targets.
        :type targets: list

        :param network: NetworkX graph.
        :type network: networkx.Graph

        :return: List of centralities, or None if no node is found in the network.
        :rtype: list or None
        """
        if network is None:
            network = self.network
        filtered_targets = self.filter_targets_in_network(network=network, targets=targets)

        if len(filtered_targets) > 0:
            centralities = []
            for node in filtered_targets:
                centralities.append(nx.closeness_centrality(network, node))

            return centralities
        else:
            return None

    def calculate_centrality_node(self, node_to_target=None, network=None, NumRandom=10000):
        """
        Calculate closeness centralities for nodes in the network.

        :param node_to_target: Map of the wanted nodes in which the centrality should be calculated. If None the all nodes are used.
        :type node_to_target: dict

        :param network: The network which should be used.
        :type network: networkx.Graph

        :param NumRandom: Number of random nodes to be sampled.
        :type NumRandom: int

        :return: Dictionary containing the centralities for each node.
        :rtype: dict
        """
        if network is None:
            network = self.network

        if node_to_target is None:
            node_to_target = self.node_to_target

        # create a random distribution of centralities on the network
        if self.verbose:
            print('Create Random Distribution')

        random_distances = []
        for i in range(0, NumRandom):
            node = rd.sample(network.nodes(), 1)[0]
            random_distances.append(nx.closeness_centrality(network, node))

        # calculate the centralities of the targets for an individual drug
        node_centrality_results = {}
        for c in node_to_target.keys():
            if len(node_to_target[c]) == 0:
                continue

            # get the centralities for all targets for a specific node
            node_centralities = self.calculate_centrality(node_to_target[c], network=network)

            # calculate the PValue,MeanCentrality and Glass' Delta if targets found in the network
            if node_centralities is not None:
                p_value = mannwhitneyu(node_centralities, random_distances)[1]
                # FoldChange = d_d / np.mean(random_distances)
                glass_delta = (np.mean(node_centralities) - np.mean(random_distances)) / np.std(random_distances)

                node_centrality_results[c] = {'MeanCentrality': np.mean(node_centralities), 'PValue': p_value,
                                              'FoldChange': glass_delta}

        # Save the results
        self.node_centrality_results = node_centrality_results
        return node_centrality_results

    def get_lcc_size(self, targets, network=None):
        """
        Return the size of the largest connected component (LCC) of given targets in the network.

        :param network: The network which should be used.
        :type network: networkx.Graph

        :param targets: Targets to be analyzed.
        :type targets: list

        :return: The size of the LCC.
        :rtype: int
        """
        if network is None:
            network = self.network

        sub_graph = nx.subgraph(network, targets)
        if len(sub_graph.nodes) > 0:
            components = (sub_graph.subgraph(c) for c in nx.connected_components(sub_graph))
            lcc = max([len(x.nodes()) for x in components])

        return lcc

    def compare_lcc_size_against_random(self, network=None, targets=None, RandomIterations=10000):
        """
        Compare the actual LCC size of a node against LCCs from randomly drawn nodes.

        Randomly drawn nodes have the same size as the number of targets in the actual node.

        :param network: A NetworkX graph.
        :type network: networkx.Graph

        :param targets: Dictionary containing the nodes to be analyzed and their targets.
        :type targets: dict

        :param RandomIterations: Number of random iterations to be performed. Default: 10000.
        :type RandomIterations: int

        :return: Dictionary containing the results (LCC, Zscore, #Genes, RelativeSize).
        :rtype: dict
        """
        if network is None:
            network = self.network
        if targets is None:
            targets = self.node_to_target

        # create a random distribution of centralities on the network
        if self.verbose:
            print('Create Random nodes...')

        node_lccs_results = {}
        for c in targets.keys():
            # Get specific drug LCC
            if len(targets[c]) == 0:
                continue
            lcc = self.get_lcc_size(targets[c], network=network)
            number_of_genes = len(targets[c])
            relative_lcc_size = float(lcc) / number_of_genes

            # Create random distribution of LCCs
            random_distribution = []
            for i in range(0, RandomIterations):
                rand_genes = rd.sample(network.nodes(), len(targets[c]))
                rand_lcc = self.get_lcc_size(rand_genes, network=network)
                random_distribution.append(rand_lcc)

            # Calculate z_score
            z_score = (lcc - np.mean(random_distribution)) / np.std(random_distribution)

            node_lccs_results[c] = {'lcc': lcc, 'Zscore': z_score, '#Genes': number_of_genes,
                                    'RelativeSize': relative_lcc_size}
        # Save the results
        self.node_lcc_results = node_lccs_results
        return node_lccs_results

    def check_shortest_path_between_targets_against_random(self, network=None, targets=None, RandomIterations=10000):
        """
        Calculate the intra-drug module distance, i.e., shortest distances between any targets of a given drug.

        The result is the mean over all such distances.

        :param network: A NetworkX graph.
        :type network: networkx.Graph

        :param targets: Dictionary containing the nodes to be analyzed and their targets.
        :type targets: dict

        :param RandomIterations: Number of random iterations to be performed. Default: 10000.
        :type RandomIterations: int

        :return: Dictionary containing the results.
        :rtype: dict
        """
        if network is None:
            network = self.network
        if targets is None:
            targets = self.node_to_target

        # create a random distribution of centralities on the network
        if self.verbose:
            print('Create Random targets...')

        # Create random distribution for randomly drawn proteins on the PPI
        random_distances = []
        for i in range(0, RandomIterations):
            targets_rand = rd.sample(network.nodes(), 2)
            if nx.has_path(network, targets_rand[0], targets_rand[1]):
                random_distances.append(len(nx.shortest_path(network, targets_rand[0], targets_rand[1])) - 1)

        # Calculate the module network diameter
        if self.verbose:
            print('Calculate Shortest Paths...')
        cloud_ShortestPaths_results = {}
        for c in targets.keys():
            if self.verbose:
                print(c)
            if len(targets[c]) == 0:
                continue

            # Extract min distances
            d_d, min_paths = self.get_shortest_mean_distance(targets[c], network=network)

            if d_d is None:
                continue

            # Calculate p_value by comparing with random Distribution
            p_value = mannwhitneyu(min_paths, random_distances)[1] # TODO: change to t-test

            # Calculate fold change (Glass' Delta)
            glass_delta = (d_d - np.mean(random_distances)) / np.std(random_distances)

            # Save Result
            cloud_ShortestPaths_results[c] = {'MeanShortestPath': d_d, 'PValue': p_value, 'FoldChange': glass_delta}

        # Save the results
        self.cloud_ShortestPaths_results = cloud_ShortestPaths_results
        return cloud_ShortestPaths_results

    def calculate_overlap_between_nodes(self, nodes=None):
        """
        Calculate the overlap between nodes.

        :param nodes: Dictionary containing the nodes to be analyzed and their targets.
        :type nodes: dict

        :return: Dictionary containing the results.
        :rtype: dict
        """

        if nodes is None:
            nodes = self.node_to_target
        within_Distances = {}
        for c in nodes:
            # print c
            if len(nodes[c]) == 0:
                continue

            d_d, min_paths = self.get_shortest_mean_distance(network=self.network,
                                                                                  targets=nodes[c])

            if d_d == None:
                continue
            else:
                within_Distances[c] = d_d

        df = pd.DataFrame(columns=["Node1", "Node2", "D_Node1", "D_Node2", "D_D1_D2", "S"])
        clouds = within_Distances.keys()
        for c in clouds:
            d_A = within_Distances[c]
            targets1 = nodes[c]
            for c2 in clouds:
                d_B = within_Distances[c2]
                targets2 = nodes[c2]
                distances = self.get_pathlengths_for_two_sets(targets1, targets2, network=self.network)
                distances1 = [min(d.values()) for k, d in distances.items() if k in targets1]
                distances2 = [min(d.values()) for k, d in distances.items() if k in targets2]
                # Dab
                if type(distances2[0]) == int and type(distances1[0]) == int:
                    between_Distance = (sum(distances1) + sum(distances2)) / float((len(distances1) + len(distances2)))

                    # Sab
                    separation = between_Distance - (d_A + d_B) / 2.0

                    df.loc[len(df)] = ([c, c2, d_A, d_B, between_Distance,
                               separation])
        self.overlap_between_nodes = df
        return df

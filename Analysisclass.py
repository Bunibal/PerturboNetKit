import pandas as pd
import os
import numpy as np
import logging  # TODO: Ask mathilde if this is necessary
import networkx as nx
import random as rd
from scipy.stats import mannwhitneyu
from Distancesclass import HubDistances

flatten = lambda *n: (e for a in n
                      for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))


class Analysis(HubDistances):
    """
    This class is used to perform network analysis on a given network.

    """

    def __init__(self, network, targets_hub_file=None, hub_to_target=None, verbose=False):
        self.cloud_ShortestPaths_results = None
        self.hub_lcc_results = None
        self.hub_centrality_results = None
        self.targets_hub_file = targets_hub_file  # A csv file containing the targets and their corresponding hubs
        self.verbose = verbose  # A boolean indicating if the class should print information
        super().__init__(network, hub_to_target)
        self.check_input()  # Check if the input is correct

    # function to check the input
    def check_input(self):
        """
        This function checks if the network is a networkx graph or not.
        It also checks if files exist that you have given the class.
        """
        if not isinstance(self.network, nx.Graph):
            raise TypeError('The network is not a networkx graph.')
        if self.targets_hub_file is not None:
            if not os.path.exists(self.targets_hub_file):
                raise FileNotFoundError('The provided targets hub file does not exist.')

    def get_targets_hub(self):
        """
        This function reads a csv file containing the targets and their corresponding hubs and returns a dictionary
        :return: dict
        """
        self.hub_to_target = {}
        df = pd.read_csv(self.targets_hub_file, sep=',', header=True, index_col=0)
        for index, row in df.iterrows():
            self.hub_to_target[index] = row.split(';')
        if self.verbose:
            print('Number of different targets: %d' % len(self.hub_to_target))
        return self.hub_to_target

    def get_degree_distribution(self, network=None):
        """
        This function returns the degree distribution of the network.
        :param network: The network to be analyzed in networkx format
        :return: degrees(list), degree_distribution(list), cumulative_degree_distribution(list)
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

    def get_subnetwork_hub(self):
        """
        This function returns a subnetwork containing only the targets.
        :return:
        """
        if self.hub_to_target is None:
            raise ValueError('No hub to target dictionary provided.')
        targets = list(set(flatten(list(self.hub_to_target.values()))))
        subnetwork_hub = self.network.subgraph(targets)
        return subnetwork_hub

    def calculate_centrality(self, targets, network=None):
        """
        Calculates closeness centralities for nodes in the network.
        :param targets: list of targets
        :param network: networkx graph
        :return list of centralities or None if no node is found in the network.
        """
        if network is None:
            network = self.network
        filtered_targets = self.filter_targets_in_network(network, targets)

        if len(filtered_targets) > 0:
            centralities = []
            for node in filtered_targets:
                centralities.append(nx.closeness_centrality(network, node))

            return centralities
        else:
            return None

    def calculate_centrality_hub(self, targets=None, network=None, NumRandom=10000):
        """
        Calculates closeness centralities for nodes in the network.
        :param targets: list of targets
        :param network: networkx graph
        :param NumRandom: Number of random nodes to be sampled
        :return: dictionary containing the centralities for each hub
        """
        if network is None:
            network = self.network

        if targets is None:
            targets = self.hub_to_target

        # create a random distribution of centralities on the network
        if self.verbose:
            print('Create Random Distribution')

        random_distances = []
        for i in range(0, NumRandom):
            node = rd.sample(network.nodes(), 1)[0]
            random_distances.append(nx.closeness_centrality(network, node))

        # calculate the centralities of the targets for an individual drug
        hub_centrality_results = {}
        for c in targets.keys():
            print(c)
            if len(targets[c]) == 0:
                continue

            # get the centralities for all targets for a specific hub
            hub_centralities = self.calculate_centrality(network, targets[c])

            # calculate the PValue,MeanCentrality and Glass' Delta if targets found in the network
            if hub_centralities is not None:
                p_value = mannwhitneyu(hub_centralities, random_distances)[1]
                # FoldChange = d_d / np.mean(random_distances)
                glass_delta = (np.mean(hub_centralities) - np.mean(random_distances)) / np.std(random_distances)

                hub_centrality_results[c] = {'MeanCentrality': np.mean(hub_centralities), 'PValue': p_value,
                                             'FoldChange': glass_delta}

        # Save the results
        self.hub_centrality_results = hub_centrality_results
        return hub_centrality_results

    def get_lcc_size(self, targets, network=None):
        """
        This function returns the lcc size of given targets in the network.
        :param network: a networkx graph
        :param targets: targets to be analyzed
        :return: the size of the lcc
        """
        if network is None:
            network = self.network

        lcc = 1
        sub_graph = nx.subgraph(network, targets)
        if len(sub_graph.nodes) > 0:
            components = (sub_graph.subgraph(c) for c in nx.connected_components(sub_graph))
            lcc = max([len(x.nodes()) for x in components])

        return lcc

    def compare_lcc_size_against_random(self, network=None, targets=None, RandomIterations=10000):
        """
        This function compares the actual lcc size of a hub against a lcc from randomly drawn hubs, which have the
        same size as the number of targets in the actual hub.
        :param network: a networkx graph
        :param targets: dictionary containing the hubs, to be analyzed and their targets
        :param RandomIterations: number of random iterations to be performed. Default: 10000
        :return: dictionary containing the results(lcc, Zscore, #Genes, RelativeSize)
        """
        if network is None:
            network = self.network
        if targets is None:
            targets = self.hub_to_target

        # create a random distribution of centralities on the network
        if self.verbose:
            print('Create Random hubs...')

        hub_lccs_results = {}
        for c in targets.keys():
            print(c)

            # Get specific drug LCC
            if len(targets[c]) == 0:
                continue
            lcc = self.get_lcc_size(network, targets[c])
            number_of_genes = len(targets[c])
            relative_lcc_size = float(lcc) / number_of_genes

            # Create random distribution of LCCs
            random_distribution = []
            for i in range(0, RandomIterations):
                rand_genes = rd.sample(network.nodes(), len(targets[c]))
                rand_lcc = self.get_lcc_size(network, rand_genes)
                random_distribution.append(rand_lcc)

            # Calculate z_score
            z_score = (lcc - np.mean(random_distribution)) / np.std(random_distribution)

            hub_lccs_results[c] = {'lcc': lcc, 'Zscore': z_score, '#Genes': number_of_genes,
                                   'RelativeSize': relative_lcc_size}
        # Save the results
        self.hub_lcc_results = hub_lccs_results
        return hub_lccs_results

    def get_shortest_distances(self, targets, network=None):
        """
        This function calculates the shortest distances between all targets in the network.
        This is always the minimum path between one target and any other target of the other set.
        This function uses only one set hence calculates the intra hub distance or drug_module diameter

        :param network: a networkx graph
        :param targets: targets to be analyzed
        :return: Mean of all paths (d_d) as well as paths (min_paths)
        """

        if network is None:
            network = self.network

        filtered_targets = self.filter_targets_in_network(targets, network)

        min_paths = []
        if len(filtered_targets) > 1:
            try:
                for t1 in filtered_targets:
                    min_distances = []
                    for t2 in filtered_targets:
                        if t1 != t2:

                            if nx.has_path(network, t1, t2):
                                min_distances.append(len(nx.shortest_path(network, t1, t2)) - 1)
                    min_paths.append(min(min_distances))
                d_d = sum(min_paths) / float(len(filtered_targets))

                return d_d, min_paths
            except:
                return None, None
        else:
            return None, None

    def check_shortest_path_between_targets_against_random(self, network=None, targets=None, RandomIterations=10000):
        """
        Function to calculate the intra_drug module distance e.g. shortest distances of any target of a given drug to any other target of this drug (mean over all)
        :param network: a networkx graph
        :param targets: dictionary containing the hubs, to be analyzed and their targets
        :param RandomIterations: number of random iterations to be performed. Default: 10000
        :return: dictionary containing the results
        """
        if network is None:
            network = self.network
        if targets is None:
            targets = self.hub_to_target

        # create a random distribution of centralities on the network
        if self.verbose:
            print('Create Random targets...')

        # Create random distribution for randomly drawn proteins on the PPI
        random_distances = []
        for i in range(0, RandomIterations):
            targets_rand = rd.sample(network.nodes(), 2)
            if nx.has_path(network, targets_rand[0], targets_rand[1]):
                random_distances.append(len(nx.shortest_path(network, targets[0], targets[1])) - 1)

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
            d_d, min_paths = self.get_shortest_distances(network, targets[c])

            if d_d is None:
                continue

            # Calculate p_value by comparing with random Distribution
            p_value = mannwhitneyu(min_paths, random_distances)[1]

            # Calculate fold change (Glass' Delta)
            glass_delta = (d_d - np.mean(random_distances)) / np.std(random_distances)

            # Save Result
            cloud_ShortestPaths_results[c] = {'MeanShortestPath': d_d, 'PValue': p_value, 'FoldChange': glass_delta}

        # Save the results
        self.cloud_ShortestPaths_results = cloud_ShortestPaths_results
        return cloud_ShortestPaths_results

    def calculate_overlap_between_hubs(self, hubs=None): # TODO: rewrite
        """
        This function calculates the overlap between the hubs
        :param hubs: dictionary containing the hubs, to be analyzed and their targets
        :return: dictionary containing the results
        """
        within_Distances = {}
        for c in cloud_targets:
            # print c
            if len(cloud_targets[c]) == 0:
                continue

            d_d, min_paths = Check_Shortest_Distances(PPI, cloud_targets[c])

            if d_d == None:
                continue
            else:
                within_Distances[c] = d_d

        clouds = within_Distances.keys()
        clouds = sorted(clouds)
        for c in clouds:
            print(c)
            d_A = within_Distances[c]
            targets1 = cloud_targets[c]

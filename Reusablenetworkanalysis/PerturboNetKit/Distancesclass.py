import networkx as nx
import numpy as np


class NodeDistances:
    """
    This class calculates the distances between hubs or targets in a network.

    Attributes:
    -----------
    network : networkx.Graph
        The network to calculate distances on.

    node_to_target : dict
        Dictionary containing the nodes and their corresponding targets.

    """

    def __init__(self, network, node_to_target):
        """
        Class constructor initializing attributes.

        :param network: the network which should be used.
        :type network: networkx.Graph

        :param node_to_target: Dictionary containing the nodes and their corresponding targets.
        :type node_to_target: dict
        """
        self.network = network
        self.node_to_target = node_to_target

    def filter_targets_in_network(self, network=None, targets=None):
        """
        Filter targets that are not in the network.

        :param network: the network which should be used.
        :type network: networkx.Graph

        :param targets: Dictionary containing the hubs and their corresponding targets.
        :type targets: dict or list

        :return: List of targets that are in the network.
        :rtype: list
        """

        if network is None:
            network = self.network
        if targets is None:
            targets = self.node_to_target

        filtered_targets = []
        for t in targets:
            if network.has_node(t):
                filtered_targets.append(t)
        return filtered_targets

    def get_pathlengths_for_single_set(self, network, node):

        """
       Calculate the shortest paths for a given set of genes in a network.

       The results are stored in a dictionary of dictionaries, where each pair of genes is represented as `gene1` and `gene2` with `gene1 < gene2`. This ensures that each pair is stored only once.

       :param network: The network.
       :type network: networkx.Graph

       :param node: The set of genes for which paths should be computed.
       :type node: set

       :return: A dictionary of dictionaries representing the shortest paths.
                - `all_path_lengths[gene1][gene2] = l` for all pairs of genes with `gene1 < gene2`,
                  where `l` is the shortest path length.
       :rtype: dict

        """
        given_gene_set = set(self.node_to_target[node])
        # remove all nodes that are not in the network
        all_genes_in_network = set(network.nodes())
        gene_set = given_gene_set & all_genes_in_network

        all_path_lenghts = {}

        # calculate the distance of all possible pairs
        for gene1 in gene_set:
            if gene1 not in all_path_lenghts.keys():
                all_path_lenghts[gene1] = {}
            for gene2 in gene_set:
                if gene1 < gene2:
                    try:
                        l = nx.shortest_path_length(network, source=gene1, target=gene2)
                        all_path_lenghts[gene1][gene2] = l
                    except:
                        continue

        return all_path_lenghts

    def get_pathlengths_for_two_sets(self, given_gene_set1,
                                     given_gene_set2, network=None):

        """
        Extract the minimum path between target_sets.

        This function calculates the minimum path between each target in the first set (targets1)
        and any other target in the second set (targets2) in the provided network.

        :param given_gene_set1: List of targets, first set
        :type given_gene_set1: list

        :param given_gene_set2: List of targets, second set
        :type given_gene_set2: list

        :param network: Networkx graph
            The network to be analyzed in networkx format.
        :type network: nx.Graph

        :return: A dictionary of dictionaries representing the shortest paths.
        :rtype: dict
        """
        if network is None:
            network = self.network

        # remove all nodes that are not in the network
        gene_set1 = self.filter_targets_in_network(network=network, targets=given_gene_set1)
        gene_set2 = self.filter_targets_in_network(network=network, targets=given_gene_set2)

        all_path_lenghts = {gene: {} for gene in gene_set1 + gene_set2}

        # calculate the distance of all possible pairs
        for gene1 in gene_set1:
            for gene2 in gene_set2:
                if gene1 != gene2:
                    try:
                        all_path_lenghts[gene1][gene2] = nx.shortest_path_length(network, source=gene1, target=gene2)
                        all_path_lenghts[gene2][gene1] = all_path_lenghts[gene1][gene2]
                    except:
                        continue
                else:
                    all_path_lenghts[gene1][gene2] = 0
                    all_path_lenghts[gene2][gene1] = 0

        return all_path_lenghts

    def get_shortest_mean_distance(self, targets, network=None):
        """
        Calculate the shortest distances between all targets in the network.

        This function computes the minimum path between one target and any other target
        of the same set, resulting in the intra-node distance or Node module diameter.

        :param network: A NetworkX graph.
        :type network: networkx.Graph

        :param targets: Targets to be analyzed.
        :type targets: list

        :return: Mean of all paths (d_d) and paths (min_paths).
        :rtype: tuple
        """

        if network is None:
            network = self.network

        filtered_targets = self.filter_targets_in_network(network=network, targets=targets)

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

    def calc_set_pair_distances(self, network, given_gene_set1, given_gene_set2):
        """
        Calculates the mean shortest distance between two sets of genes on a given network.

        :param network: The network in networkx format.
        :type network: networkx.Graph

        :param given_gene_set1: The first set of genes for which the distance will be computed.
        :type given_gene_set1: iterable

        :param given_gene_set2: The second set of genes for which the distance will be computed.
        :type given_gene_set2: iterable

        :return: Mean shortest distance between gene_set1 and gene_set2.
        :rtype: float

        """

        # remove all nodes that are not in the network
        all_genes_in_network = set(network.nodes())
        gene_set1 = self.filter_targets_in_network(given_gene_set1)
        gene_set2 = self.filter_targets_in_network(given_gene_set2)

        # get the network distances for all gene pairs:
        all_path_lenghts = self.get_pathlengths_for_two_sets(gene_set1, gene_set2, network)

        all_distances = []

        # going through all pairs starting from set 1
        for geneA in gene_set1:
            all_distances_A = []
            all_distances_B = []
            for geneB in gene_set2:

                # the genes are the same, so their distance is 0
                if geneA == geneB:
                    all_distances_A.append(0)
                    all_distances_B.append(0)
                else:
                    if geneA in all_path_lenghts.keys():
                        if geneB in all_path_lenghts.keys():
                            all_distances_A.append(all_path_lenghts[geneA][geneB])
                            all_distances_B.append(all_path_lenghts[geneB][geneA])
            if len(all_distances_A) > 0:
                l_min = min(all_distances_A)
                all_distances.append(l_min)
            if len(all_distances_B) > 0:
                all_distances.append(min(all_distances_B))

        # calculate mean shortest distance
        mean_shortest_distance = np.mean(all_distances)

        return mean_shortest_distance

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

    def get_shortest_distance_between_target_vs_targets(self, network=None, targets=None):
        """
        Calculate the shortest distances between all targets in the network.

        This function computes the minimum path between one target and any other target
        of the same set, resulting in the intra-hub distance or drug module diameter.

        :param network: the network which should be used.
        :type network: networkx.Graph

        :param targets: Targets to be analyzed.
        :type targets: list

        :return: Mean of all paths (d_d) and paths (min_paths).
        :rtype: tuple
        """

        if network is None:
            network = self.network

        filtered_targets = self.filter_targets_in_network(targets=targets, network=network)

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

    def get_shortest_distances_nodes(self, targets1, targets2, network=None):
        """
        Extract the minimum path between targets.

        This function calculates the minimum path between each target in the first set (targets1)
        and any other target in the second set (targets2) in the provided network.

        :param targets1: List of targets, first set
            :type targets1: list

        :param targets2: List of targets, second set
            :type targets2: list

        :param network: Networkx graph
            The network to be analyzed in networkx format.
            :type network: nx.Graph

        :return: Tuple containing the mean of all paths (d_d) and a list of lists representing the minimum paths
            :rtype: tuple[float, List[List[Any]]]
        """
        if network is None:
            network = self.network

        filtered_targets = self.filter_targets_in_network(targets=targets1, network=network)

        filtered_targets2 = self.filter_targets_in_network(targets=targets2, network=network)

        min_paths = []
        if len(filtered_targets) >= 1 and len(filtered_targets2) >= 1:
            try:
                for t1 in filtered_targets:
                    min_distances = []
                    for t2 in filtered_targets2:
                        if nx.has_path(network, t1, t2):
                            min_distances.append(len(nx.shortest_path(network, t1, t2)) - 1)  # TODO:if no path
                    min_paths.append(min(min_distances))
                return min_paths
            except (nx.exception.NetworkXNoPath, ValueError):
                return None, None
        else:
            return None, None

    # =============================================================================
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

    # =============================================================================
    def get_pathlengths_for_two_sets(self, network, given_gene_set1, given_gene_set2):

        """
        calculate the shortest paths between two given set of genes in a
        given network. The results are stored in a dictionary of
        dictionaries: all_path_lenghts[gene1][gene2] = l with gene1 <
        gene2, so each pair is stored only once!

        PARAMETERS:
        -----------
            - G: network
            - gene_set1/2: gene sets for which paths should be computed

        RETURNS:
        --------
            - all_path_lenghts[gene1][gene2] = l for all pairs of genes
              with gene1 < gene2

        """

        # remove all nodes that are not in the network
        all_genes_in_network = set(network.nodes())
        gene_set1 = given_gene_set1 & all_genes_in_network
        gene_set2 = given_gene_set2 & all_genes_in_network

        all_path_lenghts = {}

        # calculate the distance of all possible pairs
        for gene1 in gene_set1:
            if gene1 not in all_path_lenghts.keys():
                all_path_lenghts[gene1] = {}
            for gene2 in gene_set2:
                if gene1 != gene2:
                    try:
                        l = nx.shortest_path_length(network, source=gene1, target=gene2)
                        if gene1 < gene2:
                            all_path_lenghts[gene1][gene2] = l
                        else:
                            if gene2 not in all_path_lenghts.keys():
                                all_path_lenghts[gene2] = {}
                            all_path_lenghts[gene2][gene1] = l
                    except:
                        continue

        return all_path_lenghts

    # =============================================================================
    def calc_single_set_distance(self, network, given_gene_set):

        """
        Calculates the mean shortest distance for a set of genes on a
        given network


        PARAMETERS:
        -----------
            - G: network
            - gene_set: gene set for which distance will be computed

        RETURNS:
        --------
             - mean shortest distance

        """

        # remove all nodes that are not in the network, just to be safe
        all_genes_in_network = set(network.nodes())
        gene_set = given_gene_set & all_genes_in_network

        # get the network distances for all gene pairs:
        all_path_lenghts = self.get_pathlengths_for_single_set(network, gene_set)

        all_distances = []

        # going through all gene pairs
        for geneA in gene_set:

            all_distances_A = []
            for geneB in gene_set:

                # I have to check which gene is 'smaller' in order to know
                # where to look up the distance of that pair in the
                # all_path_lengths dict
                if geneA < geneB:
                    if geneB in all_path_lenghts[geneA].keys():
                        all_distances_A.append(all_path_lenghts[geneA][geneB])
                else:
                    if geneA in all_path_lenghts[geneB].keys():
                        all_distances_A.append(all_path_lenghts[geneB][geneA])

            if len(all_distances_A) > 0:
                l_min = min(all_distances_A)
                all_distances.append(l_min)

        # calculate mean shortest distance
        mean_shortest_distance = np.mean(all_distances)

        return mean_shortest_distance

    # =============================================================================
    def calc_set_pair_distances(self, network, given_gene_set1, given_gene_set2):

        """
        Calculates the mean shortest distance between two sets of genes on
        a given network

        PARAMETERS:
        -----------
            - G: network
            - gene_set1/2: gene sets for which distance will be computed

        RETURNS:
        --------
             - mean shortest distance

        """

        # remove all nodes that are not in the network
        all_genes_in_network = set(network.nodes())
        gene_set1 = given_gene_set1 & all_genes_in_network
        gene_set2 = given_gene_set2 & all_genes_in_network

        # get the network distances for all gene pairs:
        all_path_lenghts = self.get_pathlengths_for_two_sets(network, gene_set1, gene_set2)

        all_distances = []

        # going through all pairs starting from set 1
        for geneA in gene_set1:

            all_distances_A = []
            for geneB in gene_set2:

                # the genes are the same, so their distance is 0
                if geneA == geneB:
                    all_distances_A.append(0)

                # I have to check which gene is 'smaller' in order to know
                # where to look up the distance of that pair in the
                # all_path_lengths dict
                else:
                    if geneA < geneB:
                        try:
                            all_distances_A.append(all_path_lenghts[geneA][geneB])
                        except:
                            pass

                    else:
                        try:
                            all_distances_A.append(all_path_lenghts[geneB][geneA])
                        except:
                            pass

            if len(all_distances_A) > 0:
                l_min = min(all_distances_A)
                all_distances.append(l_min)

        # going through all pairs starting from disease B
        for geneA in gene_set2:

            all_distances_A = []
            for geneB in gene_set1:

                # the genes are the same, so their distance is 0
                if geneA == geneB:
                    all_distances_A.append(0)

                # I have to check which gene is 'smaller' in order to know
                # where to look up the distance of that pair in the
                # all_path_lengths dict
                else:
                    if geneA < geneB:
                        try:
                            all_distances_A.append(all_path_lenghts[geneA][geneB])
                        except:
                            pass

                    else:
                        try:
                            all_distances_A.append(all_path_lenghts[geneB][geneA])
                        except:
                            pass

            if len(all_distances_A) > 0:
                l_min = min(all_distances_A)
                all_distances.append(l_min)

        # calculate mean shortest distance
        mean_shortest_distance = np.mean(all_distances)

        return mean_shortest_distance

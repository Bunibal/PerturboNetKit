import networkx as nx


class HubDistances:
    """
    This class calculates the distances between hubs or targets in a network
    """

    def __init__(self, network, hub_to_target):
        """
        :param network: networkx graph
        :param hub_to_target: dictionary containing the hubs and their corresponding targets
        """
        self.network = network
        self.hub_to_target = hub_to_target

    def filter_targets_in_network(self, network=None, targets=None):
        """
        This function filters the targets that are not in the network.
        :param network: a networkx graph
        :param targets: dictionary containing the hubs and their corresponding targets
        :return: list of targets that are in the network
        """

        if network is None:
            network = self.network
        if targets is None:
            targets = self.hub_to_target

        filtered_targets = []
        for t in targets:
            if network.has_node(t):
                filtered_targets.append(t)
        return filtered_targets

    def get_shortest_distance_between_target_vs_targets(self, network=None, targets=None):
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

    def get_shortest_distances_hubs(self, targets1, targets2, network=None):
        """
        Extract the min path between targets.
        This is always the minimum path between one target and any other target of the other set.
        Returns Mean of all paths (d_d) as well as paths (min_paths)

        This function uses two sets hence calculates the inter drug distance
        :param targets1: first set of targets
        :param targets2: second set of targets
        :param network: the network to be analyzed in networkx format
        :return:
        """
        if network is None:
            network = self.network

        filtered_targets = self.filter_targets_in_network(targets1, network)

        filtered_targets2 = self.filter_targets_in_network(targets2, network)

        min_paths = []
        if len(filtered_targets) >= 1 and len(filtered_targets2) >= 1:
            try:
                for t1 in filtered_targets:
                    min_distances = []
                    for t2 in filtered_targets2:
                        if nx.has_path(network, t1, t2):
                            min_distances.append(len(nx.shortest_path(network, t1, t2)) - 1)
                    min_paths.append(min(min_distances))
                return min_paths
            except:
                return None, None
        else:
            return None, None

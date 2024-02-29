import matplotlib.pyplot as plt
import networkx as nx
from .Analysisclass import Analysis
import scipy as sc
import numpy as np


# TODO: change all outputs to pdfs?
# TODO: add checks whether there is a network and node target map available
class Networkpipeline(Analysis):
    """
    Class for analyzing networks. Contains functions for plotting network properties. Inherits from Analysis.
    """

    def plot_network(self, network=None, filename=None, only_LCC=False, node_size=50, alpha=0.5):
        """
        Plot the network.

        This method generates a plot of the network using NetworkX.

        :param network: Networkx graph, optional
            The network to be analyzed. If not provided, uses the network passed during class instantiation.
            :type network: nx.Graph or None

        :param filename: str, optional
            The filename to save the plot. If not provided, the plot will be displayed but not saved.
            :type filename: str or None

        :param only_LCC: bool, optional
            If True, plot only the largest connected component of the network. Default is False.
            :type only_LCC: bool

        :param node_size: int, optional
            The size of the nodes in the plot. Default is 50.
            :type node_size: int

        :param alpha: float, optional
            The transparency of the nodes. Should be between 0 (completely transparent) and 1 (completely opaque).
            Default is 0.5.
            :type alpha: float

        :return: None
        """
        if network is None:
            network = self.network
        if filename is None:
            filename = 'network.png'
        if only_LCC:
            Gcc = sorted(nx.connected_components(network), key=len, reverse=True)
            network = network.subgraph(Gcc[0])
        pos = nx.spring_layout(network)
        nx.draw_networkx(network, pos=pos, with_labels=False, node_size=node_size, alpha=alpha)
        plt.savefig(filename)
        plt.close()

    def plot_degree_distribution(self, network=None, filename=None):
        """
        Plot the degree distribution of the network.

        This method generates a plot of the degree distribution using NetworkX.

        :param network: Networkx graph, optional
            The network to be analyzed. If not provided, uses the network passed during class instantiation.
            :type network: nx.Graph or None

        :param filename: str, optional
            The filename to save the plot. If not provided, the plot will be displayed but not saved.
            :type filename: str or None

        :return: None
        """
        if network is None:
            network = self.network
        if filename is None:
            filename = 'degree_distribution.png'
        degrees, degree_distribution, cumulative_degree_distribution = self.get_degree_distribution(network)
        plt.title('Degree distribution')
        plt.plot(degrees, degree_distribution, 'b-', label='Degree distribution')
        plt.plot(degrees, cumulative_degree_distribution, 'r-', label='Cumulative degree distribution')
        plt.xlabel('Degree')
        plt.ylabel('Frequency')
        plt.legend()
        plt.savefig(filename)
        plt.close()

    def plot_degree_distributions(self, network1=None, network2=None, filename=None, network1_name="network1",
                                  network2_name="network2"):
        """
        Plot the degree distributions of two networks in one plot.

        This method generates a plot of the degree distributions for two networks using NetworkX.

        :param network1: Networkx graph
            The first network to be analyzed.
            :type network1: nx.Graph

        :param network2: Networkx graph
            The second network to be analyzed.
            :type network2: nx.Graph

        :param filename: str, optional
            The filename to save the plot. If not provided, the plot will be displayed but not saved.
            :type filename: str or None

        :param network1_name: str, optional
            The name of the first network (used for plot legend). If not provided, defaults to 'Network 1'.
            :type network1_name: str or None

        :param network2_name: str, optional
            The name of the second network (used for plot legend). If not provided, defaults to 'Network 2'.
            :type network2_name: str or None

        :return: None
        """

        if network1 is None:
            network1 = self.network

        if network2 is None:
            if self.node_to_target is None:
                raise ValueError(
                    'No hub to target dictionary provided. Either provide a dictionary or a second network.')
            else:
                network2 = self.get_subnetwork_nodes()

        if filename is None:
            filename = 'degree_distributions.pdf'
            if self.verbose:
                print('No filename provided. Default filename is used.')

        degrees1, degree_distribution1, cumulative_degree_distribution1 = self.get_degree_distribution(network1)
        degrees2, degree_distribution2, cumulative_degree_distribution2 = self.get_degree_distribution(network2)
        # Grid of 2x3 plots
        fig, axs = plt.subplots(2, 3, figsize=(15, 10))

        # Plot 1: Scatter Degree Distribution (Log/Log)
        axs[0, 0].scatter(degrees1, degree_distribution1, c='grey', alpha=0.4)
        axs[0, 0].scatter(degrees2, degree_distribution2, c='#40B9D4', alpha=0.4)
        axs[0, 0].legend(
            [f'{network1_name}', f'{network2_name}\nKS_pValue: {sc.stats.ks_2samp(degrees1, degrees2)[1]:.2e}'],
            frameon=False)
        axs[0, 0].set_xscale('log')
        axs[0, 0].set_yscale('log')
        axs[0, 0].set_xlabel('Degree')
        axs[0, 0].set_ylabel('P(k)')
        axs[0, 0].set_ylim(10 ** -5, 1)

        # Plot 2: Scatter Cumulative Degree Distribution (Log/Log)
        axs[0, 1].scatter(degrees1, cumulative_degree_distribution1, c='grey', alpha=0.4)
        axs[0, 1].scatter(degrees2, cumulative_degree_distribution2, c='#40B9D4', alpha=0.4)
        axs[0, 1].legend(
            [f'{network1_name}', f'{network2_name}\nKS_pValue: {sc.stats.ks_2samp(degrees1, degrees2)[1]:.2e}'],
            frameon=False)
        axs[0, 1].set_xscale('log')
        axs[0, 1].set_yscale('log')
        axs[0, 1].set_xlabel('Degree')
        axs[0, 1].set_ylabel('P(x >= k)')
        axs[0, 1].set_ylim(10 ** -4, 1)

        # Plot 3: Histogram Degree Distribution (Log/Log)
        axs[0, 2].hist(degrees1, bins='auto', alpha=0.4, color='grey')
        axs[0, 2].hist(degrees2, bins='auto', alpha=0.4, color='#40B9D4')
        axs[0, 2].legend(
            [f'{network1_name}', f'{network2_name}\nKS_pValue: {sc.stats.ks_2samp(degrees1, degrees2)[1]:.2e}'],
            frameon=False)
        axs[0, 2].set_xscale('log')
        axs[0, 2].set_yscale('log')
        axs[0, 2].set_xlabel('Degree')
        axs[0, 2].set_ylabel('P(k)')

        # Plot 4: Histogram Degree Distribution (Single Log)
        axs[1, 0].hist(degrees1, bins='auto', alpha=0.4, color='grey')
        axs[1, 0].hist(degrees2, bins='auto', alpha=0.4, color='#40B9D4')
        axs[1, 0].legend(
            [f'{network1_name}', f'{network2_name}\nKS_pValue: {sc.stats.ks_2samp(degrees1, degrees2)[1]:.2e}'],
            frameon=False)
        axs[1, 0].set_xscale('log')
        axs[1, 0].set_xlabel('Degree')
        axs[1, 0].set_ylabel('P(k)')

        # Plot 5: Normal Degree Distribution Line Plot (Log/Log)
        axs[1, 1].plot(degrees1, degree_distribution1, c='grey')
        axs[1, 1].plot(degrees2, degree_distribution2, c='#40B9D4')
        axs[1, 1].legend(
            [f'{network1_name}', f'{network2_name}\nKS_pValue: {sc.stats.ks_2samp(degrees1, degrees2)[1]:.2e}'],
            frameon=False)
        axs[1, 1].set_xscale('log')
        axs[1, 1].set_yscale('log')
        axs[1, 1].set_xlabel('Degree')
        axs[1, 1].set_ylabel('P(k)')

        # Plot 6: Normal Degree Distribution Line Plot (Single Log)
        axs[1, 2].plot(degrees1, degree_distribution1, c='grey')
        axs[1, 2].plot(degrees2, degree_distribution2, c='#40B9D4')
        axs[1, 2].legend(
            [f'{network1_name}', f'{network2_name}\nKS_pValue: {sc.stats.ks_2samp(degrees1, degrees2)[1]:.2e}'],
            frameon=False)
        axs[1, 2].set_xscale('log')
        axs[1, 2].set_xlabel('Degree')
        axs[1, 2].set_ylabel('P(k)')

        # Save the figure
        plt.savefig(filename, format='pdf')
        plt.close()

    def plot_centralities(self, filename='centrality_scatter.pdf'):
        """
        Plot the centrality distributions.

        This method generates a plot of the centrality distributions for a network using NetworkX.

        :param filename: str, optional
            The filename to save the plot. If not provided, the plot will be displayed but not saved.
            :type filename: str or None

        :return: None
        """
        # check if centrality distributions are already calculated
        if self.node_centrality_results is None:
            if self.verbose:
                print('Calculating centrality distributions...')
            self.calculate_centrality_node(network=self.network, node_to_target=self.node_to_target)

        # get foldchange and p_values
        fold_change = []
        p_values = []
        for target in self.node_centrality_results.keys():
            fold_change.append(self.node_centrality_results[target]['FoldChange'])
            p_values.append(self.node_centrality_results[target]['PValue'])

        if self.verbose:
            print('Plotting centrality scatter...')
        # Plot foldchange (=Glass' Delta) and PValues
        plt.scatter(fold_change, p_values, c='#40B9D4', alpha=0.6)
        plt.xlim([min(fold_change), max((fold_change))])
        plt.ylim([min(p_values), max((p_values))])
        plt.yscale('log')
        plt.legend([
            f'Total/Significant:  {round(len(p_values))}/{len([x for x in p_values if x < 0.05])} : {round(len([x for x in p_values if x < 0.05]) / float(len(p_values)))}'],
            frameon=False)
        plt.fill([min(fold_change), max(fold_change), max(fold_change), min(fold_change)], [1, 1, 0.05, 0.05], c='grey',
                 alpha=0.3)
        plt.xlabel("Glass' Delta")
        plt.ylabel('PValue')
        plt.savefig(filename, format='pdf')
        plt.close()

    def plot_lcc_size_results(self, filename='lcc_size_scatter.pdf'):
        """
       Plot the results of the LCC size analysis.

        This method generates a plot of the results obtained from the LCC size analysis.

        :param filename: str, optional
            The filename to save the plot. If not provided, the plot will be displayed but not saved.
            :type filename: str or None

        :return: None
        """
        # check if lcc size results are already calculated
        if self.node_lcc_results is None:
            if self.verbose:
                print('Calculating LCC size results...')
            self.compare_lcc_size_against_random()

        # get the zscores and rel_sizes
        zscores = []
        rel_sizes = []
        for target in self.node_lcc_results:
            zscores.append(self.node_lcc_results[target]['ZScore'])
            rel_sizes.append(self.node_lcc_results[target]['RelativeSize'])

        if self.verbose:
            print(f'Total: {len(zscores)}')

            print(f'Significant: {len([x for x in zscores if abs(x) > 2])}')
            print(len([x for x in zscores if abs(x) > 2]) / float(len(zscores)))

        plt.scatter(rel_sizes, zscores, alpha=0.6, c='#40B9D4')
        plt.legend(
            f'Total/Significant:  {len(zscores)}/{len([x for x in zscores if abs(x) > 2])} ({round(len([x for x in zscores if abs(x) > 2]) / float(len(zscores)))})',
            frameon=False)
        plt.fill([0, 0, 1, 1], [-2, 2, 2, -2], color='grey', alpha=0.4)
        plt.xlabel('relative module size s=Nd /S')
        plt.ylabel('z-score of module size S')
        plt.xlim([min(rel_sizes), max((rel_sizes))])
        plt.ylim([min(zscores), max((zscores))])
        plt.savefig(filename, format='pdf')

    def plot_shortest_path_between_targets(self, filename="shortest_path_between_targets.pdf"):
        """
        Plot the results of check_shortest_path_between_targets_against_random.

        This method generates a plot of the results obtained from the function check_shortest_path_between_targets_against_random.

        :param filename: str, optional
            The filename to save the plot. If not provided, the plot will be displayed but not saved.
            :type filename: str or None

        :return: None
        """

        # check if shortest path results are already calculated
        if self.node_shortest_mean_paths is None:
            if self.verbose:
                print('Calculating shortest path results...')
            self.check_shortest_path_between_targets_against_random()

        # get p_values and glass_deltas
        p_value = []
        glass_delta = []
        for target in self.node_shortest_mean_paths:
            p_value.append(self.node_shortest_mean_paths[target]['PValue'])
            glass_delta.append(self.node_shortest_mean_paths[target]['FoldChange'])

        if self.verbose:
            print("Plotting shortest path results...")
            print(f'Total: {len(p_value)}')
            print(f'Significant: {len([x for x in p_value if x < 0.05])}')
            print(f'Significant: {len([x for x in p_value if x < 0.05]) / float(len(p_value)):.2f}')

        plt.scatter(glass_delta, p_value, alpha=0.6, c='#40B9D4')
        plt.legend([
            f'Total/Significant: {len(p_value)}/{len([x for x in p_value if abs(x) < 0.05])} ({len([x for x in p_value if abs(x) < 0.05]) / len(p_value):.2f})'],
            frameon=False)
        plt.ylim(min(p_value), 1)
        plt.xlim([min(glass_delta), max(glass_delta)])
        plt.yscale('log')
        plt.xlabel("Glass' Delta")
        plt.ylabel('PValue')
        plt.savefig(filename, format='pdf')

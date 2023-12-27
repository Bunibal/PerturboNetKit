import matplotlib.pyplot as plt
import networkx as nx
from Analysisclass import Analysis
import scipy as stats
import numpy as np


class PipelineAnalysisNetwork(Analysis):
    def plot_degree_distribution(self, network=None, filename=None):
        """
        This function plots the degree distribution of the network.
        :param network: The network to be analyzed in networkx format
        :param filename: The filename of the plot
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
        This function plots the degree distributions of two networks in one plot.

        :param network1: The first network to be analyzed in networkx format
        :param network2: The second network to be analyzed in networkx format
        :param filename: The filename of the plot
        :param network1_name: The name of the first network
        :param network2_name: The name of the second network
        :return: None
        """

        if network1 is None:
            network1 = self.network

        if network2 is None:
            if self.node_to_target is None:
                raise ValueError(
                    'No hub to target dictionary provided. Either provide a dictionary or a second network.')
            else:
                network2 = self.get_subnetwork_hub()

        if filename is None:
            filename = 'degree_distribution.pdf'
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
            [f'{network1_name}', f'{network2_name}\nKS_pValue: {stats.ks_2samp(degrees1, degrees2)[1]:.2e}'],
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
            [f'{network1_name}', f'{network2_name}\nKS_pValue: {stats.ks_2samp(degrees1, degrees2)[1]:.2e}'],
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
            [f'{network1_name}', f'{network2_name}\nKS_pValue: {stats.ks_2samp(degrees1, degrees2)[1]:.2e}'],
            frameon=False)
        axs[0, 2].set_xscale('log')
        axs[0, 2].set_yscale('log')
        axs[0, 2].set_xlabel('Degree')
        axs[0, 2].set_ylabel('P(k)')

        # Plot 4: Histogram Degree Distribution (Single Log)
        axs[1, 0].hist(degrees1, bins='auto', alpha=0.4, color='grey')
        axs[1, 0].hist(degrees2, bins='auto', alpha=0.4, color='#40B9D4')
        axs[1, 0].legend(
            [f'{network1_name}', f'{network2_name}\nKS_pValue: {stats.ks_2samp(degrees1, degrees2)[1]:.2e}'],
            frameon=False)
        axs[1, 0].set_xscale('log')
        axs[1, 0].set_xlabel('Degree')
        axs[1, 0].set_ylabel('P(k)')

        # Plot 5: Normal Degree Distribution Line Plot (Log/Log)
        axs[1, 1].plot(degrees1, degree_distribution1, c='grey')
        axs[1, 1].plot(degrees2, degree_distribution2, c='#40B9D4')
        axs[1, 1].legend(
            [f'{network1_name}', f'{network2_name}\nKS_pValue: {stats.ks_2samp(degrees1, degrees2)[1]:.2e}'],
            frameon=False)
        axs[1, 1].set_xscale('log')
        axs[1, 1].set_yscale('log')
        axs[1, 1].set_xlabel('Degree')
        axs[1, 1].set_ylabel('P(k)')

        # Plot 6: Normal Degree Distribution Line Plot (Single Log)
        axs[1, 2].plot(degrees1, degree_distribution1, c='grey')
        axs[1, 2].plot(degrees2, degree_distribution2, c='#40B9D4')
        axs[1, 2].legend(
            [f'{network1_name}', f'{network2_name}\nKS_pValue: {stats.ks_2samp(degrees1, degrees2)[1]:.2e}'],
            frameon=False)
        axs[1, 2].set_xscale('log')
        axs[1, 2].set_xlabel('Degree')
        axs[1, 2].set_ylabel('P(k)')

        # Save the figure
        plt.savefig(filename, format='pdf')
        plt.close()

    def plot_centralities(self, filename='centrality_scatter.pdf'):
        """
        Plot the centrality distributions
        :param filename: Filename to save the plot
        :return: None
        """
        # check if centrality distributions are already calculated
        if self.hub_centrality_results is None:
            if self.verbose:
                print('Calculating centrality distributions...')
            self.calculate_centrality_hub(network=self.network, targets=self.node_to_target)

        # get foldchange and p_values
        fold_change = []
        p_values = []
        for target in self.hub_centrality_results.keys():
            fold_change.append(self.hub_centrality_results[target]['fold_change'])
            p_values.append(self.hub_centrality_results[target]['p_value'])

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
        Plot the results of the LCC size analysis
        :return: None
        """
        # check if lcc size results are already calculated
        if self.hub_lcc_results is None:
            if self.verbose:
                print('Calculating LCC size results...')
            self.compare_lcc_size_against_random()

        # get the zscores and rel_sizes
        zscores = []
        rel_sizes = []
        for target in self.hub_lcc_results:
            zscores.append(self.hub_lcc_results[target]['ZScore'])
            rel_sizes.append(self.hub_lcc_results[target]['RelativeSize'])

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
        Plot the results of check_shortest_path_between_targets_against_random
        :param filename:
        :return: None
        """

        # check if shortest path results are already calculated
        if self.cloud_ShortestPaths_results is None:
            if self.verbose:
                print('Calculating shortest path results...')
            self.check_shortest_path_between_targets_against_random()

        # get p_values and glass_deltas
        p_value = []
        glass_delta = []
        for target in self.cloud_ShortestPaths_results:
            p_value.append(self.cloud_ShortestPaths_results[target]['PValue'])
            glass_delta.append(self.cloud_ShortestPaths_results[target]['FoldChange'])

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

    def go_term_similarity_barplot(self, filename=""):
        """
        Plot the results of calculate_similarity_maxspecificity.
        :param filename: path and filename to save the plot
        :return: None
        """
        for go_branch in self.similarties_bin_means:
            plt.bar(range(1, max(self.similarties_bin_means[1][2]) + 1),
                    [np.mean(self.similarties_bin_means[x]) for x in
                     range(1, max(self.similarties_bin_means[1][2]) + 1)],
                    yerr=[
                        1.96 * (np.std(self.similarties_bin_means[x][0]) / float(len(self.similarties_bin_means[x][0])))
                        for x in
                        range(1, max(self.similarties_bin_means[1][2]) + 1)], align='center', alpha=0.5, ecolor='black',
                    capsize=10,
                    color='#40B9D4', zorder=2)
            plt.xticks(range(1, max(self.similarties_bin_means[1][2]) + 1),
                       ['{:2.2f}'.format(x) for x in self.similarties_bin_means[1][1]])
            plt.xlabel("Glass' Delta")
            plt.ylabel("GO Similarity (%s)\n%s Similarity (95%% CI)" % (go_branch, "MaxSpecifity"))
            # plt.show()
            plt.savefig(
                filename + go_branch + '_MaxSpecifity_BarPlot.pdf',
                format='pdf')
            plt.close()

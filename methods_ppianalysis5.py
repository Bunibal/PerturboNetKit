import networkx as nx
from matplotlib import pylab as plt
import numpy as np
from scipy import stats
import random as rand
from scipy.stats import mannwhitneyu
import gene2terms_addupstream as GO
from scipy.stats import binned_statistic
import random
import sys as sys
import os


def calculate_Centrality(PPI, targets):
    '''
    Extract centralities for nodes on the PPI. Return list of centralities or None if no node is found on the PPI
    '''
    filtered_targets = []
    for t in targets:
        if PPI.has_node(t):
            filtered_targets.append(t)

    if len(filtered_targets) > 0:
        centralities = []
        for node in filtered_targets:
            centralities.append(nx.closeness_centrality(PPI, node))

        return centralities


    else:
        return None


def CalculateCLOUD_Centralities(interactome, cloud_all, res_dir):
    '''
    Check the centralities of normal PPI nodes and nodes targeted by the CLOUD perturbers
    '''

    # Load PPI
    PPI = nx.read_gml(interactome)

    # Get all the CLOUD Targets
    cloud_targets = {}
    different_Targets = set()
    with open(cloud_all, 'r') as fp:
        next(fp)
        for line in fp:
            tmp = line.strip().split(',')
            cloud_targets[tmp[0]] = list(set(tmp[2].split(';')))
            for t in tmp[2].split(';'):
                different_Targets.add(t)
    print('Number of different targets: %d' % len(different_Targets))

    # create a random distribution of centralities on the PPI
    NumRandom = 10000
    print('Create Random Distribution')

    random_distances = []
    for i in range(0, NumRandom):
        node = rand.sample(PPI.nodes(), 1)[0]
        random_distances.append(nx.closeness_centrality(PPI, node))

    # calculate the centralities of the targets for an individual drug
    cloud_centrality_results = {}
    for c in cloud_targets.keys():
        print(c)
        if len(cloud_targets[c]) == 0:
            continue

        # get the centralities for all targets for a specific CLOUD drug
        drug_centralities = calculate_Centrality(PPI, cloud_targets[c])

        # claculate the PValue,MeanCentrality and Glass' Delta if targets found on the PPI
        if drug_centralities != None:
            pValue = mannwhitneyu(drug_centralities, random_distances)[1]
            # FoldChange = d_d / np.mean(random_distances)
            GlassDelta = (np.mean(drug_centralities) - np.mean(random_distances)) / np.std(random_distances)

            cloud_centrality_results[c] = {'MeanCentrality': np.mean(drug_centralities), 'PValue': pValue,
                                           'FoldChange': GlassDelta}

    # Save the results
    with open(res_dir + '/Centralities.csv', 'w') as fp_out:
        fp_out = open(res_dir + '/Centralities.csv', 'w')
        fp_out.write('CLOUD,MeanCentrality,PValue,FoldChange\n')
        clouds = cloud_centrality_results.keys()
        clouds = sorted(clouds)

        for c in clouds:
            fp_out.write(c + ',' + str(cloud_centrality_results[c]['MeanCentrality']) + ',' + str(
                cloud_centrality_results[c]['PValue']) + ',' + str(cloud_centrality_results[c]['FoldChange']) + '\n')


def Plot_Centralities(res_dir, graph_dir):
    '''
    Function to plot the resulting Glass' Deltas of centralities
    '''

    # Get the calculated centrality results
    pValues = []
    FoldChange = []
    with open(res_dir + '/Centralities.csv') as fp:
        next(fp)
        for line in fp:
            tmp = line.strip().split(',')

            pValues.append(float(tmp[2]))
            FoldChange.append(float(tmp[3]))

    # Plot foldchange (=Glass' Delta) and PValues
    plt.scatter(FoldChange, pValues, c='#40B9D4', alpha=0.6)
    plt.xlim([min(FoldChange), max((FoldChange))])
    plt.ylim([min(pValues), max((pValues))])
    plt.yscale('log')
    plt.legend(['Total/Significant:  %d/%d (%.2f)' % (
        len(pValues), len([x for x in pValues if x < 0.05]),
        len([x for x in pValues if x < 0.05]) / float(len(pValues)))],
               frameon=False)
    plt.fill([min(FoldChange), max(FoldChange), max(FoldChange), min(FoldChange)], [1, 1, 0.05, 0.05], c='grey',
             alpha=0.3)
    plt.xlabel("Glass' Delta")
    plt.ylabel('PValue')
    # plt.show()
    plt.savefig(graph_dir + '/Centrality_Scatter.pdf', format='pdf')


def Check_LCC_Size(PPI, targets):
    '''
    Check the LCC size of given targets on the PPI

    '''
    LCC = 1

    SubGraph = nx.subgraph(PPI, targets)
    if len(SubGraph.nodes) > 0:
        components = (SubGraph.subgraph(c) for c in nx.connected_components(SubGraph))
        LCC = max([len(x.nodes()) for x in components])

    return LCC


def create_LCC_Results(interactome, cloud_all, res_dir):
    # Load PPI
    PPI = nx.read_gml(interactome)

    # Get all the CLOUD Targets
    cloud_targets = {}
    different_Targets = set()
    with open(cloud_all, 'r') as fp:
        next(fp)
        for line in fp:
            tmp = line.strip().split(',')
            cloud_targets[tmp[0]] = list(set(tmp[2].split(';')))
            for t in tmp[2].split(';'):
                different_Targets.add(t)

    # Calculate LCC sizes and compare to randomly drawn LCC on the PPI of same size as the number of targets of the specific drug
    cloud_LCCs_results = {}
    for c in cloud_targets.keys():
        print(c)

        # Get specific drug LCC
        if len(cloud_targets[c]) == 0:
            continue
        lcc = Check_LCC_Size(PPI, cloud_targets[c])
        number_of_Genes = len(cloud_targets[c])
        relative_LCC_Size = float(lcc) / number_of_Genes

        # Create random distribution of LCCs
        random_distribution = []
        for i in range(0, 10000):
            randGenes = rand.sample(PPI.nodes(), len(cloud_targets[c]))
            rand_LCC = Check_LCC_Size(PPI, randGenes)
            random_distribution.append(rand_LCC)

        # Calculate ZScore
        ZScore = (lcc - np.mean(random_distribution)) / np.std(random_distribution)

        cloud_LCCs_results[c] = {'lcc': lcc, 'ZScore': ZScore, '#Genes': number_of_Genes,
                                 'RelativeSize': relative_LCC_Size}

    # Save Output
    with open(res_dir + '/LCC_Sizes.csv', 'w') as fp_out:
        fp_out.write('CLOUD,LCC,ZScore,#Genes,RelativeSize\n')
        clouds = cloud_LCCs_results.keys()
        clouds = sorted(clouds)

        for c in clouds:
            fp_out.write(
                c + ',' + str(cloud_LCCs_results[c]['lcc']) + ',' + str(cloud_LCCs_results[c]['ZScore']) + ',' + str(
                    cloud_LCCs_results[c]['#Genes']) + ',' + str(cloud_LCCs_results[c]['RelativeSize']) + '\n')


def Plot_LCC_Results(res_dir, graph_dir):
    '''
    Use the LCC output file and create a plot

    :return:
    '''

    with open(res_dir + '/LCC_Sizes.csv') as fp:
        zscores = []
        rel_sizes = []
        next(fp)
        for line in fp:
            tmp = line.strip().split(',')
            if tmp[2] != 'nan' and int(tmp[3]) > 1:
                zscores.append(float(tmp[2]))
                rel_sizes.append(float(tmp[4]))

    print('Total: %d' % len(zscores))

    print('Significant: %d' % len([x for x in zscores if abs(x) > 2]))
    print(len([x for x in zscores if abs(x) > 2]) / float(len(zscores)))

    # Plot relative size (e.g. percent of targets within LCC and the ZScore)
    plt.scatter(rel_sizes, zscores, alpha=0.6, c='#40B9D4')
    plt.legend(['Total/Significant:  %d/%d (%.2f)' % (len(zscores), len([x for x in zscores if abs(x) > 2]),
                                                      len([x for x in zscores if abs(x) > 2]) / float(len(zscores)))],
               frameon=False)
    plt.fill([0, 0, 1, 1], [-2, 2, 2, -2], color='grey', alpha=0.4)
    plt.xlabel('relative module size s=Nd /S')
    plt.ylabel('z-score of module size S')
    plt.xlim([min(rel_sizes), max((rel_sizes))])
    plt.ylim([min(zscores), max((zscores))])
    plt.savefig(graph_dir + '/LCC.pdf', format='pdf')


def Check_Shortest_Distances(PPI, targets):
    '''
    Extract the min path between targets.
    This is always the minimum path between one target and any other target of the other set.
    Returns Mean of all paths (d_d) as well as paths (min_paths)

    This function uses only one set hence calulcates the intra drug distance or drug_module diamter
    '''
    filtered_targets = []
    for t in targets:
        if PPI.has_node(t):
            filtered_targets.append(t)

    min_paths = []
    if len(filtered_targets) > 1:
        try:
            for t1 in filtered_targets:
                min_distances = []
                for t2 in filtered_targets:
                    if t1 != t2:

                        if nx.has_path(PPI, t1, t2):
                            min_distances.append(len(nx.shortest_path(PPI, t1, t2)) - 1)
                min_paths.append(min(min_distances))
            d_d = sum(min_paths) / float(len(filtered_targets))

            return d_d, min_paths
        except:
            return None, None
    else:
        return None, None


def Create_Shortest_Distances_Output(interactome, cloud_all, res_dir):
    '''
    Function to calculate the intra_drug module distance e.g. shortest distances of any target of a given drug to any other target of this drug (mean over all)

    '''

    # Load PPI
    PPI = nx.read_gml(interactome)

    # Get all the CLOUD Targets
    cloud_targets = {}
    different_Targets = set()
    with open(cloud_all) as fp:
        next(fp)
        for line in fp:
            tmp = line.strip().split(',')
            cloud_targets[tmp[0]] = list(set(tmp[2].split(';')))
            for t in tmp[2].split(';'):
                different_Targets.add(t)

    # Create random distribution for randomly drawn proteins on the PPI
    NumRandom = 10000
    print('Create Random Distribution')
    random_distances = []
    for i in range(0, NumRandom):
        targets = rand.sample(PPI.nodes(), 2)
        if nx.has_path(PPI, targets[0], targets[1]):
            random_distances.append(len(nx.shortest_path(PPI, targets[0], targets[1])) - 1)

    # Calculate the drug module PPI diameter
    cloud_ShortestPaths_results = {}
    for c in cloud_targets.keys():
        print(c)
        if len(cloud_targets[c]) == 0:
            continue

        # Extract min distances
        d_d, min_paths = Check_Shortest_Distances(PPI, cloud_targets[c])

        if d_d == None:
            continue

        # Calculate pValue by comparing with random Distribution
        pValue = mannwhitneyu(min_paths, random_distances)[1]

        # Calculate fold change (Glass' Delata)
        GlassDelta = (d_d - np.mean(random_distances)) / np.std(random_distances)

        # Save Result
        cloud_ShortestPaths_results[c] = {'MeanShortestPath': d_d, 'PValue': pValue, 'FoldChange': GlassDelta}

    # Save result to output
    with open(res_dir + '/ShortestPaths.csv', 'w') as fp_out:
        fp_out.write('CLOUD,MeanShortestPath,PValue,FoldChange\n')
        clouds = cloud_ShortestPaths_results.keys()
        clouds = sorted(clouds)

        for c in clouds:
            fp_out.write(c + ',' + str(cloud_ShortestPaths_results[c]['MeanShortestPath']) + ',' + str(
                cloud_ShortestPaths_results[c]['PValue']) + ',' + str(
                cloud_ShortestPaths_results[c]['FoldChange']) + '\n')


def Plot_SPath_Results(results_dir, graph_dir):
    '''
    Use the drug module diameter output file and create a plot

    :return:
    '''

    # Load the drug module diamter results
    with open(results_dir + '/ShortestPaths.csv') as fp:
        next(fp)
        pValue = []
        GlassDelta = []
        for line in fp:
            tmp = line.strip().split(',')
            pValue.append(float(tmp[2]))
            GlassDelta.append(float(tmp[3]))

    print('Total: %d' % len(pValue))
    print('Significant: %d' % len([x for x in pValue if x < 0.05]))
    print('Significant: %.2f' % (len([x for x in pValue if x < 0.05]) / float(len(pValue))))

    # Plot output again with foldchange (glass delta) and Significance (pValue)
    plt.scatter(GlassDelta, pValue, alpha=0.6, c='#40B9D4')
    plt.legend(['Total/Significant:  %d/%d  (%.2f)' % (len(pValue), len([x for x in pValue if abs(x) < 0.05]),
                                                       len([x for x in pValue if abs(x) < 0.05]) / float(
                                                           len(pValue)))], frameon=False)
    plt.ylim(min(pValue), 1)
    plt.xlim([min(GlassDelta), max(GlassDelta)])
    plt.yscale('log')
    plt.xlabel("Glass' Delta")
    plt.ylabel('PValue')
    plt.savefig(graph_dir + '/ShortestDistance.pdf', format='pdf')


def Check_Shortest_Distances(PPI, targets):
    '''
    Extract the min path between targets.
    This is always the minimum path between one target and any other target of the other set.
    Returns Mean of all paths (d_d) as well as paths (min_paths)

    This function uses only one set hence calulcates the intra drug distance or drug_module diamter

    '''
    filtered_targets = []
    for t in targets:
        if PPI.has_node(t):
            filtered_targets.append(t)

    min_paths = []
    if len(filtered_targets) > 1:
        try:
            for t1 in filtered_targets:
                min_distances = []
                for t2 in filtered_targets:
                    if t1 != t2:
                        if nx.has_path(PPI, t1, t2):
                            min_distances.append(len(nx.shortest_path(PPI, t1, t2)) - 1)
                min_paths.append(min(min_distances))
            d_d = sum(min_paths) / float(len(filtered_targets))

            return d_d, min_paths
        except:
            return None, None
    else:
        return None, None


def Check_Shortest_DistancesBetween(PPI, targets1, targets2):
    '''
    Extract the min path between targets.
    This is always the minimum path between one target and any other target of the other set.
    Returns Mean of all paths (d_d) as well as paths (min_paths)

    This function uses two sets hence calulcates the inter drug distance

    '''
    filtered_targets = []
    for t in targets1:
        if PPI.has_node(t):
            filtered_targets.append(t)

    filtered_targets2 = []
    for t in targets2:
        if PPI.has_node(t):
            filtered_targets2.append(t)

    min_paths = []
    if len(filtered_targets) >= 1 and len(filtered_targets2) >= 1:
        try:
            for t1 in filtered_targets:
                min_distances = []
                for t2 in filtered_targets2:
                    if nx.has_path(PPI, t1, t2):
                        min_distances.append(len(nx.shortest_path(PPI, t1, t2)) - 1)
                min_paths.append(min(min_distances))
            return min_paths
        except:
            return None, None
    else:
        return None, None



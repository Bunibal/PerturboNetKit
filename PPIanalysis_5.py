import os.path

from methods_ppianalysis5 import *

args = sys.argv
ppi = args[1]
cloud_targetssplit = args[2]
cloud_all = args[3]
cloud_off = args[4]
if len(sys.argv) > 4:
    graph_dir = args[5]
    results_dir = args[6]
else:
    graph_dir = 'graphs'
    results_dir = 'results'

if len(sys.argv) < 5:
    raise Exception('Please provide the path to the PPI network as an argument.')
if not os.path.exists(ppi):
    raise FileNotFoundError('The provided PPI network does not exist.')
if not os.path.exists(cloud_targetssplit):
    raise FileNotFoundError('The provided cloud target split file does not exist.')
if not os.path.exists(cloud_all):
    raise FileNotFoundError('The provided cloud all file does not exist.')
if not os.path.exists(cloud_off):
    raise FileNotFoundError('The provided cloud off file does not exist.')
if not os.path.isdir(graph_dir):
    os.makedirs(graph_dir)
    print('Created directory %s' % graph_dir)
if not os.path.isdir(results_dir):
    os.makedirs(results_dir)
    print('Created directory %s' % results_dir)

PPI = nx.read_gml(args[1])

# get cloud targets
cloud_targetssplit = {}
different_Targets = set()
with open(args[2], 'r') as fp:
    for line in fp:
        tmp = line.strip().split(',')
        cloud_targetssplit[tmp[0]] = list(set(tmp[1].split(';')))
        for t in tmp[1].split(';'):
            different_Targets.add(t)

print('Number of different targets: %d' % len(different_Targets))

# Extract the degree distribution of the whole PPI
degrees_PPI = [x[1] for x in nx.degree(PPI)]
degrees_PPI_unique = list(set(degrees_PPI))
degrees_PPI_unique.sort()
degreesPPI = []
degreeDistributionPPI = []
cumulativedegreeDistributionPPI = []
for degree in degrees_PPI_unique:
    degreesPPI.append(degree)
    degreeDistributionPPI.append(degrees_PPI.count(degree) / float(len(degrees_PPI)))
    cumulativedegreeDistributionPPI.append(len([x for x in degrees_PPI if x >= degree]) / float(len(degrees_PPI)))

# Extract the degree distribution of the subpart of the PPI containing targeted proteins
degrees_Drugs = [x[1] for x in PPI.degree(different_Targets)]
degrees_Drugs_unique = list(set(degrees_Drugs))
degrees_Drugs_unique.sort()
degreesDrugs = []
degreeDistributionDrugs = []
cumulativedegreeDistributionDrugs = []
for degree in degrees_Drugs_unique:
    degreesDrugs.append(degree)
    degreeDistributionDrugs.append(degrees_Drugs.count(degree) / float(len(degrees_Drugs)))
    cumulativedegreeDistributionDrugs.append(len([x for x in degrees_Drugs if x >= degree]) / float(len(degrees_Drugs)))

print('Mean PPI degree: %.2f' % np.mean(degrees_PPI))
print('Mean Drug degree: %.2f' % np.mean(degrees_Drugs))

# Plot the normal degree distribution (log/log)
plt.scatter(degreesPPI, degreeDistributionPPI, c='grey', alpha=0.4)
plt.scatter(degreesDrugs, degreeDistributionDrugs, c='#40B9D4', alpha=0.4)
plt.legend(['PPI', 'CLOUD\nKS_pValue: %.2e' % stats.ks_2samp(degrees_PPI, degrees_Drugs)[1]], frameon=False)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Degree')
plt.ylabel('P(k)')
plt.ylim(10 ** -5, 1)
plt.savefig(graph_dir + '/Scatter_DegreeDistribution_LogLog.pdf', format='pdf')
plt.close()

# Plot the cumulative degree distribution (log/log)
plt.scatter(degreesPPI, cumulativedegreeDistributionPPI, c='grey', alpha=0.4)
plt.scatter(degreesDrugs, cumulativedegreeDistributionDrugs, c='#40B9D4', alpha=0.4)
plt.legend(['PPI', 'CLOUD\nKS_pValue: %.2e' % stats.ks_2samp(degrees_PPI, degrees_Drugs)[1]], frameon=False)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Degree')
plt.ylabel('P(x >= k)')
plt.ylim(10 ** -4, 1)
plt.savefig(graph_dir + '/Scatter_CumulativeDegreeDistribution_LogLog.pdf', format='pdf')
plt.close()

# Plot a histogram degree distribution (log/log)
PPI_Bins = plt.hist(degrees_PPI, bins='auto', alpha=0.4, color='grey')
Drug_Bins = plt.hist(degrees_Drugs, bins='auto', alpha=0.4, color='#40B9D4')
plt.legend(['PPI', 'CLOUD\nKS_pValue: %.2e' % stats.ks_2samp(degrees_PPI, degrees_Drugs)[1]], frameon=False)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Degree')
plt.ylabel('P(k)')
plt.savefig(graph_dir + '/Histogram_LogLog.pdf', format='pdf')
plt.close()

plt.hist(degrees_PPI, bins='auto', alpha=0.4, color='grey')
plt.hist(degrees_Drugs, bins='auto', alpha=0.4, color='#40B9D4')
plt.legend(['PPI', 'CLOUD\nKS_pValue: %.2e' % stats.ks_2samp(degrees_PPI, degrees_Drugs)[1]], frameon=False)
plt.xscale('log')
plt.xlabel('Degree')
plt.ylabel('P(k)')
plt.savefig(graph_dir + '/Histogram_Log.pdf', format='pdf')
plt.close()

# Plot a noraml degree distribution line plot (log/log)
plt.plot(degreesPPI, degreeDistributionPPI, c='grey')
plt.plot(degreesDrugs, degreeDistributionDrugs, c='#40B9D4')
plt.legend(['PPI', 'CLOUD\nKS_pValue: %.2e' % stats.ks_2samp(degrees_PPI, degrees_Drugs)[1]], frameon=False)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Degree')
plt.ylabel('P(k)')
plt.savefig(graph_dir + '/LinePlot_LogLog.pdf', format='pdf')
plt.close()

# Plot a normal degree distribution line plot (single log)
plt.plot(degreesPPI, degreeDistributionPPI, c='grey')
plt.plot(degreesDrugs, degreeDistributionDrugs, c='#40B9D4')
plt.legend(['PPI', 'CLOUD\nKS_pValue: %.2e' % stats.ks_2samp(degrees_PPI, degrees_Drugs)[1]], frameon=False)
plt.xscale('log')
plt.xlabel('Degree')
plt.ylabel('P(k)')
plt.savefig(graph_dir + '/LinePlot_Log.pdf', format='pdf')
plt.close()

CalculateCLOUD_Centralities(interactome=ppi, cloud_all=cloud_all, res_dir=results_dir)
Plot_Centralities(res_dir=results_dir, graph_dir=graph_dir)

create_LCC_Results(interactome=ppi, cloud_all=cloud_all, res_dir=results_dir)
Plot_LCC_Results(res_dir=results_dir, graph_dir=graph_dir)

# 4. Shortest Path

Create_Shortest_Distances_Output(interactome=ppi, cloud_all=cloud_all, res_dir=results_dir)
Plot_SPath_Results(results_dir=results_dir, graph_dir=graph_dir)

# 5. Calculate Overlap Between two Drugs

cloud_targets = {}
different_Targets = set()

# Calculate Sab based on only Targets (remove transporters and enzymes)
with open('../data/PPI_Analysis/CLOUD_to_TargetsSplit.csv', 'r') as fp:
    next(fp)
    for line in fp:
        tmp = line.strip().split(',')
        cloud_targets[tmp[0]] = list(set(tmp[1].split(';')))
        for t in tmp[1].split(';'):
            different_Targets.add(t)

# Calcualte within distances
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

# Calculate separation between two drugs Sab
with open(results_dir + '/Separation_TargetsOnly.csv', 'w') as fp_out:
    fp_out.write('Drug1,Drug2,D_Drug1,D_Drug2,D_D1_D2,S\n')
    for c in clouds:
        print(c)
        d_A = within_Distances[c]
        targets1 = cloud_targets[c]

        for c2 in clouds:
            # if c > c2:
            d_B = within_Distances[c2]
            targets2 = cloud_targets[c2]
            distances1 = Check_Shortest_DistancesBetween(PPI, targets1, targets2)
            distances2 = Check_Shortest_DistancesBetween(PPI, targets2, targets1)

            # Dab
            between_Distance = (sum(distances1) + sum(distances2)) / float((len(distances1) + len(distances2)))

            # Sab
            separation = between_Distance - (d_A + d_B) / 2.0

            fp_out.write(
                c + ',' + c2 + ',' + str(d_A) + ',' + str(d_B) + ',' + str(between_Distance) + ',' + str(
                    separation) + '\n')


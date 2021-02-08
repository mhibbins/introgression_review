import sys
from ete3 import Tree
import numpy as np
import scipy.stats

trees = []

with open(sys.argv[1], 'r') as treefile:
    for line in treefile:
        if line.startswith("("):
            trees.append(Tree(line))

def calc_D3(treeset):

    D_53_vals, D_43_vals = 0, 0

    for tree in treeset: 
        dist_53 = tree.get_distance("6", "3", topology_only = True)
        dist_45 = tree.get_distance("6", "5", topology_only = True)
        dist_43 = tree.get_distance("5", "3", topology_only = True)

        if dist_53 <= dist_45 and dist_53 <= dist_43:
            D_53_vals += 1
        elif dist_43 <= dist_45 and dist_43 <= dist_53:
            D_43_vals += 1

    try:
        return (D_53_vals-D_43_vals)/float(D_53_vals+D_43_vals)
    except: 
        return "No discordant trees"


def D3_bootstrap(treeset, n_replicates):

    D3_replicates = []

    for i in range(n_replicates):

        bootstrap_treeset = []

        for j in range(len(treeset)):
            index = np.random.randint(0, len(treeset)-1)
            bootstrap_treeset.append(treeset[index])
        
        if calc_D3(bootstrap_treeset) == "No discordant trees":
            return "No discordance"
        else:    
            D3_replicates.append(calc_D3(bootstrap_treeset))

    mean_D3 = sum(D3_replicates)/float(len(D3_replicates))

    D3_stdev = ((sum([(x - mean_D3)**2 for x in D3_replicates]))/float((len(D3_replicates)-1)))**(0.5)
    D3_se = D3_stdev/float(len(D3_replicates)**(0.5))
    CI_lower = mean_D3 - (1.96)*float(D3_stdev)
    CI_upper = mean_D3 + (1.96)*float(D3_stdev)

    #if mean_D3 > 0:
    #    pval = sum(1 for i in D3_replicates if i <= 0)/float(len(D3_replicates))
    #elif mean_D3 < 0:
    #    pval = sum(1 for i in D3_replicates if i >= 0)/float(len(D3_replicates))

    z_score = mean_D3/float(D3_stdev)
    pval = scipy.stats.norm.sf(abs(z_score))*2

    return mean_D3, D3_stdev, D3_se, CI_lower, CI_upper, pval

if D3_bootstrap(trees, 100) == "No discordance":
    print("No discordance, delta = 0")
else:
    mean_D3, D3_stdev, D3_se, CI_lower, CI_upper, pval = D3_bootstrap(trees, 1000)
    print(mean_D3, D3_stdev, D3_se, CI_lower, CI_upper, pval)

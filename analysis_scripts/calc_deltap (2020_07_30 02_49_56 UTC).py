import sys
from ete3 import Tree

def parse_results(resultsfile):

    pvals = []

    with open(resultsfile, "r") as results:
        for line in results:
            pvals.append(line.split()[5])

    return pvals


pvals = parse_results(sys.argv[1])

def parse_seqfile(i, seqfile):

    if float(pvals[i-1]) < 0.05:

        seqs = []

        with open(seqfile, "r") as seq_file:
            for line in seq_file:
                if "5 1000" not in line:
                    seqs.append(str.strip(line))

        genome_5, genome_4, genome_3, genome_1 = [], [], [], []

        for i in range(len(seqs)):
            if "1" in seqs[i]:
                genome_1.append(seqs[i].split()[1])
            elif "3" in seqs[i]:
                genome_3.append(seqs[i].split()[1])
            elif "4" in seqs[i]:
                genome_4.append(seqs[i].split()[1])
            elif "5" in seqs[i]:
                genome_5.append(seqs[i].split()[1])

        return ''.join(genome_5), ''.join(genome_4), ''.join(genome_3), ''.join(genome_1)

    else:

        return "Test statistic not significant"

def parse_treefile(treefile):
    
    trees = []
    
    with open(treefile, "r") as tree_file:
        for line in tree_file:
            if line.startswith("("):
                trees.append(Tree(line))

    return trees

if parse_seqfile(int(sys.argv[2]), sys.argv[3]) == "Test statistic not significant":
    print("Test statistic not significant")
    sys.exit()
else:
    P1, P2, P3, outgroup = parse_seqfile(int(sys.argv[2]), sys.argv[3])
    trees = parse_treefile(sys.argv[4])

#def calc_DP(P1, P2, P3, out):
#
#    ABBA_count, BABA_count, BBAA_count = 0, 0, 0
#
#    for i in range(len(P1)):
#        if P2[i] == P3[i] and P1[i] != P3[i] and P1[i] == out[i]:
#            ABBA_count += 1
#        elif P1[i] == P3[i] and P2[i] != P3[i] and P2[i] == out[i]:
#            BABA_count += 1
#        elif P1[i] == P2[i] and P1[i] != P3[i] and P3[i] == out[i]:
#            BBAA_count += 1

#    DP = (float(ABBA_count) - float(BABA_count)) / (float(ABBA_count) + float(BABA_count) + float(BBAA_count))

#    return abs(DP) 

def calc_deltap(treeset):

    D_54_vals, D_53_vals, D_43_vals = 0, 0, 0

    for tree in treeset:
        dist_54 = tree.get_distance("5", "4", topology_only = True)
        dist_53 = tree.get_distance("5", "3", topology_only = True)
        dist_43 = tree.get_distance("4", "3", topology_only = True)

        if dist_54 <= dist_53 and dist_54 <= dist_43:
            D_54_vals += 1
        elif dist_53 <= dist_54 and dist_53 <= dist_43:
            D_53_vals += 1
        elif dist_43 <= dist_54 and dist_43 <= dist_53:
            D_43_vals += 1

    deltap = (D_53_vals-D_43_vals)/float(D_54_vals+D_53_vals+D_43_vals)

    return abs(deltap)

deltap = calc_deltap(trees)

print(deltap)



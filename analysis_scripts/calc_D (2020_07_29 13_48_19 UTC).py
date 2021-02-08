import sys
import random
import scipy.stats

def parse_seqfile(seqfile):

    seqs = []

    with open(seqfile, 'r') as seq_file:
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


def calc_D(P1, P2, P3, out):

    ABBA_count, BABA_count = 0, 0

    for i in range(len(P1)):
        if P2[i] == P3[i] and P1[i] != P3[i] and P1[i] == out[i]:
            ABBA_count += 1
        elif P1[i] == P3[i] and P2[i] != P3[i] and P2[i] == out[i]:
            BABA_count += 1

    D = (float(ABBA_count) - float(BABA_count)) / (float(ABBA_count) + float(BABA_count))

    return ABBA_count, BABA_count, D

def D_bootstrap(P1, P2, P3, outgroup, n_replicates):

    D_estimates = []

    for i in range(n_replicates):

        P1_bootstrapped, P2_bootstrapped, P3_bootstrapped, outgroup_bootstrapped = [], [], [], []
        
        for j in range(3000):

            sampling_index = random.randint(0, len(P1)-1000)

            P1_bootstrapped.append(P1[sampling_index:sampling_index+1000])
            P2_bootstrapped.append(P2[sampling_index:sampling_index+1000])
            P3_bootstrapped.append(P3[sampling_index:sampling_index+1000])
            outgroup_bootstrapped.append(outgroup[sampling_index:sampling_index+1000])

        P1_bootstrapped = ''.join(P1_bootstrapped)
        P2_bootstrapped = ''.join(P2_bootstrapped)
        P3_bootstrapped = ''.join(P3_bootstrapped)
        outgroup_bootstrapped = ''.join(outgroup_bootstrapped)

        ABBA, BABA, D = calc_D(P1_bootstrapped, P2_bootstrapped, P3_bootstrapped, outgroup_bootstrapped)
        D_estimates.append(D)

    mean_D = sum(D_estimates)/float(len(D_estimates))

    D_stdev = ((sum([(x - mean_D)**2 for x in D_estimates]))/float((len(D_estimates)-1)))**(0.5)
    D_se = D_stdev/float(len(D_estimates)**(0.5))
    CI_lower = mean_D - (1.96)*float(D_stdev)
    CI_upper = mean_D + (1.96)*float(D_stdev)

    #if mean_D > 0:
    #    pval = sum(1 for i in D_estimates if i <= 0)/float(len(D_estimates))
    #elif mean_D < 0:
    #    pval = sum(1 for i in D_estimates if i >= 0)/float(len(D_estimates))
    
    z_score = mean_D/float(D_stdev)
    pval = scipy.stats.norm.sf(abs(z_score))*2

    return mean_D, D_stdev, D_se, CI_lower, CI_upper, pval

P1, P2, P3, outgroup = parse_seqfile(sys.argv[1])
mean_D, D_stdev, D_se, CI_lower, CI_upper, pval = D_bootstrap(P1, P2, P3, outgroup, 1000)
print(mean_D, D_stdev, D_se, CI_lower, CI_upper, pval)



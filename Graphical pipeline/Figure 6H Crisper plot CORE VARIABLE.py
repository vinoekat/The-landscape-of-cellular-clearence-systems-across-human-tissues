import csv
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statistics

''''-----------------------------------This code generates plot for figure 5E----------------------------------------'''
'''To obtain the file needed to run this code please run "Figure 5E CRISPR_score_min_CORE_VARIABLE.py"'''

def glob_chap_PC_ks(lists):
    '''Calculates KS test for 2 groups of genes and creates a cumulative graph'''
    cor, arr2, prot = lists
    cor.sort()
    arr2.sort()
    prot.sort()
    joined = cor + arr2 + prot
    joined.sort()


    ### Pick the colors you want. change the next section in plt.hist
    # core_color = "#4a8499"
    # variable_color = "#69b5d0"
    # protein_coding = '#495f6d'

    #AUT
    core_color = "#79815f"
    variable_color = "#c0cd97"
    protein_coding = '#495f6d'



    #This controls the bins of the graph (for smoother graph, increase the number 100)
    bins = np.append(np.linspace(min(joined), max(joined), 1000000), [np.inf])
    plt.hist(cor, bins = bins, density=1, histtype='step', cumulative=True,color = core_color,linewidth=3)
    plt.hist(arr2, bins = bins, density=1, histtype='step', cumulative=True,color = variable_color,linewidth=3)
    plt.hist(prot, bins=bins, density=1, histtype='step', cumulative=True, color=protein_coding, linewidth=3)

    # statistics results
    print('cor/var', scipy.stats.ks_2samp(cor, arr2, alternative= 'greater'))
    print('cor/prot', scipy.stats.ks_2samp(cor, prot, alternative='greater'))
    print('prot/var', scipy.stats.ks_2samp(arr2, prot, alternative='greater'))
    print(scipy.stats.mannwhitneyu(cor, arr2), 'mw')
    print(scipy.stats.mannwhitneyu(cor, prot), 'mw')
    print(scipy.stats.mannwhitneyu(arr2, prot), 'mw')
    # number of chaperones, number of protein coding genes in the analysis
    print('Number of core = ',len(cor))
    print('Number of variable = ',len(arr2))
    print(len(prot))

    #print(all(isinstance(x, float) for x in arr1))

    #If the X axis is to wide, you can change its limits here
    plt.xlim(min(joined), max(joined))
    plt.ylim(0, 1)

    ### save figure --> choose path
    #plt.savefig('DepMap/pierre/Crisper plot.png')
    #plt.legend()
    plt.show()

gene_set='AUT'
#Choose the file to take data from
depmap = pd.read_csv('DepMap_CRISPR_min_score_core_var_dif_V8_' + gene_set + ' with prot coding consortium.csv',sep=',')


cor = []
arr2 = []
prot=[]

#Choose conditions for groups to compare
for rowNum, row in depmap.iterrows():
        if row['core'] == "Yes":
            cor.append(float(row['Min']))
        if row['core'] == 'prot coding':
            prot.append(float(row['Min']))
        if row['core'] == 'No':
            arr2.append(float(row['Min']))

#The actual function
glob_chap_PC_ks(lists=(cor,arr2, prot))

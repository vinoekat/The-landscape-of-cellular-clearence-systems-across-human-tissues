import csv
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statistics

''''-----------------------------------This code generates plot for figure 7F----------------------------------------'''
'''To obtain the file needed to run this code please run "Figure 7F CRISPR_score_min.py"'''

def glob_chap_PC_ks(lists):
    '''Calculates KS test for 2 groups of genes and creates a cumulative graph'''
    arr1, arr2, arr3= lists
    arr1.sort()
    arr2.sort()
    arr3.sort()

    joined = arr1 + arr2 + arr3
    joined.sort()


    ### Pick the colors you want. change the next section in plt.hist
    deg_color = "#a7b07b"
    fold_color = '#fce597'
    syn_color = "#8b0000"



    #This controls the bins of the graph (for smoother graph, increase the number 100)
    bins = np.append(np.linspace(min(joined), max(joined), 10000000), [np.inf])
    plt.hist(arr1, bins=bins, density=1, histtype='step', cumulative=True, color=deg_color, linewidth=3)
    plt.hist(arr2, bins=bins, density=1, histtype='step', cumulative=True, color=fold_color, linewidth=3)
    plt.hist(arr3, bins=bins, density=1, histtype='step', cumulative=True, color=syn_color, linewidth=3)


    # statistics results
    print(scipy.stats.ks_2samp(arr1, arr2, alternative='less'), 'deg vs fold less')
    print(scipy.stats.ks_2samp(arr1, arr3, alternative='less'), 'deg vs syn less')
    print(scipy.stats.ks_2samp(arr3, arr2, alternative='greater'), 'syn vs fold less')

    # number of genes in the analysis
    print('Number of deg = ', len(arr1))
    print('Number of fold = ', len(arr2))
    print('Number of syn = ', len(arr3))


    #If the X axis is to wide, you can change its limits here
    plt.xlim(min(joined), max(joined))
    plt.ylim(0, 1)

    ### save figure --> choose path
    #plt.savefig('DepMap/pierre/Crisper plot.png')
    plt.legend()
    plt.show()


#Choose the file to take data from
depmap = pd.read_csv("DepMap_CRISPR_min_score_deg_vs_fold.csv",sep=',')


arr1 = []
arr2 = []
arr3 = []


#Choose conditions for groups to compare
for rowNum, row in depmap.iterrows():
        if row['Degradation'] == "Yes":
            arr1.append(float(row['Min']))
        if row['Degradation'] == "No":
                arr2.append(float(row['Min']))
        if row['Degradation'] == "SYN":
                arr3.append(float(row['Min']))

#The actual function
glob_chap_PC_ks(lists=(arr1, arr2, arr3))

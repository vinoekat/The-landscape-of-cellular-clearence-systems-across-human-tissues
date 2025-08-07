import csv
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statistics

''''-----------------------------------This code generates plot for figure 1F----------------------------------------'''
'''To obtain the file needed to run this code please run "Figure 1F CRISPR_score_min.py"'''

def glob_chap_PC_ks(lists):
    '''Calculates KS test for 2 groups of genes and creates a cumulative graph'''
    arr1, arr2, arr3,= lists
    arr1.sort()
    arr2.sort()
    arr3.sort()
    joined = arr1 + arr2 + arr3
    joined.sort()


    ### Pick the colors you want. change the next section in plt.hist
    aut_color = "#1395ba"
    ups_color = '#a2b86c'
    pro_color = "#0d3c55"


    #This controls the bins of the graph (for smoother graph, increase the number 100)
    bins = np.append(np.linspace(min(joined), max(joined), 10000000), [np.inf])
    plt.hist(arr1, bins=bins, density=1, histtype='step', cumulative=True, color=aut_color, linewidth=3)
    plt.hist(arr3, bins=bins, density=1, histtype='step', cumulative=True, color=ups_color, linewidth=3)
    plt.hist(arr2, bins=bins, density=1, histtype='step', cumulative=True, color=pro_color, linewidth=3)

    # statistics results
    print(scipy.stats.ks_2samp(arr1, arr2, alternative='greater'), 'aut vs prot')
    print(scipy.stats.ks_2samp(arr3, arr2,alternative='greater'), 'ups vs prot')
    print(scipy.stats.ks_2samp(arr3, arr1,alternative='greater'), 'ups vs aut')
    # number of genes in the analysis
    print('Number of aut = ', len(arr1))
    print('Number of ups = ', len(arr3))
    print('Number of protein-coding genes = ', len(arr2))

    #If the X axis is to wide, you can change its limits here
    plt.xlim(min(joined), max(joined))
    plt.ylim(0, 1)

    ### save figure --> choose path
    #plt.savefig('DepMap/pierre/Crisper plot.png')
    plt.legend()
    plt.show()


#Choose the file to take data from
depmap = pd.read_csv("DepMap_CRISPR_min_score_consortium.csv",sep=',')


arr1 = []
arr2 = []
arr3 = []

#Choose conditions for groups to compare
for rowNum, row in depmap.iterrows():
        if row['Autophage'] == "Yes":
            arr1.append(float(row['Min']))
        if row['Autophage'] == "UPS":
            arr3.append(float(row['Min']))
        if row['Autophage'] == "No":
                arr2.append(float(row['Min']))

#The actual function
glob_chap_PC_ks(lists=(arr1, arr2, arr3))

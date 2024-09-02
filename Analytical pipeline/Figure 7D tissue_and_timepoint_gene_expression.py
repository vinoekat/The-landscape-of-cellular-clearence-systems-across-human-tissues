import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.stats.multitest as stm

'''------------------------------------Source files required to run this code----------------------------------------'''
gene_set = 'Deg vs folding'
aut_file = pd.read_csv('all GTex names conslist ENSG UniprotID AUT no dups w our genes.csv')
ups_file = pd.read_csv('all GTex names conslist ENSG UniprotID UPS no dups w our genes.csv')
chap_file = pd.read_csv('ENSG new gene list filtered to relevant chaperones.csv')
syn_file = pd.read_csv('consortium file with ENSG based on UniprotID SYN no dup.csv')
dev_data = pd.read_csv('Devo1_Supp_hum1.csv', sep=',')#to obtain this file please refer to the supplementary table of the
#"Gene expression across mammalian organ development" article, PMID: 31243369,Supplementary table 6.



'''------------------------------Calculate generality score for each gene in 3 systems-------------------------------'''
deg_dict_ENSG = {}
fold_dict_ENSG = {}
syn_dict_ENSG = {}

tissue_list = ['TissueTau', 'BrainTau', 'CerebellumTau', 'HeartTau', 'KidneyTau', 'LiverTau', 'OvaryTau', 'TestisTau']

for lineNum, line in aut_file.iterrows():
    line = line.to_dict()
    deg_dict_ENSG[line['ENSG']] = line['Name']
print(deg_dict_ENSG)

for lineNum, line in ups_file.iterrows():
    line = line.to_dict()
    deg_dict_ENSG[line['ENSG']] = line['Name']


for lineNum, line in chap_file.iterrows():
    line = line.to_dict()
    fold_dict_ENSG[line['ENSG']] = line['Name']

for lineNum, line in syn_file.iterrows():
    line = line.to_dict()
    syn_dict_ENSG[line['ENSG']] = line['Name']

fd = ['Tissue', 'ENSG', 'Name', 'System', 'Score']

writer = csv.DictWriter(
    open('dev_tissue_and_time_specificity_' + gene_set + '.csv', 'w', newline=''),
    fieldnames=fd)
writer.writeheader()
write_row = {}

for lineNum, line in dev_data.iterrows():
    line = line.to_dict()
    for item in line:
        if item in tissue_list:
            if line[item] == line[item]:
                if line['Human_ID'] in deg_dict_ENSG:
                    write_row['System'] = 'Degradation'
                    write_row['Name'] = deg_dict_ENSG[line['Human_ID']]
                    write_row['ENSG'] = line['Human_ID']
                    write_row['Tissue'] = item
                    write_row['Score'] = 1 - float(line[item])
                    writer.writerow(write_row)

                if line['Human_ID'] in fold_dict_ENSG:
                    write_row['System'] = 'Folding'
                    write_row['Name'] = fold_dict_ENSG[line['Human_ID']]
                    write_row['ENSG'] = line['Human_ID']
                    write_row['Tissue'] = item
                    write_row['Score'] = 1 - float(line[item])

                    writer.writerow(write_row)

                if line['Human_ID'] in syn_dict_ENSG:
                    write_row['System'] = 'Synthesis'
                    write_row['Name'] = syn_dict_ENSG[line['Human_ID']]
                    write_row['ENSG'] = line['Human_ID']
                    write_row['Tissue'] = item
                    write_row['Score'] = 1 - float(line[item])

                    writer.writerow(write_row)


'''------------------------------------------Statistical analysis----------------------------------------------------'''


def MannWhitney():

    tissue_and_time = pd.read_csv('dev_tissue_and_time_specificity_' + gene_set + '.csv')
    tissue_list2 = ['TissueTau', 'BrainTau', 'CerebellumTau', 'HeartTau', 'KidneyTau', 'LiverTau', 'OvaryTau',
                   'TestisTau']

    fd = ['Tissue', 'MW Deg vs Fold', 'BH Deg vs Fold', 'MW Syn vs Fold', 'BH Syn vs Fold', 'MW Deg vs Syn', 'BH Deg vs Syn']

    writer = csv.DictWriter(
        open('dev_tissue_and_time_specificity_MannWhitney_' + gene_set + '.csv', 'w', newline=''),
        fieldnames=fd)
    writer.writeheader()
    write_row = {}

    p_values_dict = {}

    bh_dict = {}
    k1_list = []
    k2_list = []
    k3_list = []

    for item in tissue_list2:
        deg_temp = []
        fold_temp = []
        syn_temp = []
        for rowNum, row in tissue_and_time.iterrows():
            row = row.to_dict()
            if item == row['Tissue'] and row['System'] == 'Degradation':
               deg_temp.append(float(row['Score']))
            if item == row['Tissue'] and row['System'] == 'Folding':
                fold_temp.append(float(row['Score']))

            if item == row['Tissue'] and row['System'] == 'Synthesis':
                syn_temp.append(float(row['Score']))

        #print(len(deg_temp), deg_temp)

        k1 = stats.mannwhitneyu(deg_temp, fold_temp)
        k2 = stats.mannwhitneyu(deg_temp, syn_temp)
        k3 = stats.mannwhitneyu(syn_temp, fold_temp)

        p_values_dict[item] = [k1[1], k2[1], k3[1]]
        bh_dict[item] = [0,0,0]
        k1_list.append(k1[1])
        k2_list.append(k2[1])
        k3_list.append(k3[1])

    BH_k1 = stm.multipletests(k1_list, alpha=0.05, method='fdr_bh', is_sorted=False,
                               returnsorted=False)
    BH_k2 = stm.multipletests(k2_list, alpha=0.05, method='fdr_bh', is_sorted=False,
                              returnsorted=False)
    BH_k3 = stm.multipletests(k3_list, alpha=0.05, method='fdr_bh', is_sorted=False,
                              returnsorted=False)


    i=0
    for tissue in p_values_dict:
        write_row['Tissue'] = tissue
        write_row['MW Deg vs Fold'] = p_values_dict[tissue][0]
        write_row['MW Deg vs Syn'] = p_values_dict[tissue][1]
        write_row['MW Syn vs Fold'] = p_values_dict[tissue][2]
        write_row['BH Deg vs Fold'] = BH_k1[1][i]
        write_row['BH Deg vs Syn'] = BH_k2[1][i]
        write_row['BH Syn vs Fold'] = BH_k3[1][i]

        writer.writerow(write_row)
        i += 1



MannWhitney()


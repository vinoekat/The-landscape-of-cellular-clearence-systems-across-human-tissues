import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.stats.multitest as stm


'''-------------------------------------Files needed to run this code------------------------------------------------'''
gene_set = 'UPS' #change to AUT for autophagy (ALP in article)
core_var_file = pd.read_csv('core_variable_differential_V8_' + gene_set + '_5TPM_core_consortium_trimmed.csv')
dev_data = pd.read_csv('Devo1_Supp_hum1.csv', sep=',')
other_protein_coding = pd.read_csv('protein coding 2023.csv')
'''------------------------------------------------------------------------------------------------------------------'''



def MannWhitney():
    '''function uses the file with generality score and calculates statistical significance using Mann-WHitney U test and
    Benjamini-Hochberg correction '''
    tissue_and_time = pd.read_csv('developmental_tissue_and_time_specificity_core_variable_V8_'+ gene_set + ' consortium.csv')
    tissue_list2 = ['TissueTau', 'BrainTau', 'CerebellumTau', 'HeartTau', 'KidneyTau', 'LiverTau', 'OvaryTau',
                   'TestisTau']

    fd = ['Tissue', 'MW Core vs Var', 'MW Core vs Prot', 'MW var vs Prot', 'BH Core vs Var', 'BH Core vs Prot',\
          'BH var vs Prot']

    writer = csv.DictWriter(
        open('developmental_tissue_and_time_specificity_core_variable_V8_MannWhitney_file_' + gene_set + ' consortium.csv', 'w', newline=''),
        fieldnames=fd)
    writer.writeheader()
    write_row = {}

    p_values_dict = {}

    bh_dict = {}
    k1_list = []
    k2_list = []
    k3_list = []

    for item in tissue_list2:
        core_temp = []
        variable_temp = []
        protein_temp = []
        for rowNum, row in tissue_and_time.iterrows():
            row = row.to_dict()
            if item == row['Tissue'] and row['Group'] == 'Core':
               core_temp.append(float(row['Score']))
            elif item == row['Tissue'] and row['Group'] == 'dVariable':
                variable_temp.append(float(row['Score']))
            elif item == row['Tissue'] and row['Group'] == 'Protein-coding':
                protein_temp.append(float(row['Score']))
        print(len(core_temp), core_temp)

        k1 = stats.mannwhitneyu(core_temp,variable_temp)
        k2 = stats.mannwhitneyu(core_temp,protein_temp)
        k3 = stats.mannwhitneyu(variable_temp,protein_temp)

        p_values_dict[item] = [k1[1], k2[1], k3[1]]
        bh_dict[item] = [0, 0, 0]
        k1_list.append(k1[1]) #core var
        k2_list.append(k2[1]) #core prot
        k3_list.append(k3[1]) #var prot

    BH_k1 = stm.multipletests(k1_list, alpha=0.05, method='fdr_bh', is_sorted=False,
                              returnsorted=False)
    BH_k2 = stm.multipletests(k2_list, alpha=0.05, method='fdr_bh', is_sorted=False,
                              returnsorted=False)
    BH_k3 = stm.multipletests(k3_list, alpha=0.05, method='fdr_bh', is_sorted=False,
                              returnsorted=False)

    i = 0
    for tissue in p_values_dict:
        write_row['Tissue'] = tissue
        write_row['MW Core vs Var'] = p_values_dict[tissue][0]
        write_row['MW Core vs Prot'] = p_values_dict[tissue][1]
        write_row['MW var vs Prot'] = p_values_dict[tissue][2]

        write_row['BH Core vs Var'] = BH_k1[1][i]
        write_row['BH Core vs Prot'] = BH_k2[1][i]
        write_row['BH var vs Prot'] = BH_k3[1][i]

        writer.writerow(write_row)
        i += 1


'''this section calculates generality score based on the external data and writes this to file'''


'''get gene ENSGs for each group'''
core_list_ENSG = []
variable_list_ENSG = []
protein_coding_ENSG = []
tissue_list = ['TissueTau', 'BrainTau', 'CerebellumTau', 'HeartTau', 'KidneyTau', 'LiverTau', 'OvaryTau', 'TestisTau']

for lineNum, line in core_var_file.iterrows():
    line = line.to_dict()
    if line['core/var'] == 'core' and line['ENSG'] not in core_list_ENSG:
        core_list_ENSG.append(line['ENSG'])
    elif line['core/var'] == 'variable' and line['ENSG'] not in variable_list_ENSG:
        variable_list_ENSG.append(line['ENSG'])


for lineNum, line in other_protein_coding.iterrows():
    line = line.to_dict()
    if line['ENSG'] not in protein_coding_ENSG:
        protein_coding_ENSG.append(line['ENSG'])

'''pull the data from the data set, calculate generality score and write it to file '''
fd = ['Tissue','ENSG', 'Group', 'Score']

writer = csv.DictWriter(
    open('developmental_tissue_and_time_specificity_core_variable_V8_'+ gene_set+ ' consortium.csv', 'w', newline=''),
    fieldnames=fd)
writer.writeheader()
write_row = {}

for lineNum, line in dev_data.iterrows():
    line = line.to_dict()
    for item in line:
        if item in tissue_list:
            if line[item] == line[item]:
                write_row['ENSG'] = line['Human_ID']
                write_row['Tissue'] = item
                write_row['Score'] = 1- float(line[item])
                if line['Human_ID'] in core_list_ENSG:
                    write_row['Group'] = 'Core'
                elif line['Human_ID'] in variable_list_ENSG:
                    write_row['Group'] = 'dVariable'
                elif line['Human_ID'] in protein_coding_ENSG:
                    write_row['Group'] = 'Protein-coding'

                writer.writerow(write_row)

MannWhitney()
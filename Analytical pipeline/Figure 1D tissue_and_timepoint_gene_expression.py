import numpy as np
import csv
import pandas as pd
from scipy import stats
import statsmodels.stats.multitest as stm

'''-----------------------------------Source files required to run this code-----------------------------------------'''
gene_set = 'ALL'
aut_file = pd.read_csv('all GTex names conslist ENSG UniprotID AUT no dups w our genes.csv')
ups_file = pd.read_csv('all GTex names conslist ENSG UniprotID UPS no dups w our genes.csv')
dev_data = pd.read_csv('Devo1_Supp_hum1.csv', sep=',') #to obtain this file please refer to the supplementary table of the
#"Gene expression across mammalian organ development" article, PMID: 31243369,Supplementary table 6.
other_protein_coding = pd.read_csv('protein coding 2023.csv')



'''------------------------------Calculate generality score for each gene in 3 systems-------------------------------'''
aut_dict_ENSG = {}
ups_dict_ENSG = {}
protein_coding_ENSG = {}
tissue_list = ['TissueTau', 'BrainTau', 'CerebellumTau', 'HeartTau', 'KidneyTau', 'LiverTau', 'OvaryTau', 'TestisTau']

for lineNum, line in aut_file.iterrows():
    line = line.to_dict()
    aut_dict_ENSG[line['ENSG']] = line['Name']
#print(aut_dict_ENSG)

for lineNum, line in ups_file.iterrows():
    line = line.to_dict()
    ups_dict_ENSG[line['ENSG']] = line['Name']


for lineNum, line in other_protein_coding.iterrows():
    line = line.to_dict()
    protein_coding_ENSG[line['ENSG']] = line['Name']

print(len(protein_coding_ENSG))


fd = ['Tissue', 'ENSG', 'Name', 'System', 'Score']

writer = csv.DictWriter(
    open('developmental_tissue_and_time_specificity_core_variable_V8_' + gene_set + '.csv', 'w', newline=''),
    fieldnames=fd)
writer.writeheader()
write_row = {}

for lineNum, line in dev_data.iterrows():
    line = line.to_dict()
    for item in line:
        if item in tissue_list:
            if line[item] == line[item]:
                # write_row['ENSG'] = line['Human_ID']
                # write_row['Tissue'] = item
                # write_row['Score'] = 1 - float(line[item])
                if line['Human_ID'] in aut_dict_ENSG:
                    write_row['System'] = 'AUT'
                    write_row['Name'] = aut_dict_ENSG[line['Human_ID']]
                    write_row['ENSG'] = line['Human_ID']
                    write_row['Tissue'] = item
                    write_row['Score'] = 1 - float(line[item])
                    writer.writerow(write_row)

                if line['Human_ID'] in ups_dict_ENSG:
                    write_row['System'] = 'abUPS'
                    write_row['Name'] = ups_dict_ENSG[line['Human_ID']]
                    write_row['ENSG'] = line['Human_ID']
                    write_row['Tissue'] = item
                    write_row['Score'] = 1 - float(line[item])
                    writer.writerow(write_row)

                if line['Human_ID'] in protein_coding_ENSG and line['Human_ID'] not in aut_dict_ENSG\
                        and line['Human_ID'] not in ups_dict_ENSG:
                    write_row['System'] = 'Protein-coding'
                    write_row['Name'] = protein_coding_ENSG[line['Human_ID']]
                    write_row['ENSG'] = line['Human_ID']
                    write_row['Tissue'] = item
                    write_row['Score'] = 1 - float(line[item])
                    writer.writerow(write_row)



'''------------------------------------------Statistical analysis----------------------------------------------------'''
def MannWhitney():

    tissue_and_time = pd.read_csv('developmental_tissue_and_time_specificity_core_variable_V8_' + gene_set + '.csv')
    tissue_list2 = ['TissueTau', 'BrainTau', 'CerebellumTau', 'HeartTau', 'KidneyTau', 'LiverTau', 'OvaryTau',
                   'TestisTau']

    fd = ['Tissue', 'MW AUT vs UPS', 'MW AUT vs Prot', 'MW UPS vs Prot', \
          'BH AUT vs UPS', 'BH AUT vs Prot', 'BH UPS vs Prot']

    writer = csv.DictWriter(
        open('developmental_tissue_and_time_specificity_core_variable_V8_MannWhitney_file_' + gene_set + '.csv', 'w', newline=''),
        fieldnames=fd)
    writer.writeheader()
    write_row = {}

    p_values_dict = {}

    bh_dict = {}
    k1_list = []
    k2_list = []
    k3_list = []
    k4_list = []
    k5_list = []
    k6_list = []
    for item in tissue_list2:
        aut_temp = []
        ups_temp = []
        protein_temp = []
        for rowNum, row in tissue_and_time.iterrows():
            row = row.to_dict()
            if item == row['Tissue'] and row['System'] == 'AUT':
               aut_temp.append(float(row['Score']))
            if item == row['Tissue'] and row['System'] == 'abUPS':
                ups_temp.append(float(row['Score']))
            if item == row['Tissue'] and row['System'] == 'Protein-coding':
                protein_temp.append(float(row['Score']))
        print(len(aut_temp), aut_temp)

        k1 = stats.mannwhitneyu(aut_temp, ups_temp)
        k4 = stats.mannwhitneyu(aut_temp, protein_temp)
        k5 = stats.mannwhitneyu(ups_temp, protein_temp)


        '''dict: tissue:['MW AUT vs UPS', 'MW AUT vs Prot', 'MW UPS vs Prot']'''
        p_values_dict[item] = [k1[1], k4[1], k5[1]]
        bh_dict[item] = [0, 0, 0]
        k1_list.append(k1[1])

        k4_list.append(k4[1])
        k5_list.append(k5[1])



    BH_k1 = stm.multipletests(k1_list, alpha=0.05, method='fdr_bh', is_sorted=False,
                               returnsorted=False)

    BH_k4 = stm.multipletests(k4_list, alpha=0.05, method='fdr_bh', is_sorted=False,
                              returnsorted=False)
    BH_k5 = stm.multipletests(k5_list, alpha=0.05, method='fdr_bh', is_sorted=False,
                              returnsorted=False)


    i=0
    for tissue in p_values_dict:
        write_row['Tissue'] = tissue
        write_row['MW AUT vs UPS'] = p_values_dict[tissue][0]

        write_row['MW AUT vs Prot'] = p_values_dict[tissue][1]
        write_row['MW UPS vs Prot'] = p_values_dict[tissue][2]

        write_row['BH AUT vs UPS'] = BH_k1[1][i]

        write_row['BH AUT vs Prot'] = BH_k4[1][i]
        write_row['BH UPS vs Prot'] = BH_k5[1][i]

        writer.writerow(write_row)
        i+=1



MannWhitney()
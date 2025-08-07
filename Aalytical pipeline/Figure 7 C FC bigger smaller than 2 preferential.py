import csv
import pickle
import scipy.stats as stats
import statistics
import pandas as pd

'''--------------------------------Source files required to run this code--------------------------------------------'''
aut_file = pd.read_csv('core_variable_differential_V8_AUT_5TPM_core_consortium_trimmed.csv')
ups_file = pd.read_csv('core_variable_differential_V8_UPS_5TPM_core_consortium_trimmed.csv')
chap_file = pd.read_csv('core_variable_differential_V8_Chap_5TPM_core_consortium_trimmed.csv')
syn_file = pd.read_csv('core_variable_differential_V8_SYN_5TPM_core_consortium_trimmed.csv')

dif_score_file = pd.read_csv('TCGA prefernetial expression.csv')

'''------------------------------------------------------------------------------------------------------------------'''

'''Script calculates percent of genes expressed above FC = 1, i.e. upregulated. To obtain down regulated change condition to -2 '''
condition = 2

aut_dict = {}
ups_dict = {}
chap_dict = {}
syn_dict = {}
combined_list = []
protein_dict = {}

for rowNUm, row in aut_file.iterrows():
    row = row.to_dict()
    aut_dict[row['ENSG']] = []
    aut_dict[row['ENSG']].append(row['Name'])
    aut_dict[row['ENSG']].append(row['Class'])
    aut_dict[row['ENSG']].append(row['core/var'])
    combined_list.append(row['ENSG'])

for rowNUm, row in ups_file.iterrows():
    row = row.to_dict()
    ups_dict[row['ENSG']] = []
    ups_dict[row['ENSG']].append(row['Name'])
    ups_dict[row['ENSG']].append(row['Class'])
    ups_dict[row['ENSG']].append(row['core/var'])
    combined_list.append(row['ENSG'])

for rowNUm, row in chap_file.iterrows():
    row = row.to_dict()
    if row['ENSG'] == row['ENSG']:
        chap_dict[row['ENSG']] = []
        chap_dict[row['ENSG']].append(row['Name'])
        # ups_dict[row['ENSG']].append(row['category'])
        # ups_dict[row['ENSG']].append(row['group'])
        combined_list.append(row['ENSG'])


for rowNUm, row in syn_file.iterrows():
    row = row.to_dict()
    syn_dict[row['ENSG']] = []
    syn_dict[row['ENSG']].append(row['Name'])
    syn_dict[row['ENSG']].append(row['Class'])
    syn_dict[row['ENSG']].append(row['core/var'])
    combined_list.append(row['ENSG'])


def percent_FC_bigger_than_1():
    tissues = dif_score_file.columns.tolist()
    tissues.pop(0)
    tissues.pop(0)

    fc_dict = {}
    for tissue in tissues:
        fc_dict[tissue] = {'3AUT': [0, 0], '4UPS': [0, 0], '2Chap': [0, 0], '1SYN': [0, 0]}

    for lineNum, line in dif_score_file.iterrows():
        line = line.to_dict()
        for item in line:
            if item not in ['ENSG', 'Name']:
                if line['ENSG'] in aut_dict:
                    fc_dict[item]['3AUT'][1] += 1
                    if float(line[item]) <= condition:
                        fc_dict[item]['3AUT'][0] += 1

                if line['ENSG'] in ups_dict:
                    fc_dict[item]['4UPS'][1] += 1
                    if float(line[item]) <= condition:
                        fc_dict[item]['4UPS'][0] += 1

                if line['ENSG'] in chap_dict:
                    fc_dict[item]['2Chap'][1] += 1
                    if float(line[item]) <= condition:
                        fc_dict[item]['2Chap'][0] += 1

                if line['ENSG'] in syn_dict:
                    fc_dict[item]['1SYN'][1] += 1
                    if float(line[item]) <= condition:
                        fc_dict[item]['1SYN'][0] += 1


    fd = ['Tissue', 'Percent FC above 1', 'System']
    writer = csv.DictWriter(open('FC above 2 per tissue in AUT UPS Chap and SYN consortium.csv', 'w', newline=''), fieldnames=fd)
    writer.writeheader()
    write_row = {}

    for tissue in fc_dict:

        for system in fc_dict[tissue]:
            write_row['Tissue'] = tissue
            write_row['System'] = system
            write_row['Percent FC above 2'] = (fc_dict[tissue][system][0]/fc_dict[tissue][system][1])*100
            writer.writerow(write_row)

percent_FC_bigger_than_1()
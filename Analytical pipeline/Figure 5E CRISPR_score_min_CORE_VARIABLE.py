import csv
import pandas as pd

'''------------------------------------Source files required to run this code----------------------------------------'''
autophagy = pd.read_csv('core_variable_differential_V8_AUT_5TPM_core_consortium_trimmed.csv', sep=',')
chap_genes = pd.read_csv('core_variable_differential_V8_Chap_5TPM_core_consortium_trimmed.csv', sep = ',')
crispr = pd.read_csv('DepMap_2022.csv',sep = ',')
ups_genes = pd.read_csv('core_variable_differential_V8_UPS_5TPM_core_consortium_trimmed.csv')
syn_genes = pd.read_csv('core_variable_differential_V8_SYN_5TPM_core_consortium_trimmed.csv')
'''------------------------------------------------------------------------------------------------------------------'''


'''-------------This function extracts minimal CRISPR score for each gene from array of cancer cultures--------------'''

def mininmal_growth():
    names = {}
    deg = []
    fold = []
    syn = []

    for linNum, line in autophagy.iterrows():
        names[line['Name'].strip()] = []
        if line['Name'] not in deg:
            deg.append(line['Name'].strip())

    for lineNum, line in ups_genes.iterrows():
        if line['Name'] not in deg:
            deg.append(line['Name'].strip())

    for lineNum, line in chap_genes.iterrows():
        if line['Name'] not in fold:
            fold.append(line['Name'].strip())

    for lineNum, line in syn_genes.iterrows():
        if line['Name'] not in fold:
            syn.append(line['Name'].strip())

    fieldnames = ['Gene name', 'Degradation', 'Min']

    print(deg, 'Degradation')
    print(fold, 'Folding')
    print(names, 'names')
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        write_row = {}

        for rowNum, row in crispr.iterrows():
            write_row['Gene name'] = row['Name']

            write_row['Degradation'] = ''
            if row['Name'] in deg:
                write_row['Degradation'] = 'Yes'

            if row['Name'] in fold:
                write_row['Degradation'] = 'No'

            if row['Name'] in syn:
                write_row['Degradation'] = 'SYN'

            if write_row['Degradation'] != '':

                val_arr = []
                for item in row.to_dict():
                    if item not in ['Name']:
                        if row[item] == row[item]:
                            val_arr.append(float(row[item]))
                write_row['Min'] = min(val_arr)
                writer.writerow(write_row)





output_file = 'DepMap_CRISPR_min_score_deg_vs_fold.csv'

mininmal_growth()

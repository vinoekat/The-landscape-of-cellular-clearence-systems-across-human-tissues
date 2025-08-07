import csv
import pandas as pd

'''------------------------------------Source files required to run this code----------------------------------------'''
autophagy = pd.read_csv('core_variable_differential_V8_AUT_5TPM_core_consortium_trimmed.csv', sep=',')
protein_coding = pd.read_csv('Protein coding expressing tissues 5TPM.csv', sep = ',')
crispr = pd.read_csv('DepMap_2022.csv',sep = ',')
ups_genes = pd.read_csv('core_variable_differential_V8_UPS_5TPM_core_consortium_trimmed.csv')
'''------------------------------------------------------------------------------------------------------------------'''


'''-------------This function extracts minimal CRISPR score for each gene from array of cancer cultures--------------'''
def mininmal_growth():
    names = {}
    autophages = []
    UPS = []
    prot_coding = []

    for linNum, line in autophagy.iterrows():
        names[line['Name'].strip()] = []
        autophages.append(line['Name'].strip())

    for lineNum, line in ups_genes.iterrows():
        if line['Name'] not in UPS:
            UPS.append(line['Name'])

    for lineNum, line in protein_coding.iterrows():
        names[line['Protein Name'].strip()] = []
        prot_coding.append(line['Protein Name'].strip())

    fieldnames = ['Gene name', 'Autophage', 'Min']

    print(autophages,'autophages')
    print(prot_coding,'prot')
    print(names,'names')
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        write_row={}

        for rowNum, row in crispr.iterrows():
            write_row['Gene name'] = row['Name']

            write_row['Autophage'] = ''
            if row['Name'] in autophages:
                write_row['Autophage'] = 'Yes'

            if row['Name'] in prot_coding and row['Name'] not in autophages\
                    and row['Name'] not in UPS:
                write_row['Autophage'] = 'No'

            if row['Name'] in UPS:
                write_row['Autophage'] = 'UPS'

            if write_row['Autophage'] != '':

                val_arr = []
                for item in row.to_dict():
                    if item not in ['Name']:
                        if row[item] == row[item]:
                            val_arr.append(float(row[item]))
                write_row['Min'] = min(val_arr)
                writer.writerow(write_row)


output_file = 'DepMap_CRISPR_min_score_consortium.csv'

mininmal_growth()

import csv
import pandas as pd
import numpy as np
import statistics
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import math

'''--------------------------------Source files required to run this code--------------------------------------------'''

# use UPS or AUT for the respective plots in figure5D
system = 'UPS'
core_variable_tissue_count = pd.read_csv('core_variable_differential_V8_' + system +'_5TPM_core_consortium_trimmed.csv')
output_file = 'median_expression_core_var_V8_' + system +'_bellow_5TPM_dropped_consortium.csv'
median_expression_file = pd.read_csv('median_expression_per_tissue_GTEx_2023.csv', delimiter=',')
'''------------------------------------------------------------------------------------------------------------------'''

'''-----------this script extracts median expression of core and variable genes in all tissues-----------------------'''

core_list = []
variable_list = []

fd = ['gene name', 'median expression', 'group']

writer = csv.DictWriter(
    open(output_file, 'w', newline=''), fieldnames=fd)
writer.writeheader()
write_row = {}

for rowNum, row in core_variable_tissue_count.iterrows():
    row = row.to_dict()
    if row['core/var'] == 'core':
        core_list.append(row['ENSG'])
    else:
        variable_list.append(row['ENSG'])

print(len(core_list), len(variable_list))

for rowNum, row in median_expression_file.iterrows():
    row = row.to_dict()
    if row['ENSG'].strip()[0:15] in core_list:
        temp_list = []
        for item in row:

            if item not in ['ENSG', 'Name']:
                if float(row[item]) >= 5:
                    # if float(row[item]) == 0:
                    #     write_row['median expression'] = 0
                    if float(row[item]) > 100:
                        write_row['median expression'] = 100
                    else:
                        write_row['median expression'] = row[item]

                    write_row['gene name'] = row['Name']
                    write_row['group'] = 'core'
                    writer.writerow(write_row)
    if row['ENSG'].strip()[0:15] in variable_list:
        temp_list=[]
        for item in row:

            if item not in ['ENSG', 'Name']:
                if float(row[item]) >= 5:
                    # if float(row[item]) == 0:
                    #     write_row['median expression'] = 0
                    if float(row[item]) > 100:
                        write_row['median expression'] = 100
                    else:

                        write_row['median expression'] = row[item]

                    write_row['gene name'] = row['Name']
                    write_row['group'] = 'variable'
                    writer.writerow(write_row)


def mann_whitney():

    log_median_core_variable = pd.read_csv(output_file, delimiter=',')

    array1 = []
    array2 = []


    for rowNUM, row in log_median_core_variable.iterrows():
        row = row.to_dict()
        #print(row)
        if row['group'] == 'core':
            array1.append(row['median expression'])
        if row['group'] == 'variable':
            array2.append(row['median expression'])
    print(len(array1), array1)
    print(len(array2), array2)

    print(stats.mannwhitneyu(array1, array2, alternative='two-sided'))


    print(stats.ks_2samp(array1, array2,alternative = 'two-sided'))


mann_whitney()

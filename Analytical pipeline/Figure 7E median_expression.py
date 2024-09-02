import csv
import scipy.stats as stats
import statistics
import pandas as pd

'''--------------------------------Source files required to run this code--------------------------------------------'''
autophagy_file = pd.read_csv('core_variable_differential_V8_AUT_5TPM_core_consortium_trimmed.csv')
ups_file = pd.read_csv('core_variable_differential_V8_UPS_5TPM_core_consortium_trimmed.csv')
chap_file = pd.read_csv('core_variable_differential_V8_Chap_5TPM_core_consortium_trimmed.csv')
syn_file = pd.read_csv('core_variable_differential_V8_SYN_5TPM_core_consortium_trimmed.csv')
medians_file = pd.read_csv('median_expression_per_tissue_GTEx_2023.csv')
'''------------------------------------------------------------------------------------------------------------------'''

'''-----------this script extracts median expression of core and variable genes in all tissues-----------------------'''


degradation_list = []
folding_list = []
synthesis_list = []

for rowNUm, row in autophagy_file.iterrows():
    row = row.to_dict()
    if row['ENSG'] not in degradation_list:
        degradation_list.append(row['ENSG'])

for rowNUm, row in ups_file.iterrows():
    row = row.to_dict()
    if row['ENSG'] not in degradation_list:
        degradation_list.append(row['ENSG'])

for rowNUm, row in chap_file.iterrows():
    row = row.to_dict()
    if row['ENSG'] not in folding_list:
        folding_list.append(row['ENSG'])

for rowNUm, row in syn_file.iterrows():
    row = row.to_dict()
    if row['ENSG'] not in synthesis_list:
        synthesis_list.append(row['ENSG'])

tissue_medians = {}
for item in medians_file.columns.tolist():
    if item not in ['ENSG', 'Name']:
        tissue_medians[item.strip()] = [[], [], []]


gene_values_deg = [] # 416*32 values
gene_values_fold = []
gene_values_syn = []


median_values_deg = [] # 32 values (median of median), system per tissue
median_values_fold = []
median_values_syn = []



'''write file with median of a gene per tissue ENSG/Name/Median value/System'''
fd = ['ENSG', 'Median value', 'Tissue', 'System']
writer = csv.DictWriter(open('Median RNA expression per gene deg vs fold.csv', 'w', newline=''), fieldnames=fd)
writer.writeheader()
write_row = {}

for lineNum, line in medians_file.iterrows():
    line = line.to_dict()
    if line['ENSG'].strip()[0:15] in degradation_list:
        temp_list = []
        for item in line:
            if item not in ['ENSG', 'Name'] and float(line[item]) >= 5:

                write_row['ENSG'] = line['ENSG'].strip()[0:15]
                write_row['Median value'] = line[item]
                write_row['Tissue'] = item
                write_row['System'] = 'Degradation'
                writer.writerow(write_row)

                temp_list.append(float(line[item]))

                gene_values_deg.append(float(line[item]))
                tissue_medians[item][0].append(float(line[item]))

    if line['ENSG'].strip()[0:15] in folding_list:
        temp_list = []
        for item in line:
            if item not in ['ENSG', 'Name'] and float(line[item]) >= 5:

                write_row['ENSG'] = line['ENSG'].strip()[0:15]
                write_row['Median value'] = line[item]
                write_row['Tissue'] = item
                write_row['System'] = 'Folding'
                writer.writerow(write_row)

                temp_list.append(float(line[item]))

                gene_values_fold.append(float(line[item]))
                tissue_medians[item][1].append(float(line[item]))

    if line['ENSG'].strip()[0:15] in synthesis_list:
        temp_list = []
        for item in line:
            if item not in ['ENSG', 'Name'] and float(line[item]) >= 5:

                write_row['ENSG'] = line['ENSG'].strip()[0:15]
                write_row['Median value'] = line[item]
                write_row['Tissue'] = item
                write_row['System'] = 'Synthesis'
                writer.writerow(write_row)

                temp_list.append(float(line[item]))

                gene_values_syn.append(float(line[item]))
                tissue_medians[item][2].append(float(line[item]))



fd = ['Tissue', 'Median value', 'System']
writer = csv.DictWriter(open('Median RNA expression of a system per tissue.csv', 'w', newline=''), fieldnames=fd)
writer.writeheader()
write_row = {}

systems_list = ['Degradation', 'Folding', 'Synthesis']

for system in systems_list:
    if system == 'Degradation':
        for tissue in tissue_medians:
            write_row['Tissue'] = tissue
            write_row['Median value'] = statistics.median(tissue_medians[tissue][0])
            write_row['System'] = 'Degradation'
            writer.writerow(write_row)
            median_values_deg.append(statistics.median(tissue_medians[tissue][0]))
    if system == 'Folding':
        for tissue in tissue_medians:
            write_row['Tissue'] = tissue
            write_row['Median value'] = statistics.median(tissue_medians[tissue][1])
            write_row['System'] = 'Folding'
            writer.writerow(write_row)
            median_values_fold.append(statistics.median(tissue_medians[tissue][1]))
    if system == 'Synthesis':
        for tissue in tissue_medians:
            write_row['Tissue'] = tissue
            write_row['Median value'] = statistics.median(tissue_medians[tissue][2])
            write_row['System'] = 'Synthesis'
            writer.writerow(write_row)
            median_values_syn.append(statistics.median(tissue_medians[tissue][2]))

print('all values deg vs fold', stats.mannwhitneyu(gene_values_deg, gene_values_fold), len(gene_values_deg), len(gene_values_fold))


print('tissue medians deg vs fold', stats.mannwhitneyu(median_values_deg, median_values_fold), len(median_values_deg), len(median_values_fold))


print('all values deg vs syn', stats.mannwhitneyu(gene_values_deg, gene_values_syn), len(gene_values_deg), len(gene_values_syn))
print('tissue medians deg vs syn', stats.mannwhitneyu(median_values_deg, median_values_syn), len(median_values_deg), len(median_values_syn))


print('all values syn vs fold', stats.mannwhitneyu(gene_values_syn, gene_values_fold), len(gene_values_syn), len(gene_values_fold))
print('tissue medians syn vs fold', stats.mannwhitneyu(median_values_syn, median_values_fold), len(median_values_syn), len(median_values_fold))


print(median_values_deg)


deg_wilcoxon = []
fold_wilcoxon = []


for tissue in tissue_medians:

    deg_wilcoxon.append(statistics.median(tissue_medians[tissue][0]))
    fold_wilcoxon.append(statistics.median(tissue_medians[tissue][1]))



print('tissue medians deg vs fold', stats.wilcoxon(deg_wilcoxon, fold_wilcoxon))

print('tissue medians deg vs fold', stats.wilcoxon(median_values_deg, median_values_fold))

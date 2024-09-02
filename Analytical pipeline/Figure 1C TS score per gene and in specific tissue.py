import csv
import pickle
import scipy.stats as stats
import statistics
import pandas as pd

'''-----------------------------------Files required to run this code------------------------------------------------'''

aut_file = pd.read_csv('core_variable_differential_V8_AUT_5TPM_core_consortium_trimmed.csv')
ups_file = pd.read_csv('core_variable_differential_V8_UPS_5TPM_core_consortium_trimmed.csv')
protein_file = pd.read_csv('protein coding 2023.csv')
ts_score_file = pd.read_csv('TS scores.csv') # to obtain this file please refer to the suplementary tables section of the
# "A Quantitative Proteome Map of the Human Body" PMID: 32916130 - article, tableS2G.

'''------------------------------------------------------------------------------------------------------------------'''
''' If the TS score  ≥  4 and there were no other tissues with scores in the interval (2.5, 4), the protein was
 classified as tissue-specific (in total 1,595). If there was at least one TS score  ≥  2.5 but the protein was
  not tissue-specific, the protein was classified as tissue-enriched-but-not-specific (in total 3,967). If a protein
   was observed in at least one sample before imputation in each tissue type and all its TS scores were non-NAs and 
   less than 2, the protein was classified as house-keeping (in total 1,565). The rest proteins belonged to the “others” 
   category (in total 5,500). As a requirement, at least three non-NA tissue scores were needed for a protein
    to be categorized as tissue-enriched The protein enrichment information is summarized in Table S2. snyder, cell 2020'''
'''------------------------------------------------------------------------------------------------------------------'''

list_of_files = [aut_file, ups_file, protein_file]

aut_dict = {}
ups_dict = {}
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

for rowNUm, row in protein_file.iterrows():
    row = row.to_dict()
    protein_dict[row['ENSG']] = []
    protein_dict[row['ENSG']].append(row['Name'])


'''Extract TS scores per gene in each system'''
def TS_score_per_gene(tissue):

    ts_score_dict = {}
    '''loop creates ts_score_dict that has the following built:
    'ENSG':[median TS score, number of scored tissues, type - housekeeping/enriched/tissue-specific/other]
    dict contains values that had only one repeat (i.e. '-4.44;NA_one_rep_in_raw') because of Snyder's definition
     of housekeeping'''
    for lineNum, line in ts_score_file.iterrows():
        line = line.to_dict()
        #print(line)
        '''get the values per gene'''
        temp_list = []
        for item in line:
            if item not in ['ensembl_id', 'entrez_id', 'hgnc_name', 'hgnc_symbol'] and line[item] == line[item]:
                temp_value = line[item].strip().split(';')
                #print(len(temp_value))
                temp_list.append(float(temp_value[0]))
                # if len(temp_value) == 1:
                #     temp_list.append(float(temp_value[0]))
                #print(temp_value)

        '''tissue specific if has only 1 score above 2.5 and that score is above 4'''
        '''define the type of values'''
        tissue_specific = False
        housekeeping = False
        enriched = False

        high_scores = []
        low_scores = []
        for score in temp_list:
            if score >= 2.5:
                high_scores.append(score)
            if score < 2:
                low_scores.append(score)

        #print(line[tissue])
        if line[tissue] == line[tissue]:
            tissue_value = line[tissue].split(';')
            if len(high_scores) == 1 and high_scores[0] >= 4 and float(tissue_value[0].strip()) >= 4:
                tissue_specific = True
            if len(temp_list) >= 3:
                if len(high_scores) == 1 and high_scores[0] < 4 and float(tissue_value[0].strip()) >= 2.5:
                    enriched = True
                if len(temp_list) >= 3 and len(high_scores) > 1 and float(tissue_value[0].strip()) >= 2.5:
                    enriched = True
            if len(low_scores) == 32 and float(tissue_value[0].strip()) < 2:
                housekeeping = True

            '''fill the dict'''
            if len(temp_list) != 0:
                ts_score_dict[line['ensembl_id']] = []
                #ts_score_dict[line['ensembl_id']].append(statistics.median(temp_list))  # median ts score
                ts_score_dict[line['ensembl_id']].append(line[tissue])  # ts score in the tissue
                ts_score_dict[line['ensembl_id']].append(len(temp_list)) # number of scored tissues
                if tissue_specific and not enriched:
                    ts_score_dict[line['ensembl_id']].append('Tissue-specific')
                if enriched:
                    ts_score_dict[line['ensembl_id']].append('Tissue-enriched')
                if housekeeping:
                    ts_score_dict[line['ensembl_id']].append('Housekeeping')
                if not tissue_specific and not enriched and not housekeeping:
                    ts_score_dict[line['ensembl_id']].append('Other')

        # print(temp_list)
        #print(ts_score_dict[line['ensembl_id']])

    fd = ['ENSG', 'Name', 'TS score in tissue', 'Num of scored tissues', 'core/var', 'Class', 'TS type', 'System']
    writer = csv.DictWriter(open('TS score per gene AUT and UPS in ' + tissue + '.csv', 'w', newline=''), fieldnames=fd)
    writer.writeheader()
    write_row = {}

    for ensg in ts_score_dict:
        if ensg in aut_dict:
            write_row['ENSG'] = ensg
            write_row['Name'] = aut_dict[ensg][0]
            write_row['TS score in tissue'] = ts_score_dict[ensg][0]
            write_row['Num of scored tissues'] = ts_score_dict[ensg][1]
            write_row['TS type'] = ts_score_dict[ensg][2]

            write_row['Class'] = aut_dict[ensg][1]
            write_row['core/var'] = aut_dict[ensg][2]
            write_row['System'] = 'AUT'
            writer.writerow(write_row)

        if ensg in ups_dict:
            write_row['ENSG'] = ensg
            write_row['Name'] = ups_dict[ensg][0]
            write_row['TS score in tissue'] = ts_score_dict[ensg][0]
            write_row['Num of scored tissues'] = ts_score_dict[ensg][1]
            write_row['TS type'] = ts_score_dict[ensg][2]

            write_row['Class'] = ups_dict[ensg][1]
            write_row['core/var'] = ups_dict[ensg][2]
            write_row['System'] = 'abUPS'
            writer.writerow(write_row)

        if ensg in protein_dict and ensg not in combined_list:
            write_row['ENSG'] = ensg
            write_row['Name'] = protein_dict[ensg][0]
            write_row['TS score in tissue'] = ts_score_dict[ensg][0]
            write_row['Num of scored tissues'] = ts_score_dict[ensg][1]
            write_row['TS type'] = ts_score_dict[ensg][2]

            write_row['Class'] = 'Protein-coding'
            write_row['core/var'] = 'Protein-coding'
            write_row['System'] = 'Protein-coding'
            writer.writerow(write_row)

'''Calculate percent of each protein type per tissue in 3 systems'''
def percent_TS_score(tissue):

    ts_score_systems_file = pd.read_csv('TS score per gene AUT and UPS in ' + tissue + '.csv')

    ts_percent_dict = {'AUT': {'Tissue-specific': 0, 'Tissue-enriched': 0, 'Housekeeping': 0, 'Other': 0, 'Total':0},
                       'abUPS': {'Tissue-specific': 0, 'Tissue-enriched': 0, 'Housekeeping': 0, 'Other': 0, 'Total': 0},
                       'Protein-coding': {'Tissue-specific': 0, 'Tissue-enriched': 0, 'Housekeeping': 0, 'Other': 0, 'Total': 0}}

    for lineNum, line in ts_score_systems_file.iterrows():
        line = line.to_dict()
        #print(ts_percent_dict[line['System']][line['Type']])
        ts_percent_dict[line['System']][line['TS type']] += 1
        if line['System'] == 'AUT':
            ts_percent_dict['AUT']['Total'] += 1
        if line['System'] == 'abUPS':
            ts_percent_dict['abUPS']['Total'] += 1
        if line['System'] == 'Protein-coding':
            ts_percent_dict['Protein-coding']['Total'] += 1

    fd = ['System', 'TS score percent', 'TS type']
    writer = csv.DictWriter(open('AUT and UPS TS score percent in ' + tissue + '.csv', 'w', newline=''), fieldnames=fd)
    writer.writeheader()
    write_row = {}

    for system in ts_percent_dict:
        for type in ts_percent_dict[system]:
            if type != 'Total':
                write_row['System'] = system
                write_row['TS score percent'] = (ts_percent_dict[system][type]/ts_percent_dict[system]['Total'])*100
                write_row['TS type'] = type
                writer.writerow(write_row)


tissue_list = ts_score_file.columns.tolist()
tissue_list.pop(0)
tissue_list.pop(0)
tissue_list.pop(0)
tissue_list.pop(0)

print(len(tissue_list))
for tis in tissue_list:
    TS_score_per_gene(tis)
    percent_TS_score(tis)


'''Combine tissue files into 1 file for box plot'''
fd = ['System', 'TS score percent', 'TS type', 'Tissue']
writer = csv.DictWriter(open('TS score percent per tissue AUT UPS and protein coding for box plot.csv', 'w', newline=''), fieldnames=fd)
writer.writeheader()
write_row = {}
systems = ['AUT', 'abUPS', 'Protein-coding']
for tissue in tissue_list:
    file = pd.read_csv('AUT and UPS TS score percent in ' + tissue + '.csv')
    for system in systems:
        temp_sum = 0
        for rowNum, row in file.iterrows():
            row = row.to_dict()
            if row['System'] == system:
                if row['TS type'] == 'Tissue-specific' or row['TS type'] == 'Tissue-enriched':
                    temp_sum += float(row['TS score percent'])
                if row['TS type'] == 'Housekeeping':
                    write_row['System'] = row['System']
                    write_row['TS type'] = row['TS type']
                    write_row['TS score percent'] = row['TS score percent']
                    write_row['Tissue'] = tissue
                    writer.writerow(write_row)


        write_row['System'] = system
        write_row['TS type'] = 'Tissue-specific/Tissue-enriched'
        write_row['TS score percent'] = temp_sum
        write_row['Tissue'] = tissue
        writer.writerow(write_row)
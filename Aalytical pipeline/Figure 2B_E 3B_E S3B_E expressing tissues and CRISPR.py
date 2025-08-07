import csv
import pandas as pd
import statistics
from scipy import stats

'''-----------------------------------Files required to run this code------------------------------------------------'''

#Files needed for figure 2B/E
#manually_curated_genes_to_stages = pd.read_csv('consortium core var AUT w DC w dups ESCRT no lys_ref.csv')

#files needed for figure 3B/E
#manually_curated_genes_to_stages = pd.read_csv('consortium core var UPS with diseases with dups Class Ub and UBL unified.csv')

#files needed for figure S2B/E
#manually_curated_genes_to_stages = pd.read_csv('core_variable_diff_V8_UPS PROTEASOME_5TPM.csv')


# Common files
median_expression = pd.read_csv('median_expression_per_tissue_GTEx_2023.csv')
CRISPR = pd.read_csv('DepMap_2022.csv')
protein_coding = pd.read_csv('protein coding 2023.csv')


'''---------------------------------------------Conditions-----------------------------------------------------------'''
# # use this to run for figure 2 B/E
system = 'UPS'
comment = 'Ub and UBL changed names'

# # use this to run for figure 3 B/E
# system = 'AUT'
# comment = 'ESCRT changed names'

# # use this to run for figure S2 B/E
# system = 'UPS'
# comment = 'PROTEASOME changed names'

#common conditions
TPM = 5
Number_of_tissue = 34
'''------------------------------------------------------------------------------------------------------------------'''


'''----------------------------------------Obtain gene-category dict-------------------------------------------------'''
def group_dict():
    autophagy_category_gene_dict = {}
    for rowNum, row in file_to_dict.iterrows():
        row = row.to_dict()

        if row[category] in autophagy_category_gene_dict:
            autophagy_category_gene_dict[row[category]].append(row[gene_symbol])
        if row[category] not in autophagy_category_gene_dict:
            autophagy_category_gene_dict[row[category]] = [row[gene_symbol]]

    return autophagy_category_gene_dict


'''--------------Calculate minimal CRISPR score and percent expressing tissues for each gene in each group-----------'''
def median_gene_expression_tissue_count():

    fd = ['gene name', 'category', 'abs number of expressing tissues', '% abs number of expressing tissues']

    writer = csv.DictWriter(open('Expressing tissues '\
                                 + str(TPM) + 'TPM (ABC) ' + system +' ' + category + ' consortium '\
                                 + comment +'.csv', 'w', newline=""),
                            fieldnames=fd)
    writer.writeheader()
    write_row = {}

    for rowNum, row in median_expression.iterrows():
        row = row.to_dict()
        for stage in autophagy_stage_dict:
            if row['Name'].strip() in autophagy_stage_dict[stage]:
                gene_expression = 0
                for item in row:
                    if item not in ['ENSG', 'Name']:
                        if float(row[item]) >= TPM:
                            gene_expression += 1

                write_row['gene name'] = row['Name']
                write_row['category'] = stage
                write_row['abs number of expressing tissues'] = gene_expression
                write_row['% abs number of expressing tissues'] = (gene_expression/Number_of_tissue) * 100
                writer.writerow(write_row)


def min_CRISPR_score_tissue_count():
    fd = ['gene name', 'category', 'min_CRISPR_score']

    writer = csv.DictWriter(open('min CRISPR score '\
                                 + str(TPM) + 'TPM (ABC) ' + system +' ' + category + ' consortium '\
                                 + comment +'.csv', 'w', newline=""),
                            fieldnames=fd)
    writer.writeheader()
    write_row = {}

    for rowNum, row in CRISPR.iterrows():
        row = row.to_dict()
        for stage in autophagy_stage_dict:
            if row['Name'] in autophagy_stage_dict[stage]:
                temp_line_list = []
                for item in row:
                    if item not in ['Name']:
                        if row[item] == row[item]:
                            temp_line_list.append(float(row[item]))
                min_CRISPER = min(temp_line_list)
                write_row['gene name'] = row['Name']
                write_row['category'] = stage
                write_row['min_CRISPR_score'] = min_CRISPER
                writer.writerow(write_row)


def mann_whitney():


    CRISPR = pd.read_csv('min CRISPR score '\
                                 + str(TPM) + 'TPM (ABC) ' + system +' ' + category + ' consortium '\
                                 + comment +'.csv', delimiter=',')
    tissue_expression = pd.read_csv('Expressing tissues '\
                                 + str(TPM) + 'TPM (ABC) ' + system +' ' + category + ' consortium '\
                                 + comment +'.csv')


    pvalue_array = []
    group_name_dict = {}
    for itemNum, item in tissue_expression.iterrows():
        item = item.to_dict()
        group_name_dict[item['category']] = []
    print('name dict', group_name_dict)

    for item in group_name_dict:
        tis_exp_main = []
        tis_exp_rest = []
        for rowN, row in tissue_expression.iterrows():
            row = row.to_dict()

            if row['category'] == item:
                tis_exp_main.append(row['abs number of expressing tissues'])
            elif row['category'] != item:
                tis_exp_rest.append(row['abs number of expressing tissues'])

        print('expression', item, stats.mannwhitneyu(tis_exp_main, tis_exp_rest))
        k = stats.mannwhitneyu(tis_exp_main, tis_exp_rest)
        temp_list = []
        temp_list.append('expression')
        temp_list.append(item)
        temp_list.append(str(k[1]))
        pvalue_array.append(temp_list)
        print(len(tis_exp_main), len(tis_exp_rest))

    for item in group_name_dict:
        tis_main = []
        tis_rest = []
        for rowN, row in CRISPR.iterrows():
            row = row.to_dict()

            if row['category'] == item:
                tis_main.append(row['min_CRISPR_score'])
            elif row['category'] != item:
                tis_rest.append(row['min_CRISPR_score'])
        print('CRISPR', item, stats.mannwhitneyu(tis_main, tis_rest))
        k = stats.mannwhitneyu(tis_main, tis_rest)
        temp_list = []
        temp_list.append('CRISPR')
        temp_list.append(item)
        temp_list.append(str(k[1]))
        pvalue_array.append(temp_list)
        print(len(tis_main), len(tis_rest))

    with open('pvalues new both ' + system + '.csv', 'w', newline = '') as out_object:
        out = csv.writer(out_object)

        for item in pvalue_array:
            out.writerow(item)


# data structure
# to get file for stages change key_column to 'stage'
file_to_dict = manually_curated_genes_to_stages
category = 'Class'
gene_symbol = 'Name'
autophagy_stage_dict = group_dict()

# functions
median_gene_expression_tissue_count()
differential_gene_expression_tissue_count()
min_CRISPR_score_tissue_count()
mann_whitney()


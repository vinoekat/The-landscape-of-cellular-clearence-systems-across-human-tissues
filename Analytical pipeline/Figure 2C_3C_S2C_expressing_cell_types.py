import csv
import scipy.stats as stats
import statistics
import pandas as pd

# this code runs for several figures
# use sections indicated for the specific figure and silence sections indicated for other figures
# note there are section at the botto of the code as well
'''-----------------------------------Files required to run this code------------------------------------------------'''
# use this for figure 2C
class_ups_file = pd.read_csv('all GTex names consortium file UPS with dup Class UB and UBL unified.csv')

# use this for figure 3C
class_aut_file = pd.read_csv('all ENSG GTex names consortium file AUT with dup ESCRT.csv')

# use this for figure S2C
class_ups_prot_file = pd.read_csv('all GTex names consortium file UPS with dup Class PROTEASOME.csv')



threshold = 0.15
tissue_list = ['Bladder', 'Blood', 'Bone_Marrow', 'endothelial', 'epithelial', 'Eye', 'Fat', 'Heart', 'immune', 'Kidney',\
               'Large_Intestine', 'Liver', 'Lung', 'Lymph_Node', 'Mammary', 'Muscle', 'Pancreas', 'Prostate', 'Salivary_Gland',\
               'Skin', 'Small_intestine', 'Spleen', 'stromal', 'Thymus', 'Tongue', 'Trachea', 'Uterus', 'Vasculature']


def expression_score_per_tissue(file, system, comment):
    '''Function receives a list of single cell files for various tissues, it calculates percent expressing cell types
    above the threshold in each tissue, then summarizes the result in a file'''

    gene_ensg_dict = {}

    for lineNum, line in file.iterrows():
        line = line.to_dict()
        if line['ENSG'] in gene_ensg_dict:
            gene_ensg_dict[line['ENSG']].append([line['Name'], line['Class']])
        if line['ENSG'] not in gene_ensg_dict:
            gene_ensg_dict[line['ENSG']] = [[line['Name'], line['Class']]]


    fd = ['ENSG', 'Name', 'number of expressing subsets', 'Expression score', 'Class', '% Expression score', 'Tissue']
    writer = csv.DictWriter(open(system + '_expression_score_all_tissues_with_dups_' + comment + '.csv', 'w',\
                                 newline=''), fieldnames=fd)

    writer.writeheader()
    write_row = {}

    for tissue in tissue_list:
        single_cell_file = pd.read_csv('ENSG normalized_Expression_' + tissue + '.csv')
        number_of_subsets = len(single_cell_file.columns.tolist()) - 2

        for lineNum, line in single_cell_file.iterrows():
            line = line.to_dict()
            if line['ENSG'] == line['ENSG'] and line['ENSG'] in gene_ensg_dict:
                subset_expression = 0
                for item in line:
                    if item not in ['ENSG', 'Name']:
                        if line[item] == line[item] and float(line[item]) >= threshold:
                            subset_expression += 1

                for name_class_list in gene_ensg_dict[line['ENSG']]:
                    write_row['ENSG'] = line['ENSG']
                    write_row['Name'] = line['Name']
                    write_row['number of expressing subsets'] = subset_expression
                    write_row['Expression score'] = int(subset_expression) / int(number_of_subsets)
                    write_row['Class'] = name_class_list[1]
                    write_row['% Expression score'] = (int(subset_expression) / int(number_of_subsets)) * 100
                    write_row['Tissue'] = tissue
                    writer.writerow(write_row)


# use this for figure 2C
expression_score_per_tissue(class_ups_file, system='UPS', comment='Ub UBL unified')

# use this for figure 3C
expression_score_per_tissue(class_aut_file, system='AUT', comment='ESCRT')

# use this for figure S2C
expression_score_per_tissue(class_ups_prot_file, system='UPS', comment='PROTEASOME')

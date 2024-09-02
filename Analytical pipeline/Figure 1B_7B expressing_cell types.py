import csv
import scipy.stats as stats
import statistics
import pandas as pd


'''This code produces data for 2 figures (Fig1B, Fig7B), in order to get the results for specific figure silence the 
section relating to the other figure. # - sections for silencing. To run the code copy source files into the same folder
 as the code. 
 
 Files needed to run the code:
 all GTex names conslist ENSG UniprotID AUT no dups w our genes.csv
 all GTex names conslist ENSG UniprotID UPS no dups w our genes.csv
 consortium file with ENSG based on UniprotID SYN no dup.csv
 ENSG new gene list filtered to relevant chaperones.csv
 protein coding 2023.csv

And the 28 files found in the "Single cell source files folder" named "ENSG normalized_Expression_(tissue name).csv". 
 '''



threshold = 0.15
tissue_list = ['Bladder', 'Blood', 'Bone_Marrow', 'endothelial', 'epithelial', 'Eye', 'Fat', 'Heart', 'immune', 'Kidney',\
               'Large_Intestine', 'Liver', 'Lung', 'Lymph_Node', 'Mammary', 'Muscle', 'Pancreas', 'Prostate', 'Salivary_Gland',\
               'Skin', 'Small_intestine', 'Spleen', 'stromal', 'Thymus', 'Tongue', 'Trachea', 'Uterus', 'Vasculature']
autophagy_file = pd.read_csv('all GTex names conslist ENSG UniprotID AUT no dups w our genes.csv')
ups_file = pd.read_csv('all GTex names conslist ENSG UniprotID UPS no dups w our genes.csv')
chaperone_file = pd.read_csv('ENSG new gene list filtered to relevant chaperones.csv')
synthesis_file = pd.read_csv('consortium file with ENSG based on UniprotID SYN no dup.csv')
protein_coding = pd.read_csv('protein coding 2023.csv')



def system_ensg_name_dict(file):
    '''function receives file and returns the following dict:
    receives: line - ENSG, Name
    returns: ENSG:Name
             Name: ENSG'''

    dict = {}
    for rowNum, row in file.iterrows():
        row = row.to_dict()
        dict[row['Name']] = row['ENSG']
        if row['ENSG'] == row['ENSG'] and row['ENSG'] != '':
            dict[row["ENSG"]] = row['Name']

    return dict

aut_dict = system_ensg_name_dict(autophagy_file)
ups_dict = system_ensg_name_dict(ups_file)
chap_dict = system_ensg_name_dict(chaperone_file)
syn_dict = system_ensg_name_dict(synthesis_file)
prot_dict = system_ensg_name_dict(protein_coding)

degradation_dict = {}

for item in aut_dict:
    degradation_dict[item] = aut_dict[item]
for item in ups_dict:
    degradation_dict[item] = ups_dict[item]


def expression_score(gene_dict1, gene_dict2, gene_dict3, system1, system2, system3):
    fd = ['ENSG', 'Name', 'number of expressing subsets', 'expression score', '% expression score', 'tissue', 'system']
    writer = csv.DictWriter(open('Expression_score_all_tissues_no_dups_' + system1 + '_' + system2 + '_' + system3 + '.csv',
                                 'w', newline=''), fieldnames=fd)
    writer.writeheader()
    write_row = {}

    for tissue in tissue_list:
        single_cell_file = pd.read_csv('ENSG normalized_Expression_' + tissue + '.csv')
        number_of_subsets = len(single_cell_file.columns.tolist()) - 2
        for lineNum, line in single_cell_file.iterrows():
            line = line.to_dict()
            expressing_subsets = 0
            for item in line:
                if item not in ['ENSG', 'Name']:
                    if line[item] == line[item] and float(line[item]) >= threshold:
                        expressing_subsets += 1

            if line['ENSG'] in gene_dict1:
                write_row['ENSG'] = line['ENSG']
                write_row['Name'] = gene_dict1[line['ENSG']]
                write_row['number of expressing subsets'] = expressing_subsets
                write_row['expression score'] = int(expressing_subsets) / int(number_of_subsets)
                write_row['tissue'] = tissue
                write_row['system'] = system1
                write_row['% expression score'] = (int(expressing_subsets) / int(number_of_subsets)) * 100
                writer.writerow(write_row)

            if line['ENSG'] in gene_dict2:
                write_row['ENSG'] = line['ENSG']
                write_row['Name'] = gene_dict2[line['ENSG']]
                write_row['number of expressing subsets'] = expressing_subsets
                write_row['expression score'] = int(expressing_subsets) / int(number_of_subsets)
                write_row['tissue'] = tissue
                write_row['system'] = system2
                write_row['% expression score'] = (int(expressing_subsets) / int(number_of_subsets)) * 100
                writer.writerow(write_row)



                if line['ENSG'] in gene_dict3:
                    if system3 =='Protein-coding':
                        if line['ENSG'] not in aut_dict and line['ENSG'] not in ups_dict and line['ENSG'] in prot_dict:
                            write_row['ENSG'] = line['ENSG']
                            write_row['Name'] = prot_dict[line['ENSG']]
                            write_row['number of expressing subsets'] = expressing_subsets
                            write_row['expression score'] = int(expressing_subsets) / int(number_of_subsets)
                            write_row['tissue'] = tissue
                            write_row['system'] = 'Protein-coding'
                            write_row['% expression score'] = (int(expressing_subsets) / int(number_of_subsets)) * 100
                            writer.writerow(write_row)
                    else:
                        write_row['ENSG'] = line['ENSG']
                        write_row['Name'] = gene_dict3[line['ENSG']]
                        write_row['number of expressing subsets'] = expressing_subsets
                        write_row['expression score'] = int(expressing_subsets) / int(number_of_subsets)
                        write_row['tissue'] = tissue
                        write_row['system'] = system3
                        write_row['% expression score'] = (int(expressing_subsets) / int(number_of_subsets)) * 100
                        writer.writerow(write_row)


'''Use protein coding dict as 3rd dict only'''
# #run this for figure1B, silence if running for figure 7
# expression_score(aut_dict, ups_dict, prot_dict, system1='AUT', system2='UPS', system3='Protein-coding')
# #run this for figure7B, silence if running for figure 1
expression_score(degradation_dict, chap_dict, syn_dict, system1='Degradation', system2='Folding', system3='Synthesis')






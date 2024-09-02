import csv
import scipy.stats as stats
import statistics
import pandas as pd


'''This code produces data for 2 figures (Fig1A, Fig7A), in order to get the results for specific figure silence the section relating
to the other figure. To run the code copy source files into the same folder with the code. 
 Files needed to run the code:
 median_expression_per_tissue_GTEx_2023.csv
 all GTex names conslist ENSG UniprotID AUT no dups w our genes.csv
 all GTex names conslist ENSG UniprotID UPS no dups w our genes.csv
 consortium file with ENSG based on UniprotID SYN no dup.csv
 ENSG new gene list filtered to relevant chaperones.csv
 protein coding apr 2023.txt
 '''

'''--------------------------------------Conditions------------------------------------------------------------------'''
threshold = 5  # genes are considered expressed if they are expressed above this threshold
numbers_of_tissues = 34
exclusion_threshold = 5  # genes must be expressed above this threshold at least in one tissue



'''------------------------------------Common source files and output-----------------------------------------------------------'''
median_expression = pd.read_csv(('median_expression_per_tissue_GTEx_2023.csv'), delimiter=',')
combined_output_file = 'A final_expression_in_tissues_per_gene_'+str(threshold)+'TPM_boxplot.csv'


'''------------------------------------Figure-specific source files and output---------------------------------------'''
'''Figure 1A, SILENCE THIS IF RUNNING FOR FIG7A'''
# autophagy_file = pd.read_csv('all GTex names conslist ENSG UniprotID AUT no dups w our genes.csv')
# ups_file = pd.read_csv('all GTex names conslist ENSG UniprotID UPS no dups w our genes.csv')
# protein_coding = 'protein coding apr 2023.txt'
# autophagy_output_file = 'Autophagy expressing tissues '+str(threshold)+'TPM.csv'
# protein_coding_output_file = 'Protein coding expressing tissues ' + str(threshold) + 'TPM.csv'
# ups_output_file = 'UPS expressing tissues '+str(threshold)+'TPM.csv'


'''Figure 7A, SILENCE THIS IF RUNNING FOR FIG1A'''
autophagy_file = pd.read_csv('all GTex names conslist ENSG UniprotID AUT no dups w our genes.csv')
ups_file = pd.read_csv('all GTex names conslist ENSG UniprotID UPS no dups w our genes.csv')
syn_file = pd.read_csv('consortium file with ENSG based on UniprotID SYN no dup.csv')
folding_file = pd.read_csv('ENSG new gene list filtered to relevant chaperones.csv')
folding_output_file = 'Folding expressing tissues '+str(threshold)+'TPM.csv'
degradation_output_file = 'Degradation expressing tissues '+str(threshold)+'TPM.csv'
synthesis_output_file = 'Synthesis expressing tissues '+str(threshold)+'TPM.csv'


def csv_tissue_expression(file, output_file, gene_column, system_name ):
    '''function counts number of expressing tissues for the gene set from the file and writes it to file'''

    gene_list = []
    if system_name == 'Degradation':
        for rowNum, row in autophagy_file.iterrows():
            if row['ENSG'] not in gene_list:
                gene_list.append(row['ENSG'])

        for rowNum, row in ups_file.iterrows():
            if row['ENSG'] not in gene_list:
                gene_list.append(row['ENSG'])
    else:
        for rowNum, row in file.iterrows():
            gene_list.append(row[gene_column])


    genes_above_exclusion_threshold = []
    gene_names_above_exclusion_threshold = []

    for rowNum, row in median_expression.iterrows():
        row = row.to_dict()
        if row['ENSG'].strip()[0:15] in gene_list:
            temp_list = []
            for item in row:
                if item not in ['ENSG', 'Name'] and row[item] == row[item]:
                    temp_list.append(float(row[item]))

            #Gene must be expressed above the exclusion threshold in at least one tissue
            wrong_TPM = all(count < exclusion_threshold for count in temp_list)
            if not wrong_TPM:
                # print(wrong_TPM)
                genes_above_exclusion_threshold.append(row['ENSG'].strip()[0:15])
                gene_names_above_exclusion_threshold.append(row['Name'])
    print(len(genes_above_exclusion_threshold), system_name + ' above 5TPM')
    print(len(gene_names_above_exclusion_threshold), system_name + ' names')

    fd = [system_name + ' ENSG', system_name + ' Name', 'number of expressed tissues']
    writer = csv.DictWriter(open(output_file, 'w', newline=""), fieldnames=fd)
    writer.writeheader()
    write_row = {}

    for rowNum, row in median_expression.iterrows():
        row = row.to_dict()
        if row['ENSG'].strip()[0:15] in genes_above_exclusion_threshold:
            gene_expression = 0
            for item in row:
                if item not in ['ENSG', 'Name', '']:
                    if float(row[item]) >= threshold:
                        gene_expression += 1

            write_row[system_name + ' ENSG'] = row['ENSG'].strip()[0:15]
            write_row[system_name + ' Name'] = row['Name']
            write_row['number of expressed tissues'] = gene_expression
            writer.writerow(write_row)


def protein_tissue_expression():
    '''function counts number of expressing tissues for the protein coding genes and writes it to file'''
    proteins_list = []
    with open(protein_coding, 'r', newline='') as protein_coding_object:
        protein_coding_file = protein_coding_object.readlines()

        for line in protein_coding_file:
            line = line.strip().split('\t')
            proteins_list.append(line[0])

    autophages_list = []
    for rowNum, row in autophagy_file.iterrows():
        autophages_list.append(row['ENSG'])

    ups_list = []
    for rowNum, row in ups_file.iterrows():
        autophages_list.append(row['ENSG'])


    protein_genes_above_exclusion_threshold = []
    protein_gene_names_above_exclusion_threshold = []

    for rowNum, row in median_expression.iterrows():
        row = row.to_dict()
        if row['ENSG'].strip()[0:15] in proteins_list:
            temp_list = []
            for item in row:
                if item not in ['ENSG', 'Name'] and row[item] == row[item]:
                    temp_list.append(float(row[item]))

            wrong_TPM = all(count < exclusion_threshold for count in temp_list)

            if not wrong_TPM:
                protein_genes_above_exclusion_threshold.append(row['ENSG'].strip()[0:15])
                protein_gene_names_above_exclusion_threshold.append(row['Name'])
    print(len(protein_genes_above_exclusion_threshold), 'protein-coding above 5TPM')
    print(len(protein_gene_names_above_exclusion_threshold), 'protein-coding names')


    fd = ['Protein ENSG', 'Protein Name', 'number of expressed tissues']
    writer = csv.DictWriter(open(protein_coding_output_file, 'w', newline=""),
                            fieldnames=fd)
    writer.writeheader()
    write_row = {}


    for rowNum, row in median_expression.iterrows():
        row = row.to_dict()
        if row['ENSG'].strip()[:15] in protein_genes_above_exclusion_threshold and\
                row['Name'].strip() not in autophages_list and row['Name'].strip() not in ups_list:
            gene_expression = 0
            for item in row:
                if item not in ['ENSG', 'Name']:
                    if float(row[item]) >= threshold:
                        gene_expression += 1

            write_row['Protein ENSG'] = row['ENSG'].strip()[0:15]
            write_row['Protein Name'] = row['Name']
            write_row['number of expressed tissues'] = gene_expression
            writer.writerow(write_row)


def combined_percent_expressing_tissues(file1, file2, file3, system1, system2, system3):
    '''function calculates the percent of expressing tissues per each gene for all the systems and writes it to file'''
    file1_tissue_expression = pd.read_csv(file1)
    file2_tissue_expression = pd.read_csv(file2)
    file3_tissue_expression = pd.read_csv(file3)

    print(system1)
    print(file1_tissue_expression.columns.tolist())
    print(file3_tissue_expression.columns.tolist())


    list_of_files = [file1_tissue_expression, file2_tissue_expression, file3_tissue_expression]

    fd = ['ENSG', 'Name', '% Tissues', 'System']
    writer = csv.DictWriter(open(combined_output_file, 'w', newline=""), fieldnames=fd)
    writer.writeheader()
    write_row = {}

    file1_array = []
    file2_array = []
    file3_array = []

    annotation_dict = {'Autophagy': 'AUT', 'UPS': 'UPS', 'Protein': 'Protein-coding', 'Synthesis': 'Synthesis', \
                       'Folding': 'Folding', 'Degradation':'Degradation'}

    for file in list_of_files:
        for lineNum, line in file.iterrows():
            line = line.to_dict()
            print (system1 + ' ENSG')
            if system1 + ' ENSG' in line.keys():
                print(line['number of expressed tissues'])
                write_row['ENSG'] = line[system1 + ' ENSG']
                write_row['Name'] = line[system1 + ' Name']
                write_row['% Tissues'] = (int(line['number of expressed tissues'])/34)*100
                write_row['System'] = annotation_dict[system1]
                writer.writerow(write_row)

                file1_array.append(int(line['number of expressed tissues']))

            if system2 + ' ENSG' in line.keys():
                print(line['number of expressed tissues'])
                write_row['ENSG'] = line[system2 + ' ENSG']
                write_row['Name'] = line[system2 + ' Name']
                write_row['% Tissues'] = (int(line['number of expressed tissues'])/34)*100
                write_row['System'] = annotation_dict[system2]
                writer.writerow(write_row)

                file2_array.append(int(line['number of expressed tissues']))


            if system3 + ' ENSG' in line.keys():
                print(line['number of expressed tissues'])
                write_row['ENSG'] = line[system3 + ' ENSG']
                write_row['Name'] = line[system3 + ' Name']
                write_row['% Tissues'] = (int(line['number of expressed tissues'])/34)*100
                write_row['System'] = annotation_dict[system3]
                writer.writerow(write_row)

                file3_array.append(int(line['number of expressed tissues']))



    print(system1 + ' vs ' + system3, stats.mannwhitneyu(file1_array, file3_array))
    print(system1 + ' vs ' + system2, stats.mannwhitneyu(file1_array, file2_array))
    print(system2 + ' vs ' + system3, stats.mannwhitneyu(file2_array, file3_array))

'''-------------------------------------Run function for the correct figure------------------------------------------'''
# Parameters system1-3 are comments specifying the system of the corresponding file, if file1 is autophagy system1= "Autophagy"
# Use "UPS" for UPS file, "Protein" for protein-coding file, "Synthesis" for synthesis file, "Folding" for folding file
# Use protein-coding file as the third file only.


'''Run this for figure1 A, SILENCE IF YOU RUN FOR FIG7A'''
# csv_tissue_expression(autophagy_file, autophagy_output_file, gene_column = 'ENSG', system_name='Autophagy')
# csv_tissue_expression(ups_file, ups_output_file, gene_column='ENSG', system_name='UPS')
# protein_tissue_expression()
# combined_percent_expressing_tissues(autophagy_output_file, ups_output_file, protein_coding_output_file,
#                      system1= 'Autophagy', system2= 'UPS', system3= 'Protein')
'''Run this for figure7 A, SILENCE IF YOU RUN FOR FIG1A'''
csv_tissue_expression(autophagy_file, degradation_output_file, gene_column='ENSG', system_name='Degradation')
csv_tissue_expression(folding_file, folding_output_file, gene_column='ENSG', system_name='Folding')
csv_tissue_expression(syn_file, synthesis_output_file, gene_column='ENSG', system_name='Synthesis')
combined_percent_expressing_tissues(degradation_output_file, folding_output_file, synthesis_output_file,
                          system1='Degradation', system2='Folding', system3='Synthesis')
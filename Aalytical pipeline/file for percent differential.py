import csv
import scipy.stats as stats
import statistics
import pandas as pd

'''------------------------Source files required to run this code for figure 1G and 7G-------------------------------'''
# '''--------------------------------------Use this for figure 1G,  silence for 7G-----------------------------------'''
file = pd.read_csv('consortium core var ALL with diseases combined.csv')

# '''--------------------------------------Use this for figure 7G,  silence for 1G-----------------------------------'''
file = pd.read_csv('consortium core var ALL with diseases combined.csv')


# aut_select_gene_file = pd.read_csv('consortium core var AUT with diseases with dups SELECTIVE.csv')
# ups_proteasome_gene_file = pd.read_csv('core_variable_differential_V8_UPS PROTEASOME_5TPM_core_consortium_trimmed.csv')

'''------------------------------------Obtain gene - system dictionary-----------------------------------------------'''
def gene_category_dict(file):
    dictionary = {}

    for rowNum, row in file.iterrows():
        row = row.to_dict()
        if row['ENSG'] in dictionary:
            dictionary[row['ENSG']].append(row['system'])
        else:
            dictionary[row['ENSG']] = [row['system']]

    return dictionary

'''--------------------------------------Count the disease associated genes------------------------------------------'''
def count_disease_causing_per_system(file, output_file):

    diff_file = pd.read_csv(file)
    dictionary = {'yes': {}, 'no': {}}

    for lineNum, line in diff_file.iterrows():
        line = line.to_dict()
        for state in dictionary:
            if line['system'] not in dictionary[state]:
                dictionary[state][line['system']] = []


    for lineNum, line in diff_file.iterrows():
        line = line.to_dict()
        dictionary[line['disease causing']][line['system']].append(line['ENSG'])

    for item in dictionary:
        for i in dictionary[item]:
            print(len(dictionary[item][i]))

    fd = ['system', 'disease causing', 'Number of genes']


    writer = csv.DictWriter(open(output_file, 'w', newline=''), fieldnames=fd)
    writer.writeheader()
    write_row = {}

    for state in dictionary:
        for category in dictionary[state]:
            write_row['system'] = category
            write_row['disease causing'] = state
            write_row['Number of genes'] = len(dictionary[state][category])
            writer.writerow(write_row)




'''--------------------------------------Use this for figure 1G,  silence for 7G-------------------------------------'''
count_disease_causing_per_system('consortium core var ALL with diseases combined.csv', 'All system DC percent graph.csv')

'''--------------------------------------Use this for figure 7G,  silence for 1G-------------------------------------'''
count_disease_causing_per_system('consortium core var deg folding syn with diseases combined.csv', 'All system DC percent graph.csv')




import csv
import pickle
import scipy.stats as stats
import statistics
import pandas as pd
'''-----------------------------------Files required to run this code------------------------------------------------'''
ts_score_file = pd.read_csv('TS scores.csv') # to obtain this file please refer to the suplementary tables section of the
# "A Quantitative Proteome Map of the Human Body" PMID: 32916130 - article, tableS2G.

'''-------------------This script calculates percent TS/TE  proteins in each system in all tissues-------------------'''



tissue_list = ts_score_file.columns.tolist()
tissue_list.pop(0)
tissue_list.pop(0)
tissue_list.pop(0)
tissue_list.pop(0)

fd = ['Tissue', 'TS score percent within system', 'Type', 'System']
#writer = csv.DictWriter(open('TS score per tissue in AUT UPS and Chap all tissues in all 3.csv', 'w', newline=''), fieldnames=fd)
writer = csv.DictWriter(open('TS score per tissue in AUT UPS Chap and SYN all tissues.csv', 'w', newline=''), fieldnames=fd)
writer.writeheader()
write_row = {}
system_list = ['3AUT', '4UPS', '2Chap', 'S1YN']

for tissue in tissue_list:
    #file = pd.read_csv('AUT UPS and Chap TS score percent in ' + tissue + ' in all 3.csv')
    file = pd.read_csv('AUT UPS Chap and SYN TS score percent in ' + tissue + '.csv')
    system_ts_te_dict = {'3AUT': 0, '4UPS': 0, '2Chap': 0, '1SYN':0}
    for rowNum, row in file.iterrows():
        row = row.to_dict()
        if row['Type'] == 'Other' or row['Type'] == 'Housekeeping':
            write_row['Tissue'] = tissue
            write_row['System'] = row['System']
            write_row['Type'] = row['Type']
            write_row['TS score percent within system'] = row['TS score percent']
            writer.writerow(write_row)

        if row['Type'] == 'Tissue-enriched' or row['Type'] == 'Tissue-specific':
            system_ts_te_dict[row['System']] += row['TS score percent']

    for system in system_ts_te_dict:
        write_row['Tissue'] = tissue
        write_row['System'] = system
        write_row['Type'] = 'TS/TE'
        write_row['TS score percent within system'] = system_ts_te_dict[system]
        writer.writerow(write_row)






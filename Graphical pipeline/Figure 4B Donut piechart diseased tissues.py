import csv
import pandas as pd
import statistics
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns

'''----------------------------Use these source files to obtain graphs found in figure 4B-----------------------------'''
disease_tissue_file = pd.read_csv('consortium disease to tissue file UPS global dev to multi.csv') # for UPS
#disease_tissue_file = pd.read_csv('consortium disease to tissue file AUT global dev to multi.csv') # for AUTOPHAGY


tissue_list = []

tissue_disease_dict = {'Multiple': 0, 'Global developmental': 0}


for lineNum, line in disease_tissue_file.iterrows():
    line = line.to_dict()
    if line['Tissue'] == line['Tissue']:
        temp_tissues = line['Tissue'].split(',')
        print(temp_tissues)
        '''add single tissues to the dict and count'''
        if 'Multisystem' not in temp_tissues and 'Global developmental delay' not in temp_tissues:
            if len(temp_tissues) == 1 and temp_tissues[0].strip() not in tissue_disease_dict:
                tissue_disease_dict[temp_tissues[0].strip()] = 0
            if len(temp_tissues) == 1 and temp_tissues[0].strip() in tissue_disease_dict:
                tissue_disease_dict[temp_tissues[0].strip()] += 1

            '''count up to 3 tissues as single tissues, means some genes will be counted in 2 or 3 different tissues'''
            if 1 < len(temp_tissues) <= 3:
                for tissue in temp_tissues:
                    if tissue.strip() in tissue_disease_dict:
                        tissue_disease_dict[tissue.strip()] += 1
                    if tissue.strip() not in tissue_disease_dict:
                        tissue_disease_dict[tissue.strip()] = 1


        if 'Multisystem' in temp_tissues or 'Global developmental' in temp_tissues or len(temp_tissues) > 1:
            tissue_disease_dict['Multiple'] += 1

        if 'Global developmental' in temp_tissues:
            tissue_disease_dict['Global developmental'] += 1

        for tissue in temp_tissues:
            if tissue.strip() not in tissue_list:
                tissue_list.append(tissue.strip())




def donut_pie_chart_onelayered(dictionary):

    #sort the keys so that the two files will have the same color coded tissues
    keys = list(dictionary.keys())
    keys.sort()
    sorted_dict = {i: dictionary[i] for i in keys}

    group_lables = []
    for key in sorted_dict:
        group_lables.append(key)

    group_sizes = []
    for key in sorted_dict:
        group_sizes.append(sorted_dict[key])

    brbg_full_color_palette = ["#543005", "#8C510A", "#BF812D", "#DFC27D" ,"#F6E8C3", "#F5F5F5" ,"#C7EAE5","#80CDC1",
                               "#35978F", "#01665E", "#003C30"]
    brbg_short_color_palette = ["#DFC27D", "#F6E8C3", "#C7EAE5", "#80CDC1",
                               "#35978F", "#01665E", "#799a97","#8C510A",  "#BF812D",]
    brbg_short_color_palette_aut = ["#DFC27D", "#F6E8C3",  "#80CDC1",
                                "#35978F", "#01665E", "#799a97", "#8C510A", "#BF812D", ] ##003C30
    plt.pie(group_sizes, labels=group_lables, colors=brbg_short_color_palette,  startangle=90, frame=False) #colors=colors_groups,
    # plt.pie(sizes_dis, colors=colors_dis, radius=0.75, startangle=90)  # labels=lables_dis,
    centre_circle = plt.Circle((0, 0), 0.5, color='black', fc='white', linewidth=0, )
    fig = plt.gcf()
    fig.gca().add_artist(centre_circle)

    plt.axis('equal')
    plt.tight_layout()
    plt.show()

tissue_disease_dict_concise = {'Less than 5 genes':0}
for tissue in tissue_disease_dict:
    if tissue_disease_dict[tissue] >= 5:
        tissue_disease_dict_concise[tissue] = tissue_disease_dict[tissue]
    else:
        tissue_disease_dict_concise['Less than 5 genes'] += 1

#donut_pie_chart_onelayered(tissue_disease_dict)
donut_pie_chart_onelayered(tissue_disease_dict_concise)

#donut_pie_chart_onelayered(brain_muscle_rest_dict)
print(len(tissue_disease_dict))
print(tissue_disease_dict)
print(sorted(tissue_disease_dict_concise))

print(sns.color_palette('mako'))




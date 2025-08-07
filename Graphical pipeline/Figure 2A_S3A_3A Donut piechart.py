import csv
import pandas as pd
import statistics
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

#Note there are comment section in this code that specify which section refferes to which figure
'''-----------------------------------------Source files needed to run this code-------------------------------------'''
# Run this for figure 2A
#file = pd.read_csv('consortium core var UPS with diseases with dups Class Ub and UBL unified.csv')

#Run this for figure 3A
file = pd.read_csv('consortium core var AUT with diseases with dups ESCRT no lys_ref.csv')

#Run this for figure S2A
file = pd.read_csv('core_variable_differential_V8_UPS PROTEASOME_5TPM_core_consortium_trimmed.csv')
'''------------------------------------------------------------------------------------------------------------------'''

labels_colunm = 'Class'
labels_and_sizes_dict = {}
for rowNum, row in file.iterrows():
    row = row.to_dict()
    if row[labels_colunm] in labels_and_sizes_dict:
        labels_and_sizes_dict[row[labels_colunm]] += 1
    if row[labels_colunm] not in labels_and_sizes_dict:
        labels_and_sizes_dict[row[labels_colunm]] = 1

labels_groups = []
label_group_sizes = []
for label in labels_and_sizes_dict:
    print(label, labels_and_sizes_dict[label])
    labels_groups.append(label)
    label_group_sizes.append(labels_and_sizes_dict[label])


#colors_groups = ['#706778','#798996','#ffe699','#ecb87f', '#8cc7cb','#a6b07c'] # Figure 2A

colors_groups = ['#706778','#798996','#8cc7cb','#b4857d','#a9645f','#ecb87f','#ffe699',
                              '#a6b07c','#598264'] # Figure 3A

print(labels_groups)



# lables_groups = ['E2', 'E3', 'Other', "Proteasome"]
# sizes_groups = [36, 378, 134, 54]
# colors_groups = ['#798996','#b4857d','#ecb87f',
#                              '#a6b07c']


# plot
plt.pie(label_group_sizes, labels=labels_groups, colors=colors_groups, startangle=90, frame=False)
# plt.pie(sizes_dis, colors=colors_dis, radius=0.75, startangle=90)  # labels=lables_dis,
centre_circle = plt.Circle((0, 0), 0.5, color='black', fc='white', linewidth=0 )
fig = plt.gcf()
fig.gca().add_artist(centre_circle)

plt.axis('equal')
plt.tight_layout()
plt.show()




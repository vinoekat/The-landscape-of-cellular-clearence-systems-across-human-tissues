import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

### upload your data.
data = pd.read_csv('differential_AUT_heatmap_consortium.csv', sep=',') #use for autophagy - Figure 4B
# for UPS use 'differential_UPS_heatmap_consortium.csv' - Figure 4A


### change the index to whatever your Y in the plot is.
#data = data.rename(index=data['Name'])
genes = data.pop('ENSG')
genes = data.pop('Name')
#genes = data.pop('Pair name')

# Remove not included columns in the heatmap
#data = data.drop(columns=['DepMap_ID', 'CCLE_ID', 'ploidy'])

## font
#sns.set(font_scale=10)

## heatmap
## If you want to remove the labels of y-axis or x-axis, channge yticklabels and xticklabels to 0.
## The default distance function of the heatmap is euclidean. you can change it by adding the parameter metric = ... .
## The distance function can be ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘cityblock’, ‘correlation’,
# ‘cosine’, ‘dice’, ‘euclidean’, ‘hamming’, ‘jaccard’, ‘jensenshannon’, ‘kulsinski’, ‘mahalanobis’,
# ‘matching’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’.

sns.set(font_scale=0.6)
g = sns.clustermap(data, metric='cosine', xticklabels=1, yticklabels= 0, dendrogram_ratio=0.10,
                   cmap="bwr", figsize = (2,4),
                   cbar_pos=(0, 0.85, 0.01, 0.08), center=0)
#if you want to add numbers within the clusters use -> ,annot=True,annot_kws={"size": 1} in the clustermap specifications


plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), fontsize=10, rotation=45, ha="right",
         rotation_mode="anchor")
plt.savefig('heatmap all AUT or UPS.png')
plt.show()

#bwr
#coolwarm
# Expression of ALP and UPS across human tissues
This repository contains data, code, and analysis as described in the paper "Human clearance systems have a layered architecture across tissues and cell
types that is shaped by tissue-specific proteome needs" by Ekaterina Vinogradov-Talyah, Bar Edri, Lior Ravkaie, Or Lazarescu, Fadi
Gharra, Juman Jubran, Anat Ben-Zvi, Esti Yeger-Lotem
# Article results and pipeline
## Pipeline code 
The main code was written in python. Correlations, differential analysis and graphical representation was performed using R. The code used to generate article figures can be found in the "Analytical pipeline" folder. To utilize scripts download the desired script with the relevant source files. The list of the required files is detailed at the top of the script.
## Source files and code
The source files needed for all scripts and source file code can be found on Zenode - https://zenodo.org/records/13626365. 
## Graphical pipeline
Scripts used to create each type of plot can be found in the "Graphical pipeline" folder.
## Output files
The output files can be found in the "Article results" folder. Each file *.xlsx file summarizes data that generates plots for each figure. Sheet name and letters correspont to the figure letters.
## Script-to-figure details
|Figure|Panel|Scripts - analytical pipeline|Scripts - graphical pipeline|
|:-----|:-----|:------|:------|
|Figure1/8|A|Figure 1A_8A expression per gene.py|Boxplot.R|
|Figure1/8|B|Figure 1B_8B expressing_cell types.py|Boxplot.R|
|Figure1|C|Figure 1C TS score per gene and in specific tissue.py|Boxplot.R|
|Figure1|D|Figure 1D tissue_and_timepoint_gene_expression.py|Boxplot.R|
|Figure1|F|Figure 1F CRISPR_score_min.py|Figure 1F Crisper plot.py|
|Figure1G/8I|G/I|Figure 1G_8I disease percent.py|Percent stacked bar plot.R|
| | | | |
|Figure2/3/S3|A||Figure 2A_S3A_3A Donut piechart.py|
|Figure2/3/S3|B/E|Figure 2B_E 3B_E S3_E expressing tissues and CRISPR.py|Boxplot.R|
|Figure2/3/S3|C|Figure 2C_3C_S3C_expressing_cell_types.py|Boxplot.R|
|Figure2/3/S3|D|Source code - TS score AUT and UPS by category per tissue w dups.py (run "write_general_file" function)|Figure 2D_3D_S3D Bubble plot discrete.R|
|Figure2/3/S3|F||Chord diagram.R|
|||||
|Figure4|A|Differential analysis as described in article's methods|Heatmap script differential.py|
|Figure4|C|Manual curation|Figure 4C Donut piechart diseased tissues.py|
|||||
|Figure6|A||Chord diagram.R|
|Figure6|B|Figure 6B expressing cell types CORE VARIABLE.py|Boxplot.R|
|Figure6|C|Figure 6C tissue_and_timepoint_gene_expression _CORE_VARIABLE.py|Boxplot.R|
|Figure6|D|Figure 6D median_expression_CORE_VARIABLE.py|Boxplot.R|
|Figure6|H|Figure 6H CRISPR_score_min_CORE_VARIABLE.py|Figure 6H Crisper plot CORE VARIABLE.py|
|||||
|Figure7|A|Figure 7A FC bigger smaller than 1 differential.py|Figure 7 bubble plot discrete.R, Correlation plot ggplot.R|
|Figure7|B|Folder - Figure7 B TS_TE score|Figure 7 bubble plot discrete.R, Correlation plot ggplot.R|
|Figure7|C|Figure 7 C FC bigger smaller than 2 preferential.py|Figure 7 bubble plot discrete.R, Correlation plot ggplot.R|
|||||
|Figure8|A|mentioned above||
|Figure8|B|mentioned above||
|Figure8|C|Figure 8C TS score per gene and in specific tissue.py|Boxplot.R|
|Figure8|D|Figure 8D tissue_and_timepoint_gene_expression.py|Boxplot.R|
|Figure8|E|Figure 8E median_expression.py|Boxplot.R|
|Figure8|H|Figure 8H CRISPR_score_min.py|Figure 8H Crisper plot.py|
|Figure8|I|mentioned above||
## Contact
Esti Yeger-Lotem - estiyl@bgu.ac.il, Ekaterina Vinogradov - vinoekat@gmail.com 



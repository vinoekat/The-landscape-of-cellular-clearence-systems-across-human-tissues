# Expression of ALP and UPS across human tissues
This repository contains data, code, and analysis as described in the paper "The landscape of cellular clearance systems across human tissues and cell types is shaped by tissue-specific proteome needs" by Ekaterina Vinogradov, Lior Ravkaie, Bar Edri, Juman Jubran, Anat Ben-Zvi and Esti Yeger-Lotem
# Article results
## Pipeline code 
The main code was written in python. Correlations, differential analysis and graphical representation was performed using R. The code used to generate article figures can be found in the "Analytical pipeline" folder. To utilize scripts download the desired script with the relevant source files. The list of the required files is detailed at the top of the script.
## Source files and code
The source files needed for all scripts can be found in the "Source files" folder. Single cell and protein files needed for the analysis can be found in their respective folders in the "Source files" folder. The code used to create the source files can be found in the "Source file code" folder. Conditions needed to run the source file code are detailed at the top of each script.
## Graphical pipeline
Scripts used to create each type of plot can be found in the "Graphical pipeline" folder.
## Output files
The output files can be found in the "Article results" folder. Each file *.xlsx file summarizes data that generates plots for each figure. Sheet name and letters correspont to the figure letters.
## Script-to-figure details
|Figure|Panel|Scripts - analytical pipeline|Scripts - graphical pipeline|
|:-----|:-----|:------|:------|
|Figure1/7|A|Figure 1A_7A expression per gene.py|Boxplot.R|
|Figure1/7|B|Figure 1B_7B expressing_cell types.py|Boxplot.R|
|Figure1|C|Figure 1C TS score per gene and in specific tissue.py|Boxplot.R|
|Figure1|D|Figure 1D tissue_and_timepoint_gene_expression.py|Boxplot.R|
|Figure1|F|Figure 1F CRISPR_score_min.py|Figure 1F Crisper plot.py|
|Figure1/7|G|Figure 1G_7G disease percent.py|Percent stacked bar plot.R|
| | | | |
|Figure2/3/S2|A||Figure 2A_S2A_3A Donut piechart.py|
|Figure2/3/S2|B/E|Figure 2B_E 3B_E S2B_E expressing tissues and CRISPR.py|Boxplot.R|
|Figure2/3/S2|C|Figure 2C_3C_S2C_expressing_cell_types.py|Boxplot.R|
|Figure2/3/S2|D|Source code - TS score AUT and UPS by category per tissue w dups.py (run "write_general_file" function)|Figure 2D_3D_S2D Bubble plot discrete.R|
|Figure2/3/S2|F||Chord diagram.R|
|||||
|Figure4|A|Differential analysis as described in article's methods|Heatmap script differential.py|
|Figure4|B|Manual curation|Figure 4B Donut piechart diseased tissues.py|
|||||
|Figure5|A||Chord diagram.R|
|Figure5|B|Figure 5B expressing cell types CORE VARIABLE.py|Boxplot.R|
|Figure5|C|Figure 5C tissue_and_timepoint_gene_expression _CORE_VARIABLE.py|Boxplot.R|
|Figure5|D|Figure 5D median_expression_CORE_VARIABLE.py|Boxplot.R|
|Figure5|E|Figure 5E CRISPR_score_min_CORE_VARIABLE.py|Figure 5E Crisper plot CORE VARIABLE.py|
|||||
|Figure6|A|Figure 6A FC bigger smaller than 1 differential.py|Figure 6 bubble plot discrete.R, Correlation plot ggplot.R|
|Figure6|B|Folder - Figure6 B TS_TE score|Figure 6 bubble plot discrete.R, Correlation plot ggplot.R|
|Figure6|C|Figure 6 C FC bigger smaller than 2 preferential.py|Figure 6 bubble plot discrete.R, Correlation plot ggplot.R|
|||||
|Figure7|A|mentioned above||
|Figure7|B|mentioned above||
|Figure7|C|Figure 7C TS score per gene and in specific tissue.py|Boxplot.R|
|Figure7|D|Figure 7D tissue_and_timepoint_gene_expression.py|Boxplot.R|
|Figure7|E|Figure 7E median_expression.py|Boxplot.R|
|Figure7|F|Figure 7F CRISPR_score_min.py|Boxplot.R|
|Figure7|G|mentioned above||
## Contact
Esti Yeger-Lotem - estiyl@bgu.ac.il, Ekaterina Vinogradov - vinoekat@gmail.com 



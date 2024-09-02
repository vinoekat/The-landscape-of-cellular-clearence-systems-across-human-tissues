library(data.table)
library(ggplot2)
library(corrplot)
library(colorRampPalette)

# The folowing code builts correlation graphs for figure 6 A/B/C using following source files:
#Figure 6A - FC above 1 per tissue in AUT UPS Chap and SYN consortium WIDE pval.csv
#FC bellow 1 per tissue in AUT UPS Chap and SYN consortium WIDE pval.csv

#Figure 6B - TS score per tissue in AUT UPS Chap and SYN all tissues TS TE WIDE.csv

#Figure 6C - preferential FC above 2 per tissue in AUT UPS Chap and SYN cell line names WIDE.csv
#preferential FC bellow 2 per tissue in AUT UPS Chap and SYN cell line names WIDE.csv


data =fread("file path and name")
names(data)
df <- data.frame(data, row.names=1)
col<-colorRampPalette(c("blue","grey","red"))(20)
              
#corrplot(cor(df, method = "pearson"), type="lower",method="square",
         #col=col, col.lim = c(0.0, 1), is.corr = FALSE)

corrplot.mixed(cor(df, method = 'pearson'), upper = "square", 
               lower = "number", 
               addgrid.col = "black",
               tl.col = "black",col.lim = c(0, 1),
               is.corr = FALSE)

#corrplot(cor(df), type="lower",method="color", col=col, col.lim = c(0.4, 1), is.corr = FALSE)
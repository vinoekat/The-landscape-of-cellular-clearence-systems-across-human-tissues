library(ggplot2)
library(data.table)
library(dplyr)
library(viridis)
library(hrbrthemes)

# The folowing code builts graphs for figure 6 A/B/C using following source files:
#Figure 6A - FC above 1 per tissue in AUT UPS Chap and SYN consortium.csv
#FC bellow 1 per tissue in AUT UPS Chap and SYN consortium.csv

#Figure 6B - TS score per tissue in AUT UPS Chap and SYN all tissues TS TE.csv

#Figure 6C - preferential FC above 2 per tissue in AUT UPS Chap and SYN cell line names.csv
#preferential FC bellow 2 per tissue in AUT UPS Chap and SYN cell line names.csv

data =fread("file path and name")
df <- data.frame(x_axis = data$'Tissue', percent_within_category = data$'Use percent category column name', group = data$'System') 
names(df) <- c("x_axis", "Use percent category column name", "group")
ggplot(df, aes(x=x_axis, y = group, size = percent_within_category))+

  geom_point(aes(color = percent_within_category))+
  #scale_color_gradient(low="black", high='red')+
  scale_fill_manual( "#453854")+
  scale_size(range=c(0, 10), name = 'Type')+
  
  theme_set(theme_bw() + theme(legend.key=element_blank()))+
  theme(panel.grid.major = element_line(color = "#8ccde3",
                                        size = 0.5,
                                        linetype = 2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="right",
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.key.size = unit(.5, "cm"),
        text = element_text(size=15),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle=45, hjust=1 , vjust=1 ,
                                 size = 10, color = 'black'))+
  theme(axis.text.y=element_text(size = 10, color = 'black'))+
xlab("")+
ylab("")

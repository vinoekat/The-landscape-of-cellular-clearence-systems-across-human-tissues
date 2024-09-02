library(ggplot2)
library(data.table)

# File name in fread - figure S2D
data =fread("C:/file path/TS score percent per category UPS PROTEASOME changed names TS TE united1.csv")
# file name for figure 2D - TS score percent per category UPS Ub and UBL SELECT TS TE united1.csv
# file name for figure 3D -TS score percent per category AUT Ub and UBL SELECT TS TE united1.csv
df <- data.frame(x_axis = data$'Class', percent_within_category = data$'TS score percent per class', 
                 group = data$'TS type') 
names(df) <- c("x_axis", "percent_within_category", "group")
ggplot(df, aes(x=x_axis, y = group, size = percent_within_category))+

  geom_point(aes(color = percent_within_category))+
  scale_color_gradient(low="black", high='red')+
  scale_size(range=c(0, 12), name = 'Type')+
  
  theme_set(theme_bw() + theme(legend.key=element_blank()))+
  theme(panel.grid.major = element_line(color = "#8ccde3",
                                        size = 0.5,
                                        linetype = 2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="right",
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.key.size = unit(.5, "cm"),
        text = element_text(size=15),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle=45, hjust=1 , vjust=1 ,
                                 size = 15, color = 'black'))+
  theme(axis.text.y=element_text(size = 15, color = 'black'))+
xlab("")+
ylab("")
  
 
library(data.table)
library(ggplot2)
# insert values based on the file you use, refer to article results folder for files.
# copy needed sheet from excel file into csv format.
data =fread("File path/Long file  format.csv")
names(data)
df <- data.frame(x_axis = data$'x axis values coulmn name', 
                 y_axis = data$'Score/Percent column name', 
                 group = data$'Group/System column name')
names(df) <- c("x_axis", "y_axis","group")
ggplot(df, aes(x=x_axis, y=y_axis, fill=group)) +
  xlab("")+
  ylab('y axis name')+
  ggtitle("")+
  theme( axis.title=element_text(size=10))+
  theme(axis.text.x=element_text(angle=45, hjust=1 , vjust=1 , size = 15, 
                                 color = 'black'))+
  theme(axis.text.y=element_text(size = 15,color = 'black'))+
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),panel.border = element_blank())+
  #geom_jitter(width = 0.3, size =1, colour='black',shape=101)+
  
  geom_boxplot(alpha = 0.7, width=0.6, color='black',
               position = position_dodge(width=0.8), 
               outlier.size = 3, outlier.fill ='white', 
               outlier.colour = 'black', outlier.shape = 21, 
               outlier.stroke = 1  )+
  scale_fill_manual(values=c('#1395ba','#a2b86c','#0d3c55'))+
  theme(legend.position = 'top',legend.text = element_text(size=15))

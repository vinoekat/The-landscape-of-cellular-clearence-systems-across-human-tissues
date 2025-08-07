library(data.table)
library(ggplot2)
# insert values based on the file you use, refer to article results folder for files.
# copy needed sheet from excel file into csv format.
data = fread("File path/Long file  format.csv")
names(data)
df<- data.frame(x_axis = data$'System/Group column name', 
                y_axis = data$'Number of genes containing column', 
                group = data$'Gene number column name')
names(df)<-c("x_axis","y_axis","group")
ggplot(df, aes(x=x_axis, y = y_axis, fill = group))+
  geom_bar(stat="identity",color="white", position="fill")+
  #use 3 colors for plots with 3 percent groups
  scale_fill_manual(values = c('white','red'))+
  
  xlab("categories") + ylab("percent group")+
  scale_y_continuous(expand = c(0,0))+
  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="top",
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.key.size = unit(1, "cm"),
        text = element_text(size=15),
        axis.line = element_line(colour = "black"))+
theme(axis.text.x=element_text(angle=90, hjust=1 , vjust=1 ,
                               size = 10, color = 'black'))



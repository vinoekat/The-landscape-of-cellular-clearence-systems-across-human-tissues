library(data.table)
library(ggplot2)
library(ggalluvial)
library(circlize)

#File containing 2 groups for x axis(either disease causing yes/no or core/variable)
#y axis - classes
# these files can be found in the source files folder:
#consortium core var UPS with diseases PROTEASOME.csv - Figure S2F
#consortium core var AUT with diseases with dups ESCRT no lys_ref.csv - Figure 3F,5A
#consortium core var UPS with diseases with dups Class Ub and UBL unified.csv - Figure 2F

#x axis either 'core/var'(fig5) or 'disease causing(fig2,3,S2)
data =fread("file path.csv")
names(data)
df <- data.frame(x_axis = data$'disease causing', y_axis = data$'Class') 
names(df) <- c("x_axis", "y_axis")

# for AUT - figure 3, 5
#chordDiagram(df, annotationTrack = "grid", preAllocateTracks = 1, grid.col = c('#fdee72','#8a41bb', "#317576", "#598264", "#a6b07c","#ffe699","#ecb87f", "#a9645f","#b4857d", "#8cc7cb"), transparency = 0.5)

#for UPS - figure 2,5
chordDiagram(df, annotationTrack = "grid", preAllocateTracks = 1, grid.col = c('grey','red', "#317576", "#598264", "#a6b07c","#ffe699","#ecb87f", "#a9645f", 
                                                                               "#b4857d", "#8cc7cb"), transparency = 0.5)
circos.trackPlotRegion(track.index = 1, 
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         sector.name = get.cell.meta.data("sector.index")
                         circos.text(mean(xlim), ylim[1] + .1, 
                                     sector.name, facing = "clockwise",
                                     niceFacing = TRUE, adj = c(-0.1, 0.5))
                         #circos.axis( h= "top", labels.cex = 0.5, major.tick = 0.2,
                         #sector.index = sector.name, track.index = 2) ticks and numbers
                       }, bg.border = NA)

circos.clear()
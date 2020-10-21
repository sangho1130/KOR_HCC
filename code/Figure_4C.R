library(RGtk2)
library(rsgcc)
library(RColorBrewer)

inputCibersort <- '../data/Figure 4C input.txt'


data <- read.delim(inputCibersort, sep = '\t', header = TRUE, check.names = FALSE, row.names = 1)
data <- data[,c(2,5,8,11)]
head(data)
my_palette <- colorRampPalette(c('grey90', "steelblue", 'grey30', "red2", 'grey90'))(n = 12)

heatmap.2(as.matrix(data), trace = 'none', na.color = 'grey90',
          Rowv = FALSE, Colv = FALSE, dendrogram = 'none', scale = 'none',
          breaks = c(-0.12, -0.1, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12), 
          col = my_palette, density.info = 'none', main = title, margins = c(6,12), 
          labCol = c('Partial', 'Total', 'Partial', 'Total'))
dev.off()

library(RGtk2)
library(rsgcc)
library(RColorBrewer)
library(dendsort)

input <- '../data/Figure 6D input.txt'
title <- 'Cytokines/chemokines and receptors'

data <- read.delim(input, sep = '\t', header = TRUE, check.names = FALSE, row.names = 1)
data$HR <- NULL
data$logRank <- NULL
data$HR_pretrx <- NULL
data$logRank_pretrx <- NULL

my_palette <- colorRampPalette(c("#ffffcc", "#ff6600"))(n = 20)
heatmap.2(as.matrix(t(data)), trace = 'none', na.color = 'grey95',
          Rowv = F, Colv = F, dendrogram = 'none', scale = 'none',
          breaks = c(-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
          col = my_palette, density.info = 'none', main = title)
dev.off()


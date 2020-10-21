library(RGtk2)
library(rsgcc)
library(RColorBrewer)

inputCibersort <- '../data/Figure 7B input exhausted average.txt'; title <- 'Exhausted T cell'
data <- read.delim(inputCibersort, sep = '\t', header = TRUE, check.names = FALSE, row.names = 1)
my_palette <- colorRampPalette(c('grey90', "red2"))(n = 5)

#outputPdf <- paste(c(strsplit(inputCibersort, split = '.txt')[1], '.pdf'), collapse = '')
pdf('../results/Figure 7B input exhausted average.pdf')
heatmap.2(as.matrix(data), trace = 'none', na.color = 'grey95', 
          Rowv = FALSE, Colv = FALSE,
          dendrogram = 'none', scale = 'row',
          col = my_palette, density.info = 'none', main = title, margins = c(8,8))
dev.off()


inputCibersort <- '../data/Figure 7B input mait average.txt'; title <- 'MAIT T cell'
data <- read.delim(inputCibersort, sep = '\t', header = TRUE, check.names = FALSE, row.names = 1)
my_palette <- colorRampPalette(c('grey90', "red2"))(n = 5)

#outputPdf <- paste(c(strsplit(inputCibersort, split = '.txt')[1], '.pdf'), collapse = '')
pdf('../data/Figure 7B input mait average.pdf')
heatmap.2(as.matrix(data), trace = 'none', na.color = 'grey95', 
          Rowv = FALSE, Colv = FALSE,
          dendrogram = 'none', scale = 'row',
          col = my_palette, density.info = 'none', main = title, margins = c(8,8))
dev.off()


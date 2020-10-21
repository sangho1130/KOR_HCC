#library(gplots)
library(devtools)
library(dendsort)
library(plyr)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

inputCibersort <- '../data/Figure 4A input partial hepatectomy.txt'; title <- 'Partial hepatectomy (n=19)'
data <- read.delim(inputCibersort, sep = '\t', header = TRUE, check.names = FALSE, row.names = 1)
colColors <- data[,18:27]
data <- data[,-c(18:27)]

head(colColors)
colColors$Treg <- mapvalues(colColors$Treg, from = c('Absent', 'Present'), to = c('grey', '#f9a40d'))
colColors$Angioinvasion <- mapvalues(colColors$Angioinvasion, from = c('Absent', 'Present'), to = c('grey', '#199332'))
colColors$Recurrence <- mapvalues(colColors$Recurrence, from = c('No', 'Yes'), to = c('grey', 'red2'))
colColors$Pretreatment <- mapvalues(colColors$Pretreatment, from = c('No', 'NotAvail', 'Yes'), to = c('grey', 'black', '#e81686'))
colColors$Tstage <- mapvalues(colColors$Tstage, from = c('T1', 'T2', 'T3/4'), to = c('#f0b688', '#ea6c63', '#f2160c'))
colColors$Grade <- mapvalues(colColors$Grade, from = c('G1', 'G2', 'G3/4'), to = c('#9eb1c2', '#6c87bc', '#2e5b9f'))
colColors$Deceased <- mapvalues(colColors$Deceased, from = c('Alive', 'Dead'), to = c('grey', 'red2'))
colColors$NMF <- mapvalues(colColors$NMF, from = c('C1', 'C2', 'C3'), to = c('#dda119', '#1fa2c0', '#ed6ead'))
colColors$Hoshida <- mapvalues(colColors$Hoshida, from = c('S1', 'S2', 'S3'), to = c('red2', 'royalblue', 'goldenrod1'))
colColors$Chiang <- mapvalues(colColors$Chiang, from = c('CTNNB1', 'Proliferation', 'Interferon', 'Poly7', 'Unannotated'), to = c('#0883d3', '#007901', '#9e017d', '#fe3e09', '#b1b1b1'))

#outputPdf <- paste(c(strsplit(inputCibersort, split = '.txt')[1], '.pdf'), collapse = '')
pdf('../results/Figure 4A input partial hepatectomy.pdf')
heatmap.3(as.matrix(t(data)), trace = 'none', na.color = 'grey95', 
          Rowv = dendsort(as.dendrogram(hclust(dist(t(data), method = 'euclidean'), method = 'complete'))),
          Colv = dendsort(as.dendrogram(hclust(dist(data, method = 'euclidean'), method = 'complete'))), 
          dendrogram = 'both', scale = 'row', ColSideColors = as.matrix(colColors), ColSideColorsSize = 10, 
          col = colorRampPalette(c('#333333', '#545454', '#00b2ee'))(n=10), density.info = 'none', main = title)
dev.off()



inputCibersort <- '../data/Figure 4A input total hepatectomy.txt'; title <- 'Total hepatectomy (n=35)'
data <- read.delim(inputCibersort, sep = '\t', header = TRUE, check.names = FALSE, row.names = 1)
colColors <- data[,19:28]
data <- data[,-c(19:28)]

head(colColors)
colColors$Treg <- mapvalues(colColors$Treg, from = c('Absent', 'Present'), to = c('grey', '#f9a40d'))
colColors$Angioinvasion <- mapvalues(colColors$Angioinvasion, from = c('Absent', 'Present'), to = c('grey', '#199332'))
colColors$Recurrence <- mapvalues(colColors$Recurrence, from = c('No', 'Yes'), to = c('grey', 'red2'))
colColors$Pretreatment <- mapvalues(colColors$Pretreatment, from = c('No', 'NotAvail', 'Yes'), to = c('grey', 'black', '#e81686'))
colColors$Tstage <- mapvalues(colColors$Tstage, from = c('T1', 'T2', 'T3/4'), to = c('#f0b688', '#ea6c63', '#f2160c'))
colColors$Grade <- mapvalues(colColors$Grade, from = c('G1', 'G2', 'G3/4'), to = c('#9eb1c2', '#6c87bc', '#2e5b9f'))
colColors$Deceased <- mapvalues(colColors$Deceased, from = c('Alive', 'Dead'), to = c('grey', 'red2'))
colColors$NMF <- mapvalues(colColors$NMF, from = c('C1', 'C2', 'C3'), to = c('#dda119', '#1fa2c0', '#ed6ead'))
colColors$Hoshida <- mapvalues(colColors$Hoshida, from = c('S1', 'S2', 'S3'), to = c('red2', 'royalblue', 'goldenrod1'))
colColors$Chiang <- mapvalues(colColors$Chiang, from = c('CTNNB1', 'Proliferation', 'Interferon', 'Poly7', 'Unannotated'), to = c('#0883d3', '#007901', '#9e017d', '#fe3e09', '#b1b1b1'))

#outputPdf <- paste(c(strsplit(inputCibersort, split = '.txt')[1], '.pdf'), collapse = '')
pdf('../results/Figure 4A input total hepatectomy.pdf')
heatmap.3(as.matrix(t(data)), trace = 'none', na.color = 'grey95', 
          Rowv = dendsort(as.dendrogram(hclust(dist(t(data), method = 'euclidean'), method = 'complete'))),
          Colv = dendsort(as.dendrogram(hclust(dist(data, method = 'euclidean'), method = 'complete'))), 
          dendrogram = 'both', scale = 'row', ColSideColors = as.matrix(colColors), ColSideColorsSize = 10, 
          col = colorRampPalette(c('#333333', '#545454', '#00b2ee'))(n=10), density.info = 'none', main = title)
dev.off()

library(Rtsne)
library(ggplot2)
#library(rgl)
library(plyr)

table <- '../data/Figure 2D input.txt'
data <- read.delim(table, header = TRUE, sep = '\t', row.names = 1, check.names = FALSE)

dataGrade <- data.frame(grade = data$Grade)
data$Grade <- NULL

set.seed(0011) #cibersort
tsne <- Rtsne(as.matrix(data), check_duplicates = F, theta = 0.5, pca_scale = T, dims = 3)#4)
scores <- as.data.frame(tsne$Y)
rownames(scores) <- row.names(data)
colnames(scores) <- c('tSNE1', 'tSNE2', 'tSNE3')#, 'tSNE4')

scores <- cbind(scores, grade = dataGrade$grade)
scores$disease <- mapvalues(scores$grade,
                            from = c('Catholic_Normal', 'GTEx',
                                     'Catholic_CS', 'Catholic_FH', 'Catholic_FL', 'Catholic_DL', 'Catholic_DH', 
                                     'Catholic_T1', 'Catholic_T2', 'Catholic_T3T4', 'Catholic_Mixed',
                                     'GSE77509_Nontumor', 'GSE77509_Tumor', 'GSE77509_PVTT', 
                                     'Riken_Nontumor', 'Riken_T1', 'Riken_T2', 'Riken_T3T4',
                                     'TCGA_Nontumor', 'TCGA_T1', 'TCGA_T2', 'TCGA_T3T4'), 
                            to = c('Normal', 'Normal',
                                   'FibCS', 'FibCS', 'FibCS', 'DN', 'DN', 
                                   'Tumor', 'Tumor', 'Tumor', 'Tumor',
                                   'Nontumor', 'Tumor', 'Tumor', 
                                   'Nontumor', 'Tumor', 'Tumor', 'Tumor',
                                   'Nontumor', 'Tumor', 'Tumor', 'Tumor'))
scores$disease <- factor(scores$disease, c('Normal', 'FibCS', 'DN', 'Nontumor', 'Tumor'))

plt <- ggplot(scores, aes(x = tSNE1, y = tSNE2, colour = disease)) + 
  geom_point(size = .5, alpha=0.8) + 
  scale_colour_manual(name="", values = c('Normal' = '#b3b2b2', 'FibCS' = '#528828', 'DN' = '#f2ae0e', 'Nontumor' = '#f5c2cf', 'Tumor' = '#e8495b')) +
  theme_bw(base_size = 7) + 
  theme(axis.text = element_text(colour = 'black'),
        panel.grid = element_blank()) +
  labs(title =  '')
plt
ggsave('../results/Figure 2D 2dim.pdf', units = 'cm', width = 8, height = 6)
scores.w <- data.frame(SampleID = rownames(scores), scores)
write.table(scores.w, '../results/Figure 2D tSNE projection.txt', quote = F, row.names = F, col.names = T, sep = '\t')

### 3d plot uses rgl 
#scores <- read.delim('../data/Figure 2D tSNE projection.txt', header = TRUE, sep = '\t', row.names = 1, check.names = FALSE)

#clustercolors <- mapvalues(scores$disease, from = c(c('Normal', 'FibCS', 'DN', 'Nontumor', 'Tumor')), to = c('#b3b2b2', '#528828', '#f2ae0e', '#f5c2cf', '#e8495b'))

#plot3d(scores[, 1:3], col = clustercolors, type = 'p', box = FALSE, size = 3)
#decorate3d(main = '', box = FALSE, xlab = '', ylab = '', zlab = '')
#legend3d("right", legend = c('Normal', 'FibCS', 'DN', 'Nontumor', 'Tumor'),
#         col = c('#b3b2b2', '#528828', '#f2ae0e', '#f5c2cf', '#e8495b'), pch = 16, cex=1.5, inset=c(0.01))
#view3d(theta = 280, phi = 10, zoom = .9)
#rgl.postscript(filename = "../results/Figure 2D 3dim.pdf", fmt = "pdf", drawText = TRUE)
#rgl.close()

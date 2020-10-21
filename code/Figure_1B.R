library(FactoMineR)
library(ggplot2)
library(factoextra)

table <- '../data/Figure 1B input.txt'
data <- read.delim(table, header = TRUE, sep = '\t', check.names = FALSE, row.names = 1)

dataGrade <- data.frame(row.names = rownames(data), grade = data$Grade)
data$Grade <- NULL
unique(dataGrade$grade)
dataGrade$grade = factor(dataGrade$grade, c("Normal", 'FL', 'FH', 'CS', 'DL', 'DH', 'T1', 'T2', 'T3-4', 'Mixed'))

pca <- PCA(data)
max(pca$var$contrib[,1:2])
contribution <- as.data.frame(pca$var$contrib)
colnames(contribution) <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5')
contribution <- cbind(gene = rownames(contribution), contribution)

write.table(contribution, '../results/Figure 1B PCA contribution.txt', row.names = F, col.names = T, quote = F, sep = '\t')

fviz_pca_var(pca, col.var = "contrib")+
  scale_color_gradient2(low="white", mid="steelblue2",high="red2", midpoint=4) +
  theme_bw()

eigenvalues <- pca$eig
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "Variances", xlab = "Principal Components", ylab = "Percentage of variances", col ="steelblue")
lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], type="b", pch=19, col = "red")

pcaScores <- as.data.frame(pca$ind$coord)
colnames(pcaScores) <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5')
pcaScores$Grade <- dataGrade$grade

plt <- ggplot(pcaScores, aes(x = PC1, y = PC2, colour = Grade)) + 
  geom_point(size = 1.5) + 
  scale_colour_manual(name='',
                      values = c("Normal" = "#bebdbd", 
                                 "FL" = "#bbe165", 'FH' = '#6e8a3c', 'CS' = '#546a2e', 
                                 "DL" = "#f1c055", 'DH' = '#eb8919', 
                                 "T1" = '#f69693', 'T2' = '#f7474e', 'T3-4' = '#aa0c0b', 'Mixed' = '#570a08')) + 
  theme_bw(base_size = 7) + 
  theme(axis.text = element_text(colour = 'black'), 
        axis.ticks = element_line(colour = 'black'),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())
plt
ggsave('../results/Figure 1B.pdf', plt, units = 'cm', width = 8, height = 6)

head(pcaScores)
writeTable <- data.frame(SampleID = rownames(pcaScores), pcaScores)
write.table(writeTable, '../results/Figure 1B PCA scores.txt', col.names = T, row.names = F, quote = F, sep = '\t')

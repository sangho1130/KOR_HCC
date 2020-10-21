library(reshape2)
library(ggplot2)
library(RColorBrewer)

data <- read.table('../data/Figure 2A input relative score.txt', sep = "\t",header = TRUE, check.names = FALSE)
data$Status <- factor(data$Status, c("Normal", "FL", "FH", "CS", "DL", "DH", "Nontumor", "Tumor"))
melt.data <- melt(data, variable.name = 'Cell', value.name = 'Relative')
p <- ggplot(melt.data ,aes(x = Status, y = Relative, fill = Cell)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(22)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank()) + 
  labs(x = '', y = 'Relative fraction')
p
ggsave('../results/Figure 2A input relative score.pdf', p)

###
data <- read.table('../data/Figure 2A input absolute score.txt', sep = "\t",header = TRUE, check.names = FALSE)
head(data)
data$Status <- factor(data$Status, c("Normal", "FL", "FH", "CS", "DL", "DH", "Nontumor", "Tumor"))
melt.data <- melt(data, variable.name = 'Cell', value.name = 'Absolute')
p <- ggplot(melt.data ,aes(x = Status, y = Absolute, fill = Cell)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(22)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank()) + 
  labs(x = '', y = 'Absolute fraction')
p
ggsave('../results/Figure 2A input absolute score.pdf', p)

library(ggplot2)

data <- read.delim('../data/Figure 4B input.txt', header = TRUE, sep = '\t', row.names = 1, check.names = FALSE)
head(data)

clusteringLeft <- subset(data$AbsoluteScore, data$Clustering == 'Left')
clusteringRight <- subset(data$AbsoluteScore, data$Clustering == 'Right')

t.test(clusteringLeft, clusteringRight)

plt <- ggplot(data, aes(x = Clustering, y = AbsoluteScore, col = as.character(Recurr))) +
  geom_jitter(width = 0.2, size = 1) + 
  scale_color_manual(values = c('#3f82b0', '#ea160a')) +
  theme_classic() + 
  theme(axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(color = 'black'), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = 'none') + 
  labs(title = '', y = 'Absolute Score', x = 'Clustering') + ylim(0, 1)
plt
ggsave('../results/Figure 4B.pdf', plt, height = 5, width = 4, units = 'cm')

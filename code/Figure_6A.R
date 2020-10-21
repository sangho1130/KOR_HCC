library(ggplot2)
library(Barnard)

data <- read.delim('../data/Figure 6A input.txt', header = TRUE, sep = '\t', row.names = 1, check.names = FALSE)

pval <- barnard.test(n1 = 10, n2 = 11, n3 = 2, n4 = 11)$p.value[1]
#     pretx Yes   No
# Treg       10     2
# No treg    11    11 

data$Pretrx <- factor(data$Pretrx, c('NotAvail', 'NO', 'YES'))
plt <- ggplot(data, aes(x = Pretrx, y = Tregs, col = Pretrx)) + 
  geom_jitter(width = 0.2, size = 1) + 
  scale_color_manual(values = c('#bebbbb', '#646262', '#e8137f')) +
  theme_classic() + 
  theme(axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(color = 'black'), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = 'none') + 
  labs(title = '', y = 'Absolute Score', x = 'Pretreatment')
plt
ggsave('../results/Figure 6A input.pdf', plt, units = 'cm', width = 10, height = 6)


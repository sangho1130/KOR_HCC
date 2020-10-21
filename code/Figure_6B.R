library(ggplot2)

data <- read.delim('../data/Figure 6B input.txt', header = TRUE, sep = '\t', row.names = 1, check.names = FALSE)
head(data)

pval <- t.test(subset(data$Tregs, data$Recurrence == 'Recurr'), subset(data$Tregs, data$Recurrence == 'NoRecurr'))$p.value

data$Recurrence <- factor(data$Recurrence, c('NoRecurr', 'Recurr'))

plt <- ggplot(data, aes(x = Recurrence, y = Tregs, col = Recurrence)) + 
  geom_jitter(width = 0.2, size = 1) + 
  scale_color_manual(values = c('royalblue', 'red2')) +
  theme_classic() + 
  theme(axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(color = 'black'), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = 'none') + 
  labs(title = 'Tregs (Pretreatment, n=21)', y = 'Absolute Score', x = 'Recurrence')
plt
ggsave('../results/Figure 6B input.pdf', plt, units = 'cm', width = 10, height = 6)


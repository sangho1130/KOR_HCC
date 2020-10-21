library(ggplot2)
library(survival)
library(GGally)

inputClinical <- '/data/Figure 1C bottom input.txt'
data <- read.delim(inputClinical, sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)

model <- Surv(data$`DFS (month)`, data$`Recurrence status`) ~ data$`Operational method`
data.surv <- survfit(model, data)
data.diff <- survdiff(model)
pval <- pchisq(data.diff$chisq, length(data.diff$n)-1, lower.tail = FALSE)
pval <- format(signif(pval, digits = 3), scientific = TRUE)
pval <- paste(c('P <', pval), collapse = ' ')

plt <- ggsurv(data.surv, surv.col = c('#ee756d', '#1aa6b8'), plot.cens = FALSE) + 
  labs(title = 'Operational methods', x = 'Months', y = 'Disease-free survival') + 
  guides(linetype = F) + 
  theme_classic(base_size = 7) + 
  theme(axis.text = element_text(colour="black"), 
        axis.ticks = element_line(colour='black'),
        plot.title = element_text(hjust=0.5, colour='black'), 
        legend.title = element_blank(), 
        legend.position=c(0.7,0.9),
        legend.background = element_blank()) + 
  annotate("text", x = 100, y = 0.1, label = pval, size = 2) + 
  ylim(0,1)
plt

ggsave('../results/Figure 1C bottom.pdf', plt, width = 6, height = 6, units = 'cm')

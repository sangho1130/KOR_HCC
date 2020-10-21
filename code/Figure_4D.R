library(survival)
library(GGally)
library(survminer)

data <- read.delim('../data/Figure 4D input.txt', sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)

model <- Surv(data$month_DFS, data$Recurrence) ~ data$Tregs
data.durv <- survfit(model, data)
data.diff <- survdiff(model)
pval <- pchisq(data.diff$chisq, length(data.diff$n)-1, lower.tail = FALSE)
pval <- as.numeric(pval)
print (pval)

cox <- coxph(Surv(month_DFS, Recurrence) ~ Tregs + Tstage + Grade + Age + AFP + Angioinvasion, 
             ties = 'efron', data = data)
cox
coxSumm <- summary(cox)
coxzph <- cox.zph(cox)
ggforest(cox, data = data)

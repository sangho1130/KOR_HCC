
data <- read.delim('../data/Figure 2B input.txt', header = T, row.names = 1, check.names = F)
colnames(data)
header <- c('CD3 T cells (execpt Tregs)', 'CD8 cytotoxic T cells', 'CD45RO memory T cells',
            'Monocytes+Macrophages', 'CD20 B cells', 'CD20 B cells + PCs', 
            'CD45 leukocytes', 'Foxp3 Tregs', 'C-kit mast cells',
            'MUM1 plasma cells', 'MPO neutrophils', 'CD163 macrophages M2')

all <- c()
all_pval <- c()
for (i in c(1:12)) {
  tmp <- na.omit(data[,c(1, i*2, i*2+1)])
  colnames(tmp) <- c('Group', 'V1', 'V2')
  tmpCorr <- cor.test(tmp$V1, tmp$V2, method = 'spearman')
  all <- append(all, c(data.frame(tmpCorr$estimate)[,1]))
  all_pval <- append(all_pval, c(tmpCorr$p.value))
}

writetable <- data.frame(SpearmanCorr = header, coef = all, pval = all_pval)
write.table(writetable, '../results/Figure 2B correlation.txt', row.names = F, col.names = T, sep = '\t', quote = F)

###
library(ggplot2)
data <- read.delim('../data/Figure 2B input.txt', header = T, row.names = 1, check.names = F)
colnames(data)
idx <- c(2, 4, 6, 8, 24, 20, 22)
celltypes <- c('CD3 T cells (except for Tregs)', 'CD8 T cells', 'Memory T cells', 
               'Macrophages', 'Macrophages M2', 'Plasma cells', 'Neutrophils')

for (i in c(1:length(idx))) {
  tmpdata <- data[, c(1, idx[i], idx[i]+1)]
  head(tmpdata)
  colnames(tmpdata) <- c('Group', 'IHC', 'CIBERSORT')
  tmpdata$log10IHC <- log10(tmpdata$IHC + 1)
  tmpdata <- na.omit(tmpdata)
  head(tmpdata)
  
  spearman <- cor.test(tmpdata$IHC, tmpdata$CIBERSORT, method = 'spearman') 
  pval <- round(spearman$p.value, 3)
  coef <- data.frame(spearman$estimate)
  coef <- coef$spearman.estimate
  coef <- round(coef, 3)
  text <- paste(c("Spearman Rho ", coef, '\n', 'P-value ', pval), collapse = '')
  text
  
  plt <- ggplot(tmpdata, aes(CIBERSORT, log10IHC)) + 
    geom_point(size=1) +
    geom_smooth(method = 'lm', se = FALSE, col = 'grey70') + 
    labs(title=celltypes[i], y = 'log10 # of cell/mm2', x = 'CIBERSORT absolute score') + 
    theme_bw(base_size = 7) +
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank()) +
    annotate("text", x = (max(tmpdata$CIBERSORT)+min(tmpdata$CIBERSORT))/2, 
             y = max(tmpdata$log10IHC)*0.95, label = text, size=2) + 
    xlim(0,max(tmpdata$CIBERSORT)) + ylim(0,max(tmpdata$log10IHC)); plt
  
  outputPdf <- paste(c('../results/corrPlot.', paste(unlist(strsplit(celltypes[i], split = ' ')), collapse=''), 
                       '.small.log10trans.pdf'), collapse = '')
  ggsave(outputPdf, plt, units = 'cm', height = 5, width = 4.5)
}


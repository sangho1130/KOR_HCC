
library(ggplot2)
library(reshape2)
library(plyr)

cibersort <- read.delim('../data/Figure 3 new CIBERSORT matrix.txt', row.names = 1, check.names = F)
head(cibersort); nrow(cibersort)
summary(cibersort$stage)

cbst_nt <- subset(cibersort, stage == 'Normal'); nrow(cbst_nt) # 15
cbst_nt <- cbst_nt[, c(1:22)]

cbst_fc <- subset(cibersort, stage == 'FibCS'); nrow(cbst_fc) # 30
cbst_fc <- cbst_fc[, c(1:22)]

cbst_dn <- subset(cibersort, stage == 'DN'); nrow(cbst_dn) # 16
cbst_dn <- cbst_dn[, c(1:22)]

cbst_t <- subset(cibersort, stage %in% c('T1', 'T2', 'T34')); nrow(cbst_t) # 54
cbst_t <- cbst_t[, c(1:22)]

cbst_t1_th <- subset(cibersort, stage == 'T1' & opmethod == 'Total_hepatectomy'); nrow(cbst_t1_th) # 9
cbst_t1_th <- cbst_t1_th[, c(1:22)]

cbst_t2_th <- subset(cibersort, stage == 'T2' & opmethod == 'Total_hepatectomy'); nrow(cbst_t2_th) # 20
cbst_t2_th <- cbst_t2_th[, c(1:22)]

cbst_t3_th <- subset(cibersort, stage == 'T34' & opmethod == 'Total_hepatectomy'); nrow(cbst_t3_th) # 6
cbst_t3_th <- cbst_t3_th[, c(1:22)]

cbst_t1_ph <- subset(cibersort, stage == 'T1' & opmethod == 'Partial_hepatectomy'); nrow(cbst_t1_ph) # 8
cbst_t1_ph <- cbst_t1_ph[, c(1:22)]

cbst_t2_ph <- subset(cibersort, stage == 'T2' & opmethod == 'Partial_hepatectomy'); nrow(cbst_t2_ph) # 6
cbst_t2_ph <- cbst_t2_ph[, c(1:22)]

cbst_t3_ph <- subset(cibersort, stage == 'T34' & opmethod == 'Partial_hepatectomy'); nrow(cbst_t3_ph) # 5
cbst_t3_ph <- cbst_t3_ph[, c(1:22)]


### Calculate median changes ###
medchange <- function(startArg, endArg, stageArg) {
  outputArg <- data.frame(matrix(ncol = 5, nrow = 22))
  colnames(outputArg) <- c('CellType', 'Stage', 'MedianChange', 'Direction', 'Pval')
  
  for (i in c(1:22)) {
    startmed <- median(startArg[, i])
    endmed <- median(endArg[, i])
    
    med_change <- endmed - startmed
    if (med_change < 0) {
      direction <- 'Decrease'
    } else if (med_change > 0) {
      direction <- 'Increase'
    } else if (med_change == 0) {
      direction <- 'Nochange'
    } 
    med_change <- abs(med_change)
    
    ttest <- t.test(startArg[, i], endArg[, i])
    if (is.nan(ttest$p.value)) {
      pval <- NA
    } else {
      pval <- ttest$p.value
    }

    now_celltype <- paste0(c(unlist(strsplit(colnames(startArg)[i], split = ' '))), collapse = '')
    
    outputArg[i, ] <- c(now_celltype, stageArg, med_change, direction, pval)
  }
  return (outputArg)
}
###

nt_fc <- medchange(cbst_nt, cbst_fc, 'NT-FibCS')
fc_dn <- medchange(cbst_fc, cbst_dn, 'FibCS-DN')
dn_t <- medchange(cbst_dn, cbst_t, 'DN-Tumor')

dn_t1_th <- medchange(cbst_dn, cbst_t1_th, 'DN-T1 TH')
t1_t2_th <- medchange(cbst_t1_th, cbst_t2_th, 'T1-T2 TH')
t2_t3_th <- medchange(cbst_t2_th, cbst_t3_th, 'T2-T3/4 TH')

dn_t1_ph <- medchange(cbst_dn, cbst_t1_ph, 'DN-T1 PH')
t1_t2_ph <- medchange(cbst_t1_ph, cbst_t2_ph, 'T1-T2 PH')
t2_t3_ph <- medchange(cbst_t2_ph, cbst_t3_ph, 'T2-T3/4 PH')


### Adaptive ###
merged_adapt <- rbind(nt_fc[1:10, ], fc_dn[1:10, ], dn_t[1:10, ], 
                      dn_t1_th[1:10, ], dn_t1_ph[1:10, ],
                      t1_t2_th[1:10, ], t1_t2_ph[1:10, ],
                      t2_t3_th[1:10, ], t2_t3_ph[1:10, ])
merged_adapt$CellType <- factor(merged_adapt$CellType, levels = c("Bcellnaive", "Bcellmemory", "PlasmaCells", "CD8", "CD4Naive", "CD4MemRest", "CD4MemAct", "Tfh", "Tregs", "Tgd"))
merged_adapt$Stage <- factor(merged_adapt$Stage, levels = c("NT-FibCS", "FibCS-DN", "DN-Tumor",
                                                            "DN-T1 TH", "T1-T2 TH", "T2-T3/4 TH",
                                                            "DN-T1 PH", "T1-T2 PH", "T2-T3/4 PH"))
merged_adapt$Direction <- factor(merged_adapt$Direction, levels = c('Decrease', 'Increase', 'Nochange'))
merged_adapt$MedianChange <- as.numeric(merged_adapt$MedianChange)
summary(merged_adapt)

plt <- ggplot(merged_adapt, aes(CellType, Stage, size = MedianChange, colour = Direction)) + 
  geom_point(alpha = .8) + 
  coord_flip() + 
  scale_color_manual(values=c('royalblue2','firebrick2', 'grey80')) + 
  scale_size_continuous(range = c(.1, 7)) + 
  theme_bw(base_size = 9) + 
  theme(text = element_text(colour = 'black', size = 9),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.ticks = element_line(colour = 'black'),
        panel.grid = element_blank()) + 
  labs(x = '', y = ''); plt
#ggsave('Figure 3 adaptive new.pdf', plt, units = 'cm', height = 8, width = 11)


### Innate ###
merged_innate <- rbind(nt_fc[11:22, ], fc_dn[11:22, ], dn_t[11:22, ], 
                      dn_t1_th[11:22, ], dn_t1_ph[11:22, ],
                      t1_t2_th[11:22, ], t1_t2_ph[11:22, ],
                      t2_t3_th[11:22, ], t2_t3_ph[11:22, ])
merged_innate$CellType <- factor(merged_innate$CellType, levels = c("NKResting", "NKActivated", "Monocytes", "MacrophagesM0", "MacrophagesM1", "MacrophagesM2",
                                                                    "DCResting", "DCActivated", "Mastcellresting", "Mastcellactivated","Eosinophils", "Neutrophils"))
merged_innate$Stage <- factor(merged_innate$Stage, levels = c("NT-FibCS", "FibCS-DN", "DN-Tumor",
                                                            "DN-T1 TH", "T1-T2 TH", "T2-T3/4 TH",
                                                            "DN-T1 PH", "T1-T2 PH", "T2-T3/4 PH"))
merged_innate$Direction <- factor(merged_innate$Direction, levels = c('Decrease', 'Increase', 'Nochange'))
merged_innate$MedianChange <- as.numeric(merged_innate$MedianChange)
summary(merged_innate)

plt <- ggplot(merged_innate, aes(CellType, Stage, size = MedianChange, colour = Direction)) + 
  geom_point(alpha = .8) + 
  coord_flip() + 
  scale_color_manual(values=c('royalblue2','firebrick2', 'grey80')) + 
  scale_size_continuous(range = c(.1, 7)) + 
  theme_bw(base_size = 9) + 
  theme(text = element_text(colour = 'black', size = 9),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.ticks = element_line(colour = 'black'),
        panel.grid = element_blank()) + 
  labs(x = '', y = ''); plt
#ggsave('Figure 3 innate new.pdf', plt, units = 'cm', height = 10, width = 11)


### Merged ###
merged_all <- rbind(nt_fc, fc_dn, dn_t, 
                       dn_t1_th, dn_t1_ph,
                       t1_t2_th, t1_t2_ph,
                       t2_t3_th, t2_t3_ph)
merged_all$CellType <- factor(merged_all$CellType, levels = c("Bcellnaive", "Bcellmemory", "PlasmaCells", "CD8", "CD4Naive", "CD4MemRest", "CD4MemAct", "Tfh", "Tregs", "Tgd",
                                                                    "NKResting", "NKActivated", "Monocytes", "MacrophagesM0", "MacrophagesM1", "MacrophagesM2",
                                                                    "DCResting", "DCActivated", "Mastcellresting", "Mastcellactivated","Eosinophils", "Neutrophils"))
merged_all$Stage <- factor(merged_all$Stage, levels = c("NT-FibCS", "FibCS-DN", "DN-Tumor",
                                                              "DN-T1 TH", "T1-T2 TH", "T2-T3/4 TH",
                                                              "DN-T1 PH", "T1-T2 PH", "T2-T3/4 PH"))
merged_all$Direction <- factor(merged_all$Direction, levels = c('Decrease', 'Increase', 'Nochange'))
merged_all$MedianChange <- as.numeric(merged_all$MedianChange)
summary(merged_all)

plt <- ggplot(merged_all, aes(CellType, Stage, size = MedianChange, colour = Direction)) + 
  geom_point(alpha = .8) + 
  coord_flip() + 
  scale_color_manual(values=c('royalblue2','firebrick2', 'grey80')) + 
  scale_size_continuous(range = c(.1, 7)) + 
  theme_bw(base_size = 9) + 
  theme(text = element_text(colour = 'black', size = 9),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.ticks = element_line(colour = 'black'),
        panel.grid = element_blank()) + 
  labs(x = '', y = ''); plt
#ggsave('Figure 3 merged new.pdf', plt, units = 'cm', height = 12, width = 10)

merged_all$padj <- p.adjust(merged_all$Pval, method = 'fdr', n = length(na.omit(merged_all$Pval)))
subset(merged_all, Pval <= 0.05)
subset(merged_all, padj <= 0.05)
#write.table(merged_all, 'Figure 3 merged new.txt', row.names = F, col.names = T, quote = F, sep = '\t')

  
  
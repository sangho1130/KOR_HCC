library(monocle)
library(plyr)

early_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    progenitorState <- table(pData(cds)$State, pData(cds)$Grade)[,"Normal"]
    return(as.numeric(names(progenitorState)[which (progenitorState == max(progenitorState))]))
  } else {
    return (1)
  }
}

label <- read.delim('../data/monocle_label.txt', sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
head(label)

expr <- read.delim('../data/monocle_expr.txt', sep = '\t', header = TRUE, row.names = 2, check.names = FALSE)
expr$ID <- NULL
identical(sort(rownames(label)), sort(colnames(expr)))


fd <- data.frame(row.names=rownames(expr), gene_short_name=rownames(expr)); head(fd)
monocleObj <- newCellDataSet(as.matrix(expr), phenoData = new("AnnotatedDataFrame", data = label), featureData=new("AnnotatedDataFrame", data = fd), expressionFamily=negbinomial.size())

pData(monocleObj)$Grade <- factor(pData(monocleObj)$Grade, levels = c('Normal', 'FibCS', 'DN', 'T1', 'T2', 'T3-4', 'Mixed'))
pData(monocleObj)$site <- factor(pData(monocleObj)$site, levels = c('Normal', 'FibCS', 'DN', 'Tumor'))
pData(monocleObj)$pretx <- factor(pData(monocleObj)$pretx, levels = c('Yes', 'No', 'Not applicable'))
pData(monocleObj)$op <- factor(pData(monocleObj)$op, levels = c('Total_hepatectomy', 'Partial_hepatectomy', 'Not applicable'))
head(pData(monocleObj))

monocleObj <- estimateSizeFactors(monocleObj)
monocleObj <- estimateDispersions(monocleObj)
monocleObj <- detectGenes(monocleObj, min_expr = 1)

head(fData(monocleObj))
expressed_genes <- row.names(subset(fData(monocleObj), num_cells_expressed >= 10))
length(expressed_genes)  # 17814 genes

###
#DEG_genes_site <- differentialGeneTest(monocleObj[expressed_genes,], fullModelFormulaStr = '~site', cores = 3)
#saveRDS(DEG_genes_site, 'tmp/DEG_genes_site.Rds')
#DEG_genes_grade <- differentialGeneTest(monocleObj[expressed_genes,], fullModelFormulaStr = '~Grade', cores = 3)
#saveRDS(DEG_genes_grade, 'tmp/DEG_genes_grade.Rds')
#DEG_genes_state <- differentialGeneTest(monocleObj[expressed_genes,], fullModelFormulaStr = '~State', cores = 3)
#saveRDS(DEG_genes_state, 'tmp/DEG_genes_state.Rds')
#DEG_genes_pseudotime <- differentialGeneTest(monocleObj[expressed_genes,], fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 3)
#saveRDS(DEG_genes_pseudotime, 'tmp/DEG_genes_pseudotime.Rds')
disp_table <- dispersionTable(monocleObj); rownames(disp_table) <- disp_table$gene_id
#saveRDS(disp_table, 'tmp/disp_table.Rds')


### Site DEGs ###
#monocleObj_ordering_genes <- row.names(DEG_genes_site)[order(DEG_genes_site$qval)][1:500]
#monocleObj_ordering_genes <- row.names(DEG_genes_site)[order(DEG_genes_site$qval)][1:1000]

### Grade DEGs ###
#monocleObj_ordering_genes <- row.names(DEG_genes_grade)[order(DEG_genes_grade$qval)][1:500]
#monocleObj_ordering_genes <- row.names(DEG_genes_grade)[order(DEG_genes_grade$qval)][1:1000]

### Dispersion genes ###
monocleObj_ordering_genes <- subset(disp_table, mean_expression >= 10 & dispersion_empirical >= 2 * dispersion_fit)$gene_id; length(monocleObj_ordering_genes) # 2127
disp_table$cols <- 'no'; disp_table[as.character(monocleObj_ordering_genes), 'cols'] <- 'in use'
ggplot(disp_table, aes(log10(mean_expression+1), dispersion_empirical, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'grey80'))


monocleObj <- setOrderingFilter(monocleObj, ordering_genes = monocleObj_ordering_genes)
monocleObj <- reduceDimension(monocleObj, method = 'DDRTree')
monocleObj <- orderCells(monocleObj)
monocleObj <- orderCells(monocleObj, root_state = early_state(monocleObj))
head(pData(monocleObj))
#saveRDS(monocleObj, 'tmp/monocleObj.Rds')


plot_cell_trajectory(monocleObj, color_by = 'Grade') + scale_color_manual(values = c('Normal'='#bdc1c8', 'FibCS'='#67814f', 'DN'='#f2ae2e', 'T1'='#eca9a4', 'T2'='#e65251', 'T3-4'='#963e38', 'Mixed'='#7a41c4'))
ggsave('../results/disp.grade.pdf', units = 'cm', width = 9, height = 11)

plot_cell_trajectory(monocleObj, color_by = 'site') + scale_color_manual(values = c('Normal'='#bdc1c8', 'FibCS'='#67814f', 'DN'='#f2ae2e', 'Tumor'='#e65251'))
ggsave('../results/disp.site.orig.pdf', units = 'cm', width = 9, height = 10)

plot_cell_trajectory(monocleObj, color_by = 'op') + scale_color_manual(values = c('Total_hepatectomy'='#00add3', 'Partial_hepatectomy'='#e68d83', 'Not applicable'='grey90'))
ggsave('../results/disp.operation.pdf', units = 'cm', width = 9, height = 10)

plot_cell_trajectory(monocleObj, color_by = 'Tregs') + scale_color_gradient(low = 'grey90', high = 'red2')
ggsave('../results/disp.Tregs.pdf', units = 'cm', width = 9, height = 10)

plot_cell_trajectory(monocleObj, color_by = 'AbsoluteScore') + scale_color_gradient(low = 'grey90', high = 'red2')
ggsave('../results/disp.AbsoluteScore.pdf', units = 'cm', width = 9, height = 10)

plot_cell_trajectory(monocleObj, color_by = 'Pseudotime')
ggsave('../results/disp.pseudotime.pdf', units = 'cm', width = 9, height = 10)


### below takes serious amount of time ###

#BEAM_branch1 <- BEAM(monocleObj, branch_point = 1, cores = 3)
#BEAM_branch1 <- BEAM_branch1[order(BEAM_branch1$qval),]
#BEAM_branch1 <- BEAM_branch1[,c("gene_short_name", "pval", "qval")]
head(BEAM_branch1)
#saveRDS(BEAM_branch1, 'tmp/BEAM_branch1.Rds')

#BEAM_branch2 <- BEAM(monocleObj, branch_point = 2, cores = 3)
#BEAM_branch2 <- BEAM_branch2[order(BEAM_branch2$qval),]
#BEAM_branch2 <- BEAM_branch2[,c("gene_short_name", "pval", "qval")]
head(BEAM_branch2)
#saveRDS(BEAM_branch2, 'tmp/BEAM_branch2.Rds')

# Branch 1 
usegenes <- row.names(subset(BEAM_branch1, qval < 1e-5)); length(usegenes)
BEAM_branch1_res <- plot_genes_branched_heatmap(monocleObj[usegenes,], cores = 2, branch_point = 1, #num_clusters = ,
                                                use_gene_short_name = T, show_rownames = T, return_heatmap = T)
BEAM_branch1_genes <- data.frame(gene = rownames(BEAM_branch1_res$annotation_row), BEAM_branch1_res$annotation_row)
#write.table(BEAM_branch1_genes, 'clustering/BEAM_branch1_genes_qval.txt', sep = '\t', quote = F, row.names = F, col.names = T)
#pdf('clustering/BEAM_branch1_pheatmap_qval.pdf')
BEAM_branch1_res$ph_res
dev.off()
BEAM_branch1_res <- plot_genes_branched_heatmap(monocleObj[usegenes,], cores = 2, branch_point = 1, num_clusters = 2,
                                                use_gene_short_name = T, show_rownames = T, return_heatmap = T)
BEAM_branch1_genes <- data.frame(gene = rownames(BEAM_branch1_res$annotation_row), BEAM_branch1_res$annotation_row)
#write.table(BEAM_branch1_genes, 'clustering/BEAM_branch1_genes_qval_v2.txt', sep = '\t', quote = F, row.names = F, col.names = T)
#pdf('clustering/BEAM_branch1_pheatmap_qval_v2.pdf')
BEAM_branch1_res$ph_res
dev.off()

usegenes <- row.names(BEAM_branch1)[1:100]
BEAM_branch1_res <- plot_genes_branched_heatmap(monocleObj[usegenes,], cores = 2, branch_point = 1, #num_clusters = ,
                                                use_gene_short_name = T, show_rownames = T, return_heatmap = T)
BEAM_branch1_genes <- data.frame(gene = rownames(BEAM_branch1_res$annotation_row), BEAM_branch1_res$annotation_row)
#write.table(BEAM_branch1_genes, 'clustering/BEAM_branch1_genes.txt', sep = '\t', quote = F, row.names = F, col.names = T)
#pdf('clustering/BEAM_branch1_pheatmap.pdf')
BEAM_branch1_res$ph_res
dev.off()

plot_cell_trajectory(monocleObj, markers = 'SRPX', use_color_gradient = T) # fate 1
plot_cell_trajectory(monocleObj, markers = 'KRT19', use_color_gradient = T) # fate 1
plot_cell_trajectory(monocleObj, markers = 'CTCFL', use_color_gradient = T) # fate 2
plot_cell_trajectory(monocleObj, markers = 'TERT', use_color_gradient = T) # fate 2


# Branch 2
usegenes <- row.names(subset(BEAM_branch2, qval < 1e-5)); length(usegenes)
BEAM_branch2_res <- plot_genes_branched_heatmap(monocleObj[usegenes,], cores = 2, branch_point = 2, #num_clusters = ,
                                                use_gene_short_name = T, show_rownames = T, return_heatmap = T)
BEAM_branch2_genes <- data.frame(gene = rownames(BEAM_branch2_res$annotation_row), BEAM_branch2_res$annotation_row)
#write.table(BEAM_branch2_genes, 'clustering/BEAM_branch2_genes_qval.txt', sep = '\t', quote = F, row.names = F, col.names = T)
#pdf('clustering/BEAM_branch2_pheatmap_qval.pdf')
BEAM_branch2_res$ph_res
dev.off()

BEAM_branch2_res <- plot_genes_branched_heatmap(monocleObj[usegenes,], cores = 2, branch_point = 2, num_clusters = 4,
                                                use_gene_short_name = T, show_rownames = T, return_heatmap = T)
BEAM_branch2_genes <- data.frame(gene = rownames(BEAM_branch2_res$annotation_row), BEAM_branch2_res$annotation_row)
#write.table(BEAM_branch2_genes, 'clustering/BEAM_branch2_genes_qval_v2.txt', sep = '\t', quote = F, row.names = F, col.names = T)
#pdf('clustering/BEAM_branch2_pheatmap_qval_v2.pdf')
BEAM_branch2_res$ph_res
dev.off()

usegenes <- row.names(BEAM_branch2)[1:100]
BEAM_branch2_res <- plot_genes_branched_heatmap(monocleObj[usegenes,], cores = 2, branch_point = 2, #num_clusters = ,
                                                use_gene_short_name = T, show_rownames = T, return_heatmap = T)
BEAM_branch2_genes <- data.frame(gene = rownames(BEAM_branch2_res$annotation_row), BEAM_branch2_res$annotation_row)
#write.table(BEAM_branch2_genes, 'clustering/BEAM_branch2_genes.txt', sep = '\t', quote = F, row.names = F, col.names = T)
#pdf('clustering/BEAM_branch2_pheatmap.pdf')
BEAM_branch2_res$ph_res
dev.off()

plot_cell_trajectory(monocleObj, markers = 'PAGE2', use_color_gradient = T) # fate 1
plot_cell_trajectory(monocleObj, markers = 'TERT', use_color_gradient = T) # fate 1
plot_cell_trajectory(monocleObj, markers = 'MUC6', use_color_gradient = T) # fate 2
plot_cell_trajectory(monocleObj, markers = 'MZB1', use_color_gradient = T) # fate 2
plot_cell_trajectory(monocleObj, markers = 'C7', use_color_gradient = T) # fate 2



###
head(DEG_genes_pseudotime)
sig_gene_names <- row.names(subset(DEG_genes_pseudotime, qval < 1e-5)); length(sig_gene_names)
DEG_genes_pseudotime_res <- plot_pseudotime_heatmap(monocleObj[sig_gene_names,], #num_clusters = 3, 
                                                    cores = 2, show_rownames = T, return_heatmap = T)
DEG_genes_pseudotime_res
DEG_genes_pseudotime_genes <- data.frame(gene = rownames(DEG_genes_pseudotime_res$annotation_row), DEG_genes_pseudotime_res$annotation_row)
#write.table(DEG_genes_pseudotime_genes, 'clustering/DEG_genes_pseudotime_genes.txt', sep = '\t', quote = F, row.names = F, col.names = T)
#pdf('clustering/DEG_genes_pseudotime_pheatmap.pdf')
DEG_genes_pseudotime_res$ph_res
dev.off()

DEG_genes_pseudotime_res$annotation_row

#save.image('monocle2.Rdata')
#load('monocle2.Rdata')

head(t(monocleObj@reducedDimS))
test <- data.frame(t(monocleObj@reducedDimS))
colnames(test) <- c('Component 1', 'Component 2')
test$ImmuneCluster <- pData(monocleObj)$ImmuneCluster
test$op <- pData(monocleObj)$op
ggplot(test, aes(`Component 1`, `Component 2`, shape = pData(monocleObj)$ImmuneCluster, col = pData(monocleObj)$op)) +
  geom_point() +
  scale_color_manual(values = c('Total_hepatectomy'='#00add3', 'Partial_hepatectomy'='#e68d83', 'Not applicable'='grey90')) +
  scale_shape_manual(values = c('C1' = 20, 'C2' = 17, 'None' = 20)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = 'None')
# 
plot_cell_trajectory(monocleObj, color_by = 'op', cell_size = 0) + 
  geom_point(aes(shape = pData(monocleObj)$ImmuneCluster, col = pData(monocleObj)$op)) +
  scale_color_manual(values = c('Total_hepatectomy'='#00add3', 'Partial_hepatectomy'='#e68d83', 'Not applicable'='grey90')) +
  scale_shape_manual(values = c('C1' = 20, 'C2' = 17, 'None' = 20))
#ggsave('trajectory/disp.ImmuneCluster.pdf', units = 'cm', width = 9, height = 10)

label_w <- pData(monocleObj)
label_w <- data.frame(Sample = rownames(label_w), label_w)
head(label_w)
#write.table(label_w, 'monocle2.txt', quote = F, sep = '\t', row.names = F, col.names = T)

plot_cell_trajectory(monocleObj, color_by = 'State')
#ggsave('trajectory/disp.State.pdf', units = 'cm', width = 9, height = 10)


head(pData(monocleObj))

t.test(subset(pData(monocleObj), State == '4')$Tregs,
       subset(pData(monocleObj), State != '4')$Tregs) # 0.01269
t.test(subset(pData(monocleObj), State == '4' & op =='Total_hepatectomy')$Tregs,
       subset(pData(monocleObj), State != '4' & op =='Total_hepatectomy')$Tregs) # 0.05815

t.test(subset(pData(monocleObj), State == '4')$AbsoluteScore,
       subset(pData(monocleObj), State == '3')$AbsoluteScore) # 0.0005227
t.test(subset(pData(monocleObj), State == '4')$AbsoluteScore,
       subset(pData(monocleObj), State != '4')$AbsoluteScore) # 0.00000000000000022
t.test(subset(pData(monocleObj), State == '4' & op =='Total_hepatectomy')$AbsoluteScore,
       subset(pData(monocleObj), State != '4' & op =='Total_hepatectomy')$AbsoluteScore) # 0.00001821


###
test <- subset(pData(monocleObj), State %in% c('1', '5'));head(test)
test <- subset(test, Grade %in% c('FibCS', 'DN'));head(test)

ggplot(test, aes(State, AbsoluteScore, fill = State)) +
  geom_boxplot() +
  scale_fill_manual(values = c('1' = '#f7756d', '5' = '#e76bf3')) +
  facet_wrap(~Grade) +
  labs(x = 'State', y = 'Absolute score') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'None',
        axis.text = element_text(colour = 'black'))
#ggsave('trajectory/absolutescore_state1_state5.pdf', units = 'cm', width = 8, height = 6)

t.test(subset(test, State == '1' & Grade == 'FibCS')$AbsoluteScore,
       subset(test, State == '5' & Grade == 'FibCS')$AbsoluteScore) # 0.001076
t.test(subset(test, State == '1' & Grade == 'DN')$AbsoluteScore,
       subset(test, State == '5' & Grade == 'DN')$AbsoluteScore) # 0.009992


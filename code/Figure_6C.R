library(RGtk2)
library(rsgcc)
library(RColorBrewer)
library(ggplot2)
library(scales)

drPiechart <- function(columnNames, Values, Colors, outputPdf){
  data <- data.frame(
    group = columnNames,
    value = Values
  )
  data$group <- factor(data$group, columnNames)
  
  pie <- ggplot(data, aes(x="", y=value, fill=factor(group))) +
    geom_bar(width = 1, stat = "identity") + 
    coord_polar("y", start=0) + 
    scale_fill_manual(values = Colors) + 
    coord_polar(theta = "y", direction = -1) +
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank()) 
  
  ggsave(outputPdf, units = 'cm', height = 8, width = 16)
}

#drPiechart(c(), c(), colors = c(), 'output.pdf')
drPiechart(c('Yes', 'No', 'NA'), 
           c(21, 13, 1), 
           c('#e8137f', '#646262', '#bebbbb'),
           '../results/Figure 5C left.pdf')


inputCibersort <- '../data/Figure 6C right input.txt'
data <- read.delim(inputCibersort, sep = '\t', header = TRUE, check.names = FALSE, row.names = 1)
data <- data[, -c(1, 3)]
head(data)
my_palette <- colorRampPalette(c('grey90', "steelblue", 'grey30', "red2", 'grey90'))(n = 12)

heatmap.2(as.matrix(data), trace = 'none', na.color = 'grey90',
          Rowv = FALSE, Colv = FALSE, dendrogram = 'none', scale = 'none', cexCol = 1,
          breaks = c(-0.12, -0.1, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12), 
          col = my_palette, density.info = 'none', margins = c(4,18), 
          labCol = c('Abs.', 'Rel.'))
dev.off()

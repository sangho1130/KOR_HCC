  
drPiechart <- function(columnNames, Values, Colors, outputPdf){
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
  
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

  ggsave(outputPdf, pie, units = 'cm', height = 8, width = 16)
}

#drPiechart(c(), c(), colors = c(), 'output.pdf')
drPiechart(c('Normal', 'FL', 'FH', 'CS', 'DL', 'DH', 'T1', 'T2', 'T3', 'Mixed'), 
           c(16, 10, 10, 10, 10, 7, 17, 30, 11, 7), 
           c('#bebdbd', '#bbe165', '#6e8a3c', '#546a2e', 
             '#f1c055', '#eb8919', '#f69693', '#f7474e', '#aa0c0b', '#570a08'), 
           '../results/Figure 1A.pdf')



  
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
drPiechart(c('Total hepatectomy', 'Partial hepatectomy'), 
           c(35, 19), 
           c('#1aa6b8', '#ee756d'),
           '../results/Figure 1C.pdf')



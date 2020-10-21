#install.packages("extrafont")
library(ggplot2)
library(reshape2)

### Adaptive
inFile <- '../data/Figure 3 input adaptive.txt'; plotTitle <- 'Adaptive Immunity'

data <- read.delim(inFile, sep = '\t', header = TRUE, check.names = FALSE)
data$Stage <- factor(data$Stage, levels = unique(data$Stage))
data$CellType <- factor(data$CellType, levels = unique(data$CellType))

plt <- ggplot(data, aes(CellType, Stage, size = MedianChange, colour = Direction)) + 
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
  labs(title=plotTitle, x = '', y = '')
plt

ggsave('../results/Figure 3 input adaptive.pdf', plt, units = 'cm', height = 15, width = 11)


### Innate
inFile <- '../data/Figure 3 input innate.txt'; plotTitle <- 'Innate Immunity'

data <- read.delim(inFile, sep = '\t', header = TRUE, check.names = FALSE)
data$Stage <- factor(data$Stage, levels = unique(data$Stage))
data$CellType <- factor(data$CellType, levels = unique(data$CellType))

plt <- ggplot(data, aes(CellType, Stage, size = MedianChange, colour = Direction)) + 
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
  labs(title=plotTitle, x = '', y = '')
plt

ggsave('../results/Figure 3 input innate.pdf', plt, units = 'cm', height = 15, width = 11)


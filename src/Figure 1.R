library(dplyr)
library(readr)
library(paletteer)

data <- read_tsv('data/flowsom_clustered.txt')
panel <- colnames(data)[1:49]
d <- data %>% filter(subtype1 != 'T cell')
q99 = quantile(as.matrix(d[,1:length(panel)]), 0.99)
pheatmap::pheatmap(as.matrix(d[,1:length(panel)]),
                   color = paletteer_d("rcartocolor::PurpOr", n = 100, type = 'continuous'),
                   border_color = 'white',
                   scale = 'none',
                   labels_row = paste(d$cell_cluster,' (', round(d$percentage,1), '%', ')',sep = ''),
                   display_numbers = FALSE,
                   angle_col = 45,
                   breaks = seq(0, q99, q99/90),
                   main = "Median marker expression per cluster",
                   #filename = 'figures/flowSOM_annotated_remove_Tcells.pdf',
                   height = 15, width = 18)



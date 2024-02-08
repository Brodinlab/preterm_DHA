library(dplyr)
library(readr)
library(paletteer)
library(vite)
library(ggraph)
library(igraph)
library(ggpubr)

# Figure 1b ----------------
res_bind = read.table('data/flowsom_clustered.txt', header = TRUE, sep = "\t", check.names = FALSE)
res_bind_Tcell = read.table('data/flowsom_clustered_Tcell.txt', header = TRUE, sep = "\t", check.names = FALSE)
common_marker = colnames(res_bind)[1:49]
data <- rbind(
  res_bind %>% filter(subtype1 != 'T cell') %>% 
    mutate(lineage = subtype1) %>%
    select(all_of(common_marker), n, cell_cluster, lineage, clusters),
  res_bind_Tcell %>% 
    filter(subtype3 != 'Exclude') %>% 
    mutate(cell_cluster = paste0(subtype3, '_', clusters),
           lineage = if_else(subtype3 == 'gdT', 'gdT', 'abT')) %>% 
    select(common_marker, n, cell_cluster, lineage, clusters)
) %>%
  mutate(percentage = n/sum(n) * 100)

set.seed(824)
G <- NULL
G <- build_graph(data, common_marker, filtering_T = 5)
for (i in names(data)) G <- igraph::set.vertex.attribute(G, 
                                                         name = i, value = data[, i])
cc <- igraph::multilevel.community(G)
V(G)$community_id <- as.character(cc$membership)
l_g = create_layout(G, layout="fr")
l_g$x <- V(G)$x
l_g$y <- V(G)$y
ggraph(l_g) +
  geom_edge_link(alpha=.1) + 
  geom_node_point(aes(size= percentage, fill = lineage), shape=21,  alpha = 1) +
  paletteer::scale_fill_paletteer_d("ggthemes::stata_s2color", direction = 1) + ## New palette
  scale_size_continuous(range = c(10, 20)) +
  geom_node_text(aes(label = clusters), size = 4, color = 'white', show.legend = FALSE, repel = F) ## Show label within the node
ggsave(filename = 'figures/Figure 1/1b.pdf', width = 16, height = 15)

# Figure 1c ----------------
data <- read_tsv('data/flowsom_clustered_Tcell.txt')
panel <- colnames(data)[1:49]
q99 <- quantile(as.matrix(data[,1:length(panel)]), 0.99)

pheatmap::pheatmap(as.matrix(data[,1:length(panel)]),
                   color = paletteer::paletteer_d("rcartocolor::PurpOr", n = 100, type = 'continuous'), # New color
                   border_color = 'white',
                   scale = 'none',
                   labels_row = paste(data$subtype3, '_', data$clusters,' (', round(data$percentage,1), '%', ')',sep = ''),
                   display_numbers = FALSE,
                   angle_col = 45,
                   breaks = seq(0, q99, q99/90),
                   main = "Median marker expression per cluster (Z-score of all subtypes)",
                   filename = 'figures/Figure 1/1c.pdf',
                   height = 8, width = 10)


# Figure 1d ----------------
data <- read_tsv('data/flowsom_clustered.txt')
panel <- colnames(data)[1:49]
data <- data %>% filter(subtype1 != 'T cell')
q99 = quantile(as.matrix(data[,1:length(panel)]), 0.99)
pheatmap::pheatmap(as.matrix(data[,1:length(panel)]),
                   color = paletteer_d("rcartocolor::PurpOr", n = 100, type = 'continuous'),
                   border_color = 'white',
                   scale = 'none',
                   labels_row = paste(data$cell_cluster,' (', round(data$percentage,1), '%', ')',sep = ''),
                   display_numbers = FALSE,
                   angle_col = 45,
                   breaks = seq(0, q99, q99/90),
                   main = "Median marker expression per cluster",
                   filename = 'figures/Figure 1/1d.pdf',
                   height = 15, width = 18)

# Figure 1e,f ----------------
data <- read_csv('data/flowSOM_frequency_table.csv')
m <- data %>% 
  tidyr::pivot_wider(id_cols = ID_unique, names_from = merge2, values_from = freq, values_fill = NA)
mds_result <- cmdscale(dist(m[,-1]), k=4) %>% as.data.frame()
mds_result$ID_unique <- m$ID_unique
df <- left_join(mds_result, 
                data %>% select(ID_unique, Timepoint_string) %>% distinct()) %>%
  mutate(Timepoint_string = factor(Timepoint_string, levels = c('Day 1', 'Day 3', 'Day 7', 'Day 14', 'Day 28', 'PMA 40 weeks')))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Oranges"))
ggplot(df, aes(x=V1, y=V2, color=Timepoint_string)) + geom_point() + scale_color_manual(values = getPalette(6)) + theme_bw() + xlab('MDS1') + ylab('MDS2')
ggsave('figures/Figure 1/1e.pdf', width = 5, height=4)

subpop_list <- c('Neutrophils', 'Monocytes', 'NK', 'B cell', 'CD4 T', 'CD8 T', 'gdT', 'DPT', 'Basophils', 'Eosinophil', 'mDC', 'platelets', 'LinNeg')
m_fill0 <- data %>% 
  tidyr::pivot_wider(id_cols = ID_unique, names_from = merge2, values_from = freq, values_fill = 0)
cor_m <- df %>% left_join(m_fill0) %>% select(V1, all_of(subpop_list)) %>% cor()
g1 <- data.frame(name=rownames(cor_m), cor = cor_m[,'V1']) %>%
  filter(name != 'V1') %>%
  arrange(cor) %>%
  mutate(name = factor(name, levels=.$name)) %>%
  ggplot(aes(x=cor, y=name)) +
  geom_segment(aes(x=0, xend=cor, y=name, yend=name), color="grey") +
  geom_point(color="orange", size=4) +
  ggtitle('MDS1') + 
  theme_bw()

cor_m <- df %>% left_join(m_fill0) %>% select(V2, all_of(subpop_list)) %>% cor()
g2 <- data.frame(name=rownames(cor_m), cor = cor_m[,'V2']) %>%
  filter(name != 'V2') %>%
  arrange(cor) %>%
  mutate(name = factor(name, levels=.$name)) %>%
  ggplot(aes(x=cor, y=name)) +
  geom_segment(aes(x=0, xend=cor, y=name, yend=name), color="grey") +
  geom_point(color="orange", size=4) +
  ggtitle('MDS2') + 
  theme_bw()
ggarrange(g1, g2)
ggsave('figures/Figure 1/1f.pdf', width = 8, height=8)







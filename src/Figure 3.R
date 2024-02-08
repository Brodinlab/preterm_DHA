library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggrepel)
library(lme4)
library(merTools)

olink <- read_csv('data/Olink_72_inds.csv')
olink$visit_name = factor(olink$visit_name, levels = c('Day 1', 'Day 3', 'Day 7', 'Day 14', 'Day 28', 'PMA 40 weeks'))
olink <- olink %>% mutate(Breastmilk_group=case_when(Breastmilk_group %in% c('mixed', 'maternal_milk') ~ 'maternal milk or mixed',
                                                     TRUE ~ Breastmilk_group))


# Figure 3a ----------------
olink_panel <- colnames(olink[,-c(1,2,3)])
mds_result <- cmdscale(dist(olink[,olink_panel]), k=4) %>% as.data.frame()
mds_result$visit_name <- olink$visit_name
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
ggplot(mds_result, aes(x=V1, y=V2, color=visit_name)) + geom_point() + scale_color_manual(values = getPalette(6)) + theme_bw() + 
  ggtitle('Olink protein level') + xlab('MDS1') + ylab('MDS2')
ggsave('figures/Figure 3/3a.pdf', width = 5, height=4)

# Figure 3b ----------------
ggtmp <- function(df){
  rbind(df %>% slice_max(cor, n=30),
        df %>% slice_min(cor, n=30)) %>%
    ggplot(aes(x=cor, y=name)) +
    geom_segment(aes(x=0, xend=cor, y=name, yend=name), color="grey") +
    geom_point(color="lightblue", size=4) +
    theme_pubclean()
}
cor_m1 <- olink %>% mutate(MDS1 = mds_result$V1) %>% dplyr::select(MDS1, all_of(olink_panel)) %>% cor(use = 'complete.obs')
cor_df1 <- data.frame(name=rownames(cor_m1), cor = cor_m1[,'MDS1']) %>%
  filter(name != 'MDS1') %>%
  arrange(cor) %>%
  mutate(name = factor(name, levels=.$name))
g1 <- ggtmp(cor_df1) + ggtitle('Correlation with MDS1')

cor_m2 <- olink %>% mutate(MDS2 = mds_result$V2) %>% dplyr::select(MDS2, all_of(olink_panel)) %>% cor(use = 'complete.obs')
cor_df2 <- data.frame(name=rownames(cor_m2), cor = cor_m2[,'MDS2']) %>%
  filter(name != 'MDS2') %>%
  arrange(cor) %>%
  mutate(name = factor(name, levels=.$name))
g2 <- ggtmp(cor_df2) + ggtitle('Correlation with MDS2')
ggarrange(g1, g2)
ggsave('figures/Figure 3/3b.pdf', width = 12, height = 12)

# Figure 3c ----------------
d <- lapply(unique(olink$visit_name), function(x){
  d <- olink %>% filter(visit_name == x) %>% dplyr::select(all_of(olink_panel)) %>% dist() %>% as.matrix()
  return(data.frame(visit_name = x, pdist = d[upper.tri(d)]))
}) %>% do.call(rbind, .)

d.stats <- d %>% t_test(pdist ~ visit_name, p.adjust.method = 'fdr') %>% add_xy_position(x = 'visit_name')

ggplot(d, aes(x=visit_name, y=pdist, fill=visit_name))  + geom_violin() +  
  stat_summary(fun = "mean",
               geom = "point",
               color = "black") + 
  scale_fill_manual(values = getPalette(6)) + 
  theme_bw() + 
  ggtitle('Interindividual distances of plasma protein level')
ggsave('figures/Figure 3/3c.pdf', width = 5, height=4)  

# Figure 3d ----------------
olink_mean <- olink %>% group_by(visit_name) %>% summarise(across(all_of(olink_panel), median))
olink_fc <- t(olink_mean[,-c(1)]) %>% as_tibble()
colnames(olink_fc) <- olink_mean$visit_name
olink_fc$name = olink_panel
olink_fc <- olink_fc %>% as_tibble%>% mutate(FC = log2(`PMA 40 weeks`+1) - log2(`Day 1`+1)) %>%
  arrange(abs(FC)) %>%
  mutate(name = factor(name, levels=.$name))
label_name <- rbind(olink_fc %>% slice_max(FC, n=21),
                    olink_fc %>% slice_min(FC, n=21))$name
olink_fc <- olink_fc %>% mutate(label = if_else(name %in% label_name, name, NA))
ggplot(olink_fc, aes(x=name, y = FC)) + geom_point() + 
  geom_label_repel(aes(label=label), box.padding = 0.5, max.overlaps = Inf) + 
  theme_pubclean() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle('Plasma protein logFC of Day 100 vs Day 1')
ggsave('figures/Figure 3/3d.pdf', width = 7.5, height = 7.5)

# Figure 3f ----------------
freq_anno <- read_csv('data/flowSOM_frequency_table.csv')
freq_wide <- tidyr::pivot_wider(freq_anno, id_cols = 'ID_unique', names_from = 'merge2', values_from = 'freq', values_fill = 0)
ID_link <- freq_anno %>% dplyr::select(Subject_ID, Timepoint_string, ID_unique, PN_days) %>% distinct() %>% mutate(unique_link = paste(Subject_ID, Timepoint_string, sep = '_'))
freq_names <- unique(freq_anno$merge2)

data <- olink %>% mutate(unique_link = paste(Subject_ID, visit_name, sep = '_')) %>% dplyr::select(-c(Subject_ID, visit_name)) %>%
  left_join(ID_link, by='unique_link', suffix = c('', '.y')) %>%
  left_join(freq_wide, by='ID_unique', suffix = c('', '.y')) %>%
  filter(!is.na(Neutrophils)) # filter out those we don't have CyTOF data

data.scale <- cbind(scale(data[,olink_panel]), 
                    scale(data[,freq_names]), 
                    data %>% dplyr::select(PN_days, Timepoint_string, Subject_ID))

raw_cor <- cor(data[,olink_panel], data[,freq_names], method = 'pearson', use = 'complete.obs')
cor <- raw_cor[apply(abs(raw_cor) > 0.3, 1, any), c('B cell', 'NK', 'Neutrophils', 'Monocytes', 'CD4 T', 'CD8 T')] # filter out cytokine with no correlation at all

lmer_list <- lapply(rownames(cor), function(x){
  lmer(as.formula(paste0('`', x, '` ~ `B cell` + `CD4 T` + `CD8 T` + `NK` + `Neutrophils` + `Monocytes` + (1|PN_days) + (1|Subject_ID)')), data.scale)
})
names(lmer_list) <- rownames(cor)
lmer_coeffsig <- lapply(lmer_list, function(m) {
  tmp <- Anova(m)
  tmp$subpop <- rownames(tmp)
  return(tmp %>% mutate(name = (m@call$formula %>% as.character)[2]) %>% as_tibble %>% left_join(FEsim(m, 1000) %>% rename(subpop=term), by='subpop'))
}) %>% do.call(rbind, .)
lmer_coeffsig$padj <- p.adjust(lmer_coeffsig$`Pr(>Chisq)`, method = 'BH')

df <- lmer_coeffsig %>% mutate(padj_sig = if_else(padj<0.05, padj, NA)) %>%
  mutate(name = gsub('`', '', name)) %>% filter(!is.na(padj_sig))
ggplot(df, aes(x=subpop, y=name, size=-log2(padj_sig), color=median)) + geom_point() + theme_bw() + ggtitle('Mixed effect model') + 
  scale_color_gradient2(low = scales::muted("blue"),
                        mid = "white",
                        high = scales::muted("red"))+
  labs(color='coefficient')
ggsave('figures/Figure 3/3f.pdf')




library(dplyr)
library(readr)
library(ggpubr)
library(ggrepel)
library(wesanderson)
library(rstatix)
library(igraph)

subpop_list <- c('Neutrophils', 'Monocytes', 'NK', 'B cell', 'CD4 T', 'CD8 T', 'gdT', 'DPT', 'Basophils', 'Eosinophil', 'mDC', 'platelets', 'LinNeg')
majorpop_list <- c('Neutrophils', 'Monocytes', 'NK', 'B cell', 'CD4 T', 'CD8 T')

data <- read_csv('data/flowSOM_frequency_table.csv')
cdata.raw <- read_csv('data/Olin 2018 frequency table refined.csv')
cdata_subpop_list <- c('B cell', 'CD4 T', 'CD8 T', 'Monocytes', 'NK', 'Neutrophils', 'Plasmablasts', 'Basophils', 'pDC', 'Endothelial_cells', 'mDC')
cdata <- cdata.raw %>% 
  select(Subject_ID, ID_unique, Timepoint_string, GA_days_at_birth, Relation, class, PN_days, Gender, 
         GA_days, all_of(cdata_subpop_list)) %>%
  tidyr::pivot_longer(cols = all_of(cdata_subpop_list), names_to = 'merge2', values_to = 'freq') %>%
  filter(
    class=='Control',
    Relation == 'Child') %>%
  mutate(
    class=if_else(class=='Control', 'Term', class),
    group='Term control',
    Breastmilk_group='Term control')
data.merge <- plyr::rbind.fill(data,cdata) %>% 
  filter(Breastmilk_group != 'Unknown') %>%
  mutate(milk_group = case_when(
    Breastmilk_group == 'donated' ~ 'donated',
    Breastmilk_group %in% c('mixed', 'maternal_milk') ~ 'maternal milk or mixed',
    Breastmilk_group == 'Term control' ~ 'Term control'
  ))

# Figure 5b ----------------
lapply(majorpop_list, function(x) {
  data.merge %>% 
    filter(merge2 == x) %>%
    ggplot(aes(x=PN_days, y=freq, color=milk_group)) +
    geom_point(size=0.3, alpha=0.3)+
    geom_line(aes(group = Subject_ID), alpha=0.3) +
    geom_smooth(method = 'loess', span=1.5) + 
    scale_x_continuous(trans='log1p',
                       breaks=c(0,1,3,7,14,28,100)) +
    scale_colour_manual(values=wes_palette("Darjeeling1")) +
    ggtitle(x)+
    theme_bw()
}) %>% ggarrange(plotlist = ., ncol = 3, nrow = 5, common.legend = T) %>%
  ggexport(filename='figures/Figure 5/5b.pdf', width = 10, height = 12)

# Figure 5c ----------------
lapply(majorpop_list, function(x) {
  d <- data.merge %>% 
    filter(merge2 == x, PN_days >80) 
  d.stats <- d %>%
    t_test(freq ~ milk_group, p.adjust.method = "fdr") %>%
    add_xy_position(x = "milk_group")
  ggplot(d, aes(x=milk_group, y=freq, color=milk_group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size=0.3)+
    stat_pvalue_manual(d.stats, tip.length = 0.01, label='p')+
    scale_colour_manual(values=wes_palette("Darjeeling1")) +
    ggtitle(x)+
    theme_pubr()
}) %>% ggarrange(plotlist = ., ncol = 2, nrow = 3, common.legend = T) %>%
  ggexport(filename='figures/Figure 5/5c.pdf', width = 10, height = 12)

# Figure 5d ----------------
df <- data.merge %>% 
  filter(PN_days>80,
         merge2 %in% c('Neutrophils', 'Monocytes', 'NK', 'B cell', 'CD4 T', 'CD8 T'),
         class != 'Term') %>%
  select(merge2, ID_unique, freq, group, class, milk_group, PN_days, Breastmilk_group) %>%
  tidyr::pivot_wider(names_from = merge2, values_from = freq, values_fill = 0)
subfreq <- df %>% select(all_of(c('Neutrophils', 'Monocytes', 'NK', 'B cell', 'CD4 T', 'CD8 T')))
set.seed(42)
permanova.res <- vegan::adonis2(formula=subfreq ~ milk_group * class, data=df, method='bray', permutations = 30000)
data.frame(factor=rownames(permanova.res), R2 = permanova.res$R2, sig=round(permanova.res$`Pr(>F)`, digits = 4)) %>%
  filter(factor!='Total') %>%
  ggplot(aes(x = '100 days', y=R2, fill=factor, label=sig)) +
  geom_bar(stat = 'identity') +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_bw()
ggsave('figures/Figure 5/5d.pdf', width = 4, height = 8)

# Figure 5e ----------------
data.wide <- data.merge %>%
  tidyr::pivot_wider(names_from = merge2, values_from = freq, values_fill = 0)
d <- data.wide[,c('Neutrophils', 'Monocytes', 'NK', 'B cell', 'CD4 T', 'CD8 T')]
G1 <- cccd::nng(d, k=4, method='Bray')
set.seed(100)
l <- layout_with_fr(G1)
pal <- wes_palette("Zissou1", 100, type = "continuous")
range1.100 <- function(x){1 + 99*(x-min(x))/(max(x)-min(x))}
ticks<- c(0,1,3,7,14,28,100)
pdf('figures/Figure 5/5e.pdf', width = 10, height = 4)
par(mfrow=c(1,3))
lapply(unique(data.wide$milk_group), function(x){
  data.wide.color <- data.wide %>% 
    mutate(v_color=case_when(
      milk_group == x ~ pal[round(range1.100(log2(PN_days+1)))],
      TRUE ~ '#A9A9A9'
    ))
  V(G1)$color=data.wide.color$v_color
  plot(G1, layout=l,
       vertex.size=8,
       vertex.label=NA,
       vertex.frame.color=NA,
       edge.arrow.size=0)
  fields::image.plot(zlim=range(log2(ticks+1)), legend.only=T, col=pal, axis.args=list(at=log2(ticks+1), labels=ticks))
  title(x)
})
dev.off()

# Figure 5f ----------------
olink <- read_csv('data/Olink_72_inds.csv') %>%
  mutate(Breastmilk_group=case_when(Breastmilk_group %in% c('mixed', 'maternal_milk') ~ 'maternal milk or mixed',
                                    TRUE ~ Breastmilk_group))
olink$visit_name = factor(olink$visit_name, levels = c('Day 1', 'Day 3', 'Day 7', 'Day 14', 'Day 28', 'PMA 40 weeks'))
sample_info <- read_csv(file.path('data', 'sample_info_refined.csv'))%>%
  mutate(Timepoint_string = factor(Timepoint_string, levels = c('Day 1', 'Day 3', 'Day 7', 'Day 14', 'Day 28', 'PMA 40 weeks')))

df <- left_join(olink %>% filter(visit_name == 'PMA 40 weeks'),
                sample_info %>% select(Subject_ID, preterm_class) %>% distinct(), by='Subject_ID') %>%
  tidyr::pivot_longer(cols = -c(Subject_ID, visit_name, Breastmilk_group, preterm_class))

p_stats <- df %>%
  group_by(name) %>%
  wilcox_test(value ~ Breastmilk_group) %>%
  adjust_pvalue(method = "BH")


logFC = df %>% group_by(name, Breastmilk_group) %>% summarise(value=mean(value, na.rm=T)) %>% group_by(name) %>%
  summarise(fc = value[2]-value[1]) %>%
  arrange(fc)
top_fc_list <- rbind(slice_max(logFC, fc, n=20),
                     slice_min(logFC, fc, n=20)) %>% select(name)
tmp <- logFC %>% mutate(label=if_else(name %in% unlist(top_fc_list), name, ''))
ggplot(tmp, aes(x=factor(name, levels = name), y=fc)) + 
  geom_point() + 
  geom_text_repel(aes(label=label), box.padding = 0.5, max.overlaps = Inf) + 
  theme_pubclean() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylab('log fold change (maternal milk/mixed vs donated)')
ggsave('figures/Figure 5/5f.pdf', width = 14, height = 9)




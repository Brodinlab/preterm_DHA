library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(rstatix)

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
data.merge <- plyr::rbind.fill(data,cdata)

# Figure 2a ----------------
lapply(majorpop_list, function(x) {
  tmp <- data.merge %>% filter(merge2 == x, GA_days - GA_days_at_birth <2) %>% select(freq)
  yhigh <- max(tmp)*1.1
  ylow <- min(tmp)*0.9
  data.merge %>% filter(merge2 == x,GA_days - GA_days_at_birth <2, class != 'Term') %>%
    mutate(GA_weeks = GA_days_at_birth/7) %>%
    ggplot(aes(x=GA_weeks, y=freq)) + 
    geom_point() + 
    geom_smooth(method = 'lm') + 
    stat_cor(method = "spearman") +
    theme_bw() + 
    ylim(ylow, yhigh) + 
    ggtitle(x)
}) %>% ggarrange(plotlist = ., ncol = 2, nrow = 3) %>%
  ggexport(filename='figures/Figure 2/2a.pdf', width = 10, height = 12)

# Figure 2b ----------------
lapply(majorpop_list, function(x) {
  tmp <- data.merge %>% filter(merge2 == x, GA_days - GA_days_at_birth <2) %>% select(freq)
  yhigh <- max(tmp)*1.1
  ylow <- min(tmp)*0.9
  data.merge %>% filter(merge2 == x,GA_days - GA_days_at_birth <2, class == 'Term') %>%
    mutate(GA_weeks = GA_days_at_birth/7) %>%
    ggplot(aes(x=GA_weeks, y=freq)) + 
    geom_point() + 
    theme_bw() + 
    ylim(ylow, yhigh) + 
    ggtitle(x)
}) %>% ggarrange(plotlist = ., ncol = 2, nrow = 3) %>%
  ggexport(filename='figures/Figure 2/2b.pdf', width = 4, height = 12)

# Figure 2c ----------------
lapply(majorpop_list, function(x) {
  data.merge %>% 
    filter(merge2 == x, class != 'Term') %>%
    ggplot(aes(x=PN_days, y=freq)) +
    geom_point(size=0.3, alpha=0.3)+
    geom_line(aes(group = Subject_ID), alpha=0.3) +
    geom_smooth(method = 'loess', color='red') + 
    scale_x_continuous(trans='log1p',
                       breaks=c(0,1,3,7,14,28,100)) +
    ggtitle(x)+
    theme_bw()
}) %>% ggarrange(plotlist = ., ncol = 3, nrow = 5, common.legend = T) %>%
  ggexport(filename='figures/Figure 2/2c.pdf', width = 10, height = 12)


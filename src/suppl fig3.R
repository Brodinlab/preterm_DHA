library(dplyr)
library(readr)
library(ggpubr)
library(ggrepel)
library(rstatix)
library(wesanderson)

subpop_list <- c("Neutrophils", "Monocytes", "NK", "B cell", "CD4 T", "CD8 T", "gdT", "DPT", "Basophils", "Eosinophil", "mDC", "platelets", "LinNeg")
data <- read_csv("data/flowSOM_frequency_table.csv")
cdata.raw <- read_csv("data/Olin 2018 frequency table refined.csv")
cdata_subpop_list <- c("B cell", "CD4 T", "CD8 T", "Monocytes", "NK", "Neutrophils", "Plasmablasts", "Basophils", "pDC", "Endothelial_cells", "mDC")
cdata <- cdata.raw %>%
    select(
        Subject_ID, ID_unique, Timepoint_string, GA_days_at_birth, Relation, class, PN_days, Gender,
        GA_days, all_of(cdata_subpop_list)
    ) %>%
    tidyr::pivot_longer(cols = all_of(cdata_subpop_list), names_to = "merge2", values_to = "freq") %>%
    filter(
        class == "Control",
        Relation == "Child"
    ) %>%
    mutate(
        class = if_else(class == "Control", "Term", class),
        group = "Term control",
        Breastmilk_group = "Term control"
    )
frac_breastmilk <- read_delim("data/frac_breastmilk.txt") %>%
    select(Infant_ID, frac_breastmilk_own) %>%
    rename(Subject_ID = Infant_ID) %>%
    mutate(Subject_ID = as.character(Subject_ID))

data.merge <- plyr::rbind.fill(data, cdata) %>%
    left_join(frac_breastmilk, by = "Subject_ID")

data.merge.breastmilk <- data.merge %>%
    filter(Breastmilk_group != "Term control", Breastmilk_group != "Unknown") %>%
    mutate(Breastmilk_group = factor(Breastmilk_group, levels = c("donated", "mixed", "maternal_milk")))

lapply(subpop_list, function(x) {
    d <- data.merge.breastmilk %>%
        filter(merge2 == x, PN_days > 80, Breastmilk_group != "Term control", Breastmilk_group != "Unknown")
    d.stats <- d %>%
        t_test(freq ~ Breastmilk_group, p.adjust.method = "fdr") %>%
        add_xy_position(x = "Breastmilk_group")
    ggplot(d, aes(x = Breastmilk_group, y = freq, color = Breastmilk_group)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(size = 0.3) +
        stat_pvalue_manual(d.stats, tip.length = 0.01, label = "p") +
        scale_colour_manual(values = wes_palette("BottleRocket2")) +
        ggtitle(x) +
        theme_pubr()
}) %>%
    ggarrange(plotlist = ., ncol = 2, nrow = 7, common.legend = TRUE) %>%
    ggexport(filename = "figures/suppl figure 3/c.pdf", width = 8, height = 18)


lapply(subpop_list, function(x) {
    data.merge.breastmilk %>%
        filter(merge2 == x, PN_days > 80, Breastmilk_group != "Term control", Breastmilk_group != "Unknown") %>%
        ggplot(aes(x = frac_breastmilk_own, y = freq)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        stat_cor(method = "spearman") +
        ggtitle(x) +
        theme_pubr()
}) %>%
    ggarrange(plotlist = ., ncol = 2, nrow = 7, common.legend = TRUE) %>%
    ggexport(filename = "figures/suppl figure 3/d.pdf", width = 8, height = 18)



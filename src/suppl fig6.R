library(dplyr)
library(readr)
library(ggpubr)
library(ggrepel)
library(wesanderson)
library(rstatix)
library(igraph)

subpop_list <- c("Neutrophils", "Monocytes", "NK", "B cell", "CD4 T", "CD8 T", "gdT", "DPT", "Basophils", "Eosinophil", "mDC", "platelets", "LinNeg")
majorpop_list <- c("Neutrophils", "Monocytes", "NK", "B cell", "CD4 T", "CD8 T")

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
data.merge <- plyr::rbind.fill(data, cdata) %>%
    filter(Breastmilk_group != "Unknown") %>%
    mutate(milk_group = case_when(
        Breastmilk_group == "donated" ~ "donated",
        Breastmilk_group %in% c("mixed", "maternal_milk") ~ "maternal milk or mixed",
        Breastmilk_group == "Term control" ~ "Term control"
    ))

# supplemental figure ----------------
lapply(majorpop_list, function(x) {
    d <- data.merge %>%
        filter(merge2 == x, PN_days > 80)

    # # Stats for NEC comparison within class
    # d.stats.nec <- d %>%
    #   filter(!is.na(NEC)) %>%
    #   group_by(milk_group) %>%
    #   t_test(freq ~ NEC, p.adjust.method = 'fdr') %>%
    #   add_xy_position(x = "milk_group", dodge = 0.75)

    # Stats for class comparison within NEC=0
    d.stats.milk <- d %>%
        filter(!is.na(NEC), NEC == 0) %>%
        t_test(freq ~ milk_group, p.adjust.method = "fdr") %>%
        add_xy_position(x = "milk_group")

    # Adjust xmin and xmax for d.stats.class to align with left boxplots
    d.stats.milk <- d.stats.milk %>%
        mutate(
            xmin = xmin - 0.1875, # Shift left by half the dodge width
            xmax = xmax - 0.1875
        ) # Shift left by half the dodge width

    ggplot(d %>% filter(!is.na(NEC)), aes(x = milk_group, y = freq, color = as.factor(NEC))) +
        geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1) +
        # stat_pvalue_manual(d.stats.nec, tip.length = 0.01, label = 'p') +
        stat_pvalue_manual(d.stats.milk,
            y.position = max(d$freq, na.rm = TRUE) * 1.1,
            tip.length = 0,
            label = "p",
            color = "black"
        ) +
        scale_colour_manual(values = wes_palette("Darjeeling1")) +
        ggtitle(x) +
        theme_bw()
}) %>%
    ggarrange(plotlist = ., ncol = 2, nrow = 3, common.legend = T) %>%
    ggexport(filename = "figures/suppl figure 6/6b.pdf", width = 10, height = 12)

# supplemental figure ----------------
lapply(majorpop_list, function(x) {
    d <- data.merge %>%
        filter(merge2 == x, PN_days > 80)

    # Stats for BPD comparison within class
    d.stats.bpd <- d %>%
        filter(!is.na(BPD)) %>%
        group_by(milk_group) %>%
        t_test(freq ~ BPD, p.adjust.method = "fdr") %>%
        add_xy_position(x = "milk_group", dodge = 0.75)

    # Stats for class comparison within BPD=0
    d.stats.milk <- d %>%
        filter(!is.na(BPD), BPD == 0) %>%
        t_test(freq ~ milk_group, p.adjust.method = "fdr") %>%
        add_xy_position(x = "milk_group")

    # Adjust xmin and xmax for d.stats.class to align with left boxplots
    d.stats.milk <- d.stats.milk %>%
        mutate(
            xmin = xmin - 0.1875, # Shift left by half the dodge width
            xmax = xmax - 0.1875
        ) # Shift left by half the dodge width

    ggplot(d %>% filter(!is.na(BPD)), aes(x = milk_group, y = freq, color = as.factor(BPD))) +
        geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1) +
        stat_pvalue_manual(d.stats.bpd, tip.length = 0.01, label = "p") +
        stat_pvalue_manual(d.stats.milk,
            y.position = max(d$freq, na.rm = TRUE) * 1.1,
            tip.length = 0,
            label = "p",
            color = "black"
        ) +
        scale_colour_manual(values = wes_palette("Darjeeling1")) +
        ggtitle(x) +
        theme_bw()
}) %>%
    ggarrange(plotlist = ., ncol = 2, nrow = 3, common.legend = T) %>%
    ggexport(filename = "figures/suppl figure 6/6a.pdf", width = 10, height = 12)

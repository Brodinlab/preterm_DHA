library(dplyr)
library(ggplot2)
library(readr)

fileInfo <- read.csv("data/sample_info_refined.csv") %>% mutate(
    Timepoint_string = factor(Timepoint_string, levels = c("Day 1", "Day 3", "Day 7", "Day 14", "Day 28", "PMA 40 weeks")),
    gender = if_else(SEX == 1, "boy", "girl"),
    preterm_class = factor(preterm_class, levels = c("25-28w", "22-25w"))
)

fileInfo %>% filter(Timepoint_string == "PMA 40 weeks") %>% pull(Subject_ID) %>% unique() %>% length()

ggplot(
    fileInfo %>% mutate(Subject_ID = factor(Subject_ID, levels = unique((fileInfo %>% arrange(preterm_class, gender))$Subject_ID))),
    aes(x = PN_days, y = Subject_ID)
) +
    lemon::geom_pointline(aes(group = Subject_ID, color = preterm_class, shape = gender)) +
    geom_vline(xintercept = c(1.5, 5.5, 12, 20, 45), color = "blue", size = 0.5, alpha = 0.5) +
    theme_classic() +
    xlab("PN days") +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )
ggsave(filename = "figures/suppl figure 1/suppl_fig1.pdf", width = 8, height = 8)


fileInfo_term <- read_csv("data/Olin 2018 frequency table refined.csv") %>%
    mutate(
        gender = if_else(Gender == 1, "boy", "girl"),
        gender = if_else(is.na(gender), "unknown", gender)
    ) %>%
    filter(
        class == "Control",
        Relation == "Child"
    )

ggplot(
    fileInfo_term %>% mutate(Subject_ID = factor(Subject_ID, levels = unique((fileInfo_term %>% arrange(gender))$Subject_ID))),
    aes(x = PN_days, y = Subject_ID)
) +
    lemon::geom_pointline(aes(group = Subject_ID, shape = gender)) +
    geom_vline(xintercept = c(1.5, 5.5, 12, 20, 45), color = "blue", size = 0.5, alpha = 0.5) +
    theme_classic() +
    xlab("PN days") +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )
ggsave(filename = "figures/suppl figure 1/suppl_fig1_term.pdf", width = 8, height = 8)

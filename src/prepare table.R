library(dplyr)
library(readr)


raw <- read_csv(file.path("data", "sample_info_refined.csv"))
bpd <- read_csv(file.path("data", "PetterData_BPD.csv")) %>%
    select(ID, BPD_yn) %>%
    rename(Subject_ID = ID)
nec <- read_delim(file.path("data", "NEC_info.csv"), delim = ";") %>%
    rename(Subject_ID = Inf_ID)

sample_info <- raw %>%
    left_join(bpd, by = "Subject_ID") %>%
    left_join(nec, by = "Subject_ID")

write_csv(sample_info, file.path("data", "sample_info_refined.csv"))

raw2 <- read_csv(file.path("data", "flowSOM_frequency_table.csv"))
df <- raw2 %>%
    mutate(NEC = ifelse(is.na(NEC), 0, NEC)) %>%
    mutate(BPD = ifelse(is.na(BPD_yn), 0, BPD_yn))
write_csv(df, file.path("data", "flowSOM_frequency_table.csv"))

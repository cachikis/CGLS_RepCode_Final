###############################################
# Filename: bds.R
# Author: Craig A. Chikis
# Date: 08/02/2023
# Note(s):
###############################################
rm(list = ls(all.names = TRUE))
gc()

library(tidyverse)


bds = read_csv("Data/bds2020_fa (3).csv")
entry_rate = read_csv("Data/bds2020 (2).csv")


entry_rate_summ = entry_rate %>%
    filter(year >= 1960 & year <= 2019) %>%
    summarise(median = median(estabs_entry_rate, na.rm = TRUE)) 

bds1 = bds %>%
    select(year, fage, emp) %>%
    mutate(emp = as.numeric(emp)) %>%
    group_by(year) %>%
    mutate(share = emp / sum(emp, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(fage %in% c("a) 0", "b) 1")) %>%
    group_by(year) %>%
    summarise(share = sum(share, na.rm = TRUE)) %>%
    ungroup() %>%
    summarise(median = median(share, na.rm = TRUE)) 


bds3 = bds %>%
    select(year, fage, emp) %>%
    mutate(emp = as.numeric(emp)) %>%
    group_by(year) %>%
    mutate(share = emp / sum(emp, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(fage %in% c("a) 0", "b) 1", "c) 2", "d) 3")) %>%
    group_by(year) %>%
    summarise(share = sum(share, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(year >= 1980) %>%
    summarise(median = median(share, na.rm = TRUE)) 

write_csv(bds1, "Output/Store_Data/bds1.csv")
write_csv(bds3, "Output/Store_Data/bds3.csv")
write_csv(entry_rate_summ, "Output/Store_Data/entry_rate_summ.csv")

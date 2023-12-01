#######################################################
# Filename: rates_of_return.R
# Author: Craig A. Chikis
# Date: 06/27/2023
# Note(s):
######################################################
rm(list = ls(all.names = TRUE))
gc()

library(quantmod)
library(fixest)
library(DescTools)
library(tidyverse)


args = commandArgs()
wd = args[str_detect(args, "CGLS_RepCode_Final")]

setwd(wd)



romme_data = read_csv("Output/Store_Data/usdata.csv")

romme_replicate = romme_data %>%
    filter(lubridate::year(DATE) >= 1984 & lubridate::year(DATE) <= 2000) %>%
    summarise(all_pretax = mean(return_all_capital_pre_tax, na.rm = TRUE),
              all_pretax_nogain = mean(return_all_capital_pre_tax_no_gain, na.rm = TRUE),
              bus_pretax = mean(return_business_capital_pre_tax_adj, na.rm = TRUE),
              bus_pretax_nogain = mean(return_business_capital_pre_tax_no_gain_adj, na.rm = TRUE))

romme_data = romme_data %>%
    mutate(year = lubridate::year(DATE)) %>%
    arrange(year, desc(DATE)) %>%
    distinct(year, .keep_all = TRUE) %>%
    select(year, return_all_capital_pre_tax, return_business_capital_pre_tax_adj, 
           return_all_capital_after_tax, return_business_capital_after_tax_adj) %>%
    pivot_longer(-year, names_to = "series", values_to = "values")

dgs1 = read_csv("Output/Store_Data/dgs1.csv")

dgs1_summ = dgs1 %>%
    inner_join(romme_data, by = "year") 
    
dgs1_summ = dgs1_summ %>%
    mutate(chi = (values/100 - real_rf) / 2) %>%
    group_by(series) %>%
    summarise(mean = 100*mean(chi),
              median = 100*median(chi),
              mean_rk =mean(values),
              median_rk = median(values),
              mean_rf = mean(real_rf),
              median_rf = median(real_rf)) %>%
    ungroup()

write_csv(dgs1_summ, "Output/Store_Data/chi_rk.csv")
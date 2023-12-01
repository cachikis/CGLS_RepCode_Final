############################################################
# Filename: gomme_process_IA.R
# Author: Craig A. Chikis
# Date: 10/09/2023
# Notes():
############################################################
rm(list = ls(all.names = TRUE))
gc()


library(viridis)
library(quantmod)
library(tidyverse)


args = commandArgs()
wd = args[str_detect(args, "CGLS_RepCode_Final")]

setwd(wd)

getSymbols(c("DGS5", "DGS10", "DGS1"), src = "FRED")

DGS1 = as.data.frame(DGS1) %>%
    rownames_to_column() %>% 
    as_tibble() %>%
    mutate(date = lubridate::ymd(rowname)) %>% 
    select(date, DGS1) %>% 
    filter(!is.na(date) & !is.na(DGS1)) %>% 
    mutate(month = floor_date(date, "month")) %>% 
    arrange(month, desc(date)) %>% 
    distinct(month, .keep_all = TRUE) %>% 
    select(month, yield = DGS1) %>% 
    mutate(month = ceiling_date(month, "month") - days(1)) %>%
    mutate(series = "One year")

DGS10 = as.data.frame(DGS10) %>%
    rownames_to_column() %>% 
    as_tibble() %>%
    mutate(date = lubridate::ymd(rowname)) %>% 
    select(date, DGS10) %>% 
    filter(!is.na(date) & !is.na(DGS10)) %>% 
    mutate(month = floor_date(date, "month")) %>% 
    arrange(month, desc(date)) %>% 
    distinct(month, .keep_all = TRUE) %>% 
    select(month, yield = DGS10) %>% 
    mutate(month = ceiling_date(month, "month") - days(1)) %>%
    mutate(series = "Ten year")

DGS5 = as.data.frame(DGS5) %>%
    rownames_to_column() %>% 
    as_tibble() %>%
    mutate(date = lubridate::ymd(rowname)) %>% 
    select(date, DGS5) %>% 
    filter(!is.na(date) & !is.na(DGS5)) %>% 
    mutate(month = floor_date(date, "month")) %>% 
    arrange(month, desc(date)) %>% 
    distinct(month, .keep_all = TRUE) %>% 
    select(month, yield = DGS5) %>% 
    mutate(month = ceiling_date(month, "month") - days(1)) %>%
    mutate(series = "Five year")

yields = bind_rows(DGS1, DGS10, DGS5) %>% 
    mutate(series = factor(series, levels = c("One year", "Five year", "Ten year")))
rm(DGS1, DGS10, DGS5)


mich1 = read_csv("Data/michigan_1_inflation_monthly.csv", skip = 1) %>%
    mutate(Month = as.character(Month),
           Month = if_else(nchar(Month) == 1, paste0("0", Month), Month)) %>% 
    mutate(date = lubridate::ymd(paste0(Year, "-", Month, "-", "01"))) %>% 
    select(month = date, mean = Mean, median = Median) %>% 
    mutate(month = ceiling_date(month, "month") - days(1)) %>% 
    mutate(series = "One year")

mich510 = read_csv("Data/michigan_510_inflation_monthly.csv", skip = 1) %>%
    mutate(Month = as.character(Month),
           Month = if_else(nchar(Month) == 1, paste0("0", Month), Month)) %>% 
    mutate(date = lubridate::ymd(paste0(Year, "-", Month, "-", "01"))) %>% 
    select(month = date, mean = Mean, median = Median) %>% 
    mutate(month = ceiling_date(month, "month") - days(1)) %>% 
    mutate(series = "Five year")

mich510 = mich510 %>% 
    bind_rows(mutate(mich510, series = "Ten year")) %>% 
    bind_rows(mich1) %>%
    mutate(series = factor(series, levels = c("One year", "Five year", "Ten year")))

yields = left_join(yields, mich510, by = c("month", "series")) %>%
    arrange(series, month)
rm(mich1, mich510)

yields = yields %>% 
    filter(!is.na(mean) & !is.na(yield)) %>% 
    mutate(real_yield = yield - median)


plotyields = ggplot(data = yields, mapping = aes(x = month, y = real_yield/100, color = series)) + 
    geom_line(linewidth = 1.25) + 
    geom_hline(mapping = aes(yintercept = 0)) + 
    scale_color_viridis(discrete = TRUE) +
    labs(x = "Year", y = "Real yield", color = "") + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")


gomme = read_csv("Data/gomme.csv") %>% 
    select(date = DATE, return_business_capital_after_tax_no_gain) %>% 
    na.omit() %>% 
    mutate(date = ceiling_date(date, "quarter") - days(1))

returns = ggplot(data = filter(gomme, lubridate::year(date) >= 1980), mapping = aes(x = date, y = return_business_capital_after_tax_no_gain/100)) + 
    geom_line(linewidth = 1.25, color = "black") + 
    scale_y_continuous(limits = c(0, 0.11), labels = scales::percent,
                       breaks = seq(0, 0.11, by = 0.01)) + 
    scale_color_viridis(discrete = TRUE) + 
    labs(x = "", y = "%") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")

gomme = gomme %>% 
    mutate(series = "Real return, after tax, no capital gain") %>% 
    mutate(series2 = "Return on capital") %>% 
    dplyr::rename(real_yield = return_business_capital_after_tax_no_gain,
                  month = date)

yields = yields %>% 
    mutate(series2 = "Risk-free rate")

yields = bind_rows(yields, gomme)


yields_er = yields %>% 
    filter(series == "One year") %>%
    mutate(quarter = floor_date(month, "quarter")) %>% 
    arrange(quarter, desc(month)) %>% 
    distinct(quarter, .keep_all = TRUE) %>% 
    left_join(select(gomme, month, rk = real_yield), by = c("month")) %>% 
    mutate(excess_return = rk - real_yield) %>%
    mutate(series = "Excess return",
           series2 = series) %>%
    select(month, series, series2, real_yield = excess_return)

yields = bind_rows(yields, yields_er) 

write_csv(yields, "Output/Store_Data/yields.csv")


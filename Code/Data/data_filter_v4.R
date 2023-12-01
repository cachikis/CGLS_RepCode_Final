#################################################
# Filename: data_filter_v4.R
# Author: Craig A. Chikis
# Date: 08/11/2022
# Note(s):
#################################################
rm(list = ls(all.names = TRUE))
gc()

set.seed(20220719)


packages = c("DescTools", "lfe", "data.table", "haven", "lubridate", "tidyverse")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))
}

library(DescTools)
library(lfe)
library(data.table)
library(haven)
library(lubridate)
library(tidyverse)


args = commandArgs()
wd = args[str_detect(args, "CGLS_RepCode_Final")]

setwd(wd) 



cstat_a0 = read_csv("Data/cstat_a.csv", guess_max = 1e6)
names(cstat_a0) = tolower(names(cstat_a0))

cstat_q0 = read_dta("Data/cstat_q.dta")
names(cstat_q0) = tolower(names(cstat_q0))

kpss0 = haven::read_dta("Data/MergedPatentData.dta")

kpss0 = kpss0 %>%
    mutate(permco = as.integer(permco)) %>%
    mutate(fcit010 = exp(lfcit010) - 1,
            fcit05 = exp(lfcit05) - 1,
            fcit210 = fcit25 + fcit610,
            fcit620 = fcit610 + fcit1120,
            fcit020 = fcit05 + fcit610 + fcit1120, 
            fcit2A = fcitALL - fcit01,
            fcit6A = fcitALL - fcit01 - fcit25,
            fcit11A = fcitALL - fcit01 - fcit25 - fcit610,
            lfcit2A  = log(1+fcitALL - fcit01),
            lfcit6A  = log(1+fcitALL - fcit01 - fcit25),
            lfcit11A = log(1+fcitALL - fcit01 - fcit25 - fcit610)) %>%
    select(contains("myclass"), permco, wku, fcit020, issue_year, filed_year, everything())

nber = map(.x = c(citeNBER = "Data/cite76_06.dta", assgNBER = "Data/pat76_06_assg.dta", ipcNBER = "Data/pat76_06_ipc.dta"), 
           .f = read_dta)
nberlink = read_dta("Data/pdpcohdr.dta")
dynamiclink = read_dta("Data/dynass.dta")

kogan0 = read_csv("Data/KPSS_2019_Public.csv")
permno_permco = read_csv("Data/Match_patent_permco_permno_2022.csv.zip") %>%
    distinct(patent_num, permco)

kogan0 = kogan0 %>%
    mutate(issue_date = lubridate::mdy(issue_date),
           patent_num = as.integer(patent_num),
           permco = as.integer(permco)) %>%
    mutate(y = lubridate::year(issue_date)) %>%
    arrange(permco, y)


highmarkup = read_csv("Data/deloecker_etal_targets.csv")


ntile_modify = function(x, b) {

    outvec = rep(NA, length(x))
    ntile_res = ntile(na.omit(x), b)
    outvec[!is.na(x)] = ntile_res

    return(outvec)
}

bin_quintiles = function(df0, bucketvar, nbucket = 5) {
    bucketvar_tm1 = paste0(bucketvar, "_tm1")
    bucketvar_t1 = paste0(bucketvar, "_t1")

    df = df0 %>%
        group_by(year) %>%
        mutate(bucket = ntile_modify(desc(.data[[bucketvar]]), nbucket),
               bucket_tm1 = ntile_modify(desc(.data[[bucketvar_tm1]]), nbucket),
               bucket_t1 = ntile_modify(desc(.data[[bucketvar_t1]]), nbucket)) %>%
        ungroup()

    return(df)

}


leadlag_var = function(df0, horizon) {

    name_list = list(num = c(1,2,3,4,5,-1,-2,-3,-4,-5),
                     label = c(paste0("_t", 1:5), paste0("_tm", 1:5)))

    idx = which(name_list$num == horizon) 

    salename = paste0("sale", name_list$label[idx])
    profitname = paste0("profit", name_list$label[idx])

    sales_t1 = df0 %>%
        select(lpermco, year, saleholder = sale, profitholder = profit) %>%
        mutate(year = year - horizon) 

    names(sales_t1)[names(sales_t1) == "saleholder"] = salename
    names(sales_t1)[names(sales_t1) == "profitholder"] = profitname

    df0 = left_join(df0, sales_t1, by = c("lpermco", "year"))

    return(df0)

}

tabulate_by_quintile = function(df0, col, type = "") {

  by_quintile_t = df0 %>%
        filter(!is.na(bucket)) %>% 
        group_by(bucket, year) %>%
        summarise(sd = sd(.data[[col]], na.rm = TRUE),
                  p50 = median(.data[[col]], na.rm = TRUE)) %>%
        ungroup() %>%
        group_by(bucket) %>%
        summarise(sd = mean(sd, na.rm = TRUE),
                  p50 = mean(p50, na.rm = TRUE)) %>%
        ungroup()

    by_quintile_t_uncond = df0 %>%
        group_by(year) %>%
        summarise(sd = sd(.data[[col]], na.rm = TRUE),
                  p50 = median(.data[[col]], na.rm = TRUE)) %>%
        ungroup() %>%
        summarise(sd = mean(sd, na.rm = TRUE),
                  p50 = mean(p50, na.rm = TRUE)) %>%
        mutate(bucket = 0) %>%
        ungroup()

    by_quintile_t = bind_rows(by_quintile_t, by_quintile_t_uncond)

    if (type == "sd") { 
        by_quintile_t = by_quintile_t %>%
            select(bucket, sd)
    } else if (type == "median") { 
        by_quintile_t = by_quintile_t %>%
            select(bucket, p50)
    }

    return(by_quintile_t)

}

salesmatfunc = function(bucket, bucket_t1, df0, nbucket = "quintile") {
    if (nbucket == "quintile") {
        df = df0 %>% 
            summarise(`11` = sum(!!sym(bucket) == 1 & !!sym(bucket_t1) == 1, na.rm = TRUE) / 
                             sum(!!sym(bucket) == 1, na.rm = TRUE),
                      `21` = sum(!!sym(bucket) == 2 & !!sym(bucket_t1) == 1, na.rm = TRUE) / 
                             sum(!!sym(bucket) == 2, na.rm = TRUE),
                      `31` = sum(!!sym(bucket) == 3 & !!sym(bucket_t1) == 1, na.rm = TRUE) / 
                             sum(!!sym(bucket) == 3, na.rm = TRUE),
                      `41` = sum(!!sym(bucket) == 4 & !!sym(bucket_t1) == 1, na.rm = TRUE) / 
                             sum(!!sym(bucket) == 4, na.rm = TRUE),
                      `51` = sum(!!sym(bucket) == 5 & !!sym(bucket_t1) == 1, na.rm = TRUE) / 
                             sum(!!sym(bucket) == 5, na.rm = TRUE))

        transmat = rbind(c(df[["11"]]), 
                         c(df[["21"]]), 
                         c(df[["31"]]), 
                         c(df[["41"]]), 
                         c(df[["51"]]) 
                         )
    } else if (nbucket == "halves") {
        df = df0 %>%
            summarise(`11` = sum(!!sym(bucket) == 1 & !!sym(bucket_t1) == 1, na.rm = TRUE),
                    `12` = sum(!!sym(bucket) == 1 & !!sym(bucket_t1) == 2, na.rm = TRUE),
                    `21` = sum(!!sym(bucket) == 2 & !!sym(bucket_t1) == 1, na.rm = TRUE),
                    `22` = sum(!!sym(bucket) == 2 & !!sym(bucket_t1) == 2, na.rm = TRUE)) %>%
            mutate(TOT1 = `11` + `12`, 
                   TOT2 = `21` + `22`) %>%
            mutate(`11` = `11`/TOT1,
                    `12` = `12`/TOT1,
                    `21` = `21`/TOT2,
                    `22` = `22`/TOT2)

        transmat = rbind(c(df[["11"]], df[["12"]]),
                         c(df[["21"]], df[["22"]])
                         )
        
    } else if (nbucket == "top10/bottom90") {

        df = df0 %>% 
            summarise(`11` = sum(!!sym(bucket) & !!sym(bucket_t1), na.rm = TRUE),
                    `12` = sum(!!sym(bucket) & !(!!sym(bucket_t1)), na.rm = TRUE),
                    `21` = sum(!(!!sym(bucket)) & !!sym(bucket_t1), na.rm = TRUE),
                    `22` = sum(!(!!sym(bucket)) & !(!!sym(bucket_t1)), na.rm = TRUE)) %>%
            mutate(TOT1 = `11` + `12`, 
                   TOT2 = `21` + `22`) %>%
            mutate(`11` = `11`/TOT1,
                    `12` = `12`/TOT1,
                    `21` = `21`/TOT2,
                    `22` = `22`/TOT2)

        transmat = rbind(c(df[["11"]], df[["12"]]),
                         c(df[["21"]], df[["22"]])
                         )
                         


    }

    return(transmat)
}


nqtr_index = function(cstat_q0, numqtrs) {

    cstat_q = cstat_q0 %>%
        select(datadate, lpermco, saleq, saley, cogsq, cogsy, everything()) %>%
        mutate(year = lubridate::year(datadate),
               lpermco = as.integer(lpermco)) 

    cstat_q = cstat_q %>%
        mutate(saleq = ifelse(saleq < 0, NA, saleq),
               cogsq = ifelse(cogsq < 0, NA, cogsq)) %>%
        mutate(quarter = floor_date(datadate, "quarter"),
               profit = saleq - cogsq) %>%
        arrange(lpermco, year, quarter)

    cstat_q = cstat_q %>% 
        distinct(lpermco, quarter, .keep_all = TRUE)

    cstat_q_t1 = cstat_q %>%
        select(lpermco, quarter, profit_t1 = profit) %>%
        mutate(quarter = floor_date(quarter - months(3), "quarter"))

    cstat_q_tm1 = cstat_q %>%
        select(lpermco, quarter, profit_tm1 = profit) %>%
        mutate(quarter = floor_date(quarter + months(3), "quarter"))

    cstat_q = left_join(cstat_q, cstat_q_t1, by = c("lpermco", "quarter")) %>%
        left_join(cstat_q_tm1, by = c("lpermco", "quarter")) %>%
        select(lpermco, year, quarter, profit, profit_t1, profit_tm1, everything()) %>%
        arrange(lpermco, year, quarter)

    if (is.numeric(numqtrs)) {
        idx = cstat_q %>%
            group_by(lpermco, year) %>%
            summarise(count = n(),
                      numpos = sum(profit > 0 & !is.na(profit)) >= numqtrs) %>%
            ungroup()
    } 

    if (any(idx$count > 4)) {
        stop("More than 4 quarters of data for a given permco.")
    }

    findissue = idx %>%
        filter(count > 4)

    idx = idx %>%
        select(lpermco, year, numpos) %>%
        distinct(lpermco, year, .keep_all = TRUE)

    return(idx)

}
restrict_qtr_cstat = map(.x = list(three = 3),
                         .f = ~nqtr_index(cstat_q0, .x))

# https://github.com/tidymodels/butcher/issues/141
small_felm = function(m) {
    m$fe <- NULL
    m$residuals <- NULL
    m$r.residuals <- NULL
    m$fitted.values <- NULL
    m$response <- NULL
    m$cfactor <- NULL
    m$inv <- NULL
    m$STATS$promotion <- NULL
    m$clustervar <- NULL
    return(m)
}

clean_cstat = function(cstat_a0, profit_restrict, uniformsample = TRUE, winsor_drop = c(-Inf, Inf),
                       winsor_vec = c(0.01, 0.99),
                       ficUSA = TRUE, nbucket_cstat = 5, nbucket_patent = 5, add_filters = "profit", 
                       remove_sic = list(sic2 = c(), sic4 = c()), 
                       fe = "sic2",
                       posRD = TRUE) {
    
    cstat_a = cstat_a0 %>%
        select(datadate, lpermco, sic, at, xrd, cogs, sale, dlc, dltt, pstk, mkvalt, fic, gvkey) %>%
        mutate(datadate = lubridate::ymd(datadate),
               year = lubridate::year(datadate),
               lpermco = as.integer(lpermco),
               sic2 = floor(as.integer(sic)/100)) 

    cstat_a = cstat_a %>%
        mutate(xrd = ifelse(xrd < 0, NA, xrd),
               sale = ifelse(sale < 0, NA, sale),
               cogs = ifelse(cogs < 0, NA, cogs))

    cstat_a = cstat_a %>%
        select(lpermco, year, everything()) %>%
        arrange(lpermco, year) %>%
        mutate(profit = sale - cogs,
               ev = dlc + dltt + pstk + mkvalt,
               glsc_include = (xrd > 0) & !is.na(xrd))

    
    cstat_a = cstat_a %>%
        distinct(lpermco, year, .keep_all = TRUE)

    cstat_a = leadlag_var(cstat_a, -1)
    cstat_a = leadlag_var(cstat_a, 1)

    cstat_a = cstat_a %>%
        filter(!is.na(datadate)) %>%
        filter(year >= 1960 & year <= 2019)


    cstat_a = cstat_a %>%
        mutate(sic = as.integer(sic)) %>%
        filter(!(sic >= 6000 & sic <= 6799)) %>%
        filter(!(sic >= 4900 & sic <= 4949)) %>%
        filter(!is.na(at) | at < 0) %>%
        filter(!is.na(sic)) %>%
        mutate(permco = lpermco)
    
    cstat_a = cstat_a %>%
        filter(!(floor(sic/100) %in% remove_sic$sic2)) %>%
        filter(!(sic %in% remove_sic$sic4)) 

    if (ficUSA) {
        cstat_a = cstat_a %>%
            filter(fic == "USA")
    }

    cstat_a = cstat_a %>%
        mutate(rdsales = xrd/sale)

    cstat_a = cstat_a %>%
        mutate(rdsales = ifelse(is.infinite(rdsales), NA, rdsales)) 
   
    cstat_a = cstat_a %>%
        arrange(lpermco, year) %>%
        select(lpermco, year, sic, sic2, rdsales, xrd, sale, sale_tm1, sale_t1,
               profit, profit_tm1, profit_t1, at, glsc_include, gvkey)

    cstat_a = cstat_a %>%
        mutate(salegrowth = sale_t1/sale - 1,
               profitgrowth = profit_t1/profit - 1) %>%
        mutate(profitgrowth = ifelse(profit < 0, NA, profitgrowth),
               salegrowth = ifelse(sale < 0, NA, salegrowth))

    cstat_a = cstat_a %>%               
            mutate(profitgrowth = ifelse(is.infinite(profitgrowth), NA, profitgrowth),
                   salegrowth = ifelse(is.infinite(salegrowth), NA, salegrowth))

    cstat_a = cstat_a %>%
            mutate(profitgrowth = ifelse(profitgrowth < winsor_drop[1] | profitgrowth > winsor_drop[2], NA, profitgrowth),
                   salegrowth = ifelse(salegrowth < winsor_drop[1] | profitgrowth > winsor_drop[2], NA, salegrowth),
                   rdsales = ifelse(rdsales < winsor_drop[1] | rdsales > winsor_drop[2], NA, rdsales))
    

    if (uniformsample) {
        if (posRD) {
            cstat_a = cstat_a %>%
                filter(xrd > 0) %>%
                filter(!is.na(rdsales))
        }

        cstat_a = cstat_a %>%
                filter(sale > 0 & profit > 0) %>% 
                filter(!is.na(profitgrowth) & !is.na(salegrowth)) 
                
        if (add_filters == "profit") { 
            cstat_a = cstat_a  %>% 
                filter(!is.na(profit) & !is.na(profit_tm1) & !is.na(profit_t1))

        } 

        if (nrow(profit_restrict) > 0) {
            cstat_a_pos = inner_join(cstat_a, profit_restrict, by = c("lpermco", "year")) %>%  
                filter(numpos)
        } else {
            cstat_a_pos = cstat_a 
        }

        cstat_a_noprof_restrict = cstat_a
        rm(cstat_a)
      
    } else {
        cstat_a_pos = cstat_a
        cstat_a_noprof_restrict = cstat_a
        rm(cstat_a) 
    }

    iter = list(cstat_a_pos = cstat_a_pos, cstat_a_noprof_restrict = cstat_a_noprof_restrict)

    bin_func = function(cstat_a_pos, add_filters, nbucket_cstat, nbucket_patent) {

        cstat_a_pos = bin_quintiles(cstat_a_pos, bucketvar = add_filters, nbucket = nbucket_cstat)
        cstat_a_pos_patent = bin_quintiles(cstat_a_pos, bucketvar = add_filters, nbucket = nbucket_patent) %>%
            dplyr::rename(bucket_patent = bucket, bucket_tm1_patent = bucket_tm1, bucket_t1_patent = bucket_t1) %>%
            select(lpermco, year, bucket_patent, bucket_tm1_patent, bucket_t1_patent)
        cstat_a_pos_halves = bin_quintiles(cstat_a_pos, bucketvar = add_filters, nbucket = 2) %>%
            dplyr::rename(bucket_2 = bucket, bucket_tm1_2 = bucket_tm1, bucket_t1_2 = bucket_t1) %>%
            select(lpermco, year, bucket_2, bucket_tm1_2, bucket_t1_2)

        cstat_a_pos = left_join(cstat_a_pos, cstat_a_pos_halves, by = c("lpermco", "year")) %>%
            left_join(cstat_a_pos_patent, by = c("lpermco", "year"))

        return(cstat_a_pos)

    }
    iter = map(.x = iter, 
               .f = bin_func, 
               add_filters = add_filters, nbucket_cstat = nbucket_cstat, nbucket_patent = nbucket_patent)

    cstat_a_pos = iter$cstat_a_pos
    cstat_a_noprof_restrict = iter$cstat_a_noprof_restrict
    rm(iter)

    cstat_a_pos = cstat_a_pos %>%
        group_by(year) %>%
        mutate(p90_t = quantile(.data[[add_filters]], probs = 0.9, na.rm = TRUE),
               p90_t1 = quantile(.data[[paste0(add_filters, "_t1")]], probs = 0.9, na.rm = TRUE),
               p80_t = quantile(.data[[add_filters]], probs = 0.8, na.rm = TRUE),
               p80_t1 = quantile(.data[[paste0(add_filters, "_t1")]], probs = 0.8, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(top10_t = .data[[add_filters]] > p90_t,
               top10_t1 = .data[[paste0(add_filters, "_t1")]] > p90_t1,
               top20_t = .data[[add_filters]] > p80_t,
               top20_t1 = .data[[paste0(add_filters, "_t1")]] > p80_t1) 

    cstat_a_pos = cstat_a_pos %>%
        group_by(year) %>%
        mutate(rdsales_w = Winsorize(rdsales, probs = winsor_vec, na.rm = TRUE),
               profitgrowth_w = Winsorize(profitgrowth, probs = winsor_vec, na.rm = TRUE)) %>%
        ungroup()

    rdsales_quintile_t = tabulate_by_quintile(cstat_a_pos, "rdsales_w", "median")
    profitvol_t = tabulate_by_quintile(cstat_a_pos, "profitgrowth_w", "sd")

    transmatout = salesmatfunc("bucket", "bucket_t1", df0 = cstat_a_pos, "quintile")
    transmatbottom90 = salesmatfunc("top10_t", "top10_t1", df0 = cstat_a_pos, "top10/bottom90")
    transmattop10top20 = salesmatfunc("top20_t", "top20_t1", df0 = cstat_a_pos, "top10/bottom90")
    transmathalves = salesmatfunc("bucket_2", "bucket_t1_2", df0 = cstat_a_pos, "halves")

    regdata = cstat_a_pos %>%
        mutate(bucket = as.factor(bucket),
               bucket_tm1 = as.factor(bucket_tm1),
               year = as.factor(year),
               sic = as.factor(sic),
               sic2 = as.factor(sic2), 
               lpermco = as.factor(lpermco))

    RHSiter = function(RHS) {
        if (RHS == "bucket") {
            RHS_tm1 = "bucket_tm1"
        } else if (RHS == "bucket_2") {
            RHS_tm1 = "bucket_tm1_2"
        }


        reg_profitgrowth_tm1 = felm(formula(paste0("profitgrowth_w ~ ", RHS_tm1, " | year:", fe, " | 0 | lpermco")), regdata,
                                    keepX = FALSE, keepCX = FALSE)
   
        reg_rdsales_tm1 = felm(formula(paste0("rdsales_w ~ ", RHS_tm1, " | year:", fe, " | 0 | lpermco")), regdata,
                               keepX = FALSE, keepCX = FALSE)
    

        reglist = list(reg_profitgrowth_tm1 = reg_profitgrowth_tm1,
                       reg_rdsales_tm1 = reg_rdsales_tm1)
    }
    regsiter = map(.x = c(nbucket = "bucket", two = "bucket_2"), .f = RHSiter)
    reglist = regsiter$nbucket
    reglist_2 = regsiter$two

    if (!uniformsample) {
        reglist = list()
    }

    return(list(cstat_a_pos = cstat_a_pos, 
                cstat_a_noprof_restrict = cstat_a_noprof_restrict, 
                regs = reglist,
                regs_2 = reglist_2, 
                rdsales_quintile_t = rdsales_quintile_t, 
                profitvol_t = profitvol_t,
                transmatout = transmatout,
                transmathalves = transmathalves,
                transmatoutbottom90 = transmatbottom90, 
                transmattop10top20 = transmattop10top20))
}
cstat_out_dropm1 = clean_cstat(cstat_a0, profit_restrict = restrict_qtr_cstat$three, uniformsample = TRUE, winsor_drop = c(-1, Inf), 
                        winsor_vec = c(0.01, 0.99), ficUSA = TRUE, nbucket_cstat = 5, nbucket_patent = 5, add_filters = "profit", 
                        remove_sic = list(sic2 = c(), sic4 = c()), fe = "sic2",
                        posRD = TRUE)

cstat_out_dropm1_noposRD = clean_cstat(cstat_a0, profit_restrict = restrict_qtr_cstat$three, uniformsample = TRUE, winsor_drop = c(-1, Inf), 
                        winsor_vec = c(0.01, 0.99), ficUSA = TRUE, nbucket_cstat = 5, nbucket_patent = 5, add_filters = "profit", 
                        remove_sic = list(sic2 = c(), sic4 = c()), fe = "sic2",
                        posRD = FALSE)

cstat_out_dropm1_nothreeqtr = clean_cstat(cstat_a0, profit_restrict = data.frame(), uniformsample = TRUE, winsor_drop = c(-1, Inf), 
                        winsor_vec = c(0.01, 0.99), ficUSA = TRUE, nbucket_cstat = 5, nbucket_patent = 5, add_filters = "profit", 
                        remove_sic = list(sic2 = c(), sic4 = c()), fe = "sic2",
                        posRD = TRUE)

cstat_out_full = clean_cstat(cstat_a0, profit_restrict = data.frame(), uniformsample = FALSE, winsor_drop = c(-Inf, Inf), 
                             winsor_vec = c(0.01, 0.99), ficUSA = FALSE, nbucket_cstat = 5, nbucket_patent = 5, add_filters = "profit", 
                             remove_sic = list(sic2 = c(), sic4 = c()), fe = "sic2", 
                             posRD = TRUE)

write_cstat = function(cstat_out, name) {

    cstat_write = list()
    cstat_write$regs = list()
    cstat_write$regs$reg_profitgrowth_tm1 = list()
    cstat_write$regs$reg_rdsales_tm1 = list()


    cstat_write$regs$reg_profitgrowth_tm1$coefficients = cstat_out$regs$reg_profitgrowth_tm1$coefficients
    cstat_write$regs$reg_profitgrowth_tm1$cse = cstat_out$regs$reg_profitgrowth_tm1$cse
    cstat_write$regs$reg_profitgrowth_tm1$N = cstat_out$regs$reg_profitgrowth_tm1$N

    cstat_write$regs$reg_rdsales_tm1$coefficients = cstat_out$regs$reg_rdsales_tm1$coefficients
    cstat_write$regs$reg_rdsales_tm1$cse = cstat_out$regs$reg_rdsales_tm1$cse
    cstat_write$regs$reg_rdsales_tm1$N = cstat_out$regs$reg_rdsales_tm1$N


    write_csv(cstat_out$rdsales_quintile_t, paste0("Output/Store_Data/rdsales_t_targets_", name, ".csv"))

    reg_targets = data.frame(dummies = 2:5,
                            profitgrowth_tm1 = coefficients(cstat_out$regs$reg_profitgrowth_tm1),
                            rdsales_tm1 = coefficients(cstat_out$regs$reg_rdsales_tm1))

    write_csv(cstat_out$profitvol_t, paste0("Output/Store_Data/profitvol_t_", name, ".csv"))
    
    return(list(cstat_write = list(cstat_write,  paste0("Output/Store_Data/cstat_out_", name, ".rds")),
                reg_targets = list(reg_targets, paste0("Output/Store_Data/reg_targets_", name, ".csv")),
                transmatout = list(as.data.frame(cstat_out$transmatout), paste0("Output/Store_Data/transmat_targets_", name, ".csv")),
                transmatoutbottom90 = list(as.data.frame(cstat_out$transmatoutbottom90), paste0("Output/Store_Data/transmat_targets_bottom90_", name, ".csv")),
                transmattop10top20 = list(as.data.frame(cstat_out$transmattop10top20),  paste0("Output/Store_Data/transmat_targets_top10top20_", name, ".csv")),
                transmathalves = list(as.data.frame(cstat_out$transmathalves), paste0("Output/Store_Data/transmat_targets_halves_", name, ".csv"))))


               
}
write_cstat_list = list(three_dropm1 = cstat_out_dropm1,
                        three_dropm1_noposRD = cstat_out_dropm1_noposRD,
                        three_dropm1_nothreeqtr = cstat_out_dropm1_nothreeqtr)

write_cstat_list = write_cstat_list %>%
    map2(.y = names(write_cstat_list), 
         .f = write_cstat)

write_cstat_apply = function(x) {
    if (str_detect(x[[2]], ".rds")) {
        saveRDS(x[[1]], x[[2]])
    } else if (str_detect(x[[2]], ".csv")) {
        write_csv(x[[1]], x[[2]])
    }

    return(invisible(NULL))
}

walk(.x = write_cstat_list[[1]], .f = write_cstat_apply)




nberpatent_munge = function(nber, nberlink, dynamiclink, cstat_a, merge_yr, 
                            yearcutoff = 20) {

    issueyear = nber$assgNBER %>%
        mutate(filed_year = as.integer(appyear),
               issue_year = as.integer(gyear),
               patent = as.integer(patent)) %>%
        select(wku = patent, issue_year) %>%
        distinct()

    citdata = nber$citeNBER %>%
        mutate_all(as.integer) %>%
        full_join(issueyear, by = c("cited" = "wku")) %>%
        dplyr::rename(cited_year = issue_year) %>%
        full_join(issueyear, by = c("citing" = "wku")) %>%
        dplyr::rename(citing_year = issue_year)


     citdata_tmp = citdata %>%
        filter(citing_year >= cited_year) %>%
        mutate(numyears = citing_year - cited_year)
               
    citsummary = citdata_tmp %>%
        group_by(cited) %>%
        summarise(fcit020 = length(unique(citing[(!is.na(citing)) & (numyears <= yearcutoff) ]))) %>%
        ungroup()

    tmp3 = nber$assgNBER %>%
        distinct(patent, .keep_all = TRUE) %>%
        filter(allcites == 0) %>%
        dplyr::select(cited = patent, fcit020 = allcites)

    citsummary = bind_rows(citsummary, tmp3)

    tmp_merge = nber$ipcNBER %>%
        distinct(patent, pdpass, icl_class, .keep_all = TRUE) %>%
        filter(!is.na(pdpass)) 

    ipcnber = inner_join(tmp_merge, citsummary, by = c("patent" = "cited"), relationship = "many-to-many") %>%
        select(wku = patent, issue_year = gyear, filed_year = appyear, CPC = icl_class, fcit020, everything()) %>%
        mutate(wku = as.integer(wku),
               issue_year = as.integer(issue_year),
               filed_year = as.integer(filed_year),
               CPC3 = substr(CPC, 1, 3)) 

    ipcnber = ipcnber %>%
        filter(tolower(substr(CPC,1,1)) %in% letters) %>%
        filter(!is.na(pdpass)) 
        
    ipcnber = ipcnber %>% 
        select(-CPC) %>%
        dplyr::rename(CPC = CPC3) %>%
        mutate(CPC = as.integer(as.factor(CPC))) %>% 
        distinct(wku, CPC, pdpass, .keep_all = TRUE)
    
    pdpdynass = dynamiclink %>% 
        select(c("pdpass", "pdpco1", "pdpco2", "pdpco3", "pdpco4", "pdpco5")) %>%
        pivot_longer(c("pdpco1", "pdpco2", "pdpco3", "pdpco4", "pdpco5"), names_to = "seq", values_to = "pdpco") %>%
        mutate(seq = as.integer(str_remove(seq, "pdpco")),
            pdpco = as.integer(pdpco),
            pdpass = as.integer(pdpass)) %>%
        filter(!is.na(pdpco))

    gvkeydynass = dynamiclink %>% 
        select(c("pdpass", "gvkey1", "gvkey2", "gvkey3", "gvkey4", "gvkey5")) %>%
        pivot_longer(c("gvkey1", "gvkey2", "gvkey3", "gvkey4", "gvkey5"), names_to = "seq", values_to = "gvkey") %>%
        mutate(seq = as.integer(str_remove(seq, "gvkey")),
            gvkey = as.character(gvkey),
            pdpass = as.integer(pdpass),
            gvkey = case_when(str_length(gvkey) == 1 ~ paste0("00000", gvkey),
                              str_length(gvkey) == 2 ~ paste0("0000", gvkey),
                              str_length(gvkey) == 3 ~ paste0("000", gvkey),
                              str_length(gvkey) == 4 ~ paste0("00", gvkey),
                              str_length(gvkey) == 5 ~ paste0("0", gvkey),
                              str_length(gvkey) == 6 ~ gvkey)) %>%
        filter(!is.na(gvkey))

    nberlink_merge = nberlink %>%
        mutate(gvkey = as.character(gvkey),
               gvkey = case_when(str_length(gvkey) == 1 ~ paste0("00000", gvkey),
                              str_length(gvkey) == 2 ~ paste0("0000", gvkey),
                              str_length(gvkey) == 3 ~ paste0("000", gvkey),
                              str_length(gvkey) == 4 ~ paste0("00", gvkey),
                              str_length(gvkey) == 5 ~ paste0("0", gvkey),
                              str_length(gvkey) == 6 ~ gvkey)) 

    begyrdynass = dynamiclink %>% 
        select(c("pdpass", "begyr1", "begyr2", "begyr3", "begyr4", "begyr5")) %>%
        pivot_longer(c("begyr1", "begyr2", "begyr3", "begyr4", "begyr5"), names_to = "seq", values_to = "begyr") %>%
        mutate(seq = as.integer(str_remove(seq, "begyr")),
            begyr = as.integer(begyr),
            pdpass = as.integer(pdpass)) %>%
        filter(!is.na(begyr))

    endyrdynass = dynamiclink %>% 
        select(c("pdpass", "endyr1", "endyr2", "endyr3", "endyr4", "endyr5")) %>%
        pivot_longer(c("endyr1", "endyr2", "endyr3", "endyr4", "endyr5"), names_to = "seq", values_to = "endyr") %>%
        mutate(seq = as.integer(str_remove(seq, "endyr")),
            endyr = as.integer(endyr),
            pdpass = as.integer(pdpass)) %>%
        filter(!is.na(endyr))

    dynass = full_join(pdpdynass, gvkeydynass, by = c("pdpass", "seq")) %>%
        full_join(begyrdynass, by = c("pdpass", "seq")) %>%
        full_join(endyrdynass, by = c("pdpass", "seq"))

    join1 = inner_join(dynass, ipcnber, by = "pdpass", relationship = "many-to-many") %>%
            arrange(gvkey, pdpass, wku, issue_year, CPC) %>%
            select(gvkey, pdpass, wku, issue_year, CPC, everything())

    join1 = join1 %>%
        filter(!!sym(merge_yr) >= begyr & !!sym(merge_yr) <= endyr) 


    join3 = join1 %>%
        full_join(cstat_a %>%
                    distinct(year, gvkey, .keep_all = TRUE) %>%
                    select(gvkey, year, lpermco) %>%
                    setNames(c("gvkey", merge_yr, "lpermco")), by = c("gvkey", merge_yr))

    join3 = join3 %>%
        left_join(distinct(nberlink_merge, gvkey, match), by = c("gvkey")) %>%
        arrange(gvkey, .data[[merge_yr]]) %>%
        mutate(matchflag = !is.na(match)) %>%
        filter((matchflag) | (!is.na(wku)))

    join3 = join3 %>%
        arrange(lpermco, !!sym(merge_yr), wku) %>%
        select(lpermco, !!sym(merge_yr), wku, everything())


    tmp = join3 %>%
        group_by(wku, lpermco, CPC) %>%
        mutate(N = n()) %>%
        ungroup() %>%
        arrange(desc(N), fcit020, lpermco, !!sym(merge_yr), wku, CPC) %>%
        select(N, fcit020, gvkey, lpermco, pdpass, !!sym(merge_yr), wku, CPC, everything())
   
    join3 = tmp %>%
        arrange(lpermco, CPC, wku) %>% 
        distinct(lpermco, CPC, wku, .keep_all = TRUE)

    dummy_cols = join3 %>%
        filter(!is.na(CPC)) %>%
        mutate(CPC = as.factor(CPC)) %>%
        model.matrix(~ CPC + 0, .) %>%
        as_tibble()

    join3 = bind_cols(filter(join3, !is.na(CPC)), dummy_cols) %>%
        mutate(permco = lpermco) %>%
        select(-CPC) %>%
        distinct(permco, wku, .keep_all = TRUE) 

    names(join3)[str_detect(names(join3), "CPC")] = paste0("myclass", str_remove(names(join3)[str_detect(names(join3), "CPC")], "CPC"))

    join3 = join3 %>%
        dplyr::rename(fcitALL = fcit020)
 
    return(join3)

}


nber_clean = nberpatent_munge(nber = nber, nberlink = nberlink, dynamiclink = dynamiclink, cstat_a = cstat_out_full$cstat_a_noprof_restrict, 
                              merge_yr = "issue_year", yearcutoff = 10) 


    
patent_clean = function(kpss0, cstat_use, year_type, 
                        fe, ds, prctile_inp, fcitcount = 10,
                        winsor_vec = c(NA, NA), 
                        halletal2001 = FALSE) { 

 

    kpss = kpss0 %>%
        filter(!is.na(permco))
    
    if (is.finite(fcitcount)) { 
        kpss = kpss %>% 
            filter((.data[[year_type]] >= 1960) & (.data[[year_type]] <= max(.data[[year_type]], na.rm = TRUE) - fcitcount + 1 - 5))
    } else {
        kpss = kpss %>%
            filter((.data[[year_type]] >= 1960) & (.data[[year_type]] <= max(.data[[year_type]], na.rm = TRUE) - 0 + 1 - 5))
    }

    kpss = kpss %>%
            group_by(wku, permco) %>%
            mutate(N = n()) %>%
            ungroup() %>%
            arrange(wku, desc(N)) %>%
            distinct(wku, .keep_all = TRUE) %>%
            select(-N)

    

    append_rows = function(df) {
        df = df %>%
            pivot_longer(contains("myclass"), names_to = "CPC", values_to = "logical") %>%
            mutate(logical = as.logical(logical)) %>%
            filter(logical) %>%
            mutate(CPC = str_remove(CPC, "myclass"),
                   CPC = as.integer(CPC)) %>%
            select(-logical)
    }
    wku_cpc = kpss %>%
        distinct(wku, .keep_all = TRUE) %>%
        dplyr::select(wku, contains("myclass"))

    wku_cpc_split = split(wku_cpc, rep(1:ceiling(nrow(wku_cpc)/10000), each = 10000, length.out = nrow(wku_cpc))) 
    wku_cpc = map(.x = wku_cpc_split, .f = append_rows)

    wku_cpc = rbindlist(wku_cpc) %>%
        as_tibble() %>%
        mutate(wku = as.integer(wku)) %>%
        arrange(wku, CPC)

    wku_cpc_distinct = inner_join(wku_cpc, kpss %>% select(-contains("myclass")), by = "wku")
    wku_cpc_distinct2 = distinct(wku_cpc_distinct, wku, CPC)


    offenders = wku_cpc_distinct %>%
        group_by(wku, CPC) %>%
        mutate(N = n()) %>%
        ungroup() %>%
        arrange(desc(N), wku, CPC, permco) %>%
        filter(N > 1)


    if (ds == "kpst") {
        kpss1 = kpss %>%
            select(-firm_id, -FIRMID, -assignee_ocr, -contains("myclass")) %>%
            distinct() %>%
            inner_join(wku_cpc, by = "wku") %>%
            group_by(wku, CPC) %>%
            mutate(N = n()) %>%
            ungroup()

        if (max(kpss1$N) > 1) {
            stop("Still duplicates.")
        }
    } else if (ds == "nber") {
        kpss1 = distinct(wku_cpc_distinct, wku, CPC, .keep_all = TRUE)
    }

    if (ds == "kpst") {

        if (fcitcount == 10) {
            kpss1= kpss1 %>%
                select(-fcitALL) %>%
                dplyr::rename(fcitALL = fcit010)
        } else if (fcitcount == 20) {
            kpss1= kpss1 %>%
                select(-fcitALL) %>%
                dplyr::rename(fcitALL = fcit020)
        } else if (fcitcount == 5) {
            kpss1 = kpss1 %>%
                select(-fcitALL) %>%
                dplyr::rename(fcitALL = fcit05)
        }
    } 
    if (halletal2001) {
        trunc = kpss1 %>%
                distinct(wku, .keep_all = TRUE) 
                
        if (!any(is.na(winsor_vec))) {

            trunc1 = trunc %>%
                group_by(.data[[year_type]]) %>%
                mutate(fcitALL = DescTools::Winsorize(fcitALL, na.rm = TRUE, probs = winsor_vec, type = 5)) %>%
                ungroup() 
        }
                
        trunc = trunc1 %>%
                group_by(.data[[year_type]]) %>%
                summarise(cit_year_mean = mean(fcitALL, na.rm = TRUE)) %>%
                ungroup() 

        kpss1 = kpss1 %>%
            select(-fcitALL) %>%
            left_join(select(distinct(trunc1, wku, .keep_all = TRUE), wku, fcitALL), by = "wku") %>%
            inner_join(trunc, by = year_type)
                
        normalization = trunc$cit_year_mean[trunc[[year_type]] == 1990]
            
        kpss1 = kpss1 %>%
                mutate(fcit05_adj = normalization * fcitALL / cit_year_mean)
    } else {
        kpss1 = kpss1 %>%
            mutate(fcit05_adj = fcitALL)

        normalization = NA
    }

    kpss1 = kpss1 %>%
        select(permco, filed_year, issue_year, wku, CPC, fcit05_adj) %>%
        arrange(permco, filed_year, CPC, fcit05_adj, wku) 
        
    kpss1 = kpss1 %>%
        mutate(CPC_save = CPC)
 
  
    kpss1_p90 = kpss1 %>%
        distinct(CPC_save, !!sym(year_type), wku, .keep_all = TRUE) %>%
        group_by(CPC_save, !!sym(year_type)) %>%
        summarise(p90 = quantile(fcit05_adj, probs = prctile_inp, na.rm = TRUE, type = 5)) %>%
        ungroup()

    quality1 = left_join(kpss1, kpss1_p90, by = c("CPC_save", year_type)) %>%
        mutate(CPC = as.integer(CPC),
               CPC_save = as.integer(CPC_save), 
               issue_year = as.integer(issue_year),
               filed_year = as.integer(filed_year),
               permco = as.integer(permco),
               major = fcit05_adj > p90)
               
    quality = quality1 %>%
        distinct(CPC_save, permco, !!sym(year_type), wku, .keep_all = TRUE) %>%
        filter(!is.na(fcit05_adj)) %>% 
        group_by(CPC, permco, y = !!sym(year_type)) %>%
        summarise(N = sum(!is.na(fcit05_adj)),
                  quality = sum(major)/N,
                  mean_indicator = mean(major)) %>%
        ungroup()

    alt_quality = quality1 %>%
        distinct(permco, !!sym(year_type), wku, .keep_all = TRUE) %>%
        filter(!is.na(fcit05_adj)) %>%
        group_by(permco, y = !!sym(year_type)) %>%
        summarise(num_patents = sum(!is.na(fcit05_adj)),
                  total_cit = sum(fcit05_adj),
                  avg_cit = total_cit/num_patents) %>%
        ungroup()


    cstat_a = cstat_use$cstat_a_noprof_restrict %>%
        mutate(year = as.integer(year),
               lpermco = as.integer(lpermco))

   
    nbucket_patent = max(cstat_use$cstat_a_noprof_restrict$bucket_patent, na.rm = TRUE)
    merge1 = inner_join(quality, cstat_a, by = c("permco" = "lpermco", "y" = "year")) 
   

    merge1 = merge1 %>%
        mutate(sic_FE = as.integer(sic2)) 

    merge1 = merge1 %>%
        mutate(year = y)
    
    drop_tm1 = "profit_tm1"
  
    merge1_tm1 = merge1 %>%
        filter(!is.na(!!sym(drop_tm1)))
    
    merge1_tm1 = merge1_tm1 %>%
        mutate(bucket = as.factor(bucket_patent),
               bucket_tm1 = as.factor(bucket_tm1_patent),
               y = as.factor(y),
               sic_FE = as.factor(sic_FE),
               permco=  as.factor(permco))
    
    merge1_tm1 = tryCatch(merge1_tm1 %>% mutate(CPC = as.factor(CPC)), error = function(e) merge1_tm1) 

    reg2profdummy_tm1 = felm(formula(paste0("quality ~ bucket_tm1 | ", fe, " | 0 | permco")), merge1_tm1)
    reg2profdummy_tm1_2 = felm(formula(paste0("quality ~ bucket_tm1_2 | ", fe, " | 0 | permco")), merge1_tm1)

    cdf_uncond = data.frame()

  
    
    return(list(kpss1 = kpss1, 
                normalization = normalization, 
                regs = list(reg2profdummy_tm1 = reg2profdummy_tm1,
                            reg2profdummy_tm1_2 = reg2profdummy_tm1_2),
                cdfs = list(cdf_uncond = cdf_uncond)))
}
kpss_nber = function(kpss0, ds) {
    patent_data_out_p90 = patent_clean(kpss0, cstat_use = cstat_out_dropm1, year_type = "issue_year", 
                                       fe = "y:CPC", 
                                       ds = ds, prctile_inp = 0.9, fcitcount = 10, 
                                       winsor_vec = c(NA, NA),
                                       halletal2001 = FALSE)

    patent_data_out_p80 = patent_clean(kpss0, cstat_use = cstat_out_dropm1, year_type = "issue_year", 
                                       fe = "y:CPC", 
                                       ds = ds, prctile_inp = 0.8, fcitcount = 10,
                                       winsor_vec = c(NA, NA), 
                                       halletal2001 = FALSE)

    patent_data_out_p75 = patent_clean(kpss0, cstat_use = cstat_out_dropm1, year_type = "issue_year", 
                                       fe = "y:CPC", 
                                       ds = ds, prctile_inp = 0.75, fcitcount = 10,
                                       winsor_vec = c(NA, NA), 
                                       halletal2001 = FALSE)

    return(list(patent_data_out_p90 = patent_data_out_p90, patent_data_out_p80 = patent_data_out_p80,
                patent_data_out_p75 = patent_data_out_p75))
}
patent_list = list(kpst = kpss0, nber = nber_clean)
patent_run = map2(.x = patent_list, .y = names(patent_list), .f = kpss_nber)

patent_data_out_p90 = patent_run$kpst$patent_data_out_p90
patent_data_out_p80 = patent_run$kpst$patent_data_out_p80

patent_data_out_p90_nber = patent_run$nber$patent_data_out_p90
patent_data_out_p80_nber = patent_run$nber$patent_data_out_p80



write_csv(patent_data_out_p90$cdfs$cdf_uncond, "Output/Store_Data/cdf_current.csv")
write_csv(patent_data_out_p90$kpss1, "Output/Store_Data/kpss1.csv")

write_csv(patent_data_out_p90_nber$kpss1, "Output/Store_Data/nber_cit.csv")


write_patent = function(patent_out, name, ds_name) {
      
        walk2(.x = patent_out,
            .y = names(patent_out), 
            .f = ~saveRDS(list(coefficients = map(.x = .x$regs, .f = ~.x$coefficients),
                                cse = map(.x = .x$regs, .f = ~.x$cse),
                                cdf = .x$cdfs,
                                N = map(.x = .x$regs, .f = ~.x$N)),
                            paste0("Output/Store_Data/patent_out_", name, "_", .y, "_", ds_name, ".rds")))

        return(invisible(NULL))
}
write_patent(list(patent_data_out_p90 = patent_data_out_p90), "patent_data_out_p90", ds_name = "kpst")
write_patent(list(patent_data_out_p80 = patent_data_out_p80), "patent_data_out_p80", ds_name = "kpst")
write_patent(list(patent_data_out_p90_nber = patent_data_out_p90_nber), "patent_data_out_p90", ds_name = "nber")
write_patent(list(patent_data_out_p80_nber = patent_data_out_p80_nber), "patent_data_out_p80", ds_name = "nber")

kogan_targets = function(kogan0, cstat_a0, winsor_vec, glsc_include) {

    kogan = kogan0 %>%
        group_by(permco, y) %>%
        summarise(count_patent = sum(!is.na(qje_nominal_xi)),
                  qje_nominal_xi = sum(qje_nominal_xi),
                  jor_measure_t = qje_nominal_xi/count_patent,
                  jor_measure_t_nodiv = qje_nominal_xi,
                  nominal_xi = sum(nominal_xi)) %>%
        ungroup()
    
    cstat_a = cstat_a0 %>%
        filter(year <= max(kogan$y))

    mergedf = left_join(cstat_a, kogan, by = c("year" = "y", "lpermco" = "permco")) %>%
        dplyr::rename(permco = lpermco) 

    dups = distinct(mergedf, permco, year) 
    if (nrow(dups) != nrow(mergedf)) {
        stop("duplicate obs.")
    } 

    mergedf = mergedf %>%
        arrange(permco, year) %>%
        mutate(qje_nominal_xi_NA = qje_nominal_xi,
               qje_nominal_xi = ifelse(is.na(qje_nominal_xi), 0, qje_nominal_xi), 
               theta = qje_nominal_xi/at)
               
    mergedf = mergedf %>%
        group_by(year) %>%
        mutate(theta_w = Winsorize(theta, probs= winsor_vec, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(kogan_years = (year >= 1950 & year <= 2010 & !is.na(year)))

    if (glsc_include) {
        mergedf = mergedf %>%
            filter(glsc_include)
    }

    
    check_uncond = mergedf %>%
        summarise(mean = 100*mean(theta_w, na.rm = TRUE),
                  sd = 100*sd(theta_w, na.rm = TRUE),
                  p1 = 100*quantile(theta_w, na.rm = TRUE, probs = 0.01),
                  p5 = 100*quantile(theta_w, na.rm = TRUE, probs = 0.05),
                  p10 = 100*quantile(theta_w, na.rm = TRUE, probs = 0.1),
                  p25 = 100*quantile(theta_w, na.rm = TRUE, probs = 0.25),
                  p50 = 100*quantile(theta_w, na.rm = TRUE, probs = 0.5),
                  p75 = 100*quantile(theta_w, na.rm = TRUE, probs = 0.75),
                  p90 = 100*quantile(theta_w, na.rm = TRUE, probs = 0.9),
                  p95 = 100*quantile(theta_w, na.rm = TRUE, probs = 0.95), 
                  p99 = 100*quantile(theta_w, na.rm = TRUE, probs = 0.99)) %>%
        ungroup()

    return(list(check_uncond = check_uncond))

}
kogan_targets_out = kogan_targets(kogan0, cstat_a0 = cstat_out_full$cstat_a_pos, winsor_vec = c(0.01, 0.99), 
                                  glsc_include = TRUE)
kogan_targets_out_FALSE = kogan_targets(kogan0, cstat_a0 = cstat_out_full$cstat_a_pos, winsor_vec = c(0.01, 0.99), 
                                        glsc_include = FALSE)
write_csv(kogan_targets_out$check_uncond, "Output/Store_Data/kogan_targets.csv")
write_csv(kogan_targets_out_FALSE$check_uncond, "Output/Store_Data/kogan_targets_glscincludeFALSE.csv")

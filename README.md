# Financial discounts and creative-destruction dynamics

This is the replication archive for Chikis, Goldberg, and Lopez-Salido (2023). 

*Abstract*: We develop a growth model of strategic firm innovation and financial discounts.  A central focus is that firms' discount rate exceeds the risk-free rate.  Constraints on creative destruction, firm financing, and knowledge diffusion determine, by disadvantaging market laggards relative to leaders, whether a high required excess return fosters or restrains productivity growth and market competition.  The estimated model matches many cross-sectional patterns capturing R&D, patent quality, profits, and reallocation.  Variation in the excess return implies low safe interest rates are associated with weak productivity growth, competition, and entry, consistent with macroeconomic trends.  Our results have novel implications for monetary policy.  


## Getting started
To run this code, you will need MATLAB and R. To run the empirical portion of the code, you will need a subscription to WRDS. Note that you can run the model portion of the replication code without bothering with the empirical portion. 


### Software 
* MATLAB
* R

### External data
Following the below instructions, you will need to move all datasets into `Data/`

* Download Compustat Annual from WRDS with the following options selected
    * Get Data
    * CRSP
    * CRSP/Compustat Merged
    * Fundamentals Annual
        * Select full time series
        * gvkey 
        * Search the entire database
        * Uncheck CAD
        * Select: LC, LU
        * Select: datadate, lpermco, sic, at, xrd, cogs, sale, dlc, dltt, pstk, mkvalt, fic, gvkey
        * Export as CSV and save as `Data/cstat_a.csv`
    * Fundamentals Quarterly
        * Select full time series
        * gvkey 
        * Search the entire database
        * Uncheck CAD
        * Select: LC, LU
        * Select: datadate, lpermco, saleq, saley, cogsq, cogsy
        * Export as dta and save as `Data/cstat_q.dta`
* Download patent data from [Kogan et al. (2021)](https://github.com/KPSS2017/Measuring-Technological-Innovation-Over-the-Long-Run-Replication-Kit)
    * Download their repository and follow their instructions
    * You will need to move the file `MergedPatentData.dta` to `Data/MergedPatentData.dta` 
* Download the patent valuation data from [Kogan et al. (2017)](https://github.com/KPSS2017/Technological-Innovation-Resource-Allocation-and-Growth-Extended-Data)
    * Navigate to branch `update_2020`
    * Navigate to commit history
    * Download the the files `KPSS_2019_Public.csv` and `Match_patent_permco_permno_2022.csv.zip` from the commit from Jun 25, 2020 and save them in `Data/`
* Download the NBER patent data from [NBER](https://sites.google.com/site/patentdataproject/Home/downloads?authuser=0)
    * `cite76_06.dta`
    * `pat76_06_assg.dta`
    * `pat76_06_ipc.dta`
    * `pdpcohdr.dta`
    * `pat76_06_ipc.dta`
    * `dynass.dta`
* Download the BDS data from [Census](https://www.census.gov/data/datasets/time-series/econ/bds/bds-datasets.html)
    * Economy-wide (save as `Data/bds2020 (2).csv`) 
    * Firm age (save as `Data/bds2020_fa (3).csv`)
* Download [Michigan consumer survey data](https://data.sca.isr.umich.edu/)
    * Table 32 (saved as `Data/michigan_1_inflation_monthly.csv`)
    * Table 33 (saved as `Data/Data/michigan_510_inflation_monthly.csv`)
* Download data from [Paul Gomme's website](https://paulgomme.github.io/#data)
    * "comma-separated values (CSV) file."
    * Save as `Data/gomme.csv`
    * Code that generates a time-series of intangible-adjusted rates of return on capital is modified from code provided on Gomme's website and is available upon request from the authors


### Executing program

We provide multiple options for running the code. If you are working in a PBS/Slurm environment, and you have downloaded the external data specified above, you can run the full suite of programs by navigating to wherever you saved the repository and submitting to the terminal something akin to:
```bash
qsub Code/Main/run_main_cgls.sh
```

Note you should first modify the working directory part of the shell script
```bash
cd "where/you/saved/CGLS_RepCode_Final/"
workingdirectory=where/you/saved/CGLS_RepCode_Final/
```
in addition to the output for the logs of your run
```bash
#PBS -o where/you/saved/CGLS_RepCode_Final/Logs/output.txt
```

Alternatively, you can open up MATLAB; navigate to the directory you stored the archive in and run `Code/Main/main_cgls.m`. Again, you will want to change the working directory:
```matlab
% Fill this in per your local file setup
cd("where/you/saved/CGLS_RepCode_Final/")
```

You can also choose to run without the external data via a job scheduler
```bash
qsub Code/Main/run_main_cgls_no_data.sh
```
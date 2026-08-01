# GrayWhaleCalfProduction
This repository contains R code to extract and analyze gray whale calf production data that are collected at Piedras Blancas, CA, annually. 
Necessary functions are stored in GrayWhaleCalfProduction_fcns.R. Some of the functions in the file may be obsolete and unused. The analytical
part is conducted via JAGS (Just Another Gibbs Sampler, which can be downloaded from https://sourceforge.net/projects/mcmc-jags/files/. In 2026, 
an updated analytical approach was introduced because of the failure of the previous model when the counts for the 2026 season declined and they 
were spread over several weeks. The new approach uses the Negative Binomial model instead of the Poisson and is coded in Stan (https://mc-stan.org/). 

It is assumed that this repository is downloaded to your local computer via RStudio and a project is created. Please create three empty folders
within the project folder and name them "data," "figures," and "RData." 

In 2023, new data extraction code was developed due to inconsistencies found in the previous version (Collating and Formatting Passdown.R).
The comparison of the two extraction methods is provided in a notebook file named "compare_data_extraction.Rmd." In the past report, I did not 
use the new extraction method. I appended the new estimate for the 2023 season to the 2022 report. For the 2023 report, I did use the new 
extraction but the table was replaced with the old version. These reports can be found in Report_calf_production_2022.Rmd and Report_calf_production_2023.Rmd.

The order of the data processing and analysis should be: (THIS NEEDS TO BE UPDATED 2026-08-01)
1. Raw data should be saved in the "data/All data" directory. Folder names are C_C YYYY, where YYYY indicates a four-digit year. 
2. Edit list.sheet.names.inshore, list.col.types.inshore, and list.col.names.inshore in GrayWhaleCalfProduction_fcns.R (if you pull the most recent version, this should be done) 
3. Run Extract Excel Data.Rmd
4. Run calf_production_jags.R
6. Edit Report_calf_production_2023.Rmd to create a new report.

Previous approach:
In Stewart and Weller (2021), the mean ($\lambda$) of the number of whales passing by the survey area is assumed to be the same within a week. 
There is no assumption about the mean except it is bounded between 0 and 40 (i.e., $ \lambda \sim UNIF(0,40)$). The weekly mean assumption is a bit arbitrary. It was also assumed that the true number of whales in the survey area per 3-hr period was an independent Poisson random deviate (i.e., $n_{true} \sim POI(\lambda)$). The observed number of whales per shift is assumed to be a binomial random deviate (i.e., $n_{obs} \sim BIN(n_{true}, p_{obs})$), where $p_{obs}$ was estimated through a calibration study ($p_{obs} \sim N(0.889, 0.06375^2)$). 

Updated approach:


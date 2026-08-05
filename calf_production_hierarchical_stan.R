# This R script file conducts all hierarchical analyses in the tech memo using Stan.
# Results can be saved when the save.files switch is TRUE. 
# The data extraction is best with "v3", which is the most up-to-date data extraction
# version.
# Analysis model can be selected by choosing either one of "nb_fixedyear", 
# "nb_shiftingpeak", or "nb_shiftingpeak_centered". The "nb_fixedyear" model is
# selected to produce final estimates for 2026. To fit the Negative-Binomial model
# to each year independently, use 'calf_production_single_year_stan.R'. 

rm(list = ls())
# Load required library
library(tidyverse)
library(cmdstanr)
library(loo)
library(ggplot2)
library(bayesplot)

source("GrayWhaleCalfProduction_fcns_v2.R")

save.files <- FALSE #TRUE #

data.ext <- "v3" # "v2"  # or v2
#model <- "nb"
#model <- "nb_sens_sigma_year"
model <- "nb_fixedyear"
#model <- "nb_shiftingpeak"
#model <- "nb_shiftingpeak_centered"

# Change this to use results from a previous run
run.date <- "2026-07-09" #Sys.Date() #"2026-06-26"    

data.path <- paste0("data//Formatted Annual Data Combined ", 
                    data.ext, "//")
FILES <- list.files(path = data.path,
                    pattern = paste0(data.ext, ".csv"))

# get data
all.data <- count.obs <- effort <- week <- n.obs <- n.weeks <- list()
years <- vector(mode = "numeric", length = length(FILES))

for(i in 1:length(FILES)){
  years[i] <- as.numeric(str_split(FILES[i], " Formatted_combined_v3.csv")[[1]][1]) 
  
    data <- read.csv(paste0(data.path, FILES[i]))
    data$Effort[is.na(data$Effort)] <- 0
    data$Effort[data$Effort > 3] <- 3
    data$Sightings[data$Effort == 0] <- 0  # no effort, no sightings
    
    all.data[[i]] <- data.frame(count_obs = as.numeric(data$Sightings),
                                effort_hours = data$Effort,
                                week_num = data$Week,
                                calendar_year = year(as.Date(data$Date)),
                                iso_wk = isoweek(as.Date(data$Date)),
                                Month = month(as.Date(data$Date)),
                                Day = day(as.Date(data$Date)),
                                wday = wday(as.Date(data$Date)),
                                Date = as.Date(data$Date))
    
}

calf_data <- do.call("rbind", all.data) 
calf_data %>%
  mutate(
    year_idx   = as.integer(as.factor(calendar_year)),
    
    # This aligns every year to the exact same biological timeline:
    week_idx   = as.integer(iso_wk - min(iso_wk) + 1),
    
    log_offset = ifelse(effort_hours > 0, log(effort_hours / 3.0), 0)
  ) %>%
  arrange(year_idx, week_idx) -> prepped_data_clean

all.years <- unique(prepped_data_clean$calendar_year)
model.file <- paste0("models//GWCalfCount_", model, ".stan")
out.file <- paste0("GWCalfCount_Stan_", model, "_", run.date)
stan_data <- list(
  n_obs      = nrow(prepped_data_clean),
  n_years    = max(prepped_data_clean$year_idx),
  n_weeks    = max(prepped_data_clean$week_idx),
  count_obs  = prepped_data_clean$count_obs,
  effort     = prepped_data_clean$effort_hours, # Passing raw effort to Stan
  log_offset = prepped_data_clean$log_offset, 
  year_idx   = prepped_data_clean$year_idx,
  week_idx   = prepped_data_clean$week_idx
)

# This only works when running smooth-centered peak model
if (model == "nb_shiftingpeak_centered"){
  init_fn <- function() list(
    p_obs = 0.889, 
    year_eff = rep(0, stan_data$n_years),
    amp = 2, width = 2, 
    mu_peak = (stan_data$n_weeks + 1)/2,
    sigma_peak = 0.5, 
    peak_raw = rep(0, stan_data$n_years), phi = 2
  )
  
} else {
  init_fn <- function() list(p_obs = 0.889,
                             phi = 2)
}


if (!file.exists(paste0("RData//", out.file, ".rds"))){
  # 1. Compile and Fit the Negative Binomial Model
  # In your R script where you load your model, add force_recompile = TRUE just once:
  mod_ <- cmdstan_model(stan_file = model.file,
                        #force_recompile = TRUE,
                        cpp_options = list(stan_threads = TRUE, 
                                           O = 3))
  tic <- Sys.time()
  fit_ <- mod_$sample(
    data            = stan_data,
    seed            = 123,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 1000,
    threads_per_chain = 2,
    iter_sampling   = 2000,
    #init            = 0.1,
    adapt_delta     = 0.99,
    max_treedepth   = 12,
    init            = init_fn,
    refresh         = 500
  )
  toc <- Sys.time() -  tic
  
  fit_$save_object(file = paste0("RData//", out.file, ".rds"))
  
  out.list <- list(out.filename = paste0("RData//", out.file, ".rds"),
                   data.ext = data.ext, # "v2"  # or v2
                   data.path = data.path,
                   model = model.file,
                   pre.stan.data = prepped_data_clean,
                   stan.data = stan_data,
                   Run.time = toc,
                   Run.Date = Sys.Date(),
                   System = Sys.getenv())
  
  saveRDS(out.list, 
          file = paste0("RData//", out.file, ".info"))
} else {
  
  fit_ <- readRDS(paste0("RData//", out.file, ".rds"))
  out.list <- readRDS(paste0("RData//", out.file, ".info"))
}

# 1. Extract the raw log-likelihood matrix [8000 iterations x 17008 columns]
raw_log_lik <- fit_$draws("log_lik", format = "matrix")

# 2. Identify the columns that actually represent real survey days
real_data_cols <- which(stan_data$effort > 0)

# 3. Subset the matrix to ONLY include those 6,120 real days
clean_log_lik_mat <- raw_log_lik[, real_data_cols]

# 4. Run LOO on the cleaned matrix
loo_nb <- loo(clean_log_lik_mat, cores = 4)
#print(loo_nb_fixed)

# Summarize global parameters and hierarchical standard deviations
if (model == "nb" | model == "nb_fixedyear"){
  global_summary <- fit_$summary(c("p_obs", "sigma_week", "year_eff",
                                   "phi",  "week_eff"),
                                 posterior::default_summary_measures(), 
                                 posterior::default_convergence_measures(),
                                 extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)))
  
} else if (model == "nb_shiftingpeak" | model == "nb_shiftingpeak_centered"){
  global_summary <- fit_$summary(c("p_obs", "year_eff",
                                   "phi",  "amp", "width",
                                   "mu_peak", "sigma_peak", "peak_week"),
                                 posterior::default_summary_measures(), 
                                 posterior::default_convergence_measures(),
                                 extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)))
  
}

stan.out <- stan.post.process(stan.fit = fit_, 
                              pre.stan.data = prepped_data_clean,
                              stan.data = stan_data, 
                              out.file.name = out.file,
                              save.file = save.files, 
                              save.fig = F)
  
PPC.out <- PPC_counts(stan.fit = fit_,
                      stan.data = stan_data)

diag.summary <- fit_$diagnostic_summary()

all.estimates <- stan.out$abundance

loo.out <- list(raw.log.lik = raw_log_lik,
                real.data.cols = real_data_cols,
                clean.log.lik.mat = clean_log_lik_mat,
                loo = loo_nb)

if (save.files){
  write.csv(all.estimates,
            file = paste0("data//all_estimates_", model, ".csv"))
  write.csv(global_summary,
            file = paste0("data//global_summary_", model, ".csv"))
  saveRDS(loo.out,
          file = paste0("RData//loo_out_", model, ".rds"))
  saveRDS(stan.out,
          file = paste0("RData//stan_out_", model, ".rds"))
  saveRDS(PPC.out,
          file = paste0("RData//PPC_out_", model, ".rds"))
}


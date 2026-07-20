# Converting the year-specific data to a hierarchical model that contains shared
# parameters over years:

rm(list = ls())
# Load required library
library(tidyverse)
library(cmdstanr)
library(loo)
library(ggplot2)
library(bayesplot)

source("GrayWhaleCalfProduction_fcns_v2.R")

data.ext <- "v3" # "v2"  # or v2
#model <- "nb"
model <- "nb_singleyear"
run.date <- Sys.Date() # "2026-06-30"    #  Change this to use results from a previous run

if (data.ext == "v1"){
  data.path <- "data/Formatted Annual Data/"
  FILES <- list.files(path = data.path, 
                      pattern = "Formatted.csv")
  
} else {
  data.path <- paste0("data/Formatted Annual Data Combined ", data.ext, "/")
  FILES <- list.files(path = data.path,
                      pattern = paste0(data.ext, ".csv"))
  
}

# get data
all.data <- count.obs <- effort <- week <- n.obs <- n.weeks <- list()
years <- vector(mode = "numeric", length = length(FILES))

for(i in 1:length(FILES)){
  if (data.ext == "v1"){
    years[i] <- as.numeric(str_split(FILES[i], " Formatted.csv")[[1]][1])
  } else {
    years[i] <- as.numeric(str_split(FILES[i], " Formatted_combined_v3.csv")[[1]][1]) 
  }
  
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
  arrange(year_idx, week_num) -> prepped_data_clean

    
all.years <- unique(prepped_data_clean$calendar_year)
model.file <- paste0("models/GWCalfCount_", model, ".stan")
stan.out <- PPC.out <- global_summary <- diag.summary <- list()
for (y in 1:max(prepped_data_clean$year_idx)){
  out.file <- paste0("GWCalfCount_Stan_", model, "_", 
                     all.years[y], "_", run.date)
  data.1year <- prepped_data_clean %>%
    filter(year_idx == y)
  stan_data <- list(
    n_obs      = nrow(data.1year),
    n_years    = max(data.1year$year_idx),
    n_weeks    = max(data.1year$week_idx),
    count_obs  = data.1year$count_obs,
    effort     = data.1year$effort_hours, # Passing raw effort to Stan
    log_offset = data.1year$log_offset, 
    year_idx   = data.1year$year_idx,
    week_idx   = data.1year$week_num
  )
  if (!file.exists(paste0("Rdata/", out.file, ".rds"))){
    # 1. Compile and Fit the Negative Binomial Model
    # In your R script where you load your model, add force_recompile = TRUE just once:
    mod_ <- cmdstan_model(stan_file = model.file,
                          #force_recompile = TRUE,   # This is useful when moving to a new R version
                          cpp_options = list(stan_threads = TRUE, 
                                             O = 3))
    tic <- Sys.time()
    fit_ <- mod_$sample(
      data            = stan_data,
      seed            = 12,
      chains          = 4,
      parallel_chains = 4,
      iter_warmup     = 1000,,
      threads_per_chain = 2,
      iter_sampling   = 2000,
      init            = 0.1,
      adapt_delta     = 0.99,
      refresh         = 500
    )
    toc <- Sys.time() -  tic
    
    fit_$save_object(file = paste0("Rdata/", out.file, ".rds"))
    
    out.list <- list(out.filename = paste0("Rdata/", out.file, ".rds"),
                     data.ext = data.ext, # "v2"  # or v2
                     data.path = data.path,
                     model = model.file,
                     pre.stan.data = prepped_data_clean,
                     stan.data = stan_data,
                     Run.time = toc,
                     Run.Date = Sys.Date(),
                     System = Sys.getenv())
    
    saveRDS(out.list, 
            file = paste0("Rdata/", out.file, ".info"))
  } else {
    
    fit_ <- readRDS(paste0("Rdata/", out.file, ".rds"))
    out.list <- readRDS(paste0("Rdata/", out.file, ".info"))
  }  
  
  # 1. Extract the raw log-likelihood matrix [8000 iterations x 17008 columns]
  raw_log_lik <- fit_$draws("log_lik", format = "matrix")
  
  # 2. Identify the columns that actually represent real survey days
  real_data_cols <- which(stan_data$effort > 0)
  
  # 3. Subset the matrix to ONLY include those 6,120 real days
  clean_log_lik_mat <- raw_log_lik[, real_data_cols]
  
  # 4. Run LOO on the cleaned matrix
  loo_nb_fixed <- loo(clean_log_lik_mat, cores = 4)
  #print(loo_nb_fixed)
  
  # Summarize global parameters and hierarchical standard deviations
  global_summary[[y]] <- fit_$summary(c("p_obs", "sigma_week", 
                                        "phi", "beta0", "week_eff"))
  
  
  #print(global_summary)
  stan.out[[y]] <- stan.post.process.1year(
    stan.fit = fit_, 
    pre.stan.data = data.1year,
    stan.data = stan_data, 
    out.file.name = out.file,
    save.file = F)
  
  PPC.out[[y]] <- PPC_counts(stan.fit = fit_,
                             stan.data = stan_data)
  
  diag.summary[[y]] <- fit_$diagnostic_summary()
}

all.estimates.list <- list() 
for (y in 1:length(stan.out)){
  all.estimates.list[[y]] <- stan.out[[y]]$abundance
}

all.estimates <- do.call("rbind", all.estimates.list)

rownames(all.estimates) <- all.years
all.estimates %>%
  rownames_to_column("Year") -> all.estimates 

write.csv(all.estimates,
          file = paste0("Data/all_estimates_", model, ".csv"))

nrows <- lapply(global_summary, FUN = nrow) %>%
  unlist() 

global.summary.df <- do.call("rbind", global_summary)
global.summary.df %>%
  mutate(Year = rep(all.years, times = nrows)) -> global.summary.df
write.csv(global.summary.df,
          file = "data/global_summary_single_year.csv")

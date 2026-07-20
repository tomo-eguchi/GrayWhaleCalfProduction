
# This script is a modified version of "Running Individual Year Model.R"
# I inherited it from Josh Stewart in early 2022.

# In this version, all models can be run by specifying the model name (v1 - vy) at the beginning.

#V7 IS INCOMPLETE. NEEDS TO BE FIXED. 2023-05-19

rm(list=ls())
library(jagsUI)
library(tidyverse)
library(lubridate)
library(bayesplot)

source("GrayWhaleCalfProduction_fcns.R")

save.file <- F
data.ext <- "v3" # "v2"  # or v2
model <- "v1" # "v1b" #"v8" #  "v1a"

#FILES <- list.files(pattern = ".csv$")
#
if (data.ext == "v1"){
  data.path <- "data/Formatted Annual Data/"
  FILES <- list.files(path = data.path, 
                      pattern = "Formatted.csv")
  
} else {
  data.path <- paste0("data/Formatted Annual Data ", data.ext, "/")
  FILES <- list.files(path = data.path,
                      pattern = paste0("Formatted_inshore_", data.ext, ".csv"))
  
}


MCMC.params <- list(n.samples = 100000,
                    n.thin = 100,
                    n.burnin = 50000,
                    n.chains = 5)

n.samples <- MCMC.params$n.chains * ((MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin)

jags.params <- c("count.true",
                 "lambda", "psi", "alpha",
                 "beta1", "beta2", "eps",
                 "p.obs.corr",
                 "p.obs",
                 "Total.Calves",
                 "loglik")

# get data
count.obs <- effort <- week <- n.obs <- n.weeks <- list()
years <- vector(mode = "numeric", length = length(FILES))
jm.out <- list()
for(i in 1:length(FILES)){
  if (data.ext == "v1"){
    years[i] <- as.numeric(str_split(FILES[i], " Formatted.csv")[[1]][1])
  } else {
    years[i] <- as.numeric(str_split(FILES[i], " Formatted_inshore")[[1]][1])
    
  }
  
  out.file.name <- paste0("RData/calf_estimates_", data.ext,
                          "_M", model, "_", years[i], ".rds") 
  
  if (file.exists(out.file.name)){
    jm.out[[i]] <- readRDS(out.file.name)
    
  } else {
    
    data <- read.csv(paste0(data.path, FILES[i]))
    data$Effort[is.na(data$Effort)] <- 0
    data$Effort[data$Effort > 3] <- 3
    data$Sightings[data$Effort == 0] <- 0  # no effort, no sightings
    
    jags.data <- list(count.obs = data$Sightings,
                      effort = data$Effort,
                      week = data$Week,
                      n.obs = length(data$Sightings),
                      n.weeks = max(data$Week),
                      weekly.max = data %>% 
                        group_by(Week) %>% 
                        summarize(weekly.max = max(Sightings)) %>% 
                        select (weekly.max) %>% 
                        as.vector() %>%
                        unlist() %>% 
                        unname())
    

    if (model == "v1b"){
      inits <- function() {list(z = rep(0, jags.data$n.obs) )}
    } else {
      inits <- NULL
    }
    
    
    jm <- jags(jags.data,
               inits = inits,
               parameters.to.save= jags.params,
               paste0("models/GWCalfCount_", model, ".jags"), 
               n.chains = MCMC.params$n.chains,
               n.burnin = MCMC.params$n.burnin,
               n.thin = MCMC.params$n.thin,
               n.iter = MCMC.params$n.samples,
               DIC = T, parallel=T)
    
    # This function is in GrayWhaleCalfProduction_fcns.R
    # params.to.monitor is a string of regular expression.
    # e.g., "^BF\\.Fixed|^K\\["
    params <- "^count\\.true\\[|^lambda\\[|^p\\.obs\\.corr\\[|^p\\.obs|^Total\\.Calves\\["
    jm.MCMC <- MCMC.diag(jm = jm, 
                         MCMC.params = MCMC.params,
                         params.to.monitor = params)
    
    jm.out[[i]] <- list(jm = jm,
                        MCMC.diag = jm.MCMC,
                        jags.data = jags.data,
                        MCMC.params = MCMC.params,
                        run.date = Sys.Date())     
    
    saveRDS(jm.out[[i]],
            file = out.file.name)
  }
  
}

posterior.summary <- list()
for (k in 1:length(jm.out)){
  posterior.summary[[k]] <- posterior::summarize_draws(jm.out[[k]]$jm$samples) %>%
    mutate(Year = years[k])
}

posterior.summary.df <- do.call(rbind, posterior.summary) 

posterior.summary.df %>%
  group_by(Year) %>%
  summarise(Max.Rhat = max(rhat, na.rm = T),
            Min.ESS.bulk = min(ess_bulk, na.rm = T),
            Min.ESS.tail = min(ess_tail, na.rm = T)) %>%
  mutate(n.obs = lapply(jm.out, FUN = function(x) x$jags.data$n.obs) %>% unlist(),
         effort = lapply(jm.out, FUN = function(x) x$jags.data$effort %>% sum()) %>% unlist(),
         n.whales = lapply(jm.out, FUN = function(x) x$jags.data$count.obs %>% sum()) %>% unlist()) -> Rhat.ESS.summary

saveRDS(Rhat.ESS.summary,
        file = "RData/Rhat_ESS_Summary_Jags2026.rds")

stats.total.calves <- lapply(jm.out, 
                             FUN = function(x){
  Mean <- x$jm$mean$Total.Calves
  Median <- x$jm$q50$Total.Calves
  LCL <- x$jm$q2.5$Total.Calves
  UCL <- x$jm$q97.5$Total.Calves
  
  return(data.frame(Mean = Mean,
                    Median = Median,
                    LCL = LCL,
                    UCL = UCL))
})

Estimates <- do.call(rbind, stats.total.calves)
Estimates$Year <- years
Estimates$Method <- model
#Estimates$Sys_env <- Sys.getenv()

if (save.file)
  write.csv(Estimates,
            paste0("data/Calf_Estimates_", data.ext, "_M", model, "_", Sys.Date(), ".csv"),
            row.names = F)


max.Rhat.list <- lapply(jm.out, FUN = function(x) x$MCMC.diag$max.Rhat)
max.Rhat <- do.call(rbind, max.Rhat.list) %>%
  as.data.frame() %>%
  mutate(Year = years)

# look at just 2025/2026
dat.2025 <- data.frame(Year = 2025,
                       week = jm.out[[31]]$jags.data$week,
                       n.obs = jm.out[[31]]$jags.data$n.obs,
                       count = jm.out[[31]]$jags.data$count.obs,
                       effort = jm.out[[31]]$jags.data$effort)

dat.2026 <- data.frame(Year = 2026,
                       week = jm.out[[32]]$jags.data$week,
                       n.obs = jm.out[[32]]$jags.data$n.obs,
                       count = jm.out[[32]]$jags.data$count.obs,
                       effort = jm.out[[32]]$jags.data$effort)

data.2025.2026 <- rbind(dat.2025, dat.2026)

data.2025.2026 %>%
  group_by(Year, week) %>%
  summarise(n.obs = sum(n.obs),
            count = sum(count),
            effort = sum(effort),
            count.per.effort = count/effort) -> summary.2025.2026

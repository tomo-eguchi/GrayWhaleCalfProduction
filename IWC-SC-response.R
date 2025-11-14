# A script to create plots for a request by Aimee Lang for the IWC Scientific Committee's
# requests on the calf production analysis.
# 

rm(list = ls())
library(tidyverse)
library(ggplot2)

save.fig <- F

# Comparing scaling probability and scaling counts:
prob.results <- read.csv("Data/Calf_Estimates_v3_Mv1_2025-11-13.csv")
count.results <- read.csv("Data/Calf_Estimates_v3_Mv1b_2025-11-13.csv")

all.results <- rbind(prob.results, count.results)
p.prob.count <- ggplot(all.results) +
  geom_point(aes(x = Year, y = Median, color = Method)) +
  geom_errorbar(aes(x = Year, ymin = LCL, ymax = UCL, color = Method))

if (save.fig)
  ggsave(filename = "figures/prob-vs-counts.png",
         plot = p.prob.count,
         device = "png",
         dpi = 600)

# find which years have huge estimates from counts:
count.results %>% 
  left_join(prob.results, by = "Year") %>% #-> tmp
  mutate(d.Mean = Mean.x - Mean.y,
         d.Median = Median.x - Median.y) %>%
  rename(p.Median = Median.y,
         c.Median = Median.x,
         p.Mean = Mean.y,
         c.Mean = Mean.x) %>%
  select(Year, p.Median, c.Median, p.Mean, c.Mean,
         d.Mean, d.Median) -> d.results
  
model <- "v1b"
all.years <- c(1994:2019, 2021:2025)
all.results <- list()
pareto.results <- list()
Rhats.results <- list()
dir.name <- "RData/"

for (k in 1:length(all.years)){
  tmp <- readRDS(file = paste0(dir.name, "calf_estimates_v3_M",
                               model, "_", all.years[k], ".rds"))
  
  count.obs <- tmp$jags.data$count.obs
  MCMC.diag <- tmp$MCMC.diag
  pareto <- MCMC.diag$loo.out$diagnostics$pareto_k
  pareto[is.infinite(pareto)] <- 5
  week <- tmp$jags.data$week
  
  pareto.results[[k]] <- data.frame(count = count.obs,
                                    Pareto = pareto,
                                    week = week,
                                    Year = all.years[k])
  
  Rhats.results[[k]] <- MCMC.diag$max.Rhat
}

pareto.results.df <- do.call(rbind, pareto.results)

p.pareto <- ggplot(pareto.results.df) +
  geom_point(aes(x = week, y = Pareto, color = count),
             alpha = 0.5) + 
  geom_hline(yintercept = 0.7, color = "yellow", size = 1) +
  geom_hline(yintercept = 1.0, color = "red", size = 1.2) +
  facet_wrap(~ Year)

if (save.fig)
  ggsave(filename = paste0("figures/Pareto-k_", model, ".png"),
         plot = p.pareto,
         device = "png", 
         dpi = 600)

Rhats.results.df <- do.call(rbind, Rhats.results)
colnames(Rhats.results.df) <- names(MCMC.diag$max.Rhat)



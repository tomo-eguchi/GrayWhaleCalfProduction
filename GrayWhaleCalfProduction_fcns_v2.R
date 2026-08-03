#This includes updated functions from TE downloaded June 2026
#Excel file definitions (that need to be updated each year) moved to "ExcelFileParams.R"
library(rstan)
library(loo)
library(readxl)

source("ExcelFileParams.R")


# Some functions that were made by Claude:
# 
# =====================================================================
# Fitted seasonal curves for the NB fixed-year and shifting-peak models
# =====================================================================
# Returns, for a chosen year, the fitted per-interval passage curve across
# weeks with a credible ribbon (computed per posterior draw, then summarised).
#
# Curve definition (both models), on the TRUE-passage scale:
#     mu(y, w) = exp( year_eff[y] + shape(y, w) )
# i.e. expected pairs per FULLY-OBSERVED 3-h interval, BEFORE detection/effort
# thinning. (For the expected *observed* count at full effort, multiply by p_obs.)
# The two models differ only in shape(y, w):
#     fixed-year      : shape = week_eff[w]                 (discrete, shared)
#     shifting-peak   : shape = amp * (g(w) - gbar[y])      (smooth, year-specific peak)
#
# include_level = TRUE  -> mu on the per-interval mean scale (level included)
# include_level = FALSE -> exp(shape): the seasonal MULTIPLIER only (level stripped),
#                          for comparing curve SHAPE independent of abundance.
#
# Requires: posterior. (dplyr/ggplot2 only for the example at the bottom.)
# =====================================================================

## Coerce any fit/draws object to a draws_matrix once; index columns by name.
##   M <- as_draws_matrix(fit$draws())     # cmdstanr
as_M <- function(draws) posterior::as_draws_matrix(draws)

# ------------------------------------------------------------------
# FIXED-YEAR: shape is the shared, discrete week_eff.
# Evaluated at INTEGER weeks only (week_eff is a free per-week value; there is
# no smooth interpolation in the model, so don't invent one).
# ------------------------------------------------------------------
# M        : draws_matrix
# year_i   : integer index matching year_eff[year_i]  (see year_index() helper)
# n_weeks  : global number of weeks (dat$n_weeks); week_eff[1..n_weeks]
# weeks    : integer weeks to evaluate (default 1:n_weeks)
curve_fixedyear <- function(M, year_i, n_weeks, 
                            weeks = 1:n_weeks,
                            include_level = TRUE, 
                            probs = c(0.05, 0.5, 0.95)) {
  ye <- if (include_level) 
    as.numeric(M[, sprintf("year_eff[%d]", year_i)]) else 0
  
  rows <- lapply(weeks, function(w) {
    we <- as.numeric(M[, sprintf("week_eff[%d]", w)])
    q  <- quantile(exp(ye + we), probs)
    data.frame(week = w, 
               lower = q[[1]], 
               median = q[[2]], 
               upper = q[[3]],
               model = "fixed-year", 
               row.names = NULL)
  })
  
  return(do.call(rbind, rows))
}

# ------------------------------------------------------------------
# SHIFTING-PEAK: shape is the mean-centered Gaussian pulse.
# CRITICAL: gbar is the mean of the raw Gaussian over the INTEGER week grid
# 1:n_weeks, exactly as in the Stan model. Computing it over the fine plotting
# grid instead would shift the curve vertically and break comparability with
# year_eff (and hence with the fixed-year curve). Do not change this.
# ------------------------------------------------------------------
# grid : fine week grid for a smooth curve (default step 0.1)
# Uses a single shared `width`. For the skew variant, replace `width` below with
# width_l / width_r selected by sign of (w - peak) in BOTH the gbar loop and g.
curve_shiftpeak <- function(M, year_i, 
                            n_weeks, 
                            grid = seq(1, n_weeks, by = 0.1),
                            include_level = TRUE, 
                            probs = c(0.05, 0.5, 0.95)) {
  ye   <- if (include_level) 
    as.numeric(M[, sprintf("year_eff[%d]", year_i)]) else 0
  amp  <- as.numeric(M[, "amp"])
  wid  <- as.numeric(M[, "width"])
  peak <- as.numeric(M[, sprintf("peak_week[%d]", year_i)])
  
  # gbar[draw] = mean over integer weeks 1:n_weeks of exp(-0.5*((k-peak)/wid)^2)
  g_int <- vapply(1:n_weeks,
                  function(k) exp(-0.5 * ((k - peak) / wid)^2),
                  numeric(length(peak)))          # draws x n_weeks
  gbar  <- rowMeans(g_int)                          # draws-length
  
  rows <- lapply(grid, function(w) {
    g    <- exp(-0.5 * ((w - peak) / wid)^2)        # draws-length
    seas <- amp * (g - gbar)                        # mean-centred pulse
    q    <- quantile(exp(ye + seas), probs)
    data.frame(week = w, 
               lower = q[[1]], 
               median = q[[2]], upper = q[[3]],
               model = "shifting-peak", 
               row.names = NULL)
  })
  return(do.call(rbind, rows))
}

# ------------------------------------------------------------------
# Thin dispatcher, so you can call one function with model = "..."
# ------------------------------------------------------------------
season_curve <- function(M, 
                         model = c("fixedyear", "shiftingpeak"),
                         year_i, n_weeks, 
                         include_level = TRUE, ...) {
  model <- match.arg(model)
  if (model == "fixedyear")
    curve_fixedyear(M, year_i, 
                    n_weeks, 
                    include_level = include_level, ...)
  else
    curve_shiftpeak(M, year_i, 
                    n_weeks, 
                    include_level = include_level, ...)
}

# ------------------------------------------------------------------
# Helper: map a calendar year to its year_eff index (2020 is missing, so the
# index is NOT year - 1993). Pass the same year vector used to build the model.
#   year_levels <- sort(unique(raw$Year))
#   year_index(2015, year_levels)  ->  the integer for year_eff[.]
# ------------------------------------------------------------------
year_index <- function(year, 
                       year_levels) match(year, year_levels)


# =====================================================================
# EXAMPLE — the 2014/2015 divergence figure
# =====================================================================
# Overlays the two fitted curves for one year plus the observed passage timing.
# Requires the fixed-year fit (M_fix) and the shifting-peak fit (M_shift), each
# from the SAME data (same n_weeks, same year ordering).
#
# library(dplyr); library(ggplot2)
#
# M_fix   <- as_M(fit_fixedyear$draws())
# M_shift <- as_M(fit_shiftpeak$draws())
# year_levels <- sort(unique(raw$Year))
# n_weeks <- dat$n_weeks
# yr <- 2015; yi <- year_index(yr, year_levels)
#
# cf <- curve_fixedyear(M_fix,   yi, n_weeks)                 # integer weeks
# cs <- curve_shiftpeak(M_shift, yi, n_weeks)                 # smooth grid
# curves <- rbind(cf, cs)
#
# ## observed passage timing for the same year: count per effort by week.
# ## This is ~proportional to mu * p_obs, so it shows WHERE the peak is, not its
# ## absolute height -> put it on a secondary axis rather than forcing one scale.
# obs <- raw %>% filter(Year == yr, effort > 0) %>%
#   group_by(week) %>% summarise(cpe = sum(count) / sum(effort), .groups = "drop")
# scl <- max(curves$median) / max(obs$cpe)                    # visual scale only
#
# ggplot(curves, aes(week)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = .20) +
#   geom_line(aes(y = median, colour = model), linewidth = .8) +
#   geom_point(data = obs, aes(x = week, y = cpe * scl),
#              inherit.aes = FALSE, colour = "black") +
#   scale_y_continuous(
#     name = "fitted pairs per interval (true passage)",
#     sec.axis = sec_axis(~ . / scl, name = "observed count / effort")) +
#   labs(x = "week (model index)", title = paste("Seasonal curve —", yr),
#        subtitle = "fixed-year peaks at the common week; shifting-peak follows this year's data") +
#   theme_minimal()
#
# What to look for: the shifting-peak median (smooth) peaking at a DIFFERENT week
# than the fixed-year median (integer-week line), with the observed points
# supporting the shifted position -> the mechanism behind the 2015 disagreement.
#
# Tip: to compare SHAPE only (ignore that 2015 was a high-abundance year), call
# both curve fns with include_level = FALSE and plot the seasonal multiplier.

# Added 2026-06-30 TE
# Using the stan output calf abundance estimates are computed using posterior samples for single year
# The Stan model is GWCalfCount_nb_singleyear.stan
# The model is ran in calf_production_single_year_stan.R
# =======================================================
# POST-PROCESSING: Reconstruct True Counts via Conjugacy
# =======================================================
stan.post.process.1year <- function(stan.fit, stan.data, pre.stan.data, out.file.name, save.file = F){
  fit_ <- stan.fit
  
  mu_true_samples <- fit_$draws("mu_true", format = "matrix")
  p_obs_samples   <- fit_$draws("p_obs_corr", format = "matrix")
  phi_samples     <- as.vector(fit_$draws("phi", format = "matrix"))
  
  n_iterations <- nrow(mu_true_samples)
  n_obs        <- stan.data$n_obs
  
  # Create matrix equivalents for element-wise array operations
  count_obs_matrix <- matrix(rep(stan.data$count_obs, 
                                 each = n_iterations), 
                             nrow = n_iterations)
  phi_matrix       <- matrix(rep(phi_samples, times = n_obs), 
                             nrow = n_iterations)
  
  # Gamma Conjugacy Parameters
  shape_matrix <- phi_matrix + count_obs_matrix
  rate_matrix  <- (phi_matrix / mu_true_samples) + p_obs_samples
  
  # Draw the continuous underlying latent rates (lambda)
  lambda_samples <- matrix(
    rgamma(length(shape_matrix), 
           shape = shape_matrix, 
           rate = rate_matrix),
    nrow = n_iterations, ncol = n_obs
  )
  
  # Draw missed calves from the unobserved portion of the process
  mu_missed <- lambda_samples * (1 - p_obs_samples)
  count_missed_samples <- matrix(
    rpois(length(mu_missed), lambda = mu_missed),
    nrow = n_iterations, ncol = n_obs
  )
  
  # True Integer Counts = Observed + Simulated Missed
  count_true_samples <- count_obs_matrix + count_missed_samples
  
  # 4. Summarize Annual Trajectories
  summary_stats <- list()
  annual_posterior_samples <- rowSums(count_true_samples)
    
  summary_stats <- data.frame(
    Mean   = mean(annual_posterior_samples),
    Median = median(annual_posterior_samples),
    SE = sqrt(var(annual_posterior_samples)),
    Lower_95_CI = quantile(annual_posterior_samples, 
                           probs = 0.025),
    Upper_95_CI = quantile(annual_posterior_samples, 
                           probs = 0.975)
  )
  
  #print(annual_abundance)
  if (save.file)
    write.csv(summary_stats,
              file = paste0("data/", out.file.name, ".csv"),
              quote = FALSE)
  
  return(list(abundance = summary_stats))
}


# Added by TE 2026-06-26
# Using the stan output calf abudnance estimates are computed using posterior samples
# The Stan model is GWCalfCount_nb.stan
# The model is ran in calf_production_stan_2026_v4TE.R
# =======================================================
# POST-PROCESSING: Reconstruct True Counts via Conjugacy
# =======================================================
stan.post.process <- function(stan.fit, stan.data, pre.stan.data, out.file.name, save.file = F, save.fig = F){
  fit_ <- stan.fit
  
  # pull out posterior samples from the stan.fit object
  mu_true_samples <- fit_$draws("mu_true", format = "matrix")
  p_obs_samples   <- fit_$draws("p_obs_corr", format = "matrix")
  phi_samples     <- as.vector(fit_$draws("phi", format = "matrix"))
  
  n_iterations <- nrow(mu_true_samples)
  n_obs        <- stan.data$n_obs
  
  # Create matrix equivalents for element-wise array operations
  count_obs_matrix <- matrix(rep(stan.data$count_obs, 
                                 each = n_iterations), 
                             nrow = n_iterations)
  phi_matrix       <- matrix(rep(phi_samples, times = n_obs), 
                             nrow = n_iterations)
  
  # Gamma Conjugacy Parameters
  shape_matrix <- phi_matrix + count_obs_matrix
  rate_matrix  <- (phi_matrix / mu_true_samples) + p_obs_samples
  
  # Draw the continuous underlying latent rates (lambda)
  lambda_samples <- matrix(
    rgamma(length(shape_matrix), 
           shape = shape_matrix, 
           rate = rate_matrix),
    nrow = n_iterations, ncol = n_obs
  )
  
  # Draw missed calves from the unobserved portion of the process
  mu_missed <- lambda_samples * (1 - p_obs_samples)
  count_missed_samples <- matrix(
    rpois(length(mu_missed), lambda = mu_missed),
    nrow = n_iterations, ncol = n_obs
  )
  
  # True Integer Counts = Observed + Simulated Missed
  count_true_samples <- count_obs_matrix + count_missed_samples
  
  # 4. Summarize Annual Trajectories
  summary_stats <- list()
  annual_posterior_samples <- list()
  for (y in 1:stan.data$n_years) {
    matching_rows <- which(pre.stan.data$year_idx == y)
    annual_posterior_samples[[y]] <- rowSums(count_true_samples[, matching_rows])
    calendar_year_label <- unique(pre.stan.data$calendar_year[pre.stan.data$year_idx == y])
    
    summary_stats[[y]] <- data.frame(
      Year        = calendar_year_label,
      Mean   = mean(annual_posterior_samples[[y]]),
      Median = median(annual_posterior_samples[[y]]),
      SE = sqrt(var(annual_posterior_samples[[y]])),
      Lower_95_CI = quantile(annual_posterior_samples[[y]], 
                             probs = 0.025),
      Upper_95_CI = quantile(annual_posterior_samples[[y]], 
                             probs = 0.975)
    )
  }
  
  annual_abundance <- do.call(rbind, summary_stats)
  rownames(annual_abundance) <- NULL
  
  #print(annual_abundance)
  if (save.file)
    write.csv(annual_abundance,
              file = paste0("data/", out.file.name, ".csv"),
              quote = FALSE)
  
  p.plot <- ggplot(annual_abundance, aes(x = Year, y = Median)) +
    geom_ribbon(aes(ymin = Lower_95_CI, ymax = Upper_95_CI),
                fill = "skyblue", alpha = 0.4) +
    geom_line(color = "blue", linewidth = 1) +
    geom_point(color = "darkblue", size = 2) +
    labs(
      title = "Estimated Annual Calf Abundance",
      x = "Year",
      y = "Total Estimated Calves"
    ) +
    theme_minimal()
  
  if (save.fig)
    ggsave(filename = "figures/calf_abundance.png",
           plot = p.plot,
           device = "png",
           dpi = 600)
  
  return(list(abundance = annual_abundance,
              plot = p.plot,
              posterior.samples = annual_posterior_samples))
}

# Added by TE 2026-06-26
# Posterior Prdictive Check of the Stan output from GWCalfCount_nb.stan
# Use calf_production_stan_2026_v4TE.R to run the Stan model. 
PPC_counts <- function(stan.fit, stan.data, save.fig = F){
  # 1. Extract the simulated replication matrix [8000 iterations x 17008 columns]
  y_rep <- stan.fit$draws("count_rep", format = "matrix")
  
  # 2. Identify the active survey days
  real_days <- which(stan.data$effort > 0)
  
  # 3. Subset both your actual data and your simulated data to ONLY include active days
  y_observed_clean <- stan.data$count_obs[real_days]
  y_rep_clean      <- y_rep[, real_days]
  
  # 4. Plot Density Overlay
  p.ppc.bars <- ppc_bars(y = y_observed_clean, 
                         yrep = y_rep_clean[1:2000, ]) +
    coord_cartesian(xlim = c(-0.5, 6.5)) + 
    labs(title = "Posterior Predictive Check: Discrete Sighting Counts",
         x = "Number of mother-calf pairs observed",
         y = "Total number of 3-hr intervals")
  
  # 5. Plot Summary Statistics Check
  # This evaluates if your model successfully matches the proportion of zeros 
  # observed in the field on active watch days.
  prop_zero <- function(x) mean(x == 0)
  p.ppc.stat <- ppc_stat(y = y_observed_clean, 
                         yrep = y_rep_clean, stat = "prop_zero") +
    labs(title = "Proportion of Zeros: Observed vs. Simulated")
  
  return(list(bars = p.ppc.bars,
              stat = p.ppc.stat,
              y.rep = y_rep_clean,
              y.obs = y_observed_clean))
  
  if (save.fig){
    ggsave(p.ppc.bars,
           filename = "figures/ppc_bars.png",
           device = "png",
           dpi = 600)
    ggsave(p.ppc.stat,
           filename = "figures/ppc_stat.png",
           device = "png",
           dpi = 600)
  }
}

# Removed compute.LOOIC as this can be done easily using the posterior package
# Removed rank.normalized.R.hat as this can be done easily using the posterior package

get.all.data <- function(sheet.name.inshore, sheet.name.offshore,
                         col.types.inshore, col.types.offshore,
                         col.names.inshore, col.names.offshore,
                         Year, 
                         xls.file.name.inshore, xls.file.name.offshore,
                         start.time){
  
  data.inshore <- get.data(Year = Year, 
                           xls.file.name = xls.file.name.inshore, 
                           sheet.name = sheet.name.inshore,
                           col.types = col.types.inshore, 
                           col.names = col.names.inshore, 
                           start.time = start.time) %>%
    mutate(Area = "in") %>%
    extract.all.vars()
  
  data.offshore <- get.data(Year = Year, 
                            xls.file.name = xls.file.name.offshore, 
                            sheet.name = sheet.name.offshore,
                            col.types = col.types.offshore, 
                            col.names = col.names.offshore, 
                            start.time = start.time)%>%
    mutate(Area = "off") %>%
    extract.all.vars()
  
  # data.offshore %>% 
  #   filter(Mother_Calf > 0) %>%
  #   rbind(data.inshore) %>%
  #   arrange(Date.date, Minutes_since_0000) %>% 
  #   extract.all.vars() -> data.all
  
  #shift.offshore <- find.effort(data.offshore, T0 = start.time) 
  #shift.inshore <- find.effort(data.inshore, T0 = start.time) 
  
  shift.offshore <- find.effort.2(data.offshore) 
  shift.inshore <- find.effort.2(data.inshore) 
  
  # formatted.offshore <- format.output(shift.offshore$out.df, 
  #                                     max.shift = max(shift.offshore$out.df$Shift))
  # formatted.inshore <- format.output(shift.inshore$out.df, 
  #                                    max.shift = max(shift.inshore$out.df$Shift))
  
  formatted.offshore <- format.output.2(shift.offshore$out.df)
  formatted.inshore <- format.output.2(shift.inshore$out.df)

  return(list(data.offshore = data.offshore,
              data.inshore = data.inshore,
              shift.data.offshore = shift.offshore$out.df,
              shift.data.inshore = shift.inshore$out.df,
              formatted.data.offshore = formatted.offshore,
              formatted.data.inshore = formatted.inshore,
              shift.all.offshore = shift.offshore$shift.df,
              shift.all.inshore = shift.inshore$shift.df))
}


# col.names has to be in the same order of:
# c("Date", "Event", "Time", "Obs. Code", 
# "Sea State", "Vis. IN", "Cow  / Calf")
# The exact names may be different in each file. Make sure to 
# match them. 

# T0 is the beginning of the first shift (start.time)

get.data <- function(Year, xls.file.name, sheet.name,
                     col.types, col.names, start.time){
  
  if (length(col.names) < 8){
    data.out <- read_excel(xls.file.name,
                           sheet = sheet.name,
                           col_types = col.types,
                           col_names = TRUE) %>% 
      dplyr::select(all_of(col.names)) %>%
      transmute(Date = .data[[col.names[[1]]]],
                Event = .data[[col.names[[2]]]],
                Time = .data[[col.names[[3]]]],
                Minutes_since_T0 = sapply(Time, FUN = char_time2min, start.time),
                Minutes_since_0000 = sapply(Time, FUN = char_time2min),
                #Shift = sapply(Minutes_since_0000, FUN = find.shift),
                Shift = sapply(Minutes_since_0000, FUN = find.shift.2),
                Obs = .data[[col.names[[4]]]],
                SeaState = .data[[col.names[[5]]]],
                Vis = .data[[col.names[[6]]]],
                Mother_Calf = .data[[col.names[[7]]]],
                Year = Year,
                Area = "off")
     
  } else {
    data.out <- read_excel(xls.file.name,
                           sheet = sheet.name,
                           col_types = col.types,
                           col_names = TRUE) %>% 
      dplyr::select(all_of(col.names)) %>%
      transmute(Date = .data[[col.names[[1]]]],
                Event = .data[[col.names[[2]]]],
                Time = .data[[col.names[[3]]]],
                Minutes_since_T0 = sapply(Time, FUN = char_time2min, start.time),
                Minutes_since_0000 = sapply(Time, FUN = char_time2min),
                #Shift = sapply(Minutes_since_0000, FUN = find.shift),
                Shift = sapply(Minutes_since_0000, FUN = find.shift.2),
                Obs = .data[[col.names[[4]]]],
                SeaState = .data[[col.names[[5]]]],
                Vis = .data[[col.names[[6]]]],
                Mother_Calf = .data[[col.names[[7]]]],
                Year = Year,
                Area =  .data[[col.names[[8]]]])
  }
  # Some files contain wrong years...
  years <- year(data.out$Date)
  dif.years <- sum(years - Year, na.rm = T)
  if (dif.years != 0){
    year(data.out$Date) <- Year
  } 

  # filter out rows with NAs in certain fields, then order them using
  # Date and time since 0000 .
  # Create Date fields with character and date format - may not be necessary?
  data.out %>% 
    filter(!is.na(Date)) %>%
    filter(!is.na(Event)) %>%
    filter(!is.na(Shift)) %>%
    arrange(Date, Minutes_since_0000) %>%
    mutate(Date.date = as.Date(Date),
           Date.char = as.character(Date)) %>%
    dplyr::select(-Date) %>%
    relocate(Date.date) -> data.out
  
  # fill in sea state and visibility
  sea.state <- data.out[1, "SeaState"]
  vis <- data.out[1, "SeaState"]
  for (k6 in 2:nrow(data.out)){
    if (is.na(data.out[k6, "SeaState"])){
      data.out[k6, "SeaState"] <- sea.state
    } else {
      sea.state <- data.out[k6, "SeaState"]
    }
    
    if (is.na(data.out[k6, "Vis"])){
      data.out[k6, "Vis"] <- vis
    } else {
      vis <- data.out[k6, "Vis"]
    }
    
  }
  
  return(data.out)  
}

get.data.inshore.only <- function(sheet.name, 
                                  col.types, 
                                  col.names, 
                                  Year, 
                                  xls.file.name, 
                                  start.time){
  
  data.inshore <- get.data(Year = Year,
                           xls.file.name = xls.file.name,
                           sheet.name = sheet.name,
                           col.types = col.types,
                           col.names = col.names,
                           start.time = start.time) %>%
    extract.all.vars()
  
  # calculate effort and other statistics per 3-hr shift, using find.effort
  # function in Piedras_Blancas_fcns.R
  #data.shift <- find.effort(data.inshore, T0 = start.time) #%>% extract.shift.vars()
  
  data.shift <- find.effort.2(data.inshore) 
  
  # formatted.inshore <- format.output(data.shift$out.df, 
  #                                    max.shift = max(data.shift$out.df$Shift))
  
  formatted.inshore <- format.output.2(data.shift$out.df)
  out.list <- list(data.inshore = data.inshore,
                   shift.data.inshore = data.shift$out.df,
                   formatted.data.inshore = formatted.inshore)
  return(out.list)
  
}

# x is a data.frame with at least three fields: Date (character), 
# Minutes_since_0000
# and Shift - find.shift is used to come up with this (or find.shift.2)
# 1 = START EFFORT
#2 = CHANGE OBSERVERS
#3 = CHANGE SIGHTING CONDITIONS
#4 = GRAY WHALE SIGHTING
#5 = END EFFORT
#6 = OTHER SPECIES SIGHTING
# provide a data frame that came back from get.data function and the start
# time of the first shift (T0) in character.
# find.effort <- function(x, T0){
# 2023-06-01 A HA MOMENT. RATHER THAN MAKING ARBITRARY STARTING TIME OF 0630 OR 0700, 
# CREATE 3-HR SHIFTS STARTING FROM 0100. This has not been implemented. The comparison
# between Ver1 and Ver2 needs to be completed first. 

# T0 is start time in character, e.g., "0630", "0700)
find.effort <- function(x, T0){
  
  # NEED TO ADD SHIFT INDEX. TO DO THAT, I NEED TO KNOW WHAT TIME THE OBSERVATION
  # STARTS EVERYDAY AND THE MAXIMUM NUMBER OF SHIFTS PER DAY, WHICH ARE NOT CONSISTENT
  # AMONG YEARS. I NEED TO THINK ABOUT THIS A BIT MORE... 2023-05-01
  # 
  
  # THE FUNCTION ABOVE find.shift DEFINES SHIFTS. ANYTIME BEFORE 1000 IS SHIFT 1, ETC.
  # get.data USES find.shift AND RETURNS ASSIGNED SHIFTS IN THE OUTPUT.
  # ALTERNATIVELY, find.shift.2 DEFINES SHIFTS DIFFERENTLY. IT USES A DIFFERENT START TIME
  # WHERE 0100 IS THE BEGINNING OF EACH DAY. THIS WAY, OCCASIONAL START OF THE FIRST SHIFT
  # BEFORE 0700 CAN BE ASSIGNED TO THE CORRECT SHIFT BIN. IN THIS WAY, START TIME (T0)
  # SHOULD NOT MATTER ANY MORE. 2023-06-02 - SEE find.effort.2 FOR THIS. 
  
  # REGARDLESS OF STARTING TIME, THEY SHOULD BE SEPARATED INTO 3-HR BLOCKS
  
  # 2023-05-10 THIS IS STILL NOT WORKING PROPERLY. EFFORT IS NOT SUMMED CORRECTLY, PARTIALLY
  # DUE TO NON-MATCHING 1/5 EVENT CODE (ON AND OFF EFFORT MARKERS). SEA STATE AND VISIBILITY
  # FILTERS NEED TO BE MORE EFFICIENT. THIS HAS BEEN FIXED.
  
  # this turns dates in to character
  all.dates <- unique(x$Date)
  
  # Start time of shift 1:
  T0.minutes <- char_time2min(T0)
  
  # End of shift time in minutes since 0100. 600 is 1000hrs
  shift.ends <- c(600, 780, 960, 1140, 1320)
  
  shift.begins <- c(T0.minutes, shift.ends)
  
  out.list <- list()
  shift.list <- list()
  d <- k1 <- c <- k2 <- 1
  d <- 13
  k4 <- 1
  for (d in 1:length(all.dates)){
    # pick just one day's worth of data
    one.day <- filter(x, Date == as.Date(all.dates[d])) %>%
      arrange(Minutes_since_0000)
    
    # Must have at least one Event = 1
    one.day %>% filter(Event == 1) -> one.day.event.1
    
    # Entries when there were no shift at all (e.g., 1998-04-28) need to be removed
    # also, an entry before effort starts (e.g., 1998-03-23 offshore)
    if (nrow(one.day) > 1 & one.day$Event[1] != 5 & nrow(one.day.event.1) > 0){
      # Find start and end time for the day:
      one.day %>%
        arrange(by = Minutes_since_0000) -> one.day
      
      # Visibility is coded in characters. This needs to be converted into double in order
      # to find maximum value within each shift. 2023-05-10
      all.sea.state <- one.day$SeaState
      idx.slash <- grep("/", all.sea.state)
      if (length(idx.slash) == 0){
        one.day$SeaState.num <- as.numeric(all.sea.state)
      } else {
        # Should I take the first one, or maximum? I go with first one for now. 2023-05-11
        one.day$SeaState.num <- lapply(str_split(all.sea.state, "/"), 
                                       FUN = function(x) x[1]) %>% 
          unlist() %>% 
          as.numeric()
      }
      
      all.vis <- one.day$Vis
      idx.slash <- grep("/", all.vis)
      if (length(idx.slash) == 0){
        one.day$Vis.num <- as.numeric(all.vis)
      } else {
        one.day$Vis.num <- lapply(str_split(all.vis, "/"), 
                                  FUN = function(x) x[1]) %>% 
          unlist() %>% 
          as.numeric()
      }
      
      # Check to see if the line is the end of one and the beginning of the next
      shift <- one.day$Shift[1]
      if (length(grep("/", shift)) == 0){
        shift.first <- as.numeric(shift)
      } else {
        shift.first <- strsplit(shift, "/")[[1]][2] %>% as.numeric()
      }
      
      shift <- one.day$Shift[nrow(one.day)]
      if (length(grep("/", shift)) == 0){
        shift.last <- as.numeric(shift)
      } else {
        shift.last <- strsplit(shift, "/")[[1]][1] %>% as.numeric()
      }
      
      # go through one shift at a time but remove the 5th shift
      if (shift.last > 4) shift.last <- 4
      
      for (k4 in shift.first:shift.last){
        if (k4 == 1){
          # for the first shift, sometimes they put a comment in the first row
          # with Event = 6, rather than Event = 1. Those need to be removed.
          # take all Event == 1
          tmp <- one.day %>% 
            filter(Event == 1)
          
          # Some years contain "shifts" with no observations but just comments
          # Event = 6. E.g., 1997-03-24. These need to be dealt with. 
          #if (nrow(tmp) == 0){   # no start events
            #START HERE NEXT 2023-05-23
          #}
          T0000.begin <- min(tmp$Minutes_since_0000, na.rm = T) 
          
          # if (nrow(tmp) == 0) stop("error")
        } else {
          T0000.begin <- shift.begins[k4]  
        }
        
        T0000.end <- shift.ends[k4]
        
        if (k4 == 1){
          one.day %>% filter(Minutes_since_0000 >= (T0000.begin),
                             Minutes_since_0000 <= T0000.end) -> one.shift  
        } else {
          one.day %>% filter(Minutes_since_0000 > (T0000.begin),
                             Minutes_since_0000 <= T0000.end) -> one.shift
          
        }
        
        # This will pick up the shift plus (shift)/(shift+1) 
        # shift.idx <- grep(as.character(k4), one.day$Shift)
        # 
        # one.shift <- one.day[shift.idx,]
        
        if (nrow(one.shift) > 0){
          # Sometimes the beginning of one shift and the end of the previous 
          # shift is shared in one line with Shift = x/y. When that and a sighting
          # happens simultaneously, the sighting gets double counted between the
          # two shifts. So, the sighting has to be placed in one or the other.
          line.bottom <- one.shift[nrow(one.shift),]
          if (line.bottom$Event == 4 & str_detect(line.bottom$Shift, "/")){
            one.shift <- one.shift[1:(nrow(one.shift)-1),]
            
          }
          
          # add one row at the top and end of one.shift, so that it has
          # Event 1 at the top, and 5 at the bottom if they are not there:
          if (one.shift[1,"Event"] != 1){
            one.shift.eft <- rbind(one.shift[1,], one.shift)
            one.shift.eft[1, "Event"] <- 1
            one.shift.eft[1, "Time"] <- minutes2time_char(shift.begins[k4])
            one.shift.eft[1, "Minutes_since_T0"] <- shift.begins[k4] - char_time2min(T0)
            one.shift.eft[1, "Minutes_since_0000"] <- shift.begins[k4]
            one.shift.eft[1, "Mother_Calf"] <- NA

          } else {
            one.shift.eft <- one.shift
          }
          
          if (one.shift[nrow(one.shift), "Event"] != 5){
            one.shift.eft <- rbind(one.shift.eft, one.shift[nrow(one.shift),])
            one.shift.eft[nrow(one.shift.eft), "Event"] <- 5
            one.shift.eft[nrow(one.shift.eft), "Time"] <- minutes2time_char(shift.ends[k4])
            one.shift.eft[nrow(one.shift.eft), "Minutes_since_T0"] <- shift.ends[k4] - char_time2min(T0)
            one.shift.eft[nrow(one.shift.eft), "Minutes_since_0000"] <- shift.ends[k4]
            one.shift.eft[nrow(one.shift.eft), "Mother_Calf"] <- NA 
          }
          
          # Figure out off effort time due to visibility and sea state (> 4)
          one.shift.eft %>% 
            mutate(Effort = ifelse((Vis.num > 4 | SeaState.num > 4),
                                   "off", "on")) -> one.shift.eft
          
          # # find out how many on/off effort existed
          row.1 <- which(one.shift.eft$Event == 1)
          row.5 <- which(one.shift.eft$Event == 5)
          
          # sometimes too many event == 1 and event == 5
          # 
          if (length(row.1) > 1 & length(row.5) > 1){
            for (k3 in 2:length(row.1)){
              if (row.1[k3] < row.5[k3-1]){
                stop()
                row.1[k3] <- NA
              }
              
            }
            row.1 <- row.1[!is.na(row.1)]
            
          }
          
          # if they don't match, adjust accordingly.
          if (length(row.1) != length(row.5)){
            nrow <- min(length(row.1), length(row.5))
            row.1.1 <- vector(mode = "numeric", length = nrow)
            row.5.1 <- vector(mode = "numeric", length = nrow)
            for (k2 in 1:nrow){
              if (k2 == 1){
                row.1.1[k2] <- row.1[k2]
                row.5.1[k2] <- row.5[k2]
              } else {
                row.1.1[k2] <- first(row.1[row.1 > row.5.1[k2-1]])
                row.5.1[k2] <- first(row.5[row.5 > row.1.1[k2]])
              }
              
            }
            row.1 <- row.1.1[!is.na(row.1.1)]
            row.5 <- row.5.1[!is.na(row.5.1)]
          }
          
          # These are probably unnecessary?
          one.shift.eft[row.1, "Event"] <- 1
          one.shift.eft[row.5, "Event"] <- 5
          
          one.shift.eft[row.1, "Mother_Calf"] <- NA
          one.shift.eft[row.5, "Mother_Calf"] <- NA
          
          one.shift.eft <- one.shift.eft[!is.na(one.shift.eft$Effort),]
          # calculate effort for each "on" period per shift
          #tmp.eft <- 0
          idx.off <- which(one.shift.eft$Effort == "off")
          idx.on <- which(one.shift.eft$Effort == "on")
          
          if (length(idx.on) > 0){
            for (k2 in min(idx.on):nrow(one.shift.eft)){
              if (k2 == min(idx.on)){
                tmp.eft <- 0
                time.0 <- one.shift.eft[k2, "Minutes_since_0000"]
                effort <- "on"
              } else {
                if (one.shift.eft[k2, "Effort"] == "off" & effort == "on"){
                  d.time <- one.shift.eft[k2, "Minutes_since_0000"] - time.0
                  tmp.eft <- tmp.eft + d.time
                  effort <- "off"
                } else if (one.shift.eft[k2, "Effort"] == "on" & effort == "off"){
                  time.0 <- one.shift.eft[k2, "Minutes_since_0000"]
                  effort <- "on"
                } else if (one.shift.eft[k2, "Effort"] == "on" & effort == "on"){
                  d.time <- one.shift.eft[k2, "Minutes_since_0000"] - time.0
                  tmp.eft <- tmp.eft + d.time
                  effort <- "on"
                  time.0 <- one.shift.eft[k2, "Minutes_since_0000"]
                }
              }
              
            }
            
          } else {
            tmp.eft <- 0
          }
          
          # Sometimes, there is no Vis (or also sea state) update in an entire
          # shift. Use one from previous shift in those cases
          # These should not happen any longer but keep it anyways. 2023-05-10
          if (sum(!is.na(one.shift$SeaState.num)) > 0){
            max.sea.state <- max(one.shift$SeaState.num, na.rm = T)          
          } else {
            max.sea.state <- max.sea.state
          }
          
          if (sum(!is.na(one.shift$Vis.num)) > 0){
            max.Vis <- max(one.shift$Vis.num, na.rm = T) 
          } else {
            max.Vis <- max.Vis
          }
          
          out.list[[c]] <- data.frame(Date = all.dates[d],
                                      Minutes_since_0000 = one.shift.eft[1, "Minutes_since_0000"],
                                      Shift = k4,
                                      Effort = unname(tmp.eft),
                                      Mother_Calf = one.shift.eft %>%
                                        filter(Effort == "on") %>%
                                        summarize(MandC = sum(Mother_Calf, na.rm = TRUE)) %>%
                                        pull(),
                                      Sea_State = ifelse(!is.infinite(max.sea.state),
                                                         max.sea.state, NA),
                                      Vis = ifelse(!is.infinite(max.Vis),
                                                   max.Vis, NA))
          
          shift.list[[c]] <- one.shift.eft
          c <- c + 1
          
        }
      }  
    }
    
  }
  
  out.df <- do.call(rbind, out.list)
  shift.df <- do.call(rbind, shift.list)
  # out.df %>% 
  #   mutate(Date = as.Date(Date.char, 
  #                              format = "%Y-%m-%d")) -> out.df.1
  return(list(out.df = out.df,
              shift.df = shift.df))
}

# x is a data.frame with at least three fields: Date (character), 
# Minutes_since_0000
# and Shift - find.shift is used to come up with this (or find.shift.2)
# 1 = START EFFORT
#2 = CHANGE OBSERVERS
#3 = CHANGE SIGHTING CONDITIONS
#4 = GRAY WHALE SIGHTING
#5 = END EFFORT
#6 = OTHER SPECIES SIGHTING
# provide a data frame that came back from get.data function
# Unlike find.effort, this one does not need the start time of the first
# shift. The first shift is considered to start at 0100 This way,
# Shift 1: 0100-0400
# Shift 2: 0400-0700
# Shift 3: 0700-1000
# Shift 4: 1000-1300
# Shift 5: 1300-1600
# Shift 6: 1600-1900
# Shift 7: 1900-2200
# Shift 8: 2200-0100
find.effort.2 <- function(x){
  
  # this turns dates in to character
  all.dates <- unique(x$Date)
  
  # Start time of shift 2:
  T0.minutes <- char_time2min("0400")
  
  # End of shift time in minutes since 0100. 600 is 1000hrs
  shift.ends <- c(240, 420, 600, 780, 960, 1140)
  
  shift.begins <- c(T0.minutes, shift.ends)
  
  out.list <- list()
  shift.list <- list()
  d <- k1 <- c <- k2 <- 1
  d <- 13
  k4 <- 1
  for (d in 1:length(all.dates)){
    # pick just one day's worth of data
    one.day <- filter(x, Date == as.Date(all.dates[d])) %>% 
      arrange(Minutes_since_0000)
    
    # Must have at least one Event = 1
    one.day %>% filter(Event == 1) -> one.day.event.1
    
    # Entries when there were no shift at all (e.g., 1998-04-28) need to be removed
    # also, an entry before effort starts (e.g., 1998-03-23 offshore)
    if (nrow(one.day) > 1 & one.day$Event[1] != 5 & nrow(one.day.event.1) > 0){
      # Find start and end time for the day:
      one.day %>%
        arrange(by = Minutes_since_0000) -> one.day
      
      # Visibility is coded in characters. This needs to be converted into double in order
      # to find maximum value within each shift. 2023-05-10
      all.sea.state <- one.day$SeaState
      idx.slash <- grep("/", all.sea.state)
      if (length(idx.slash) == 0){
        one.day$SeaState.num <- as.numeric(all.sea.state)
      } else {
        # Should I take the first one, or maximum? I go with first one for now. 2023-05-11
        one.day$SeaState.num <- lapply(str_split(all.sea.state, "/"), 
                                       FUN = function(x) x[1]) %>% 
          unlist() %>% 
          as.numeric()
      }
      
      all.vis <- one.day$Vis
      idx.slash <- grep("/", all.vis)
      if (length(idx.slash) == 0){
        one.day$Vis.num <- as.numeric(all.vis)
      } else {
        one.day$Vis.num <- lapply(str_split(all.vis, "/"), 
                                  FUN = function(x) x[1]) %>% 
          unlist() %>% 
          as.numeric()
      }
      
      # Check to see if the line is the end of one and the beginning of the next
      shift <- one.day$Shift[1]
      if (length(grep("/", shift)) == 0){
        shift.first <- as.numeric(shift)
      } else {
        shift.first <- strsplit(shift, "/")[[1]][2] %>% as.numeric()
      }
      
      shift <- one.day$Shift[nrow(one.day)]
      if (length(grep("/", shift)) == 0){
        shift.last <- as.numeric(shift)
      } else {
        shift.last <- strsplit(shift, "/")[[1]][1] %>% as.numeric()
      }
      
      for (k4 in shift.first:shift.last){
        if (k4 == shift.first){
          # for the first shift, sometimes they put a comment in the first row
          # with Event = 6, rather than Event = 1. Those need to be removed.
          # take all Event == 1
          tmp <- one.day %>% 
            filter(Event == 1)
          
          # Some years contain "shifts" with no observations but just comments
          # Event = 6. E.g., 1997-03-24. These need to be dealt with. 
          #if (nrow(tmp) == 0){   # no start events
          #START HERE NEXT 2023-05-23
          #}
          T0000.begin <- min(tmp$Minutes_since_0000, na.rm = T) 
          
          # if (nrow(tmp) == 0) stop("error")
        } else {
          T0000.begin <- shift.begins[k4]  
        }
        
        T0000.end <- shift.ends[k4]
        
        if (k4 == shift.first){
          one.day %>% filter(Minutes_since_0000 >= (T0000.begin),
                             Minutes_since_0000 <= T0000.end) -> one.shift  
        } else {
          one.day %>% filter(Minutes_since_0000 > (T0000.begin),
                             Minutes_since_0000 <= T0000.end) -> one.shift
          
        }
        
        # This will pick up the shift plus (shift)/(shift+1) 
        # shift.idx <- grep(as.character(k4), one.day$Shift)
        # 
        # one.shift <- one.day[shift.idx,]
        
        if (nrow(one.shift) > 0){
          # Sometimes the beginning of one shift and the end of the previous 
          # shift is shared in one line with Shift = x/y. When that and a sighting
          # happens simultaneously, the sighting gets double counted between the
          # two shifts. So, the sighting has to be placed in one or the other.
          line.bottom <- one.shift[nrow(one.shift),]
          if (line.bottom$Event == 4 & str_detect(line.bottom$Shift, "/")){
            one.shift <- one.shift[1:(nrow(one.shift)-1),]
            
          }
          
          # add one row at the top and end of one.shift, so that it has
          # Event 1 at the top, and 5 at the bottom if they are not there already:
          if (one.shift[1,"Event"] != 1){
            one.shift.eft <- rbind(one.shift[1,], one.shift)
            one.shift.eft[1, "Event"] <- 1
            one.shift.eft[1, "Time"] <- minutes2time_char(shift.begins[k4])
            one.shift.eft[1, "Minutes_since_T0"] <- shift.begins[k4] - char_time2min("100")
            one.shift.eft[1, "Minutes_since_0000"] <- shift.begins[k4]
            one.shift.eft[1, "Mother_Calf"] <- NA
            
          } else {
            one.shift.eft <- one.shift
          }
          
          if (one.shift[nrow(one.shift), "Event"] != 5){
            one.shift.eft <- rbind(one.shift.eft, one.shift[nrow(one.shift),])
            one.shift.eft[nrow(one.shift.eft), "Event"] <- 5
            one.shift.eft[nrow(one.shift.eft), "Time"] <- minutes2time_char(shift.ends[k4])
            one.shift.eft[nrow(one.shift.eft), "Minutes_since_T0"] <- shift.ends[k4] - char_time2min("100")
            one.shift.eft[nrow(one.shift.eft), "Minutes_since_0000"] <- shift.ends[k4]
            one.shift.eft[nrow(one.shift.eft), "Mother_Calf"] <- NA 
          }
          
          # Figure out off effort time due to visibility and sea state (> 4)
          one.shift.eft %>% 
            mutate(Effort = ifelse((Vis.num > 4 | SeaState.num > 4),
                                   "off", "on")) -> one.shift.eft
          
          # # find out how many on/off effort existed
          row.1 <- which(one.shift.eft$Event == 1)
          row.5 <- which(one.shift.eft$Event == 5)
          
          # sometimes too many event == 1 and event == 5
          # 
          if (length(row.1) > 1 & length(row.5) > 1){
            for (k3 in 2:length(row.1)){
              if (row.1[k3] < row.5[k3-1]){
                stop()
                row.1[k3] <- NA
              }
              
            }
            row.1 <- row.1[!is.na(row.1)]
            
          }
          
          # if they don't match, adjust accordingly.
          if (length(row.1) != length(row.5)){
            nrow <- min(length(row.1), length(row.5))
            row.1.1 <- vector(mode = "numeric", length = nrow)
            row.5.1 <- vector(mode = "numeric", length = nrow)
            for (k2 in 1:nrow){
              if (k2 == 1){
                row.1.1[k2] <- row.1[k2]
                row.5.1[k2] <- row.5[k2]
              } else {
                row.1.1[k2] <- first(row.1[row.1 > row.5.1[k2-1]])
                row.5.1[k2] <- first(row.5[row.5 > row.1.1[k2]])
              }
              
            }
            row.1 <- row.1.1[!is.na(row.1.1)]
            row.5 <- row.5.1[!is.na(row.5.1)]
          }
          
          # These are probably unnecessary?
          one.shift.eft[row.1, "Event"] <- 1
          one.shift.eft[row.5, "Event"] <- 5
          
          one.shift.eft[row.1, "Mother_Calf"] <- NA
          one.shift.eft[row.5, "Mother_Calf"] <- NA
          
          one.shift.eft <- one.shift.eft[!is.na(one.shift.eft$Effort),]
          # calculate effort for each "on" period per shift
          #tmp.eft <- 0
          idx.off <- which(one.shift.eft$Effort == "off")
          idx.on <- which(one.shift.eft$Effort == "on")
          
          if (length(idx.on) > 0){
            for (k2 in min(idx.on):nrow(one.shift.eft)){
              if (k2 == min(idx.on)){
                tmp.eft <- 0
                time.0 <- one.shift.eft[k2, "Minutes_since_0000"]
                effort <- "on"
              } else {
                if (one.shift.eft[k2, "Effort"] == "off" & effort == "on"){
                  d.time <- one.shift.eft[k2, "Minutes_since_0000"] - time.0
                  tmp.eft <- tmp.eft + d.time
                  effort <- "off"
                } else if (one.shift.eft[k2, "Effort"] == "on" & effort == "off"){
                  time.0 <- one.shift.eft[k2, "Minutes_since_0000"]
                  effort <- "on"
                } else if (one.shift.eft[k2, "Effort"] == "on" & effort == "on"){
                  d.time <- one.shift.eft[k2, "Minutes_since_0000"] - time.0
                  tmp.eft <- tmp.eft + d.time
                  effort <- "on"
                  time.0 <- one.shift.eft[k2, "Minutes_since_0000"]
                }
              }
              
            }
            
          } else {
            tmp.eft <- 0
          }
          
          # Sometimes, there is no Vis (or also sea state) update in an entire
          # shift. Use one from previous shift in those cases
          # These should not happen any longer but keep it anyways. 2023-05-10
          if (sum(!is.na(one.shift$SeaState.num)) > 0){
            max.sea.state <- max(one.shift$SeaState.num, na.rm = T)          
          } else {
            max.sea.state <- max.sea.state
          }
          
          if (sum(!is.na(one.shift$Vis.num)) > 0){
            max.Vis <- max(one.shift$Vis.num, na.rm = T) 
          } else {
            max.Vis <- max.Vis
          }
          
          out.list[[c]] <- data.frame(Date = all.dates[d],
                                      Minutes_since_0000 = one.shift.eft[1, "Minutes_since_0000"],
                                      Shift = k4,
                                      Effort = unname(tmp.eft),
                                      Mother_Calf = one.shift.eft %>%
                                        filter(Effort == "on") %>%
                                        summarize(MandC = sum(Mother_Calf, na.rm = TRUE)) %>%
                                        pull(),
                                      Sea_State = ifelse(!is.infinite(max.sea.state),
                                                         max.sea.state, NA),
                                      Vis = ifelse(!is.infinite(max.Vis),
                                                   max.Vis, NA))
          
          shift.list[[c]] <- one.shift.eft
          c <- c + 1
          
        }
      }  
    }
    
  }
  
  out.df <- do.call(rbind, out.list)
  shift.df <- do.call(rbind, shift.list)
  # out.df %>% 
  #   mutate(Date = as.Date(Date.char, 
  #                              format = "%Y-%m-%d")) -> out.df.1
  return(list(out.df = out.df,
              shift.df = shift.df))
}

# Removed MCMC.diag because they can be done easily with the posterior package

# Function to convert character time (e.g., 1349) to minutes from 
# a particular start time of the day (e.g., 0700). Default is midnight (0000)
# Use as char_time2min(x, origin = "0000"). 
char_time2min <- function(x, origin = "0000"){
  chars.origin <- strsplit(origin, split = "") %>% unlist
  if (length(chars.origin) == 3){
    hr <- as.numeric(chars.origin[1])
    m <- as.numeric(c(chars.origin[2:3]))
    M.0 <- hr * 60 + m[1] * 10 + m[2]
  } else {
    hr <- as.numeric(chars.origin[1:2])
    m <- as.numeric(c(chars.origin[3:4]))
    M.0 <- (hr[1] * 10 + hr[2]) * 60 + m[1] * 10 + m[2]
  }
  
  if (!is.na(x)){
    chars <- strsplit(x, split = "") %>% unlist()
    if (length(chars) == 3){
      hr <- as.numeric(chars[1])
      m <- as.numeric(c(chars[2:3]))
      M.1 <- (hr * 60 + m[1] * 10 + m[2])
    } else {
      hr <- as.numeric(chars[1:2])
      m <- as.numeric(c(chars[3:4]))
      M.1 <- (hr[1] * 10 + hr[2]) * 60 + m[1] * 10 + m[2]
    }
  } else {
    M.1 <- NA
  }
  
  if (!is.na(M.1)){
    out <- M.1 - M.0
    # if (M.1 >= M.0){
    #   out <- M.1 - M.0
    # } else {
    #   out <- 1440 - (M.0 - M.1)
    # }
    
  } else {
    out <- NA
  }
  
  return(out)
}

# converts minutes since 0000 in integer to time in character
# minute2time_char(420) will return "700"
minutes2time_char <- function(min_0000){
  H <- trunc(min_0000/60)
  M <- formatC(((min_0000/60 - H) * 60), width = 2, flag = "0")
  
  out <- paste0(H, M)
  return(out)
}

# change max.shift if needed. max.shift can be up to 5
find.shift <- function(x0){
  x <- as.numeric(x0)
  if (!is.na(x)){
    if (x < 600){
      shift <- "1"
    } else if (x == 600) {
      shift <- "1/2"
    } else if (x > 600 & x < 780) {
      shift <- "2"
    } else if (x == 780){
      shift <- "2/3"
    } else if (x > 780 & x < 960){
      shift <- "3"
    } else if (x == 960){
      shift <- "3/4"
    } else if (x > 960 & x < 1140){
      shift <- "4"
    } else if (x == 1140){
      shift <- "4/5"        
    } else if (x > 1140 & x <= 1320){ 
      shift <- "5"
    } else {
      shift <- NA
    }
  } else {
    shift <- NA
  }
  return(shift)
}

# A different way to look at shifts. Rather than starting at 0600 or 0700,
# start at 0100, so that occasional extra 30 minutes can be assigned to the 
# correct shift
find.shift.2 <- function(x0){
  x <- as.numeric(x0)
  if (!is.na(x)){
    if (x > 240 & x < 420){
      shift <- "2"
    } else if (x == 420) {
      shift <- "2/3"
    } else if (x > 420 & x < 600) {
      shift <- "3"
    } else if (x == 600){
      shift <- "3/4"
    } else if (x > 600 & x < 780){
      shift <- "4"
    } else if (x == 780){
      shift <- "4/5"
    } else if (x > 780 & x < 960){
      shift <- "5"
    } else if (x == 960){
      shift <- "5/6"        
    } else if (x > 960 & x <= 1140){ 
      shift <- "6"
    } else {
      shift <- NA
    }
  } else {
    shift <- NA
  }
  return(shift)
}


format.output <- function(data.shift, max.shift){
  # create a data frame with the full set of date/shift for dates with observations
  day.shifts <- data.frame(Date = rep(unique(data.shift$Date),
                                           each = max.shift),
                           Shift = rep(c(1,2,3,4), 
                                       times = length(unique(data.shift$Date))))
  
  # Combine the full set of date/shift with the observed - some shifts are NAs because
  # they were not in the dataset
  day.shifts %>% 
    left_join(data.shift, by = c("Date", "Shift")) -> all.day.shifts
  
  
  # create a vector with sequential weeks
  all.weeks <- seq.Date(as.Date(min(data.shift$Date)),
                        as.Date(max(data.shift$Date)),
                        by = "week")
  
  # select necessary data from per-shift data frame, mutate the column names
  # then create a new data frame
  data.shift %>% 
    select(Shift, Date, Effort, Mother_Calf) %>%
    mutate(Date = Date,
           Effort = Effort/60,
           Sightings = Mother_Calf) %>%
    select(Date, Shift, Effort, Sightings) %>%
    #mutate(Date = as.Date(Date)) %>%
    data.frame() -> raw.data
  
  # Create all dates, including weekends
  all.dates <- seq.Date(as.Date(min(data.shift$Date)),
                        as.Date(max(data.shift$Date)),
                        by = "day")
  
  # create all shifts, including nights
  all.shifts <- data.frame(Date = rep(all.dates,
                                      each = 8),
                           Shift = rep(c(1:8), 
                                       times = length(all.dates),
                                       by = "day"))
  
  all.shifts %>% 
    left_join(raw.data, by = c("Date", "Shift")) %>%
    mutate(Week = ceiling(difftime(Date, min(all.dates), 
                                   units = "weeks")) %>%
             as.numeric(),
           Date = Date) %>%
    select(Week, Date, Shift, Effort, Sightings) -> formatted.all.data
  
  # need to change week = 0 to week = 1...
  formatted.all.data[formatted.all.data$Week == 0, "Week"] <- 1
  
  # Change NAs to 0s
  formatted.all.data[is.na(formatted.all.data)] <- 0
  
  return(formatted.all.data)
}

format.output.2 <- function(data.shift){
  # create a data frame with the full set of date/shift for dates with observations
  # WITH NEW DEFINITIONS OF SHIFTS, THE FOLLOWING DATA FRAME NEEDS TO BE DEFINED
  # DIFFERENTLY. THE FIRST SHIFT STARTS AT 0100. SO MAX.SHIFT CANNOT BE USED TO 
  # REPLICATE DATE. 
  data.all.shifts <- data.frame(Date = rep(unique(data.shift$Date),
                                           each = 8),
                                Shift = rep(c(1:8), 
                                            times = length(unique(data.shift$Date))))
  
  # Combine the full set of date/shift with the observed - some shifts are NAs because
  # they were not in the dataset
  data.all.shifts %>% 
    left_join(data.shift, by = c("Date", "Shift")) -> all.data.shifts
  
  
  # create a vector with sequential weeks
  all.weeks <- seq.Date(as.Date(min(data.shift$Date)),
                        as.Date(max(data.shift$Date)),
                        by = "week")
  
  # select necessary data from per-shift data frame, mutate the column names
  # then create a new data frame
  all.data.shifts %>% 
    select(Shift, Date, Effort, Mother_Calf) %>%
    mutate(Date = Date,
           Effort = Effort/60,
           Sightings = Mother_Calf) %>%
    select(Date, Shift, Effort, Sightings) %>%
    #mutate(Date = as.Date(Date)) %>%
    data.frame() -> raw.data
  
  # Create all dates, including weekends
  all.dates <- seq.Date(as.Date(min(data.shift$Date)),
                        as.Date(max(data.shift$Date)),
                        by = "day")
  
  # create all shifts, including nights
  all.shifts <- data.frame(Date = rep(all.dates,
                                      each = 8),
                           Shift = rep(c(1:8), 
                                       times = length(all.dates),
                                       by = "day"))
  
  all.shifts %>% 
    left_join(raw.data, by = c("Date", "Shift")) %>%
    mutate(Week = ceiling(difftime(Date, min(all.dates), 
                                   units = "weeks")) %>%
             as.numeric(),
           Date = Date) %>%
    select(Week, Date, Shift, Effort, Sightings) -> formatted.all.data
  
  # need to change week = 0 to week = 1...
  formatted.all.data[formatted.all.data$Week == 0, "Week"] <- 1
  
  # Change NAs to 0s
  formatted.all.data[is.na(formatted.all.data)] <- 0
  
  return(formatted.all.data)
}



extract.all.vars <- function(x) {
  x %>%
    mutate(Date = Date.date) %>%
    select(-c(Date.date, Date.char)) %>%
    relocate(Date) -> x
  
  return(x)
}

extract.shift.vars <- function(x) {
  x %>%
    select(Date.char, Shift, Effort, Mother_Calf, Sea_State, Vis) %>%
    transmute(Date = Date.char,
              Shift = Shift,
              Sea_State = Sea_State,
              Vis = Vis,
              Effort = Effort, 
              Mother_Calf = Mother_Calf) -> x.out
  return(x.out)
}

# ver is the data extraction version. v2 is my original one
# v3 is one with the first shift starting at 0100.
file.names <- function(out.dir, ver, Year, out.list){
  out.file.name.inshore <- paste0(out.dir, "Processed_inshore_data_",
                              Year, "_", ver, ".csv")
  out.file.name.offshore <- paste0(out.dir, "Processed_offshore_data_",
                                  Year, "_", ver, ".csv")
  out.file.name.shift.inshore <- paste0(out.dir, "Processed_by_shift_inshore_data_", 
                                Year, "_", ver, ".csv")
  out.file.name.shift.offshore <- paste0(out.dir, "Processed_by_shift_offshore_data_", 
                                        Year, "_", ver, ".csv")
  
  out.file.name.formatted.inshore <- paste0(out.formatted.dir,
                                    Year, " Formatted_inshore_", ver, ".csv")
  out.file.name.formatted.offshore <- paste0(out.formatted.dir,
                                            Year, " Formatted_offshore_", 
                                            ver, ".csv")
  
  files <- list(inshore = out.file.name.inshore,
                offshore = out.file.name.offshore,
                shift.inshore = out.file.name.shift.inshore,
                shift.offshore = out.file.name.shift.offshore,
                formatted.inshore = out.file.name.formatted.inshore,
                formatted.offshore = out.file.name.formatted.offshore)
  
  #if (!file.exists(files$inshore))
    write_csv(out.list$data.inshore,
              file = files$inshore,
              quote = "none")
  
  # NEWER DATASETS DO NOT HAVE OFFSHORE LOGS (>2016). 
  if (!is.null(out.list$data.offshore))
    write_csv(out.list$data.offshore,
              file = files$offshore,
              quote = "none")
  
  #if (!file.exists(files$shift.inshore))
    write_csv(out.list$shift.data.inshore,
              file = files$shift.inshore, 
              quote = "none")
  
  if (!is.null(out.list$shift.data.offshore))
    write_csv(out.list$shift.data.offshore,
              file = files$shift.offshore, 
              quote = "none")
  
  #if (!file.exists(files$formatted.inshore))
    write_csv(out.list$formatted.data.inshore,
              file = files$formatted.inshore, 
              quote = "none",
              na = "")
  
  if (!is.null(out.list$formatted.data.offshore))
    write_csv(out.list$formatted.data.offshore,
              file = files$formatted.offshore, 
              quote = "none",
              na = "")
  
  return(files)
}


find.effort.dif <- function(Y, daily.summary.list, out.list){
  daily.summary.list[[which(years == Y)]]$daily.summary.1.2 %>%
    filter(abs(dif.effort) > 0.05 ) %>%
    dplyr::select(Date.char) %>%
    pull() -> date.dif.effort
  
  raw.data.all <- shift.dif <- data.2.dif <- data.1.dif <- list(length(date.dif.effort))
  
  k <- 13
  for (k in 1:length(date.dif.effort)){
    # For v1 extraction, extract the first four rows
    daily.summary.list[[which(years == Y)]]$data.1 %>%
      filter(Date == as.Date(date.dif.effort[k])) %>%
      slice(1:4) -> tmp 
    
    # Add one row at the top with 0 effort and 0 sightings
    # unelegant way of doing this but it works
    tmp <- rbind(tmp[1,], tmp) 
    tmp[1, "Effort"] <- 0
    tmp[1, "Sightings"] <- 0
    data.1.dif[[k]] <- tmp
    
    # for v3 extraction, extract rows 2-6, where rows 3-6 correspond
    # to v1 extraction.
    daily.summary.list[[which(years == Y)]]$data.2 %>%
      filter(Date == as.Date(date.dif.effort[k])) %>%
      slice(2:6) -> data.2.dif[[k]]
    
    # Absolute difference in effort is greater than 0.05 hr (3 min)
    # shift.dif is defined with the v3 shift index.
    shift.dif[[k]] <- data.2.dif[[k]]$Shift[which(abs(data.1.dif[[k]]$Effort - data.2.dif[[k]]$Effort) > 0.05)]
    
    if (length(shift.dif[[k]]) > 0){
      raw.data <- list(length(shift.dif[[k]]))
      k2 <- 1
      for (k2 in 1:length(shift.dif[[k]])){
        tmp.data <- out.list$shift.all.inshore %>% 
          filter(Date == as.Date(date.dif.effort[k]))
        #tmp.data <- out.list$data.inshore %>% 
        #  filter(Date == as.Date(date.dif.effort[k]))
        
        #NEED TO FILTER SO THAT THE CHUNK STARTS WITH EVENT = 1 AND ENDS WITH EVENT = 5
        tmp <- tmp.data[grep(shift.dif[[k]][k2], tmp.data$Shift),]
        
        if (nrow(tmp) == 0){
          raw.data[[k2]] <- NULL
        } else {
          raw.data[[k2]] <- tmp[which(tmp$Event == 1)[1]:max(which(tmp$Event == 5)),]          
        }

        
        #raw.data[[k2]] <- tmp.data[grep(shift.dif[[k]][k2], tmp.data$Shift),]
      }
      raw.data.all[[k]] <- raw.data
      
    }
    
  }
  
  return(list(date.dif = date.dif.effort,
              dif.1 = data.1.dif,
              dif.2 = data.2.dif,
              dif.shift = shift.dif,
              raw.data = raw.data.all))  
}

find.sightings.dif <- function(Y, daily.summary.list, out.list){
  daily.summary.list[[which(years == Y)]]$daily.summary.1.2 %>%
    filter(dif.sightings != 0) %>%
    dplyr::select(Date.char) %>%
    pull() -> date.dif.sightings
  
  raw.data.all <- shift.dif <- data.2.dif <- data.1.dif <- list(length(date.dif.sightings))
  
  k <- 1
  for (k in 1:length(date.dif.sightings)){
    daily.summary.list[[which(years == Y)]]$data.1 %>%
      filter(Date == as.Date(date.dif.sightings[k])) %>%
      slice(1:4) -> tmp 
    
    # Add one row at the top with 0 effort and 0 sightings
    # unelegant way of doing this but it works
    tmp <- rbind(tmp[1,], tmp) 
    tmp[1, "Effort"] <- 0
    tmp[1, "Sightings"] <- 0
    data.1.dif[[k]] <- tmp
    
    daily.summary.list[[which(years == Y)]]$data.2 %>%
      filter(Date == as.Date(date.dif.sightings[k])) %>%
      slice(2:6) -> data.2.dif[[k]]
    
    # Absolute difference in sightings is greater than 0
    shift.dif[[k]] <- data.2.dif[[k]]$Shift[which(abs(data.1.dif[[k]]$Sightings - data.2.dif[[k]]$Sightings) != 0)]
    
    if (length(shift.dif[[k]]) > 0){
      raw.data <- list(length(shift.dif[[k]]))
      k2 <- 1
      for (k2 in 1:length(shift.dif[[k]])){
        tmp.data <- out.list$shift.all.inshore %>% 
          filter(Date == as.Date(date.dif.sightings[k]))
        
        #NEED TO FILTER SO THAT THE CHUNK STARTS WITH EVENT = 1 AND ENDS WITH EVENT = 5
        tmp <- tmp.data[grep(shift.dif[[k]][k2], tmp.data$Shift),]
        
        if (nrow(tmp) == 0){
          raw.data[[k2]] <- NULL
        } else {
          raw.data[[k2]] <- tmp[which(tmp$Event == 1)[1]:max(which(tmp$Event == 5)),]          
        }
        
        # raw.data[[k2]] <- out.list$shift.all.inshore %>% 
        #   filter(Date == as.Date(date.dif.sightings[k])) %>%
        #   filter(Shift == shift.dif[[k]][k2])
        
      }
      raw.data.all[[k]] <- raw.data
      
    }
    
  }
  
  return(list(date.dif = date.dif.sightings,
              dif.1 = data.1.dif,
              dif.2 = data.2.dif,
              dif.shift = shift.dif,
              raw.data = raw.data.all))  
}





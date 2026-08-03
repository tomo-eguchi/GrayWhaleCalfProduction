
stop("Use GrayWhaleCalfProduction_fcns_v2.R instead.")

library(rstan)
library(loo)
library(readxl)


compute.LOOIC <- function(loglik.array, MCMC.params){
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  n.samples <- dim(loglik.array)[1]
  
  # Create an empty list to hold the log likelihood values
  #ll.list <- vector("list", length = n.samples)
  
  loglik.list <- lapply(1:n.samples, function(s) {
    # Extract the slice for the current sample 's'
    slice <- loglik.array[s, , , ]
    # Return the non-NA values from that slice
    slice[!is.na(slice)]
  })
  
  #loglik.list <- apply(loglik.array, 1, function(slice) slice[!is.na(slice)])
  
  # loop through each MCMC sample and collect log likelihood values
  # for (s in 1:n.samples){
  #   ll.cube <- loglik.array[s,,,]
  #   ll.list[[s]] <- ll.cube[!is.na(loglik.array[s,,,])]
  # }
  
  loglik.mat <- do.call(rbind, loglik.list)
  
  Reff <- relative_eff(exp(loglik.mat),
                       chain_id = rep(1:MCMC.params$n.chains,
                                      each = n.per.chain),
                       cores = MCMC.params$n.chains)
  
  loo.out <- rstanarm::loo(loglik.mat, 
                           r_eff = Reff, 
                           cores = MCMC.params$n.chains, 
                           k_threshold = 0.7)
  
  out.list <- list(Reff = Reff,
                   loo.out = loo.out)
  
  return(out.list)  
}


# Compute the "rank-normalized R-hat" by Vehtari et al. (2021) from jagsUI
# output.
# Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2021). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. Bayesian Analysis, 16(2), 667–718.
# https://doi.org/10.1214/20-BA1221
# 
# The first input is MCMC samples. If the jagsUI output is jm, this is jm$samples.
# The second input is a string of regular expression. This is a bit
# complicated. For example, to select "BF.Fixed" and all K parameters, which 
# are indexed, Use "^BF\\.Fixed|^K\\[" A '^' specifies that the following letter
# is the beginning of a string. '\\.' specifies a literal period, which needs
# to be "escaped" by two backslashes (\\). A square bracket needs to be escaped
# with two backslashes as well. The pipe (|) indicates 'or'.  
# 
rank.normalized.R.hat <- function(samples, params, MCMC.params){
  library(posterior)
  library(coda)
  
  col.names <- grep(params, varnames(samples), value = TRUE, perl = TRUE)
  subset.mcmc.samples <- samples[, col.names]
  
  subset.mcmc.array <- as_draws_array(subset.mcmc.samples, 
                                      .nchains = MCMC.params$n.chains)
  
  rhat.values <- apply(subset.mcmc.array, MARGIN = 3, FUN = posterior::rhat)
  
  return(rhat.values)
}


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


# extract MCMC diagnostic statistics, including Rhat, loglikelihood, DIC, and LOOIC
MCMC.diag <- function(jm, MCMC.params, params.to.monitor){
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  
  Rank.norm.Rhat <- rank.normalized.R.hat(jm$samples, 
                                          params = params.to.monitor)
  
  max.Rhat <- unlist(lapply(jm$Rhat, FUN = max, na.rm = T))
  
  loglik <- jm$sims.list$loglik
  
  Reff <- relative_eff(exp(loglik), 
                       chain_id = rep(1:MCMC.params$n.chains, 
                                      each = n.per.chain),
                       cores = MCMC.params$n.chains)
  loo.out <- loo(loglik, 
                 r_eff = Reff, 
                 cores = MCMC.params$n.chains)
  
  return(list(DIC = jm$DIC,
              loglik.obs = loglik,
              Reff = Reff,
              Rank.norm.Rhat = Rank.norm.Rhat,
              max.Rhat = max.Rhat,
              loo.out = loo.out))
  
}


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


# Excel file definitions
# Sheet names for inshore logs
list.sheet.names.inshore <- list(Y1994 = "INSHORE LOG", Y1995 = "Inshore Data", Y1996 = "data",
                                 Y1997 = "LOGS", Y1998 = "LOG", Y1999 = "DATA",
                                 Y2000 = "DATA", Y2001 = "LOGS", Y2002 = "LOGS",
                                 Y2003 = "Logs", Y2004 = "DATA", Y2005 = "LOGS",
                                 Y2006 = "LOGS", Y2007 = "log", Y2008 = "LOG",
                                 Y2009 = "LOGS", Y2010 = "LOGS", Y2011 = "LOGS",
                                 Y2012 = "LOGS", Y2013 = "LOGS", Y2014 = "LOGS",
                                 Y2015 = "LOGS", Y2016 = "LOGS", Y2017 = "LOGS",
                                 Y2018 = "LOGS", Y2019 = "INSHORE LOGS", Y2021 = "INSHORE LOGS",
                                 Y2022 = "MomCalf Logs 2022 postQC",
                                 Y2023 = "PiedrasBlancasCalfProductionDat",
                                 Y2024 = "MOTHER-CALF LOG 2024",
                                 Y2025 = "MOTHER-CALF LOG 2025")

# Sheet names for offshore logs. Starting 2017, offshore logs were eliminated.
list.sheet.names.offshore <- list(Y1994 = NULL, Y1995 = "Offshore Data", Y1996 = "Sheet1",
                                  Y1997 = "LOGS", Y1998 = "LOGS", Y1999 = "DATA",
                                  Y2000 = "LOGS", Y2001 = "LOGS", Y2002 = "LOG",
                                  Y2003 = "logs", Y2004 = "DATA", Y2005 = "LOGS",
                                  Y2006 = "LOGS", Y2007 = "log", Y2008 = "LOGS",
                                  Y2009 = "LOGS", Y2010 = "LOGS", Y2011 = "LOGS",
                                  Y2012 = "LOGS", Y2013 = "LOGS", Y2014 = "LOGS",
                                  Y2015 = "LOGS", Y2016 = NULL)

# Column types for inshore logs. Some years contain more columns (typically one) with
# no identifier, which return "New Name" when the file is read.
list.col.types.inshore <- list(Y1994 = c("date", "numeric", "text", "numeric", "text",
                                         "text", "numeric", "text", "text"),
                               Y1995 = c("date", "numeric", "text", "numeric", "text", 
                                         "text", "numeric", "numeric", "text", "text", 
                                         "text", "text", "text", "text", "text", "text"),
                               Y1996 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "numeric", "text", "text", 
                                         "text", "text", "text"),
                               Y1997 = c("date", "numeric", "text", "numeric", "numeric", 
                                         "numeric", "text", "text", "text", 
                                         "text", "numeric", "text", "text", 
                                         "text"),
                               Y1998 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "numeric", "text", "text", "text", 
                                         "text", "text"),
                               Y1999 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "numeric", "text", "text", "text",
                                         "text", "text", "text"),
                               Y2000 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text",
                                         "text", "text"),
                               Y2001 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text",
                                         "text", "text"),
                               Y2002 =  c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text", "text", "text",
                                          "text", "text"),
                               Y2003 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text",
                                         "text", "text"),
                               Y2004 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text",
                                         "text", "text"),
                               Y2005 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text", 
                                         "text", "text"),
                               Y2006 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text", 
                                         "text", "text", "text"),
                               Y2007 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text", 
                                         "text", "text", "text", "text"),
                               Y2008 = c("date", "numeric", "text", "numeric", "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text", "text", "text", "text"),
                               Y2009 = c("date", "numeric", "text", "numeric", "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text", "text", "text", "text"),
                               Y2010 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text", 
                                         "text", "text", "text"),
                               Y2011 =  c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text", "text", "text", 
                                          "text", "text", "text"),
                               Y2012 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text", 
                                         "text", "text", "text"),
                               Y2013 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text", 
                                         "text", "text", "text"),
                               Y2014 =  c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text", "text", "text", 
                                          "text", "text", "text"),
                               Y2015 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text", 
                                         "text", "text", "text"),
                               Y2016 = c("date", "numeric", "text", "numeric", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text", 
                                         "text", "text", "text"),
                               Y2017 = c("date", "numeric", "text", "text", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text", 
                                         "text", "text", "text", "text",
                                         "text", "text"),
                               Y2018 =  c("numeric", "date", "numeric", "text", "text", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text", "text", "text", 
                                          "text", "text", "text", "text",
                                          "text", "text", "text", "text"),
                               Y2019 =  c("date", "numeric", "numeric", "text", "text", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text", "text", "text", 
                                          "text", "text", "text", "text",
                                          "text", "text", "text"),
                               Y2021 = c("date", "numeric", "numeric", "text", "text", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "text", "text", "text", 
                                         "text", "text", "text", "text",
                                         "text", "text", "text", "text", "text"),
                               Y2022 = c("date", "numeric", "numeric", "text", "text", 
                                         "text", "text", "numeric", "numeric", "text", 
                                         "text", "text", "text", "numeric", "text", 
                                         "text", "text", "text"),
                               Y2023 = c("date", "numeric", "numeric", "text", "text",
                                         "text", "text","text", "numeric", "text", "text", 
                                         "text", "text", "text", "text", "text", 
                                         "text", "text"),
                               Y2024 = c("date", "numeric", "numeric", "text", "numeric",
                                         "numeric", "numeric","numeric", "numeric", "numeric", 
                                         "numeric", "text", "text", "text", "text", "numeric", 
                                         "text", "text", "text", "text"),
                               Y2025 = c("date", "numeric", "text", "numeric", "numeric",
                                         "numeric", "numeric", "text", "numeric", "numeric", "text"))

# Column types for offshore logs. Some years contain more columns (typically one) with
# no identifier, which return "New Name" when the file is read. No offshore files starting
# in 2017.
list.col.types.offshore <- list(Y1994 = c("date", "numeric", "text", "numeric", "text",
                                          "text", "numeric", "numeric", "text", "text"),
                                Y1995 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "numeric", "text", "numeric", "text"),
                                Y1996 =  c("date", "numeric", "text", "numeric", 
                                           "text", "text", "numeric", "numeric", 
                                           "numeric", "text",  "text"),
                                Y1997 = c("date", "numeric", "text", "numeric", 
                                          "numeric", "text", "text", "numeric", 
                                          "text", "text",  "numeric", "text"),
                                Y1998 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text",  "text", "text"),
                                Y1999 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text",  "text", "text"),
                                Y2000 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text",  "text"),
                                Y2001 =  c("date", "numeric", "text", "numeric", 
                                           "text", "text", "numeric", "numeric", 
                                           "text", "text",  "text", "text"),
                                Y2002 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text",  "text", "text"),
                                Y2003 =  c("date", "numeric", "text", "numeric", 
                                           "text", "text", "numeric", "numeric", 
                                           "text", "text",  "text", "text"),
                                Y2004 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text",  "text", "text"),
                                Y2005 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text",  "text", "text"),
                                Y2006 =  c("date", "numeric", "text", "numeric", 
                                           "text", "text", "numeric", "numeric", 
                                           "text", "text",  "text", "text"),
                                Y2007 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text",  "text", "text"),
                                Y2008 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text",  "text", "text"),
                                Y2009 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text",  "text", "text"),
                                Y2010 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text",  "text", "text"),
                                Y2011 =  c("date", "numeric", "text", "numeric", 
                                           "text", "text", "numeric", "numeric", 
                                           "text", "text",  "text", "text"),
                                Y2012 =  c("date", "numeric", "text", "numeric", 
                                           "text", "text", "numeric", "numeric", 
                                           "text", "text",  "text", "text"),
                                Y2013 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text",  "text", "text"),
                                Y2014 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text",  "text", "text"),
                                Y2015 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text",  "text", "text"),
                                Y2016 = c("date", "numeric", "text", "numeric", 
                                          "text", "text", "numeric", "numeric", 
                                          "text", "text",  "text", "text"))

# Column names that are extracted. These are not all the column names in each Excel spreadsheet.
# These are the column names from each Excel spreadsheet that were extracted to be used. 
# Columns are fixed - see get.data above. For recent years, columns need to be the following
# 1: Date
# 2: Event
# 3: Time
# 4: Observer
# 5: Sea state
# 6: Visibility
# 7: Mother/Calf count
# 8: Area
list.col.names.inshore <- list(Y1994 = c("Date", "Event", "Time", "Obs. Code", "Sea State",
                                         "Vis. IN", "Cow  / Calf"),
                               Y1995 = c("Date",	"Event Code",	"Time",	"Observer Code",	
                                         "Sea State",	"Vis In",	"C/C Pairs"),
                               Y1996 = c("DATE",	"EVENT CODE",	"TIME", "OBS.",	
                                         "SEA STATE",	"VIS (IN/OUT)",	"C/C"),
                               Y1997 = c("DATE",	"EVENT",	"TIME", "OBSER", 
                                         "SEA STATE", "VIS I/O", "C/C"),
                               Y1998 = c("Date",	"Event Code",	"Time",	"Obs",
                                         "Vis (I/O)", "Sea State",	"c/c"	),
                               Y1999 = c("Date",	"Event Code",	"Time", 	"Obs. Code",
                                         "Vis", "Sea State", "C/C"),
                               Y2000 = c("DATE",	"EVENT",	"TIME",	"OBS.",	"VIS (IN/OUT)",	"SEA STATE",	"C/C"),
                               Y2001 = c("DATE",	"EVENT",	"TIME",	"OBS.", "VIS I/O",	"SEA ST",	"C/C"	),
                               Y2002 = c("DATE",	"EVENT",	"TIME", 	"OBS", "VIS (I/O)",	"SS",	"C/C"),
                               Y2003 = c("DATE",	"EVENT",	"TIME", 	"OBS", "VIS (I/O)",	"SS",	"C/C"),
                               Y2004 = c("DATE",	"EVENT",	"TIME", 	"OBS", "VIS (I/O)",	"SS",	"C/C"),
                               Y2005 = c("DATE",	"EVENT",	"TIME", 	"OBS", "VIS (I/O)",	"SS",	"C/C"),
                               Y2006 =  c("DATE",	"EVENT", 	"TIME",	"OBS#", "VIS IN",	"SEA ST",	"C/C"),
                               Y2007 = c("DATE",	"EVENT", 	"TIME",	"OBS#", "VIS IN",	"SEA ST",	"C/C"),
                               Y2008 = c("DATE",	"EVENT", 	"TIME",	"OBS#", "VIS IN",	"SEA ST",	"C/C"),
                               Y2009 =  c("DATE",	"EVENT", 	"TIME",	"OBS#", "VIS IN",	"SEA ST",	"C/C"),
                               Y2010 = c("DATE",	"EVENT", 	"TIME",	"OBS#", "VIS IN",	"SEA ST",	"C/C"),
                               Y2011 = c("DATE",	"EVENT", 	"TIME",	"OBS#", "VIS IN",	"SEA ST",	"C/C"),
                               Y2012 = c("DATE",	"EVENT", 	"TIME",	"OBS#", "VIS IN",	"SEA ST",	"C/C"),
                               Y2013 = c("DATE",	"EVENT", 	"TIME",	"OBS#", "VIS IN",	"SEA ST",	"C/C"),
                               Y2014 = c("DATE",	"EVENT", 	"TIME",	"OBS#", "VIS IN",	"SEA ST",	"C/C", 
                                         "IN/OFF"),
                               Y2015 = c("DATE",	"EVENT", 	"TIME",	"OBS#", "VIS IN",	"SEA ST",	"C/C", 
                                         "IN/OFF"),
                               Y2016 =  c("DATE",	"EVENT", 	"TIME",	"OBS#", "VIS IN",	"SEA ST",	"C/C", 
                                          "IN/OFF"),
                               Y2017 = c("DATE",	"EVENT", 	"TIME",	"OBS#", "VIS IN",	"SEA ST",	"C/C", 
                                         "IN/OFF"),
                               Y2018 = c("DATE",	"EVENT", 	"TIME",	"OBS#", "VIS IN",	"SEA ST",	"C/C", 
                                         "IN/OFF"),
                               Y2019 = c("DATE",	"EVENT", 	"TIME",	"OBS #", "VIS IN",	"SEA ST",	"C/C", 
                                         "IN/OFF"),
                               Y2021 = c("DATE",	"EVENT", 	"TIME",	"OBS #", "VIS IN",	"SEA ST",	"C/C",
                                         "IN/OFF"),
                               Y2022 = c("DATE", 	"EVENT", 	"TIME",	"OBS #", "VIS IN", "SEA ST", "M/C",
                                         "LOCATION (IN/OFF)"),
                               Y2023 = c("DATE", "EVENT", "TIME", "OBS #", "SEA ST", "VIS IN", "M/C",
                                         "LOCATION (IN/OFF)"),
                               Y2024 = c("DATE", "EVENT", "TIME", "OBS #", "SEA STATE INSHORE", 
                                         "VIS INSHORE IN", "M/C", "LOCATION (IN/OFF)"),
                               Y2025 = c("Date", "Event Code", "Time", "Obs.1", "Sea State", "Vis",
                                         "M/C", "Location (Inshore/Offshore of PB Rocks)",
                                         "Time at Rocks", "Time Out", "Notes"))

# Same for offfshore spreadsheets
list.col.names.offshore <- list(Y1994 = c("Date", "Event", "Time", "Obs.", "Sea State",
                                          "Vis. Code", "Cow/Calf"),
                                Y1995 = c("Date",	"Event",	"Time",	"OBS.", 
                                          "Sea", "Vis.", 	"Cow/Calf"),
                                Y1996 = c("DATE",	"EVENT CODE",	"TIME", 	"OBS",
                                          "SEA ST.", "VIS",	"C/C"),
                                Y1997 = c("DATE",	"EVENT CODE",	"TIME", "OBS",
                                          "S STATE", "VIS", "COW/CA"),
                                Y1998 = c("Date",	"Event Code",	"Time",	"Obs",	
                                          "Vis",	"SeaState",	"c/c"),
                                Y1999 = c("Date",	"Event Code",	"Time",	"Observer",	
                                          "Vis. Code",	"Sea State",	"C/C"),
                                Y2000 = c("DATE",	"EVENT",	"TIME",	"OBS.",
                                          "VIS CODE",	"SEA STATE",	"C/C"),
                                Y2001 =  c("DATE",	"EVENT",	"TIME",	"OBS",
                                           "VIS",	"SEA STATE",	"C/C"),
                                Y2002 = c("DATE",	"EVENT",	"TIME",	"OBS", "VIS",	"SS",	"C/C"),
                                Y2003 = c("DATE",	"EVENT",	"TIME",	"OBS", "VIS",	"SS",	"C/C"),
                                Y2004 = c("DATE",	"EVENT",	"TIME",	"OBS", "VIS",	"SS",	"C/C"),
                                Y2005 = c("DATE",	"EVENT",	"TIME",	"OBS", "VIS",	"SS",	"C/C"),
                                Y2006 = c("DATE",	"EVENT",	"TIME",	"OBS", "VIS",	"SEA ST",	"C/C"),
                                Y2007 = c("DATE",	"EVENT",	"TIME",	"OBS", "VIS",	"SEA ST",	"C/C"),
                                Y2008 = c("DATE",	"EVENT",	"TIME",	"OBS", "VIS",	"SEA ST",	"C/C"),
                                Y2009 = c("DATE",	"EVENT",	"TIME",	"OBS", "VIS",	"SEA ST",	"C/C"),
                                Y2010 = c("DATE",	"EVENT",	"TIME",	"OBS", "VIS",	"SEA ST",	"C/C"),
                                Y2011 = c("DATE",	"EVENT",	"TIME",	"OBS", "VIS",	"SEA ST",	"C/C"),
                                Y2012 =  c("DATE", "EVENT",	"TIME",	"OBS", "VIS",	"SEA ST",	"C/C"),
                                Y2013 = c("DATE",	"EVENT",	"TIME",	"OBS", "VIS",	"SEA ST",	"C/C"),
                                Y2014 = c("DATE",	"EVENT",	"TIME",	"OBS", "VIS",	"SEA ST",	"C/C"),
                                Y2015 = c("DATE",	"EVENT",	"TIME",	"OBS", "VIS",	"SEA ST",	"C/C"),
                                Y2016 = c("DATE",	"EVENT",	"TIME",	"OBS",	
                                          "VIS",	"SEA ST",	"C/C"))




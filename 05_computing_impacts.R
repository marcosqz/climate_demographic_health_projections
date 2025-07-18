################################################################################

# This scripts estimates the health impacts of projected climate exposures by
# combining the latest climate projections with future demographic trends. 
# Specifically, it calculates heat-related mortality projection under the 
# SSP2-4.5 scenario for London. The analysis separates the contribution of
# demographic change and climate change, and quantifies the uncertainty
# associated with both sources.

#### LOAD LIBRARIES ############################################################

library(lubridate) # year
library(dlnm) # onebasis

#### LOAD DATA #################################################################

# Load model parameters
load("indata/processed/study_parameters.RData")

# Load output from the epi model
load("outdata/file/01_epi_model/arglag.RData")
load("outdata/file/01_epi_model/argvar.RData")
load("outdata/file/01_epi_model/dlnm_var.RData")
load("outdata/file/01_epi_model/coefsim_age.RData")
load("outdata/file/01_epi_model/mmt_age.RData")

# Load calibrated mortality and population projections
load(paste0(
  "outdata/file/02_calibrated_demographic_projections/",
  "data_proj_mort_popu_calibrated_ssp2.RData"))

# Load mortality and population observations
load("indata/processed/data_obs_temp_mort.RData")
load("indata/processed/data_obs_popu.RData")

# Load bias-corrected temperature projections
load(paste0(
  "outdata/file/03_calibrated_climate_projections/", 
  "data_proj_temp_biascorrection_ssp245.RData"))
load(paste0(
  "outdata/file/03_calibrated_climate_projections/",
  "data_proj_temp_constant_ssp245.RData"))

# Load warming levels
load("indata/processed/data_global_warming_levels_ssp245.RData")

#### EXPAND DEMOGRAPHIC PROJECTIONS TO DAILY VALUES ############################

# Repeat each row 5 times to create one row per year (within each 5-year block)
proj_mortpopu_expanded <- 
  proj_mortpopu_cal[rep(1:(nrow(proj_mortpopu_cal)), each = 5),]

# Replace repeated year values with consecutive years
# (1950, 1950, 1950, 1950, 1950, 1955 => 1950, 1951, 1952, 1953, 1954, 1955)
proj_mortpopu_expanded$year <- 
  c(sapply(proj_mortpopu_cal$year, function(y) {y + 0:4}))
rm(proj_mortpopu_cal)

# WE WANT TO TRANSALTE THE OBSERVED SEASONAL PATTERN OF THE MORTALITY TO
# THE MORTALITY PROJECTION SERIES

# Calculate the seasonality of the mortality in the observed data
weights_seas <- lapply(study_param$age_groups, function(i_age) { # loop age-groups
  
  # Average mortality over day of year (there are some with 365, other 366)
  seas_leap <- tapply(
    data_tempmort[[paste0("mort.", i_age)]], 
    as.numeric(format(data_tempmort$date, "%j")), mean)
  seas <- seas_leap[seq(365)] # this in a easy solution, method can be refined
  
  # Return the proportion of deaths for each day in leap and non-leap years
  return(list(
    weights_noleap = seas/sum(seas),
    weights_leap = seas_leap/sum(seas_leap)))
  
})
names(weights_seas) <- study_param$age_groups

# Put the weights in a dataset with the same days as in the demographics 
# projections
all_dates <- lapply(proj_mortpopu_expanded$year, function(i_year) { # loop projection years
  
  # Create a data frame with all dates
  dates <- data.frame(
    year = i_year,
    date = seq(as.Date(paste0(i_year, "-01-01")),
                as.Date(paste0(i_year, "-12-31")),
                by = "day"))
  
  # Add the mortality weights for each day in leap and non-leap years
  if(nrow(dates) == 365) {
    for(i_age in study_param$age_groups) { # loop age-groups
      dates[[paste0("weights.", i_age)]] <- 
        weights_seas[[i_age]]$weights_noleap
    }
  } else if(nrow(dates) == 366) {
    for(i_age in study_param$age_groups) { # loop age-groups
      dates[[paste0("weights.", i_age)]] <- 
        weights_seas[[i_age]]$weights_leap
    }
  }
  
  return(dates)
  
})
all_dates <- do.call(rbind, all_dates)
rm(weights_seas)

# Expand demographics projections from years to days
proj_mortpopu_daily <- merge(all_dates, proj_mortpopu_expanded)
rm(proj_mortpopu_expanded)

# Compute seasonal daily mortality by multiplying deaths by the weights
for(i_age in study_param$age_groups) { # loop age-groups
  
  # In each day multiply yealy value by the daily seasonal weight
  proj_mortpopu_daily[[paste0("mort.", i_age)]] <-
    proj_mortpopu_daily[[paste0("mort.", i_age)]] * 
    proj_mortpopu_daily[[paste0("weights.", i_age)]]
  
  # Remove the column with the weight
  proj_mortpopu_daily[[paste0("weights.", i_age)]] <- NULL
}; rm(i_age)

#### CREATE CONSTANT DEMOGRAPHIC PROJECTIONS ###################################

# EXTRACT THE MEAN OF THE LAST 5 YEARS OF OBSERVATIONS (SIMILAR TO THE CONSTANT
# CLIMATE PROJECTIONS CALCULATED BEFORE)

# Create dataset with yearly mortality observations
data_mort <- data_tempmort[, c("year", "mort.00_74", "mort.75plus")]
rm(data_tempmort) 
data_mort <- aggregate(. ~ year, data = data_mort, FUN = sum)
# Remove row for 2012 (we don't have mortality data for the whole year)
data_mort <- subset(data_mort, year != 2012)

# Calculate the mean of last 5 years
data_mort_baseline <- colMeans(tail(data_mort, 5))
data_popu_baseline <- colMeans(tail(data_popu, 5))
rm(data_mort, data_popu)

# Maintain the data.frame structure
data_mort_baseline <- data.frame(t(data_mort_baseline))
data_popu_baseline <- data.frame(t(data_popu_baseline))

# Repeat the mean of last 5 years for the whole projection period
years_proj <-  unique(all_dates$year)
data_mort_baseline <- data_mort_baseline[rep(1, length(years_proj)),]
data_popu_baseline <- data_popu_baseline[rep(1, length(years_proj)),]
data_mort_baseline$year <- years_proj
data_popu_baseline$year <- years_proj
rm(years_proj)

# Expand data_mort_baseline and popu from year to daily
data_mort_baseline <- merge(all_dates, data_mort_baseline)
data_popu_baseline <- merge(all_dates[, c("year", "date")], data_popu_baseline)
rm(all_dates)

# Compute seasonal daily mortality by multiplying deaths by the weights
for(i_age in study_param$age_groups) {
  
  # In each day multiply yealy value by the daily seasonal weight
  data_mort_baseline[[paste0("mort.", i_age)]] <- 
    data_mort_baseline[[paste0("mort.", i_age)]] * 
    data_mort_baseline[[paste0("weights.", i_age)]]
  
  # Remove the column with the weight
  data_mort_baseline[[paste0("weights.", i_age)]] <- NULL
}; rm(i_age)

# Bind constant mortality and population
proj_mortpopu_constant <- merge(data_mort_baseline, data_popu_baseline)
rm(data_mort_baseline, data_popu_baseline)

#### PREPARE TEMPERATURE PROJECTION DATASETS ###################################

# Create a function to transform temperature projections from wide to long
transform_tolong_tempproj <- function(x) {
  proj_temp_long <- reshape(
    x,
    varying = list(colnames(x[, !(names(x) %in% "date")])),  # all columns except 'date'
    v.names = "temperature",
    timevar = "gcm",
    times = colnames(x[, !(names(x) %in% "date")]),
    idvar = "date",
    direction = "long")
  
  proj_temp_long <- proj_temp_long[order(proj_temp_long$date), ]
  proj_temp_long$gcm <- gsub("temp\\.", "", proj_temp_long$gcm)
  rownames(proj_temp_long) <- NULL
  proj_temp_long
}

# Temperature projections from wide to long
proj_temp_bc <- transform_tolong_tempproj(proj_temp_bc)
proj_temp_constant <- transform_tolong_tempproj(proj_temp_constant)

#### CALCULATE HEAT-RELATED IMPACTS ############################################

# Build a function for the calculation of the impacts
calculate_health_impacts <- function(constant_temp,
                                     constant_mortpopu) {
  
  # Select between original or constant climate projections
  if(!constant_temp){
    f.temp <- proj_temp_bc
  } else {
    f.temp <- proj_temp_constant
  }
  
  # Select between original or constant demographic projections
  if(!constant_mortpopu){
    f.mortpopu <- proj_mortpopu_daily
  } else {
    f.mortpopu <- proj_mortpopu_constant
  }
  
  # Loop age-groups to compute heat-related mortality
  impacts_age <- lapply(study_param$age_groups, function(i_age){
    
    print(paste0("Calculation heat-related mortality: ", i_age))
    
    # CALCULATE ATTRIBUTABLE FRACTIONS FOR EACH DAY, GCM AND SIMULATION
    
    # Exposure-response basis at the age-specific mmt
    cenvec <- onebasis(x = mmt_age[i_age], 
                       fun = argvar$fun, 
                       knots = argvar$knots,
                       Boundary.knots = argvar$Bound)
    
    # Exposure-response basis at the projected temperatures centered at the
    # age-specific mmt
    bcen <- scale(onebasis(x = f.temp$temperature,
                           fun = argvar$fun, 
                           knots = argvar$knots,
                           Boundary.knots = argvar$Bound),
                  center = cenvec, scale = FALSE)
    
    # Relative risks and attributable fractions calculation
    rrsim <- exp(bcen %*% coefsim_age[[i_age]])
    afsim <- (rrsim - 1) / rrsim # AF = (RR-1)/RR
    
    # Calculate heat-related deaths by setting AF=0 to temperatures below mmt
    ind_heat <- f.temp$temperature > mmt_age[i_age]
    afsim[!ind_heat,] <- 0
    
    # TRANSFORM AF IN A LONG DATA-FRAME TO ENABLE EASY CALCULATION OF OTHER
    # IMPACT MEASUERES
    
    # Each column as a simulation of the epidemiological model
    colnames(afsim) <- paste0("sim", 1:study_param$n_sim)
    impacts <- data.frame(afsim); rm(afsim)
    
    # Add columns with date and gcm (datasets need to be ordered)
    impacts$date <- f.temp$date
    impacts$gcm <- f.temp$gcm
    
    # Transform afsim from wide to long 
    # (this loop works slightly faster than reshape function)
    impacts_long <- lapply(1:study_param$n_sim, function(i) { # loop simulations
      
      impacts_subset <- impacts[,c("date", "gcm", paste0("sim", i))]
      
      # Add column with the epidemiological simulation number
      impacts_subset$simulation <- i
      
      # Colum paste0("sim", i) stores the AF
      colnames(impacts_subset)[colnames(impacts_subset) == paste0("sim", i)] <-
        paste0("af.", i_age)
      
      return(impacts_subset)
      
    })
    impacts_long <- do.call(rbind, impacts_long)
    
    # MERGE IMPACTS DATAFRAME WITH DEMOGRAPHIC PROJECTIONS
    impacts_long <- merge(
      impacts_long,
      f.mortpopu[, c("year", "date", 
                     paste0("mort.", i_age), paste0("popu.", i_age))],
      by = "date")
    
    # CALCULATE ATTRIBUTABLE NUMBER (AN = AF * MORT)
    impacts_long[[paste0("an.", i_age)]] <- # AN
      impacts_long[[paste0("af.", i_age)]] * # AF
      impacts_long[[paste0("mort.", i_age)]] # MORT 

    # CALCULATE ATTRIBUTABLE NUMBER RATE (AN_RATE = AN / POPULATION * 100000)
    impacts_long[[paste0("an_rate.", i_age)]] <- # AN RATE
      impacts_long[[paste0("an.", i_age)]] / # AN
      impacts_long[[paste0("popu.", i_age)]] * 100000 # POPULATION * 100000
    
    return(impacts_long[, c("year", "date", "gcm", "simulation", 
                            paste0(c("af.", "an.", "an_rate."), i_age))])
    
  })
  # Bind by columns the age-specific impact datasets
  impacts <- do.call(cbind, impacts_age)

  # Remove the duplicated rows (date, gcm, sim) created by binding datasets
  impacts <- impacts[, !duplicated(as.list(impacts))]
  
  # CALCULATE TOTAL AN BY SUMMING AGE-SPECIFIC ANs
  impacts$an <- rowSums(impacts[, paste0("an.", study_param$age_groups)])
  
  return(impacts)
  
}

# CALCULATE THE IMPACTS FOR PROJECTIONS (FULL, CLIMATE AND DEMOGRAPHIC 
# CONTRIBUTIONS)
# Combined climate and demographic contributions
impacts_full <- calculate_health_impacts(
  constant_temp = FALSE, 
  constant_mortpopu = FALSE)
# Demographic contribution (temperatures are held constant)
impacts_exclude_temp <- calculate_health_impacts(
  constant_temp = TRUE, 
  constant_mortpopu = FALSE)
# Climate contribution (demographics are held constant)
impacts_exclude_demo <- calculate_health_impacts(
  constant_temp = FALSE, 
  constant_mortpopu = TRUE)

#### TEMPORAL AGGREGATION OF IMPACTS (DAILY TO YEARLY) #########################

# Create function to aggregate impact datasets from daily to yearly
transform_impacts_to_yearly <- function(impacts) {
  
  # Subset to relevant columns
  impacts_year <- 
    impacts[ , c("year", "gcm", "simulation", "an", 
                 paste0("an.", study_param$age_groups))]
  
  # Response variables to be aggregated (total and age-specific AN)
  response_vars <- c("an", paste0("an.", study_param$age_groups))
  
  # Construct the formula to aggregate response by year, gcm, and simulation
  agg_formula <- paste("cbind(", paste(response_vars, collapse = ", "), 
                       ") ~ year + gcm + simulation")
  agg_formula <- as.formula(agg_formula)
  
  # Aggregate 
  impacts_year <- aggregate(agg_formula, data = impacts_year, FUN = sum)
  
  return(impacts_year)
  
}

# TEMPORAL AGGREGATION TO YEARLY IMPACTS
an_year_full <- 
  transform_impacts_to_yearly(impacts = impacts_full)
an_year_exclude_temp <- 
  transform_impacts_to_yearly(impacts = impacts_exclude_temp)
an_year_exclude_demo <- 
  transform_impacts_to_yearly(impacts = impacts_exclude_demo)

#### TEMPORAL AGGREGATION OF IMPACTS (BY PERIOD AND GWLs) ######################

# Create function to aggregate yearly impact datasets to longer periods
compute_an_warming <- function(impacts_year){
  
  an_warming <- lapply(c("gwl", "end_century"), function(method) {
    
    # Loop for each gcm (GLW periods depend on the GCM)
    an_warming <- lapply(study_param$selected_gcms, function(i_gcm) {
      
      if(method == "gwl") {
        # Extract the start and end of the 21-year period for that specific gcm 
        # and level of warming
        years_warming <- subset(
          data_gwl, 
          (gcm == i_gcm) & (warming_level == study_param$selected_warming),
          select = c("year1", "year2"))
        years_warming <- unlist(years_warming)
      } else if(method == "end_century") {
        years_warming <- c(2079, 2099)
      }
      
      # Subset for the gcm and the years defining the period
      subset_data <- subset(
        impacts_year, 
        (gcm == i_gcm) & (year %in% years_warming[1]:years_warming[2]),
        select = c(simulation, an))
      
      # Sum AN for the 21-year period
      subset_data <- aggregate(an ~ simulation, subset_data, FUN = sum)
      subset_data$gcm <- i_gcm
      subset_data[, c("gcm", "simulation", "an")]
      
      return(subset_data)
      
    })
    an_warming <- do.call(rbind, an_warming)
    
    return(an_warming)
    
  })
  names(an_warming) <- c("gwl", "end_century")
  
  return(an_warming)
  
}

# TEMPORAL AGGREGATION TO GWL (2C)
an_period_full <- compute_an_warming(an_year_full)
an_period_exclude_temp <- compute_an_warming(an_year_exclude_temp)
an_period_exclude_demo <- compute_an_warming(an_year_exclude_demo)

#### SAVE OUTPUTS ##############################################################

# Save heat-related mortality by year
save(an_year_full, file = paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_year_full_",
  study_param$ssp_rcp_scenario,".RData"))
save(an_year_exclude_temp, file = paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_year_exclude_temp_",
  study_param$ssp_rcp_scenario,".RData"))
save(an_year_exclude_demo, file = paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_year_exclude_demo_",
  study_param$ssp_rcp_scenario,".RData"))

# Save heat-related mortality by 21-year periods
save(an_period_full, file = paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_period_full_",
  study_param$ssp_rcp_scenario,".RData"))
save(an_period_exclude_temp, file = paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_period_exclude_temp_",
  study_param$ssp_rcp_scenario,".RData"))
save(an_period_exclude_demo, file = paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_period_exclude_demo_",
  study_param$ssp_rcp_scenario,".RData"))
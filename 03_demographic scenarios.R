################################################################################

# Calibration of national-level SSP2 mortality and population projection for the
# United Kingdom and Northern Ireland to adjust to London's absolute numbers.
# This scripts computes city-specific correction factors by comparing
# national-level projections with observed values for London during overlapping
# historical periods. The entire projection series is the scaled by the ratio of
# observed to projected means in the baseline period.

#### LOAD LIBRARIES ############################################################

#### LOAD DATA #################################################################

load("indata/processed/data_obs_temp_mort.RData")
load("indata/processed/data_obs_popu.RData")
load("indata/processed/data_proj_mort_popu_ssp2.RData")
load("indata/processed/study_parameters.RData")

#### CALIBRATION OF DEMOGRAPHIC PROJECTIONS ####################################

# TRANSFORM DAILY MORTALITY OBSERVATIONS INTO YEARLY MORTALITY OBSERVATIONS

# Extract mortality from the temperature-mortality dataset
data_mort <- subset(
  data_tempmort, 
  select = c("year",  paste0("mort.", study_param$age_groups)))

# Create formula for aggregation: cbind(mort.00_74, mort.75plus) ~ year
formula_agg <- as.formula(paste0(
  "cbind(", 
  paste(paste0("mort.", study_param$age_groups), collapse = ", "),
  ") ~ year"))

# Aggregate by year
data_mort <- aggregate(formula_agg, data = data_mort, FUN = sum)
  
# Remove row for 2012 (we don't have mortality data for the whole year)
data_mort <- subset(data_mort, year != 2012)

# LOOP MORTALITY AND POPULATION CALIBRATION (IT IS THE SAME PROCEDURE)
proj_mortpopu_cal <- lapply(c("mort", "popu"), function(var) {
  
  # SELECT MORT/POPU OBSERVATIONS DATASET
  if(var == "mort") {
    data <- data_mort
  } else if (var == "popu") {
    data <- data_popu
  }
  
  # SUBSET PROJECTIONS DATASET WITH THE CORRESPONDING VARIABLE AND AGE GROUP
  proj <- subset(
    proj_mortpopu,
    select = c("year", paste0(var, ".",  study_param$age_groups)))
  
  # MERGE OBSERVATIONS AND PROJECTIONS
  merged_data <- merge(data, proj, by = "year")
  
  # CALCULATE MEAN OF OBSERVATIONS AND PROJECTIONS IN THE COINCIDING YEARS
  merged_data <- subset(merged_data, select = -year)
  merged_data <- colMeans(merged_data)
  
  # LOOP AGE GROUPS
  for(i_age in study_param$age_groups) {
    
    # CORRECTION AS THE RATIO BETWEEN COINCIDING OBSERVATIONS AND PROJECTIONS
    correction <- 
      merged_data[[paste0(var, ".", i_age, ".x")]] / # Observations (London)
      merged_data[[paste0(var, ".", i_age, ".y")]]   # Projections (UK and NI)
    
    # APPLY THE CORRECTION TO ALL THE PROJECTIONS
    proj[[paste0(var, ".", i_age)]] <- 
      proj[[paste0(var, ".", i_age)]] * correction
    
    rm(correction)
    
  }; rm(i_age)
  
  return(proj)
  
})

# MERGE CALIBRATED PROJECTIONS OF MORTALITY AND POPULATION
proj_mortpopu_cal <- do.call(merge, proj_mortpopu_cal)

#### TEMPORAL OF DEMOGRAPHIC PROJECTIONS #######################################

# EXPAND DEMOGRAPHIC PROJECTIONS TO DAILY VALUES

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
data_mort <- subset(
  data_tempmort,
  select = c("year", paste0("mort.", study_param$age_groups)))
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
  
  # Remove the column with the mortality weight
  data_mort_baseline[[paste0("weights.", i_age)]] <- NULL
}; rm(i_age)

# Bind constant mortality and population
proj_mortpopu_constant <- merge(data_mort_baseline, data_popu_baseline)
rm(data_mort_baseline, data_popu_baseline)

#### SAVE OUTPUTS ##############################################################

save(proj_mortpopu_daily, file = paste0(
  "outdata/file/02_calibrated_demographic_projections/",
  "data_proj_mort_popu_calibrated_daily_ssp", 
  study_param$ssp_scenario,
  ".RData"))
save(proj_mortpopu_constant, file = paste0(
  "outdata/file/02_calibrated_demographic_projections/",
  "data_proj_mort_popu_calibrated_constant_ssp", 
  study_param$ssp_scenario,
  ".RData"))
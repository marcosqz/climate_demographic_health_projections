################################################################################

# This script performs spatial and temporal calibration of national-level 
# quinquennial SSP2 mortality projections to generate city-specific and daily 
# calibrated series. It then computes the attributable number (AN) of 
# heat-related deaths by combining the calibrated demographic projections with 
# the climate-driven attributable fractions (AF) obtained in the previous steps.

rm(list = ls())

#### LOAD DATA #################################################################

load("indata/processed/data_obs_temp_mort.RData")
load("indata/processed/data_obs_popu.RData")
load("indata/processed/data_proj_mort_ssp2.RData")
load("indata/processed/study_parameters.RData")

# Load attributable fractions
load(paste0(
  "outdata/file/02_calibrated_climate_projections/attributable_fraction_",
  study_param$ssp_rcp_scenario,".RData"))

#### CALIBRATION OF DEMOGRAPHIC PROJECTIONS ####################################

# ---- SPATIAL CALIBRATION OF MORTALITY PROJECTIONS ----
# Adjust national-level mortality projections (UK & NI) to reflect 
# city-specific (London) demographic levels using observed data.

# Extract mortality from the temperature-mortality dataset
data_mort <- subset(
  data_tempmort, 
  select = c("year",  paste0("mort.", study_param$age_groups)))

# Create formula for aggregate mortality over years
formula_agg <- as.formula(paste0(
  "cbind(", 
  paste(paste0("mort.", study_param$age_groups), collapse = ", "),
  ") ~ year"))

# Aggregate observed mortality by year
data_mort <- aggregate(formula_agg, data = data_mort, FUN = sum)
  
# Remove row for 2012 (we don't have mortality data for the whole year)
data_mort <- subset(data_mort, year != 2012)
  
# Merge observation and projections (it only keeps matching years)
merged_data <- merge(data_mort, proj_mort, by = "year")
  
# Calculate mean of observations and projections in the coinciding years
merged_data <- subset(merged_data, select = -year)
merged_data <- colMeans(merged_data)

# Calculate correction factors
corr_factor <- sapply(study_param$age_groups, function(i_age) {
  merged_data[[paste0("mort.", i_age, ".x")]] / # Observations (London)
    merged_data[[paste0("mort.", i_age, ".y")]] # Projections (UK and NI)
})

# Create a new dataset with the spatially calibrated projection
proj_mort_ldn <- proj_mort
for(i_age in study_param$age_groups) {
  proj_mort_ldn[[paste0("mort.", i_age)]] <-
    proj_mort_ldn[[paste0("mort.", i_age)]] * corr_factor[[i_age]]
}; rm(i_age)

# ---- TEMPORAL CALIBRATION OF DEMOGRAPHIC PROJECTIONS ----
# Disaggregate annual mortality projections into daily values 
# by applying the observed seasonal pattern of mortality.

# Calculate the seasonal pattern of within-year, age-specific observed mortality
weights_seas_doy <- lapply(study_param$age_groups, function(i_age) {
  
  # Calculate average counts for each day of the year (doy) from observed series 
  deathdoy <- tapply(
    data_tempmort[[paste0("mort.", i_age)]], 
    as.numeric(format(data_tempmort$date, "%j")), 
    mean, na.rm = TRUE)[seq(365)]
  
  # Compute relative weights by dividing each dow’s mean mortality by the 
  # annual total
  weights_seas <- deathdoy/sum(deathdoy)
  
  return(weights_seas)
  
}); names(weights_seas_doy) <- study_param$age_groups

# Create a full-period dataset of daily seasonal weights for each
# projection year
weights_seas_period <- lapply(proj_mort$year, function(i_year) { # loop projection years
  
  # Create a data frame with all dates
  dates <- data.frame(
    year = i_year,
    date = seq(as.Date(paste0(i_year, "-01-01")),
               as.Date(paste0(i_year, "-12-31")),
               by = "day"))

  # Remove leap days
  if(lubridate::leap_year(i_year)){
    dates <- subset(dates, date != as.Date(paste0(i_year, "-02-29")))
  }
  
  # Add the mortality weights for each day to the dates dataset
  for(i_age in study_param$age_groups) { # loop age-groups
    dates[[paste0("weights.", i_age)]] <- weights_seas_doy[[i_age]]
  }
  
  return(dates)
  
})
weights_seas_period <- do.call(rbind, weights_seas_period)
rm(weights_seas_doy)

# Expand demographics projections from years to days
proj_mort_ldn_daily <- merge(weights_seas_period, proj_mort_ldn)

# Compute seasonal daily mortality by multiplying yearly deaths by the weights
# of the day of the year
for(i_age in study_param$age_groups) { # loop age-groups
  
  # Multiply annual mortality totals by daily weights to obtain daily series
  proj_mort_ldn_daily[[paste0("mort.", i_age)]] <-
    proj_mort_ldn_daily[[paste0("mort.", i_age)]] * 
    proj_mort_ldn_daily[[paste0("weights.", i_age)]]
  
  # Remove the column with the weights
  proj_mort_ldn_daily[[paste0("weights.", i_age)]] <- NULL
  
}; rm(i_age)

#### ESTIMATION OF THE ATTRIBUTABLE NUMBER OF DEATHS ###########################
# Compute daily heat-related deaths (AN) by multiplying the attributable 
# fraction (AF) by the calibrated daily mortality for each age group and GCM.

# Loop age groups
an <- lapply(study_param$age_groups, function(i_age) {
  
  # Loop GCMs
  an <- lapply(study_param$selected_gcms, function(i_gcm) {
    
    # AN = AF * mort
    an <- af[[i_age]][[i_gcm]] * # Attributable fraction
      proj_mort_ldn_daily[[paste0("mort.", i_age)]] # Mortality
    
    return(an)
    
  }); names(an) <- study_param$selected_gcms
  
  return(an)
  
}); names(an) <- study_param$age_groups

#### SAVE OUTPUTS ##############################################################

save(corr_factor, file = paste0(
  "outdata/file/03_calibrated_demographic_projections/",
  "correction_factor_spatialcal_ssp", 
  study_param$ssp_scenario,
  ".RData"))
save(proj_mort_ldn, file = paste0(
  "outdata/file/03_calibrated_demographic_projections/",
  "data_proj_mort_spatialcal_ssp", 
  study_param$ssp_scenario,
  ".RData"))
save(proj_mort_ldn_daily, file = paste0(
  "outdata/file/03_calibrated_demographic_projections/",
  "data_proj_mort_spatialcal_tempcal_ssp", 
  study_param$ssp_scenario,
  ".RData"))
save(an, file = paste0(
  "outdata/file/03_calibrated_demographic_projections/attributable_number_", 
  study_param$ssp_rcp_scenario, 
  ".RData"))
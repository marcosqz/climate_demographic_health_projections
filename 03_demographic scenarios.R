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

# ARRANGE MORTALITY DATASET AND AGGREGATE BY YEAR (SAME AS PROJECTIONS)
data_mort <- data_tempmort[, c("year", 
                               paste0("mort.", study_param$age_groups))]
rm(data_tempmort)

# Create formula for aggregation: cbind(mort.00_74, mort.75plus) ~ year
formula_agg <- as.formula(paste0(
  "cbind(", 
  paste(paste0("mort.", study_param$age_groups), collapse = ", "),
  ") ~ year"))

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
  
  # PROJECTIONS DATASET
  proj <- proj_mortpopu[, c("year", paste0(var, ".",  study_param$age_groups))] 
  
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

#### SAVE OUTPUTS ##############################################################

save(proj_mortpopu_cal, file = paste0(
  "outdata/file/02_calibrated_demographic_projections/",
  "data_proj_mort_popu_calibrated_ssp", study_param$ssp_scenario,".RData"))
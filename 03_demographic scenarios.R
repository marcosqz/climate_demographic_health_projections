#### 3.2. CALIBRATION OF MORTALITY AND POPULATION PROJECTIONS ####

# Load mortality observations
load("../../outdata/file/data_tempmort.RData")
# Load population observations
load("../../outdata/file/data_popu.RData")
# Load mortality and population projections
load("../../outdata/file/data_projection_mortality_population_ssp2.RData")

# Keep only the mortality columns on the observations
data_mort <- data_tempmort[, c("year", "mort.00_74", "mort.75plus")]
rm(data_tempmort)

# Aggregate mortality observations by years (same as mortality projections)
data_mort <- aggregate(cbind(mort.00_74, mort.75plus) ~ year, 
                       data = data_mort, FUN = sum) # TODO: 2012 is not complete, choose between 1) account for that or 2) use only complete years

# Calibrate projections to respect to observations for mortalities and temperatures
proj_mortpopu_cal <- lapply(c("mort", "popu"), function(var) {
  
  if(var == "mort") {
    data <- data_mort
  } else if (var == "popu") {
    data <- data_popu
  }
  
  proj <- proj_mortpopu[, c("year", paste0(var, ".", c("00_74", "75plus")))] # Keep columns for morttalities or populations
  
  # MERGE OBSERVATIONS AND PROJECTIONS
  merged_data <- merge(data, proj, by = "year")
  
  # CALCULATE MEAN OF OBSERVATIONS AND PROJECTIONS IN THE COINCIDING YEARS
  merged_data <- subset(merged_data, select = -year)
  merged_data <- colMeans(merged_data)
  
  for(i_age in c("00_74", "75plus")) { # LOOP AGE GROUPS
    
    # CORRECTION AS THE RATIO BETWEEN OBSERVATIONS AND OBSVATION
    correction <- 
      merged_data[[paste0(var, ".", i_age, ".x")]] / # Observations (London)
      merged_data[[paste0(var, ".", i_age, ".y")]]   # Projections (UK and NI)
    
    # APPLY THE CORRECTION FOR THE CORRESPONDING AGE GROUP
    proj[[paste0(var, ".", i_age)]] <- 
      proj[[paste0(var, ".", i_age)]] * correction
    
    rm(correction)
    
  }; rm(i_age)
  
  proj # return the calibrated projections
  
})
# Merge calibrated projections of mortality and population
proj_mortpopu_cal <- do.call(merge, proj_mortpopu_cal) 

# Save dataset
save(proj_mortpopu_cal, 
     file = "../../outdata/file/data_projection_temperature_calibrated_ssp2rcp45.RData")
#### TODO: ADD DESCRIPTION OF THIS SCRIPT

# Calibration of mortality and population projection (SS2) for United Kingdom
# and Northern Ireland to respect to the London observations

#### LOAD LIBRARIES ############################################################

#### LOAD DATA #################################################################

load("outdata/file/data_tempmort.RData")
load("outdata/file/data_popu.RData")
load("outdata/file/data_projection_mortality_population_ssp2.RData") # TODO: change ssp2

#### CALIBRATION OF DEMOGRAPHIC PROJECTIONS ####################################

# ARRANGE MORTALITY DATASET AND AGGREGATE BY YEAR (SAME AS PROJECTIONS)
data_mort <- data_tempmort[, c("year", "mort.00_74", "mort.75plus")]
rm(data_tempmort)
data_mort <- aggregate(cbind(mort.00_74, mort.75plus) ~ year, 
                       data = data_mort, FUN = sum)
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
  proj <- proj_mortpopu[, c("year", paste0(var, ".", c("00_74", "75plus")))] 
  
  # MERGE OBSERVATIONS AND PROJECTIONS
  merged_data <- merge(data, proj, by = "year")
  
  # CALCULATE MEAN OF OBSERVATIONS AND PROJECTIONS IN THE COINCIDING YEARS
  merged_data <- subset(merged_data, select = -year)
  merged_data <- colMeans(merged_data)
  
  # LOOP AGE GROUPS
  for(i_age in c("00_74", "75plus")) {
    
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

save(proj_mortpopu_cal, 
     file = "outdata/file/data_projection_mortpopu_calibrated_ssp2rcp45.RData") # TODO: change ssp2rcp45
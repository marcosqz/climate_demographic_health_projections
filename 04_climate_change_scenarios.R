#### LOAD LIBRARIES ####

library(lubridate) # year, month

#### 4.3 CALIBRATION OF TEMPERATURE PROJECTIONS ####

# Load selected parameters
load("outdata/file/study_parameters.RData")

# Load temperature observations
load("outdata/file/data_tempmort.RData")

# Load temperature projections
load("outdata/file/data_projection_temperature_ssp2rcp45.RData")

# Load ISIMIP v3 method for bias correction
source("isimip3_from_masselot.R")

# Keep only the temperature columns on the observations
data_temp <- data_tempmort[, c("date", "tmean")]
rm(data_tempmort)

# Function to apply bias correction to GCM temperature projections
Run_BC_Temperature_Projections <- function(
    data_obshist,     # Observed historical temperature data
    data_simhist,     # Historical modeled temperature from GCMs
    data_simfut,      # Future modeled temperature from GCMs
    period_obshist = NULL,  # Period to use from observed historical data
    period_simhist = NULL,  # Period to use from historical simulations
    period_simfut = NULL,   # Period to use from future projections
    var_obshist,      # Column name or 'use_gcm' for temp. gcm column
    var_simhist,      # Column name or 'use_gcm' for temp. gcm column
    var_simfut,       # Column name or 'use_gcm' for temp. gcm column
    calperiod         # List of calibration periods
    ){
  
  # Loop through each selected GCM model
  proj_temp_bc <- lapply(study_param$selected_gcm, function(i_gcm) {
    
    # Dynamically generate variable names if 'use_gcm' is specified
    if (var_obshist == "use_gcm") {var_obshist <- paste0("temp.", i_gcm)}
    if (var_simhist == "use_gcm") {var_simhist <- paste0("temp.", i_gcm)}
    if (var_simfut == "use_gcm") {var_simfut <- paste0("temp.", i_gcm)}
    
    proj_temp_bc <- lapply(calperiod, function(i_calperiod) {
      
      print(paste0(
        i_gcm, ": ", i_calperiod[1], "-", i_calperiod[length(i_calperiod)]))
      
      # Default to current calibration period if none provided
      if (is.null(period_obshist)) {period_obshist <- i_calperiod}
      if (is.null(period_simhist)) {period_simhist <- i_calperiod}
      if (is.null(period_simfut)) {period_simfut <- i_calperiod}
      
      # Loop over each calendar month
      proj_temp_bc <- lapply(1:12, function(i_month){
        
        # Subset datasets for given month and year period
        subset_obshist <- subset(
          data_obshist, 
          (year(date) %in% period_obshist) & (month(date) == i_month))
        subset_simhist <- subset(
          data_simhist, 
          (year(date) %in% period_simhist) & (month(date) == i_month))
        subset_simfut <- subset(
          data_simfut, 
          (year(date) %in% period_simfut) & (month(date) == i_month))
        
        # Extract temperature values
        obshist <- subset_obshist[[var_obshist]]
        simhist <- subset_simhist[[var_simhist]]
        simfut <- subset_simfut[[var_simfut]]
        
        # Extract corresponding years
        yearobshist <- year(subset_obshist[["date"]])
        yearsimhist <- year(subset_simhist[["date"]])
        yearsimfut <- year(subset_simfut[["date"]])
        
        # Apply ISIMIP3 bias correction method
        temp_bc <- isimip3(
          obshist = obshist,
          simhist = simhist,
          simfut = simfut,
          yearobshist = yearobshist,
          yearsimhist = yearsimhist,
          yearsimfut = yearsimfut)
        
        proj_temp_bc <- data.frame(
          date = subset_simfut[["date"]],
          temp = temp_bc[,1])
        
        proj_temp_bc
        
      })
      
      # Combine all months for this calibration period
      proj_temp_bc <- do.call(rbind, proj_temp_bc)
      proj_temp_bc
      
    })
    
    # Combine all calibration periods for this GCM
    proj_temp_bc <- do.call(rbind, proj_temp_bc)
    colnames(proj_temp_bc)[colnames(proj_temp_bc) == "temp"] <- 
      paste0("temp.", i_gcm)
    proj_temp_bc
    
  })
  proj_temp_bc <- 
    Reduce(function(x, y) merge(x, y, by = "date", all = TRUE), proj_temp_bc)
  
  proj_temp_bc
  
}

# Apply bias correction to full future projections
proj_temp_bc <- Run_BC_Temperature_Projections(
  data_obshist = data_temp,
  data_simhist = proj_temp,
  data_simfut = proj_temp,
  period_obshist = 1990:2011,
  period_simhist = 1990:2011,
  period_simfut = NULL,
  var_obshist = "tmean",
  var_simhist = "use_gcm",
  var_simfut = "use_gcm",
  calperiod = list(
    1950:2011,
    2012:2029,
    2030:2039,
    2040:2049, 
    2050:2059,
    2060:2069, 
    2070:2079,
    2080:2089, 
    2090:2099))

# Apply bias correction to demo future projections (TODO: IMPROVE THIS COMMENT)
proj_temp_demo <- Run_BC_Temperature_Projections(
  data_obshist = proj_temp_bc,
  data_simhist = proj_temp_bc,
  data_simfut = proj_temp_bc,
  period_obshist = 2007:2011,
  period_simhist = NULL,
  period_simfut = NULL,
  var_obshist = "use_gcm",
  var_simhist = "use_gcm",
  var_simfut = "use_gcm",
  calperiod = lapply(0:17, function(i){2010:2014+5*i}))

# Save the final bias-corrected temperature projections
save(proj_temp_bc, file = "outdata/file/data_projection_temperature_biascorrection_ssp2rcp45.RData")
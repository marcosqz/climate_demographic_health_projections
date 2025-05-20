#### TODO: ADD DESCRIPTION OF THIS SCRIPT

# Bias-correction of temperature projections (SSP2-4.5) for London

#### LOAD LIBRARIES ############################################################

library(lubridate) # year, month

#### LOAD DATA #################################################################

load("outdata/file/study_parameters.RData")
load("outdata/file/data_tempmort.RData")
load("outdata/file/data_projection_temperature_ssp2rcp45.RData")
# Load ISIMIP v3 method for bias correction
source("isimip3_from_masselot.R")

#### BIAS-CORRECTION TEMPERATURE PROJECTIONS ###################################

# KEEP ONLY TEMPEARTURE OBSERVATIONS FROM THE MORT/TEMP DATASET
data_temp <- data_tempmort[, c("date", "tmean")]
rm(data_tempmort)

# BUILD FUNCTION TO APPLY BIAS CORRECTION TO TEMPERATURE PROJECTIONS
Run_BC_Temperature_Projections <- function(
    data_obshist,     # Observed historical temperature data
    data_simhist,     # Historical modeled temperature from GCMs
    data_simfut,      # Future modeled temperature from GCMs
    period_obshist = NULL,  # Period to use from observed historical data
    period_simhist = NULL,  # Period to use from historical simulations
    period_simfut = NULL,   # Period to use from future projections
    var_obshist,      # Var name or if 'use_gcm' takes paste0("temp.", i_gcm)
    var_simhist,      # Var name or if 'use_gcm' takes paste0("temp.", i_gcm)
    var_simfut,       # Var name or if 'use_gcm' takes paste0("temp.", i_gcm)
    calperiod         # Loop over some calibration periods
    ){
  
  # LOOP THROUGH EACH SELECTED GCM
  proj_temp_bc <- lapply(study_param$selected_gcm, function(i_gcm) {
    
    print(paste0("Run bias-correction: GCM ",  
                 i_gcm, " (", which(i_gcm == study_param$selected_gcm), "/",
                 length(study_param$selected_gcm), ")"))
    
    # GENERATE VARIABLES NAMES IF "use_gcm" is specified
    if (var_obshist == "use_gcm") {var_obshist <- paste0("temp.", i_gcm)}
    if (var_simhist == "use_gcm") {var_simhist <- paste0("temp.", i_gcm)}
    if (var_simfut == "use_gcm") {var_simfut <- paste0("temp.", i_gcm)}
    
    # LOOP THROUGH EACH CALIBRATION PERIOD
    proj_temp_bc <- lapply(calperiod, function(i_calperiod) {
      
      # DEFAUL TO CURRENT CALIBRATION PERIOD IF NONE PROVIDED
      if (is.null(period_obshist)) {period_obshist <- i_calperiod}
      if (is.null(period_simhist)) {period_simhist <- i_calperiod}
      if (is.null(period_simfut)) {period_simfut <- i_calperiod}
      
      # LOOP OVER EACH CALENDAR MONTH
      proj_temp_bc <- lapply(1:12, function(i_month){
        
        # SUBSET DATASETS FOR GIVEN MONTH AND YEAR PERIOD
        subset_obshist <- subset(
          data_obshist, 
          (year(date) %in% period_obshist) & (month(date) == i_month))
        subset_simhist <- subset(
          data_simhist, 
          (year(date) %in% period_simhist) & (month(date) == i_month))
        subset_simfut <- subset(
          data_simfut, 
          (year(date) %in% period_simfut) & (month(date) == i_month))
        
        # EXTRACT TEMPERATURE VALUES
        obshist <- subset_obshist[[var_obshist]]
        simhist <- subset_simhist[[var_simhist]]
        simfut <- subset_simfut[[var_simfut]]
        
        # EXTRACT YEARS FOR ALL THREE DATASETS
        yearobshist <- year(subset_obshist[["date"]])
        yearsimhist <- year(subset_simhist[["date"]])
        yearsimfut <- year(subset_simfut[["date"]])
        
        # APPLY ISIMIP3 BIAS CORRECTION METHOD
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
        
        return(proj_temp_bc)
        
      })
      
      # COMBINE ALL MONTHS FOR THIS CALIBRATION PERIOD
      proj_temp_bc <- do.call(rbind, proj_temp_bc)
      
      return(proj_temp_bc)
      
    })
    
    # COMBINE ALL CALIBRATION PERIODS FOR THIS GCM
    proj_temp_bc <- do.call(rbind, proj_temp_bc)
    colnames(proj_temp_bc)[colnames(proj_temp_bc) == "temp"] <- 
      paste0("temp.", i_gcm)
    
    return(proj_temp_bc)
    
  })
  # MERGE THE DATASETS FOR ALL GCMs IN A WIDE FORMAT
  proj_temp_bc <- 
    Reduce(function(x, y) merge(x, y, by = "date", all = TRUE), proj_temp_bc)
  
  return(proj_temp_bc)
  
}

# APPLY BIAS CORRECTION TO FUTURE PROJECTIONS
proj_temp_bc <- Run_BC_Temperature_Projections(
  data_obshist = data_temp,
  data_simhist = proj_temp,
  data_simfut = proj_temp,
  period_obshist = 1990:2011,
  period_simhist = 1990:2011,
  period_simfut = NULL, # NULL to take calibration period from calperiod
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

# APPLY BIAS CORRECTION TO HELD TEMPERATURE CONSTANT OVER THE PROJECTION PERIOD
proj_temp_demo <- Run_BC_Temperature_Projections(
  data_obshist = proj_temp_bc,
  data_simhist = proj_temp_bc,
  data_simfut = proj_temp_bc,
  period_obshist = 2007:2011,
  period_simhist = NULL, # NULL to take calibration period from calperiod
  period_simfut = NULL, # NULL to take calibration period from calperiod
  var_obshist = "use_gcm",
  var_simhist = "use_gcm",
  var_simfut = "use_gcm",
  calperiod = lapply(0:17, function(i){2010:2014+5*i}))

#### SAVE OUTPUTS ##############################################################

save(
  proj_temp_bc, 
  file = "outdata/file/data_projection_temperature_biascorrection_ssp2rcp45.RData")
save(
  proj_temp_demo, 
  file = "outdata/file/data_projection_temperature_biascorrection_constant.RData")
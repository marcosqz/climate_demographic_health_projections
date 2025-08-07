################################################################################

# TODO: Investigate and ask about the source of the isimip3 function (change name of the file)
# TODO: Masselot uses decade calibration periods. Is it correct to use directly
#       1950:2099

# This script aligns projected and observed temperatures, calibrating GCM 
# outputs using the ISIMIP3BASD method, with observed temperatures at the 
# 1990–2011 historical period as reference. It also creates an alternative
# projection series with static temperature distributions by remapping future
# decades to match the 2007-2011 observed (TODO???) distribution using 

#### LOAD LIBRARIES ############################################################

library(lubridate) # year, month

#### LOAD DATA #################################################################

load("indata/processed/study_parameters.RData")
load("indata/processed/data_obs_temp_mort.RData")
load("indata/processed/data_proj_temp_ssp245.RData")
source("isimip3_from_masselot.R")

#### BIAS-CORRECTION TEMPERATURE PROJECTIONS ###################################

# KEEP ONLY TEMPEARTURE OBSERVATIONS FROM THE MORT/TEMP DATASET
data_temp <- subset(data_tempmort, select = c("date", "tmean"))
rm(data_tempmort)

# DEFINE PERIODS
period_obshist <- 1990:2011 # Observational historical period (reference)
period_simhist <- 1990:2011 # Simulations historical period (used to calculate the bias to respect to the observations)
period_simfut <- 1950:2099 # Simulations future period (where the bc is applied)

# LOOP THROUGH EACH SELECTED GCM
proj_temp_bc <- lapply(study_param$selected_gcm, function(i_gcm) {
  
  print(paste0("Run bias-correction: GCM ",  
               i_gcm, " (", which(i_gcm == study_param$selected_gcm), "/",
               length(study_param$selected_gcm), ")"))
  
  # LOOP OVER EACH CALENDAR MONTH
  proj_temp_bc <- lapply(1:12, function(i_month){
    
    # SUBSET DATASETS FOR GIVEN MONTH AND YEAR PERIOD
    subset_obshist <- subset(
      data_temp, 
      (year(date) %in% period_obshist) & (month(date) == i_month))
    subset_simhist <- subset(
      proj_temp, 
      (year(date) %in% period_simhist) & (month(date) == i_month))
    subset_simfut <- subset(
      proj_temp, 
      (year(date) %in% period_simfut) & (month(date) == i_month))
    
    # EXTRACT TEMPERATURE VALUES
    obshist <- subset_obshist$tmean
    simhist <- subset_simhist[[paste0("temp.", i_gcm)]]
    simfut <- subset_simfut[[paste0("temp.", i_gcm)]]
    
    # EXTRACT YEARS FOR ALL THREE DATASETS
    yearobshist <- year(subset_obshist$date)
    yearsimhist <- year(subset_simhist$date)
    yearsimfut <- year(subset_simfut$date)
    
    # APPLY ISIMIP3 BIAS CORRECTION METHOD
    temp_bc <- isimip3(
      obshist = obshist,
      simhist = simhist,
      simfut = simfut,
      yearobshist = yearobshist,
      yearsimhist = yearsimhist,
      yearsimfut = yearsimfut)
    
    proj_temp_bc <- data.frame(
      date = subset_simfut$date,
      temp = temp_bc[,1])
    
    return(proj_temp_bc)
    
  })
  
  # COMBINE ALL MONTHS
  proj_temp_bc <- do.call(rbind, proj_temp_bc)
  colnames(proj_temp_bc)[colnames(proj_temp_bc) == "temp"] <- 
    paste0("temp.", i_gcm)
  
  return(proj_temp_bc)
  
})
# MERGE THE DATASETS FOR ALL GCMs IN A WIDE FORMAT
proj_temp_bc <- 
  Reduce(function(x, y) merge(x, y, by = "date", all = TRUE), proj_temp_bc)

#### SAVE OUTPUTS ##############################################################
save(proj_temp_bc, file = paste0(
  "outdata/file/03_calibrated_climate_projections/data_proj_temp_biascorrection_",
  study_param$ssp_rcp_scenario,".RData"))


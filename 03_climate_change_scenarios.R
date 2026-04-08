################################################################################

# This script aligns projected and observed daily temperatures by 
# bias-correcting General Circulation Model (GCM) outputs using the ISIMIP3BASD 
# method, with observed temperatures from the 1990–2011 period as the reference.
# It then computes the daily attributable fraction (AF) of heat-related deaths 
# for each age group, GCM, and epidemiological curve by combining the 
# exposure-response functions with the bias-corrected temperature projections.

rm(list = ls())

#### LOAD LIBRARIES ############################################################

library(lubridate) # year, month
library(dlnm) # onebasis

#### LOAD DATA #################################################################

load("indata/processed/study_parameters.RData")
load("indata/processed/data_obs_temp_mort.RData")
load("indata/processed/data_proj_temp_ssp245.RData")

# Load ISIMIP3 bias-correction function
source("isimip3.R")

# Load output from the epidemiological model
load("outdata/file/01_epi_model/arglag.RData")
load("outdata/file/01_epi_model/argvar.RData")
load("outdata/file/01_epi_model/coef_age.RData")
load("outdata/file/01_epi_model/mmt_age.RData")

#### BIAS-CORRECTION TEMPERATURE PROJECTIONS ###################################
# Calibrate raw GCM temperature projections against observed temperatures using 
# the ISIMIP3 method, to correct systematic model biases relative to the 
# 1990–2011 observed reference period.

# Extract observed daily mean temperatures from the temperature-mortality
# dataset
data_temp <- subset(data_tempmort, select = c("date", "tmean"))
rm(data_tempmort)

# Define periods
period_obshist <- 1990:2011 # Observational historical period (reference)
period_simhist <- 1990:2011 # Simulations historical period (used to calculate the bias to respect to the observations)
period_simfut <- 1950:2099 # Simulations future period (where the bc is applied)

# Loop GCMs
proj_temp_bc <- lapply(study_param$selected_gcms, function(i_gcm) {
  
  print(paste0("Run bias-correction: GCM ",  
               i_gcm, " (", which(i_gcm == study_param$selected_gcms), "/",
               length(study_param$selected_gcms), ")"))
  
  # Loop over calendar months (month-wise bias correction to account for
  # seasonal differences)
  proj_temp_bc <- lapply(1:12, function(i_month){
    
    # Subset observed, simulated historical and simulated future data
    subset_obshist <- subset(
      data_temp, 
      (year(date) %in% period_obshist) & (month(date) == i_month))
    subset_simhist <- subset(
      proj_temp, 
      (year(date) %in% period_simhist) & (month(date) == i_month))
    subset_simfut <- subset(
      proj_temp, 
      (year(date) %in% period_simfut) & (month(date) == i_month))
    
    # Extract temperature values
    obshist <- subset_obshist$tmean
    simhist <- subset_simhist[[paste0("temp.", i_gcm)]]
    simfut <- subset_simfut[[paste0("temp.", i_gcm)]]
    
    # Extract the years for each dataset
    yearobshist <- year(subset_obshist$date)
    yearsimhist <- year(subset_simhist$date)
    yearsimfut <- year(subset_simfut$date)
    
    # Apply the ISIMIP method to bias-correct simulated future temperatures
    temp_bc <- isimip3(
      obshist = obshist,
      simhist = simhist,
      simfut = simfut,
      yearobshist = yearobshist,
      yearsimhist = yearsimhist,
      yearsimfut = yearsimfut)
    
    # Create output dataframe with dates and bias-corrected temperatures
    proj_temp_bc <- data.frame(
      date = subset_simfut$date,
      temp = temp_bc[,1])
    
    return(proj_temp_bc)
    
  })
  
  # Combine monthly bias-corrected datasets
  proj_temp_bc <- do.call(rbind, proj_temp_bc)
  colnames(proj_temp_bc)[colnames(proj_temp_bc) == "temp"] <- 
    paste0("temp.", i_gcm)
  
  return(proj_temp_bc)
  
})

# Merge all bias-corrected GCM datasets into a single wide dataframe
proj_temp_bc <- 
  Reduce(function(x, y) merge(x, y, by = "date", all = TRUE), proj_temp_bc)

#### ESTIMATION OF THE ATTRIBUTABLE FRACTION OF DEATHS #########################
# Calculate the daily attributable fraction (AF) of deaths due to heat exposure 
# for each age group and GCM, combining the bias-corrected temperature 
# projections with the age-specific exposure–response relationships.

# Loop age groups
af <- lapply(study_param$age_groups, function(i_age) {
  
  # Loop GCMs
  af <- lapply(study_param$selected_gcms, function(i_gcm) {
    
    # Build exposure-response basis for the age-specific mmt
    cenvec <- onebasis(
      x = mmt_age[i_age], 
      fun = argvar$fun,
      knots = argvar$knots,
      Boundary.knots = argvar$Bound)
    
    # Build exposure-response basis for projected temperatures, centered at the
    # age-specific mmt
    bcen <- scale(
      onebasis(
        x = proj_temp_bc[[paste0("temp.", i_gcm)]],
        fun = argvar$fun, 
        knots = argvar$knots,
        Boundary.knots = argvar$Bound),
      center = cenvec, 
      scale = FALSE)
    
    # Estimate relative risks using estimated and sampled coefficients
    rr <- exp(bcen %*% coef_age[[i_age]])
    
    # Calculate attributable fraction
    af <- (rr - 1) / rr # AF = (RR-1)/RR
    
    # Set AF = 0 for days with temperatures below the MMT (non-heat days)
    ind_heat <- proj_temp_bc[[paste0("temp.", i_gcm)]] > mmt_age[i_age]
    af[!ind_heat,] <- 0
    
    return(af)
    
  }); names(af) <- study_param$selected_gcms
  
  return(af)
  
}); names(af) <- study_param$age_groups
  
#### SAVE OUTPUTS ##############################################################
save(proj_temp_bc, file = paste0(
  "outdata/file/02_calibrated_climate_projections/data_proj_temp_biascorrection_",
  study_param$ssp_rcp_scenario,".RData"))
save(af, file = paste0(
  "outdata/file/02_calibrated_climate_projections/attributable_fraction_",
  study_param$ssp_rcp_scenario,".RData"))
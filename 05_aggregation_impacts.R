################################################################################

# This script performs the final steps to summarise the future heat-related 
# mortality in London over the 21st century, calculating the points estimates 
# and 95% empirical confidence intervals that combine climate and 
# epidemiological uncertainty for aggregations over time and space.

#### LOAD DATA #################################################################

load("indata/processed/study_parameters.RData")

# Load mortality projections
load(paste0(
  "outdata/file/03_calibrated_demographic_projections/",
  "data_proj_mort_spatialcal_tempcal_ssp", 
  study_param$ssp_scenario,
  ".RData"))
# Load attributable numbers
load(file = paste0(
  "outdata/file/03_calibrated_demographic_projections/",
  "attributable_number_", 
  study_param$ssp_rcp_scenario,
  ".RData"))
# Load years for global warming levels
load("indata/processed/data_global_warming_levels_ssp245.RData")

#### COMPUTE DAILY TOTAL SUMMARY HEALTH IMPACTS ################################

# Build a data frame with dates, corresponding years and decades (for time-based
# grouping)
data_time <- data.frame(
  date = proj_mort_ldn_daily$date,
  year = lubridate::year(proj_mort_ldn_daily$date),
  decade = lubridate::year(proj_mort_ldn_daily$date) %/% 10 * 10)

# ---- 1.SUMMARISE DECADAL TOTAL AND AGE-SPECIFIC HEALTH IMPACTS ----

# AGGREGATE OVER DEMOGRAPHIC DIMENSIONS (sum across age groups)
an[["total"]] <- lapply(study_param$selected_gcms, function(i_gcm) {
  
  # Sum ANs for both age groups for each GCM and epidemiological association
  list_age <- lapply(study_param$age_groups, function(i_age) {
    an[[i_age]][[i_gcm]]
  })
  list_age <- Reduce("+", list_age)
  
  return(list_age)
  
}); names(an[["total"]]) <- study_param$selected_gcms

# Loop over total and age groups
an_decade_summary_age <- lapply(names(an), function(i_age) {
  
  # AGGREGATE OVER TIME (sum over decades)
  an_decade <- lapply(study_param$selected_gcms, function(i_gcm) { # Loop GCMs
    
    # Split daily data by decade
    data_decades <- split(data.frame(an[[i_age]][[i_gcm]]), data_time$decade)
    
    # Compute decadal sums for each GCM
    data_decades <- sapply(data_decades, function(x) colSums(x))
    
    return(t(data_decades))
    
  }); names(an_decade) <- study_param$selected_gcms
  
  # DERIVE POINT ESTIMATE AND 95% CONFIDENCE INTERVAL
  
  # Compute point estimate as the mean across GCMs of ANs using only the
  # estimated coefficients of the ERF
  an_decade_pe <- rowMeans(sapply(an_decade, function(x) x[,"est"]))
  
  # Compute 95% CI from AN distribution across GCM and sampled coefficients
  # of the ERF
  an_decade_ci <- lapply(an_decade, function(x) {subset(x, select = -est)})
  an_decade_ci <- do.call(cbind, an_decade_ci)
  an_decade_ci <- t(apply(an_decade_ci, 1, quantile, c(0.025, 0.975)))
  
  # Store decadal summary results in a data frame
  an_decade_summary <- data.frame(
    time = unique(data_time$decade),
    fit = an_decade_pe,
    low = an_decade_ci[,"2.5%"],
    high = an_decade_ci[,"97.5%"])
  
  return(an_decade_summary)

}); names(an_decade_summary_age) <- names(an)

# ---- 2.TEMPORAL AGGREGATION OF THE IMPACTS BY RELEVANT TEMPORAL WINDOWS ----
# Aggregate heat-related mortality by specific temporal windows:
# (a) global warming level (GWL) of 2°C, and
# (b) fixed future period (end-of-century, 2079-2099).

# 2a. global warming level (GWL) of 2°C
an_gwl_summary_age <- lapply(names(an), function(i_age) {
  
  # AGGREGATE OVER TIME (sum years for GCM-specific GWL periods)
  an_gwl <- t(sapply(study_param$selected_gcms, function(i_gcm) {
    
    # Extract start and end years defining the GWL period for each GCM
    years_gwl <- unlist(subset(
      data_gwl, 
      (warming_level == study_param$selected_warming) & (gcm == i_gcm),
      select = c("year1", "year2")))
    
    # Sum daily ANs across the selected GWL period for all estimated and 
    # simulated erfs
    an_period <- colSums(
      an[[i_age]][[i_gcm]][data_time$year %in% c(years_gwl[1]:years_gwl[2]),])
    
    return(an_period)
    
  }))
  
  # DERIVE POINT ESTIMATE AND 95% CONFIDENCE INTERVAL
  
  # Compute point estimate as the mean across GCMs of ANs using only the
  # estimated coefficients of the ERF
  an_gwl_pe <- mean(an_gwl[,"est"])
  
  # Compute 95% CI from AN distribution across GCM and sampled coefficients
  # of the ERF
  an_gwl_ci <- quantile(an_gwl[,colnames(an_gwl)!="est"], c(0.025, 0.975))
  
  # Store GWL summary results in a data frame
  an_gwl_summary <- data.frame(
    time = "gwl",
    fit = an_gwl_pe,
    low = an_gwl_ci["2.5%"],
    high = an_gwl_ci["97.5%"])
  
  return(an_gwl_summary)
  
}); names(an_gwl_summary_age) <- names(an)

# 2b. fixed future period (end-of-century, 2079-2099).
an_endcentury_summary_age <- lapply(names(an), function(i_age) {
  
  # AGGREGATE OVER TIME (sum years for end-of-century period)
  an_endcentury <- t(
    sapply(study_param$selected_gcms, function(i_gcm) {colSums(
      an[[i_age]][[i_gcm]][data_time$year %in% c(2079:2099),])
    }))
  
  # DERIVE POINT ESTIMATE AND 95% CONFIDENCE INTERVAL
  
  # Compute point estimate as the mean across GCMs of ANs using only the
  # estimated coefficients of the ERF
  an_endcentury_pe <- mean(an_endcentury[,"est"])
  
  # Compute 95% CI from AN distribution across GCM and sampled coefficients
  # of the ERF
  an_endcentury_ci <- 
    quantile(an_endcentury[,colnames(an_endcentury)!="est"], c(0.025, 0.975))
  
  # Store end-of-century summary results in a data frame
  an_endcentury_summary <- data.frame(
    time = "endcentury",
    fit = an_endcentury_pe,
    low = an_endcentury_ci["2.5%"],
    high = an_endcentury_ci["97.5%"])
  
  return(an_endcentury_summary)
  
}); names(an_endcentury_summary_age) <- names(an)

#### SAVE OUTPUTS ##############################################################
save(an_decade_summary_age, file = paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_decades_",
  study_param$ssp_rcp_scenario,".RData"))
save(an_endcentury_summary_age, file = paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_endcentury_",
  study_param$ssp_rcp_scenario,".RData"))
save(an_gwl_summary_age, file = paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_gwl_",
  study_param$ssp_rcp_scenario,".RData"))
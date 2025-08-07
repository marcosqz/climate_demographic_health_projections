################################################################################

# This scripts estimates the health impacts of projected climate exposures by
# combining the latest climate projections with future demographic trends. 
# Specifically, it calculates heat-related mortality projections under the 
# SSP2-4.5 scenario for London. The analysis separates the contribution of
# climate change only, and quantifies the uncertainty associated with both 
# sources.

#### LOAD LIBRARIES ############################################################

library(lubridate) # year
library(dlnm) # onebasis

#### LOAD DATA #################################################################

# Load model parameters
load("indata/processed/study_parameters.RData")

# Load output from the epidemiological model
load("outdata/file/01_epi_model/arglag.RData")
load("outdata/file/01_epi_model/argvar.RData")
load("outdata/file/01_epi_model/dlnm_var.RData")
load("outdata/file/01_epi_model/coef_age.RData")
load("outdata/file/01_epi_model/coefsim_age.RData")
load("outdata/file/01_epi_model/mmt_age.RData")

# Load calibrated mortality and population projections
load(paste0(
  "outdata/file/02_calibrated_demographic_projections/",
  "data_proj_mort_popu_calibrated_ssp2.RData"))
load(paste0(
  "outdata/file/02_calibrated_demographic_projections/",
  "data_proj_mort_popu_calibrated_daily_ssp2.RData"))
load(paste0(
  "outdata/file/02_calibrated_demographic_projections/",
  "data_proj_mort_popu_calibrated_constant_ssp2.RData"))

# Load bias-corrected temperature projections
load(paste0(
  "outdata/file/03_calibrated_climate_projections/", 
  "data_proj_temp_biascorrection_ssp245.RData"))

# Load warming levels
load("indata/processed/data_global_warming_levels_ssp245.RData")

#### PREPARE TEMPERATURE PROJECTION DATASETS ###################################

# Create a function to transform temperature projections from wide to long
transform_tolong_tempproj <- function(data) {
  
  # Use reshape function to transform dataset from wide to long
  proj_temp_long <- reshape(
    data,
    varying = list(names(subset(data, select = -date))),  # all columns except 'date'
    v.names = "temperature",
    timevar = "gcm",
    times = names(subset(data, select = -date)),
    idvar = "date",
    direction = "long")
  
  # Arrange the new dataset
  proj_temp_long <- proj_temp_long[order(proj_temp_long$date), ]
  proj_temp_long$gcm <- gsub("temp\\.", "", proj_temp_long$gcm)
  rownames(proj_temp_long) <- NULL
  
  return(proj_temp_long)
  
}

# Temperature projections from wide to long
proj_temp_bc <- transform_tolong_tempproj(proj_temp_bc)

#### CALCULATE HEAT-RELATED IMPACTS ############################################

# Build a function for the calculation of the impacts
calculate_health_impacts <- function(
    epidemiological, # "point_estimate" or "simulations"
    climate, # "ensemble_mean" or "gcms"
    constant_mortpopu # TRUE or FALSE
    ) {
  
  if(epidemiological == "point_estimate") {
    # Use the estimated coefficients
    f.coef <- coef_age
  } else if (epidemiological == "simulations") {
    # Use the simulated coefficients from the Monte Carlo simulations
    f.coef <- coefsim_age
  }
  
  if(climate == "ensemble_mean") {
    # Use the ensemble mean of the selected GCMs
    f.temp <- aggregate(temperature ~ date, data = proj_temp_bc, FUN = "mean")
  } else if (climate == "gcms") {
    # Use all GCMs 
    f.temp <- proj_temp_bc
  }
  
  # Select between original or constant demographic projections
  if(!constant_mortpopu){
    # Original demographic projections
    f.mortpopu <- proj_mortpopu_daily
  } else {
    # Constant demographic projections
    f.mortpopu <- proj_mortpopu_constant
  }
  
  # Loop age-groups to compute heat-related mortality
  impacts_age <- lapply(study_param$age_groups, function(i_age){
    
    # Exposure-response basis at the age-specific mmt
    cenvec <- onebasis(
      x = mmt_age[i_age], 
      fun = argvar$fun,
      knots = argvar$knots,
      Boundary.knots = argvar$Bound)
    
    # Exposure-response basis at the projected temperatures centered at the
    # age-specific mmt
    bcen <- scale(
      onebasis(
        x = f.temp$temperature,
        fun = argvar$fun, 
        knots = argvar$knots,
        Boundary.knots = argvar$Bound),
      center = cenvec, 
      scale = FALSE)
    
    # Relative risks and attributable fractions calculation
    rr <- exp(bcen %*% f.coef[[i_age]])
    af <- (rr - 1) / rr # AF = (RR-1)/RR
    
    # Calculate heat-related deaths by setting AF=0 to temperatures below mmt
    ind_heat <- f.temp$temperature > mmt_age[i_age]
    af[!ind_heat,] <- 0
    
    # Create a data.frame to store all the impacts
    colnames(af) <- paste0("af", 1:ncol(af))
    impacts <- data.frame(af); rm(af)
    impacts$date <- f.temp$date
    if (climate == "gcms") {
      impacts$gcm <- f.temp$gcm
    }
    
    # IMPACTS FROM WIDE TO LONG (single simulation column into multiple 
    # simulation columns)
    
    if(epidemiological == "simulations") {
      
      # Loop simulations (this loop works slightly faster than reshape function)
      impacts <- lapply(1:study_param$n_sim, function(i_sim) {
        
        if(climate == "ensemble_mean"){
          
          impacts_subset <- subset(
            impacts, 
            select = c("date", paste0("af", i_sim)))
          
        } else if (climate == "gcms") {
          
          impacts_subset <- subset(
            impacts,
            select = c("date", "gcm", paste0("af", i_sim)))
          
        }
        
        # Add column with the epidemiological simulation number
        impacts_subset$simulation <- i_sim
        names(impacts_subset)[names(impacts_subset) == paste0("af", i_sim)] <-
          paste0("af.", i_age)
        
        return(impacts_subset)
        
      })
      
      impacts <- do.call(rbind, impacts) 
      
    } else if(epidemiological == "point_estimate") {
      
      names(impacts)[names(impacts) == "af1"] <- paste0("af.", i_age)
      
    }
    
    # MERGE IMPACTS DATAFRAME WITH DEMOGRAPHIC PROJECTIONS
    impacts <- merge(
      x = impacts,
      y = subset(
        f.mortpopu, 
        select = c("year", "date", 
                   paste0("mort.", i_age), 
                   paste0("popu.", i_age))),
      by = "date")
    
    # CALCULATE ATTRIBUTABLE NUMBER (AN = AF * MORT)
    impacts[[paste0("an.", i_age)]] <- # AN
      impacts[[paste0("af.", i_age)]] * # AF
      impacts[[paste0("mort.", i_age)]] # MORT

    # CALCULATE ATTRIBUTABLE NUMBER RATE (AN_RATE = AN / POPULATION * 100000)
    impacts[[paste0("an_rate.", i_age)]] <- # AN RATE
      impacts[[paste0("an.", i_age)]] / # AN
      impacts[[paste0("popu.", i_age)]] * 100000 # POPULATION * 100000
    
    # REMOVE COLUMNS WITH THE DEMOGRAPHIC DATA
    impacts[, paste0(c("mort.", "popu."), i_age)] <- NULL
    impacts <- impacts[, c("year", setdiff(names(impacts), "year"))]
    
    return(impacts)
    
  })
  # Bind by columns the age-specific impact datasets
  impacts <- do.call(cbind, impacts_age)

  # Remove the duplicated columns (date, gcm, sim) created by binding datasets
  impacts <- impacts[, !duplicated(as.list(impacts))]
  
  # CALCULATE TOTAL AN BY SUMMING AGE-SPECIFIC ANs
  impacts$an <- rowSums(
    subset(impacts, 
           select = paste0("an.", study_param$age_groups)))
  
  return(impacts)
  
}

# SET DIFFERENT COMBINATIONS FOR THE CALCULATION OF HEALTH-IMPACT PROJECTIONS
opt_epi <- c("epi_est" = "point_estimate", "epi_sim" = "simulations")
opt_clim <- c("clim_ensmean" = "ensemble_mean", "clim_gcm" = "gcms")
opt_proj <- c("full_democlim" = FALSE, "only_clim" = TRUE)

# Initialize impacts object
impacts <- list()

# Loop different combinations
for(i_epi in opt_epi) { 
  for(i_clim in opt_clim) { 
    for (i_proj in opt_proj) {
      
      print(paste0(
        "Run: epidemiological-", i_epi, 
        " / climate-", i_clim, 
        " / constant_population-", i_proj))
      
      # Define the name indexing the list of the specific combination
      name_list <- paste(names(opt_epi)[i_epi == opt_epi],
                         names(opt_clim)[i_clim == opt_clim],
                         names(opt_proj)[i_proj == opt_proj], sep = ".")
      
      # Calculate the impacts and store it in the list
      impacts[[name_list]] <- calculate_health_impacts(
        epidemiological = i_epi,
        climate = i_clim,
        constant_mortpopu = i_proj)
      
      rm(name_list)
      
}}}; rm(i_epi, i_clim, i_proj)

#### TEMPORAL AGGREGATION OF IMPACTS (DAILY TO YEARLY) #########################

# NOTE: From now we are only keeping attributable numbers (AN), but any of the
# other impacts attributable fractions (AF) or AN rates could be used in a
# similar way

an_daily <- lapply(impacts, function(data) {
  
  # List of names with "af" and "an_rate"
  col_eliminate <- c(
    paste0("af.", study_param$age_groups),
    paste0("an_rate.", study_param$age_groups))
  col_keep <- names(data)[!(names(data) %in% col_eliminate)]
  
  # Keep columns different than columns_eliminate
  data <- subset(data, select = col_keep)
  
  return(data)
  
})

# Create function to aggregate impact datasets from daily to yearly
transform_impacts_to_yearly <- function(data) {

  # Remove date column
  data_year <- subset(data, select = -date)
  
  # Response variables to be aggregated (total and age-specific AN)
  response_vars <- c("an", paste0("an.", study_param$age_groups))
  
  # Construct the formula to aggregate response by year, gcm, and simulation
  agg_formula <- paste(
    "cbind(", 
    paste(response_vars, collapse = ", "), 
    ") ~",
    paste0(colnames(data_year)[!(colnames(data_year) %in% response_vars)],
           collapse = "+ "))
  agg_formula <- as.formula(agg_formula)

  # Aggregate AN by year and the other variable (e.g. gcm, simulations)
  data_year <- aggregate(agg_formula, data = data_year, FUN = sum)
  
  # With only simulations and gcms could be done directly
  # agg_formula <- as.formula(paste0(
  #   "cbind(", 
  #   paste(response_vars, collapse = ", "), 
  #   ") ~ year + gcm + simulation"))
  # data_year <- aggregate(agg_formula, data = data, FUN = sum)
  
  return(data_year)
  
}

# Run the yearly impacts for all datasets
an_year <- lapply(an_daily, function(x) {transform_impacts_to_yearly(x)})

# CALCULATE THE ONLY DEMOGRAPHIC SCENARIOS

calculate_only_demo <- function(data1, data2) {
  
  # Find columns that are different of any of the an columns
  # simpler version: c("year", "gcm", "simulation")
  names_diff <- names(data1)[!(names(data1) %in%
                               c("an", paste0("an.", study_param$age_groups)))]
  
  # Merge the full demographic and climate dataset with the only climate
  data_merged <- merge(
    data1,
    data2,
    by = names_diff, 
    suffixes = c("_full", "_clim"))
  
  # DEMO = FULL - CLIM
  an_demo <- sapply(
    c("", paste0(".", study_param$age_groups)), function(var) {
      data_merged[[paste0("an", var, "_full")]] - 
        data_merged[[paste0("an", var, "_clim")]]
    })
  colnames(an_demo) <- c("an", paste0("an.", study_param$age_groups))
  
  # Bind the an_demo with the year, gcm and simulation columns
  an_demo <- cbind(
    subset(data_merged, select = names_diff),
    an_demo)
  
  return(an_demo)
  
}

# Epi estimates and clim ensemble mean
an_year$epi_est.clim_ensmean.only_demo <- calculate_only_demo(
  data1 = an_year$epi_est.clim_ensmean.full_democlim, 
  data2 = an_year$epi_est.clim_ensmean.only_clim)

# Epi est and clim gcms
an_year$epi_est.clim_gcm.only_demo <- calculate_only_demo(
  data1 = an_year$epi_est.clim_gcm.full_democlim,
  data2 = an_year$epi_est.clim_gcm.only_clim)

# Epi sim and clim ensmean
an_year$epi_sim.clim_ensmean.only_demo <- calculate_only_demo(
  data1 = an_year$epi_sim.clim_ensmean.full_democlim,
  data2 = an_year$epi_sim.clim_ensmean.only_clim)

# Epi sim and clim gcm
an_year$epi_sim.clim_gcm.only_demo <- calculate_only_demo(
  data1 = an_year$epi_sim.clim_gcm.full_democlim,
  data2 = an_year$epi_sim.clim_gcm.only_clim)

#### TEMPORAL AGGREGATION OF IMPACTS (BY GWLs AND FIXED PERIODS) ###############

# Initialize an_period object
an_period <- list()

# TEMPORAL AGGREGATION TO GWL (2C) PERIID

# Function to compute the AN in the GWL 21-year period
# (NOTE: It could be simplified if we only use simulation and gcms)
aggregation_gwl <- function(data, method) {
  
  # Merge dataset with yearly ANs with the 21-year period of GWL
  data <- merge(
    data,
    subset(data_gwl, 
           warming_level == study_param$selected_warming,
           select = c("gcm", "year1", "year2")))
  
  # Select only years within the GWL
  data <- subset(data, (year >= year1) & (year <= year2))
  
  # Aggregate ANs by the gcms (or simulation for uncertainty)
  if(method == "point_estimate") {
    # Point estimate for GWL differ by GCMs
    formula <- cbind(an, an.00_74, an.75plus) ~ gcm
  } else if (method == "uncertainty") {
    # For the uncertainty we also consider the simulations of the epi models
    formula <- cbind(an, an.00_74, an.75plus) ~ gcm + simulation
  }
  
  # Aggregate data
  data <- aggregate(formula, data = data, FUN = "sum")
  
  if(method == "point_estimate") {
    # For the point estimate give the mean of the point estimate for the
    # selected gcms
    data <- colMeans(subset(data, select = -gcm))
    
  } else if(method == "uncertainty") {
    # For the uncertainty extract quantile 2.5 and 97.5% of the ensemble
    # of the combination of simulations of the epi models and the gcms
    data <- apply(
      subset(data, select = c("an", paste0("an.", study_param$age_groups))),
      2, quantile, c(0.5, 0.025, 0.975))
  }

  return(data)
  
}

# Full demographic and climate
an_period$gwl$full_democlim <- 
  aggregation_gwl(data = an_year$epi_sim.clim_gcm.full_democlim,
                  method = "uncertainty")

# Only climate
an_period$gwl$only_clim <- 
  aggregation_gwl(data = an_year$epi_sim.clim_gcm.only_clim,
                  method = "uncertainty")

# Only demographic
an_period$gwl$only_demo <- 
  aggregation_gwl(data = an_year$epi_sim.clim_gcm.only_demo,
                  method = "uncertainty")

# TEMPORAL AGGREGATION END-OF-CENTURY 2079-2099

# Function to compute the AN in fixed 21-year period
# (NOTE: This function is different from the previous one, because in the
# previous function we are forced to use different GCMs in the point estimate, 
# because each GCM has different periods. For the fixed we can use directly
# the point estimate. If we use only point estimate and all gcms for everything,
# this part will be simplified)
aggregation_fixed_period <- function(data, method, year_period = 2079:2099) {
  
  # Subset data to a 21-year period for the end of century
  data <- subset(
    data,
    year %in% year_period, 
    select = -year)
  
  if (method == "point_estimate") {
    
    data <- colSums(data)
    
  } else if (method == "uncertainty") {
    
    data <- aggregate(
      cbind(an, an.00_74, an.75plus) ~ gcm + simulation, 
      data = data, 
      FUN = "sum")
    
    data <- apply(
      subset(data, select = c("an", paste0("an.", study_param$age_groups))), 
      2, quantile, c(0.5, 0.025, 0.975))
  }
  
  return(data)
  
}

# Full demographic and climate
an_period$end_century$full_democlim <- 
  aggregation_fixed_period(data = an_year$epi_sim.clim_gcm.full_democlim,
                           method = "uncertainty")

# Only climate
an_period$end_century$only_clim <- 
  aggregation_fixed_period(data = an_year$epi_sim.clim_gcm.only_clim,
                           method = "uncertainty")

# Only demographic
an_period$end_century$only_demo <- 
  aggregation_fixed_period(data = an_year$epi_sim.clim_gcm.only_demo,
                           method = "uncertainty")

#### SAVE OUTPUTS ##############################################################

# Save heat-related mortality by year
save(an_year, file = paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_year_",
  study_param$ssp_rcp_scenario,".RData"))

# Save heat-related mortality by 21-year periods
save(an_period, file = paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_period_",
  study_param$ssp_rcp_scenario,".RData"))
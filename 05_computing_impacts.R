#### TODO: ADD DESCRIPTION OF THIS SCRIPT

# Calculate heat-related mortality projections in London under SSP2-4.5: 
# age-specific temperature-mortality associations (<75 and ≥75 years), 
# age-specific calibrated demographic projections under SSP2, and bias-corrected
# temperature projections under SSP2-4.5

#### LOAD LIBRARIES ############################################################

library(lubridate) # year
library(dlnm) # onebasis

#### LOAD DATA #################################################################

# Load model parameters
load("outdata/file/study_parameters.RData")

# Load output from the epi model
load("outdata/file/epi_model_argvar.RData")
load("outdata/file/epi_model_arglag.RData")
load("outdata/file/epi_model_dlnm_varibles.RData")
load("outdata/file/epi_model_coefsimage.RData")
load("outdata/file/epi_model_mmt.RData")

# Load calibrated mortality and population projections
load("outdata/file/data_projection_mortpopu_calibrated_ssp2rcp45.RData")

# Load mortality and population observations
load("outdata/file/data_tempmort.RData")
load("outdata/file/data_popu.RData")

# Load bias-corrected temperature projections
load("outdata/file/data_projection_temperature_biascorrection_ssp2rcp45.RData")
load("outdata/file/data_projection_temperature_biascorrection_constant.RData")

# Load warming levels
load("outdata/file/data_warming_level_window_ssp2rcp45.RData")

#### EXPAND DEMOGRAPHIC PROJECTIONS TO ANNUAL VALUES ###########################

# Repeat each row 5 times to create one row per year (within each 5-year block)
proj_mortpopu_expanded <- 
  proj_mortpopu_cal[rep(1:(nrow(proj_mortpopu_cal)), each = 5),]

# Replace repeated year values with consecutive years
# (1950, 1950, 1950, 1950, 1950, 1955 => 1950, 1951, 1952, 1953, 1954, 1955)
proj_mortpopu_expanded$year <- 
  c(sapply(proj_mortpopu_cal$year, function(y) {y + 0:4}))

# Divide the yearly mortality values to get daily mortality values
proj_mortpopu_expanded[, paste0("mort.", c("00_74", "75plus"))] <-
  proj_mortpopu_expanded[, paste0("mort.", c("00_74", "75plus"))] / 365

#### CREATE A DUMMY CONSTANT MORTALITY AND PROJECTIONS #########################

# Create dataset with yearly mortality observations
data_mort <- data_tempmort[, c("year", "mort.00_74", "mort.75plus")]
rm(data_tempmort)
data_mort <- aggregate(. ~ year, data = data_mort, FUN = sum)
# Remove row for 2012 (we don't have mortality data for the whole year)
data_mort <- subset(data_mort, year != 2012)

# Calculate the mean of last 5 years
data_mort <- colMeans(tail(data_mort, 5))
data_popu <- colMeans(tail(data_popu, 5))

# Maintain the data.frame structure
data_mort <- data.frame(t(data_mort))
data_popu <- data.frame(t(data_popu))

# Repeat the mean of last 5 years for the whole projection period
years_proj <-  unique(year(proj_temp_bc$date))
data_mort <- data_mort[rep(1, length(years_proj)),]
data_popu <- data_popu[rep(1, length(years_proj)),]
data_mort$year <- years_proj
data_popu$year <- years_proj

# Transform yearly mortalities to daily
data_mort$mort.00_74 <- data_mort$mort.00_74 / 365
data_mort$mort.75plus <- data_mort$mort.75plus / 365

# Bind constant mortality and population
proj_mortpopu_constant <- 
  cbind(data_mort, data_popu[, c("popu.00_74", "popu.75plus")])

#### PREPARE TEMPERATURE PROJECTION DATASETS ###################################

# Create a function to transform temperature projections from wide to long
transform_tolong_tempproj <- function(x) {
  proj_temp_long <- reshape(
    x,
    varying = list(names(x)[-1]),  # all columns except 'date'
    v.names = "temperature",
    timevar = "gcm",
    times = names(x)[-1],
    idvar = "date",
    direction = "long")
  
  proj_temp_long <- proj_temp_long[order(proj_temp_long$date), ]
  proj_temp_long$gcm <- gsub("temp\\.", "", proj_temp_long$gcm)
  rownames(proj_temp_long) <- NULL
  proj_temp_long
}

# Temperature projections from wide to long
proj_temp_bc <- transform_tolong_tempproj(proj_temp_bc)
proj_temp_demo <- transform_tolong_tempproj(proj_temp_demo)

#### CALCULATE HEAT-RELATED IMPACTS ############################################

# Build a function for the calculation of the impacts
calculate_health_impacts <- function(constant_temp,
                                     constant_mortpopu) {
  
  # Select between projection temperatures or constant temperatures
  if(!constant_temp){
    f.temp <- proj_temp_bc
  } else {
    f.temp <- proj_temp_demo
  }
  
  # Select between demographic projections or constant mortality/population
  if(!constant_mortpopu){
    f.mortpopu <- proj_mortpopu_expanded
  } else {
    f.mortpopu <- proj_mortpopu_constant
  }
  
  # Loop age-groups to compute heat-related mortality
  impacts_age <- lapply(study_param$age_groups, function(i_age){
    
    print(i_age)
    
    # Exposure-response basis at the age-specific mmt
    cenvec <- onebasis(mmt_age[i_age], fun = argvar$fun, knots = argvar$knots,
                       Boundary.knots = argvar$Bound)
    
    # Centered exposure-response basis at each daily projected temperature
    # TODO: Consider that this will not work with interaction b-dlnm
    bcen <- scale(onebasis(f.temp$temperature,
                           fun = argvar$fun, knots = argvar$knots,
                           Boundary.knots = argvar$Bound),
                  center = cenvec, scale = FALSE)
    
    # Relative risks and attributable fractions calculation
    rrsim <- exp(bcen %*% coefsim_age[[i_age]])
    afsim <- (rrsim - 1) / rrsim # AF = (RR-1)/RR
    
    # Calculate heat-related deaths by setting AF=0 to temperatures below mmt
    ind_heat <- f.temp$temperature > mmt_age[i_age]
    afsim[!ind_heat,] <- 0
    
    
    # STRUCUTRE IMPACTS AS A DATA.FRAME
    # Each column as a simulation of the epidemiological model
    colnames(afsim) <- paste0("sim", 1:100)
    impacts <- data.frame(afsim); rm(afsim)
    # Add columns with date and gcm (datasets need to be ordered)
    impacts$date <- f.temp$date
    impacts$gcm <- f.temp$gcm
    # Transform afsim from wide to long (slighter faster than reshape)
    impacts_long <- lapply(1:100, function(i) {
      
      impacts_subset <- impacts[,c("date", "gcm", paste0("sim", i))]
      # Add column with epidemiological simulation
      impacts_subset$simulation <- i
      colnames(impacts_subset)[colnames(impacts_subset) == paste0("sim", i)] <-
        paste0("af.", i_age)
      impacts_subset
      
    })
    impacts_long <- do.call(rbind, impacts_long)
    # Merge the impacts data.frame with mortality and population projections
    impacts_long$year <- year(impacts_long$date)
    impacts_long <- merge(
      impacts_long,
      f.mortpopu[, c("year", paste0("mort.", i_age), paste0("popu.", i_age))],
      by = "year")
    
    # ATTRIBUTABLE NUMBER (AN = AF * MORT)
    impacts_long[[paste0("an.", i_age)]] <- # AN
      impacts_long[[paste0("af.", i_age)]] * # AF
      impacts_long[[paste0("mort.", i_age)]] # MORT 
    # TODO: Mortality is kept constant trough all the period, check if I want to consider the past seasonality of mortality)
    
    # ATTRIBUTABLE NUMBER RATE (AN_RATE = AN / POPULATION * 100000)
    impacts_long[[paste0("an_rate.", i_age)]] <- # AN RATE
      impacts_long[[paste0("an.", i_age)]] / # AN
      impacts_long[[paste0("popu.", i_age)]] * 100000 # POPULATION * 100000
    
    return(impacts_long[, c("date", "gcm", "simulation", 
                            paste0(c("af.", "an.", "an_rate."), i_age))])
    
  })
    
  # Bind by columns the age-specific impact datasets
  impacts <- do.call(cbind, impacts_age)
  
  # Remove the duplicated rows (date, gcm, sim) created by binding datasets
  impacts <- impacts[, !duplicated(as.list(impacts))]
  impacts$year <- year(impacts$date)
  
  # Total AN by summing age-specific ANs for day, sim, gcm
  impacts$an <- rowSums(impacts[, paste0("an.", study_param$age_groups)])
  
  return(impacts)
  
}

# CALCULATE THE IMPACTS FOR PROJECTIONS (FULL, CLIMATE AND DEMOGRAPHIC 
# CONTRIBUTIONS)
# Combined climate and demographic contributions
impacts_full <- 
  calculate_health_impacts(constant_temp = FALSE, constant_mortpopu = FALSE)
# Demographic contribution (temperatures are held constant)
impacts_exclude_temp <- 
  calculate_health_impacts(constant_temp = TRUE, constant_mortpopu = FALSE)
# Climate contribution (demographics are held constant)
impacts_exclude_demo <- 
  calculate_health_impacts(constant_temp = FALSE, constant_mortpopu = TRUE)

#### TEMPORAL AGGREGATION OF IMPACTS (DAILY TO YEARLY) #########################

# Create function to aggregate impact datasets from daily to yearly
transform_impacts_to_yearly <- function(impacts) {
  
  # Subset only AN (heat-related mortality)
  impacts_year <- 
    impacts[ , c("year", "gcm", "simulation", "an", 
                 paste0("an.", study_param$age_groups))]
  
  # Response variables to be aggregated
  response_vars <- c("an", paste0("an.", study_param$age_groups))
  # Construct the formula to aggregate response by year (and gcm+simulation)
  agg_formula <- paste("cbind(", paste(response_vars, collapse = ", "), 
                       ") ~ year + gcm + simulation")
  agg_formula <- as.formula(agg_formula)
  
  # Aggregate 
  impacts_year <- aggregate(agg_formula, data = impacts_year, FUN = sum)
  
  return(impacts_year)
  
}

# TEMPORAL AGGREGATION TO YEARLY IMPACTS
impacts_year_full <- 
  transform_impacts_to_yearly(impacts = impacts_full)
impacts_year_exclude_temp <- 
  transform_impacts_to_yearly(impacts = impacts_exclude_temp)
impacts_year_exclude_demo <- 
  transform_impacts_to_yearly(impacts = impacts_exclude_demo)

#### TEMPORAL AGGREGATION OF IMPACTS (BY PERIOD - GWLs) ########################

# Create function to aggregate yearly impact datasets to longer periods
compute_an_warming <- function(impacts_year){

  # Loop for each gcm (GLW periods depend on the GCM)
  an_warming <- lapply(study_param$selected_gcms, function(i_gcm) {
    
    # Extract the start and end of the 21-year period for that specific gcm and
    # level of warming
    years_warming <-
      data_warming[data_warming$model_run == i_gcm, c("year1", "year2")]
    
    # Subset for the gcm and the warming period
    subset_data <- subset(
      impacts_year, (gcm == i_gcm) & 
        (year %in% years_warming$year1:years_warming$year2),
      select = c(simulation, an))
    
    # Sum AN for the 21-year period
    subset_data <- aggregate(an ~ simulation, subset_data, FUN = sum)
    subset_data$gcm <- i_gcm
    subset_data[, c("gcm", "simulation", "an")]
    
  })
  an_warming <- do.call(rbind, an_warming)
  
  return(an_warming)
}

# TEMPORAL AGGREGATION TO GWL (2C)
an_warming_full <- compute_an_warming(impacts_year_full)
an_warming_exclude_temp <- compute_an_warming(impacts_year_exclude_temp)
an_warming_exclude_demo <- compute_an_warming(impacts_year_exclude_demo)

#### SAVE OUTPUTS ##############################################################

save(an_warming_full, 
     file = "outdata/file/attributable_number_warming_2C_full.RData")
save(an_warming_exclude_temp, 
     file = "outdata/file/attributable_number_warming_2C_exclude_temp.RData")
save(an_warming_exclude_demo, 
     file = "outdata/file/attributable_number_warming_2C_exclude_demo.RData")
save(impacts_year_full, 
     file = "outdata/file/attributable_number_warming_years_full.RData")
save(impacts_year_exclude_temp, 
     file = "outdata/file/attributable_number_warming_years_exclude_temp.RData")
save(impacts_year_exclude_demo, 
     file = "outdata/file/attributable_number_warming_years_exclude_demo.RData")
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
    constant_mortpopu # TRUE or FALSE
    ) {
  
  # Use the estimated coefficients
  f.coef <- coef_age
  # Use the simulated coefficients from the Monte Carlo simulations
  f.coefsim <- coefsim_age
  
  # Use the bias-corrected temperature projections
  f.temp <- proj_temp_bc
  
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
    
    # Relative risks with the coefficients of the epidemiological association 
    rrest <- exp(bcen %*% f.coef[[i_age]])
    
    # Relative risks with the simulations of the epidemiological association 
    rrsim <- exp(bcen %*% f.coefsim[[i_age]])
    
    # Bind estimated and simulations (change names to differentiate them)
    colnames(rrest) <- "est"
    colnames(rrsim) <- paste0(
      "sim", formatC(1:study_param$n_sim, 
                     width = (log10(study_param$n_sim)+1), 
                     flag = "0"))
    rr <- cbind(rrest, rrsim)
    
    # Calculate attributable fraction
    af <- (rr - 1) / rr # AF = (RR-1)/RR
    
    # Calculate heat-related deaths by setting AF=0 to temperatures below mmt
    ind_heat <- f.temp$temperature > mmt_age[i_age]
    af[!ind_heat,] <- 0
    
    # Create a data.frame to store all the impacts
    impacts <- data.frame(
      date = f.temp$date,
      gcm = f.temp$gcm,
      af); rm(af)
    
    # IMPACTS FROM WIDE TO LONG (single simulation column into multiple 
    # simulation columns)
      
    # Loop simulations (this loop works slightly faster than reshape function)
    impacts <- lapply(
      colnames(impacts)[!(colnames(impacts) %in% c("date", "gcm"))], 
      function(i_col) {
      
      # Subset the data for a specific column with AF values                    
      impacts_subset <- subset(
        impacts,
        select = c("date", "gcm", i_col))
      
      # Add column with the epidemiological estimate/simulation number
      impacts_subset$epi <- i_col
      names(impacts_subset)[names(impacts_subset) == i_col] <-
        paste0("af.", i_age)
      
      return(impacts_subset)
    
    })
    impacts <- do.call(rbind, impacts)
    # Reorder columns
    impacts <- impacts[, c("date", "gcm", "epi", paste0("af.", i_age))]
        
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
  
  # Reorder rows
  impacts <- impacts[
    order(impacts$date, impacts$gcm, impacts$epi),]
  
  return(impacts)
  
}


# SET DIFFERENT COMBINATIONS FOR THE CALCULATION OF HEALTH-IMPACT PROJECTIONS
opt_proj <- c("full_democlim" = FALSE, "only_clim" = TRUE)

# Initialize impacts object
impacts <- list()

# Loop different combinations
for (i_proj in opt_proj) {
      
  print(paste0(
    "Run daily impacts: constant_population-", i_proj))
  
  # Define the name indexing the list of the specific combination
  name_list <- names(opt_proj)[i_proj == opt_proj]
  
  # Calculate the impacts and store it in the list
  impacts[[name_list]] <- calculate_health_impacts(
    constant_mortpopu = i_proj)
  
  rm(name_list)
      
}; rm(i_proj)

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
    ") ~ year + gcm + epi")
  agg_formula <- as.formula(agg_formula)

  # Aggregate AN by year and the other variable (e.g. gcm, simulations)
  data_year <- aggregate(agg_formula, data = data_year, FUN = sum)
  
  return(data_year)
  
}

# Run the yearly impacts for all datasets
an_year <- lapply(an_daily, function(x) {transform_impacts_to_yearly(x)})

# CALCULATE THE ONLY DEMOGRAPHIC SCENARIOS

# Check that full and climate scnarion have the same order
if(all.equal(an_year$full_democlim[, c("year", "gcm", "epi")],
             an_year$only_clim[, c("year", "gcm", "epi")])) {
  
  # DEMO = FULL-CLIM
  an_year$only_demo <- cbind(
    an_year$full_democlim[, c("year", "gcm", "epi")],
    an_year$full_democlim[, c("an", paste0("an.", study_param$age_groups))] -
      an_year$only_clim[, c("an", paste0("an.", study_param$age_groups))])
  
}

#### TEMPORAL AGGREGATION OF IMPACTS (BY GWLs AND FIXED PERIODS) ###############

# Initialize an_period object
an_period <- list()

# TEMPORAL AGGREGATION TO GWL (2C) PERIID

# Function to compute the AN in the GWL 21-year period
# (NOTE: It could be simplified if we only use simulation and gcms)
aggregation_period <- function(
    data, 
    method, # "gwl" for using the 21-year global warming periods by GCM or
            # "fixed" for a unique fixed period for all GCMs
    period = NULL # specify the first and last year of the fixed period for
                  # method = "fixed"
    ) {
  
  if(method == "gwl") {
    
    # Merge dataset with yearly ANs with the 21-year period of GWL
    data <- merge(
      data,
      subset(data_gwl, 
             warming_level == study_param$selected_warming,
             select = c("gcm", "year1", "year2")))
    
    # Select only years within the GWL
    data <- subset(data, (year >= year1) & (year <= year2))
    
  } else if (method == "fixed") {
    # Select only years within the 21-year fixed period
    data <- subset(data, (year >= period[1]) & (year <= period[2]))
  } 
    
  # Response variables to be aggregated (total and age-specific AN)
  response_vars <- c("an", paste0("an.", study_param$age_groups))
  
  # Construct the formula to aggregate response by gcm and epi uncertainty
  agg_formula <- paste(
    "cbind(", paste(response_vars, collapse = ", "), ") ~ gcm + epi")
  agg_formula <- as.formula(agg_formula)
  
  # Aggregate data
  data <- aggregate(agg_formula, data = data, FUN = "sum")
  
  # Calculate point estimate
  estimate <- colMeans(subset(data, epi == "est", select = response_vars))
  
  # Calculate confidence intervals
  confidence_interval <- apply(
    subset(data, epi != "est", select = response_vars), 2, 
    quantile, c(0.025, 0.975))
  
  output <- rbind(
    estimate,
    confidence_interval)
  rownames(output) <- c("est", "ci1", "ci2")
  
  return(output)
  
}

an_period$gwl <- lapply(an_year, function(x) {
  aggregation_period(data = x, method = "gwl")})
an_period$end_century <- lapply(an_year, function(x) {
  aggregation_period(data = x, method = "fixed", period = c(2079, 2099))})

#### SAVE OUTPUTS ##############################################################

# Save heat-related mortality by year
save(an_year, file = paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_year_",
  study_param$ssp_rcp_scenario,".RData"))

# Save heat-related mortality by 21-year periods
save(an_period, file = paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_period_",
  study_param$ssp_rcp_scenario,".RData"))
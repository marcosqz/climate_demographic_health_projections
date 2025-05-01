#### LOAD LIBRARIES ####

library(lubridate) # year
library(dlnm) # onebasis

#### LOAD DATA ####

# Load model parameters
load("outdata/file/study_parameters.RData")

# Load output from the epi model
load("outdata/file/epi_model_argvar.RData")
load("outdata/file/epi_model_arglag.RData")
load("outdata/file/epi_model_dlnm_varibles.RData")
load("outdata/file/epi_model_coefsimage.RData")
load("outdata/file/epi_model_mmt.RData")

# Load calibrated mortality and population projections
load("outdata/file/data_projection_temperature_calibrated_ssp2rcp45.RData")

# Load bias-corrected temperature projections
load("outdata/file/data_projection_temperature_biascorrection_ssp2rcp45.RData")

# Load warming levels
load("outdata/file/data_warming_level_window_ssp2rcp45.RData")

#### 6. COMPUTING CLIMATE CHANGE IMPACTS ACCOUNTING FOR ALL SCENARIOS ####

# EXPAND THE MORTALITY AND POPULATION PROJECTIONS TO ANNUAL VALUES INSTEAD OF 5-YEARLY VALUES
# Repeat each row 5 times to create one row per year (within each 5-year block)
proj_mortpopu_expanded <- 
  proj_mortpopu_cal[rep(1:(nrow(proj_mortpopu_cal)), each = 5),]

# Replace repeated year values with consecutive years
# (e.g., 1950, 1950, 1950, 1950, 1950, 1955 => 1950, 1951, 1952, 1953, 1954, 1955)
proj_mortpopu_expanded$year <- 
  c(sapply(proj_mortpopu_cal$year, function(y) {y + 0:4}))

# Divide the yearly mortality values to get daily mortality values
proj_mortpopu_expanded[, paste0("mort.", c("00_74", "75plus"))] <-
  proj_mortpopu_expanded[, paste0("mort.", c("00_74", "75plus"))] / 365

# TEMPERATURE PROJECTIONS DATASETS FROM WIDE TO LONG
proj_temp_bc <- reshape(
  proj_temp_bc,
  varying = list(names(proj_temp_bc)[-1]),  # all columns except 'date'
  v.names = "temperature",
  timevar = "gcm",
  times = names(proj_temp_bc)[-1],
  idvar = "date",
  direction = "long"
)
proj_temp_bc <- proj_temp_bc[order(proj_temp_bc$date), ]
proj_temp_bc$gcm <- gsub("temp\\.", "", proj_temp_bc$gcm)
rownames(proj_temp_bc) <- NULL

# LOOP FOR AGE GROUPS TO COMPUTE THE AGE-SPECIFIC ATTRIBUTABLE MEASURES
impacts_age <- lapply(study_param$age_groups, function(i_age){

  print(i_age)

  # EXPOSURE RESPONSE BASIS CENTERED AT THE AGE-SPECIFIC MMT
  cenvec <- onebasis(mmt_age[i_age], fun = argvar$fun, knots = argvar$knots,
                     Boundary.knots = argvar$Bound)

  # PREDICT THE EXPOSURE-RESPONSE ASSOCIATION CENTERED AT THE MMT FOR THE
  # SELECTED TEMPERATURE PROJECTIONS
  bcen <- scale(onebasis(proj_temp_bc$temperature,
                         fun = argvar$fun, knots = argvar$knots,
                         Boundary.knots = argvar$Bound),
                center = cenvec, scale = FALSE)

  # CALCULATE THE SAMPLES OF RELATIVE RISKS AND ATTRIBUTABLE FRACTIONS
  rrsim <- exp(bcen %*% coefsim_age[[i_age]])
  afsim <- (rrsim - 1) / rrsim # AF = (RR-1)/RR
  
  # STRUCTURE THE ATTRIBUTABLE FRACTION DATA AS A DATA FRAME WITH DATES AND GCMS
  colnames(afsim) <- paste0("sim", 1:100)
  impacts <- data.frame(afsim); rm(afsim)
  # AFSIM AND PROJ_TEMP_BC ARE SORTED THE SAME WAY TO MATCH ROWS CORRECTLY
  impacts$date <- proj_temp_bc$date
  impacts$gcm <- proj_temp_bc$gcm

  # TRANSFORM AFSIM FROM WIDE TO LONG (SLIGHTER FASTER THAN RESHAPE)
  impacts_long <- lapply(1:100, function(i) {

    impacts_subset <- impacts[,c("date", "gcm", paste0("sim", i))]
    impacts_subset$simulation <- i
    colnames(impacts_subset)[colnames(impacts_subset) == paste0("sim", i)] <-
      paste0("af.", i_age)
    impacts_subset

  })
  impacts_long <- do.call(rbind, impacts_long)
  
  # TODO: REMOVE
  # afsim_long <- reshape(
  #   afsim,
  #   varying = list(paste0("sim", 1:100)),
  #   v.names = "value",
  #   timevar = "simulation",
  #   times = paste0("sim", 1:100),
  #   direction = "long")

  # MERGE THE ATTRIBUTABLE FRACTION DATASET WITH THE MORTALITY AND POPULATION
  # PROJECTIONS
  impacts_long$year <- year(impacts_long$date)
  impacts_long <- merge(
    impacts_long,
    proj_mortpopu_expanded[, c("year",
                               paste0("mort.", i_age),
                               paste0("popu.", i_age))],
    by = "year")

  # CALCULATE THE ATTRIBUTABLE NUMBER (AN = AF * MORT)
  impacts_long[[paste0("an.", i_age)]] <- # AN
    impacts_long[[paste0("af.", i_age)]] * # AF
    impacts_long[[paste0("mort.", i_age)]] # MORT (TODO: Mortality is kept constant trough all the period, check if I want to consider the past seasonality of mortality)

  # CALCULATE THE ATTRIBUTABLE NUMBER RATE (AN_RATE = AN / POPULATION * 100000)
  impacts_long[[paste0("an_rate.", i_age)]] <- # AN RATE
    impacts_long[[paste0("an.", i_age)]] / # AN
    impacts_long[[paste0("popu.", i_age)]] * 100000 # POPULATION * 100000

  impacts_long[, c("date", "gcm", "simulation", paste0(c("af.", "an.", "an_rate."), i_age))]

})

# BIND BY COLUMNS THE AGE-SPECIFIC IMPACT DATASETS  
impacts <- do.call(cbind, impacts_age)

# REMOVE THE DUPLICATED ROWS (DATE, GCM, SIMULATION) CREATED BY BINDIN THE DATASETS
impacts <- impacts[, !duplicated(as.list(impacts))]
impacts$year <- year(impacts$date)

# CREATE A COLUMN WITH THE TOTAL ATTRIBUTABLE NUMBER BY SUMMING THE AN FOR THE DIFFRENT AGE GROUPS
impacts$an <- rowSums(impacts[, paste0("an.", study_param$age_groups)])

# WE WILL KEEP ANNUAL IMPACTS
# SELECT ONLY COLUMNS WITH THE AN IMPACT
impacts_year <- impacts[, 
  c("year", "gcm", "simulation", "an", paste0("an.", study_param$age_groups))]

# AGGREGATE ANs BY YEAR, GCM AND SIMULATION
response_vars <- c("an", paste0("an.", study_param$age_groups))
# Construct the formula string
formula_str <- paste("cbind(", paste(response_vars, collapse = ", "), 
                     ") ~ year + gcm + simulation")

# Convert to a formula
agg_formula <- as.formula(formula_str)

# Apply aggregate
impacts_year <- aggregate(agg_formula, data = impacts_year, FUN = sum)

# CALCULTE THE IMPACT FOR THE WARMING PERIOD

# Loop for each gcm
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

# SAVE FINAL RESULTS
save(an_warming, file = "outdata/file/attributable_number_warming_2C.RData")
save(impacts_year, file = "outdata/file/attributable_number_warming_years.RData")
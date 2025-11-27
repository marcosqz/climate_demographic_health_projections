################################################################################
# This script prepares and harmonises input datasets for an illustrative study 
# projecting future heat-related mortality in London over the 21st century. 
#
# The data sources include:
#   (1) Observed temperature and mortality
#   (2) Observed population data
#   (3) National population and survival ratio projections
#   (4) Climate model temperature projections (GCMs)
#   (5) Global Warming Level (GWL) periods

#### LOAD LIBRARIES ############################################################

library(lubridate) # dow
library(eurostat) # get_eurostat
library(wcde) # get_wcde0
library(sf) # st_read / st_transform
library(raster) # brick
library(exactextractr) # exact_extract

#### SELECT STUDY PARAMETERS ###################################################

# Define parameters for the illustrative health impact projection study
study_param <- list(
  age_groups = c("00_74", "75plus"), # <75 and ≥75 years
  n_sim = 100, # Number of simulations for the epidemiological models
  ssp_scenario = 2, # Shared Socioeconomic Pathway 2
  ssp_rcp_scenario = "ssp245", # SSP2-4.5 / choose between c(ssp126, ssp245, ssp370, ssp585)
  selected_gcms = c("ACCESS-CM2", "BCC-CSM2-MR", "CESM2"), # GCM model
  selected_warming = "2") # Global warming level of 2°C / choose between (1.5, 2, 3, 4)

#### DATASET 1: OBSERVED TEMPERATURE AND MORTALITY ##############

# Load daily temperature and mortality data for London (1990–2012)
# (Raw source file also available locally at indata/raw/lndn_obs.csv)
url_data <- "https://raw.githubusercontent.com/gasparrini/2019_vicedo-cabrera_Epidem_Rcodedata/refs/heads/master/lndn_obs.csv" 
data_tempmort <- read.csv(url_data)
data_tempmort$date <- as.Date(data_tempmort$date, format = "%d/%m/%Y")

# Re-group age categories (<75 and ≥75 years)
data_tempmort$mort.00_74 <- data_tempmort$all_0_64 + data_tempmort$all_65_74
data_tempmort$mort.75plus <- data_tempmort$all_75_84 + data_tempmort$all_85plus

# Keep relevant variables
data_tempmort <- subset(
  data_tempmort, 
  select = c("date", "year", "dow", "tmean", "mort.00_74", "mort.75plus"))

rm(url_data)

#### DATASET 2: OBSERVED POPULATION ############################################

# Download and aggregate annual population estimates for London (region code 
# "UKI") from Eurostat
data_popu <- get_eurostat(
  id = "demo_r_d2jan", 
  filters = list(geo = "UKI", sex = "T"))
data_popu <- subset(data_popu, select = c("age", "time", "values"))
data_popu$year <- year(data_popu$time) # e.g "2005-01-01" to "2005"
data_popu$time <- NULL

# Aggregate into two age groups (<75 and ≥75)
data_popu_00_74 <- subset(data_popu, age %in% c("Y_LT1", paste0("Y", 1:74)))
data_popu_75plus <- subset(data_popu, age %in% c(paste0("Y", 75:99), "Y_OPEN"))

# Aggregate counts by year and merge both age groups
data_popu_00_74 <- setNames(
  aggregate(values ~ year, data = data_popu_00_74, FUN = sum),
  c("year", "popu.00_74"))
data_popu_75plus <- setNames(
  aggregate(values ~ year, data = data_popu_75plus, FUN = sum),
  c("year", "popu.75plus"))
data_popu <- merge(data_popu_00_74, data_popu_75plus)
rm(data_popu_00_74, data_popu_75plus)

#### DATASET 3: POPULATION AND MORTALITY PROJECTIONS FOR UK AND NI #############

# Retrieve national-level demographic projections (population and survival 
# ratios) from the Wittgenstein Centre Human Capital Data Explorer (WCDE v2)

# Population size (thousands), 1 value every 5 year: 1950, 1955, 1960, ...
proj_popu <- get_wcde(
  indicator = "pop", # Population Size (000's)
  scenario = study_param$ssp_scenario,
  pop_age =  "all",
  pop_sex = "both",
  version = "wcde-v3",
  country_name = "United Kingdom of Great Britain and Northern Ireland")
colnames(proj_popu)[colnames(proj_popu) == "pop"] <- "popu"

# Age-Specific Survival Ratios (ASSR), 1 value per 5 year-period intervals: 
# 1950-1955, 1955-1960, ...
proj_mort <- get_wcde(
  indicator = "assr", # Age-Specific Survival Ratio
  scenario = study_param$ssp_scenario,
  pop_age =  "all",
  version = "wcde-v3",
  country_name = "United Kingdom of Great Britain and Northern Ireland")
proj_mort <- subset(proj_mort, age != "Newborn") # Exclude "Newborn" age group

# Convert 5-year intervals (e.g., “2000–2005”) to start year (e.g., 2000)
proj_mort$year <- as.numeric(substr(proj_mort$period, 1, 4))

# Merge projections, scale to absolute counts, and compute annual deaths
proj_mort <- merge(proj_popu, proj_mort)
rm(proj_popu)
proj_mort$popu <- proj_mort$popu * 1000
proj_mort$mort <- proj_mort$popu * (1 - proj_mort$assr)

# We only keep the mortality projections, but keep population projections too 
# if AN rates are one of the health outcome of interest
proj_mort <- subset(proj_mort, select = c("age", "sex", "year", "mort"))

# Aggregate by age and year, combining both sexes
proj_mort <- aggregate(mort ~ age + year, data = proj_mort, FUN = sum)

# Define age groups and aggregate into <75 and ≥75
proj_mort$age_group <- ifelse(
  test = proj_mort$age %in% 
    c("75--79", "80--84", "85--89", "90--94", "95--99", "100+"),
  yes = "75plus",
  no = "00_74")
proj_mort$age <- NULL
proj_mort <- aggregate(mort ~ age_group + year, data = proj_mort, FUN = sum)

# Reshape long to wide and rescale annual deaths
proj_mort <- reshape(
  data = proj_mort, 
  timevar = "age_group", 
  idvar = "year", 
  direction = "wide")
proj_mort[,paste0("mort.", study_param$age_groups)] <-
  proj_mort[,paste0("mort.", study_param$age_groups)]/5

# Expand to one row per year within each 5-year block
years_demodata <- proj_mort$year
proj_mort <- proj_mort[rep(1:(nrow(proj_mort)), each = 5),]
proj_mort$year <- c(sapply(years_demodata, function(y) {y + 0:4}))
rm(years_demodata)

#### DATASET 4: TEMPERATURE PROJECTIONS ############

# Download temperature projections for SSP2-RCP4.5
# This section downloads gridded daily temperature projections from the 
# NEX-GDDP-CMIP6 dataset hosted by NASA. The download fetches data for selected
# GCMs and years (1950-2100), covering both historical and future periods.

# IMPORTANT NOTES:
# - The download may take a long time depending on your internet connection,
#   NASA's server availability, and file size (~33.5 MB in this case). It may 
#   take up to 2-3 hours.
# - This step only needs to be run once. Set "run_download <- FALSE" to skip
#   it if the data have already been downloaded.
# - If you're only running this code for tutorial or learning purposes, the
#   original temperature projection data is downloaded in indata/raw/

# Direct download like this depend on evolving infrastructure. For obtaining 
# URLs for new GCMs of scenarions, check NASA's NEX-GDDP-CMIP6:
# 1) https://ds.nccs.nasa.gov/thredds/ncss/grid/AMES/NEX/GDDP-CMIP6/CESM2/ssp245/r4i1p1f1/tas/tas_day_CESM2_ssp245_r4i1p1f1_gn_2100_v2.0.nc/dataset.html
# 2) https://www.nccs.nasa.gov/services/data-collections/land-based-products/nex-gddp-cmip6

run_download <- FALSE

# Loop through GCMs and years to download daily projections (optional)
if(run_download) {
  
  # Bounding box for London
  lon <- c(-1, 1)
  lat <- c(51, 52)
  
  # Create directories
  dir.create("indata/raw/temperature_projections/", recursive = TRUE)
  
  lapply(study_param$selected_gcms, function(i_gcm) {
    
    # Create directories
    dir.create(
      paste0("indata/raw/temperature_projections/", i_gcm, "/"), 
      recursive = TRUE)
    dir.create(
      paste0("indata/raw/temperature_projections/", i_gcm, "/historical/"), 
      recursive = TRUE)
    lapply(study_param$ssp_rcp_scenario, function(i_scenario) {
      dir.create(
        paste0("indata/raw/temperature_projections/", i_gcm, "/", i_scenario), 
        recursive = TRUE)
    })
    
    lapply(1950:2100, function(i_year) {
      
      # Select scenario: "historical" BEFORE 2015, "ssp245" (or defined) afterward
      if(i_year < 2015) {
        i_scenario <- "historical"
      } else {
        i_scenario <- study_param$ssp_rcp_scenario
      }
      
      # Set correct ensemble member and time based on GCM (differ by model)
      if(i_gcm %in% c("ACCESS-CM2", "BCC-CSM2-MR")) {
        diff_url <- "r1i1p1f1"; time_url <- "12"
      } else if (i_gcm %in% c("CESM2")) {
        diff_url <- "r4i1p1f1"; time_url <- "00"
      }
      
      # Construct download URL
      url <- paste0(
        "https://ds.nccs.nasa.gov/thredds/ncss/grid/AMES/NEX/GDDP-CMIP6/", 
        i_gcm, "/", i_scenario , "/", diff_url, "/tas/tas_day_", i_gcm, "_", 
        i_scenario , "_", diff_url, "_gn_", i_year, ".nc")
      
      # Subset spatial domain (optional: omit for full global data)
      url_subset <- paste0(
        "?var=tas",
        "&north=", lat[2], 
        "&west=", lon[1], 
        "&east=", lon[2], 
        "&south=", lat[1],
        "&horizStride=1&time_start=", i_year, 
        "-01-01T", time_url, ":00:00Z&time_end=", i_year,
        "-12-31T", time_url, ":00:00Z&&&accept=netcdf3&addLatLon=true")
      
      # Define destination path
      destfile <- paste0(
        "indata/raw/temperature_projections/", i_gcm, "/", 
        i_scenario, "/proj_temp_grid_", 
        i_gcm, "_", i_scenario, "_", i_year, ".nc")
      
      # Download the NETCDF file
      print(paste0(
        "Downloading rasters: ", i_gcm, " ", i_year, 
        " (", which(i_gcm == study_param$selected_gcms), "/", 
        length(study_param$selected_gcms), ")"))
      download.file(url = paste0(url, url_subset), destfile = destfile)
      
    })
  })
  
  # (Remove temporary datasets)
  rm(lon, lat)
  
}

# Convert gridded rasters into daily mean temperature time series for London
# Load the shapefile defining the boundaries of London
# (downloaded from: https://data.london.gov.uk/dataset/statistical-gis-boundary-files-london)
shp_london <- st_read("indata/raw/shapefile_london/London_GLA_Boundary.shp")
# Transform the coordinate system of the shapefile to EPSG:4326 (WGS84) to match
# the resolution of the gridded data
shp_london <- st_transform(shp_london, 4326)

# Loop selected GCMs
proj_temp <- lapply(study_param$selected_gcms, function(i_gcm){ 
  
  print(paste0(
    "Process rasters: ", i_gcm, " (", 
    which(i_gcm == study_param$selected_gcms), 
    "/", length(study_param$selected_gcms), ")"))
  
  # Loop years
  proj_temp <- lapply(1950:2100, function(i_year) {
    
    # Load appropriate scenario (historical or future)
    if(i_year < 2015) {
      file_path <- paste0(
        "indata/raw/temperature_projections/", i_gcm, "/", "historical", 
        "/proj_temp_grid_", i_gcm, "_", "historical", "_", i_year, ".nc")
      raster_data <- raster::brick(file_path)
    } else {
      file_path <- paste0(
        "indata/raw/temperature_projections/", i_gcm, "/", 
        study_param$ssp_rcp_scenario, "/proj_temp_grid_", i_gcm, "_", 
        study_param$ssp_rcp_scenario, "_", i_year, ".nc")
      raster_data <- raster::brick(file_path)
    }
    
    # Some NetCDF files have longitudes ranging from 0 to 360 instead 
    # of -180 to 180.
    # This adjustment ensures compatibility with the London shapefile
    if (xmin(raster_data) > 180) {
      extent(raster_data) <- extent(xmin(raster_data) - 360,
                                    xmax(raster_data) - 360,
                                    ymin(raster_data),
                                    ymax(raster_data))
    }
    
    # Extract mean temperature for London
    y <- exact_extract(raster_data, shp_london, "mean")
    y <- unlist(y)
    
    # Create the dates
    dates <- seq(
      as.Date(paste0(i_year, "-01-01")),
      as.Date(paste0(i_year, "-12-31")),
      by = "day")
    
    # For leap years, we remove the values for February 29th.
    # This will not significantly affect the results and simplifies the code
    # by avoiding the need to account for different year lengths.
    if(lubridate::leap_year(i_year)){
      dates <- dates[dates != as.Date(paste0(i_year, "-02-29"))]
      if(length(y) == 366) { # Some GCMs have 365 days in leap years
        y <- y[!grepl("02.29", names(y))]
      }
    }
    
    # Create the output for 365 days
    output <- data.frame(
      date = dates,
      gcm = i_gcm,
      temp = y - 273.15 # Convert Kelvin to Celsius,
    )

    return(output)
    
  })
  
  proj_temp <- do.call(rbind, proj_temp)
  rownames(proj_temp) <- NULL
  
  return(proj_temp)
  
})
proj_temp <- do.call(rbind, proj_temp)

# (Remove temporary datasets)
rm(shp_london)

# Reshape from long to wide
proj_temp <- reshape(proj_temp, timevar = "gcm", idvar = "date", 
                     direction = "wide")
proj_temp <- proj_temp[order(proj_temp$date),]

#### DATASET 5: GLOBAL WARMING LEVELS ##########################################
# Download information to work with global warming levels (1.5°, 2°, 3°, 4°)
# under warming-levels from IPCC-WG1

# Load global warming levels (GWLs) (Raw source file also available locally at: 
# indata/raw/CMIP6_Atlas_WarmingLevels.csv)
url_data <- "https://raw.githubusercontent.com/IPCC-WG1/Atlas/main/warming-levels/CMIP6_Atlas_WarmingLevels.csv"
data_gwl <- read.csv(url_data); rm(url_data)

# Keep selected GCMs and the target SSP–RCP scenario
data_gwl <- subset(
  data_gwl,
  grepl(paste0(paste(study_param$selected_gcms, collapse = "|"), "_"), 
        model_run)) # (e.g. identify ACCESS-CM in ACCESS-CM2_r1i1p1f1)
data_gwl <- subset(
  data_gwl,
  select = grepl(study_param$ssp_rcp_scenario, names(data_gwl)) | # columns that include the name of the ssp/rcp scenario
    names(data_gwl) == "model_run") # keep also the column with the gcm models

# Clean names of GCMs (e.g. BCC-CSM2-MR_r1i1p1f1 to BCC-CSM2-MR)
data_gwl$model_run <- sub("_.*", "", data_gwl$model_run) 

# Reshape from wide to long
data_gwl <- reshape(
  data = data_gwl, 
  varying = which(names(data_gwl) != "model_run"),
  v.names = "year", 
  timevar = "warming_level", 
  times = names(data_gwl)[names(data_gwl) != "model_run"],
  direction = "long")
                    
# Clean GWLs dataset
data_gwl$warming_level <- gsub(
  paste0("X|_", study_param$ssp_rcp_scenario), "", data_gwl$warming_level)
data_gwl$id <- NULL
rownames(data_gwl) <- NULL
colnames(data_gwl)[colnames(data_gwl) == "model_run"] <- "gcm"

# Create 21-year period centered in the GCM-specific GWL years
data_gwl$year1 <- data_gwl$year - 10
data_gwl$year2 <- data_gwl$year + 10

#### SAVE OUTPUTS ##############################################################

# Save study parameters
save(study_param, file = "indata/processed/study_parameters.RData")
# Save observed temperature and mortality
save(data_tempmort, file = "indata/processed/data_obs_temp_mort.RData")
# Save observed population data
save(data_popu, file = "indata/processed/data_obs_popu.RData")
# Save temperature projections
save(proj_temp, file = paste0(
  "indata/processed/data_proj_temp_",study_param$ssp_rcp_scenario, ".RData"))
# Save mortality and population projections
save(proj_mort, file = paste0(
  "indata/processed/data_proj_mort_ssp", study_param$ssp_scenario,".RData"))
# Save Global Warming Level (GWL) periods
save(data_gwl, file = paste0(
  "indata/processed/data_global_warming_levels_",
  study_param$ssp_rcp_scenario, ".RData"))
################################################################################

# This script prepares data from multiple publicly available sources for an
# illustrative projection study of future heat-related mortality in Greater
# London over the 21st century. The data sources include:
#    1) Temperature and mortality observations
#    2) Population observations
#    3) Population and survival ratio projections
#    4) Temperature projections
#    Other) Global warming level periods

#### LOAD LIBRARIES ############################################################

library(lubridate) # dow
library(eurostat) # get_eurostat
library(wcde) # get_wcde0
library(sf) # st_read / st_transform
library(raster) # brick
library(exactextractr) # exact_extract

#### SELECT STUDY PARAMETERS ###################################################

# LIST OF PARAMETERS FOR THE ILLUSTRATIVE HEALTH IMPACT PROJECTION STUDY
study_param <- list(
  age_groups = c("00_74", "75plus"), # <75 and ≥75 years
  n_sim = 100, # Number of simulations for the epidemiological models
  ssp_scenario = 2, # Shared Socioeconomic Pathways 2
  ssp_rcp_scenario = "ssp245", # SSP2-4.5 / choose between c('ssp126, ssp245, ssp370, ssp585)
  selected_gcms = c("ACCESS-CM2", "BCC-CSM2-MR", "CESM2"), # GCM model
  selected_warming = "2") # Global warming level of 2°C / choose between (1.5, 2, 3, 4 °C)

#### DATASET 1: TEMPERATURE AND MORTALITY OBSERVATIONS FOR LONDON ##############

# LOAD TEMPERATURE AND MORTALITY OBSERVATIONS
# (Raw source file also available locally at: indata/raw/lndn_obs.csv)
data_tempmort <- read.csv(
  "https://raw.githubusercontent.com/gasparrini/2019_vicedo-cabrera_Epidem_Rcodedata/refs/heads/master/lndn_obs.csv")
data_tempmort$date<- as.Date(data_tempmort$date, format = "%d/%m/%Y")

# RE-GROUP AGE CATEGORIES (00_74 AND 75PLUS)
data_tempmort$mort.00_74 <- data_tempmort$all_0_64 + data_tempmort$all_65_74
data_tempmort$mort.75plus <- data_tempmort$all_75_84 + data_tempmort$all_85plus

# SELECT RELEVANT COLUMNS
data_tempmort <- data_tempmort[, c("date", "year", "dow", "tmean", 
                                   "mort.00_74", "mort.75plus")]

#### DATASET 2: POPULATION OBSERVATIONS FOR LONDON #############################

# DOWNLOAD POPULATION DATA FOR LONDON (REGION CODE "UKI") FROM EUROSTAT
data_popu <- get_eurostat(id = "demo_r_d2jan", 
                          filters = list(geo = "UKI", sex = "T"))

# KEEP RELEVANT COLUMNS AND EXTRACT THE YEAR FROM THE TIME VARIABLE
data_popu <- data_popu[, c("age", "time", "values")]
data_popu$year <- year(data_popu$time) # e.g "2005-01-01" to "2005"
data_popu$time <- NULL

# RE-GROUP AGE CATEGORIES (00_74 AND 75PLUS)

# (Subset age group under 75: less than 1yo "Y_LT1" and 1 "Y1" to 74 "Y74")
data_popu_00_74 <- 
  data_popu[data_popu$age %in% c("Y_LT1", paste0("Y", 1:74)), ]
# (Subset age group 75 plus: 75 "Y75" to 99 "Y99" and open-ended age class 
# "Y_OPEN")
data_popu_75plus <- 
  data_popu[data_popu$age %in% c(paste0("Y", 75:99), "Y_OPEN"), ]

# (Aggregate population counts by year within each age group)
data_popu_00_74 <- setNames(
  aggregate(values ~ year, data = data_popu_00_74, FUN = sum),
  c("year", "popu.00_74"))
data_popu_75plus <- setNames(
  aggregate(values ~ year, data = data_popu_75plus, FUN = sum),
  c("year", "popu.75plus"))

# MERGE AGGREGATED POPULATION DATA FOR THE AGE GROUPS INTO A SINGLE DATASET
data_popu <- merge(data_popu_00_74, data_popu_75plus)
# (Remove temporary datasets)
rm(data_popu_00_74, data_popu_75plus)

#### DATASET 3: POPULATION AND MORTALITY PROJECTIONS FOR UK AND NI #############

# DOWNLOAD POPULATION PROJECTIONS FROM THE WITTGENSTEIN CENTRE HUMAN CAPITAL 
# DATA EXPLORER 
# (1 value every 5 year: 1950, 1955, 1960, ...)
proj_popu <- get_wcde(
  indicator = "pop", # Population Size (000's)
  scenario = study_param$ssp_scenario,
  pop_age =  "all",
  pop_sex = "both",
  version = "wcde-v2", # we use v2, because v3 only has population from 2000
  country_name = "United Kingdom of Great Britain and Northern Ireland")
colnames(proj_popu)[colnames(proj_popu) == "pop"] <- "popu"

# DOWNLOAD SURVIVAL RATIO PROJECTIONS FROM THE WITTGENSTEIN CENTRE HUMAN CAPITAL 
# DATA EXPLORER 
# (1 value per 5 year-period intervals: 1950-1955, 1955-1960, ...)
proj_mort <- get_wcde(
  indicator = "assr", # Age-Specific Survival Ratio
  scenario = study_param$ssp_scenario,
  pop_age =  "all",
  version = "wcde-v2", # we use v2, because v3 only has population from 2000
  country_name = "United Kingdom of Great Britain and Northern Ireland")
proj_mort <- proj_mort[proj_mort$age != "Newborn",] # Exclude "Newborn" age group

# CONVERT 5-YEAR PERIOD LABELS IN MORTALITY PROJECTIONS TO A SINGLE 
# REPRESENTATIVE YEAR (e.g. convert "2000-2005" to 2000 by extracting the start
# year)
proj_mort$year <- as.numeric(substr(proj_mort$period, 1, 4))

# MERGE MORTALITY AND POPULATION PROJECTIONS
proj_mortpopu <- merge(proj_popu, proj_mort)
rm(proj_mort, proj_popu)

# SCALE POPULATIONS AND COMPUTE DEATHS
proj_mortpopu$popu <- proj_mortpopu$popu * 1000
proj_mortpopu$mort <- proj_mortpopu$popu * (1 - proj_mortpopu$assr)

# KEEP RELEVANT COLUMNS
proj_mortpopu <- proj_mortpopu[, c("age", "sex", "year", "popu", "mort")]

# AGGREGATE SEX-SPECIFIC COUNTS BY YEAR AND AGE
proj_mortpopu <- aggregate(cbind(popu, mort) ~ age + year, 
                           data = proj_mortpopu, FUN = sum)

# IDENTIFY AGE GROUPS 00_74 AND 75PLUS
proj_mortpopu$age_group <- proj_mortpopu$age %in% c(
  "75--79", "80--84", "85--89", "90--94", "95--99", "100+")
proj_mortpopu$age_group[proj_mortpopu$age_group == FALSE] <- "00_74"
proj_mortpopu$age_group[proj_mortpopu$age_group == TRUE] <- "75plus"
proj_mortpopu$age <- NULL

# AGGREGATE COUNTS BY YEAR AND AGE CATEGORIES 00_74 AND 75PLUS
proj_mortpopu <- aggregate(cbind(popu, mort) ~ age_group + year, 
                           data = proj_mortpopu, 
                           FUN = sum)

# RESPHAPE THE PROJECTION DATASET FROM LONG TO WIDE
proj_mortpopu <- reshape(proj_mortpopu, 
                         timevar = "age_group", 
                         idvar = "year", 
                         direction = "wide")

# RESCALE NUMBER OF DEATHS AS ANNUAL AVERAGE
proj_mortpopu$mort.00_74 <- proj_mortpopu$mort.00_74/5
proj_mortpopu$mort.75plus <- proj_mortpopu$mort.75plus/5

#### DATASET 4: TEMPERATURE PROJECTIONS FOR COMBINATIONS OF SSP/RCP ############

# DOWNLOAD TEMPERATURE PROJECTIONS FOR SSP2-RCP4.5 FOR THE PERIOD 2015-2100

# This section download gridded daily temperature projections from the 
# NEX-GDDP-CMIP6 dataset hosted by NASA. The download fetches data for selected
# GCMs and years (1950-2100), covering both historical and future periods.

# IMPORTANT NOTES:
# - The download may take a long time depending on your internet connection,
#   NASA's server availability, and file size (~350 MB in this case). It may 
#   take up to 2-3 hours.
# - This step only needs to be run once. Set "run_download <- FALSE" to skip
#   it if the data have already been downloaded.
# - If you're only running this code for tutorial or learning purposes, the
#   processed temperature projection data (already provided for following 
#   scripts) is sufficient, and download is not required.

# Direct download like this depend on evolving infrastructure. For obtaining 
# URLs for new GCMs of scenarions, check NASA's NEX-GDDP-CMIP6:
# 1) https://ds.nccs.nasa.gov/thredds/ncss/grid/AMES/NEX/GDDP-CMIP6/CESM2/ssp245/r4i1p1f1/tas/tas_day_CESM2_ssp245_r4i1p1f1_gn_2100_v2.0.nc/dataset.html
# 2) https://www.nccs.nasa.gov/services/data-collections/land-based-products/nex-gddp-cmip6

run_download <- TRUE

# DEFINE BOUNDING BOX AROUND LONDON (CAN BE ADJUSTED TO REDUCE FILE SIZE)
lon <- c(-5, 5)
lat <- c(50, 53)

# LOOP OVER GCMs AND YEARS TO DOWNLOAD TEMPERATURE PROJECTIONS
if(run_download) {
  
  # CREATE DIRECTORIES
  dir.create("indata/raw/temperature_projections/", recursive = TRUE)
  
  lapply(study_param$selected_gcms, function(i_gcm) {
    
    # CREATE DIRECTORIES
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
        "https://ds.nccs.nasa.gov/thredds/ncss/grid/AMES/NEX/GDDP-CMIP6/", i_gcm, 
        "/", i_scenario , "/", diff_url, "/tas/tas_day_", i_gcm, 
        "_", i_scenario , "_", diff_url, "_gn_", i_year, ".nc")
      
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
}

# (Remove temporary datasets)
rm(run_download, lon, lat)

# CONVERT GRIDDED TEMPERATURE PROJECTIONS INTO A TIME-SERIES

# Load the shapefile define the boundaries of London
# (downloaded from: https://data.london.gov.uk/dataset/statistical-gis-boundary-files-london)
shp_london <- st_read("indata/raw/shapefile_london/London_GLA_Boundary.shp")
# Transform the coordinate system of the shapefile to EPSG:4326 (WGS84) to match
# the resolution of the gridded data
shp_london <- st_transform(shp_london, 4326)

# Loop selected GCMs
proj_temp <- lapply(study_param$selected_gcms, function(i_gcm){ 
  
  print(paste0(
    "Process rasters: ", i_gcm, " (", which(i_gcm == study_param$selected_gcms), 
    "/", length(study_param$selected_gcms), ")"))
  
  # Loop years
  proj_temp <- lapply(1950:2100, function(i_year) {
    
    # Load raster data with the temperature projections
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
    
    # Create the dataset with the extracted mean temperatures
    # Handle the special case of leap years when the raster has only 365 layers:
    # this indicates that the 29th of February is missing in the data
    if(!(lubridate::leap_year(i_year) & # Check if it's a leap year
         (dim(raster_data)[3] == 365))) # And raster has only 365 days
    {
      
      # Create the dates
      dates <- seq(
        as.Date(paste0(i_year, "-01-01")),
        as.Date(paste0(i_year, "-12-31")),
        by = "day")
      
      # Create the output with an NA at the end for the 29th of Feb
      output <- data.frame(
        date = dates,
        gcm = i_gcm,
        temp = y - 273.15 # Convert Kelvin to Celsius,
      )
      
    } else {
      
      # Handle leap years when the raster is missing Feb 29
      
      # Create the dates for the leap year
      dates_leap <- seq(
        as.Date(paste0(i_year, "-01-01")),
        as.Date(paste0(i_year, "-12-31")),
        by = "day")
      
      # Identify the leap day
      leap_day <- as.Date(paste0(i_year, "-02-29"))
      
      # Move leap day to the end
      dates_leap <- c(dates_leap[dates_leap != leap_day], leap_day)
      
      # Create the output with an NA at the end for the 29th of Feb
      output <- data.frame(
        date = dates_leap,
        gcm = i_gcm,
        temp = c(y - 273.15, NA) # Convert Kelvin to Celsius,
      )
      
    }
    
    return(output)
    
  })
  
  proj_temp <- do.call(rbind, proj_temp)
  rownames(proj_temp) <- NULL
  
  return(proj_temp)
  
})
proj_temp <- do.call(rbind, proj_temp)

# (Remove temporary datasets)
rm(shp_london)

# TRANSFORM PROJECTIONS FRON LONG TO WIDE
proj_temp <- reshape(proj_temp, timevar = "gcm", idvar = "date", 
                     direction = "wide")
proj_temp <- proj_temp[order(proj_temp$date),]

# HANDLE MISSING TEMPERATURE VALUES FOR LEAP DAYS IN SOME GCM
# Some GCMs omit the 29th of February in leap years. For the missing days,
# impute the temperature as the average of the precedin (Feb 28) and following
# day (Mar 1)
for(i_gcm in study_param$selected_gcms) {
  
  # Identify missing temperature values
  ind <- is.na(proj_temp[[paste0("temp.", i_gcm)]])
  
  if(any(ind)){
    
    # Extract dates with missing values
    na_dates <- proj_temp[ind,]$date
    
    # Loop through each missing date and impute temperature
    for(i_date in na_dates){
      
      i_date <- as.Date(i_date)
      
      # Impute missing value as the mean of temperatures from the previos and
      # next day
      proj_temp[proj_temp$date == i_date, paste0("temp.", i_gcm)] <- 
        mean(
          proj_temp[proj_temp$date %in% c(i_date-1, i_date+1), 
                    paste0("temp.", i_gcm)])
    }; rm(i_date, na_dates)
    
  }; rm(ind)
  
}; rm(i_gcm)

#### DATASET 5: GLOBAL WARMING LEVELS ##########################################

# DOWNLOAD INFORMATION TO WORK WITH GLOBAL WARMING LEVELS (1.5°, 2°, 3°, 4°) 
# UNDER WARMING-LEVELS FROM IPCC-WG1

# LOAD GLOBAL WARMING LEVELS (GWLs)
# (Raw source file also available locally at: 
# indata/raw/CMIP6_Atlas_WarmingLevels.csv)
data_gwl <- read.csv(
  "https://raw.githubusercontent.com/IPCC-WG1/Atlas/main/warming-levels/CMIP6_Atlas_WarmingLevels.csv")

# KEEP SELECTED GCMs 
data_gwl <- lapply(study_param$selected_gcms, function(x) {
  data_gwl[grepl(paste0(x, "_"), data_gwl$model_run),]
})
data_gwl <- do.call(rbind, data_gwl)

# KEEP ONLY THE SELECTED SSP/RCP SCENARIO
data_gwl <- 
  data_gwl[, grepl(study_param$ssp_rcp_scenario, colnames(data_gwl)) | # columns that include the name of the ssp/rcp scenario
                 colnames(data_gwl) == "model_run"] # keep also the column with the gcm models

# CLEAN NAMES OF GCMs (e.g. BCC-CSM2-MR_r1i1p1f1 to BCC-CSM2-MR)
data_gwl$model_run <- sub("_.*", "", data_gwl$model_run) 

# GWLs DATASET FROM WIDE TO LONG
data_gwl <- reshape(data_gwl,
                    varying = which(names(data_gwl) != "model_run"),
                    v.names = "year", 
                    timevar = "warming_level", 
                    times = names(data_gwl)[names(data_gwl) != "model_run"],
                    direction = "long")

# CLEAN GWLs DATASET
data_gwl$warming_level <- 
  gsub("X|_ssp245", "", data_gwl$warming_level)
data_gwl$id <- NULL
rownames(data_gwl) <- NULL
colnames(data_gwl)[colnames(data_gwl) == "model_run"] <- "gcm"

# CREATE 21-YEARS GCM-SPECIFIC WARMING LEVEL WINDOWS
data_gwl$year1 <- data_gwl$year - 10
data_gwl$year2 <- data_gwl$year + 10

#### SAVE OUTPUTS ##############################################################

save(study_param, file = "indata/processed/study_parameters.RData")
save(data_tempmort, file = "indata/processed/data_obs_temp_mort.RData")
save(data_popu, file = "indata/processed/data_obs_popu.RData")
save(proj_mortpopu, file = paste0(
  "indata/processed/data_proj_mort_popu_ssp", 
  study_param$ssp_scenario,".RData"))
save(proj_temp, file = paste0(
  "indata/processed/data_proj_temp_",
  study_param$ssp_rcp_scenario, ".RData"))
save(data_gwl, file = paste0(
  "indata/processed/data_global_warming_levels_",
  study_param$ssp_rcp_scenario, ".RData"))
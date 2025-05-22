#### TODO: ADD DESCRIPTION OF THIS SCRIPT

# Here we downloaded both observed and projected climate and demographic 
# information from publicly available data sources to conduct a projection study
# on future heat-related mortality impacts in London during the 21st century:
#    1) Temperature and mortality observations
#    2) Population observations
#    3) Population and survival ratio projections
#    4) Temperature projections
#    5) Global warming level periods

#### LOAD LIBRARIES ############################################################

library(lubridate) # dow
library(eurostat) # get_eurostat
library(wcde) # get_wcde0
# RClimChange is an experimental package. 
# To install RClimChange use: devtools::install_github("hllauca/RClimChange")
library(RClimChange) # nex_download
library(sf) # st_read / st_transform
library(raster) # brick
library(exactextractr) # exact_extract

#### SELECT STUDY PARETERS #####################################################

# LIST OF PARAMETERS FOR THE ILLUSTRATIVE HEALTH IMPACT PROJECTION STUDY
study_param <- list(
  age_groups = c("00_74", "75plus"), # <75 and ≥75 years
  ssp_scenario = 2, # SSP2
  ssp_rcp_scenario = "ssp245", # c('ssp126, ssp245, ssp370, ssp585)
  selected_gcms = c("BCC-CSM2-MR", "MIROC6", "IPSL-CM6A-LR"), # GCM model
  selected_warming = "2") # Global warming level (1.5, 2, 3, 4 C)

#### DATASET 1: TEMPERATURE AND MORTALITY OBSERVATIONS FOR LONDON ##############

data_tempmort <- read.csv("indata/lndn_obs.csv", 
                          colClasses = c(date = "Date"))

# RE-GROUP AGE CATEGORIES (00_74 AND 75PLUS)
data_tempmort$mort.00_74 <- data_tempmort$all_0_64 + data_tempmort$all_65_74
data_tempmort$mort.75plus <- data_tempmort$all_75_84 + data_tempmort$all_85plus

# SELECT RELEVANT COLUMNS
data_tempmort <- data_tempmort[, c("date", "year", "dow", "tmean", 
                                   "mort.00_74", "mort.75plus")]

#### DATASET 2: POPULATION OBSERVATIONS FOR LONDON #############################

# DOWNLOAD POPULATION DATA FOR LONDON (UKI) FROM EUROSTAT
data_popu <- get_eurostat("demo_r_d2jan", filters = list(geo = "UKI", 
                                                         sex = "T"))

# SELECT RELEVANT COLUMNS AND EXTRACT THE YEARS
data_popu <- data_popu[, c("age", "time", "values")]
data_popu$year <- year(data_popu$time) # e.g "2005-01-01" to "2005"
data_popu$time <- NULL

# SUBSET POPULATION DATA BY AGE GROUPS
# Age group under 75 (less than 1yo and 1 to 74)
data_popu_00_74 <- 
  data_popu[data_popu$age %in% c("Y_LT1", paste0("Y", 1:74)), ]
# Age group 75 plus (75 to 99 and open-ended age class)
data_popu_75plus <- 
  data_popu[data_popu$age %in% c(paste0("Y", 75:99), "Y_OPEN"), ]

# AGGREGATE SMALLER POPULATION GROUPS FOR THE LARGER 00_74 AND 75plus
data_popu_00_74 <- setNames(
  aggregate(values ~ year, data = data_popu_00_74, FUN = sum),
  c("year", "popu.00_74"))
data_popu_75plus <- setNames(
  aggregate(values ~ year, data = data_popu_75plus, FUN = sum),
  c("year", "popu.75plus"))

# MERGE AGE-SPECIFIC DATASETS
data_popu <- merge(data_popu_00_74, data_popu_75plus)
rm(data_popu_00_74, data_popu_75plus)

#### DATASET 3: POPULATION AND MORTALITY PROJECTIONS FOR UK AND NI #############

# DOWNLOAD POPULATION PROJECTIONS (1 value every 5 year: 1950, 1955, 1960, ...)
proj_popu <- get_wcde(
  indicator = "pop", # Population Size (000's)
  scenario = study_param$ssp_scenario,
  pop_age =  "all",
  pop_sex = "both",
  version = "wcde-v2", # we use v2, because v3 only has population from 2000
  country_name = "United Kingdom of Great Britain and Northern Ireland")
colnames(proj_popu)[colnames(proj_popu) == "pop"] <- "popu"

# DOWNLOAD SURVIVAL RATIO PROJECTIONS (1 value every 5 year-period: 1950-1955, 
# 1955-1960, ...)
proj_mort <- get_wcde(
  indicator = "assr", # Age-Specific Survival Ratio
  scenario = study_param$ssp_scenario,
  pop_age =  "all",
  version = "wcde-v2", # we use v2, because v3 only has population from 2000
  country_name = "United Kingdom of Great Britain and Northern Ireland")
proj_mort <- proj_mort[proj_mort$age != "Newborn",]

# TRANSFORM THE 5-YEAR PERIODS IN THE MORTALITY PROJECTIONS INTO A UNIQUE YEAR
# (E.G. 2000-2005 TO 2000)
proj_mort$year <- as.numeric(substr(proj_mort$period, 1, 4))

# MERGE MORTALITY AND POPULATION PROJECTIONS
proj_mortpopu <- merge(proj_popu, proj_mort)
rm(proj_mort, proj_popu)

# SCALE POPULATIONS AND COMPUTE DEATHS
proj_mortpopu$popu <- proj_mortpopu$popu * 1000
proj_mortpopu$mort <- proj_mortpopu$popu * (1 - proj_mortpopu$assr)

# SELECT RELEVANT COLUMNS
proj_mortpopu <- proj_mortpopu[, c("age", "sex", "year", "popu", "mort")]

# AGGREGATE SEX GROUPS
proj_mortpopu <- aggregate(cbind(popu, mort) ~ age + year, 
                           data = proj_mortpopu, FUN = sum)

# CREATE AGE GROUPS 00_74 AND 75PLUS
proj_mortpopu$age_group <- proj_mortpopu$age %in% c(
  "75--79", "80--84", "85--89", "90--94", "95--99", "100+")
proj_mortpopu$age_group[proj_mortpopu$age_group == FALSE] <- "00_74"
proj_mortpopu$age_group[proj_mortpopu$age_group == TRUE] <- "75plus"
proj_mortpopu$age <- NULL

# AGGREGATE AGE GROUPS CONTAINING THE 00_74 AND 75PLUS
proj_mortpopu <- aggregate(cbind(popu, mort) ~ age_group + year, 
                           data = proj_mortpopu, FUN = sum)

# RESPHAPE THE PROJECTION DATASET FROM LONG TO WIDE
proj_mortpopu <- reshape(proj_mortpopu, timevar = "age_group", idvar = "year", 
                         direction = "wide")

# RESCALE NUMBER OF DEATHS AS ANNUAL AVERAGE
proj_mortpopu$mort.00_74 <- proj_mortpopu$mort.00_74/5
proj_mortpopu$mort.75plus <- proj_mortpopu$mort.75plus/5

#### DATASET 4: TEMPERATURE PROJECTIONS FOR COMBINATIONS OF SSP/RCP ############

# DOWNLOAD TEMPERATURE PROJECTIONS FOR SSP2-RCP4.5 IN THE FUTURE PERIOD  
# (2015-2100)
# Set run_download <- TRUE to download the original gridded datasets
# TODO: Specify approximate times to download
run_download <- TRUE # If already downloaded can be set to FALSE.
if(run_download == TRUE) {
  nex_download(location = getwd(),
               model = study_param$selected_gcms,
               scenario = study_param$ssp_rcp_scenario,
               variable = 'tas',
               years = 2015:2100,
               version = NULL,
               roi = c(-5, 5, 50, 53), # subset gridded dataset to a smaller area that includes London
               method = 'curl')
} 

# DOWNLOAD TEMPERATURE PROJECTIONS FOR THE HISTORICAL PERIOD (1950-2014)
if(run_download == TRUE) {
  nex_download(location = getwd(),
               model = study_param$selected_gcms,
               scenario = 'historical',
               variable = 'tas',
               years = 1950:2014,
               version = NULL,
               roi = c(-5, 5, 50, 53),
               method = 'curl')
}
rm(run_download)

# CONVERT GRIDDED TEMPERATURE PROJECTIONS INTO A TIME-SERIES
# Load the shapefile from London
# (downloaded from: https://data.london.gov.uk/dataset/statistical-gis-boundary-files-london)
shp_london <- st_read("indata/shapefile_london/London_GLA_Boundary.shp")
# Transform the coordinate system of the shapefile to EPSG:4326 (WGS84) to match
# the resolution of the gridded data
shp_london <- st_transform(shp_london, 4326)

# Define variable to find the path of the downloaded gridded dataset
var_path <- c("BCC-CSM2-MR" = "gn", "MIROC6" = "gn", "IPSL-CM6A-LR" = "gr")

# Loop selected GCMs
proj_temp <- lapply(study_param$selected_gcms, function(i_gcm){ 
  
  print(paste0(i_gcm, " (", which(i_gcm == study_param$selected_gcms), 
               "/", length(study_param$selected_gcms), ")"))
  
  # Loop years
  proj_temp <- lapply(1950:2100, function(i_year) {
    
    # Load raster data with the temperature projections
    if(i_year < 2015) {
      file_path <- paste0("../../tas/historical/", i_gcm, 
                          "/tas_day_", i_gcm, 
                          "_historical_r1i1p1f1_", var_path[i_gcm],
                          "_", i_year, ".nc")
      raster_data <- raster::brick(file_path)
    } else {
      file_path <- paste0("../../tas/ssp245/", i_gcm, 
                          "/tas_day_", i_gcm, 
                          "_ssp245_r1i1p1f1_", var_path[i_gcm],
                          "_", i_year, ".nc")
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
    
    # Extract Mean Temperature for London
    y <- exact_extract(raster_data, shp_london, "mean")
    y <- unlist(y)
    data.frame(
      date= as.Date(gsub("mean.X", "", names(y)), format = "%Y.%m.%d"),
      gcm = i_gcm,
      temp = y - 273.15 # Convert Kelvin to Celsius,
    )
    
  })
  
  proj_temp <- do.call(rbind, proj_temp)
  rownames(proj_temp) <- NULL
  
  return(proj_temp)
  
})
proj_temp <- do.call(rbind, proj_temp)
rm(shp_london)

# TRANSFORM PROJECTIONS FRON LONG TO WIDE
proj_temp <- reshape(proj_temp, timevar = "gcm", idvar = "date", 
                     direction = "wide")
proj_temp <- proj_temp[order(proj_temp$date),]

# SOME GCMs (BCC-CSM2-MR) MISS THE 31ST DECEMBER OF LEAP YEARS 
# IMPUTE THE MEAN VALUE FROM THE PRECEDING AND FOLLOWING DAYS.
na_dates <- proj_temp[is.na(proj_temp$`temp.BCC-CSM2-MR`),]$date

for(i_date in na_dates){
  
  i_date <- as.Date(i_date)
  
  proj_temp[proj_temp$date == i_date, "temp.BCC-CSM2-MR"] <- 
    mean(
      proj_temp[proj_temp$date %in% c(i_date-1, i_date+1), "temp.BCC-CSM2-MR"])
}
rm(i_date, na_dates)

#### DATASET 5: GLOBAL WARMING LEVELS ##########################################

# READ THE CSV FILE WITH GLWs DIRECTLY INTO R
data_warming <- read.csv(
  "https://raw.githubusercontent.com/IPCC-WG1/Atlas/main/warming-levels/CMIP6_Atlas_WarmingLevels.csv", 
  stringsAsFactors = FALSE)

# KEEP ONLY THE SELECTED GCMS
data_warming <- lapply(study_param$selected_gcms, function(x) {
  data_warming[grepl(x, data_warming$model_run),]
})
data_warming <- do.call(rbind, data_warming)

# KEEP ONLY THE SELECTED SSP/RCP SCENARIO
data_warming <- 
  data_warming[, grepl(study_param$ssp_rcp_scenario, colnames(data_warming)) | # columns that include the name of the ssp/rcp scenario
                 colnames(data_warming) == "model_run"] # keep also the column with the gcm models
data_warming$model_run <- sub("_.*", "", data_warming$model_run) # keep only the name of the scenario (e.g. BCC-CSM2-MR_r1i1p1f1 to BCC-CSM2-MR)

# DATASET FROM WIDE TO LONG
data_warming <- reshape(data_warming,
                        varying = which(names(data_warming) != "model_run"),
                        v.names = "year", 
                        timevar = "warming_level", 
                        times = names(data_warming)[names(data_warming) != "model_run"],
                        direction = "long")

# CLEAN THE NAMES OF THE DATASET
data_warming$warming_level <- 
  gsub("X|_ssp245", "", data_warming$warming_level)
data_warming$id <- NULL
rownames(data_warming) <- NULL

# SUBSET SELECTED WARMING LEVEL
data_warming <- subset(data_warming,
                       warming_level == study_param$selected_warming,
                       select = c("model_run", "year"))

# 21-YEARS GCM-SPECIFIC WARMING LEVEL WINDOW
data_warming$year1 <- data_warming$year - 10
data_warming$year2 <- data_warming$year + 10

#### SAVE OUTPUTS ##############################################################

save(study_param, file = "outdata/file/study_parameters.RData")
save(data_tempmort, file = "outdata/file/data_tempmort.RData")
save(data_popu, file = "outdata/file/data_popu.RData")
save(proj_mortpopu, file = "outdata/file/data_projection_mortality_population_ssp2.RData")
save(proj_temp, file = "outdata/file/data_projection_temperature_ssp2rcp45.RData")
save(data_warming, file = "outdata/file/data_warming_level_window_ssp2rcp45.RData")
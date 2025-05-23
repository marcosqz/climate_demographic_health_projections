#### TODO: ADD DESCRIPTION OF THIS SCRIPT

# Plot all figures in the tutorial

#### LOAD LIBRARIES ############################################################

library(dlnm) # onebasis
library(sf) # st_read
library(raster) # brick, xmin

#### LOAD DATA #################################################################

load("outdata/file/data_tempmort.RData")
load("outdata/file/data_popu.RData")
load("outdata/file/epi_model_argvar.RData")
load("outdata/file/epi_model_arglag.RData")
load("outdata/file/epi_model_dlnm_varibles.RData")
load("outdata/file/epi_model_coefsimage.RData")
load("outdata/file/epi_model_mmt.RData")

load("outdata/file/study_parameters.RData")

load("outdata/file/data_projection_mortality_population_ssp2.RData") # TODO: Change mortality_population to mortpopu
load("outdata/file/data_projection_mortpopu_calibrated_ssp2rcp45.RData") # TODO: Change name, remove rcp45 

load("outdata/file/data_projection_temperature_ssp2rcp45.RData")
load("outdata/file/data_projection_temperature_biascorrection_ssp2rcp45.RData")
load("outdata/file/data_projection_temperature_biascorrection_constant.RData")

load("outdata/file/attributable_number_warming_years_full.RData")
load("outdata/file/attributable_number_warming_years_exclude_temp.RData")
load("outdata/file/attributable_number_warming_years_exclude_demo.RData")

load("outdata/file/attributable_number_warming_2C_full.RData")
load("outdata/file/attributable_number_warming_2C_exclude_temp.RData")
load("outdata/file/attributable_number_warming_2C_exclude_demo.RData")

#### FIGURE 1. AGE-SPECIFIC TEMPERATURE-MORTALITY ASSOCIATIONS #################

# VISUALLITATION PARAMETERS BY AGE GROUPS
age_parameters <- data.frame(
  response = c("mort.00_74", "mort.75plus"),
  groups = c("00_74", "75plus"),
  col = c("#1B9E77", "#D95F02"))

title_plot <- c(
  "00_74" = expression("Young ("< 75*" years)"),
  "75plus" = expression("Old (">= 75*" years)"))

# Define the temperatures for the x-axis
pred_perc <- c(seq(0, 1, 0.1), 2:98, seq(99, 100, 0.1))
x_temp <- quantile(data_tempmort$tmean, pred_perc / 100)
rm(pred_perc)

# Consider also temperatures beyond the boundaries of observed temperatures
x_temp_extra <- seq(max(x_temp), max(x_temp) + 5, length = 11) 
x_temp <- c(x_temp, x_temp_extra[-1])
rm(x_temp_extra)

# PLOT PARAMETERS
ymax <- 5
perc_rr <- c(1, 5, 95, 99)/100
col_rr <- c("blue", "cyan", "orange", "red", "red4", "purple4")
# Highlight the RR in percentiles 1, 5, 95, 99, 100, 100 + 2ºC
ind_rr <- (x_temp %in% quantile(data_tempmort$tmean, perc_rr)) |
  (x_temp %in% c(max(data_tempmort$tmean), max(data_tempmort$tmean) + 2))
temp_rr <- x_temp[ind_rr]

# PLOT FIGURE WITH AGE-SPECIFIC TEMPERATURE-MORTALITY ASSOCIATIONS
nrow.fig <- 1; ncol.fig <- 2
pdf("outdata/plot/fig1_age_specific_associations.pdf",
    width = ncol.fig*3*1.5, height = nrow.fig*3*1.5)
layout(matrix(seq(ncol.fig * nrow.fig), nrow = nrow.fig, byrow = TRUE))
par(mex = 0.8, mgp = c(3, 1, 0), las = 1, oma = c(0, 0, 0, 0))

lapply(study_param$age_groups, function(i_age) {
  
  # Exposure-response basis at the age-specific mmt
  cenvec <- onebasis(mmt_age[[i_age]], 
                     fun = argvar$fun, knots = argvar$knots, 
                     Boundary.knots = argvar$Bound)
  
  # Centered exposure-response basis at each daily projected temperature
  bcen <- scale(onebasis(x_temp, fun = argvar$fun, knots = argvar$knots, 
                         Boundary.knots = argvar$Bound), 
                center = cenvec, scale = FALSE)
  
  # Relative risks samples
  rrsim <- exp(bcen %*% coefsim_age[[i_age]])
  
  # Calculate min, median and max of the sample of RR
  rrmedian <- apply(rrsim, 1, median)
  rrmin <- apply(rrsim, 1, min)
  rrmax <- apply(rrsim, 1, max)
  
  # Extract the RR that are higlighted in the plot
  focus_rr_median <- rrmedian[ind_rr]
  focus_rr_min <- rrmin[ind_rr]
  focus_rr_max <- rrmax[ind_rr]
  
  # Extract color for the correspondin age group
  col_age <- age_parameters$col[age_parameters$groups == i_age]
  
  # Plot empty plot
  plot(x_temp, rrsim[,1], 
       xlab = expression(paste("Temperature (", degree, "C)")),
       ylab = "Relative risk", ylim = c(1.000, ymax), 
       log = "y", type = "n")
  title(title_plot[i_age], line = 0.65, font.main = 1, cex.main = 1)
  
  # Plot sample of RR
  for(i in 1:100){
    lines(x_temp, rrsim[,i], col = col_age, lwd = 0.1)
  }; rm(i)
  
  # Plot median RR
  lines(x_temp, rrmedian, col = "black", lwd = 2)
  
  # Other lines in the plot
  abline(h = 1)
  abline(v = max(data_tempmort$tmean), lty = 2, col = "grey")
  
  # Plot lines and RR for the highlighted RR
  for(i in 1:length(col_rr)){
    lines(c(temp_rr[i], temp_rr[i]), c(1, focus_rr_median[i]), 
          lty = 4, col = col_rr[i], lwd = 1.5)
  }; rm(i)
  points(temp_rr, focus_rr_median, col = "black", pch = 21, 
         cex = 1.5, bg = col_rr, lwd = 2)
  
  # Plot a legend with the highlighted RR
  legend(
    "top",
    paste0("P", c(sprintf("%02d", perc_rr*100), c(100, 100)), 
           c("", "", "", "", "", "+2ºC"), ": ",
           formatC(focus_rr_median, digits = 2, format = "f", flag="0"), " (",
           formatC(focus_rr_min, digits = 2, format = "f", flag="0"), ",",
           formatC(focus_rr_max, digits = 2, format = "f", flag="0"), ")"), 
    col = col_rr, pch = 19, box.lty = 0, cex = 1.000)
  
})

# Title common to both panels
mtext("Age-specific temperature-mortality associations (London, 1990-2012)", 
      side = 3, outer = TRUE, line = -2.2, cex = 1.3, font = 2)

dev.off()

rm(age_parameters, title_plot, col_rr, ind_rr, ncol.fig, nrow.fig, perc_rr,
   temp_rr, x_temp, ymax)

#### FIGURE 2. DEMOGRAPHIC PROJECTIONS #########################################

# Transform the mortality observation dataset to yearly
data_mort_year <- 
  aggregate(cbind(mort.00_74, mort.75plus) ~ year, 
            data = data_tempmort, FUN = sum)
data_mort_year <- data_mort_year[data_mort_year$year != 2012,]

# Add total mortality and population in the projections
proj_mortpopu$mort <- proj_mortpopu$mort.00_74 + proj_mortpopu$mort.75plus
proj_mortpopu$popu <- proj_mortpopu$popu.00_74 + proj_mortpopu$popu.75plus

# VISUALLITATION PARAMETERS BY AGE GROUPS
title_plot <- c(
  "00_74" = expression("Young ("< 75*" years)"),
  "75plus" = expression("Old (">= 75*" years)"))

# Define colors for the age groups
col_plot <- c("00_74" = "#1B9E77", "75plus" = "#D95F02")
col_plot_dark <- c(
  "00_74" = adjustcolor(col_plot["00_74"], alpha.f = 1, 
                        red.f = .75, green.f = .75, blue.f = .75),
  "75plus" = adjustcolor(col_plot["75plus"], alpha.f = 1, 
                         red.f = .75, green.f = .75, blue.f = .75))
col_plot_light <- c(
  "00_74" = adjustcolor(col_plot["00_74"], alpha.f = 1, 
                        red.f = 1.25, green.f = 1.25, blue.f = 1.25),
  "75plus" = adjustcolor(col_plot["75plus"], alpha.f = 1, 
                         red.f = 1.25, green.f = 1.25, blue.f = 1.25))

# PLOT FIGURE WITH DEMOGRAPHIC PROJECTIONS
nrow.fig <- 1; ncol.fig <- 3
pdf("outdata/plot/fig2_demographic_projections.pdf",
    width = ncol.fig*2.5*1.3, height = nrow.fig*2.5*1.3)
layout(matrix(seq(ncol.fig * nrow.fig), nrow = nrow.fig, byrow = TRUE))
par(mex = 0.8, mgp = c(3, 1, 0), las = 1, oma = c(0, 0, 0, 0))

# a) PLOT PROJECTIONS OF MORTALITY
plot(x = proj_mortpopu$year, 
     y = proj_mortpopu$mort.00_74/1e5, type = "l",
     ylim = c(0, 
              max(proj_mortpopu$mort.00_74, proj_mortpopu$mort.75plus)/1e5),
     main = "a) Mortality projections \n(SSP2 - UK and Northern Ireland)",
     xlab = "Year", ylab = "Number of deahts (x100,000)", cex.lab = 1.1,
     col = col_plot[["00_74"]], lwd = 2)
lines(x = proj_mortpopu$year, 
      y = proj_mortpopu$mort.75plus/1e5, 
      lwd = 2, col = col_plot[["75plus"]])
abline(h = 0)
legend("topleft", 
       legend = c(title_plot[1], title_plot[2]), 
       fill = col_plot, 
       border = "black", bty = "o", cex = 1,
       bg = "white")

# b) PLOT PROJECTIONS OF POPULATION PERCERTANGE
plot(x = proj_mortpopu$year, 
     y = rep(0.5, length(proj_mortpopu$year)), 
     type = "n", ylim = c(0, 100),
     main = "b) Populations projections \n(SSP2 - UK and Northern Ireland)", 
     xlab = "Year", ylab = "Percentage (%)", cex.lab = 1.1)
polygon(x = c(proj_mortpopu$year, rev(proj_mortpopu$year)), 
        y = c(rep(0, length(proj_mortpopu$year)), 
              rev(proj_mortpopu$popu.00_74/proj_mortpopu$popu*100)), 
        col = col_plot[["00_74"]], border = NA)
polygon(c(proj_mortpopu$year, rev(proj_mortpopu$year)),
        c(100 - proj_mortpopu$popu.75plus/proj_mortpopu$popu*100, 
          rev(rep(100, length(proj_mortpopu$popu.75plus)))), 
        col = col_plot[["75plus"]], border = NA)
legend("bottomleft", 
       legend = c(title_plot[[1]], title_plot[[2]]), 
       fill = col_plot, 
       border = "black", bty = "o", cex = 1,
       bg = "white")

# c) PLOT CALIBRATED MORTALITY AND POPULATION PROJECTIONS

# Adjust margins: bottom, left, top, right
par(mar = c(5, 4, 4, 5) + 0.1)

# Plot population data (left y-axis)
plot(proj_mortpopu_cal$year, proj_mortpopu_cal$popu.00_74/1e5, 
     ylim = range(c(proj_mortpopu_cal$popu.00_74, 
                    proj_mortpopu_cal$popu.75plus))/1e5,
     type = "l", lwd = 2,
     col = col_plot_dark["00_74"], lty = 1,
     xlab = "Year", ylab = "Population (x100,000)", cex.lab = 1.1,
     main = "c) Calibration demographic projections \n(SSP2 - London)")
lines(proj_mortpopu_cal$year, proj_mortpopu_cal$popu.75plus/1e5, 
      lwd = 2, col = col_plot_dark["75plus"], lty = 1)
points(data_popu$year, data_popu$popu.00_74/1e5,
       pch = 4, col = col_plot_dark["00_74"], cex = 0.6)
points(data_popu$year, data_popu$popu.75plus/1e5,
       pch = 4, col = col_plot_dark["75plus"], cex = 0.6)
lines(data_popu$year, data_popu$popu.00_74/1e5,
      col = col_plot_dark["00_74"], lwd = 0.5)
lines(data_popu$year, data_popu$popu.75plus/1e5,
      col = col_plot_dark["75plus"], lwd = 0.5)
abline(h = 0)

# Add second plot with mortality data (right y-axis)
par(new = TRUE)

plot(proj_mortpopu_cal$year, proj_mortpopu_cal$mort.00_74/1e5,
     ylim = range(c(proj_mortpopu_cal$mort.00_74, proj_mortpopu_cal$mort.75plus,
                    data_mort_year$mort.00_74, data_mort_year$mort.75plus))/1e5,
     type = "l", lwd = 2, axes = FALSE, xlab = "", ylab = "",
     col = adjustcolor(col_plot_light["00_74"], alpha.f = 1), lty = 1)
lines(proj_mortpopu_cal$year, proj_mortpopu_cal$mort.75plus/1e5, 
      lwd = 2, col = adjustcolor(col_plot_light["75plus"], alpha.f = 1), lty = 1)
points(data_mort_year$year, data_mort_year$mort.00_74/1e5, pch = 4, 
       col = adjustcolor(col_plot_light["00_74"], alpha.f = 1), cex = 0.6)
points(data_mort_year$year, data_mort_year$mort.75plus/1e5, pch = 4, 
       col = adjustcolor(col_plot_light["75plus"], alpha.f = 1), cex = 0.6)
lines(data_mort_year$year, data_mort_year$mort.00_74/1e5,
      col = adjustcolor(col_plot_light["00_74"], alpha.f = 1), lwd = 0.5)
lines(data_mort_year$year, data_mort_year$mort.75plus/1e5,
      col = adjustcolor(col_plot_light["75plus"], alpha.f = 1), lwd = 0.5)

# Add right y-axis for mortality
axis(side = 4)
mtext("Number of deahts (x100,000)", side = 4, line = 3, las = 0, cex = 0.75)

# Add legend
text_legend <- c(
  expression("popu,proj,"< 75*""),
  expression("popu,proj,">= 75*""),
  expression("popu,obs,"< 75*""),
  expression("popu,obs,">= 75*""),
  expression("mort,proj,"< 75*""),
  expression("mort,proj,">= 75*""),
  expression("mort,obs,"< 75*""),
  expression("mort,obs,">= 75*""))
legend("right", 
       legend = text_legend,
       col = c(rep(col_plot_dark, 2),
               rep(col_plot_light, 2)),
       lty = c(1, 1, 1, 1, 1, 1, 1, 1),
       pch = c(NA, NA, 4, 4, NA, NA, 4, 4),
       lwd = c(2, 2, 0.5, 0.5, 2, 2, 0.5, 0.5),
       pt.cex = 1,
       cex = 0.85,
       bty = "n")

dev.off()

rm(data_mort_year, col_plot, col_plot_dark, col_plot_light, ncol.fig, nrow.fig,
   text_legend, title_plot)

#### FIGURE 3. CLIMATE PROJECTIONS #############################################

# Aggregate projection temperature to years
proj_temp_year <- 
  aggregate(.~ lubridate::year(date), data = proj_temp, FUN = mean)
proj_temp_bc_year <- 
  aggregate(.~ lubridate::year(date), data = proj_temp_bc, FUN = mean)
proj_temp_demo_year <- 
  aggregate(.~ lubridate::year(date), data = proj_temp_demo, FUN = mean)

# Arrange yearly temperature projections
proj_temp_year$date <- NULL
proj_temp_bc_year$date <- NULL
proj_temp_demo_year$date <- NULL
colnames(proj_temp_year)[1] <- "year"
colnames(proj_temp_bc_year)[1] <- "year"
colnames(proj_temp_demo_year)[1] <- "year"

# Keep only above 1990 for the plot
proj_temp_year <- subset(proj_temp_year, year >= 1990)
proj_temp_bc_year <- subset(proj_temp_bc_year, year >= 1990)
proj_temp_demo_year <- subset(proj_temp_demo_year, year >= 1990)

# Aggregate yearly temperature observations
data_temp_year <- data_tempmort[, c("date", "tmean")]
data_temp_year$year <- lubridate::year(data_temp_year$date)
data_temp_year$date <- NULL
data_temp_year <- aggregate(.~year, data = data_temp_year, FUN = mean)
data_temp_year <- data_temp_year[data_temp_year$year != 2012,]

# TODO: REMOVE this function, don't subset to specific warming at code01 to 
# load directly here the dataset
Create_data_warming <- function() {
  
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
  
  # 21-YEARS GCM-SPECIFIC WARMING LEVEL WINDOW
  data_warming$year1 <- data_warming$year - 10
  data_warming$year2 <- data_warming$year + 10
  
  return(data_warming)
  
}

data_warming <- Create_data_warming()

# PLOT PARAMETERS
col_plot_a <- c("#1F77B4", "#2CA02C", "#E73F74", "#6699CC")
col_plot_b <- c("#118DFF", "#750985", "#C83D95")
levels_gwl <- c("1.5", "2", "3")
offsets_gwl <- c(-0.15, 0, 0.15)
points_gwl <- c(15, 16, 17)

# PLOT FIGURE WITH CLIMATE MODELS
pdf(file = "outdata/plot/fig3_climate_models.pdf",
    width = 5*1.2, height = 4*1.2 + 4)
par(mfrow = c(2, 1), mar = c(3.8, 4.1, 4.1, 0.5), las = 1)

# a) TIME-SERIES TEMPERATURE PROJECTIONS FOR GCM (IPSL-CM6A-LR)
# Raw projections
plot(proj_temp_year$year, proj_temp_year$`temp.IPSL-CM6A-LR`, 
     col = col_plot_a[1], type = "l",
     ylim = range(proj_temp_year$`temp.IPSL-CM6A-LR`,
                  proj_temp_bc_year$`temp.IPSL-CM6A-LR`),
     xlab = "Year", ylab = expression(paste("Temperature (", degree, "C)")),
     main = "a) GCM IPSL-CM6A-LR (London, SSP2-4.5)",
     lwd = 2, lty = 2)
# Bias-corrected projections
lines(proj_temp_bc_year$year, proj_temp_bc_year$`temp.IPSL-CM6A-LR`, 
      col = col_plot_a[2], lwd = 2)
# Constant projections
lines(proj_temp_demo_year$year, proj_temp_demo_year$`temp.IPSL-CM6A-LR`, 
      col = col_plot_a[3], lwd = 2)
# Observations
lines(data_temp_year$year, data_temp_year$tmean,
      col = col_plot_a[4], lwd = 2)

# Add legend
legend("topleft",
       c("Observations",
         "Raw projections",
         "Bias-corrected projections",
         "Constant projections"),
       lwd = c(2, 2, 2, 2),
       lty = c(1, 2, 1, 1),
       col = col_plot_a[c(4, 1, 2, 3)],
       cex = 0.8)

# b) PERIODS OF GLOBAL WARMING LEVELS BY GCM 
par(mar = c(5, 8, 3.5, 0.5))  # (bottom, left, top, right)
plot(data_warming$year[1], 1 - 0.15,
     xlim = c(2000, 2100), ylim = c(0.5, 3.5), type = "n",
     main = "b) Global warming levels (SSP2-RCP4.5)",
     xlab = "Year", ylab = "", , yaxt = "n")
axis(2, at = 1:3, labels = unique(data_warming$model_run), las = 1) 

# Loop through each GWL
for (i in 1:length(levels_gwl)) {
  data <- subset(data_warming, warming_level == levels_gwl[i])
  
  # Plot points
  points(data$year, 1:length(data$year) + offsets_gwl[i],
         pch = points_gwl[i], col = col_plot_b[i], cex = 2)
  
  # Add segments (horizontal lines)
  segments(data$year - 10, 1:length(data$year) + offsets_gwl[i],
           data$year + 10, 1:length(data$year) + offsets_gwl[i], 
           col = col_plot_b[i], lwd = 2)
}; rm(data, i)
abline(h = seq(1.5, 18.5), col = "grey", lty = 2)

# Add legend below the title, outside the plot
legend("top", legend = paste0(levels_gwl, "ºC"), col = col_plot_b,
       pch = points_gwl, horiz = TRUE, cex = 1.2,
       bty = "n", xpd = TRUE, inset = c(0, -0.13))

dev.off()

rm(proj_temp_year, proj_temp_bc_year, proj_temp_demo_year, data_temp_year,
   Create_data_warming, offsets_gwl, levels_gwl, points_gwl,
   col_plot_a, col_plot_b)

#### FIGURE 4. HEAT-RELATED MORTALITY ##########################################

# CALCULATE SUMMARY VALUES OF HEAT-RELATED MORTALIY FOR EACH AGE GROUP

# Estimates by year
periods_an <- 1950:2099
names(periods_an) <- 1950:2099

# Get median, and CI of the AN for each group
an_summary <- list(
  
  # All age
  total = sapply(periods_an, function(i_year) {
    quantile(impacts_year_full[impacts_year_full$year %in% i_year, "an"], 
             c(0.025, 0.5, 0.975))
  }),
  
  # Young
  "00_74" = sapply(periods_an, function(i_year) {
    quantile(impacts_year_full[impacts_year_full$year %in% i_year, "an.00_74"], 
             c(0.025, 0.5, 0.975))
  }),
  
  # Old
  "75plus" = sapply(periods_an, function(i_year) {
    quantile(impacts_year_full[impacts_year_full$year %in% i_year, "an.75plus"], 
             c(0.025, 0.5, 0.975))
}))

# TEMPORAL AGGREGATION OF IMPACTS (BY PERIOD - END OF THE CENTURY)
# TODO: Move this calculation to code 05_computing_impacts

# Create function to aggregate yearly impact datasets to longer periods
compute_an_warming <- function(impacts_year){
  
  # Loop for each gcm
  an_warming <- lapply(study_param$selected_gcms, function(i_gcm) {
    
    # First and last year of the 21-period
    years_warming <- c(2079, 2099)
    
    # Subset for the gcm and the period
    subset_data <- subset(
      impacts_year, (gcm == i_gcm) & 
        (year %in% years_warming[1]:years_warming[2]),
      select = c(simulation, an))
    
    # Sum AN for the 21-year period
    subset_data <- aggregate(an ~ simulation, subset_data, FUN = sum)
    subset_data$gcm <- i_gcm
    subset_data[, c("gcm", "simulation", "an")]
    
  })
  an_warming <- do.call(rbind, an_warming)
  an_warming
}

# Temporal aggregation to end-of-century (2079-2099)
an_endcentury_full <- compute_an_warming(impacts_year_full)
an_endcentury_exclude_temp <- compute_an_warming(impacts_year_exclude_temp)
an_endcentury_exclude_demo <- compute_an_warming(impacts_year_exclude_demo)

# CALCULATE SUMMARY VALUES OF HEAT-RELATED MORTALIY FOR 21-YEAR PERIODS

# Put in a list all the datasets of 21-years period
datasets_period <- list(an_warming_exclude_demo,
                        an_warming_exclude_temp,
                        an_warming_full,
                        an_endcentury_exclude_demo,
                        an_endcentury_exclude_temp,
                        an_endcentury_full)

# Extract medians and quantiles of 21-years periods
values_period <- 
  sapply(datasets_period, function(df) median(df$an)) / 21
min_ci_period <- 
  sapply(datasets_period, function(df) quantile(df$an, 0.025)) / 21
max_ci_period <- 
  sapply(datasets_period, function(df) quantile(df$an, 0.975)) / 21

# Names of the barplots corresponding to each 21-year period dataset
names(values_period) <- rep(c("Climate", 
                              "Demographic", 
                              "Climate + \ndemographic"), 2)
  
# PLOT PARAMETERS
col_plot_a <- c("total" = "#000000", "00_74" = "#1B9E77", "75plus" = "#D95F02")
col_plot_b <- c("#C969A1", "#CE4441")
title_plot_a <- list(
  "Total (Old age)",
  expression("Young ("< 75*" years)"),
  expression("Old (">= 75*" years)"))

# PLOT HEALT-RELATED IMPACTS
nrow.fig <- 1; ncol.fig <- 2
pdf("outdata/plot/fig4_healt_impacts.pdf",
    width = ncol.fig*3*1.5, height = nrow.fig*3*1.5)
layout(matrix(seq(ncol.fig * nrow.fig), nrow = nrow.fig, byrow = TRUE))
par(mex = 1, mgp = c(3, 1, 0), las = 1, oma = c(0, 0, 0, 0))

# a) YEARLY TIME-SERIES HEAT-RELATED MORTALITY BY AGE GROUPS
plot(names(periods_an),
     an_summary$total[2,],
     ylim = c(0, max(an_summary$total[2,], na.rm = TRUE)), type = "n",
     xlab = "Year",
     ylab = "Yearly heat-related deaths",
     main = "", )
title("a) Age-specific time evolution", 
      line = 0.65, font.main = 1, cex.main = 1)

# Plot 95% CI
polygon(x = c(names(periods_an), rev(names(periods_an))),
        y = c(an_summary$total[1,], rev(an_summary$total[3,])),
        border = NA, col = adjustcolor(col_plot_a["total"], alpha.f = 0.2))
polygon(x = c(names(periods_an), rev(names(periods_an))),
        y = c(an_summary$`00_74`[1,], rev(an_summary$`00_74`[3,])),
        border = NA, col = adjustcolor(col_plot_a["00_74"], alpha.f = 0.2))
polygon(x = c(names(periods_an), rev(names(periods_an))),
        y = c(an_summary$`75plus`[1,], rev(an_summary$`75plus`[3,])),
        border = NA, col = adjustcolor(col_plot_a["75plus"], alpha.f = 0.2))

# Plot point-estimates
lines(periods_an,
      an_summary$total[2,], lty = 1, lwd = 2, col = col_plot_a["total"])
lines(periods_an,
      an_summary$`00_74`[2,], lty = 2, lwd = 2, col = col_plot_a["00_74"])
lines(periods_an,
      an_summary$`75plus`[2,], lty = 2, lwd = 2, col = col_plot_a["75plus"])

# Horizontal line at h = 0
abline(h = 0)

# Add legend
legend("topleft", 
       c(title_plot_a[[1]],
         title_plot_a[[2]],
         title_plot_a[[3]]), 
       lty = c(1, 2, 2), lwd = 2,
       col = col_plot_a,
       cex = 0.75)

# b) CLIMATE AND DEMOGRAPHIC CONTRIBUTIONS IN DIFFERENT PERIODS
par(mar = c(6, 4, 4, 2))  # Bottom, Left, Top, Right

# Create the barplot and capture the bar midpoints (x positions)
bar_centers <- barplot(values_period,
                       ylim = c(0, max(max_ci_period) + 2),
                       las = 2,
                       col = c(rep(col_plot_b[1], 3), rep(col_plot_b[2], 3)),
                       ylab = "Yearly heat-related deaths",
                       main = "")
title("b) Aggregated 21-year period impacts", 
      line = 0.65, font.main = 1, cex.main = 1)

# Add vertical error bars using arrows()
arrows(x0 = bar_centers, y0 = min_ci_period,
       x1 = bar_centers, y1 = max_ci_period,
       angle = 90, code = 3, length = 0.1, col = "black")

# Add a legend
legend("topleft", 
       legend = c(expression(paste("Global warming level (2", degree, "C)")), 
                  "End-of-century"), 
       fill = col_plot_b, bty = "n", cex = 0.75)

# COMMON LEGENDS TO BOTH PLOTS
mtext("Heat-related mortality projetions (London, SSP2-4.5)", 
      side = 3, outer = TRUE, line = -2.2, cex = 1.3, font = 2)

dev.off()

rm(bar_centers, datasets_period, title_plot_a, col_plot_a, col_plot_b,
   max_ci_period, min_ci_period, ncol.fig, nrow.fig, periods_an, values_period,
   compute_an_warming, an_summary)

#### FIGURE S1. PROCESS GRIDDED TEMPERATURES ###################################

# Load shapefile city of London
shp_london <- st_read("indata/shapefile_london/London_GLA_Boundary.shp")
shp_london <- st_transform(shp_london, 4326)

# Load the raster of temperature for one gcm and any year
file_path <- paste0("indata/tas/ssp245/", study_param$selected_gcms[1],
                    "/tas_day_", study_param$selected_gcms[1],
                    "_ssp245_r1i1p1f1_gn",
                    "_", 2100, ".nc")
raster_data <- brick(file_path)

# Some NetCDF files have longitudes ranging from 0 to 360 instead of -180 to 180.
# This adjustment ensures compatibility with the London shapefile
if (xmin(raster_data) > 180) {
  extent(raster_data) <- extent(xmin(raster_data) - 360,
                                xmax(raster_data) - 360,
                                ymin(raster_data),
                                ymax(raster_data))
}

# PLOT GRIDDED TEMPERATURES AND SHAPEFILE CITY
pdf(file = "outdata/plot/figs1_process_gridded_temperatures.pdf",
    width = 5, height = 4)
plot(raster_data[[226]]-273.15, 
     xlab = "Longitude", ylab = "Latitude",
     main = paste0("Mean temperature 2100-08-15 \n(SSP2-RCP4.5, ",
                   study_param$selected_gcms[1], ")"))
plot(shp_london$geometry, add = TRUE)
dev.off()

rm(shp_london, file_path, raster_data)
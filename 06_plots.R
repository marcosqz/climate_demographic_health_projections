################################################################################

# This scripts generates figures using outputs from previous scripts

#### LOAD LIBRARIES ############################################################

library(dlnm) # onebasis
library(sf) # st_read
library(raster) # brick, xmin

#### LOAD DATA #################################################################

# Load study parameters
load("indata/processed/study_parameters.RData")

# Load data observations
load("indata/processed/data_obs_temp_mort.RData")
load("indata/processed/data_obs_popu.RData")

# Load data with global warming level periods
load("indata/processed/data_global_warming_levels_ssp245.RData")

# Load outputs from the epidemiological models
load("outdata/file/01_epi_model/argvar.RData")
load("outdata/file/01_epi_model/arglag.RData")
load("outdata/file/01_epi_model/dlnm_var.RData")
load("outdata/file/01_epi_model/coefsim_age.RData")
load("outdata/file/01_epi_model/mmt_age.RData")

# Load raw and calibrated demographic projections
load("indata/processed/data_proj_mort_popu_ssp2.RData")
load(paste0(
  "outdata/file/02_calibrated_demographic_projections/",
  "data_proj_mort_popu_calibrated_ssp2.RData"))

# Load raw, calibrated and constant temperature projections
load("indata/processed/data_proj_temp_ssp245.RData")
load(paste0("outdata/file/03_calibrated_climate_projections/",
            "data_proj_temp_biascorrection_ssp245.RData"))
load(paste0("outdata/file/03_calibrated_climate_projections/",
            "data_proj_temp_constant_ssp245.RData"))

# Load health impact results
load(file = paste0(
  "outdata/file/04_health_impacts/",
  "heat_related_mortality_year_full_ssp245.RData"))
load(file = paste0(
  "outdata/file/04_health_impacts/", 
  "heat_related_mortality_year_exclude_temp_ssp245.RData"))
load(paste0(
  "outdata/file/04_health_impacts/",
  "heat_related_mortality_year_exclude_demo_ssp245.RData"))
load(paste0(
  "outdata/file/04_health_impacts/",
  "heat_related_mortality_period_full_ssp245.RData"))
load(paste0(
  "outdata/file/04_health_impacts/",
  "heat_related_mortality_period_exclude_temp_ssp245.RData"))
load(paste0(
  "outdata/file/04_health_impacts/",
  "heat_related_mortality_period_exclude_demo_ssp245.RData"))

#### FIGURE 1. AGE-SPECIFIC TEMPERATURE-MORTALITY ASSOCIATIONS #################

# VISUALLITATION PARAMETERS BY AGE GROUPS
age_parameters <- data.frame(
  response = c("mort.00_74", "mort.75plus"),
  groups = c("00_74", "75plus"),
  col = c("#1B9E77", "#D95F02"))

title_plot <- c(
  "00_74" = expression("a) Young ("< 75*" years)"),
  "75plus" = expression("b) Old (">= 75*" years)"))

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
       xlab = c(4, 4, 0, 0),
       ylab = "Relative risk", ylim = c(1.000, ymax), 
       log = "y", type = "n")
  title(title_plot[i_age], line = 0.65, font.main = 1, cex.main = 1)
  
  # Plot sample of RR
  for(i in 1:study_param$n_sim){
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
proj_temp_constant_year <- 
  aggregate(.~ lubridate::year(date), data = proj_temp_constant, FUN = mean)

# Arrange yearly temperature projections datasets
proj_temp_year$date <- NULL
proj_temp_bc_year$date <- NULL
proj_temp_constant_year$date <- NULL
colnames(proj_temp_year)[1] <- "year"
colnames(proj_temp_bc_year)[1] <- "year"
colnames(proj_temp_constant_year)[1] <- "year"

# Keep only years above 1990 for the plot
proj_temp_year <- subset(proj_temp_year, year >= 1990)
proj_temp_bc_year <- subset(proj_temp_bc_year, year >= 1990)
proj_temp_constant_year <- subset(proj_temp_constant_year, year >= 1990)

# Aggregate yearly temperature observations
data_temp_year <- data_tempmort[, c("date", "tmean")]
data_temp_year$year <- lubridate::year(data_temp_year$date)
data_temp_year$date <- NULL
data_temp_year <- aggregate(.~year, data = data_temp_year, FUN = mean)
data_temp_year <- data_temp_year[data_temp_year$year != 2012,]

# PLOT PARAMETERS
col_plot <- c("#1F77B4", "#2CA02C", "#E73F74", "#6699CC")

# PLOT FIGURE WITH CLIMATE MODELS
pdf(file = "outdata/plot/fig3_climate_models.pdf",
    width = 7, height = 12)
par(mfrow = c(3, 1), mar = c(3.8, 6, 4.1, 0.5), las = 1)

lapply(seq_along(study_param$selected_gcms), function(i_loop) {
  
  i_gcm <- study_param$selected_gcms[i_loop]
  
  # TIME-SERIES TEMPERATURE PROJECTIONS FOR GCM
  # Raw projections
  plot(proj_temp_year$year, proj_temp_year[[paste0("temp.", i_gcm)]], 
       col = col_plot[1], type = "l",
       ylim = range(proj_temp_year[[paste0("temp.", i_gcm)]],
                    proj_temp_bc_year[[paste0("temp.", i_gcm)]],
                    proj_temp_constant_year[[paste0("temp.", i_gcm)]],
                    data_temp_year$tmean),
       xlab = "Year", ylab = expression(paste("Temperature (", degree, "C)")),
       main = paste0(letters[i_loop], ") GCM ", i_gcm," (London, SSP2-4.5)"),
       cex.main = 2,
       cex = 1.2,
       cex.lab = 1.5,
       cex.axis = 1.5,
       lwd = 2, lty = 2)
  # Bias-corrected projections
  lines(proj_temp_bc_year$year, proj_temp_bc_year[[paste0("temp.", i_gcm)]], 
        col = col_plot[2], lwd = 2)
  # Constant projections
  lines(proj_temp_constant_year$year, 
        proj_temp_constant_year[[paste0("temp.", i_gcm)]], 
        col = col_plot[3], lwd = 2)
  # Observations
  lines(data_temp_year$year, data_temp_year$tmean,
        col = col_plot[4], lwd = 4)
  
  # Add legend
  legend("bottomright",
         c("Observations",
           "Raw projections",
           "Bias-corrected projections",
           "Constant projections"),
         lwd = c(4, 2, 2, 2),
         lty = c(1, 2, 1, 1),
         col = col_plot[c(4, 1, 2, 3)],
         bg = "white",
         cex = 1.2)
  
})

dev.off()

rm(proj_temp_year, proj_temp_bc_year, proj_temp_constant_year, data_temp_year,
   col_plot)

#### FIGURE 4. HEAT-RELATED MORTALITY ##########################################

# CALCULATE SUMMARY VALUES OF HEAT-RELATED MORTALIY FOR EACH AGE GROUP

# Define the periods to plot
periods_an <- 1950:2099
names(periods_an) <- 1950:2099

# Compute median and confidence interval of the attributable number (an) by
# age group an period
an_summary <- lapply(c("an", paste0("an.", study_param$age_groups)), function(var) {
  sapply(periods_an, function(i_year) {
    quantile(an_year_full[an_year_full$year %in% i_year, var], 
             c(0.025, 0.5, 0.975))
  })
})
names(an_summary) <- c("total", study_param$age_groups)

# Put in a list all the datasets of 21-years period
datasets_period <- list(an_period_exclude_demo$gwl,
                        an_period_exclude_temp$gwl,
                        an_period_full$gwl,
                        an_period_exclude_demo$end_century,
                        an_period_exclude_temp$end_century,
                        an_period_full$end_century)

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

# Create the plot with the point-estimates
plot(
  x = seq_along(values_period),
  y = values_period,
  xlim = c(min(seq_along(values_period)) - 0.2, 
           max(seq_along(values_period)) + 0.2),
  ylim = c(0, max(max_ci_period) + 2),
  xaxt = "n",
  pch = 21,
  cex = 2,
  bg = c(rep(col_plot_b[1], 3), rep(col_plot_b[2], 3)),
  xlab = "",
  ylab = "Yearly heat-related deaths",
  main = "")
abline(h = 0)
abline(v = mean(seq_along(values_period)))

# Add custom x-axis labels
axis(1, at = seq_along(values_period), labels = names(values_period), las = 2)

title("b) Aggregated 21-year period impacts", 
      line = 0.65, font.main = 1, cex.main = 1)

# Add vertical error bars using arrows()
arrows(x0 = seq_along(values_period), y0 = min_ci_period,
       x1 = seq_along(values_period), y1 = max_ci_period,
       angle = 90, code = 3, length = 0.1, col = "black")

# Add a legend
legend("topleft", 
       legend = c(expression(paste("Global warming level (2", degree, "C)")), 
                  "End-of-century"), 
       pch = 21,
       pt.bg = col_plot_b, bty = "n", cex = 0.75)

# COMMON LEGENDS TO BOTH PLOTS
mtext("Heat-related mortality projetions (London, SSP2-4.5)", 
      side = 3, outer = TRUE, line = -2.2, cex = 1.3, font = 2)

dev.off()

rm(datasets_period, title_plot_a, col_plot_a, col_plot_b,
   max_ci_period, min_ci_period, ncol.fig, nrow.fig, periods_an, values_period,
   an_summary)

#### FIGURE S1. PROCESS GRIDDED TEMPERATURES ###################################

# Load shapefile city of London
shp_london <- st_read(
  "indata/raw/shapefile_london/London_GLA_Boundary.shp")
shp_london <- st_transform(shp_london, 4326)

# Select target dates
target_dates <- as.Date(c("2073-07-17", "2073-07-18", "2073-07-19"))

# Extract raster data for each GCM and date
raster_data <- lapply(study_param$selected_gcms, function(i_gcm) {
  
  # Load the raster of temperature for one gcm and any year
  file_path <- paste0("indata/raw/temperature_projections/", i_gcm,
                      "/ssp245/proj_temp_grid_", i_gcm,
                      "_ssp245_", 2073, ".nc")
  raster_data <- brick(file_path)
  
  raster_data <- raster_data[[paste0("X", format(target_dates, "%Y.%m.%d"))]]
  
  # Some NetCDF files have longitudes ranging from 0 to 360 instead of -180 to 180.
  # This adjustment ensures compatibility with the London shapefile
  if (xmin(raster_data) > 180) {
    extent(raster_data) <- extent(xmin(raster_data) - 360,
                                  xmax(raster_data) - 360,
                                  ymin(raster_data),
                                  ymax(raster_data))
  }

  return(raster_data)

})
names(raster_data) <- study_param$selected_gcms

# Plot figure
nrow.fig <- length(study_param$selected_gcms); ncol.fig <- length(target_dates)
pdf(file = "outdata/plot/figs1_process_gridded_temperatures.pdf",
    width = 8, height = 6)
par(mfrow = c(nrow.fig, ncol.fig),
    mar = c(3, 3, 1, 5),   # inner margins
    oma = c(0, 3, 3, 0),   # outer margins
    mgp = c(1.5, 0.5, 0),
    xpd = NA)

lapply(study_param$selected_gcms, function(i_gcm) {
  lapply(target_dates, function(i_date) {
    
    # Format date to the name used to index raster
    date_raster <- paste0("X", format(i_date, "%Y.%m.%d"))
    
    # Plot raster
    plot(raster_data[[i_gcm]][[date_raster]] - 273.15,
         xlim = c(-3, 3),
         ylim = c(50.5, 53),
         xlab = "Longitude", ylab = "Latitude",
         main = "", legend = TRUE)
    
    # Add shapefile London
    plot(shp_london$geometry, add = TRUE)
    
    # Plot arrow
    x0 <- -0.1
    y0 <- 51.5
    x1 <- -1.6
    y1 <- 50.3
    arrows(x0, y0, x1, y1 + 0.15, length = 0.08, col = "red", lwd = 1.5)
    
    # Write mean temperature London
    temp_val <- round(
      proj_temp[proj_temp$date %in% i_date, paste0("temp.", i_gcm)], 1)
    text(x1, y1, paste0("London: ", format(temp_val, nsmall = 1), "°C"), 
         cex = 0.8)
    
  })
})

# Add column and row labels at the top and left respectively
mtext(target_dates, side = 3, line = 0.5, outer = TRUE, 
      at = (1:ncol.fig - 0.5) / 3 - 0.015)
mtext(study_param$selected_gcms, side = 2, line = 0.5, outer = TRUE, 
      at = 1 - (1:nrow.fig - 0.5) / 3)

dev.off()

rm(shp_london, target_dates, raster_data, ncol.fig, nrow.fig)

#### FIGURE S2. GLOBAL WARMING LEVEL PERIODS ####################################

# PLOT PARAMETERS

col_plot <- c("#118DFF", "#750985", "#C83D95")
levels_gwl <- c("1.5", "2", "3")
offsets_gwl <- c(-0.15, 0, 0.15)
points_gwl <- c(15, 16, 17)

# PLOT FIGURE WITH CLIMATE MODELS
pdf(file = "outdata/plot/figs2_global_warming_level_periods.pdf",
    width = 7, height = 5)
par(mfrow = c(1, 1), mar = c(3.8, 6, 4.1, 0.5), las = 1)

# b) PERIODS OF GLOBAL WARMING LEVELS BY GCM 
par(mar = c(5, 8, 3.5, 0.5))  # (bottom, left, top, right)
plot(data_gwl$year[1], 1 - 0.15,
     xlim = c(2000, 2100), ylim = c(0.5, 3.5), type = "n",
     main = "Global warming levels (SSP2-RCP4.5)",
     xlab = "Year", ylab = "", , yaxt = "n")
axis(2, at = 1:3, labels = unique(data_gwl$gcm), las = 1) 

# Loop through each GWL
for (i in 1:length(levels_gwl)) {
  data <- subset(data_gwl, warming_level == levels_gwl[i])
  
  # Plot points
  points(data$year, 1:length(data$year) + offsets_gwl[i],
         pch = points_gwl[i], col = col_plot[i], cex = 2)
  
  # Add segments (horizontal lines)
  segments(data$year - 10, 1:length(data$year) + offsets_gwl[i],
           data$year + 10, 1:length(data$year) + offsets_gwl[i], 
           col = col_plot[i], lwd = 2)
}; rm(data, i)
abline(h = seq(1.5, 18.5), col = "grey", lty = 2)

# Add legend below the title, outside the plot
legend("top", legend = paste0(levels_gwl, "ºC"), col = col_plot,
       pch = points_gwl, horiz = TRUE, cex = 1.2,
       bty = "n", xpd = TRUE, inset = c(0, -0.11))


dev.off()

rm(offsets_gwl, levels_gwl, points_gwl, col_plot)
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
load("outdata/file/01_epi_model/coef_age.RData")
load("outdata/file/01_epi_model/vcov_age.RData")
load("outdata/file/01_epi_model/coefsim_age.RData")
load("outdata/file/01_epi_model/mmt_age.RData")

# Load raw and calibrated demographic projections
load("indata/processed/data_proj_mort_popu_ssp2.RData")
load(paste0(
  "outdata/file/02_calibrated_demographic_projections/",
  "data_proj_mort_popu_calibrated_ssp2.RData"))
load(paste0(
  "outdata/file/02_calibrated_demographic_projections/",
  "data_proj_mort_popu_calibrated_daily_ssp2.RData"))
load(paste0(
  "outdata/file/02_calibrated_demographic_projections/",
  "data_proj_mort_popu_calibrated_constant_ssp2.RData"))

# Load raw, calibrated and constant temperature projections
load("indata/processed/data_proj_temp_ssp245.RData")
load(paste0("outdata/file/03_calibrated_climate_projections/",
            "data_proj_temp_biascorrection_ssp245.RData"))

# Load health impact results
load("outdata/file/04_health_impacts/heat_related_mortality_year_ssp245.RData")
load("outdata/file/04_health_impacts/heat_related_mortality_period_ssp245.RData")

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

# PLOT PARAMETERS
ymax <- 3

# PLOT FIGURE WITH AGE-SPECIFIC TEMPERATURE-MORTALITY ASSOCIATIONS
nrow.fig <- 1; ncol.fig <- 2
pdf("outdata/plot/fig1_age_specific_associations.pdf",
    width = ncol.fig*3*1.5, height = nrow.fig*3*1.5)
layout(matrix(seq(ncol.fig * nrow.fig), nrow = nrow.fig, byrow = TRUE))
par(mex = 0.8, mgp = c(3, 1, 0), las = 1, oma = c(0, 0, 0, 0))

lapply(study_param$age_groups, function(i_age) {
  
  # Exposure-response basis
  basis_temp <- onebasis(x_temp, 
                         fun = argvar$fun, 
                         knots = argvar$knots, 
                         Boundary.knots = argvar$Bound)
  
  # Predict the estimated association centred at the mmt
  cp <- crosspred(basis_temp,
                  coef = coef_age[[i_age]],
                  vcov = vcov_age[[i_age]],
                  model.link = "log",
                  at = x_temp,
                  cen = mmt_age[[i_age]])
  
  # Extract color for the corresponding age group
  col_age <- age_parameters$col[age_parameters$groups == i_age]
  
  # Plot curve
  plot(cp, lwd = 1, 
       col = "grey5", lty = 2,
       xlab = expression(paste("Temperature (", degree, "C)")),
       ylab = "Relative risk", 
       ylim = c(1.000, ymax), 
       log = "y")
  title(title_plot[i_age], line = 0.65, font.main = 1, cex.main = 1)
  
  # Highlight heat tail of the curve
  lines(
    cp$predvar[cp$predvar >=  mmt_age[[i_age]]],
    cp$allRRfit[cp$predvar >=  mmt_age[[i_age]]],
    col = col_age, lwd = 2)
  
  # Other lines in the plot
  abline(h = 1)
  abline(v = mmt_age[[i_age]], lwd = 1, lty = 3, col = "black")
  
})

# Title common to both panels
mtext("Age-specific temperature-mortality associations (London, 1990-2012)", 
      side = 3, outer = TRUE, line = -2.2, cex = 1.3, font = 2)

dev.off()

rm(age_parameters, title_plot, ncol.fig, nrow.fig, x_temp, ymax)

#### FIGURE 2. DEMOGRAPHIC PROJECTIONS #########################################

# Calculate yearly mortality observations
data_mort_year <- aggregate(
  cbind(mort.00_74, mort.75plus) ~ year, 
  data = data_tempmort, FUN = sum)
data_mort_year <- subset(data_mort_year, year != 2012)

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
nrow.fig <- 2; ncol.fig <- 3
pdf("outdata/plot/fig2_demographic_projections.pdf",
    width = ncol.fig*3.25, height = nrow.fig*3.25)
layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = nrow.fig, byrow = TRUE))
par(mex = 0.8, mgp = c(3, 1, 0), las = 1, oma = c(0, 0, 0, 0))

# a) PLOT PROJECTIONS OF MORTALITY

# Initialize plot
plot(x = 1, 
     y = 1,
     xlim = range(proj_mortpopu$year),
     ylim = c(0, max(proj_mortpopu$mort.00_74, 
                     proj_mortpopu$mort.75plus)/1e5),
     type = "n",
     main = "a) Mortality projections \n(SSP2 - UK and Northern Ireland)",
     xlab = "Time (years)", 
     ylab = "Number of deahts (x100,000)", 
     cex.lab = 1.1,
     xaxt = "n")

# Define tick positions
ticks <- seq(min(proj_mortpopu$year), max(proj_mortpopu$year), by = 5)
abline(v = ticks, col = "grey", lwd = 0.5, lty = 2)

# Add axis
axis(1, at = ticks)

# Plot mortality projection by age group
lapply(study_param$age_groups, function(i_age) {
  points(x = proj_mortpopu$year, 
        y = proj_mortpopu[[paste0("mort.", i_age)]]/1e5,
        col = "black",
        bg = col_plot[[i_age]],
        pch = 21)
})
abline(h = 0)

# Add legend
legend("topleft", 
       legend = c(title_plot[1], title_plot[2]),
       pch = 21,
       pt.bg = col_plot,
       col = "black",
       # col = "black", 
       cex = 1,
       bg = "white")

# b) PLOT PROJECTIONS OF POPULATION PERCERTANGE

# Initialize plot
plot(x = 1, 
     y = 1,
     xlim = range(proj_mortpopu$year),
     ylim = c(0, 100),
     type = "n", 
     main = "b) Populations projections \n(SSP2 - UK and Northern Ireland)", 
     xlab = "Time (years)", 
     ylab = "Percentage (%)", 
     cex.lab = 1.1,
     xaxt = "n")

# Define tick positions
ticks <- seq(min(proj_mortpopu$year), max(proj_mortpopu$year), by = 5)

# Add axis
axis(1, at = ticks)

# Add 1st age group from 0 to % population of that age group
polygon(x = c(proj_mortpopu$year, rev(proj_mortpopu$year)), 
        y = c(rep(0, length(proj_mortpopu$year)), 
              rev(proj_mortpopu$popu.00_74/proj_mortpopu$popu*100)), 
        col = col_plot[["00_74"]], 
        border = NA)

# Add last age age group from % of that age group to 100
polygon(x = c(proj_mortpopu$year, rev(proj_mortpopu$year)),
        y = c(100 - proj_mortpopu$popu.75plus/proj_mortpopu$popu*100, 
          rev(rep(100, length(proj_mortpopu$popu.75plus)))), 
        col = col_plot[["75plus"]], 
        border = NA)

# Add vertical lines to identify ticks
abline(v = ticks, col = "grey", lwd = 0.5, lty = 2)

# Add legend
legend("bottomleft", 
       legend = c(title_plot[[1]], title_plot[[2]]), 
       fill = col_plot, 
       border = "black", 
       bty = "o", 
       cex = 1,
       bg = "white")

# c) PLOT CALIBRATED MORTALITY AND POPULATION PROJECTIONS

# Adjust margins: bottom, left, top, right
par(mar = c(5, 4, 4, 5) + 0.1)

# Initialize population plot (left y-axis)
plot(x = 1,
     y = 1,
     xlim = range(proj_mortpopu_cal$year),
     ylim = range(c(proj_mortpopu_cal$popu.00_74, 
                    proj_mortpopu_cal$popu.75plus))/1e5,
     type = "n", 
     xlab = "Time (years)", 
     ylab = "Population (x100,000)", 
     cex.lab = 1.1,
     main = "c) Spatial calibration demographic \n projections (SSP2 - London)",
     xaxt = "n")

# Define tick positions
ticks <- seq(min(proj_mortpopu$year), max(proj_mortpopu$year), by = 5)
abline(v = ticks, col = "grey", lwd = 0.5, lty = 2)

# Add axis
axis(1, at = ticks)

# Loop age groups
lapply(study_param$age_groups, function(i_age) {
  
  # Plot calibrated population projections
  points(proj_mortpopu_cal$year, 
         proj_mortpopu_cal[[paste0("popu.", i_age)]]/1e5, 
         bg = col_plot_dark[i_age],
         col = "black",
         pch = 21)
  
  # Plot observed populations
  points(data_popu$year,
         data_popu[[paste0("popu.", i_age)]]/1e5,
         pch = 4,
         col = col_plot_dark[i_age],
         cex = 0.6)
  
  lines(data_popu$year, 
        data_popu[[paste0("popu.", i_age)]]/1e5,
        col = col_plot_dark[i_age], 
        lwd = 0.5)
  
})
abline(h = 0)

# Add second plot with mortality data (right y-axis)
par(new = TRUE)

# Initialize mortality plot (right y-axis)
plot(x = 1,
     y = 1,
     xlim = range(proj_mortpopu_cal$year),
     ylim = range(c(
       proj_mortpopu_cal[, paste0("mort.", study_param$age_groups)],
       data_mort_year[, paste0("mort.", study_param$age_groups)]))/1e5,
     axes = FALSE, 
     xlab = "", 
     ylab = "")

# Loop age groups

lapply(study_param$age_groups, function(i_age) {
  
  # Plot calibrated mortality projections
  points(proj_mortpopu_cal$year, 
         proj_mortpopu_cal[[paste0("mort.", i_age)]]/1e5, 
         col = "black", 
         bg = adjustcolor(col_plot_light[i_age], alpha.f = 1),
         pch = 21)
  
  # Plot mortality observations
  points(data_mort_year$year, 
         data_mort_year[[paste0("mort.", i_age)]]/1e5, 
         pch = 4, 
         col = adjustcolor(col_plot_light[i_age], alpha.f = 1), 
         cex = 0.6)
  lines(data_mort_year$year, 
        data_mort_year[[paste0("mort.", i_age)]]/1e5,
        col = adjustcolor(col_plot_light[i_age], alpha.f = 1), 
        lwd = 0.5)
  
})

# Add right y-axis for mortality
axis(side = 4)
mtext("Number of deahts (x100,000)", side = 4, line = 3, las = 0, cex = 0.75)

# Add legend
legend("right",
       c(expression(underline("Shapes")),
         "Projections",
         "Observations",
         expression(underline("Colours")),
         expression("Population "< 75*""),
         expression("Population ">= 75*""),
         expression("Mortality "< 75*""),
         expression("Mortality ">= 75*"")),
       col = c(NA, 1, 1, NA, col_plot_dark, col_plot_light),
       lty = c(NA, NA, 1, NA, 1, 1, 1, 1),
       pch = c(NA, 21, 4, NA, NA, NA, NA, NA),
       lwd = c(NA, NA, 1, NA, 7, 7, 7, 7),
       text.font = c(2, 1, 1, 2, 1, 1, 1, 1),
       pt.cex = 1,
       cex = 0.85,
       bty = "n")

# d) PLOT TEMPORAL CALIBRATION MORTALITY PROJECTIONS
par(mar = c(4.2, 4.2, 2, 1))

# Select only projection period
ind_plot <- proj_mortpopu_daily$year >= 2010

# Initialize mortality plot
plot(x = proj_mortpopu_daily$date[ind_plot],
     y = rep(1, length(proj_mortpopu_daily$date[ind_plot])),
     ylim = c(0, max(subset(proj_mortpopu_daily[ind_plot,], 
                            select = paste0("mort.", study_param$age_groups)))),
     type = "n",
     xlab = "Time (days)",
     ylab = "Number of deaths",
     main = "d) Temporal calibration mortality projections (SSP2 - London)",
     xaxs = "i",
     cex.lab = 1.1)
abline(h = 0)

# Plot daily mortality projections by age groups
lapply(study_param$age_groups, function(i_age) {
  
  # Calibrated daily mortality projections
  lines(proj_mortpopu_daily$date[ind_plot],
        proj_mortpopu_daily[[paste0("mort.", i_age)]][ind_plot],
        col = col_plot_light[i_age])
  
  # Fixed daily mortality projections
  lines(proj_mortpopu_constant$date[ind_plot],
        proj_mortpopu_constant[[paste0("mort.", i_age)]][ind_plot],
        col = col_plot_dark[i_age])
  
})

# Add legend
text_legend <- c(
  expression("Young ("< 75*" years)"),
  expression("Old (">= 75*" years)"),
  expression("Young ("< 75*" years) - fixed"),
  expression("Old (">= 75*" years) - fixed"))
legend("topleft", 
       legend = text_legend,
       col = c(col_plot_light, col_plot_dark),
       lty = 1,
       lwd = 1.5)
rm(text_legend)

dev.off()

rm(data_mort_year, col_plot, col_plot_dark, col_plot_light, ncol.fig, nrow.fig,
   title_plot, ind_plot, ticks)

#### FIGURE 3. CLIMATE PROJECTIONS #############################################

# Aggregate projection temperature to years
proj_temp_year <- 
  aggregate(.~ lubridate::year(date), data = proj_temp, FUN = mean)
proj_temp_bc_year <- 
  aggregate(.~ lubridate::year(date), data = proj_temp_bc, FUN = mean)

# Arrange yearly temperature projections datasets
proj_temp_year$date <- NULL
proj_temp_bc_year$date <- NULL
colnames(proj_temp_year)[1] <- "year"
colnames(proj_temp_bc_year)[1] <- "year"

# Calculate ensemble mean of temperature projections
proj_temp_bc_year$ens_mean <- rowMeans(
  proj_temp_bc_year[paste0("temp.", study_param$selected_gcms)])

# Keep only years above 1990 for the plot
proj_temp_year <- subset(proj_temp_year, year >= 1990)
proj_temp_bc_year <- subset(proj_temp_bc_year, year >= 1990)

# Aggregate yearly temperature observations
data_temp_year <- data_tempmort[, c("date", "tmean")]
data_temp_year$year <- lubridate::year(data_temp_year$date)
data_temp_year$date <- NULL
data_temp_year <- aggregate(.~year, data = data_temp_year, FUN = mean)
data_temp_year <- data_temp_year[data_temp_year$year != 2012,]

# PLOT PARAMETERS
col_plot <- c("#1B9E77", "#D95F02", "#7570B3")

# PLOT FIGURE WITH CLIMATE MODELS
pdf(file = "outdata/plot/fig3_climate_models.pdf",
    width = 14, height = 8)
par(mfrow = c(2, 2), 
    mar = c(3.8, 6, 4.1, 0.5),
    mgp = c(3, 1, 0),
    las = 1)

lapply(seq_along(study_param$selected_gcms), function(i_loop) {
  
  i_gcm <- study_param$selected_gcms[i_loop]
  
  # TIME-SERIES TEMPERATURE PROJECTIONS FOR GCM
  # Raw projections
  plot(x = proj_temp_year$year, 
       y = proj_temp_year[[paste0("temp.", i_gcm)]], 
       type = "l",
       col = col_plot[1],
       ylim = range(proj_temp_year[[paste0("temp.", i_gcm)]],
                    proj_temp_bc_year[[paste0("temp.", i_gcm)]],
                    data_temp_year$tmean),
       xlab = "Year", 
       ylab = expression(paste("Temperature (", degree, "C)")),
       main = "",
       lty = 2,
       lwd = 2, 
       cex = 1.2,
       cex.lab = 1.5,
       cex.axis = 1.5)
  title(paste0(letters[i_loop], ") GCM ", i_gcm), 
        line = 0.5, cex.main = 1.3)
  # Bias-corrected projections
  lines(proj_temp_bc_year$year, proj_temp_bc_year[[paste0("temp.", i_gcm)]], 
        col = col_plot[2], lwd = 2)
  # Observations
  lines(data_temp_year$year, data_temp_year$tmean,
        col = col_plot[3], lwd = 4)
  
  # Add legend
  legend("bottomright",
         c("Observations",
           "Raw projections",
           "Bias-corrected projections"),
         lwd = c(4, 2, 2),
         lty = c(1, 2, 1),
         col = col_plot[c(3, 1, 2)],
         bg = "white",
         cex = 1)
  
})

# Plot one more panel with the ensemble mean of the GCMs
plot(x = proj_temp_bc_year$year, 
     y = proj_temp_bc_year$ens_mean, 
     type = "l",
     col = col_plot[2],
     xlab = "Year", 
     ylab = expression(paste("Temperature (", degree, "C)")),
     main = "",
     lty = 1,
     lwd = 2, 
     cex = 1.2,
     cex.lab = 1.5,
     cex.axis = 1.5)
title(paste0(letters[length(study_param$selected_gcms)+1], 
             ") Ensemble mean GCMs"), 
      line = 0.5, cex.main = 1.3)

mtext("Temperature projections (London, SSP2-4.5)", 
      side = 3, outer = TRUE, line = -1.75, cex = 2, font = 2)

dev.off()

rm(proj_temp_year, proj_temp_bc_year, data_temp_year, col_plot)

#### FIGURE 4. HEAT-RELATED MORTALITY ##########################################

# CALCULATE YEARLY VALUES OF DEMOGRPAHIC AND CLIMATE PROJECTIONS
proj_mortpopu_yearly <- aggregate(
  cbind(mort.00_74, mort.75plus) ~ year, 
  data = proj_mortpopu_daily, 
  FUN = "sum")

# Aggregate projection temperature to years
proj_temp_bc_year <- 
  aggregate(.~ lubridate::year(date), data = proj_temp_bc, FUN = mean)
proj_temp_bc_year$ens_mean <- rowMeans(
  proj_temp_bc_year[paste0("temp.", study_param$selected_gcms)])

# Arrange yearly temperature projections datasets
proj_temp_bc_year$date <- NULL
colnames(proj_temp_bc_year)[1] <- "year"

# CALCULATE SUMMARY VALUES OF HEAT-RELATED MORTALIY FOR EACH AGE GROUP

# Define the periods to plot (panels a and b)
periods_an <- 1950:2099
names(periods_an) <- 1950:2099

# Compute median and confidence interval of the attributable number (an) by
# age group an period

# TODO: DECIDE AND KEEP ONLY ONE OPTION
mode_impact_est <- "median"
# mode_impact_est <- "coef_ensmean"

an_summary <- lapply(
  c("an", paste0("an.", study_param$age_groups)), function(var) {
  
  # Loop by years
      
  sapply(periods_an, function(i_year) {
    
    if(mode_impact_est == "median") {
      data_unc <- an_year$epi_sim.clim_gcm.full_democlim
      
      an_values <- quantile(
        data_unc[data_unc$year %in% i_year, var], 
        c(0.5, 0.025, 0.975))
      
    }
    
    if(mode_impact_est == "coef_ensmean") {
      data_est <- an_year$epi_est.clim_ensmean.full_democlim
      data_unc <- an_year$epi_sim.clim_gcm.full_democlim
      
      an_values <- c(
        subset(data_est, year == i_year, select = var)[,1],
        quantile(data_unc[data_unc$year %in% i_year, var], c(0.025, 0.975)))
    }
    
    return(an_values)
    
  })
})
names(an_summary) <- c("total", study_param$age_groups)

# Extract point_estimate and ci of 21-years periods

values_period <- lapply(1:3, function(i) {
  c(
    "Demographic" = an_period$gwl$only_demo[i,"an"],
    "Climate" = an_period$gwl$only_clim[i,"an"],
    "Climate + \ndemographic" = an_period$gwl$full_democlim[i,"an"],
    "Demographic" = an_period$end_century$only_demo[i,"an"],
    "Climate" = an_period$end_century$only_clim[i,"an"],
    "Climate + \ndemographic" = an_period$end_century$full_democlim[i,"an"]) / 21
})
names(values_period) <- c("point_estimate", "ci_low", "ci_high")
  
# PLOT PARAMETERS
col_plot <- list(
  gcm = c("#72190E", "#332288", "#225555"),
  age = c("total" = "#000000", "00_74" = "#1B9E77", "75plus" = "#D95F02"),
  period = c("#C969A1", "#CE4441")
)

title_plot <- list(
  a = "a) Time trends in age-specific impacts", 
  b = "b) Aggregated 21-year period impacts",
  c = "c) Demographic and climate projections")

legend_plot <- list(
  a = c(
    "Total (Old age)",
    expression("Young ("< 75*" years)"),
    expression("Old (">= 75*" years)")),
  b = c(
    expression(paste("Global warming level (2", degree, "C)")), 
    "End-of-century"),
  c = c(
    expression("Mortality (">= 75*" years)"),
    "Temperature (ensemble mean)",
    paste0("Temperature (GCM", 1:3, ")"),
    expression("Mortality ("< 75*" years)")))

# PLOT HEALT-RELATED IMPACTS
nrow.fig <- 2; ncol.fig <- 2
pdf("outdata/plot/fig4_health_impacts.pdf",
    width = ncol.fig*3*1.5, height = nrow.fig*3*1.5)
layout(matrix(c(1, 2, 3, 3), nrow = nrow.fig, byrow = TRUE))

# a) YEARLY TIME-SERIES HEAT-RELATED MORTALITY BY AGE GROUPS
par(mex = 1, mar = c(5, 4, 4, 1), 
    mgp = c(3, 1, 0),
    oma = c(0, 0, 0, 0),
    las = 1)

plot(names(periods_an),
     an_summary$total[1,],
     ylim = c(0, max(an_summary$total[1,], na.rm = TRUE)), type = "n",
     xlab = "Time (years)",
     ylab = "Yearly heat-related deaths",
     main = "")
title(title_plot$a, 
      line = 0.65, font.main = 2, cex.main = 1.3)

# Plot 95% CI
lapply(names(an_summary), function(i_age) {
  polygon(x = c(names(periods_an), rev(names(periods_an))),
          y = c(an_summary[[i_age]][2,], rev(an_summary[[i_age]][3,])),
          border = NA, 
          col = adjustcolor(col_plot$age[i_age], alpha.f = 0.2))
})

# Define the type of line for each age group
panel_lty <- c("total" = 1,
               "00_74" = 2,
               "75plus" = 2)

# Plot point-estimates
lapply(names(an_summary), function(i_age) {
  lines(x = periods_an,
        y = an_summary[[i_age]][1,], 
        lty = panel_lty[i_age], 
        lwd = 2, 
        col = col_plot$age[i_age])
})

# Horizontal line at h = 0
abline(h = 0)

# Add legend
legend("topleft", 
       legend_plot$a, 
       lty = panel_lty, 
       lwd = 2,
       col = col_plot$age,
       cex = 1)

# b) CLIMATE AND DEMOGRAPHIC CONTRIBUTIONS IN DIFFERENT PERIODS
par(mar = c(6, 4, 4, 2))  # Bottom, Left, Top, Right

# Create the plot with the point-estimates
plot(
  x = seq_along(values_period$point_estimate),
  y = values_period$point_estimate,
  xlim = c(min(seq_along(values_period$point_estimate)) - 0.2, 
           max(seq_along(values_period$point_estimate)) + 0.2),
  ylim = c(0, max(values_period$ci_high) + 2),
  xaxt = "n",
  pch = 21,
  cex = 2,
  bg = c(rep(col_plot$period[1], length(values_period$point_estimate)/2), 
         rep(col_plot$period[2], length(values_period$point_estimate)/2)),
  xlab = "",
  ylab = "Yearly heat-related deaths",
  main = "")
abline(h = 0)
abline(v = mean(seq_along(values_period$point_estimate)))

# Add custom x-axis labels
axis(1, 
     at = seq_along(values_period$point_estimate), 
     labels = names(values_period$point_estimate), las = 2)

title(title_plot$b, 
      line = 0.65, font.main = 2, cex.main = 1.3)

# Add vertical error bars using arrows()
arrows(x0 = seq_along(values_period$point_estimate), 
       y0 = values_period$ci_low,
       x1 = seq_along(values_period$point_estimate), 
       y1 = values_period$ci_high,
       angle = 90, code = 3, length = 0.1, col = "black")

# Add a legend
legend("topleft", 
       legend = legend_plot$b, 
       pch = 21,
       bg = "white",
       pt.bg = col_plot$period, 
       cex = 1)

# c) DEMOGRAPHIC AND CLIMATE PROJECTIONS

# Adjust margins: bottom, left, top, right
par(mar = c(4, 4, 4, 5) + 0.1,
    mex = 1)

# Initialize plot
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(2010, 2100),
     ylim = c(0, max(proj_mortpopu_yearly$mort.75plus)/1e5),
     xlab = "Time (years)",
     ylab = "Number of deahts (x100,000)")
title(title_plot$c, 
      line = 0.65, font.main = 2, cex.main = 1.3)
abline(h = 0)

# Plot yearly mortality projections
lapply(c("00_74", "75plus"), function(i_age) {
  lines(x = proj_mortpopu_yearly$year,
        y = proj_mortpopu_yearly[[paste0("mort.", i_age)]]/1e5,
        col = col_plot$age[i_age], 
        lwd = 3)
})

# Subset to select warming level
period_selected_gwl <- subset(
  data_gwl, warming_level == study_param$selected_warming)

# Plot arrows defining the GWL for each GCM
lapply(1:length(study_param$selected_gcms), function(i) {
  
  # Extract the years for the specific GWL period
  years <- subset(period_selected_gwl, 
                  gcm == study_param$selected_gcms[i],
                  select = c("year1", "year2"))
  
  # Draw a transparent polygon with a rectangle with the GWL period
  polygon(
    x = c(years[1], years[2], years[2], years[1]),
    y = c(-0.1, -0.1, 1, 1),
    col = adjustcolor(col_plot$gcm[i], alpha.f = 0.2),
    border = NA)
  
  # Draw horizontal lines defining the GWL period
  arrows(
    x0 = years[[1]], 
    x1 = years[[2]], 
    y0 = 0.75- i*0.02, 
    lwd = 1.5,
    code = 3, 
    length = 0.1, 
    col =  col_plot$gcm[i])
  
  # Put text specifying to which GCM correspond each GWL period
  text(
    paste0("GWL for GCM", i),
    x = 2025,
    y = 0.75-i*0.02,
    col = col_plot$gcm[i],
    cex = 0.8)
  
})

# Add second plot with climate projections (right y-axis)
par(new = TRUE)

# Initialize plot
plot(
  x = 1,
  y = 1,
  xlim = c(2010, 2100),
  ylim = range(proj_temp_bc_year[proj_temp_bc_year$year>=2010,-1]),
  type = "n", 
  axes = FALSE, 
  xlab = "", 
  ylab = "", 
  col = "black")

# Plot yearly temperature in each GCM
lapply(1:length(study_param$selected_gcms), function(i) {
  lines(
    x = proj_temp_bc_year$year,
    y = proj_temp_bc_year[[paste0("temp.",study_param$selected_gcms[i])]],
    col = col_plot$gcm[i], 
    lwd = 1.5, 
    lty = 1)
})

# Plot ensemble mean temperatures
lines(
  x = proj_temp_bc_year$year,
  y = proj_temp_bc_year$ens_mean,
  lwd = 4,
  col = "black")

# Add right y-axis for temperature projections
axis(side = 4)
mtext(
  expression(paste("Temperature (", degree, "C)")), 
  side = 4, 
  line = 3, 
  las = 0, 
  cex = 0.8)

# Add legend
legend(
  "bottomright",
  legend_plot$c,
  lty = c(1, 1, rep(1, length(study_param$selected_gcms)), 1),
  lwd = c(2.5, 2.5, rep(2.5, length(study_param$selected_gcms)), 2.5),
  col = c(col_plot$age["75plus"],
          "black",
          col_plot$gcm,
          col_plot$age["00_74"]),
  cex = 1,
  title = "Projections",
  bg = "white")

dev.off()

rm(proj_mortpopu_yearly, proj_temp_bc_year, periods_an, 
   an_summary, values_period, col_plot, title_plot, legend_plot, nrow.fig, 
   ncol.fig, panel_lty, period_selected_gwl)

#### FIGURE S1. PROCESS GRIDDED TEMPERATURES ###################################

# Load shapefile city of London
shp_london <- st_read("indata/raw/shapefile_london/London_GLA_Boundary.shp")
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
    mar = c(3, 3, 1, 5),
    oma = c(0, 3, 3, 0),
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

#### FIGURE S3. EXTRAPOLATION TEMPERATURE-MORTALTITY ASSOCIATION ###############

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

# PLOT FIGURE WITH AGE-SPECIFIC TEMPERATURE-MORTALITY ASSOCIATIONS
nrow.fig <- 1; ncol.fig <- 2
pdf("outdata/plot/figs3_extrapolation_age_specific_associations.pdf",
    width = ncol.fig*3*1.5, height = nrow.fig*3*1.5)
layout(matrix(seq(ncol.fig * nrow.fig), nrow = nrow.fig, byrow = TRUE))
par(mex = 0.8, mgp = c(3, 1, 0), las = 1, oma = c(0, 0, 0, 0))

lapply(study_param$age_groups, function(i_age) {
  
  # Exposure-response basis
  basis_temp <- onebasis(x_temp, 
                         fun = argvar$fun, 
                         knots = argvar$knots, 
                         Boundary.knots = argvar$Bound)
  
  # Predict the estimated association centred at the mmt
  cp <- crosspred(basis_temp,
                  coef = coef_age[[i_age]],
                  vcov = vcov_age[[i_age]],
                  model.link = "log",
                  at = x_temp,
                  cen = mmt_age[[i_age]])
  
  # Extract color for the corresponding age group
  col_age <- age_parameters$col[age_parameters$groups == i_age]
  
  # Plot curve
  plot(cp,
       type = "n",
       xlab = expression(paste("Temperature (", degree, "C)")),
       ylab = "Relative risk", 
       ylim = c(1.000, ymax),
       log = "y")
  title(title_plot[i_age], line = 0.65, font.main = 1, cex.main = 1)
  
  # Highlight heat tail of the curve
  
  ind_cold <- cp$predvar <= mmt_age[[i_age]]
  ind_heat_obs <- (cp$predvar >=  mmt_age[[i_age]]) &
    (cp$predvar <=  x_temp["100.0%"])
  ind_heat_proj <- (cp$predvar >=  x_temp["100.0%"])
  
  lines(
    x = cp$predvar[ind_cold],
    y = cp$allRRfit[ind_cold],
    col = "grey5", 
    lty = 2)
  lines(
    x = cp$predvar[ind_heat_obs],
    y = cp$allRRfit[ind_heat_obs],
    col = col_age,
    lwd = 2,
    lty = 1)
  lines(
    x = cp$predvar[ind_heat_proj],
    y = cp$allRRfit[ind_heat_proj],
    col = col_age,
    lwd = 2,
    lty = 2)
  
  # Other lines in the plot
  abline(h = 1)
  abline(v = x_temp["100.0%"], lwd = 1, lty = 2, col = "black")
  abline(v = mmt_age[[i_age]], lwd = 1, lty = 3, col = "black")
  
})

# Title common to both panels
mtext("Extrapolation age-specific temperature-mortality associations (London, 1990-2012)", 
      side = 3, outer = TRUE, line = -2.2, cex = 1.3, font = 2)

dev.off()

rm(age_parameters, title_plot, ncol.fig, nrow.fig, x_temp, ymax)

#### FIGURE S4. EXTRAPOLATION TEMPERATURE-MORTALTITY ASSOCIATION ###############

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

# PLOT FIGURE WITH AGE-SPECIFIC TEMPERATURE-MORTALITY ASSOCIATIONS
nrow.fig <- 1; ncol.fig <- 2
pdf("outdata/plot/figs4_simulations_age_specific_associations.pdf",
    width = ncol.fig*3*1.5, height = nrow.fig*3*1.5)
layout(matrix(seq(ncol.fig * nrow.fig), nrow = nrow.fig, byrow = TRUE))
par(mex = 0.8, mgp = c(3, 1, 0), las = 1, oma = c(0, 0, 0, 0))

lapply(study_param$age_groups, function(i_age) {
  
  # Exposure-response basis at the age-specific mmt
  cenvec <- onebasis(
    mmt_age[[i_age]],
    fun = argvar$fun,
    knots = argvar$knots,
    Boundary.knots = argvar$Bound)
  
  # Centered exposure-response basis at each daily projected temperature
  bcen <- scale(
    onebasis(x_temp, 
             fun = argvar$fun, 
             knots = argvar$knots, 
             Boundary.knots = argvar$Bound), 
    center = cenvec, 
    scale = FALSE)
  
  # Relative risks samples
  rrsim <- exp(bcen %*% coefsim_age[[i_age]])
  
  # Extract color for the correspondin age group
  col_age <- age_parameters$col[age_parameters$groups == i_age]
  
  # Plot empty plot
  plot(x_temp, 
       rrsim[,1], 
       type = "n",
       xlab = expression(paste("Temperature (", degree, "C)")),
       ylab = "Relative risk", 
       ylim = c(1.000, ymax), 
       log = "y")
  title(title_plot[i_age], line = 0.65, font.main = 1, cex.main = 1)
  
  # Plot sample of RR
  ind_cold <- x_temp <= mmt_age[[i_age]]
  ind_heat_obs <- (x_temp >=  mmt_age[[i_age]]) &
    (x_temp <=  x_temp["100.0%"])
  ind_heat_proj <- (x_temp >=  x_temp["100.0%"])
  
  for(i in 1:study_param$n_sim){
    lines(x_temp[ind_cold], 
          rrsim[ind_cold, i], 
          col = "grey5", 
          lty = 2, 
          lwd = 0.5)
    lines(x_temp[ind_heat_obs], 
          rrsim[ind_heat_obs, i], 
          col = col_age, 
          lwd = 0.5)
    lines(x_temp[ind_heat_proj], 
          rrsim[ind_heat_proj, i], 
          col = col_age,
          lty = 2,
          lwd = 0.5)
  }; rm(i)
  
  # Other lines in the plot
  abline(h = 1)
  abline(v = x_temp["100.0%"], lwd = 1, lty = 2, col = "black")
  abline(v = mmt_age[[i_age]], lwd = 1, lty = 3, col = "black")
  
})

# Title common to both panels
mtext("Simulations of age-specific temperature-mortality associations (London, 1990-2012)", 
      side = 3, outer = TRUE, line = -2.2, cex = 1.3, font = 2)

dev.off()

rm(age_parameters, title_plot, ncol.fig, nrow.fig, x_temp, ymax)

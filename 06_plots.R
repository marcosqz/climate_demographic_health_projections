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
load("outdata/file/01_epi_model/mmt_age.RData")

# Load raw, calibrated and constant temperature projections
load("indata/processed/data_proj_temp_ssp245.RData")
load(paste0("outdata/file/02_calibrated_climate_projections/",
            "data_proj_temp_biascorrection_ssp245.RData"))

# Load attributable fractions
load(paste0(
  "outdata/file/02_calibrated_climate_projections/attributable_fraction_",
  study_param$ssp_rcp_scenario,".RData"))

# Load raw and calibrated demographic projections
load("indata/processed/data_proj_mort_popu_ssp2.RData")
load(paste0(
  "outdata/file/03_calibrated_demographic_projections/",
  "data_proj_mort_popu_spatialcal_ssp2.RData"))
load(paste0(
  "outdata/file/03_calibrated_demographic_projections/",
  "data_proj_mort_popu_spatialcal_tempcal_ssp2.RData"))

# Load attributable numbers
load(paste0(
  "outdata/file/03_calibrated_demographic_projections/attributable_number_", 
  study_param$ssp_rcp_scenario, ".RData"))

# Load health impact summary results
load(paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_decades_",
  study_param$ssp_rcp_scenario, ".RData"))
load(paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_endcentury_",
  study_param$ssp_rcp_scenario,".RData"))
load(paste0(
  "outdata/file/04_health_impacts/heat_related_mortality_gwl_",
  study_param$ssp_rcp_scenario,".RData"))

#### FIGURE 1. AGE-SPECIFIC TEMPERATURE-MORTALITY ASSOCIATIONS #################

# Define the temperatures for the x-axis
pred_perc <- c(seq(0, 1, 0.1), 2:98, seq(99, 100, 0.1))
x_temp <- quantile(data_tempmort$tmean, pred_perc / 100)
rm(pred_perc)

# Consider also temperatures beyond the boundaries of observed temperatures
x_temp_extra <- seq(max(x_temp), max(x_temp) + 5, length = 11) 
x_temp <- c(x_temp, x_temp_extra[-1])
rm(x_temp_extra)

# Define plotting variables
ymax <- 5
col_plot <- c("00_74" = "#1B9E77", "75plus" = "#D95F02")
title_plot <- c(
  "00_74" = expression("a) Young ("< 75*" years)"),
  "75plus" = expression("b) Old (">= 75*" years)"))

# ---- PLOT AGE-SPECIFIC EXPOSURE-RESPONSE FUNCTIONS ----
nrow.fig <- 1; ncol.fig <- 2
pdf("outdata/plot/fig1_age_specific_associations.pdf",
    width = ncol.fig*3*1.5, height = nrow.fig*3*1.5)
layout(matrix(seq(ncol.fig * nrow.fig), nrow = nrow.fig, byrow = TRUE))
par(mex = 0.8, mgp = c(3, 1, 0), las = 1, oma = c(0, 0, 0, 0))

run_loop <- lapply(study_param$age_groups, function(i_age) {
  
  # Exposure-response basis
  basis_temp <- onebasis(x_temp, 
                         fun = argvar$fun, 
                         knots = argvar$knots, 
                         Boundary.knots = argvar$Bound)
  
  # Predict the estimated association centered at the mmt
  cp <- crosspred(basis_temp,
                  coef = coef_age[[i_age]][,"est"],
                  vcov = vcov_age[[i_age]],
                  model.link = "log",
                  at = x_temp,
                  cen = mmt_age[[i_age]])
  
  # Extract color for the corresponding age group
  col_age <- col_plot[i_age]
  
  # Plot curve
  plot(cp,
       col = rgb(1, 1, 1, 0), # line invisible
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
  
  lines(x_temp[ind_cold], 
        cp$allRRfit[ind_cold], 
        col = "grey5", 
        lty = 2, 
        lwd = 2)
  lines(x_temp[ind_heat_obs], 
        cp$allRRfit[ind_heat_obs], 
        col = col_age, 
        lwd = 2)
  lines(x_temp[ind_heat_proj], 
        cp$allRRfit[ind_heat_proj], 
        col = col_age,
        lty = 2,
        lwd = 2)
  
  # Other lines in the plot
  abline(h = 1)
  abline(v = x_temp["100.0%"], lwd = 1, lty = 2, col = "black")
  abline(v = mmt_age[[i_age]], lwd = 1, lty = 3, col = "black")
  
})

# Title common to both panels
mtext("Age-specific temperature-mortality associations (London, 1990-2012)", 
      side = 3, outer = TRUE, line = -2.2, cex = 1.3, font = 2)

dev.off()

rm(ymax, title_plot, col_plot, ncol.fig, nrow.fig, x_temp, run_loop)

#### FIGURE 2. CLIMATE PROJECTIONS AND ATTRIBUTABLE FRACTIONS ##################

# Create a list with annual mean raw and bias-corrected temperatures
proj_temp_year <- lapply(list(proj_temp, proj_temp_bc), function(data) {
  
  # Aggregate projection temperature to years
  proj_temp_year <- 
    aggregate(.~ lubridate::year(date), data = data, FUN = mean)
  
  # Arrange yearly temperature projections datasets
  proj_temp_year$date <- NULL
  colnames(proj_temp_year)[1] <- "year"
  
  # Keep only years above 1990 for the plot
  proj_temp_year <- subset(proj_temp_year, year >= 1990)
  
  return(proj_temp_year)
  
}); names(proj_temp_year) <- c("raw", "bc")

# Compute annual mean temperature observations
data_temp_year <- data_tempmort[, c("date", "tmean")]
data_temp_year$year <- lubridate::year(data_temp_year$date)
data_temp_year$date <- NULL
data_temp_year <- aggregate(.~year, data = data_temp_year, FUN = mean)
data_temp_year <- data_temp_year[data_temp_year$year != 2012,]

# Define plotting parameters
col_plot <- c("#88CCEE", "#AA4499", "#CC6677", "#332288", "#117733") 	
names(col_plot) <- c("raw", "observations", study_param$selected_gcms)

# Build a data frame with dates, corresponding years and decades (for time-based
# grouping)
data_time <- data.frame(
  date = proj_temp_bc$date,
  year = lubridate::year(proj_temp_bc$date),
  decade = lubridate::year(proj_temp_bc$date) %/% 10 * 10)

# Titles for plots with AFs
title_plot_af <- c(
  "00_74" = 
    expression(bold("d) Attributable Fraction - Young ("< 75*" years)")),
  "75plus" = 
    expression(bold("e) Attributable fraction - Old (">= 75*" years)")))

# Subset GWL dataset to the study GWL
period_selected_gwl <- subset(
  data_gwl, warming_level == study_param$selected_warming)

# ---- PLOT CLIMATE MODELS AND ATTRIBUTABLE FRACTIONS ----
pdf(file = "outdata/plot/fig2_climate_models.pdf", width = 12, height = 8)
layout_matrix <- matrix(c(
  1, 1, 2, 2, 3, 3,
  4, 4, 4, 5, 5, 5),
  nrow = 2, byrow = TRUE)
layout(layout_matrix)
par(mar = c(4.5, 6, 2.5, 0.5),
    mgp = c(3, 1, 0),
    las = 1)  

# ---- PANELS A-C: CLIMATE MODELS ----

# Define common ylim for all panels
ylim_plot <- sapply(
  proj_temp_year, function(x) x[,paste0("temp.", study_param$selected_gcms)])
ylim_plot <- range(ylim_plot, data_temp_year$tmean)

# Loop GCMs
run_loop <- lapply(seq_along(study_param$selected_gcms), function(i_loop) {
  
  i_gcm <- study_param$selected_gcms[i_loop]
  
  # Plot raw projections
  plot(x = proj_temp_year$raw$year, 
       y = proj_temp_year$raw[[paste0("temp.", i_gcm)]], 
       type = "l",
       lty = 1,
       lwd = 2, 
       cex = 1.2,
       cex.lab = 1.5,
       cex.axis = 1.5,
       col = col_plot[1],
       ylim = ylim_plot,
       xlab = "Year", 
       ylab = expression(paste("Temperature (", degree, "C)")),
       main = "")
  title(paste0(letters[i_loop], ") GCM ", i_gcm), 
        line = 1.0, cex.main = 2)
  
  # Draw a transparent polygon with a rectangle with the calibration period
  polygon(
    x = c(1990, 2011, 2011, 1990),
    y = c(-0.1, -0.1, 100, 100),
    col = adjustcolor("grey", alpha.f = 0.2),
    border = NA)

  # Draw horizontal lines defining the calibration period
  arrows(
    x0 = 1990,
    x1 = 2011,
    y0 = 15,
    lwd = 1.5,
    code = 3,
    length = 0.1,
    col = "grey")
  text(
    paste0("Calibration \n period"),
    x = 2001,
    y = 15.5,
    col = "black",
    cex = 1.25)
  
  # Plot bias-corrected projections
  lines(x = proj_temp_year$bc$year,
        y = proj_temp_year$bc[[paste0("temp.", i_gcm)]],
        col = col_plot[i_gcm],
        lwd = 2)
  
  # Plot temperature observations
  lines(x = data_temp_year$year,
        y = data_temp_year$tmean,
        col = col_plot[2],
        lwd = 4)
  
  # Extract the years for the specific GWL period
  years <- subset(period_selected_gwl,
                  gcm == i_gcm,
                  select = c("year1", "year2"))
  
  # Draw a transparent polygon with a rectangle with the GWL period
  polygon(
    x = c(years[1], years[2], years[2], years[1]),
    y = c(-0.1, -0.1, 100, 100),
    col = adjustcolor("#DDCC77", alpha.f = 0.2),
    border = NA)
  
  # Draw horizontal lines defining the GWL period
  arrows(
    x0 = years[[1]],
    x1 = years[[2]],
    y0 = 15,
    lwd = 1.5,
    code = 3,
    length = 0.1,
    col = "#DDCC77")
  text(
    paste0("GWL 2°C"),
    x = mean(c(years[[1]], years[[2]])),
    y = 15.5,
    col = "black",
    cex = 1.25)
  
  # Add legend
  legend("bottomright",
         c("Observations",
           "Raw projections",
           "Bias-corrected projections"),
         lwd = c(4, 2, 2),
         lty = c(1, 1, 1),
         col = c(col_plot[2], col_plot[1], col_plot[i_gcm]),
         bg = "white",
         cex = 1.25)
  
}); rm(run_loop, ylim_plot)

# ---- PANELS D-E: ATTRIBUTABLE FRACTIONS ----

# Loop age groups
run_loop <- lapply(study_param$age_groups, function(i_age) {
  
  # Initalizate plot
  plot(
    x = unique(data_time$decade),
    y = rep(1, length(unique(data_time$decade))),
    type = "n",
    ylim = c(0, 0.05)*100,
    xlab = "Decade",
    ylab = "Mean daily AF heat-related deaths (%)",
    xaxt = "n",
    main = "",
    cex = 1.2,
    cex.lab = 1.5,
    cex.axis = 1.5)
  title(title_plot_af[i_age], line = 1.5, cex.main = 2)
  
  # Plot x-axis labels with decades
  axis(1, 
       at = data_time$decade, 
       labels = paste0(data_time$decade, "s"), 
       cex.axis = 1.5)
  abline(h = 0)
  
  # Loop GCMs
  run_loop <- lapply(study_param$selected_gcms, function(i_gcm) {
    # Loop simulated epidemiological curves
    run_loop <- lapply(1:study_param$n_sim, function(i) {
      
      # Plot the AFs for the simulated coefficients from the epi models
      lines(
        x = unique(data_time$decade),
        y = sapply(
          split(af[[i_age]][[i_gcm]][,paste0("sim", i)], data_time$decade), 
          mean) * 100,
        col = col_plot[i_gcm],
        lty = 2,
        lwd = 0.25)
    })
    
    # Plot the AFs for the estimated coefficients from the epi models
    lines(
      x = unique(data_time$decade),
      y = sapply(split(af[[i_age]][[i_gcm]][,"est"], data_time$decade), 
                 mean) * 100,
      col = col_plot[i_gcm],
      lwd = 4)
  })
  
  # Add legend
  legend("topleft",
         c(study_param$selected_gcms,
           "Epidemiological estimates",
           "Epidemiological simulations"),
         col = c(col_plot[study_param$selected_gcms], 1, 1),
         lwd = c(3, 3, 3, 3, 1.25),
         lty = c(1, 1, 1, 1, 2),
         cex = 1.5)
  
}); rm(run_loop)

dev.off()

rm(proj_temp_year, data_temp_year, col_plot, data_time, title_plot_af, 
   period_selected_gwl, layout_matrix)

#### FIGURE 3. DEMOGRAPHIC PROJECTIONS AND ATTRIBUTABLE NUMBERS ################

# Calculate yearly mortality observations
data_mort_year <- aggregate(
  cbind(mort.00_74, mort.75plus) ~ year, 
  data = data_tempmort, FUN = sum)
data_mort_year <- subset(data_mort_year, year != 2012)

# Define colors for the age groups
col_plot <- c("#1B9E77", "#D95F02", "#CC6677", "#332288", "#117733") 	
names(col_plot) <- c(study_param$age_groups, study_param$selected_gcms)

# Build a data frame with dates, corresponding years and decades (for time-based
# grouping)
data_time <- data.frame(
  date = proj_temp_bc$date,
  year = lubridate::year(proj_temp_bc$date),
  decade = lubridate::year(proj_temp_bc$date) %/% 10 * 10)

# Title for panels with AN
title_plot_an <- c(
  "00_74" = expression(bold("d) AN - Young ("< 75*" years)")),
  "75plus" = expression(bold("e) AN - Old (">= 75*" years)")))

# ---- PLOT DEMOGRAPHIC PROJECTIONS AND ATTRIBUTABLE NUMBERS ----
pdf("outdata/plot/fig3_demographic_projections.pdf", width = 12, height = 8)

layout_matrix <- matrix(c(
  1, 1, 2, 2, 4, 4, 
  3, 3, 3, 3, 5, 5),
  nrow = 2, byrow = TRUE)
layout(layout_matrix)
par(mar = c(4.5, 6, 2.5, 0.5),
    mgp = c(3.5, 1, 0),
    las = 1)

# ---- PANELS A-B: SPATIAL CALIBRATION OF DEMOGRAPHIC PROJECTIONS ----

# Define plot-specific parameters for "mort" and "popu"
title_plot_cal <- c(
  "mort" = "a) Spatial calibration: mortality",
  "popu" = "b) Spatial cal.: population")
ylab_plot_cal <- c(
  "mort" = "Number of deahts (x100,000)",
  "popu" = "Population (x100,000)")
calperiod_plot <- list(
  "mort" = c(1990, 2011),
  "popu" = c(1992, 2019))
pos_text_cal <- list(
  "mort" = c(0.62, 0.54, 0.40, 0.15), # arrow, text_cal, cf1, cf2
  "popu" = c(40, 30, 57, 7))
text_cf <- list(
  "mort" = c("cf: x0.09", "cf: x0.10"), # young, old
  "popu" = c("cf: x0.13", "cf: x0.09"))

# Loop "mort" and "popu" variables
run_loop <- lapply(c("mort", "popu"), function(var) {
  
  if(var == "mort") {data_obs <- data_mort_year
  } else if (var == "popu") {data_obs <- data_popu}
  
  # Initialize plot
  plot(x = 1, 
       y = 1,
       xlim = range(proj_mortpopu_ldn$year),
       ylim = c(0,
                max(proj_mortpopu_ldn[[paste0(var, ".00_74")]], 
                    proj_mortpopu_ldn[[paste0(var, ".75plus")]],
                    unlist(data_obs
                           [,paste0(var, ".", study_param$age_groups)]))/1e5),
       type = "n",
       main = "",
       xlab = "Year", 
       ylab = ylab_plot_cal[var], 
       cex.lab = 1.5,
       cex.axis = 1.5,
       xaxt = "n")
  title(title_plot_cal[[var]], line = 1.0, cex.main = 2)
  
  # Draw a transparent polygon with a rectangle with the calibration period
  polygon(
    x = c(calperiod_plot[[var]][1], calperiod_plot[[var]][2], 
          calperiod_plot[[var]][2], calperiod_plot[[var]][1]),
    y = c(-1000, -1000, 1000, 1000),
    col = adjustcolor("grey", alpha.f = 0.2),
    border = NA)
  
  # Draw horizontal lines defining the calibration period
  arrows(
    x0 = calperiod_plot[[var]][1], 
    x1 = calperiod_plot[[var]][2], 
    y0 = pos_text_cal[[var]][1], 
    lwd = 1.5,
    code = 3, 
    length = 0.1, 
    col = "grey")
  text(
    paste0("Calibration \n period"),
    x = mean(calperiod_plot[[var]]),
    y = pos_text_cal[[var]][2],
    col = "black",
    cex = 1.25)
  text(
    text_cf[[var]][1],
    x = 1960,
    y = pos_text_cal[[var]][3],
    col = col_plot["00_74"],
    cex = 1.25)
  text(
    text_cf[[var]][2],
    x = 1960,
    y = pos_text_cal[[var]][4],
    col = col_plot["75plus"],
    cex = 1.25)
  
  # Add axis
  axis(1, cex.axis = 1.5)
  
  # Plot demographic projections by age group
  lapply(study_param$age_groups, function(i_age) {
    
    # Plot calibrated demographic projections
    points(proj_mortpopu_ldn$year, 
           proj_mortpopu_ldn[[paste0(var, ".", i_age)]]/1e5, 
           col = "black",
           bg = col_plot[i_age],
           pch = 21,
           cex = 2)
    
    # Plot demographic observations
    points(data_obs$year, 
           data_obs[[paste0(var, ".", i_age)]]/1e5, 
           pch = 4, 
           col = col_plot[i_age], 
           cex = 2)
    lines(data_obs$year, 
          data_obs[[paste0(var, ".", i_age)]]/1e5,
          col = col_plot[i_age], 
          lwd = 2)
    
  })
  abline(h = 0)
  
  # Add legend
  legend("right", 
         legend = c(expression("Young ("< 75*" years)"),
                    expression("Old (">= 75*" years)"), 
                    "Projections", 
                    "Observations"),
         pch = c(15, 15, 21, 4),
         pt.bg = c(NA, NA, "darkgrey", NA),
         col = c(col_plot[study_param$age_groups], "black", "darkgrey"),
         cex = 1.25)
  
})

# ---- PANEL C: TEMPORAL CALIBRATION OF DEMOGRAPHIC PROJECTIONS ----

# Select only projection period
ind_plot <- proj_mortpopu_ldn_daily$year >= 2000

# Initialize mortality plot
plot(x = proj_mortpopu_ldn_daily$date[ind_plot],
     y = rep(1, length(proj_mortpopu_ldn_daily$date[ind_plot])),
     ylim = c(0, max(subset(proj_mortpopu_ldn_daily[ind_plot,], 
                            select = paste0("mort.", study_param$age_groups)))),
     type = "n",
     xlab = "Day",
     ylab = "Number of deaths",
     main = "",
     xaxs = "i",
     cex.lab = 1.5,
     cex.axis = 1.5)
abline(h = 0)
title("c) Temporal calibration: mortality", line = 1.0, cex.main = 2)

# Plot daily mortality projections by age groups
run_loop <- lapply(study_param$age_groups, function(i_age) {
  
  # Calibrated daily mortality projections
  lines(proj_mortpopu_ldn_daily$date[ind_plot],
        proj_mortpopu_ldn_daily[[paste0("mort.", i_age)]][ind_plot],
        col = col_plot[i_age])
  
})

# Add legend
legend("topleft", 
       legend = c(
         expression("Young ("< 75*" years)"),
         expression("Old (">= 75*" years)")),
       col = col_plot[study_param$age_groups],
       lty = 1,
       lwd = 1.5,
       cex = 1.4)


# ---- PANELS D-E: ATTRIBUTABLE NUMBERS ----

# Define different ylim for the age groups
ylim <- c("00_74" = 3000, "75plus" = 30000)

# Loop age groups
run_loop <- lapply(study_param$age_groups, function(i_age) {
  
  # Initialize plot
  plot(
    x = unique(data_time$decade), 
    y = rep(1, length(unique(data_time$decade))), 
    type = "n",
    ylim = c(0, ylim[i_age]) / 1e4,
    xlab = "Decade",
    ylab = "AN heat-related deaths (x10,000)",
    main = "",
    xaxt = "n",
    cex = 1.2,
    cex.lab = 1.5,
    cex.axis = 1.5)
  title(title_plot_an[i_age], line = 1.5, cex.main = 2)
  axis(1, 
       at = data_time$decade, 
       labels = paste0(data_time$decade, "s"), 
       cex.axis = 1.5)
  abline(h = 0)
  
  # Loop GCMs
  run_loop <- lapply(study_param$selected_gcms, function(i_gcm) {
    # Loop simulated epidemiological curves
    run_loop <- lapply(1:study_param$n_sim, function(i) {
      # Plot AN
      lines(
        x = unique(data_time$decade),
        y = sapply(
          split(an[[i_age]][[i_gcm]][,paste0("sim", i)], data_time$decade), 
          sum) / 1e4, 
        col = col_plot[i_gcm],
        lty = 2,
        lwd = 0.25)
    })
    # Plot AN for estimated epidemiological curves
    lines(
      x = unique(data_time$decade),
      y = sapply(split(an[[i_age]][[i_gcm]][,"est"], data_time$decade), 
                 sum) / 1e4, 
      col = col_plot[i_gcm], 
      lwd = 3)
  })
  
  # Add legend
  legend("topleft",
         c(study_param$selected_gcms, 
           "Epidemiological estimates", 
           "Epidemiological simulations"),
         col = c(col_plot[study_param$selected_gcms], 1, 1),
         lwd = c(3, 3, 3, 3, 1.25),
         lty = c(1, 1, 1, 1, 2),
         cex = 1.25)
  
})

dev.off()

rm(data_mort_year, col_plot, data_time, title_plot_an, 
   layout_matrix, title_plot_cal, ylab_plot_cal, calperiod_plot, pos_text_cal, 
   text_cf, ind_plot, ylim)

#### FIGURE 4. SUMMARY HEAT-RELATED MORTALITY ##################################

# Calculate annual sums of demographic projections
proj_mortpopu_yearly <- aggregate(
  cbind(mort.00_74, mort.75plus) ~ year, 
  data = proj_mortpopu_ldn_daily, 
  FUN = "sum")

# Calculate annual means of annuals projections
proj_temp_bc_year <- 
  aggregate(.~ lubridate::year(date), data = proj_temp_bc, FUN = mean)
proj_temp_bc_year$date <- NULL
colnames(proj_temp_bc_year)[1] <- "year"

# Extract point_estimate and ci of 21-years periods
values_period <- list(
  gwl = sapply(an_gwl_summary_age, function(x) 
    {x[c("fit", "low", "high")]/21}),
  end_century = sapply(an_endcentury_summary_age, function(x) 
    {x[c("fit", "low", "high")]/21}))

# Define the order of the panel b
order_plot <- c("Young" = "00_74", 
                "Old" = "75plus", 
                "Total" = "total")

# Extract the values for estimation and CI ordered for panel b
values_plot <- sapply(c("fit", "low", "high"), function(value) {
  c(values_period$gwl[value, order_plot], 
    values_period$end_century[value, order_plot])
})
values_plot <- as.matrix(values_plot)

# Build a data frame with dates, corresponding years and decades (for time-based
# grouping)
data_time <- data.frame(
  date = proj_temp_bc$date,
  year = lubridate::year(proj_temp_bc$date),
  decade = lubridate::year(proj_temp_bc$date) %/% 10 * 10)

# Define plotting variables
col_plot <- list(
  age = c("total" = "#000000", "00_74" = "#1B9E77", "75plus" = "#D95F02"),
  period = c("#C969A1", "#CE4441"))

title_plot <- list(
  a = "a) Time trends in age-specific impacts", 
  b = "b) Aggregated 21-year period impacts")

panel_lty <- c("total" = 1,
               "00_74" = 2,
               "75plus" = 2)

legend_plot <- list(
  a = c(
    "Total (Old age)",
    expression("Young ("< 75*" years)"),
    expression("Old (">= 75*" years)")),
  b = c(
    expression(paste("Global warming level (2", degree, "C)")), 
    "End-of-century"))

# ---- PLOT SUMMARY HEALTH IMPACT PROJECTIONS ----
nrow.fig <- 1; ncol.fig <- 2
pdf("outdata/plot/fig4_health_impacts.pdf",
    width = ncol.fig*3*1.5, height = nrow.fig*3*1.5)
layout(matrix(c(1, 2), nrow = nrow.fig, byrow = TRUE))

# ---- PANEL A: DECADAL TIME SERIES OF HEAT-RELATED DEATHS BY AGE GROUP ----
par(mex = 1, mar = c(5, 4, 4, 1), 
    mgp = c(3, 1, 0),
    oma = c(0, 0, 0, 0),
    las = 1)

# Initializate plot
plot(x = an_decade_summary_age$total$time,
     y = an_decade_summary_age$total$fit/10,
     ylim = c(0, max(an_decade_summary_age$total$high, na.rm = TRUE))/10, 
     type = "n",
     xaxt = "n",
     xlab = "Decade",
     ylab = "Yearly heat-related deaths",
     main = "")
title(title_plot$a, line = 0.65, font.main = 2, cex.main = 1)
abline(h = 0)

# Plot x-axis
axis(1, 
     at = data_time$decade, 
     labels = paste0(data_time$decade, "s"), 
     cex.axis = 1)

# Loop age groups
run_loop <- lapply(c("total", study_param$age_groups), function(i_age) {
  
  data <- an_decade_summary_age[[i_age]]
  
  # Plot 95% CI
  polygon(x = c(data$time, rev(data$time)),
          y = c(data$low, rev(data$high))/10,
          border = NA, 
          col = adjustcolor(col_plot$age[i_age], alpha.f = 0.2))
})

# Loop age groups
run_loop <- lapply(c("total", study_param$age_groups), function(i_age) {
  
  data <- an_decade_summary_age[[i_age]]
  
  # Plot point-estimates
  lines(x = data$time,
        y = data$fit/10, 
        lty = panel_lty[i_age], 
        lwd = 2, 
        col = col_plot$age[i_age])
})

# Add legend
legend("topleft", 
       legend_plot$a, 
       lty = panel_lty, 
       lwd = 2,
       col = col_plot$age,
       cex = 1)

# ---- PANEL B: HEAT-RELATED DEATHS BY AGE GROUP IN SPECIFIC PERIODS ----

# Create the plot with the point-estimates
plot(
  x = seq_along(values_plot[,"fit"]),
  y = values_plot[,"fit"],
  xlim = c(min(seq_along(values_plot[,"fit"])) - 0.2, 
           max(seq_along(values_plot[,"fit"])) + 0.2),
  ylim = c(0, max(unlist(values_plot[,"high"])) + 2),
  xaxt = "n",
  pch = 21,
  cex = 2,
  bg = c(rep(col_plot$period[1], length(values_plot[,"fit"])/2), 
         rep(col_plot$period[2], length(values_plot[,"fit"])/2)),
  xlab = "",
  ylab = "Yearly heat-related deaths",
  main = "")
title(title_plot$b, 
      line = 0.65, font.main = 2, cex.main = 1)
abline(h = 0)
abline(v = mean(seq_along(values_plot[,"fit"])))

# Add x-axis labels
axis(1, 
     at = seq_along(values_plot[,"fit"]), 
     labels = rep(names(order_plot), 2), las = 1)

# Add 95% CI
arrows(x0 = seq_along(values_plot[,"fit"]), 
       y0 = unlist(values_plot[,"low"]),
       x1 = seq_along(values_plot[,"fit"]),
       y1 = unlist(values_plot[,"high"]),
       angle = 90, code = 3, length = 0.1, col = "black")

# Add legend
legend("topleft", 
       legend = legend_plot$b, 
       pch = 21,
       bg = "white",
       pt.bg = col_plot$period, 
       cex = 1)

# Add title for the plot
mtext("Heat-related mortality projections (London, SSP2-4.5)", 
      side = 3, outer = TRUE, line = -2.2, cex = 1.3, font = 2)

dev.off()

rm(proj_mortpopu_yearly, proj_temp_bc_year, values_period, order_plot, 
   values_plot, data_time, col_plot, title_plot, legend_plot)

#### FIGURE S1. UNCERTAINTY TEMPERATURE-MORTALTITY ASSOCIATION ###############

# Define the temperatures for the x-axis
pred_perc <- c(seq(0, 1, 0.1), 2:98, seq(99, 100, 0.1))
x_temp <- quantile(data_tempmort$tmean, pred_perc / 100)
rm(pred_perc)

# Consider also temperatures beyond the boundaries of observed temperatures
x_temp_extra <- seq(max(x_temp), max(x_temp) + 5, length = 11) 
x_temp <- c(x_temp, x_temp_extra[-1])
rm(x_temp_extra)

# Define plotting variables
ymax <- 5
col_plot <- c("00_74" = "#1B9E77", "75plus" = "#D95F02")
title_plot <- c(
  "00_74" = expression("a) Young ("< 75*" years)"),
  "75plus" = expression("b) Old (">= 75*" years)"))

# PLOT FIGURE WITH AGE-SPECIFIC TEMPERATURE-MORTALITY ASSOCIATIONS
nrow.fig <- 1; ncol.fig <- 2
pdf("outdata/plot/figs1_simulations_age_specific_associations.pdf",
    width = ncol.fig*3*1.5, height = nrow.fig*3*1.5)
layout(matrix(seq(ncol.fig * nrow.fig), nrow = nrow.fig, byrow = TRUE))
par(mex = 0.8, mgp = c(3, 1, 0), las = 1, oma = c(0, 0, 0, 0))

run_loop <- lapply(study_param$age_groups, function(i_age) {
  
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
  coefsim <- subset(coef_age[[i_age]], select = -est)
  rrsim <- exp(bcen %*% coefsim)
  
  # Extract color for the correspondin age group
  col_age <- col_plot[i_age]
  
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

rm(col_plot, title_plot, ncol.fig, nrow.fig, x_temp, ymax)

#### FIGURE S2. PROCESS GRIDDED TEMPERATURES ###################################

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
pdf(file = "outdata/plot/figs2_process_gridded_temperatures.pdf",
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
    x1 <- -0.6
    y1 <- 51.1
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

#### FIGURE S3. ANNUAL OBSERVED AGE-SPECIFIC SEASONALITY PATTERNS #### 

title_plot <- c(
  "00_74" = expression("a) Young ("< 75*" years)"),
  "75plus" = expression("b) Old (">= 75*" years)"))

nrow.fig <- 1; ncol.fig <- 2
pdf("outdata/plot/figs3_mean_annual_cycle_mortality.pdf",
    width = ncol.fig*3*1.5, height = nrow.fig*3*1.5)
layout(matrix(seq(ncol.fig * nrow.fig), nrow = nrow.fig, byrow = TRUE))
par(mex = 0.8, mgp = c(3, 1, 0), las = 1, oma = c(0, 0, 0, 0))

run_loop <- lapply(study_param$age_groups, function(i_age) {
  
  # Calculate mean of daily deaths by day of the year
  deathdoy <- tapply(
    data_tempmort[[paste0("mort.", i_age)]], 
    as.numeric(format(data_tempmort$date, "%j")), 
    mean, na.rm = TRUE)[seq(365)]
  
  # Seasonality weights by day of the year
  weights_seas <- deathdoy / sum(deathdoy)
  
  ts_weight_mort <- data.frame(
    mort = data_tempmort[[paste0("mort.", i_age)]],
    doy = as.numeric(format(data_tempmort$date, "%j")),
    year = lubridate::year(data_tempmort$date))
  # Remove 2012 for the plot because the year is not complete
  ts_weight_mort <- subset(ts_weight_mort, year != 2012)
  
  ts_weight_mort_sum <- aggregate(mort ~ year, data = ts_weight_mort, FUN = "sum")
  colnames(ts_weight_mort_sum)[colnames(ts_weight_mort_sum) == "mort"] <- "sum_mort"
  
  ts_weight_mort <- merge(ts_weight_mort, ts_weight_mort_sum, by = "year")
  
  ts_weight_mort$weight <- ts_weight_mort$mort / ts_weight_mort$sum_mort
  
  plot(ts_weight_mort$doy, 
       ts_weight_mort$weight*100, 
       pch = 16, 
       col = "grey", 
       cex = 0.5,
       xlab = "Day of the year",
       ylab = "Yearly percentage of deaths (%)")
  title(title_plot[i_age], line = 0.65, font.main = 1, cex.main = 1)
  lines(1:365, weights_seas*100, col = "blue", lwd = 2)
  
})
# Title common to both panels
mtext("Observed mortality (London, 1990-2012)", 
      side = 3, outer = TRUE, line = -2.2, cex = 1.3, font = 2)
dev.off()

rm(title_plot)

#### FIGURE S4. GLOBAL WARMING LEVEL PERIODS ####################################

# Define plotting parameters
col_plot <- c("#118DFF", "#750985", "#C83D95")
levels_gwl <- c("1.5", "2", "3")
offsets_gwl <- c(-0.15, 0, 0.15)
points_gwl <- c(15, 16, 17)

# Plot different GWL depending on the GCM
pdf(file = "outdata/plot/figs4_global_warming_level_periods.pdf",
    width = 7, height = 5)
par(mfrow = c(1, 1), mar = c(5, 8, 3.5, 0.5), las = 1)

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
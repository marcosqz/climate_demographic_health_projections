#### TODO: ADD DESCRIPTION OF THIS SCRIPT

# Extra code! It is not needed to run the full tutorial.
# In this script we estimate age-specific temperature-mortality association for 
# London using a single interaction model

#### LOAD LIBRARIES ############################################################

library(INLA) # inla
library(dlnm) # logknots
library(splines) # ns

#### LOAD DATA #################################################################

load("outdata/file/data_tempmort.RData")

#### SET VARIABLES DEFINING THE DLNM MODEL #####################################

# PARAMETERS DEFINING DLNMs
dlnm_var <- list(
  var_prc = c(10, 75, 90),
  var_fun = "ns",
  lag_fun = "ns",
  max_lag = 21,
  lagnk = 3)

# DEFINE THE EXPOSURE- AND LAG-RESPNSE FUNCTIONS
#### TODO: Why Masselot uses "bs", it is not necessary the "ns" to the log-liner extrapolation
argvar <- list(fun = dlnm_var$var_fun,
               knots = quantile(data_tempmort$tmean, 
                                dlnm_var$var_prc/100, na.rm = TRUE),
               Bound = range(data_tempmort$tmean, na.rm = TRUE))
arglag <- list(fun = dlnm_var$lag_fun, 
               knots = logknots(dlnm_var$max_lag, nk = dlnm_var$lagnk))

# BUILD THE CROSS-BASIS
cb <- crossbasis(data_tempmort$tmean, lag = dlnm_var$max_lag, argvar, arglag)

# BUILD THE BASIS WITH THE SEASONALITY TERM
seas <- ns(data_tempmort$date, 
           df = round(8 * length(data_tempmort$date) / 365.25))

# DATASET FROM WIDE TO LONG TO INCLUDE IT IN A SINGLE MODEL
data_long <- reshape(data_tempmort,
                     varying = c("mort.00_74", "mort.75plus"),
                     v.names = "mortality",
                     timevar = "age_group",
                     times = c("00_74", "75plus"),
                     direction = "long")

# BIND TO COPIES OF cb AND seas TO FOLLOW THE SAME STRUCTURE OF data_long 
cb <- rbind(cb, cb)
seas <- rbind(seas, seas)

# REMOVE NAs (inla FUNCTION DOES NOT REMOVE THEM)
ind_na <- apply(cb, 1, function(x) any(is.na(x)))

data_long <- data_long[!ind_na,]
cb <- cb[!ind_na,]
seas <- seas[!ind_na,]

# ARREGEMENTS FOR THE MODEL
colnames(cb) <- paste0("cb", 1:ncol(cb))
data_long$age_group <- factor(data_long$age_group)

# RUN INTERACTION MODEL
model <- 
  inla(mortality ~ cb*age_group + factor(dow)*age_group + seas*age_group,
       data = data_long, family = "poisson", 
       control.compute = list(config = TRUE))

# CREATE A LIST OF PART OF SAMPLE WE ARE INTERESTED (THE CB AND INTERACTION)
selection_names <- c(colnames(cb), paste0(colnames(cb), ":age_group75plus"))
selection_list <- as.list(rep(1, length(selection_names)))
names(selection_list) <- selection_names

# EXTRACT THE ENSAMBLE OF SAMPLE COEFFICIENTS FROM THE JOINT POSTERIOR 
# DISTRIBUTIONS
coef <- inla.posterior.sample(100, model, selection = selection_list)

# INITIALIZE OBJECT TO SAVE THE OUTPUTS
coefsim_age <- list()

# ARRANGE THE COEFFICIENTS OF THE CB IN THE SAME WAY AS IN CODE 02
coefsim_age[["00_74"]] <- sapply(coef, function(x) {
  sapply(colnames(cb), function(coef_name) {
    x$latent[paste0(coef_name, ":1"),]
  })
})

coefsim_age[["75plus"]] <- sapply(coef, function(x) {
  sapply(colnames(cb), function(coef_name) {
    x$latent[paste0(coef_name, ":1"),] + 
      x$latent[paste0(coef_name, ":age_group75plus:1"),]
  })
})

#### SAVE OUTPUTS ##############################################################
save(coefsim_age, 
     file = "outdata/file/extra_epi_model_coefsimage_interaction.RData")

#### PLOT AGE-SPECIFIC ASSOCIATIONS ############################################

# DEFINE VARIABLES FOR AGE GROUPS AND PLOTS
# Age groups
age_parameters <- data.frame(
  response = c("mort.00_74", "mort.75plus"),
  groups = c("00_74", "75plus"),
  col_est = c("#1B9E77", "#D95F02"))
n_groups <- 2

# Title plot
title_plot <- list(
  expression("Young ("< 75*" years)"),
  expression("Old (">= 75*" years)"))

# Plot
# TODO: I HAVE TO BETTER EXPLAIN ALL THE SMALL VARIABLES THAT ARE USED FOR THE
# PLOT (e.g. perc) AND CHANGE SOME NAME (e.g. TPLOT, RRPLOT) THAT ARE NOT 
# INFORMATIVE ENOUGH
main <- NULL
title <- NULL
perc <- c(1, 5, 95, 99)/100
col <- c("blue", "cyan", "orange", "red", "red4", "purple4")
ymax <- 5

# INITIALIZE PLOT
nrow.fig <- 1; ncol.fig <- 2
pdf("outdata/plot/fig99_extra_age_specific_associations_interaction.pdf",
    width = ncol.fig*3*1.5, height = nrow.fig*3*1.5)
layout(matrix(seq(ncol.fig * nrow.fig), nrow = nrow.fig, byrow = TRUE))
par(mex = 0.8, mgp = c(3, 1, 0), las = 1, oma = c(0, 0, 0, 0))

for(iter in 1:n_groups) {
  
  # TEMPERATURES TO PLOT IN THE X-AXIS
  predper <- c(seq(0, 1, 0.1), 2:98, seq(99, 100, 0.1))
  tper <- quantile(data_tempmort$tmean, predper / 100)
  tper_extra <- seq(max(tper), max(tper) + 5, length = 11)
  tper <- c(tper, tper_extra[-1])
  
  # GENERATE THE CROSSBASIS FOR TEMPERATURES
  # (Note, we have here the non-reduced coefficients. With the reduced we use
  # the onebasis for only the exposure. Here we use an equivalent alternative,
  # where the crossbasis is calculated with a matrix that for each exposure has
  # a series of lags that are all the same)
  cb <- crossbasis(
    matrix(rep(tper, each = dlnm_var$max_lag + 1), 
           ncol = dlnm_var$max_lag + 1, 
           byrow = TRUE), 
    lag = dlnm_var$max_lag, argvar, arglag)
  
  # CALCULATE THE SAMPLE OF RELATIVE RISKS
  rrsim <- cb %*% coefsim_age[[iter]]
  
  # LOCATE POSITION OF THE TEMPERATURE WITH MINIMUM RISK
  immt <- which.min(rowMeans(rrsim))
  coef_mmt <- rrsim[immt,] 
  
  # BUILD A MATRIX REPLICATING THE COEFFICIENTS IN THE MMT
  centring_matrix <- 
    matrix(rep(coef_mmt, length(tper)), 
           nrow = length(tper), byrow =TRUE)
  
  # CENTER RISK IN THE MMT
  rrsim <- exp(rrsim - centring_matrix)
  
  # GET MIN, MEDIAN AND MAX OF RRs TO PLOT THEM
  rrmin <- apply(rrsim, 1, min)
  rrmax <- apply(rrsim, 1, max)
  rrmedian <- apply(rrsim, 1, median)
  
  # FIND SPECIFIC TEMPERATURES AND RISKS TO PLOT IN THE LEGEND
  ind <- (tper %in% quantile(data_tempmort$tmean, perc)) |
    (tper %in% c(max(data_tempmort$tmean), max(data_tempmort$tmean) + 2))
  
  Tplot <- tper[ind]
  RRplot <- rrmedian[ind]
  RRplot_min <- rrmin[ind]
  RRplot_max <- rrmax[ind]
  
  # EMPTY PLOT
  plot(tper, rrsim[,1], 
       xlab = expression(paste("Temperature (", degree, "C)")),
       ylab = "Relative risk", ylim = c(1.000, ymax), 
       log = "y", lwd = 2,
       col = age_parameters$col_est[iter], 
       ci.arg = list(col = rgb(0.8, 0.8, 0.8)), 
       type = "n")
  title(title_plot[[iter]], line = 0.65, font.main = 1, cex.main = 1)
  
  # PLOT SAMPLE OF RR
  for(i in 1:100){
    lines(tper, rrsim[,i], col = age_parameters$col_est[iter], lwd = 0.1)
  }
  
  # PLOT MEDIAN RR
  lines(tper, rrmedian, col = "black", lwd = 2)
  
  # OTHER PLOTTING ASPECTS
  abline(h = 1)
  abline(v = max(data_tempmort$tmean), lty = 2, col = "grey")
  
  # PLOT LINES AND RR FOR THE SPECIFIC TEMPERATURES OF INTEREST
  for(i in 1:length(col)){
    lines(c(Tplot[i], Tplot[i]), c(1, RRplot[i]), 
          lty = 4, col = col[i], lwd = 1.5)
  }
  rm(i)
  
  points(Tplot, RRplot, col = "black", pch = 21, 
         cex = 1.5, bg = col, lwd = 2)
  
  legend("top", paste0("P", c(sprintf("%02d", perc*100), c(100, 100)), 
                       c("", "", "", "", "", "+1ºC"), ": ",
                       formatC(RRplot, digits = 2, format = "f", flag="0"), " (",
                       formatC(RRplot_min, digits = 2, format = "f", flag="0"), ",",
                       formatC(RRplot_max, digits = 2, format = "f", flag="0"), ")"), 
         col = col, pch = 19, box.lty = 0, cex = 1.000)
  
}

mtext("Temperature-mortality associations (London, 1990-2012)", 
      side = 3, outer = TRUE, line = -2.2, cex = 1.3, font = 2)

dev.off()
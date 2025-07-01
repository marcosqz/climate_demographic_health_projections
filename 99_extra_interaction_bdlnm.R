#### TODO: ADD DESCRIPTION OF THIS SCRIPT

# Extra code! It is not needed to run the full tutorial.
# In this script we estimate age-specific temperature-mortality association for 
# London using a single interaction model

#### LOAD LIBRARIES ############################################################

library(INLA) # inla
library(dlnm) # logknots
library(splines) # ns

#### LOAD DATA #################################################################

load("indata/processed/data_obs_temp_mort.RData")

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

# SAVE DIMENSIONS FOR REDUCING COEFFICIENTS
vx <- attr(cb, "df")[1]
vl <- attr(cb, "df")[2]

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

coefsim_age_interaction <- sapply(coef, function(x) {
  sapply(colnames(cb), function(coef_name) {
      x$latent[paste0(coef_name, ":age_group75plus:1"),]
  })
})

# REDUCE COEFFICIENTS
# (The object with the reduced coefficient has the same structure as coefsim_age
# from the frequentist stratified approach, allowing compatibility with the rest
# of the code when using the B-DLNM interaciton model)

# See https://doi.org/10.1186/1471-2288-13-1 for details on the reduction of the
# coefficients, in which the formula for the dimension-reducing matrix M for the
# cumulative overall association is M = (1^T_{L+1} * C) ⊗ I_{vx}

# Uni-dimensional basis for lags
C <- onebasis(x = 0:dlnm_var$max_lag,
              fun = dlnm_var$lag_fun,
              knots = logknots(dlnm_var$max_lag, nk = dlnm_var$lagnk),
              intercept = TRUE)

# Row vector of ones of length (L+1)
ones <- matrix(1, nrow = 1, ncol = dlnm_var$max_lag + 1)

# Matrix multiplication
product <- ones %*% C

# Identity matrix of size vx
I_vx <- diag(vx)

# Compute dimension-reducing matrix
M <- I_vx %x% product

# (Remove temporary datasets)
rm(C, ones, product, I_vx, vx, vl)

# Apply the reduction
coefsim_age[["00_74"]] <- 
  apply(coefsim_age[["00_74"]], 2, function(x) M %*% x)
coefsim_age[["75plus"]] <- 
  apply(coefsim_age[["75plus"]], 2, function(x) M %*% x)
coefsim_age_interaction <- 
  apply(coefsim_age_interaction, 2, function(x) M %*% x)

#### SAVE OUTPUTS ##############################################################
save(coefsim_age, 
     file = "outdata/file/01_epi_model/coefsim_age_interaction.RData")

#### PLOT AGE-SPECIFIC ASSOCIATIONS ############################################

# ADD THE INTERACTION TO coefsim_age OBJECT to plot it in the loop
coefsim_age$interaction <- coefsim_age_interaction
  
# DEFINE VARIABLES FOR AGE GROUPS AND PLOTS
# Age groups
age_parameters <- data.frame(
  response = c("mort.00_74", "mort.75plus", "NA"),
  groups = c("00_74", "75plus", "interaction"),
  col_est = c("#1B9E77", "#D95F02", "#7570B3"))
n_groups <- 3

# Title plot
title_plot <- list(
  expression("a) Young ("< 75*" years)"),
  expression("b) Old (">= 75*" years)"),
  expression("c) Interaction old vs young"))

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
nrow.fig <- 1; ncol.fig <- 3
pdf("outdata/plot/figs4_extra_age_specific_associations_interaction.pdf",
    width = ncol.fig*3*1.5*1, height = nrow.fig*3*1.5)
layout(matrix(seq(ncol.fig * nrow.fig), nrow = nrow.fig, byrow = TRUE))
par(mex = 1, mgp = c(3, 1, 0), las = 1, oma = c(0, 0, 0, 0))

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
  bcen <- onebasis(
    x = tper, 
    fun = argvar$fun, 
    knots = argvar$knots, 
    Boundary.knots = argvar$Bound)
  
  # CALCULATE THE SAMPLE OF RELATIVE RISKS
  rrsim <- bcen %*% coefsim_age[[iter]]
  
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
  plot(tper, rrsim[,1], type = "n",
       xlab = expression(paste("Temperature (", degree, "C)")),
       ylab = "Relative risk", ylim = c(1.000, ymax), 
       log = "y", cex.lab = 1.5, cex.axis = 1.5)
  title(title_plot[[iter]], line = 0.65, font.main = 1, cex.main = 1.5)
  
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
          lty = 4, col = col[i], lwd = 2)
  }
  rm(i)
  
  points(Tplot, RRplot, col = "black", pch = 21, 
         cex = 1.5*1.5, bg = col, lwd = 1.5*1.5)
  
  legend("top", paste0("P", c(sprintf("%02d", perc*100), c(100, 100)), 
                       c("", "", "", "", "", "+2ºC"), ": ",
                       formatC(RRplot, digits = 2, format = "f", flag="0"), " (",
                       formatC(RRplot_min, digits = 2, format = "f", flag="0"), ",",
                       formatC(RRplot_max, digits = 2, format = "f", flag="0"), ")"), 
         col = col, pch = 19, box.lty = 0, cex = 1.5)
  
}

mtext("Age-specific temperature-mortality associations with an interaction model and B-DLNM (London, 1990-2012)", 
      side = 3, outer = TRUE, line = -2.2, cex = 1.3, font = 2)

dev.off()
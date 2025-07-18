################################################################################

# This script compares the heat-related mortality aggregation across age groups,
# using identical versus independent seeds in the  Monte Carlo simulations of 
# the exposure-response functions for each group.

#### LOAD LIBRARIES ############################################################

library(dlnm) # logknots
library(splines) # ns
library(MASS) # mvrnorm

#### LOAD DATA #################################################################

load("indata/processed/data_obs_temp_mort.RData")
data_tempmort$tmean <- round(data_tempmort$tmean, 4) # TODO: See how the association changes in panel a) if I remove this row

#### SET VARIABLES DEFINING THE DLNM MODEL #####################################

# PARAMETERS DEFINING DLNMs
dlnm_var <- list(
  var_prc = c(10, 75, 90),
  var_fun = "ns",
  lag_fun = "ns",
  max_lag = 21,
  lagnk = 3)

# USE 1000 SIMULATION TO SEE MORE CLEARLY THE CORRELATIONS IN THE PLOTS
n_sim <- 1000

# DEFINE THE EXPOSURE- AND LAG-RESPNSE FUNCTIONS
argvar <- list(fun = dlnm_var$var_fun,
               knots = quantile(data_tempmort$tmean, 
                                dlnm_var$var_prc/100, na.rm = TRUE),
               Bound = range(data_tempmort$tmean, na.rm = TRUE))
arglag <- list(fun = dlnm_var$lag_fun, 
               knots = logknots(dlnm_var$max_lag, nk = dlnm_var$lagnk))

# BUILD THE CROSS-BASIS
cb <- crossbasis(data_tempmort$tmean, lag = dlnm_var$max_lag, argvar, arglag)

# DEFINE VARIABLE NAMES FOR AGE GROUPS
age_parameters <- data.frame(
  response = c("mort.00_74", "mort.75plus"),
  groups = c("00_74", "75plus"))
n_groups <- nrow(age_parameters)

# LOOP HEAT-RELATED MORTALITY AGGREGATION USING IDENTICAL VERSUS INDEPENDENT IN 
# THE MONTE CARLO SIMULATIONS
impacts_age <- lapply(c("seed", "no_seed"), function(iter) {
  
  # INITIALIZE OBJECTS TO SAVE THE OUTPUTS
  coefsim_age <- list()
  mmt_age <- c()
  
  # LOOP SIMULATIONS OF EXPOSURE-RESPONSE FUNCTION FOR EACH AGE GROUP
  for(i in 1:n_groups) {
    
    model <- glm(
      as.formula(paste0(
        age_parameters$response[i], 
        "~ cb + factor(dow) + ns(date, df = round(8 * length(date) / 365.25))")),
      data = data_tempmort, family = quasipoisson)
    
    # REDUCE COEFFICIENTS TO KEEP THE CUMULATIVE EXPOSURE-RESPONSE FUNCTION
    reduced <- crossreduce(cb, model)
    
    # FIND THE MINIMUM MORTALITY TEMPERATURE
    mmt_age[age_parameters$groups[i]] <- 
      reduced$predvar[which.min(reduced$RRfit)]
    
    # EXTRACT MODEL COEFFICIENTS AND VARIANCES
    coef <- coef(reduced)
    vcov <- vcov(reduced)
    
    # GENERATE MONTE CARLO SAMPLES TO GENERATE SIMUALATIONS OF THE
    # RESPONSE FUNCTION
    if(iter == "seed") {
      set.seed(13041975) # Important! Same random seed for each group
    }
    coefsim_age[[age_parameters$groups[i]]] <- t(mvrnorm(n_sim, coef, vcov))
    
  }
  
  # LOOP CALCULATION HEALTH IMPACTS FOR EACH AGE GROUP
  impacts_age <- lapply(age_parameters$groups, function(i_age){
    
    # EXPOSURE-RESPONSE BASIS AT THE AGE-SPECIFIC MMT
    cenvec <- onebasis(
      x = mmt_age[i_age], 
      fun = argvar$fun, 
      knots = argvar$knots,
      Boundary.knots = argvar$Bound)

    # EXPOSURE-RESPONSE BASIS AT THE PROJECTED TEMPERATURES CENTERED AT THE
    # AGE-SPECIFIC MMT
    bcen <- scale(
      onebasis(
        x = data_tempmort$tmean,
        fun = argvar$fun,
        knots = argvar$knots,
        Boundary.knots = argvar$Bound),
      center = cenvec, scale = FALSE)
    
    # RELATIVE RISKS AND ATTRIBUTABLE FRACTIONS CALCULATION
    rrsim <- exp(bcen %*% coefsim_age[[i_age]])
    afsim <- (rrsim - 1) / rrsim # AF = (RR-1)/RR
  
    # CALCULATE HEAT-RELARED DEATHS BY SETTING AF=0 TO TEMPERATURES BELOW MMT
    ind_heat <- data_tempmort$tmean > mmt_age[i_age]
    afsim[!ind_heat,] <- 0
    
    # STRUCUTRE IMPACTS AS A DATA.FRAME
    colnames(afsim) <- paste0("sim", 1:n_sim)
    impacts <- data.frame(afsim); rm(afsim)
    
    # MATRIX OF LAGGED MORTALITES (FORWARD PERSPECTIVE)
    lagged_mort <- tsModel::Lag(
      v = data_tempmort[[paste0("mort.", i_age)]],
      k = -seq(0,dlnm_var$max_lag)) # an integer vector giving lag numbers
    
    # CALCULATE THE DAILY MORTALITIES MEAN OVER THE LAG PERIOD
    lagged_mort <- rowMeans(lagged_mort)
    
    # CREATE A MATRIX BY REPLICATING THE LAGGED MORTALITY VECTOR ACROSS COLUMNS,
    # RESULTING IN A MATRIX WITH 'NSIM' COLUMNS
    matrix_mort <- matrix(rep(lagged_mort, n_sim), ncol = n_sim)
    
    # CALCULATE ATTRIBUTABLE NUMBER 
    an <- colSums(impacts * matrix_mort, na.rm = TRUE) # AN = AF * MORT
    
    return(an)
    
  })
})
names(impacts_age) <- c("seed", "no_seed") 

# ALL-AGE ANs BY SUMMING THE PAIRED AGE-SPECIFIC ANs
estimates_aggr_an <- lapply(impacts_age, function(x) {
  
  # CALCULATE MEAN AND CONFIDENCE INTERVAL OF THE SUM OF ANs
  estimates_aggr_an <- quantile(x[[1]] + x[[2]], c(0.025, 0.5, 0.975))
  
  # FORMAT IT WITH ONE DECIMAL NUMBER
  estimates_aggr_an <-format(round(estimates_aggr_an, 1), nsmall = 1)
  
  return(estimates_aggr_an)
}); names(estimates_aggr_an) <- c("seed", "no_seed") 

# TITLE PLOTS
title_plot <- c(
  "seed" = "a) Identical seeds",
  "no_seed" = "b) Independent seeds")

# PLOT CORRELATION BETWEEN PARALLEL MONTE CARLO SIMULATIONS OF ANs
nrow.fig <- 1; ncol.fig <- 2
pdf(file = "outdata/plot/figs3_correlation_an_monte_carlo.pdf",
    width = ncol.fig * 4.5, height = nrow.fig * 4.5)
par(mfrow = c(nrow.fig, ncol.fig),
    mar = c(4, 4, 2, 2),    
    oma = c(0, 0, 0, 0),
    mgp = c(2.5, 0.5, 0),
    las = 1)

par(mfrow = c(1, 2))

lapply(1:length(impacts_age), function(i) {
  plot(impacts_age[[i]][[1]], impacts_age[[i]][[2]],
       xlab = expression("Heat-related mortality (young)"), 
       ylab = expression("Heat-related mortality (old)"),
       xlim = range(impacts_age[[1]][[1]], impacts_age[[2]][[1]]),
       ylim = range(impacts_age[[1]][[2]], impacts_age[[2]][[2]]),
       main = title_plot[i])
  text(x = mean(c(par("usr")[1], par("usr")[2])), 
       y = par("usr")[4],
       labels = paste0(
         "Heat-related mortality (total): ", 
         estimates_aggr_an[[i]][2], " (", 
         estimates_aggr_an[[i]][1], ", ", 
         estimates_aggr_an[[i]][3], ")"),
       pos = 1, 
       offset = 0.15, 
       cex = 0.75, 
       col = "red")
})

dev.off()
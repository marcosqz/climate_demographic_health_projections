#### TODO: ADD DESCRIPTION OF THIS SCRIPT

# Estimation of age-specific temperature-mortality association for London using
# stratified models

#### LOAD LIBRARIES ############################################################

library(dlnm) # logknots
library(splines) # ns
library(MASS) # mvrnorm

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

# DEFINE VARIABLE NAMES FOR AGE GROUPS
age_parameters <- data.frame(
  response = c("mort.00_74", "mort.75plus"),
  groups = c("00_74", "75plus"))
n_groups <- 2

# INITIALIZE OBJECTS TO SAVE THE OUTPUTS
coefsim_age <- list()
mmt_age <- c()

# LOOP AGE-SPECIFIC TEMPERATURE-MORTALITY MODELS
for(i in 1:n_groups) {
  
  model <- glm(
    as.formula(paste0(
      age_parameters$response[i], 
      "~ cb + factor(dow) + ns(date, df = round(8 * length(date) / 365.25))")),
    data = data_tempmort, family = quasipoisson)
  
  # REDUCE COEFFICIENTS TO KEEP THE CUMULATIVE EXPOSURE-RESPONSE FUNCTION
  reduced <- crossreduce(cb, model)
  
  # FIND THE MINIMUM MORTALITY TEMPERATURE
  mmt_age[age_parameters$groups[i]] <- reduced$predvar[which.min(reduced$RRfit)]
  
  # CENTER CURVES IN THE MMT
  reduced <- crossreduce(cb, model, cen = mmt_age[age_parameters$groups[i]])
  
  # EXTRACT MODEL COEFFICIENTS AND VARIANCES
  coef <- coef(reduced)
  vcov <- vcov(reduced)
  
  # GENERATE MONTE CARLO SAMPLES TO GENERATE SIMUALATIONS OF THE
  # RESPONSE FUNCTION
  set.seed(13041975) # Important! Same random seed for each group
  coefsim_age[[age_parameters$groups[i]]] <- t(mvrnorm(100, coef, vcov))
  
}

#### SAVE OUTPUTS ##############################################################

save(argvar, file = "outdata/file/epi_model_argvar.RData")
save(arglag, file = "outdata/file/epi_model_arglag.RData")
save(dlnm_var, file = "outdata/file/epi_model_dlnm_varibles.RData")
save(coefsim_age, file = "outdata/file/epi_model_coefsimage.RData")
save(mmt_age, file = "outdata/file/epi_model_mmt.RData")
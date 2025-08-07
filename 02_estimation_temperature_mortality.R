################################################################################

# This script prepares the data to estimate the association between temperature 
# and mortality for two age groups (<75 and ≥75 years) in London during 
# 1990–2012 using stratified quasi-Poisson regression models combined with 
# distrusted lag non-linear models (DLMNs). 

#### LOAD LIBRARIES ############################################################

library(dlnm) # logknots
library(splines) # ns
library(MASS) # mvrnorm

#### LOAD DATA #################################################################

load("indata/processed/data_obs_temp_mort.RData")
load("indata/processed/study_parameters.RData")

#### SET VARIABLES DEFINING THE DLNM MODEL #####################################

# PARAMETERS DEFINING DLNMs
dlnm_var <- list(
  var_prc = c(10, 75, 90),
  var_fun = "ns",
  lag_fun = "ns",
  max_lag = 21,
  lagnk = 3)

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
  response = paste0("mort.", study_param$age_groups),
  groups = study_param$age_groups)
n_groups <- length(study_param$age_groups)

# INITIALIZE OBJECTS TO SAVE THE OUTPUTS
coef_age <- list()
vcov_age <- list()
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
  reduced <- crossreduce(cb, model, at = data_tempmort$tmean)
  
  # FIND THE MINIMUM MORTALITY TEMPERATURE
  mmt_age[age_parameters$groups[i]] <- reduced$predvar[which.min(reduced$RRfit)]
  
  # EXTRACT MODEL COEFFICIENTS AND VARIANCES
  coef_age[[age_parameters$groups[i]]] <- coef(reduced)
  vcov_age[[age_parameters$groups[i]]] <- vcov(reduced)
  
  # GENERATE MONTE CARLO SAMPLES TO GENERATE SIMUALATIONS OF THE
  # RESPONSE FUNCTION
  set.seed(13041975) # Important! Same random seed for each group
  coefsim_age[[age_parameters$groups[i]]] <- 
    t(mvrnorm(study_param$n_sim,
              coef_age[[age_parameters$groups[i]]],
              vcov_age[[age_parameters$groups[i]]]))
  
}; rm(i)

#### SAVE OUTPUTS ##############################################################

save(argvar, file = "outdata/file/01_epi_model/argvar.RData")
save(arglag, file = "outdata/file/01_epi_model/arglag.RData")
save(dlnm_var, file = "outdata/file/01_epi_model/dlnm_var.RData")
save(coef_age, file = "outdata/file/01_epi_model/coef_age.RData")
save(vcov_age, file = "outdata/file/01_epi_model/vcov_age.RData")
save(coefsim_age, file = "outdata/file/01_epi_model/coefsim_age.RData")
save(mmt_age, file = "outdata/file/01_epi_model/mmt_age.RData")
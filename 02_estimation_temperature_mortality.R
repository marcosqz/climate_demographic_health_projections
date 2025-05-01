#### LOAD LIBRARIES ####

library(dlnm) # logknots
library(splines) # ns
library(MASS) # mvrnorm

#### LOAD DATA

load("outdata/file/data_tempmort.RData")

#### 2. ESTIMATION EXPOSURE-RESPONSE FUNCTIONS ####  
#### (ACCOUNTING FOR DIFFERENT POPULATION VULNERABILITIES) ####

# SET VARIABLE DEFINING THE DLNM MODEL

dlnm_var <- list(
  var_prc = c(10, 75, 90),
  var_fun = "ns",
  lag_fun = "ns",
  max_lag = 21,
  lagnk = 3)

# DEFINE THE CROSS-BASIS TERM
argvar <- list(fun = dlnm_var$var_fun, # TODO: Why masselot uses "bs", it is not necessary the "ns" to the log-liner extrapolation
               knots = quantile(data_tempmort$tmean, 
                                dlnm_var$var_prc/100, na.rm = TRUE),
               Bound = range(data_tempmort$tmean, na.rm = TRUE))
arglag <- list(fun = dlnm_var$lag_fun, 
               knots = logknots(dlnm_var$max_lag, nk = dlnm_var$lagnk))

# BUILD THE CROSS-BASIS
cb <- crossbasis(data_tempmort$tmean, lag = dlnm_var$max_lag, argvar, arglag)

# VISUALLITATION PARAMETERS BY AGE GROUPS
age_parameters <- data.frame(
  response = c("mort.00_74", "mort.75plus"),
  groups = c("00_74", "75plus"))

coefsim_age <- list()
mmt_age <- c()

# RUN AGE-SPECIFIC MODELS
for(i in 1:2) { # AGE GROUPS
  
  model <- glm(
    as.formula(paste0(
      age_parameters$response[i], 
      "~ cb + factor(dow) + ns(date, df = round(8 * length(date) / 365.25))")),
    data = data_tempmort, family = quasipoisson)
  
  # Predict the curve
  reduced <- crossreduce(cb, model)
  
  # Find the mmt
  mmt <- reduced$predvar[which.min(reduced$RRfit)]
  mmt_age[age_parameters$groups[i]] <- mmt
  
  # Center curves in the mmt
  reduced <- crossreduce(cb, model, cen = mmt)
  
  # Extract model coefficients and variances
  coef <- coef(reduced)
  vcov <- vcov(reduced)
  
  # Create the 1000 samples of coefficients
  set.seed(13041975) # TODO: Check this! I think the set.seed it should be the same but I'm not sure
  coefsim <- t(mvrnorm(100, coef, vcov))
  coefsim_age[[age_parameters$groups[i]]] <- coefsim
  
}

# SAVE DATASETS
save(argvar, file = "outdata/file/epi_model_argvar.RData")
save(arglag, file = "outdata/file/epi_model_arglag.RData")
save(dlnm_var, file = "outdata/file/epi_model_dlnm_varibles.RData")
save(coefsim_age, file = "outdata/file/epi_model_coefsimage.RData")
save(mmt_age, file = "outdata/file/epi_model_mmt.RData")
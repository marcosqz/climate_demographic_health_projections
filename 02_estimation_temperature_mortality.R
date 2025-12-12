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

# Set parameters defining the DLNMs
dlnm_var <- list(
  var_prc = c(10, 75, 90),
  var_fun = "ns",
  lag_fun = "ns",
  max_lag = 21,
  lagnk = 3)

# Define the exposure- and lag-response functions
argvar <- list(fun = dlnm_var$var_fun,
               knots = quantile(data_tempmort$tmean, 
                                dlnm_var$var_prc/100, na.rm = TRUE),
               Bound = range(data_tempmort$tmean, na.rm = TRUE))
arglag <- list(fun = dlnm_var$lag_fun, 
               knots = logknots(dlnm_var$max_lag, nk = dlnm_var$lagnk))

# Build the crossbasis
cb <- crossbasis(data_tempmort$tmean, lag = dlnm_var$max_lag, argvar, arglag)

# Initialize objects to save the outputs
coef_age <- list()
vcov_age <- list()
mmt_age <- c()

# Loop age-specific temperature-mortality models
for(i in 1:length(study_param$age_groups)) {
  
  # Set the age group label
  i_age <- study_param$age_groups[i]
  
  # Run model
  model <- glm(
    as.formula(paste0("mort.", i_age,
      "~ cb + factor(dow) + ns(date, df = round(8 * length(date) / 365.25))")),
    data = data_tempmort, family = quasipoisson)
  
  # Reduce coefficients keeping the cumulative exposure-response function
  reduced <- crossreduce(cb, model, at = data_tempmort$tmean)
  
  # Find the minimum mortality temperature
  mmt_age[i_age] <- reduced$predvar[which.min(reduced$RRfit)]
  
  # Extract model coefficients and variances
  coef_age[[i_age]] <- coef(reduced)
  vcov_age[[i_age]] <- vcov(reduced)
  
  # Generate monte carlo samples to generate simulation of the ERF
  set.seed(13041975 + i) # Different random seed for each group
  coefsim_age <- t(mvrnorm(study_param$n_sim,
                           coef_age[[i_age]],
                           vcov_age[[i_age]]))     

  # Save estimated and sampled coefficients together
  coef_age[[i_age]] <- cbind(coef_age[[i_age]], coefsim_age)
  colnames(coef_age[[i_age]]) <- c("est", paste0("sim", 1:study_param$n_sim))
  
}; rm(i)

#### SAVE OUTPUTS ##############################################################

save(argvar, file = "outdata/file/01_epi_model/argvar.RData")
save(arglag, file = "outdata/file/01_epi_model/arglag.RData")
save(dlnm_var, file = "outdata/file/01_epi_model/dlnm_var.RData")
save(coef_age, file = "outdata/file/01_epi_model/coef_age.RData")
save(vcov_age, file = "outdata/file/01_epi_model/vcov_age.RData")
save(mmt_age, file = "outdata/file/01_epi_model/mmt_age.RData")
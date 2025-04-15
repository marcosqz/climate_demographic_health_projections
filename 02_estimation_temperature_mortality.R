#### LOAD LIBRARIES ####

library(dlnm) # logknots
library(splines) # ns
library(MASS) # mvrnorm

#### LOAD DATA

load("../../outdata/file/data_tempmort.RData")

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
argvar <- list(fun = dlnm_var$var_fun, # TODO: Why masselot's uses "bs", it is not necessary the "ns" to the log-liner extrapolation
               knots = quantile(data_tempmort$tmean, 
                                dlnm_var$var_prc/100, na.rm = TRUE),
               Bound = range(data_tempmort$tmean, na.rm = TRUE))
arglag <- list(fun = dlnm_var$lag_fun, 
               knots = logknots(dlnm_var$max_lag, nk = dlnm_var$lagnk))

# BUILD THE CROSS-BASIS
cb <- crossbasis(data_tempmort$tmean, lag = maxlag, argvar, arglag)

# VISUALLITATION PARAMETERS BY AGE GROUPS
age_parameters <- data.frame(
  response = c("mort.00_74", "mort.75plus"),
  groups = c("00_74", "75plus"),
  col_est = c("#1B9E77", "#D95F02"),
  title_plot = c("<74 yo", "+75 yo"))

coefsim_age <- list()
mmt_age <- c()

# RUN AGE-SPECIFIC MODELS
for(i in 1:2) {
  
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
  set.seed(13041975) # TODO: Check this! I think the set.see it should be the same but I'm not sure
  coefsim <- t(mvrnorm(100, coef, vcov))
  coefsim_age[[age_parameters$groups[i]]] <- coefsim
  
  # Temperature percentiles
  predper <- c(seq(0, 1, 0.1), 2:98, seq(99, 100, 0.1))
  tper <- quantile(data_tempmort$tmean, predper / 100)
  tper_extra <- seq(max(tper), max(tper) + 10, length = 50)
  tper <- c(tper, tper_extra[-1])
  
  # Exposure_response basis centred in the corresponding mmt by age groups
  cenvec <- onebasis(mmt, fun = argvar$fun, knots = argvar$knots, 
                     Boundary.knots = argvar$Bound)
  
  bcen <- scale(onebasis(tper, fun = argvar$fun, knots = argvar$knots, 
                         Boundary.knots = argvar$Bound), 
                center = cenvec, scale = FALSE)
  
  # Calculate the samples of relative risks 
  rrsim <- exp(bcen %*% coefsim)
  
  # Calculate the point-estimate of relative risks 
  rrest <- exp(bcen %*% coef)[,1]
  
  # Plot first age group
  plot(tper, rrsim[,1], type = "n",
       ylim = c(0.95, 4), log = "y",
       xlab = "Temperature (ºC)",
       ylab = "Relative risk",
       main = age_parameters$title_plot[i],
       xlim = c(min(tper), max(tper)))
  for(j in 1:100) {lines(tper, rrsim[,j], lwd = 0.5, col = "lightgrey")}
  lines(tper, rrest, lwd = 3, col = age_parameters$col_est[i])
  abline(v = max(data_tempmort$tmean), col = "black", lty = 2)
  abline(h = 1)
  text(max(data_tempmort$tmean)-0, 3.9, "observed \ntemperatures \n<-----", 
       pos = 2, srt = 0, cex = 0.7)
  text(max(data_tempmort$tmean)+0, 3.9, "previously \nunobserved \ntemperatures \n----->", 
       pos = 4, srt = 0, cex = 0.7)
  
}

# SAVE DATASETS
save(argvar, file = "../../outdata/file/argvar.RData")
save(arglag, file = "../../outdata/file/arglag.RData")
save(dlnm_var, file = "../../outdata/file/epi_model_dlnm_varibles.RData")
save(coefsim_age, file = "../../outdata/file/epi_model_coefsimage.RData")
save(mmt_age, file = "../../outdata/file/epi_model_mmt.RData")
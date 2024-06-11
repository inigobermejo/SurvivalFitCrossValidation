library(survival)
library(flexsurv)
library(readxl)
library(tidyverse)

# Breast cancer data sets used in Royston and Altman (2013)
# German Breast Cancer Study Group
nrow(gbsg)
# 686
KM.gbsg <- survfit(Surv(rfstime, status) ~ 1, data = gbsg)
plot(KM.gbsg, conf.int= T, mark.time= T, xlim = c(0,3000))

# Monoclonal gammopathy data
# Natural history of 1341 sequential patients with monoclonal gammopathy of undetermined significance (MGUS).
nrow(mgus2)
colnames(mgus2)
# 1384
KM.mgus2 <- survfit(Surv(futime, death) ~ 1, data = mgus2)
plot(KM.mgus2, conf.int= T, mark.time= T)

# Acute myeloid leukemia
# Simulated
nrow(myeloid)
# 646

# Survival times of patients with multiple myeloma
nrow(myeloma)
# 3882
colnames(myeloma)
KM.myeloma <- survfit(Surv(futime, death) ~ 1, data = myeloma)
plot(KM.myeloma, conf.int= T, mark.time= T)

# Breast cancer in the Rotterdam tumour bank
nrow(rotterdam)
# 2982
colnames(rotterdam)
KM.rotterdam <- survfit(Surv(dtime, death) ~ 1, data = rotterdam)
plot(KM.rotterdam, conf.int= T, mark.time= T)

# Subjects on a liver transplant waiting list from 1990-1999, and their disposition: received a transplant
# transplant
nrow(transplant)
# 815
colnames(transplant)

df_transplant <- transplant
df_transplant$status <- ifelse(df_transplant$event == 'ltx', 1,0)
df_transplant <- df_transplant[df_transplant$futime>0,]
df_transplant <- df_transplant[df_transplant$futime<1200,]
KM.transplant <- survfit(Surv(futime, status) ~ 1, data = df_transplant)
plot(KM.transplant, conf.int= T, mark.time= T)

setwd("~/GitHub/SurvSet/SurvSet/_datagen/output")

data.ovarian <- read.csv('dataOvarian1.csv')
data.ovarian <- data.ovarian[data.ovarian$time > 0,]
colnames(data.ovarian)
nrow(data.ovarian)
KM.data.ovarian <- survfit(Surv(time, event) ~ 1, data = data.ovarian)
plot(KM.data.ovarian, conf.int= T, mark.time= T)


setwd("~/Data/SurvCV")
df_tcga <- readxl::read_excel("TCGA-CDR-SupplementalTableS1.xlsx", sheet = 'TCGA-CDR')

df_tcga_gbm <- df_tcga[df_tcga$type=='GBM',]
df_tcga_gbm$status <- ifelse(df_tcga_gbm$vital_status == 'Dead', 1,0)
df_tcga_gbm <- df_tcga_gbm[df_tcga_gbm$OS.time>0,]
KM.tcga_gbm <- survfit(Surv(OS.time, ifelse(vital_status == 'Dead', 1,0)) ~ 1, data = df_tcga_gbm)
plot(KM.tcga_gbm, conf.int= T, mark.time= T)

#colnames(df_tcga)
#seer.breast <- read.csv('SEER Breast Cancer Dataset .csv')
#KM.seer <- survfit(Surv(Survival.Months, ifelse(Status == 'Dead', 1,0)) ~ 1, data = seer.breast)
#plot(KM.seer, conf.int= T, mark.time= T)

sample_censor <- function(df, n, time = 'time', status = 'status')
{
  df$time   <-   df[[time]]  
  df$status <-   df[[status]]
  df_selected <- df[sample(nrow(df),n),]
  KM.df <- survfit(Surv(time, status) ~ 1, data = df_selected)
  time_censor <- KM.df$time[min(which(KM.df$surv < 0.5))]
  df_selected$status <- ifelse(df_selected$time > time_censor, 0, df_selected$status)
  df_selected$time <- ifelse(df_selected$time > time_censor, time_censor, df_selected$time)
  return(df_selected)
}

surv_fit <- function (KM.data, time = 'time', status = 'status')
{
  KM.data$time   <-   KM.data[[time]]  
  KM.data$status <-   KM.data[[status]]
  

  fit.llogis    <- flexsurvreg(Surv(time, status) ~ 1, data = KM.data, dist = "llogis" )       # fit model with loglogistic distribution
  fit.weib      <- flexsurvreg(Surv(time, status) ~ 1, data = KM.data, dist = "weibull")       # fit model with Weibull distribution
  fit.lnorm     <- flexsurvreg(Surv(time, status) ~ 1, data = KM.data, dist = "lnorm"  )       # fit model with lognormal distribution
  fit.gamma     <- flexsurvreg(Surv(time, status) ~ 1, data = KM.data, dist = "gamma"  )       # fit model with gamma distribution 
  fit.exp       <- flexsurvreg(Surv(time, status) ~ 1, data = KM.data, dist = "exp"    )       # fit model with exponential distribution
  fit.gompertz  <- flexsurvreg(Surv(time, status) ~ 1, data = KM.data, dist = "gompertz" , inits=c(-1, 1/mean(KM.data$time)))    # fit model with gompertz  
  
  spline_model_k_1 <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 1) , error = function(e) NULL)
  spline_model_k_2 <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 2) , error = function(e) NULL)
  spline_model_k_3 <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 3) , error = function(e) NULL)
  spline_model_k_4 <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 4) , error = function(e) NULL)
  
  models <- list( fit.llogis, fit.weib, fit.lnorm, fit.gamma, fit.exp, fit.gompertz,
               spline_model_k_1, spline_model_k_2, spline_model_k_3, spline_model_k_4)

  modelnames <- c('loglogistic', 'weibull', 'lognormal', 'gamma', 'exponential', 'gompertz',
                  'spline k=1', 'spline k=2', 'spline k=3', 'spline k=4')
  
  AIC <- c(    AIC(fit.llogis),                                         
               AIC(fit.weib), 
               AIC(fit.lnorm), 
               AIC(fit.gamma),
               AIC(fit.exp),
               AIC(fit.gompertz),
               tryCatch(AIC(spline_model_k_1), error = function(e) Inf),
               tryCatch(AIC(spline_model_k_2), error = function(e) Inf),
               tryCatch(AIC(spline_model_k_3), error = function(e) Inf),
               tryCatch(AIC(spline_model_k_4), error = function(e) Inf))
  
  
  # compare BIC values
  BIC <- c(    BIC(fit.llogis),                                         
               BIC(fit.weib), 
               BIC(fit.lnorm), 
               BIC(fit.gamma),
               BIC(fit.exp),
               BIC(fit.gompertz),
               tryCatch(BIC(spline_model_k_1), error = function(e) Inf),
               tryCatch(BIC(spline_model_k_2), error = function(e) Inf),
               tryCatch(BIC(spline_model_k_3), error = function(e) Inf),
               tryCatch(BIC(spline_model_k_4), error = function(e) Inf))
  

  return (list (Models = models, Modelnames = modelnames, Metrics = data.frame(AIC, BIC)))
}


auc_km <- function (KM.data)
{
  time_diff <- c(0, diff(KM.data$time))
  
  # Calculate the area under the curve using the trapezoidal rule
  auc <- sum(KM.data$surv * time_diff)
  
  return(auc)
}


cross_validation <- function(df, k = 10)
{
    metrics <- list()
    # Shuffle dataset
    df <- df[sample(nrow(df)),]
    # Determine fold size
    fold_size <- floor(nrow(df) / k)
    
    for(i in 1:k)
    {
      # Separate training and testing set
      test_set <- ((i-1)*fold_size + 1):(i*fold_size)
      df_test <- df[test_set,]
      df_train <- df[-test_set,]
      
      # Fit the logistic regression
      results_fit <- surv_fit(df_train)
      metrics[[i]] <- assess_fit <- assess_fit(results_fit$Models, df_test)
    }
    return(metrics)
}

assess_fit <- function(models, data)
{
  n <- length(models)
  metrics <- data.frame(aic = rep(0, n), bic = rep(0, n))
  for(i in 1:n)
  {
    model <- models[[i]]
    predictions <- predict(model, newdata = data, times = unique(data$time), type = "survival")
    hazard <- predict(model, newdata = data, times = unique(data$time), type = "hazard")
    # Extract the .pred_survival for the exact times in test_times
    p_surv <- map2_dbl(predictions$.pred, data$time, ~ .x %>% filter(.time == .y) %>% pull(.pred_survival))
    p_haz <- map2_dbl(hazard$.pred, data$time, ~ .x %>% filter(.time == .y) %>% pull(.pred_hazard))
    log_lik <- sum(data$status * log(p_haz) + log(p_surv))
    num_params <- length(coef(model))

    metrics$aic[i] <- -2 * log_lik + 2 * num_params
    metrics$bic[i] <- -2 * log_lik + log(model$N) * num_params
  }
  return(metrics)
}

process_dataset <- function(df, time = 'time', status = 'status', n = 250, k = 100)
{
  max_time = max(df[, time])

  results_full <- surv_fit(df, time, status)
  best_model_name_aic <- results_full$Modelnames[which.min(results_full$Metrics$AIC)]
  best_model_name_bic <- results_full$Modelnames[which.min(results_full$Metrics$BIC)]
  
  best_model_aic <- results_full$Models[which.min(results_full$Metrics$AIC)][[1]]
  best_model_bic <- results_full$Models[which.min(results_full$Metrics$BIC)][[1]]
  
  rmst_best_aic_full <- summary(best_model_aic, type='rmst')[[1]]$est
  rmst_best_bic_full <- summary(best_model_bic, type='rmst')[[1]]$est

  sampling_results <- list()
  rmst_trad_aic <- c()
  rmst_trad_bic <- c()
  rmst_cv_aic <- c()
  rmst_cv_bic <- c()
  
  for(i in 1:k)
  {
    df_sample        <- sample_censor(df, n, time, status)
    results_trad     <- surv_fit(df_sample, time, status)

    best_model_name_trad_aic   <- results_trad$Modelnames[which.min(results_trad$Metrics$AIC)][[1]]
    best_model_name_trad_bic   <- results_trad$Modelnames[which.min(results_trad$Metrics$BIC)][[1]]
    
    best_model_aic   <- results_trad$Models[which.min(results_trad$Metrics$AIC)][[1]]
    best_model_bic   <- results_trad$Models[which.min(results_trad$Metrics$BIC)][[1]]
    
    rmst_trad_aic[i] <- summary(best_model_aic, type='rmst', t = max_time)[[1]]$est
    rmst_trad_bic[i] <- summary(best_model_bic, type='rmst', t = max_time)[[1]]$est
    
    results_cv       <- cross_validation(df_sample, k = 5)
    results_cv_avg   <- Reduce("+", results_cv) / length(results_cv)
    
    best_model_name_cv_aic   <- results_trad$Modelnames[which.min(results_cv_avg$aic)][[1]]
    best_model_name_cv_bic   <- results_trad$Modelnames[which.min(results_cv_avg$bic)][[1]]

    best_model_aic   <- results_trad$Models[which.min(results_cv_avg$aic)][[1]]
    best_model_bic   <- results_trad$Models[which.min(results_cv_avg$bic)][[1]]

    rmst_cv_aic[i]   <- summary(best_model_aic, type='rmst', t = max_time)[[1]]$est
    rmst_cv_bic[i]   <- summary(best_model_bic, type='rmst', t = max_time)[[1]]$est
    
    sampling_results[[i]] <- list(cv = results_cv_avg, trad = results_trad$Metrics)

    print(paste0(i,":", "AIC: full:",  best_model_name_aic,
                 "; trad:", best_model_name_trad_aic, 
                 "; cv:", best_model_name_cv_aic, 
                 "; BIC:", best_model_name_bic,
                 "; trad:", best_model_name_trad_bic, 
                 "; cv:", best_model_name_cv_bic))
    
    print(paste0(i,": AIC",
                 ": trad:", round(mean(abs(rmst_trad_aic - rmst_best_aic_full)),3),
                 "; cv:",   round(mean(abs(rmst_cv_aic - rmst_best_aic_full)),3), 
                 "; BIC",
                 ": trad:", round(mean(abs(rmst_trad_bic - rmst_best_bic_full)),3), 
                 "; cv:",   round(mean(abs(rmst_cv_bic - rmst_best_bic_full)), 3)))
  }
  
  error_aic_trad <-  mean(abs(rmst_trad_aic - rmst_best_aic_full))
  error_aic_cv   <-  mean(abs(rmst_cv_aic - rmst_best_aic_full))
  error_bic_trad <-  mean(abs(rmst_trad_bic - rmst_best_bic_full))
  error_bic_cv   <-  mean(abs(rmst_cv_bic - rmst_best_bic_full))
  
  return(list(sampling_results,rmst_best_aic_full, 
              rmst_best_bic_full, 
              rmst_trad_aic, 
              rmst_trad_bic, 
              rmst_cv_aic, 
              rmst_cv_bic,
              error_aic_trad,
              error_aic_cv,
              error_bic_trad,
              error_bic_cv))
}

plot_KMs<-function(df, time = 'time', status = 'status', n = 250, cv = TRUE, second_sample = FALSE)
{
  df$time   <-   df[,   time]  
  df$status <-   df[, status]
  times    <- seq(0, max(df$time) + 10, 1)

  KM.df <- survfit(Surv(time, status) ~ 1, data = df)
  plot(KM.df, conf.int= T, mark.time= T, xlim = c(0,(max(df$time) + 10)))
  fit_full            <- surv_fit(df, time, status)
  best_model_full     <- fit_full$Models[which.min(fit_full$Metrics$AIC)][[1]]
  lines(best_model_full,   t = times, ci = F, col = 'black')
  
  df_sample1 <- sample_censor(df, n, time, status)
  lines(survfit(Surv(time, status) ~ 1, data = df_sample1), conf.int= T, mark.time= T, col = 'red')
  fit_sample1         <- surv_fit(df_sample1, time, status)
  best_model_sample1  <- fit_sample1$Models[which.min(fit_sample1$Metrics$AIC)][[1]]
  lines(best_model_sample1,   t = times,ci = F, col = 'red')
  
  if(cv)
  {
    results_cv       <- cross_validation(df_sample1, k = 5)
    results_cv_avg   <- Reduce("+", results_cv) / length(results_cv)
    best_model_cv_aic   <- fit_sample1$Models[which.min(results_cv_avg$aic)][[1]]
    lines(best_model_cv_aic,   t = times,ci = F, col = 'blue')
  }
  
  if(second_sample)
  {
    df_sample2 <- sample_censor(df, n, time, status)
    lines(survfit(Surv(time, status) ~ 1, data = df_sample2), conf.int= T, mark.time= T, col = 'blue')
    fit_sample2         <- surv_fit(df_sample2, time, status)
    best_model_sample2  <- fit_sample2$Models[which.min(fit_sample2$Metrics$AIC)][[1]]
    lines(best_model_sample2,   t = times, ci = F, col = 'blue')
  }
  
}

results_gbsg       <- process_dataset(gbsg,          time = 'rfstime', status = 'status')
results_mgus2      <- process_dataset(mgus2,         time = 'futime',  status = 'death')
results_myeloma    <- process_dataset(myeloma,       time = 'futime',  status = 'death')
results_rotterdam  <- process_dataset(rotterdam,     time = 'dtime',   status = 'death')
results_transplant <- process_dataset(df_transplant, time = 'futime',  status = 'status')
results_ovarian    <- process_dataset(data.ovarian,  time = 'time',    status = 'event')
results_tcga_gbm   <- process_dataset(df_tcga_gbm,   time = 'OS.time', status = 'status')



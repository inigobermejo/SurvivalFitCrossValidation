library(survival)
library(flexsurv)
library(readxl)
library(tidyverse)

df_seer_selected <- read.csv('seer_selected.csv')


df_seer_breast     <- df_seer_selected %>%  subset(Site.recode.ICD.O.3.WHO.2008 == 'Breast')
df_seer_pancreas   <- df_seer_selected %>%  subset(Site.recode.ICD.O.3.WHO.2008 == 'Pancreas')
df_seer_colorectal <- df_seer_selected %>%  subset(Site.recode.ICD.O.3.WHO.2008 %in% c('Ascending Colon',
                                                                                       'Cecum',
                                                                                       'Hepatic Flexure',
                                                                                       'Rectosigmoid Junction',
                                                                                       'Sigmoid Colon',
                                                                                       'Splenic Flexure',
                                                                                       'Transverse Colon',
                                                                                       'Descending Colon',
                                                                                       'Rectum'))
df_seer_lung       <- df_seer_selected %>%  subset(Site.recode.ICD.O.3.WHO.2008 == 'Lung and Bronchus')
df_seer_SCLC       <- df_seer_lung     %>%  subset(Histologic.Type.ICD.O.3 %in% c("8041", "8042", "8043", "8044", "8045"))
df_seer_NSCLC      <- df_seer_lung     %>%  subset(!Histologic.Type.ICD.O.3 %in% c("8041", "8042", "8043", "8044", "8045"))

KM.seer_pancreas <- survfit(Surv(Survival.months, status) ~ 1, data = df_seer_pancreas)
plot(KM.seer_pancreas, conf.int= T, mark.time= T, xlim=c(0,300))


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

setwd("~/Inigo/SurvivalFitCrossValidation")

data.ovarian <- read.csv('dataOvarian1.csv')
data.ovarian <- data.ovarian[data.ovarian$time > 0,]
colnames(data.ovarian)
nrow(data.ovarian)
KM.data.ovarian <- survfit(Surv(time, event) ~ 1, data = data.ovarian)
plot(KM.data.ovarian, conf.int= T, mark.time= T)


setwd("~/Inigo/SurvivalFitCrossValidation")
df_tcga <- readxl::read_excel("TCGA-CDR-SupplementalTableS1.xlsx", sheet = 'TCGA-CDR')

df_tcga_gbm <- df_tcga[df_tcga$type=='GBM',]
df_tcga_gbm$status <- ifelse(df_tcga_gbm$vital_status == 'Dead', 1,0)
df_tcga_gbm <- df_tcga_gbm[!is.na(df_tcga_gbm$status),]
df_tcga_gbm <- df_tcga_gbm[df_tcga_gbm$OS.time>0,]
KM.tcga_gbm <- survfit(Surv(OS.time, ifelse(vital_status == 'Dead', 1,0)) ~ 1, data = df_tcga_gbm)
plot(KM.tcga_gbm, conf.int= T, mark.time= T)

sample_censor <- function(df, time = 'time', status = 'status', n = 250, cutoff = 0.5)
{
  df$time   <-   df[[time]]  
  df$status <-   df[[status]]
  df_selected <- df[sample(nrow(df),n),]
  KM.df <- survfit(Surv(time, status) ~ 1, data = df_selected)
  if(!is.infinite(min(which(KM.df$surv < cutoff))))
  {
    time_censor <- KM.df$time[min(which(KM.df$surv < cutoff))]
    df_selected$status <- ifelse(df_selected$time > time_censor, 0, df_selected$status)
    df_selected$time <- ifelse(df_selected$time > time_censor, time_censor, df_selected$time)
  }
  return(df_selected)
}

surv_fit <- function (KM.data, time = 'time', status = 'status')
{
  KM.data$time   <-   KM.data[[time]]  
  KM.data$status <-   KM.data[[status]]
  
  fit.llogis    <- flexsurvreg(Surv(time, status) ~ 1, data = KM.data, dist = "llogis" )       # fit model with loglogistic distribution
  fit.weib      <- flexsurvreg(Surv(time, status) ~ 1, data = KM.data, dist = "weibull")       # fit model with Weibull distribution
  fit.lnorm     <- flexsurvreg(Surv(time, status) ~ 1, data = KM.data, dist = "lnorm"  )       # fit model with lognormal distribution
  fit.gamma     <- tryCatch(flexsurvreg(Surv(time, status) ~ 1, data = KM.data, dist = "gamma")   , error = function(e) NULL) # fit model with gamma distribution 
  fit.exp       <- tryCatch(flexsurvreg(Surv(time, status) ~ 1, data = KM.data, dist = "exp")     , error = function(e) NULL)       # fit model with exponential distribution
  fit.gompertz  <- tryCatch(flexsurvreg(Surv(time, status) ~ 1, data = KM.data, dist = "gompertz" , inits=c(-1, 1/mean(KM.data$time))) , error = function(e) NULL)    # fit model with gompertz  
  fit.gengamma  <- tryCatch(flexsurvreg(Surv(time, status) ~ 1, data = KM.data, dist = "gengamma"), error = function(e) NULL)       # fit model with generalized gamma distribution 
  
  spline_model_k_1        <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 1) , error = function(e) NULL)
  spline_model_k_2        <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 2) , error = function(e) NULL)
  spline_model_k_3        <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 3) , error = function(e) NULL)
  spline_model_k_4        <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 4) , error = function(e) NULL)
  spline_model_k_1_odds   <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 1, scale = 'odds') , error = function(e) NULL)
  spline_model_k_2_odds   <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 2, scale = 'odds') , error = function(e) NULL)
  spline_model_k_3_odds   <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 3, scale = 'odds') , error = function(e) NULL)
  spline_model_k_4_odds   <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 4, scale = 'odds') , error = function(e) NULL)
  spline_model_k_1_normal <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 1, scale = 'normal') , error = function(e) NULL)
  spline_model_k_2_normal <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 2, scale = 'normal') , error = function(e) NULL)
  spline_model_k_3_normal <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 3, scale = 'normal') , error = function(e) NULL)
  spline_model_k_4_normal <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = KM.data, k = 4, scale = 'normal') , error = function(e) NULL)
  
  models <- list( fit.llogis, fit.weib, fit.lnorm, fit.gamma, fit.exp, fit.gompertz, fit.gengamma,
               spline_model_k_1, spline_model_k_2, spline_model_k_3, spline_model_k_4,
               spline_model_k_1_odds, spline_model_k_2_odds, spline_model_k_3_odds, spline_model_k_4_odds,
               spline_model_k_1_normal, spline_model_k_2_normal, spline_model_k_3_normal, spline_model_k_4_normal)

  modelnames <- c('loglogistic', 'weibull', 'lognormal', 'gamma', 'exponential', 'gompertz', 'generalised gamma',
                  'spline k=1 hazard', 'spline k=2 hazard', 'spline k=3 hazard', 'spline k=4 hazard',
                  'spline k=1 odds',   'spline k=2 odds',   'spline k=3 odds',   'spline k=4 odds',
                  'spline k=1 normal', 'spline k=2 normal', 'spline k=3 normal', 'spline k=4 normal')
  
  aic <- c(    AIC(fit.llogis),                                         
               AIC(fit.weib), 
               AIC(fit.lnorm), 
               tryCatch(AIC(fit.gamma),               error = function(e) Inf),
               tryCatch(AIC(fit.exp),                 error = function(e) Inf),
               tryCatch(AIC(fit.gompertz),            error = function(e) Inf),
               tryCatch(AIC(fit.gengamma),            error = function(e) Inf),
               tryCatch(AIC(spline_model_k_1),        error = function(e) Inf),
               tryCatch(AIC(spline_model_k_2),        error = function(e) Inf),
               tryCatch(AIC(spline_model_k_3),        error = function(e) Inf),
               tryCatch(AIC(spline_model_k_4),        error = function(e) Inf),
               tryCatch(AIC(spline_model_k_1_odds),   error = function(e) Inf),
               tryCatch(AIC(spline_model_k_2_odds),   error = function(e) Inf),
               tryCatch(AIC(spline_model_k_3_odds),   error = function(e) Inf),
               tryCatch(AIC(spline_model_k_4_odds),   error = function(e) Inf),
               tryCatch(AIC(spline_model_k_1_normal), error = function(e) Inf),
               tryCatch(AIC(spline_model_k_2_normal), error = function(e) Inf),
               tryCatch(AIC(spline_model_k_3_normal), error = function(e) Inf),
               tryCatch(AIC(spline_model_k_4_normal), error = function(e) Inf)
  )
  
  
  # compare BIC values
  bic <- c(    BIC(fit.llogis),                                         
               BIC(fit.weib), 
               BIC(fit.lnorm), 
               tryCatch(BIC(fit.gamma),               error = function(e) Inf),
               tryCatch(BIC(fit.exp),                 error = function(e) Inf),
               tryCatch(BIC(fit.gompertz),            error = function(e) Inf),
               tryCatch(BIC(fit.gengamma),            error = function(e) Inf),
               tryCatch(BIC(spline_model_k_1),        error = function(e) Inf),
               tryCatch(BIC(spline_model_k_2),        error = function(e) Inf),
               tryCatch(BIC(spline_model_k_3),        error = function(e) Inf),
               tryCatch(BIC(spline_model_k_4),        error = function(e) Inf),
               tryCatch(BIC(spline_model_k_1_odds),   error = function(e) Inf),
               tryCatch(BIC(spline_model_k_2_odds),   error = function(e) Inf),
               tryCatch(BIC(spline_model_k_3_odds),   error = function(e) Inf),
               tryCatch(BIC(spline_model_k_4_odds),   error = function(e) Inf),
               tryCatch(BIC(spline_model_k_1_normal), error = function(e) Inf),
               tryCatch(BIC(spline_model_k_2_normal), error = function(e) Inf),
               tryCatch(BIC(spline_model_k_3_normal), error = function(e) Inf),
               tryCatch(BIC(spline_model_k_4_normal), error = function(e) Inf)
  )
  

  return (list (Models = models, Modelnames = modelnames, Metrics = data.frame(aic, bic)))
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
      
      # Fit the curves
      results_fit <- surv_fit(df_train)
      metrics[[i]]  <- assess_fit(results_fit$Models, df_test)
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
    if(is.null(model))
    {
      #print(paste0("Warning: Model " ,i, " is null"))
      metrics$aic[i] <- Inf
      metrics$bic[i] <- Inf
      
    }else
    {
      predictions <- predict(model, newdata = data, times = unique(data$time), type = "survival")
      hazard <- predict(model, newdata = data, times = unique(data$time), type = "hazard")
      # Extract the .pred_survival for the exact times in test_times
      p_surv <- map2_dbl(predictions$.pred, data$time, ~ .x %>% filter(.eval_time == .y) %>% pull(.pred_survival))
      p_haz <- map2_dbl(hazard$.pred, data$time, ~ .x %>% filter(.eval_time == .y) %>% pull(.pred_hazard))
      log_lik <- sum(data$status * log(p_haz) + log(p_surv))
      num_params <- length(coef(model))
  
      metrics$aic[i] <- -2 * log_lik + 2 * num_params
      metrics$bic[i] <- -2 * log_lik + log(model$N) * num_params
      
    }
  }
  return(metrics)
}


plot_KMs<-function(df, time = 'time', status = 'status', n = 250, k = 5, cv = TRUE, second_sample = FALSE, cutoff = 0.5)
{
  df$time   <-   df[,   time]  
  df$status <-   df[, status]
  times    <- seq(0, max(df$time) + 10, 1)
  
  KM.df <- survfit(Surv(time, status) ~ 1, data = df)
  plot(KM.df, conf.int= T, mark.time= T, xlim = c(0,(max(df$time) + 10)))
  fit_full            <- surv_fit(df, time, status)
  best_model_full     <- fit_full$Models[which.min(fit_full$Metrics$aic)][[1]]
  lines(best_model_full,   t = times, ci = F, col = 'black')
  
  df_sample1 <- sample_censor(df, time, status, n, cutoff)
  lines(survfit(Surv(time, status) ~ 1, data = df_sample1), conf.int= T, mark.time= T, col = 'red')
  fit_sample1         <- surv_fit(df_sample1)
  best_model_sample1  <- fit_sample1$Models[which.min(fit_sample1$Metrics$aic)][[1]]
  lines(best_model_sample1,   t = times,ci = F, col = 'red')
  
  if(cv)
  {
    results_cv       <- cross_validation(df_sample1, k = k)
    results_cv_avg   <- Reduce("+", results_cv) / length(results_cv)
    best_model_cv_aic   <- fit_sample1$Models[which.min(results_cv_avg$aic)][[1]]
    lines(best_model_cv_aic,   t = times,ci = F, col = 'blue')
  }
  
  if(second_sample)
  {
    df_sample2 <- sample_censor(df, time, status, n, cutoff)
    lines(survfit(Surv(time, status) ~ 1, data = df_sample2), conf.int= T, mark.time= T, col = 'blue')
    fit_sample2         <- surv_fit(df_sample2, time, status)
    best_model_sample2  <- fit_sample2$Models[which.min(fit_sample2$Metrics$aic)][[1]]
    lines(best_model_sample2,   t = times, ci = F, col = 'blue')
  }
  
}

process_dataset <- function(df, time = 'time', status = 'status', n = 250, m = 100, k = 5, cutoff = 0.5)
{
  max_time = max(df[, time])

  results_full <- surv_fit(df, time, status)
  best_model_name_aic <- results_full$Modelnames[which.min(results_full$Metrics$aic)]
  best_model_name_bic <- results_full$Modelnames[which.min(results_full$Metrics$bic)]
  
  best_model_aic <- results_full$Models[which.min(results_full$Metrics$aic)][[1]]
  best_model_bic <- results_full$Models[which.min(results_full$Metrics$bic)][[1]]
  
  rmst_best_aic_full <- summary(best_model_aic, type='rmst')[[1]]$est
  rmst_best_bic_full <- summary(best_model_bic, type='rmst')[[1]]$est

  sampling_results <- list()
  sample_rmst <- c()
  rmst_trad_aic <- c()
  rmst_trad_bic <- c()
  rmst_cv_aic <- c()
  rmst_cv_bic <- c()
  best_model_name_trad_aic <- c()
  best_model_name_trad_bic <- c()
  best_model_name_cv_aic <- c()
  best_model_name_cv_bic <- c()
  
  for(i in 1:m)
  {
    df_sample        <- sample_censor(df, time, status, n, cutoff)
    results_trad     <- surv_fit(df_sample)
    sample_rmst[[i]] <- sapply(results_trad$Models, function(x){ return(tryCatch(summary(x, type='rmst', t = max_time)[[1]]$est, error = function(e) NA))})

    best_model_name_trad_aic[i]   <- results_trad$Modelnames[which.min(results_trad$Metrics$aic)]
    best_model_name_trad_bic[i]   <- results_trad$Modelnames[which.min(results_trad$Metrics$bic)]
    
    best_model_aic   <- results_trad$Models[which.min(results_trad$Metrics$aic)][[1]]
    best_model_bic   <- results_trad$Models[which.min(results_trad$Metrics$bic)][[1]]
    
    rmst_trad_aic[i] <- summary(best_model_aic, type='rmst', t = max_time)[[1]]$est
    rmst_trad_bic[i] <- summary(best_model_bic, type='rmst', t = max_time)[[1]]$est
    
    results_cv       <- cross_validation(df_sample, k = k)
    results_cv_avg   <- Reduce("+", results_cv) / length(results_cv)
    
    best_model_name_cv_aic[i]   <- results_trad$Modelnames[which.min(results_cv_avg$aic)]
    best_model_name_cv_bic[i]   <- results_trad$Modelnames[which.min(results_cv_avg$bic)]

    best_model_aic   <- results_trad$Models[which.min(results_cv_avg$aic)][[1]]
    best_model_bic   <- results_trad$Models[which.min(results_cv_avg$bic)][[1]]

    rmst_cv_aic[i]   <- summary(best_model_aic, type='rmst', t = max_time)[[1]]$est
    rmst_cv_bic[i]   <- summary(best_model_bic, type='rmst', t = max_time)[[1]]$est
    
    sampling_results[[i]] <- list(cv = results_cv_avg, trad = results_trad$Metrics)

    # print(paste0(i,": AIC",
    #              ": full:", best_model_name_aic,
    #              "; trad:", best_model_name_trad_aic[i], 
    #              "; cv:",   best_model_name_cv_aic[i], 
    #              "; BIC:",  best_model_name_bic,
    #              "; trad:", best_model_name_trad_bic[i], 
    #              "; cv:",   best_model_name_cv_bic[i]))
    
    print(paste0(i, ": n:", n, ";k:", k, "; AIC",
                 ": trad:", round(mean(abs(rmst_trad_aic - rmst_best_aic_full)), 2),
                 "; cv:",   round(mean(abs(rmst_cv_aic   - rmst_best_aic_full)), 2), 
                 "; BIC",
                 ": trad:", round(mean(abs(rmst_trad_bic - rmst_best_bic_full)), 2), 
                 "; cv:",   round(mean(abs(rmst_cv_bic   - rmst_best_bic_full)), 2)))
  }
  
  error_aic_trad <-  mean(abs(rmst_trad_aic - rmst_best_aic_full))
  error_aic_cv   <-  mean(abs(rmst_cv_aic - rmst_best_aic_full))
  error_bic_trad <-  mean(abs(rmst_trad_bic - rmst_best_bic_full))
  error_bic_cv   <-  mean(abs(rmst_cv_bic - rmst_best_bic_full))
  
   
   df_rmst_all           <- as.data.frame(do.call(rbind, sample_rmst))
   colnames(df_rmst_all) <- results_full$Modelnames
  
  return(list(metrics                  = sampling_results,
              rmst                     = df_rmst_all,
              rmst_best_aic_full       = rmst_best_aic_full, 
              rmst_best_bic_full       = rmst_best_bic_full, 
              rmst_trad_aic            = rmst_trad_aic, 
              rmst_trad_bic            = rmst_trad_bic, 
              rmst_cv_aic              = rmst_cv_aic, 
              rmst_cv_bic              = rmst_cv_bic,
              best_model_name_trad_aic = best_model_name_trad_aic,
              best_model_name_trad_bic = best_model_name_trad_bic,
              best_model_name_cv_aic   = best_model_name_cv_aic,
              best_model_name_cv_bic   = best_model_name_cv_bic,
              error_aic_trad           = error_aic_trad,
              error_aic_cv             = error_aic_cv,
              error_bic_trad           = error_bic_trad,
              error_bic_cv             = error_bic_cv))
}


process_datasets <- function(n = 250, m = 100, k = 5, cutoff = 0.5)
{
  print("------------------- SEER - Breast")
  results_seer_breast     <- process_dataset(df_seer_breast,     time = 'Survival.months', status = 'status', m = m, n = n, k = k, cutoff = cutoff)
  print("------------------- SEER - Pancreas") 
  results_seer_pancreas   <- process_dataset(df_seer_pancreas,   time = 'Survival.months', status = 'status', m = m, n = n, k = k, cutoff = cutoff)
  print("------------------- SEER - Colorectal")
  results_seer_colorectal <- process_dataset(df_seer_colorectal, time = 'Survival.months', status = 'status', m = m, n = n, k = k, cutoff = cutoff)
  print("------------------- SEER - SCLC")
  results_seer_sclc       <- process_dataset(df_seer_SCLC,       time = 'Survival.months', status = 'status', m = m, n = n, k = k, cutoff = cutoff)
  print("------------------- SEER - NSCLC")
  results_seer_nsclc      <- process_dataset(df_seer_NSCLC,      time = 'Survival.months', status = 'status', m = m, n = n, k = k, cutoff = cutoff)
  print("------------------- gbsg")
  results_gbsg            <- process_dataset(gbsg,               time = 'rfstime', status = 'status', m = m, n = n, k = k, cutoff = cutoff)
  print("------------------- mgus2")
  results_mgus2           <- process_dataset(mgus2,              time = 'futime',  status = 'death' , m = m, n = n, k = k, cutoff = cutoff)
  print("------------------- myeloma")
  results_myeloma         <- process_dataset(myeloma,            time = 'futime',  status = 'death' , m = m, n = n, k = k, cutoff = cutoff)
  print("------------------- rotterdam")
  results_rotterdam       <- process_dataset(rotterdam,          time = 'dtime',   status = 'death' , m = m, n = n, k = k, cutoff = cutoff)
  print("------------------- transplant")
  results_transplant      <- process_dataset(df_transplant,      time = 'futime',  status = 'status', m = m, n = n, k = k, cutoff = cutoff)
  print("------------------- ovarian")
  results_ovarian         <- process_dataset(data.ovarian,       time = 'time',    status = 'event' , m = m, n = n, k = k, cutoff = cutoff)
  print("------------------- TCGA GBM")
  results_tcga_gbm        <- process_dataset(df_tcga_gbm,        time = 'OS.time', status = 'status', m = m, n = n, k = k, cutoff = cutoff)
  
  return(list(seer_breast     = results_seer_breast,
              seer_pancreas   = results_seer_pancreas,
              seer_colorectal = results_seer_colorectal,
              seer_nsclc      = results_seer_nsclc,
              seer_sclc       = results_seer_sclc,
              gbsg            = results_gbsg,
              mgus2           = results_mgus2,
              myeloma         = results_myeloma,
              rotterdam       = results_rotterdam,
              transplant      = results_transplant,
              ovarian         = results_ovarian,
              tcga_gbm        = results_tcga_gbm))
  
}

results_n_250_k_10 <-process_datasets(n = 250, k = 10)
saveRDS(results_n_250_k_10, file = paste0('results_n_250_k_10', '.rds'))
results_n_150_k_10 <- process_datasets(n = 150, k = 10)
saveRDS(results_n_150_k_10, file = paste0('results_n_150_k_10', '.rds'))
results_n_350_k_10 <- process_datasets(n = 350, k = 10)
saveRDS(results_n_350_k_10, file = paste0('results_n_350_k_10', '.rds'))

results_n_250_k_10_co_04 <-process_datasets(n = 250, k = 10, cutoff = 0.4)
saveRDS(results_n_250_k_10_co_04, file = paste0('results_n_250_k_10_cutoff_0.4', '.rds'))
results_n_250_k_10_co_03 <-process_datasets(n = 250, k = 10, cutoff = 0.3)
saveRDS(results_n_250_k_10_co_03, file = paste0('results_n_250_k_10_cutoff_0.3', '.rds'))

results_n_150_k_5 <- process_datasets(n = 150, k = 5)
saveRDS(results_n_150_k_5, file = paste0('results_n_150_k_5', '.rds'))
results_n_250_k_5 <- process_datasets(n = 250, k = 5)
saveRDS(results_n_250_k_5, file = paste0('results_n_250_k_5', '.rds'))
results_n_350_k_5 <- process_datasets(n = 350, k = 5)
saveRDS(results_n_350_k_5, file = paste0('results_n_350_k_5', '.rds'))



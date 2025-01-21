# Load necessary libraries
library(survival)
library(rlang)
library(flexsurv)
library(ggplot2)

boostrapped_survival_curves_comparison <- function(km_fit, bootstrapped_kms, fit1_results, fit2_results)
{  

  km_matrix   <- do.call(cbind, bootstrapped_kms)
  fit1_matrix <- do.call(cbind, fit1_results)
  fit2_matrix <- do.call(cbind, fit2_results)
  
  # Calculate percentiles for confidence intervals
  km_ci_lower   <- apply(km_matrix, 1, quantile, probs = 0.025)
  km_ci_upper   <- apply(km_matrix, 1, quantile, probs = 0.975)
  fit1_ci_lower <- apply(fit1_matrix, 1, quantile, probs = 0.025)
  fit1_median   <- apply(fit1_matrix, 1, quantile, probs = 0.5)
  fit1_ci_upper <- apply(fit1_matrix, 1, quantile, probs = 0.975)
  fit2_ci_lower <- apply(fit2_matrix, 1, quantile, probs = 0.025)
  fit2_median   <- apply(fit2_matrix, 1, quantile, probs = 0.5)
  fit2_ci_upper <- apply(fit2_matrix, 1, quantile, probs = 0.975)
  
  fit_time_points <- summary(km_fit, type = "survival", t = seq(0, max(km_fit$time), length.out = 100))$time
  
  # Get the original survival curves
  km_curve <- approx(km_fit$time, km_fit$surv, xout = seq(0, max(km_fit$time), length.out = max(km_fit$time)))$y
  km_time_points <- summary(km_fit, type = "survival", t = seq(0, max(km_fit$time), length.out = max(km_fit$time)))$time
  
  # Create data frames for plotting
  km_data <- data.frame(time = km_time_points, surv = km_curve)
  km_data$surv[1] = 1
  #km_data   <- data.frame(time = fit_time_points, lower = km_ci_lower, upper = km_ci_upper)
  fit1_data <- data.frame(time = fit_time_points, surv = fit1_median, lower = fit1_ci_lower, upper = fit1_ci_upper)
  fit2_data <- data.frame(time = fit_time_points, surv = fit2_median, lower = fit2_ci_lower, upper = fit2_ci_upper)
  # Plot the survival curves with shaded confidence intervals
  return(
    ggplot() + 
      geom_step(data = km_data, aes(x = time, y = surv, color = "Kaplan-Meier", linetype = "Kaplan-Meier"), size = 1) +
      #geom_ribbon(data = km_data, aes(x = time, ymin = lower, ymax = upper), fill = "black", alpha = 0.2) +
      geom_line(data = fit1_data, aes(x = time, y = surv, color = "Traditional", linetype = "Traditional"), size = 1) +
      geom_ribbon(data = fit1_data, aes(x = time, ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
      geom_line(data = fit2_data, aes(x = time, y = surv, color = "Cross Validation", linetype = "Cross Validation"), size = 1) +
      geom_ribbon(data = fit2_data, aes(x = time, ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
      scale_color_manual(values = c("Kaplan-Meier" = "black", "Traditional" = "blue", "Cross Validation" = "red")) +
      scale_linetype_manual(values = c("Kaplan-Meier" = "solid", "Traditional" = "solid", "Cross Validation" = "dashed")) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            legend.position = "right") +
      labs(title = "",
           x = "Time", y = "Survival Probability",
           color = "Group", linetype = "Group")

  )
}

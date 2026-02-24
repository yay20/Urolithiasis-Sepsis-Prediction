# ==============================================================================
# Script 03: Model Evaluation (K-M Curve & DCA Curve)
# ==============================================================================
library(survival)
library(survminer)
library(timeROC)
library(dcurves)
library(ggplot2)
library(cowplot) 
library(dplyr)

# Create output directory for figures
if(!dir.exists("Figures")) dir.create("Figures")

# ------------------------------------------------------------------------------
# 1. Data Loading and Preprocessing
# ------------------------------------------------------------------------------
if(!exists("test_val")) stop("Error: 'test_val' dataset is not loaded.")

# Predict risk scores (Linear Predictor)
test_val$lp <- predict(cox_model, newdata = test_val, type = "lp")

# ------------------------------------------------------------------------------
# 2. Risk Stratification (Determining Cut-off)
# ------------------------------------------------------------------------------
roc_res <- timeROC(
  T = test_val$time,
  delta = test_val$event,
  marker = test_val$lp,
  cause = 1,
  times = 7, 
  weighting = "marginal"
)

# Find optimal cut-off using Youden Index
optimal_idx <- which.max(roc_res$TP[, 1] - roc_res$FP[, 1])
cutoff_val  <- roc_res$cumulative_stats[optimal_idx]

cat("Optimal Risk Score Cutoff (at Day 7):", cutoff_val, "\n")

# Create Risk Group Variable
test_val$risk_group <- ifelse(test_val$lp > cutoff_val, "High Risk", "Low Risk")
test_val$risk_group <- factor(test_val$risk_group, levels = c("Low Risk", "High Risk"))

# ------------------------------------------------------------------------------
# 3. Pre-process Subgroup Variables
# ------------------------------------------------------------------------------
# Ensure variables exist before processing to avoid errors
if("NLR" %in% names(test_val)) {
  nlr_median <- median(test_val$NLR, na.rm = TRUE)
  test_val$nlr_group <- ifelse(test_val$NLR > nlr_median, "High NLR", "Low NLR")
  test_val$nlr_group <- factor(test_val$nlr_group, levels = c("Low NLR", "High NLR"))
}

if("drainage_history" %in% names(test_val)) {
  test_val$drainage_group <- ifelse(test_val$drainage_history == 1, "With Drainage", "No Drainage")
  test_val$drainage_group <- factor(test_val$drainage_group, levels = c("No Drainage", "With Drainage"))
}

if("diabetes" %in% names(test_val)) {
  test_val$diabetes_group <- ifelse(test_val$diabetes == 1, "Diabetes", "No Diabetes")
  test_val$diabetes_group <- factor(test_val$diabetes_group, levels = c("No Diabetes", "Diabetes"))
}

if("hypertension" %in% names(test_val)) {
  test_val$htn_group <- ifelse(test_val$hypertension == 1, "Hypertension", "No Hypertension")
  test_val$htn_group <- factor(test_val$htn_group, levels = c("No Hypertension", "Hypertension"))
}

if("gender" %in% names(test_val)) {
  test_val$gender_group <- ifelse(test_val$gender == 1, "Male", "Female")
  test_val$gender_group <- factor(test_val$gender_group, levels = c("Female", "Male"))
}

# ------------------------------------------------------------------------------
# 4. Universal Function for KM Plots
# ------------------------------------------------------------------------------
plot_km_stratified <- function(data, group_var, title_text, filename) {
  
  # Check if variable exists
  if(!group_var %in% names(data)) {
    warning(paste("Variable", group_var, "not found in data. Skipping plot."))
    return(NULL)
  }
  
  form <- as.formula(paste("Surv(time, event) ~", group_var))
  fit <- survfit(form, data = data)
  
  # Calculate HR
  cox_res <- coxph(form, data = data)
  summ <- summary(cox_res)
  hr <- round(summ$coefficients[2], 2)
  ci_lower <- round(summ$conf.int[3], 2)
  ci_upper <- round(summ$conf.int[4], 2)
  pval <- round(summ$coefficients[5], 4)
  pval_txt <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", pval))
  
  # Plot
  p <- ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = FALSE,
    risk.table = TRUE, 
    risk.table.col = "strata",
    risk.table.height = 0.25,
    palette = c("#2E9FDF", "#E74C3C"),
    ggtheme = theme_classic(),
    xlab = "Time (days)",
    ylab = "Sepsis-free Survival Probability",
    title = title_text,
    xlim = c(0, 14),
    break.time.by = 2,
    legend.labs = levels(data[[group_var]])
  )
  
  # Add HR Annotation
  p$plot <- p$plot + 
    annotate("text", x = 4, y = 0.15, 
             label = paste0("HR = ", hr, " (95% CI: ", ci_lower, "-", ci_upper, ")\n", pval_txt),
             size = 4, fontface = "bold", color = "black", hjust = 0)
  
  # Save to Figures folder
  pdf_name <- paste0("Figures/KM_", filename, ".pdf")
  
  final_plot <- arrange_ggsurvplots(list(p), print = FALSE, ncol = 1, nrow = 1)
  ggsave(pdf_name, final_plot, width = 7, height = 7)
  
  message(paste("Saved KM plot:", pdf_name))
}

# ------------------------------------------------------------------------------
# 5. Execute KM Loops
# ------------------------------------------------------------------------------
km_tasks <- list(
  list(var = "risk_group",     title = "A. Risk Stratification", file = "RiskGroup"),
  list(var = "diabetes_group", title = "B. Diabetes Status",     file = "Diabetes"),
  list(var = "gender_group",   title = "C. Gender",              file = "Gender"),
  list(var = "drainage_group", title = "D. Drainage History",    file = "Drainage"),
  list(var = "htn_group",      title = "E. Hypertension",        file = "Hypertension"),
  list(var = "nlr_group",      title = "F. NLR Stratification",  file = "NLR")
)

for (task in km_tasks) {
  plot_km_stratified(test_val, task$var, task$title, task$file)
}

# ------------------------------------------------------------------------------
# 6. Decision Curve Analysis (DCA) Loops
# ------------------------------------------------------------------------------

# --- 6.1 Data Preparation ---

dca_data_raw <- data.frame(
  time = cox_model$y[, 1],
  status = cox_model$y[, 2],
  lp = cox_model$linear.predictors
)
lp_centered <- dca_data_raw$lp - mean(dca_data_raw$lp)

# Function to extract Baseline Survival (S0) from model object
get_S0_at_time <- function(model_obj, target_day) {
  times <- model_obj$y[, 1]; status <- model_obj$y[, 2]; lp <- model_obj$linear.predictors
  df <- data.frame(time = times, status = status, risk_score = exp(lp))
  df <- df[order(df$time), ]
  unique_times <- sort(unique(df$time[df$status == 1]))
  baseline_hazard <- numeric(length(unique_times)); cum_haz <- 0
  for(i in 1:length(unique_times)) {
    t_current <- unique_times[i]
    risk_set_sum <- sum(df$risk_score[df$time >= t_current])
    if(risk_set_sum > 0) cum_haz <- cum_haz + (sum(df$time == t_current & df$status == 1) / risk_set_sum)
    baseline_hazard[i] <- cum_haz
  }
  idx <- which(unique_times <= target_day)
  if(length(idx) == 0) return(1.0)
  return(exp(-baseline_hazard[max(idx)]))
}

# --- 6.2 DCA Plotting Function ---
plot_dca_custom <- function(data, time_point, ylim_max, filename) {
  
  # Calculate Risk for specific time point using 'cox_model'
  s0 <- get_S0_at_time(cox_model, time_point)
  
  risk_col_name <- paste0("risk_", time_point, "d")
  data[[risk_col_name]] <- 1 - s0 ^ exp(lp_centered)
  
  dca_form <- as.formula(paste("Surv(time, status) ~", risk_col_name))
  
  dca_res <- dca(
    dca_form,
    data = data,
    time = time_point,
    thresholds = seq(0, 0.50, by = 0.01),
    label = setNames(list("Dynamic RCS Model"), risk_col_name)
  )
  
  p <- plot(dca_res, smooth = TRUE, span = 0.3) +
    coord_cartesian(ylim = c(-0.002, ylim_max)) +
    labs(title = paste0("DCA (", time_point, "-Day Prediction)")) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = c(0.8, 0.8),
      panel.grid.minor = element_blank()
    )
  
  # Save to Figures folder
  ggsave(paste0("Figures/DCA_", filename, ".pdf"), p, width = 5, height = 5)
  
  message(paste("Saved DCA plot for Day", time_point))
}

# --- 6.3 Execute DCA Loops ---
dca_tasks <- list(
  list(day = 1,  ylim = 0.005, file = "Day1"), 
  list(day = 3,  ylim = 0.005, file = "Day3"), 
  list(day = 7,  ylim = 0.015, file = "Day7"), 
  list(day = 14, ylim = 0.015, file = "Day14") 
)

for (task in dca_tasks) {
  plot_dca_custom(dca_data_raw, task$day, task$ylim, task$file)
}

message("All evaluation plots generated successfully in 'Figures' folder!")
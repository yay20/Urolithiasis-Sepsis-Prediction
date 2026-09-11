# ==============================================================================
# Script 03: Model Evaluation (K-M Survival Analysis & Decision Curve Analysis)
# ==============================================================================

library(survival)
library(survminer)
library(timeROC)
library(dcurves)
library(ggplot2)
library(dplyr)

if (!dir.exists("Figures")) dir.create("Figures")

# ------------------------------------------------------------------------------
# 1. Data Ingestion & Linear Predictors
# ------------------------------------------------------------------------------
cox_model  <- readRDS("cox_model.rds")
train_data <- read.csv("train_data.csv")
test_data  <- read.csv("test_data.csv")

eval_horizons <- c(1, 3, 7, 14)

train_data$lp <- predict(cox_model, newdata = train_data, type = "lp")
test_data$lp  <- predict(cox_model, newdata = test_data,  type = "lp")

# ------------------------------------------------------------------------------
# 2. Cut-off Determination & Risk Stratification
# ------------------------------------------------------------------------------
roc_train <- timeROC(
  T = train_data$time, delta = train_data$event, marker = train_data$lp,
  cause = 1, times = 7, weighting = "marginal", ROC = FALSE, iid = FALSE
)

opt_idx    <- which.max(roc_train$TP[, 1] - roc_train$FP[, 1])
cutoff_val <- roc_train$cumulative_stats[opt_idx]

train_data$risk_group <- factor(ifelse(train_data$lp > cutoff_val, "High risk", "Low risk"), levels = c("Low risk", "High risk"))
test_data$risk_group  <- factor(ifelse(test_data$lp > cutoff_val, "High risk", "Low risk"), levels = c("Low risk", "High risk"))

# ------------------------------------------------------------------------------
# 3. Standard Kaplan-Meier Plotting Engine
# ------------------------------------------------------------------------------
plot_km_pipeline <- function(df, strata_var, title_label, filename_prefix) {
  fit_obj <- survfit(as.formula(paste("Surv(time, event) ~", strata_var)), data = df)
  
  km_plot <- ggsurvplot(
    fit_obj,
    data           = df,
    palette        = c("#2E9FDF", "#E74C3C"),
    censor         = TRUE,
    risk.table     = TRUE,
    risk.table.col = "strata",
    pval           = TRUE,
    pval.method    = TRUE,
    conf.int       = FALSE,
    title          = title_label,
    legend.title   = "Strata",
    legend.labs    = levels(df[[strata_var]]),
    xlab           = "Time (days)",
    ylab           = "Sepsis-free Survival Probability",
    ggtheme        = theme_classic()
  )
  
  export_obj <- arrange_ggsurvplots(list(km_plot), print = FALSE, ncol = 1, nrow = 1)
  ggsave(paste0("Figures/", filename_prefix, ".png"), export_obj, width = 7.0, height = 6.5, dpi = 300)
  ggsave(paste0("Figures/", filename_prefix, ".pdf"), export_obj, width = 7.0, height = 6.5)
}

# ------------------------------------------------------------------------------
# 4. Multivariable & Univariable K-M Evaluation
# ------------------------------------------------------------------------------
# 4.1 Primary Model-Stratified Analysis 
plot_km_pipeline(test_data, "risk_group", "Model Risk Stratification", "KM_Model_Stratification")

# 4.2 Univariable Categorical Predictor Analysis 
df_eval <- test_data

# Automatically select only categorical and binary features
cat_vars <- names(df_eval)[sapply(df_eval, function(x) is.factor(x) || length(unique(na.omit(x))) == 2)]
cat_vars <- setdiff(cat_vars, c("event", "risk_group"))

for (target in cat_vars) {
  plot_km_pipeline(df_eval, target, paste("Univariable Stratification:", target), paste0("KM_Univariable_", target))
}

# ------------------------------------------------------------------------------
# 5. Decision Curve Analysis (DCA) Across Horizons 
# ------------------------------------------------------------------------------
calc_s0 <- function(model_obj, t) {
  bh  <- basehaz(model_obj, centered = TRUE)
  idx <- which(bh$time <= t)
  if (length(idx) == 0) return(1.0)
  exp(-bh$hazard[max(idx)])
}

for (h in eval_horizons) {
  s0_h <- calc_s0(cox_model, h)
  risk_col <- paste0("risk_", h, "d")
  test_data[[risk_col]] <- 1 - s0_h^exp(test_data$lp - mean(train_data$lp))
  
  dca_res <- dca(
    as.formula(paste("Surv(time, event) ~", risk_col)),
    data       = test_data,
    time       = h,
    thresholds = seq(0, 0.50, by = 0.01),
    label      = setNames(list("Dynamic RCS Model"), risk_col)
  )
  
  p_dca <- plot(dca_res, smooth = TRUE, span = 0.3) +
    labs(title = paste0("DCA (", h, "-Day Prediction)")) +
    theme_bw() +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(paste0("Figures/DCA_Horizon_", h, "d.png"), p_dca, width = 6.0, height = 5.5, dpi = 300)
  ggsave(paste0("Figures/DCA_Horizon_", h, "d.pdf"), p_dca, width = 6.0, height = 5.5)
}
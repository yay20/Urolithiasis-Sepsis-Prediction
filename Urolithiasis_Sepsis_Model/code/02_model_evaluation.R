# ==============================================================================
# Script 02: Model Evaluation
# ==============================================================================

library(survival)
library(timeROC)
library(pROC)
library(riskRegression)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
# 1. Data Preparation
# ------------------------------------------------------------------------------
cox_model  <- readRDS("cox_model.rds")
train_data <- read.csv("train_data.csv")
val_data   <- read.csv("validation_data.csv")
test_data  <- read.csv("test_data.csv")

time_points <- c(1, 3, 7, 14)

categorical_vars <- names(test_data)[sapply(test_data, function(x) is.factor(x) || is.character(x) || length(unique(na.omit(x))) <= 2)]
categorical_vars <- setdiff(categorical_vars, c("event", "time", "pred_risk"))

train_data <- train_data %>% mutate(across(all_of(categorical_vars), as.factor))
val_data   <- val_data   %>% mutate(across(all_of(categorical_vars), as.factor))
test_data  <- test_data  %>% mutate(across(all_of(categorical_vars), as.factor))

train_data$pred_risk <- predict(cox_model, newdata = train_data, type = "risk")
val_data$pred_risk   <- predict(cox_model, newdata = val_data,   type = "risk")
test_data$pred_risk  <- predict(cox_model, newdata = test_data,  type = "risk")

cohorts_list <- list(Training = train_data, Validation = val_data, Testing = test_data)

# ------------------------------------------------------------------------------
# 2. Time-Dependent AUC Evaluation
# ------------------------------------------------------------------------------
calc_auc_ci <- function(df, times, cohort_name) {
  roc_res <- timeROC(T = df$time, delta = df$event, marker = df$pred_risk,
                     cause = 1, times = times, weighting = "marginal",
                     ROC = FALSE, iid = FALSE)
  aucs <- roc_res$AUC
  n1   <- sum(df$event == 1)
  n0   <- sum(df$event == 0)
  q1   <- aucs / (2 - aucs)
  q2   <- (2 * aucs^2) / (1 + aucs)
  se   <- sqrt((aucs * (1 - aucs) + (n1 - 1) * (q1 - aucs^2) + (n0 - 1) * (q2 - aucs^2)) / (n1 * n0))
  
  ci_low <- pmax(0, aucs - 1.96 * se)
  ci_upp <- pmin(1, aucs + 1.96 * se)
  
  data.frame(
    Time      = times,
    Cohort    = cohort_name,
    AUC       = aucs,
    Formatted = sprintf("%.6f(%.6f-%.6f)", aucs, ci_low, ci_upp)
  )
}

auc_records <- do.call(rbind, lapply(names(cohorts_list), function(c_name) {
  calc_auc_ci(cohorts_list[[c_name]], time_points, c_name)
}))

# Export summary table
auc_summary_table <- reshape(auc_records[, c("Time", "Cohort", "Formatted")],
                             idvar = "Time", timevar = "Cohort", direction = "wide")
colnames(auc_summary_table) <- gsub("Formatted.", "", colnames(auc_summary_table))
write.csv(auc_summary_table, "time_dependent_auc_summary.csv", row.names = FALSE)

# Plot time-dependent AUC curves
auc_records_fig3 <- auc_records %>%
  filter(Cohort %in% c("Validation", "Testing")) %>%
  mutate(Cohort = ifelse(Cohort == "Testing", "Test", Cohort))

p_auc <- ggplot(auc_records_fig3, aes(x = Time, y = AUC, color = Cohort, group = Cohort)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(1, 3, 7, 14), labels = c("24h", "72h", "168h", "336h")) +
  scale_y_continuous(limits = c(0.60, 0.90), breaks = seq(0.60, 0.90, 0.05)) +
  scale_color_manual(values = c("Test" = "#1E88E5", "Validation" = "#D81B60")) +
  theme_classic(base_size = 11) +
  labs(
    title = "Time-Dependent AUC Performance",
    x     = "Follow-up Time",
    y     = "Time-Dependent AUC",
    color = NULL
  ) +
  theme(
    plot.title           = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position      = "top",
    legend.justification = "center",
    panel.grid.major     = element_line(color = "gray90", linewidth = 0.3)
  )

ggsave("time_dependent_auc_curves.png", p_auc, width = 6.5, height = 4.5, dpi = 300)
ggsave("time_dependent_auc_curves.pdf", p_auc, width = 6.5, height = 4.5)

# ------------------------------------------------------------------------------
# 3. Model Calibration Assessment
# ------------------------------------------------------------------------------
calib_score <- Score(list("Model" = cox_model), formula = Surv(time, event) ~ 1,
                     data = test_data, times = time_points, plots = "calibration",
                     metrics = NULL, summary = NULL)

all_calib_df <- calib_score$Calibration$plotframe
all_calib_df$Horizon <- factor(paste0(all_calib_df$times, "-Day"),
                               levels = paste0(time_points, "-Day"))

p_calib <- ggplot(all_calib_df, aes(x = pred, y = Obs)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "loess", se = FALSE, color = "#1976D2", linewidth = 0.9) +
  facet_wrap(~ Horizon, ncol = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 11) +
  labs(x = "Predicted Risk Probability", y = "Observed Event Rate") +
  theme(
    strip.background = element_rect(fill = "#F0F0F0", color = "black"),
    strip.text       = element_text(face = "bold"),
    aspect.ratio     = 1
  )

ggsave("calibration_curves_multihorizon.png", p_calib, width = 7, height = 6.5, dpi = 300)
ggsave("calibration_curves_multihorizon.pdf", p_calib, width = 7, height = 6.5)

# ------------------------------------------------------------------------------
# 4. Classification Performance and Precision-Recall Analysis
# ------------------------------------------------------------------------------
pr_metrics_list <- list()
pr_curves_list  <- list()

for (c_name in names(cohorts_list)) {
  df <- cohorts_list[[c_name]]
  for (t in time_points) {
    valid <- !(df$time < t & df$event == 0)
    sub_y <- ifelse(df$time[valid] <= t & df$event[valid] == 1, 1, 0)
    sub_p <- df$pred_risk[valid]
    
    roc_obj  <- roc(sub_y, sub_p, quiet = TRUE)
    cutoff   <- coords(roc_obj, "best", best.method = "youden", transpose = FALSE)$threshold[1]
    pred_bin <- ifelse(sub_p >= cutoff, 1, 0)
    
    tp <- sum(pred_bin == 1 & sub_y == 1)
    tn <- sum(pred_bin == 0 & sub_y == 0)
    fp <- sum(pred_bin == 1 & sub_y == 0)
    fn <- sum(pred_bin == 0 & sub_y == 1)
    
    sen <- tp / (tp + fn)
    spe <- tn / (tn + fp)
    ppv <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
    npv <- ifelse(tn + fn > 0, tn / (tn + fn), 0)
    
    ord       <- order(sub_p, decreasing = TRUE)
    y_ord     <- sub_y[ord]
    tp_seq    <- cumsum(y_ord == 1)
    fp_seq    <- cumsum(y_ord == 0)
    rec_vec   <- c(0, tp_seq / sum(sub_y == 1))
    prec_vec  <- c(1, tp_seq / (tp_seq + fp_seq))
    pr_auc    <- sum(diff(rec_vec) * (prec_vec[-1] + prec_vec[-length(prec_vec)]) / 2)
    
    time_label <- ifelse(t == 1, "1 day", paste0(t, " days"))
    pr_metrics_list[[paste(c_name, t)]] <- data.frame(
      Time_Point  = time_label,
      Cohort      = c_name,
      Sensitivity = sprintf("%.3f", sen),
      Specificity = sprintf("%.3f", spe),
      PPV         = sprintf("%.3f", ppv),
      NPV         = sprintf("%.3f", npv),
      PR_AUC      = sprintf("%.3f", pr_auc)
    )
    
    pr_curves_list[[paste(c_name, t)]] <- data.frame(
      Recall    = rec_vec,
      Precision = prec_vec,
      Time      = paste0(t, "Days"),
      Cohort    = paste0(c_name, " Set")
    )
  }
}

classification_table <- do.call(rbind, pr_metrics_list)
write.csv(classification_table, "classification_precision_recall_metrics.csv", row.names = FALSE)

all_pr_df <- do.call(rbind, pr_curves_list)
all_pr_df$Time   <- factor(all_pr_df$Time, levels = paste0(time_points, "Days"))
all_pr_df$Cohort <- factor(all_pr_df$Cohort, levels = c("Training Set", "Validation Set", "Testing Set"))

p_pr <- ggplot(all_pr_df, aes(x = Recall, y = Precision, color = Time)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ Cohort, nrow = 1) +
  scale_color_manual(values = c("#2E7D32", "#1976D2", "#E65100", "#C2185B")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 11) +
  labs(x = "Recall", y = "Precision", color = "Time Points") +
  theme(
    strip.background = element_rect(fill = "#F5F5F5", color = "black"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "top",
    aspect.ratio     = 1
  )

ggsave("precision_recall_curves_multicohort.png", p_pr, width = 11, height = 4.2, dpi = 300)
ggsave("precision_recall_curves_multicohort.pdf", p_pr, width = 11, height = 4.2)

# ------------------------------------------------------------------------------
# 5. Missingness Sensitivity Analysis
# ------------------------------------------------------------------------------
cont_features <- setdiff(names(test_data), c(categorical_vars, "time", "event", "pred_risk"))

format_sensitivity_row <- function(df, label) {
  res <- calc_auc_ci(df, time_points, label)
  data.frame(
    Strategy = label,
    AUC_1d   = res$Formatted[res$Time == 1],
    AUC_3d   = res$Formatted[res$Time == 3],
    AUC_7d   = res$Formatted[res$Time == 7],
    AUC_14d  = res$Formatted[res$Time == 14]
  )
}

res_base <- format_sensitivity_row(test_data, "Primary(missForest)")

test_10 <- test_data
set.seed(123)
for (col in cont_features) {
  if (is.numeric(test_10[[col]])) {
    mask <- sample(nrow(test_10), size = floor(0.10 * nrow(test_10)))
    test_10[[col]][mask] <- median(test_10[[col]][-mask], na.rm = TRUE)
  }
}
test_10$pred_risk <- predict(cox_model, newdata = test_10, type = "risk")
res_10 <- format_sensitivity_row(test_10, "10% Missing (Median Imputed)")

test_20 <- test_data
set.seed(456)
for (col in cont_features) {
  if (is.numeric(test_20[[col]])) {
    mask <- sample(nrow(test_20), size = floor(0.20 * nrow(test_20)))
    test_20[[col]][mask] <- median(test_20[[col]][-mask], na.rm = TRUE)
  }
}
test_20$pred_risk <- predict(cox_model, newdata = test_20, type = "risk")
res_20 <- format_sensitivity_row(test_20, "20% Missing (Median Imputed)")

sensitivity_table <- rbind(res_base, res_10, res_20)
write.csv(sensitivity_table, "missingness_sensitivity_analysis.csv", row.names = FALSE)
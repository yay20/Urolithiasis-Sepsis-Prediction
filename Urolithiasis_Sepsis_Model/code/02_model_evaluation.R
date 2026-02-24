# ==============================================================================
# Script 02: Model Evaluation (Calibration & Time-dependent AUC)
# ==============================================================================

library(survival)
library(timeROC)
library(ggplot2)
library(dplyr)
library(riskRegression) # Excellent for calibration plots




# Load necessary data


# ------------------------------------------------------------------------------
# 1. Calibration Curves (Loop for 1, 3, 7, 14 Days)
# ------------------------------------------------------------------------------
# Define time points of interest
time_points <- c(1, 3, 7, 14)

plot_calibration_loop <- function(model, data, times) {
  
  for (t in times) {
    cat("Generating Calibration Curve for Day", t, "...\n")
    calib_res <- Score(list("Model" = model),
                       formula = Surv(time, event) ~ 1,
                       data = data,
                       times = t,
                       plots = "calibration",
                       metrics = NULL, 
                       summary = NULL)
    
    # Extract plotting data
    plot_df <- calib_res$Calibration$plotframe
    
    # Plot using ggplot2
    p <- ggplot(plot_df, aes(x = pred, y = Obs)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      geom_smooth(method = "loess", se = FALSE, color = "#2E9FDF", size = 1.2, span = 0.8) +
      # geom_point(alpha=0.5) + # Optional: show bins
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), name = "Predicted Probability") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), name = "Observed Probability") +
      labs(title = paste0("Calibration Curve (", t, "-Day)")) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 12)
      ) +
      coord_fixed()
    
    # Save plot
    ggsave(paste0("Calibration_Day", t, ".pdf"), p, width = 6, height = 6)
    print(p)
  }
}


# ------------------------------------------------------------------------------
# 2. Time-Dependent AUC (Combined Plot for Val/Test)
# ------------------------------------------------------------------------------

# Calculate Time-Dependent AUC
calc_tdAUC <- function(data, model, times) {
  lp <- predict(model, newdata = data, type = "lp")
  
  res <- timeROC(
    T = data$time,
    delta = data$event,
    marker = lp,
    cause = 1,
    times = times,
    weighting = "marginal",
    iid = TRUE 
  )
  return(res)
}

# Plot AUC for Validation and Test sets
plot_val_test_auc <- function(auc_val, auc_test, time_points) {
  
  # Prepare data
  plot_data <- rbind(
    data.frame(Time = time_points, AUC = auc_val$AUC, Dataset = "Validation"),
    data.frame(Time = time_points, AUC = auc_test$AUC, Dataset = "Test")
  )
  
  # Calculate mean AUC for horizontal lines
  mean_aucs <- plot_data %>% 
    group_by(Dataset) %>% 
    summarise(MeanAUC = mean(AUC))
  
  # Define colors: Test (Blue), Validation (Red)
  my_colors <- c("Test" = "#2196F3", "Validation" = "#E74C3C")
  
  p <- ggplot(plot_data, aes(x = Time, y = AUC, color = Dataset, group = Dataset)) +
    
    # Vertical lines at time points
    geom_vline(xintercept = time_points, color = "grey80", linewidth = 1) +
    
    # Horizontal mean lines
    geom_hline(data = mean_aucs, aes(yintercept = MeanAUC, color = Dataset), 
               linetype = "dashed", linewidth = 0.8, alpha = 0.8) +
    
    # Main AUC curves
    geom_line(linewidth = 1.2) +
    
    # Axis formatting (Days to Hours)
    scale_x_continuous(breaks = time_points, 
                       labels = paste0(time_points * 24, "h")) +
    scale_y_continuous(limits = c(0.60, 0.95), breaks = seq(0.60, 0.95, 0.05)) +
    
    scale_color_manual(values = my_colors) +
    
    labs(title = "Time-Dependent AUC Performance",
         x = "Follow-up Time",
         y = "Time-Dependent AUC",
         color = NULL) +
    
    # Theme adjustments
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "top",
      legend.text = element_text(size = 10),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      panel.grid.major = element_line(color = "grey95"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "grey80")
    )
  
  return(p)
}




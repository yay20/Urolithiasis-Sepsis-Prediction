# ==============================================================================
  # Script 01: Variable Selection and Model Construction
  # ==============================================================================
  
  library(survival)
  library(glmnet)
  library(rms)
  library(ggplot2)
  library(dplyr)
  
  # ------------------------------------------------------------------------------
  # 1. Data Preparation
  # ------------------------------------------------------------------------------
  categorical_vars <- names(df)[sapply(df, function(x) is.factor(x) || is.character(x) || length(unique(na.omit(x))) <= 2)]
  categorical_vars <- setdiff(categorical_vars, c("event", "time"))
  
  df <- df %>%
    mutate(across(all_of(categorical_vars), as.factor))
  
  df <- df[df$time > 0, ]
  
  # ------------------------------------------------------------------------------
  # 2. Variable Selection via Regularized Elastic Net
  # ------------------------------------------------------------------------------
  predictor_vars <- setdiff(names(df), c("time", "event"))
  x <- data.matrix(df[, predictor_vars])
  y <- Surv(df$time, df$event)
  
  set.seed(123)
  cv_enet <- cv.glmnet(
    x, y,
    family      = "cox",
    alpha       = 0.5,
    nfolds      = 10,
    standardize = TRUE
  )
  
  coef_enet     <- coef(cv_enet, s = "lambda.1se")
  selected_vars <- rownames(coef_enet)[which(as.numeric(coef_enet) != 0)]
  
  write.csv(data.frame(Predictor = selected_vars), "selected_predictors.csv", row.names = FALSE)
  
  # ------------------------------------------------------------------------------
  # 3. Model Formula Construction with Restricted Cubic Splines
  # ------------------------------------------------------------------------------
  selected_categorical <- intersect(selected_vars, categorical_vars)
  selected_continuous  <- setdiff(selected_vars, categorical_vars)
  
  formula_terms <- c(
    if (length(selected_continuous) > 0) paste0("rcs(", selected_continuous, ", 5)") else NULL,
    selected_categorical
  )
  
  cox_formula <- as.formula(paste("Surv(time, event) ~", paste(formula_terms, collapse = " + ")))
  
  # ------------------------------------------------------------------------------
  # 4. Multivariable Cox Model Fitting
  # ------------------------------------------------------------------------------
  dd <- datadist(df); options(datadist = "dd")
  cox_model <- cph(cox_formula, data = df, x = TRUE, y = TRUE, surv = TRUE)
  
  saveRDS(cox_model, "cox_model.rds")
  
  # ------------------------------------------------------------------------------
  # 5. Predictor Risk Quantification Across Quartiles
  # ------------------------------------------------------------------------------
  term_labels <- attr(cox_model$terms, "term.labels")
  pvars       <- as.list(attr(cox_model$terms, "predvars"))
  
  calc_rcs_basis <- function(x, k) {
    t1 <- k[1]; t4 <- k[4]; t5 <- k[5]
    cube_p <- function(u) pmax(0, u)^3
    b <- matrix(NA, nrow = length(x), ncol = 4)
    b[, 1] <- x
    b[, 2] <- (cube_p(x - k[1]) - cube_p(x - t4)*(t5 - k[1])/(t5 - t4) + cube_p(x - t5)*(t4 - k[1])/(t5 - t4)) / (t5 - t1)^2
    b[, 3] <- (cube_p(x - k[2]) - cube_p(x - t4)*(t5 - k[2])/(t5 - t4) + cube_p(x - t5)*(t4 - k[2])/(t5 - t4)) / (t5 - t1)^2
    b[, 4] <- (cube_p(x - k[3]) - cube_p(x - t4)*(t5 - k[3])/(t5 - t4) + cube_p(x - t5)*(t4 - k[3])/(t5 - t4)) / (t5 - t1)^2
    return(b)
  }
  
  table_records <- list()
  curve_records <- list()
  
  for (i in seq_along(term_labels)) {
    term     <- term_labels[i]
    coef_idx <- cox_model$assign[[term]]
    betas    <- as.numeric(cox_model$coefficients[coef_idx])
    sub_var  <- as.matrix(cox_model$var[coef_idx, coef_idx])
    
    if (grepl("^rcs\\(", term)) {
      var_name <- gsub("^rcs\\(([^,]+),.*$", "\\1", term)
      
      knots <- NULL
      for (p in pvars) {
        call_str <- paste(deparse(p), collapse = "")
        if (grepl(paste0("^rcs\\(", var_name), call_str)) {
          knots <- eval(p[[3]])
          break
        }
      }
      
      x_seq      <- seq(knots[1], knots[5], length.out = 200)
      basis      <- calc_rcs_basis(x_seq, knots)
      ref_basis  <- calc_rcs_basis(knots[3], knots)
      basis_diff <- sweep(basis, 2, ref_basis, "-")
      
      log_hr  <- basis_diff %*% betas
      se_log  <- sqrt(rowSums((basis_diff %*% sub_var) * basis_diff))
      hr_raw  <- exp(log_hr)
      
      curve_records[[var_name]] <- data.frame(
        Value    = x_seq,
        HR       = hr_raw,
        Lower    = exp(log_hr - 1.96 * se_log),
        Upper    = pmin(exp(log_hr + 1.96 * se_log), pmax(hr_raw * 3, 10)),
        Variable = var_name
      )
      
      q25_val <- knots[1] + (knots[2] - knots[1]) * (25 - 5) / (27.5 - 5)
      q75_val <- knots[4] + (knots[5] - knots[4]) * (75 - 72.5) / (95 - 72.5)
      pts     <- c(q25_val, knots[3], q75_val)
      
      pts_diff <- sweep(calc_rcs_basis(pts, knots), 2, ref_basis, "-")
      pts_log  <- pts_diff %*% betas
      pts_se   <- sqrt(rowSums((pts_diff %*% sub_var) * pts_diff))
      
      wald_stat <- as.numeric(t(betas) %*% solve(sub_var) %*% betas)
      p_val     <- 1 - pchisq(wald_stat, df = length(betas))
      p_str     <- ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val))
      
      table_records[[var_name]] <- data.frame(
        Variable = var_name,
        Quartile = c("25th", "50th", "75th"),
        Value    = sprintf("%.2f", pts),
        HR_95CI  = sprintf("%.3f (%.3f-%.3f)", exp(pts_log), exp(pts_log - 1.96 * pts_se), exp(pts_log + 1.96 * pts_se)),
        P_value  = c(p_str, "", ""),
        stringsAsFactors = FALSE
      )
    } else {
      var_name <- term
      b        <- betas[1]
      s        <- sqrt(sub_var[1, 1])
      p_val    <- 1 - pchisq((b^2) / (s^2), df = 1)
      p_str    <- ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val))
      
      table_records[[var_name]] <- data.frame(
        Variable = var_name,
        Quartile = "-",
        Value    = "-",
        HR_95CI  = sprintf("%.3f (%.3f-%.3f)", exp(b), exp(b - 1.96 * s), exp(b + 1.96 * s)),
        P_value  = p_str,
        stringsAsFactors = FALSE
      )
    }
  }
  
  hazard_ratios_df <- do.call(rbind, table_records)
  write.csv(hazard_ratios_df, "hazard_ratios_quartiles.csv", row.names = FALSE)
  
  # ------------------------------------------------------------------------------
  # 6. Non-linear Risk Trajectory Visualization
  # ------------------------------------------------------------------------------
  if (length(curve_records) > 0) {
    curves_df <- do.call(rbind, curve_records)
    
    p_trajectories <- ggplot(curves_df, aes(x = Value, y = HR)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#1976D2", alpha = 0.2) +
      geom_line(color = "#D32F2F", linewidth = 0.75) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey40") +
      facet_wrap(~ Variable, scales = "free", ncol = 5) +
      theme_bw(base_size = 10) +
      labs(x = "Biomarker Value", y = "Hazard Ratio (95% CI)") +
      theme(
        strip.background = element_rect(fill = "#F0F0F0", color = "black"),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        aspect.ratio     = 0.8
      )
    
    ggsave("continuous_risk_trajectories.png", p_trajectories, width = 14, height = 9, dpi = 300)
    ggsave("continuous_risk_trajectories.pdf", p_trajectories, width = 14, height = 9)
  }

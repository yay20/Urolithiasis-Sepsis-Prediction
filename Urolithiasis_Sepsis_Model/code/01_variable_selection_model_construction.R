# ==============================================================================
# Script 01: Variable Selection and Model Construction
# ==============================================================================

library(survival)
library(glmnet)
library(dplyr)
library(rms)

# ------------------------------------------------------------------------------
# 1. Data Preparation
# ------------------------------------------------------------------------------



# Convert categorical factors to numeric for glmnet input
df <- df %>%
  mutate(across(c(gender, hypertension, drainage_history, diabetes), 
                ~as.numeric(as.character(.x))))

# Calculate survival time
# Ensure time is positive
df$time <- df$time_stop - df$time_start
df <- df[df$time > 0, ]



# ------------------------------------------------------------------------------
# 2. Variable Selection via Elastic Net
# ------------------------------------------------------------------------------
# Prepare matrix for glmnet
predictor_vars <- names(df)[!names(df) %in% c("time", "event")]
x <- as.matrix(df[, predictor_vars])
y <- Surv(df$time, df$event)

set.seed(123)
# Alpha = 0.5 balances Lasso (L1) and Ridge (L2) penalties
cv_enet <- cv.glmnet(x, y,
                     family = "cox",
                     alpha = 0.5,
                     nfolds = 10,
                     standardize = TRUE)

# Extract non-zero coefficients at lambda.1se (more regularized model)
coef_enet <- coef(cv_enet, s = "lambda.1se")
selected_vars_enet <- rownames(coef_enet)[which(coef_enet != 0)]

cat("Selected Variables by Elastic Net:\n", selected_vars_enet, "\n")

# Save selection results
write.csv(selected_vars_enet, "selected_variables.csv", row.names = FALSE)


# 3. Optimize RCS Knots for Continuous Variables
# ------------------------------------------------------------------------------
formula_terms <- sapply(selected_vars_enet, function(var) {
  if (exists("continuous_vars_std") && var %in% continuous_vars_std) {
    return(paste0("rcs(", var, ", 5)"))
  } else {
    return(var)
  }
})

cox_formula <- as.formula(paste("Surv(time, event) ~", paste(formula_terms, collapse = " + ")))

# 4. Fit the Cox Proportional Hazards Model
dd <- datadist(df); options(datadist = "dd")
cox_model <- cph(cox_formula, data = df, x = TRUE, y = TRUE, surv = TRUE)






library(survival)
library(rms)

# 1. Load Model
risk_engine <- readRDS("clean_model_for_public.rds")
model <- risk_engine$model_obj
std_params <- risk_engine$std_params
base_s0 <- risk_engine$base_s0

# 2. Create Dummy Patient (Example)
new_patient <- data.frame(
  age=65, gender=1, diabetes=1, hypertension=1, drainage_history=0,
  wbc_value=15.0, rbc_value=3.5, hct_value=35.0, mcv_value=85.0, plt_value=80.0,
  gran1_value=12.0, gran2_value=85.0, lymph1_value=0.8, lymph2_value=10.0,
  alb_value=28.0, cr_value=150.0, ua_value=400.0, uph_value=6.0, k_value=4.5,
  dbil_value=10.0, fib_value=5.0
)
new_patient$NLR <- new_patient$gran1_value / new_patient$lymph1_value

# 3. Standardization & Clamping
std_patient <- new_patient
for(var in names(std_params)) {
  if(var %in% names(std_patient)) {
    p <- std_params[[var]]
    z <- (new_patient[[var]] - p[1]) / p[2]
    std_patient[[var]] <- max(min(z, 3.0), -3.0) # Clamp [-3, 3]
  }
}

# 4. Type Conversion
for (n in names(std_patient)) {
  if (n %in% names(model$xlevels)) {
    std_patient[[n]] <- factor(std_patient[[n]], levels = model$xlevels[[n]])
  }
}

# 5. Predict
lp_patient <- predict(model, newdata = std_patient, type = "lp")

# Calculate Reference LP (Z=0)
ref_df <- std_patient
for(var in names(std_params)) ref_df[[var]] <- 0
lp_ref <- predict(model, newdata = ref_df, type = "lp")

# Calculate Risk
hr <- exp(lp_patient - lp_ref)
risk_14d <- 1 - base_s0^hr

cat(sprintf("Predicted 14-Day Sepsis Risk: %.1f%%\n", risk_14d * 100))


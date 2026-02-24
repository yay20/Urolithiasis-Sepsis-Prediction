
# Dynamic Prediction Model for Sepsis in Urolithiasis Patients

This repository contains the pre-trained **Cox Proportional Hazards Model** and prediction scripts associated with the study:

> **Dynamic Prediction Model for Sepsis in Urolithiasis Using Routine Clinical Parameters: A Large-scale Multi-center Cohort Study**

## 1. Overview

The model was developed using a multi-center cohort of **>100,000 patients**. It utilizes **Restricted Cubic Splines (RCS)** to capture non-linear relationships between **22 routine clinical parameters** (e.g., Blood Routine, Biochemistry) and the risk of sepsis.

## 2. Repository Contents

*   `risk_model.rds`: The R object containing the trained model coefficients, RCS knot locations, and standardization parameters.
    *   *Note: Patient-level data has been strictly removed to protect privacy.*
*   `run_prediction.R`: An R script demonstrating how to load the model and calculate risk for a new patient.

## 3. Requirements

*   R (>= 4.0.0)
*   R packages: `survival`, `rms`

## 4. How to Use

1.  Download `risk_model.rds` and `run_prediction.R` to your local machine.
2.  Open `run_prediction.R` in RStudio.
3.  Install necessary packages if missing: `install.packages(c('survival', 'rms'))`.
4.  Run the script to see a demo prediction.
5.  Replace the dummy data in the script with your own local data to test the model.

## 5. Privacy & Data Availability Statement

The original training dataset (containing ~100,000 individual records) is **not publicly available** due to strict ethical and legal restrictions regarding patient privacy. 

## 6. Disclaimer

This tool is intended for **research purposes only**. It should not be used as the sole basis for clinical decision-making. The authors and their institutions assume no liability for any medical decisions made based on the predictions of this model.

## 7. License

This work is licensed under a **Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)**.


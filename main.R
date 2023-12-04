library(tidyverse)
library(survival)
library(broom)
library(dplyr)
library(tidyr)
library(forcats)

# Assuming your CSV has a similar structure to the RDS file in the example
df0 <- read_csv("Data/Breast Cancer METABRIC.csv")

# Assuming 'df' is your dataframe
df <- na.omit(df0)

# You can now view the structure of the cleaned data and the number of rows removed
str(df)
cat("Rows before cleaning: ", nrow(df0), "\nRows after cleaning: ", nrow(df), "\n")


# ---------------------------------------------- CHECK CORRELATION MATRIX -----------

# Select only numerical columns from the dataset
numerical_data <- df %>% select_if(is.numeric)

# Compute the correlation matrix for the numerical data
correlation_matrix <- cor(na.omit(numerical_data))


#----------

# Checking for missing data
#missing_data2 <- sapply(df, function(x) sum(is.na(x)))
#missing_data2

df <- df %>%
  # Convert 'Overall Survival (Months)' to numeric if it's not already
  mutate(`Overall Survival (Months)` = as.numeric(as.character(`Overall Survival (Months)`))) %>%
  # Remove rows with missing values in 'Overall Survival Status' or 'Overall Survival (Months)'
  filter(!is.na(`Overall Survival Status`), !is.na(`Overall Survival (Months)`))

df <- df %>%
  mutate(died_of_cancer = case_when(
    `Overall Survival Status` == "Living" | (`Overall Survival Status` == "Deceased" & `Patient's Vital Status` != "Died of Disease") ~ 0,
    `Patient's Vital Status` == "Died of Disease" ~ 1,
    TRUE ~ 2
  ))

df_clean <- df %>%
  filter(died_of_cancer != 2)



# -----------------------------


# Baseline model
M0 <- coxph(Surv(`Overall Survival (Months)`, died_of_cancer) ~ 1, data = df_clean)

# Model with 'Tumor Size' as a predictor
M1 <- coxph(Surv(`Overall Survival (Months)`, died_of_cancer) ~ `Tumor Size`, data = df_clean)

# Model with 'Cellularity' as a predictor
M2 <- coxph(Surv(`Overall Survival (Months)`, died_of_cancer) ~ `Cellularity`, data = df_clean)

# Model with all predictors
M3 <- coxph(Surv(`Overall Survival (Months)`, died_of_cancer) ~ `Tumor Size` + `Cellularity` + `Age at Diagnosis` + `Type of Breast Surgery` + `Cancer Type Detailed` + `Pam50 + Claudin-low subtype` + `Nottingham prognostic index` + `Mutation Count` + `Lymph nodes examined positive`, data = df_clean)
# Model with `Tumor Size` + `Cellularity` + `Age at Diagnosis`  + `Lymph nodes examined positive`+`Relapse Free Status` predictors

M4 <- coxph(Surv(`Overall Survival (Months)`, died_of_cancer) ~ `Tumor Size` + `Cellularity` + `Age at Diagnosis`  + `Lymph nodes examined positive`+`Relapse Free Status`, data = df_clean)
# Model with `Tumor Size` + ER Status +`Cellularity` + `Age at Diagnosis`  + `Lymph nodes examined positive`+`Relapse Free Status` predictors

M5 <- (coxph(Surv(`Overall Survival (Months)`, died_of_cancer) ~ `Tumor Stage`+`Chemotherapy`+`ER Status`+`Age at Diagnosis`+`Relapse Free Status`+`Lymph nodes examined positive` , data = df_clean))
# Model with `Tumor Stage`+`Chemotherapy`+`ER Status`+`Age at Diagnosis`+`Relapse Free Status` predictors

M6 <- (coxph(Surv(`Overall Survival (Months)`, died_of_cancer) ~ `Tumor Stage`+`Chemotherapy`+`ER Status`+`Age at Diagnosis`+`Relapse Free Status` , data = df_clean))
# Model with `Tumor Stage`+`Chemotherapy`+`ER Status`+`Age at Diagnosis` predictors

M7 <- (coxph(Surv(`Overall Survival (Months)`, died_of_cancer) ~ `Tumor Stage`+`Chemotherapy`+`ER Status`+`Age at Diagnosis` , data = df_clean))

summary(M3)

# Check proportional hazards assumption
cox.zph(M3)
# -------------------------

# LRT
anova(M1, M3)

# AIC
fits <- list(M1 = M1, M2 = M2, M3 = M3,M4 = M4, M5 = M5, M6 = M6, M7 = M7 )
sapply(fits, AIC)

# -------------------

# Diagnostics for M3, replace with the model of interest
residuals_diag <- cox.zph(M3)
plot(residuals_diag)

# -----------------
dat <- df_clean
set.seed(123) # for reproducibility
train_index <- sample(1:nrow(dat), 0.7 * nrow(dat))
train_dat <- dat[train_index, ]
test_dat <- dat[-train_index, ]

# Retrain your model on the training dataset
M3_train <- coxph(Surv(`Overall Survival (Months)`, died_of_cancer) ~ `Tumor Size` + `Cellularity`, data = train_dat)

# Make predictions on the testing set
test_dat$predicted_risk <- predict(M3_train, newdata = test_dat, type = "risk")

# Evaluate the model
surv_obj <- Surv(test_dat$`Overall Survival (Months)`, test_dat$died_of_cancer)
cox_fit <- coxph(surv_obj ~ predicted_risk, data = test_dat)
summary(cox_fit)

# ------------------- tumor size + cellularity
# 
# coef (0.08633): This is the estimated coefficient for predicted_risk. It indicates the log hazard ratio, which is the expected change in the log hazard for a one-unit increase in predicted_risk.
# exp(coef) (1.09017): This is the hazard ratio (HR). A HR greater than 1 suggests a higher hazard (or risk) as the predicted_risk increases. Specifically, for each one-unit increase in predicted_risk, the hazard of the event occurring is expected to increase by 9%.
# se(coef) (0.04826): The standard error of the coefficient. It measures the variability or uncertainty in the estimated coefficient.
# z (1.789): The z-value is the coefficient divided by its standard error. It's used to test the null hypothesis that the coefficient equals 0.
# Pr(>|z|) (0.0736): The p-value for the z-test. It indicates whether the predictor is statistically significant. A common threshold for significance is 0.05. Here, the p-value is just above 0.05, indicated with a dot (.), which means it's suggestive but not conventionally statistically significant.
# The confidence interval for the hazard ratio:
#   
#   lower .95 (0.9918) and upper .95 (1.198): The 95% confidence interval for the HR does not include 1, which suggests that there may be an effect, but since the p-value is greater than 0.05, we would not typically consider this to be statistically significant.
# Model diagnostics:
#   
#   Concordance (0.629): The concordance statistic is a measure of the predictive accuracy of the model and ranges from 0.5 (no predictive ability) to 1 (perfect prediction). A value of 0.629 indicates a fair predictive ability.
# Likelihood ratio test (p=0.1): This test compares the goodness of fit of the model against a null model with no predictors. A p-value of 0.1 suggests that the model is not significantly better than the null model at the 0.05 level.
# Wald test (p=0.07) and Score (logrank) test (p=0.06): These are additional tests for the significance of the predictors. Both are marginally non-significant, suggesting that there is a trend towards predicted_risk being associated with the outcome, but the evidence is not strong enough to confirm this at the conventional p < 0.05 level.
# In conclusion, the predicted_risk variable appears to have a trend towards being a predictor of survival, but this relationship is not statistically significant at the p < 0.05 level. The model has a fair predictive ability as indicated by the concordance statistic.


#---------------------------

# Coefficients table interpretation:
#   
#   coef (0.34947): This is the estimated coefficient for the variable predicted_risk. It's the log hazard ratio (HR), indicating that for each one-unit increase in predicted_risk, the hazard increases by about 35%.
# exp(coef) (1.41832): This is the hazard ratio itself. It suggests that as predicted_risk increases by one unit, the event's hazard (e.g., death, recurrence) increases by 41.8%.
# se(coef) (0.06334): The standard error of the coefficient estimate, indicating the level of uncertainty around the coefficient estimate.
# z (5.518): The z-statistic, used in testing the null hypothesis that the coefficient equals zero (no effect). It is calculated as the coefficient divided by its standard error.
# Pr(>|z|) (< 0.0000000344 or 3.44e-08): The p-value associated with the z-statistic. Since this p-value is very small (much less than 0.05), we reject the null hypothesis and conclude that predicted_risk is a statistically significant predictor of the hazard.
# Confidence interval for the hazard ratio:
#   
#   lower .95 (1.253) & upper .95 (1.606): This interval is quite narrow and does not include 1, which further supports the significance of the predicted_risk predictor.
# Model diagnostics:
#   
#   Concordance (0.638): The concordance statistic is a measure of the model's predictive accuracy, which is 0.638. This suggests the model has a moderate ability to discriminate between individuals who will experience the event versus those who will not.
# Likelihood ratio test (p=0.00009), Wald test (p=0.0000000344), and Score (logrank) test (p=0.0000000003 or 3e-11): All these tests provide evidence against the null hypothesis, indicating that the model with predicted_risk is significantly better at explaining the variation in the survival times than the null model without any predictors.
# In summary, the predicted_risk variable is a statistically significant predictor of the hazard according to this Cox model. The model has a moderate predictive ability, and all diagnostic tests suggest that the model fits the data better than a null model without predictors.


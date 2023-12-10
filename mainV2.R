library(tidyverse)
library(survival)
library(broom)
library(dplyr)
library(tidyr)
library(forcats)
library(MASS)
library(glmnet)
library(boot)
library(rms)
library(survminer)
library(cvms)
library(survivalROC)
library(groupdata2)

# Importing Dataset
df0 <- read_csv("Breast Cancer METABRIC.csv")
df <- na.omit(df0)

# structure of the cleaned data and the number of rows removed
str(df)
cat("Rows before cleaning: ", nrow(df0), "\nRows after cleaning: ", nrow(df), "\n")


# ---------------------------------------------- CHECK CORRELATION MATRIX -----------

# Select only numerical columns from the dataset
numerical_data <- df %>% select_if(is.numeric)

# Compute the correlation matrix for the numerical data
correlation_matrix <- cor(na.omit(numerical_data))
print(correlation_matrix)

#----------

# Checking for missing data
#missing_data2 <- sapply(df, function(x) sum(is.na(x)))
#missing_data2

df <- df %>%
  # Convert 'Overall Survival (Months)' to numeric if it's not already
  mutate(`Overall Survival (Months)` = as.numeric(as.character(`Overall Survival (Months)`))) %>%
  # Remove rows with missing values in 'Overall Survival Status' or 'Overall Survival (Months)'
  filter(!is.na(`Overall Survival Status`), !is.na(`Overall Survival (Months)`))

#  (`Overall Survival Status` == "Deceased" & `Patient's Vital Status` != "Died of Disease") -> if not dead of disease then event did not occur

df <- df %>%
  mutate(died_of_cancer = case_when(
    `Overall Survival Status` == "Living"  ~ 0,
    `Patient's Vital Status` == "Died of Disease" ~ 1,
    TRUE ~ 0
  ))

df_filtered <- dplyr::select(df, -c(`Patient ID`, `Overall Survival Status`, `Relapse Free Status (Months)`, 
                                    `Relapse Free Status`, `Patient's Vital Status`))


names(df_filtered) <- gsub(" ", "_", names(df_filtered))
names(df_filtered) <- gsub("\\+", "plus", names(df_filtered))
names(df_filtered) <- gsub("-", "_", names(df_filtered))
names(df_filtered) <- gsub("\\(", "", names(df_filtered))
names(df_filtered) <- gsub("\\)", "", names(df_filtered))
names(df_filtered) <- make.names(names(df_filtered), unique = TRUE)
print(names(df_filtered))

#Died_of_cancer will be our event 


#See dsitribution of values inside the dataset :
for (var in names(df_filtered)) {
  # Check if the column is numeric
  if (is.numeric(df_filtered[[var]])) {
    # Create a histogram for numeric columns
    plot <- ggplot(df_filtered, aes_string(x = var)) + 
      geom_histogram(bins = 30, fill = "blue", color = "black") + 
      ggtitle(paste("Histogram of", var)) +
      theme_minimal()
    
  } else {
    # Create a bar plot for categorical columns
    plot <- ggplot(df_filtered, aes_string(x = var)) + 
      geom_bar(fill = "orange", color = "black") + 
      ggtitle(paste("Bar Plot of", var)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x labels for readability
  }
  
  # Print the plot
  print(plot)
}

# We remove cancer tpye and sex because there is only one category
df_filtered <- dplyr::select(df_filtered, -c("Sex", "Cancer_Type"))

# ---------------------- Univariate testing -----------------
# Define the outcome and predictor variables
surv_response <- with(df_filtered, Surv(Overall_Survival_Months, died_of_cancer))
predictors <- df_filtered[, !(names(df_filtered) %in% c("Overall_Survival_Months", "died_of_cancer"))]

# Univariate Cox Regression
univariate_results <- list()
for (var in names(predictors)) {
  formula_str <- sprintf("Surv(Overall_Survival_Months, died_of_cancer) ~ %s", var)
  formula <- as.formula(formula_str)
  model <- coxph(formula, data = df_filtered)
  univariate_results[[var]] <- summary(model)
}

# Review univariate results to decide which variables to include in the multivariate model
print(univariate_results)



#------------ Find best variable combination ------------------

#According to last part these are the significant variables
initial_formula <- as.formula("Surv(Overall_Survival_Months, died_of_cancer) ~ Type_of_Breast_Surgery + Chemotherapy + Pam50_plus_Claudin_low_subtype + ER_status_measured_by_IHC + ER_Status + Neoplasm_Histologic_Grade + HER2_status_measured_by_SNP6 + HER2_Status + Lymph_nodes_examined_positive + Nottingham_prognostic_index + PR_Status + X3_Gene_classifier_subtype + Tumor_Size + Tumor_Stage")
full_model <- coxph(initial_formula, data = df_filtered)

# Apply stepwise model selection
stepwise_model <- stepAIC(full_model, direction = "both")
print(summary(stepwise_model))



#------------------Improving the model : parameter removal --------

# X3_Gene_classifier: have p-values that are not statistically significant (p > 0.05). 
# It suggests that the specific subtypes may not be significant predictors of survival in the model.
#Let's try to recompute the model without them :

initial_formula <- as.formula("Surv(Overall_Survival_Months, died_of_cancer) ~ Type_of_Breast_Surgery + Chemotherapy + Pam50_plus_Claudin_low_subtype + ER_status_measured_by_IHC + ER_Status + Neoplasm_Histologic_Grade + HER2_status_measured_by_SNP6 + HER2_Status + Lymph_nodes_examined_positive + Nottingham_prognostic_index + PR_Status + Tumor_Size + Tumor_Stage")
full_model <- coxph(initial_formula, data = df_filtered)
# Apply stepwise model selection
Model2 <- stepAIC(full_model, direction = "both")
print(summary(Model2))

#Still better model by keeping Pam50

#------------------Trying from every parameters possible --------

initial_formula <- as.formula("Surv(Overall_Survival_Months, died_of_cancer) ~ Age_at_Diagnosis + Type_of_Breast_Surgery + Cancer_Type_Detailed + Cellularity + Chemotherapy + Pam50_plus_Claudin_low_subtype + Cohort + ER_status_measured_by_IHC + ER_Status + Neoplasm_Histologic_Grade + HER2_status_measured_by_SNP6 + HER2_Status + Tumor_Other_Histologic_Subtype + Hormone_Therapy + Inferred_Menopausal_State + Integrative_Cluster + Primary_Tumor_Laterality + Lymph_nodes_examined_positive + Mutation_Count + Nottingham_prognostic_index + Oncotree_Code + PR_Status + Radio_Therapy + X3_Gene_classifier_subtype + Tumor_Size + Tumor_Stage")
full_model <- coxph(initial_formula, data = df_filtered)
# Apply stepwise model selection
FullParameterModel <- stepAIC(full_model, direction = "both")
print(summary(FullParameterModel))

# ---------------------- Model 3 :
#According to the previous model it would be better with these parameters :
#   Age_at_Diagnosis                                 0.010567  1.010623  0.004509  2.344 0.019091 *  
#   Neoplasm_Histologic_Grade                       -0.449369  0.638031  0.133163 -3.375 0.000739 ***
#   HER2_StatusPositive                              0.594275  1.811716  0.263859  2.252 0.024307 *  
#   Hormone_TherapyYes                              -0.363205  0.695444  0.127779 -2.842 0.004477 ** 
#   Nottingham_prognostic_index                      0.613990  1.847789  0.076750  8.000 1.25e-15 ***
#   Radio_TherapyYes                                -0.303352  0.738339  0.114543 -2.648 0.008088 ** 
#   Tumor_Size                                       0.010859  1.010919  0.002492  4.358 1.31e-05 ***

#X3 gene classifier subtype still gives better result despite not being much statistically significant

Model3 <- coxph(Surv(Overall_Survival_Months, died_of_cancer) ~ Age_at_Diagnosis + Neoplasm_Histologic_Grade + HER2_Status + Hormone_Therapy + Nottingham_prognostic_index + Radio_Therapy + Tumor_Size + X3_Gene_classifier_subtype, data = df_filtered)
print(summary(Model3))

#Best model so far with Concordance of 0.711

# ----------------- Cross Validation : not working
# Define the Cox model formula
cox_formula <- as.formula("Surv(Overall_Survival_Months, died_of_cancer) ~ 
    Age_at_Diagnosis + Neoplasm_Histologic_Grade + HER2_Status + 
    Hormone_Therapy + Nottingham_prognostic_index + Radio_Therapy + 
    Tumor_Size + X3_Gene_classifier_subtype")

# Create a new data frame with the fold column
df_with_folds <- df_filtered

# Add a fold column to the new data frame using the "n_dist" method
df_with_folds$.folds <- groupdata2::fold(df_with_folds, k = num_folds)

# Create the Cox model
cox_model <- coxph(cox_formula, data = df_with_folds)

# Specify the formulas for cross-validation (you can modify these as needed)
formulas <- c("Surv(Overall_Survival_Months, died_of_cancer) ~ 
    Age_at_Diagnosis + Neoplasm_Histologic_Grade + HER2_Status + 
    Hormone_Therapy + Nottingham_prognostic_index + Radio_Therapy + 
    Tumor_Size + X3_Gene_classifier_subtype")

# Perform cross-validation using cross_validate
cv_results <- cross_validate(
  data = df_with_folds,
  formulas = formulas,
  family = "gaussian",  # Update with the appropriate family
  fold_cols = as.factor(df_with_folds$.folds),  # Use the fold column you created
)

# Print or access the cross-validation results
print(cv_results)




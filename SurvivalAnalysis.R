# Load necessary libraries
library(readr)
library(dplyr)
library(ggplot2)

# Import the CSV file
data <- read_csv("Breast Cancer METABRIC.csv")

# View the first few rows of the dataset
head(data)

# Summary statistics for numeric variables
summary(data)

# Explore the structure of the dataset
str(data)

# Count of unique values for categorical variables
sapply(data, function(x) if (is.factor(x)) length(unique(x)))

# Statistical analysis examples

# Age distribution of patients
age_distribution <- ggplot(data, aes(x=`Age at Diagnosis`)) +
  geom_histogram(binwidth=1, color='black', fill='blue') +
  labs(title="Age Distribution of Patients", x="Age at Diagnosis", y="Count")
print(age_distribution)

# Tumor size by tumor stage, assuming 'Tumor Size' and 'Tumor Stage' are the exact column names
tumor_size_stage <- ggplot(data, aes(x=`Tumor Stage`, y=`Tumor Size`)) +
  geom_boxplot() +
  labs(title="Tumor Size by Tumor Stage", x="Tumor Stage", y="Tumor Size (mm)")
print(tumor_size_stage)

# Overall Survival by Age, assuming 'Overall Survival (Months)' is the exact column name
survival_age <- ggplot(data, aes(x=`Age at Diagnosis`, y=`Overall Survival (Months)`)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(title="Overall Survival by Age", x="Age at Diagnosis", y="Overall Survival (Months)")
print(survival_age)

# Checking for missing data
missing_data <- sapply(data, function(x) sum(is.na(x)))
missing_data

# Correlation between numeric variables, assuming all other numeric column names are exact
numeric_data <- data %>% select_if(is.numeric)
correlation_matrix <- cor(na.omit(numeric_data))
print(correlation_matrix)

# Save a summary table of the statistics
summary_stats <- summary(data)
write.csv(summary_stats, file="SummaryStatistics.csv")

# This is a very basic analysis. Depending on your hypothesis and goals,
# you might want to perform more specific statistical tests or predictive modeling.



# ---------------------------------------------- CHECK CORRELATION MATRIX -----------

# Select only numerical columns from the dataset
numerical_data <- data %>% select_if(is.numeric)

# Compute the correlation matrix for the numerical data
correlation_matrix <- cor(na.omit(numerical_data))

# View the correlation matrix
print(correlation_matrix)

# If you want to view the correlation matrix as a more readable formatted table, you can use the `knitr` package
library(knitr)
kable(correlation_matrix)

#---------------------------------------------




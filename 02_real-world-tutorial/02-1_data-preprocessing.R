### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 02  | real-world-tutorial
### (1) | data-preprocessing.R
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### [NOTE]  
### This is a preprocessing script for ABCD cohort data
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### written by
### ---
### Jinwoo Lee
### SNU Connectome Lab
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(fastDummies)
library(caret)

### _______________________________________________________________________ ####
### @@@ PART I. DEFAULT SETTING @@@ ####
## > I-1. Preparing Raw Data ####
# You can access the rawdata of ABCD cohort in https://nda.nih.gov/abcd/request-access. 
data.raw <- read.csv('GRF Protocol Release5.1_new.csv')

target.env.cov.list <- c("age", "site_id_l", "sex", "parent_age", "high_educ", "income", 
                         "birth_weight", "maternal_age", "bmi","history_ratio", 
                         "race_ethnicity", "married", "parent_identity")

target.vol.cov.list <- names(data.raw)[startsWith(names(data.raw), "smri_vol")]

# using only cortical/subcortical gray-matter volume features 
target.volgm.cov.list <- target.vol.cov.list[!grepl("cbwmatter|crbwmatter|cc|csf|ventricle|bstem|wmhint|vedc|intracranialv|suprateialv|aal|aar", 
                                             target.vol.cov.list, ignore.case = TRUE)] # excluding the miscellaneous regions
target.outcome <- names(data.raw)[startsWith(names(data.raw), "cbcl")]
target.treatment <- "Bullying"

target.total.vars.list <- c(target.treatment, target.outcome, target.env.cov.list, target.volgm.cov.list)


### _______________________________________________________________________ ####
### @@@ PART II. PREPROCESSING @@@ ####
## > II-1.Filtering the Variables of Interest ####
data.raw <- data.raw[, target.total.vars.list]

## > II-2. Excluding Samples with Missing Values ####
data.raw <- na.omit(data.raw) # N = 10,815 -> 9,213

## > II-3. Defining the Types of Variables ####
data.raw$Bullying <- as.factor(data.raw$Bullying)
data.raw$site_id_l <- as.factor(data.raw$site_id_l)
data.raw$sex <- as.factor(data.raw$sex)
data.raw$high_educ <- as.factor(data.raw$high_educ)
data.raw$income <- as.factor(data.raw$income)
data.raw$race_ethnicity <- as.factor(data.raw$race_ethnicity)
data.raw$married <- as.factor(data.raw$married)
data.raw$parent_identity <- as.factor(data.raw$parent_identity)

summary(data.raw) # columns with 'too much' outliers: sex, bmi

## > II-4. Outlier Detection for BMI - continuous ####
detect_outliers_iqr <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower <- Q1 - 1.5 * IQR
  upper <- Q3 + 1.5 * IQR
  return(x < lower | x > upper)
}

bmi.outliers <- detect_outliers_iqr(data.raw$bmi)
data.raw <- data.raw[!bmi.outliers, ] # excluding 433 samples (N = 8,780)

## > II-5. Outlier Detection for Sex - Categorical ####
data.raw <- data.raw[data.raw$sex %in% c("F", "M"), ] # excluding 2 samples (N = 8,778)
data.raw$sex <- factor(data.raw$sex, levels = c("M", "F"))

# saving the dataset before the normalization for Appendix B
write.csv(data.raw, "GRFProtocol_real-world-dataset_pped-by-JW_for-appendixB.csv")

## > II-6. Volume Normalization with Whole-brain Volume ####
for (col in target.volgm.cov.list) {
  norm_col <- paste0(col, "_norm")
  data.raw[[norm_col]] <- data.raw[[col]] / data.raw$smri_vol_scs_wholeb
}

summary(data.raw)

# excluding the not-normalized volume features
data.raw <- data.raw[, !colnames(data.raw) %in% target.volgm.cov.list]

# excluding the whole-brain measure
data.raw <- data.raw[, !colnames(data.raw) %in% "smri_vol_scs_wholeb_norm"]

## > II-7. Zero Variance Check ####
nearZeroVar(data.raw, saveMetrics = TRUE) # all variables do not show near-zero-variance

## > II-8. Z-transformation for Continuous Covariates ####
numeric_vars <- sapply(data.raw, is.numeric)
data.raw[, numeric_vars] <- scale(data.raw[, numeric_vars])

## > II-9. One Hot Encoding for Categorical Covariates ####
except_cat_vars <- c('Bullying', 'site_id_l')
cat_vars <- names(data.raw)[sapply(data.raw, is.factor) & !(names(data.raw) %in% except_cat_vars)]

data.dummy <- dummy_cols(data.raw, 
                         select_columns = cat_vars,
                         remove_selected_columns = TRUE,
                         remove_first_dummy = FALSE)

data.dummy <- data.dummy[, !colnames(data.dummy) %in% "sex_F"]

# save the preprocessed dataset
write.csv(data.dummy, "GRFProtocol_real-world-dataset_pped-by-JW.csv")
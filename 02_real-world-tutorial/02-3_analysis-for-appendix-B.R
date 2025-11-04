### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 02  | real-world-tutorial
### (3) | analysis-for-appendix-B.R
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### [NOTE]  
### This is a script for all analyses in Appendix B 
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### written by
### ---
### Jinwoo Lee
### SNU Connectome Lab
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(dplyr)

### _______________________________________________________________________ ####
### @@@ PART I. DEFAULT SETTING @@@ ####
## > I-1. Preparing Preprocessed Dataset ####
data <- read.csv("GRFProtocol_real-world-dataset_pped-by-JW_for-appendixB.csv",  header = TRUE)

variables <- c("Bullying", "cbcl_dsm_depression_t", "age", "sex", "parent_age",
               "high_educ", "income", "birth_weight", "maternal_age", "bmi", "history_ratio", 
               "race_ethnicity", "married", "parent_identity")

data.focus <- data[, variables]

## > I-2. Defining the Variable Types ####
data.focus$Bullying <- as.factor(data.focus$Bullying)
data.focus$sex <- as.factor(data.focus$sex)
data.focus$high_educ <- as.factor(data.focus$high_educ)
data.focus$income <- as.factor(data.focus$income)
data.focus$race_ethnicity <- as.factor(data.focus$race_ethnicity)
data.focus$married <- as.factor(data.focus$married)
data.focus$parent_identity <- as.factor(data.focus$parent_identity)

cont_vars <- names(data.focus)[sapply(data.focus, is.numeric)]
cat_vars  <- names(data.focus)[sapply(data.focus, function(x) is.factor(x) | is.character(x))]


### _______________________________________________________________________ ####
### @@@ PART II. SUMMARY STATISTICS @@@ ####
## > II-1. for Continuous Variables ####
summarize_cont <- function(data, vars) {
  data %>%
    summarise(across(all_of(vars),
                     ~ sprintf("%.1f (%.1f)", mean(.x, na.rm = TRUE), sd(.x, na.rm = TRUE)),
                     .names="{.col}"))
}

## > II-1. for Categorical Variables ####
summarize_cat <- function(data, vars) {
  res <- list()
  for (v in vars) {
    tab <- table(data[[v]], useNA="no")
    n   <- sum(tab)
    df  <- data.frame(
      Variable = v,
      Level = names(tab),
      Count = as.numeric(tab),
      Percent = round(100 * as.numeric(tab)/n, 1)
    )
    res[[v]] <- df
  }
  bind_rows(res)
}

## > II-3. Calculating Descriptive Stats ####
groups <- list(
  All = data.focus,
  Bullying1 = data.focus %>% filter(Bullying == 1),
  Bullying0 = data.focus %>% filter(Bullying == 0)
)

summary_results <- list()

for (g in names(groups)) {
  dat <- groups[[g]]
  
  cont_summary <- summarize_cont(dat, cont_vars) %>%
    mutate(Group = g)
  
  cat_summary <- summarize_cat(dat, cat_vars) %>%
    mutate(Group = g)
  
  summary_results[[g]] <- list(Continuous = cont_summary,
                               Categorical = cat_summary)
}


### _______________________________________________________________________ ####
### @@@ PART III. STATISTICAL TESTING: TREATED vs. UNTREATED @@@ ####
p_values <- list()

## > III-1. Continuous: Student's T-Test ####
cont_tests <- lapply(cont_vars, function(v) {
  test <- t.test(data.focus[[v]] ~ data.focus$Bullying)
  data.frame(Variable = v, P_value = test$p.value, Test = "t-test")
})

## > III-2. Categorical: Fisher's Exact Test ####
cat_tests <- lapply(cat_vars, function(v) {
  tbl <- table(data.focus[[v]], data.focus$Bullying)
  if (all(dim(tbl) > 1)) {
    test <- tryCatch(
      fisher.test(tbl, workspace = 2e7),
      error = function(e) fisher.test(tbl, simulate.p.value = TRUE, B = 1e5)
    )
    data.frame(Variable = v, P_value = test$p.value, Test = "Fisher's exact")
  } else {
    data.frame(Variable = v, P_value = NA, Test = "Not applicable")
  }
})
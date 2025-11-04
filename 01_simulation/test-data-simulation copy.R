### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 01   | simulation
###      | test-data-simulation.R
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### [NOTE]  
### This is a data simulation code for test dataset.
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Test-sets are used in only backward-elimination analyses.
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### written by
### ---
### Maria Pak & Jinwoo Lee
### SNU Connectome Lab
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(grf)
library(dplyr)
library(tidyr)
library(caret)
library(tibble)
library(fastDummies)

### _______________________________________________________________________ ####
### @@@ PART I. DEFAULT SETTING @@@ ####
## > I-1. Load the Workspace with Preprocessed ABCD Data ####
# You can access the rawdata of ABCD cohort in https://nda.nih.gov/abcd/request-access. 
load("full_results.RData")

full_data <- full_join(demo, smri, by='src_subject_id') %>%
  drop_na()

## > I-2. Z-transformation for Non-nominal Variables ####
not_numeric <- c('src_subject_id','demo_sex_v2','site_id_l',
                 'demo_prnt_marital_v2_1','demo_prnt_marital_v2_2','demo_prnt_marital_v2_3',
                 'demo_prnt_marital_v2_4','demo_prnt_marital_v2_5','demo_prnt_marital_v2_6',
                 'race_ethnicity_1','race_ethnicity_2','race_ethnicity_3',
                 'race_ethnicity_4','race_ethnicity_5')

numeric_cols <- names(full_data)[!(names(full_data) %in% not_numeric)]

for(i in 2:ncol(full_data)){
  if(names(full_data)[i] %in% numeric_cols){
    full_data[,i] <- as.vector(scale(full_data[,i]))
  }
}


### _______________________________________________________________________ ####
### @@@ PART II. SIMULATION OF DATASET @@@ ####
# This script creates four simulated datasets:
# 1. data_nonlin_strong
# 2. data_nonlin_weak
# 3. data_lin_strong
# 4. data_lin_weak

## > II-1. Common Setup (Run ONCE for all datasets) ####
N <- 1000
seed <- 1
sd_1 <- 0.5
sd_2 <- 0.5

# set seed to sample the same base data for all simulations
set.seed(123 * seed)
data_base <- full_data %>% sample_n(size = N)

# define features for heterogeneity (x_het) - key moderators
# using 'ind_x' to avoid confusion with 'ind1'/'ind2' below 
ind_x1 <- which(names(data_base) == 'smri_vol_scs_amygdalarh')
ind_x2 <- which(names(data_base) == 'smri_vol_cdk_cdacaterh')
ind_x3 <- which(names(data_base) == 'smri_vol_cdk_parahpallh')

x1 <- data_base[, ind_x1]
x2 <- data_base[, ind_x2]
x3 <- data_base[, ind_x3]

x_het <- data.frame(x1, x2, x3) 

# generate common noise and treatment vectors
# > we used '123' and '1234' as random seeds for noise generation to generate 
# > different noises for the train and test sets, respectively.
set.seed(1234*seed)     
noise_1 <- rnorm(N, 0, sd_1)

set.seed(1234*seed*2)   
noise_2 <- rnorm(N, 0, sd_2) 

set.seed(123*seed)
w_prop <- 1 / (1 + exp(-(x1 + x2)))
w <- rbinom(N, prob = w_prop, size = 1)

# get indices for direct effect calculation (as in your original script)
ind1 <- which(names(data_base) == 'src_subject_id')
ind2 <- which(names(data_base) == 'site_id_l')


## > II-2. Define a Function to Generate a Single Dataset ####
# This function uses the common variables defined above
create_sim_data <- function(beta_0, beta_k, gamma, gamma_0, type) {
  
  # calculate direct effect
  dir_eff <- beta_k * rowSums(data_base[, -c(ind1, ind2)])
  
  # calculate theta (CATE) based on type
  if (type == "nonlinear") {
    theta <- gamma_0
    for(i in 1:ncol(x_het)){
      theta <- theta + gamma * cos(x_het[, i]) + gamma * sin(x_het[, i])
    }
  } else if (type == "linear") {
    # this is the linear equation you provided
    theta <- gamma_0 + gamma * rowSums(x_het)
  } else {
    stop("Invalid type. Must be 'nonlinear' or 'linear'.")
  }
  
  # add common noise_1 to theta
  theta <- theta + noise_1
  
  # calculate final outcome y
  y <- beta_0 + dir_eff + theta * w + noise_2
  
  # combine and return the final dataset
  data_out <- cbind(y, theta, w, data_base)
  return(data_out)
}

## > II-3. Generate the Four Datasets by Calling the Function ####
## > In backward-simulation analyses, we only test-sets with weak heterogeneity
## > for replication analysis (see Fig 3).
# >> II-3-a. data_nonlin_weak ####
data_nonlin_weak <- create_sim_data(
  beta_0 = 0.1, 
  beta_k = 0.4, 
  gamma = 4, 
  gamma_0 = 0.4,
  type = "nonlinear"
)

# >> II-3-b. data_lin_weak ####
data_lin_weak <- create_sim_data(
  beta_0 = 0.1, 
  beta_k = 0.7, 
  gamma = 5, 
  gamma_0 = 0.4,
  type = "linear"
)

## > II-4. (Optional) Clean Up Intermediate Variables ####
rm(data_base, x1, x2, x3, x_het, noise_1, noise_2, w, 
   ind_x1, ind_x2, ind_x3, ind1, ind2,
   N, seed, sd_1, sd_2, create_sim_data)

# save the generated dataset
write.csv(data_nonlin_weak, 'data/test_nonlin_weak.csv')
write.csv(data_lin_weak, 'data/test_lin_weak.csv')
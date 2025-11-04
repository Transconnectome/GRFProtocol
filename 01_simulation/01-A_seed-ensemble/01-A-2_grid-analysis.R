### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 01-A | simulation > seed_ensemble
### (2)  | grid_analysis.R
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### [NOTE]  
### This is a data analysis script for Fig 2c.
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### written by
### ---
### Jinwoo Lee & Junghoon Park
### SNU Connectome Lab
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(grf)
library(caret)
library(dplyr)
library(scales)
library(ggplot2)

for (i in 1:4) {
  ### _____________________________________________________________________ ####
  ### @@@ PART I. INITIAL SETUP @@@ ####
  ### > I-1. Data Preparation ####
  # assume we have simulation data using ../train_data_simulation.R.
  if (i == 1) {data <- read.csv("../data/train_lin_strong.csv")}
  else if (i == 2) {data <- read.csv("../data/train_lin_weak.csv")}
  else if (i == 3) {data <- read.csv("../data/train_nonlin_strong.csv")}
  else if (i == 4) {data <- read.csv("../data/train_nonlin_weak.csv")}
  
  print(paste0("Start running ", dataset_name))
  
  ### > I-2. Variable Definition ####
  X <- data[, c(5:18, 20:150)]    # except for site_id
  C <- as.factor(data$site_id_l)  # site_id for clustering
  Y <- as.numeric(data$y)
  W <- as.numeric(data$w)
  
  ### > I-3. Random Seeds ####
  max.seed.num <- 10   # the maximum number of seed condition
  test.num <- 50       # the number of tests in each condition (i.e., each cell in grid)
  
  set.seed(1234)
  seed.list <- sample(c(1:10000), max.seed.num * test.num, replace = FALSE)
  landmark.list <- seq(from = 1, to = max.seed.num * test.num, by = max.seed.num)
  
  n.seed.list <- c(1, 2, 5, 10)
  n.tree.list <- c(1000, 2000, 5000, 10000)
  
  ### > I-4. Results Template ####
  results.df <- data.frame(matrix(nrow = 0, ncol = 10))
  total.idx <- length(n.seed.list) * length(n.tree.list) * test.num
  loop.idx <- 0
  
  ### _____________________________________________________________________ ####
  ### @@@ PART II. MODEL FITTING: SEED > TREE > TRIAL @@@ ####      
  for (n.seed in n.seed.list) {
    for (n.tree in n.tree.list) {
      for (idx.test in c(1:test.num)) {
        
        start.t <- proc.time()
        options(warn = 1)
        
        ## > II-1. Model Fitting ####
        # single seed condition
        if (n.seed == 1) {
          tmp.seed <- seed.list[landmark.list[idx.test]]
          
          tmp.cf <- causal_forest(X, Y, W,
                                  Y.hat = NULL, W.hat = NULL,
                                  tune.parameters = "all",
                                  num.trees = n.tree,
                                  num.threads = 8,
                                  clusters = C,
                                  seed = tmp.seed)
        }
        
        # seed ensemble condition
        else {
          tmp.seed.list <- seed.list[c(landmark.list[idx.test] : (landmark.list[idx.test] + n.seed - 1))]
          mini.cf.list <- list()
          
          for (seed.idx in c(1:n.seed)) {
            tmp.mini.cf <- causal_forest(X, Y, W,
                                         Y.hat = NULL, W.hat = NULL,
                                         tune.parameters = "all",
                                         num.trees = n.tree,
                                         num.threads = 8,
                                         clusters = C,
                                         seed = tmp.seed.list[seed.idx])
            
            # saving 'mini' models
            mini.cf.list[[paste0("mini.cf_", tmp.seed.list[seed.idx])]] <- tmp.mini.cf
          }
          
          # >> II-2-b. performing seed ensemble and making a 'big' model
          tmp.cf <- merge_forests(mini.cf.list, compute.oob.predictions = FALSE)
        }
        
        ## > II-2. Calibration Test ####
        tmp.calib <- test_calibration(tmp.cf)
        beta_ATE.coef <- tmp.calib[1]
        beta_ATE.p <- tmp.calib[7]
        beta_ITE.coef <- tmp.calib[2]
        beta_ITE.p <- tmp.calib[8]
        
        ## > II-3. Evaluating ITE Prediction Performance ####
        ITE.est <- predict(tmp.cf)$predictions
        ITE.true <- data$theta
        est.true.cor <- cor.test(ITE.est, ITE.true)
        recov.coef <- as.numeric(est.true.cor$estimate)
        recov.p <- as.numeric(est.true.cor$p.value)
        
        ## > II-4. Assessing Computational Costs ####
        end.t <- proc.time()
        elapsed.dur <- end.t - start.t
        
        ## > II-5. Results Archiving ####
        results.df <- rbind(results.df, c(n.seed, n.tree, idx.test,
                                          beta_ATE.coef, beta_ATE.p,
                                          beta_ITE.coef, beta_ITE.p,
                                          recov.coef, recov.p,
                                          elapsed.dur))
        
        loop.idx <- loop.idx + 1
        
        print(paste0("[", sprintf("%03d", loop.idx), "/", total.idx, "] ", idx.test, "-th analysis with ", 
                     n.seed, " seed(s) and ", n.tree, " trees was done!"))
        
      }
    }
  }
  
  results.df <- results.df[, c(1:10)] # excluding unnecessary columns
  
  colnames(results.df) <- c("n_seed", "n_tree", "idx_test", 
                            "beta_ATE_coef", "beta_ATE_p", 
                            "beta_ITE_coef", "beta_ITE_p", 
                            "recov_coef", "recov_p", "duration")
  
  ## > II-6. Saving the Final Results for Each ITE Condition ####
  write.csv(results.df, file = paste0("results/grid-analysis-results_", dataset_name, ".csv"), row.names = FALSE)
}


### _____________________________________________________________________ ####
### @@@ PART III. VISUALIZATION - FIG 2C @@@ #### 
## > III-1. Data Preparation ####
results.df <- read.csv("results/grid-analysis-results_data_nonlin_weak.csv")

results.df$beta_ATE_p_neglog <- -log(results.df$beta_ATE_p)
results.df$beta_ITE_p_neglog <- -log(results.df$beta_ITE_p)

results.df$n_seed <- as.factor(results.df$n_seed)
results.df$n_tree <- as.factor(results.df$n_tree)

seed.list <- c('1', '2', '5', '10')
tree.list <- c('1000', '2000', '5000', '10000')

total.df <- data.frame(matrix(nrow = 0, ncol = 7))

for (seed in seed.list) {
  for (tree in tree.list) {
    
    tmp.df <- results.df[results.df$n_seed == seed & results.df$n_tree == tree, ]
    
    beta_ATE_coef_SD <- sd(tmp.df[, 'beta_ATE_coef'])
    beta_ATE_p_neglog_SD <- sd(tmp.df[, 'beta_ATE_p_neglog'])
    beta_ITE_coef_SD <- sd(tmp.df[, 'beta_ITE_coef'])
    beta_ITE_p_neglog_SD <- sd(tmp.df[, 'beta_ITE_p_neglog'])
    
    beta_ATE_coef_mean <- mean(tmp.df[, 'beta_ATE_coef'])
    beta_ITE_coef_mean <- mean(tmp.df[, 'beta_ITE_coef'])
    beta_ATE_p_neglog_mean <- mean(tmp.df[, 'beta_ATE_p_neglog'])
    beta_ITE_p_neglog_mean <- mean(tmp.df[, 'beta_ITE_p_neglog'])
    
    # coefficient of variation
    beta_ATE_coef_CV <- beta_ATE_coef_SD / beta_ATE_coef_mean
    beta_ITE_coef_CV <- beta_ITE_coef_SD / beta_ITE_coef_mean
    beta_ATE_p_neglog_CV <- beta_ATE_p_neglog_SD / beta_ATE_p_neglog_mean
    beta_ITE_p_neglog_CV <- beta_ITE_p_neglog_SD / beta_ITE_p_neglog_mean
    
    # duration
    duration_M <- mean(tmp.df[, 'duration'])
    
    total.df <- rbind(total.df, c(seed, tree, 
                                  beta_ATE_coef_CV, beta_ATE_p_neglog_CV,
                                  beta_ITE_coef_CV, beta_ITE_p_neglog_CV,
                                  duration_M))
  }
}

rm(list = ls(pattern = "^beta"), duration_M, seed, tree, tmp.df)

colnames(total.df) <- c("seed", "tree", 
                        "beta_ATE_coef_CV", "beta_ATE_p_neglog_CV",
                        "beta_ITE_coef_CV", "beta_ITE_p_neglog_CV",
                        "duration_M")

for (col in c(1:ncol(total.df))) {total.df[, col] <- as.numeric(total.df[, col])}
total.df[, "seed"] <- as.factor(total.df[, "seed"])
total.df[, "tree"] <- as.factor(total.df[, "tree"])


## > III-2. GRID (1): beta_{ATE} ####
fig2c.ATE <- ggplot(total.df, aes(x = tree, y = seed, fill = beta_ATE_p_neglog_CV)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f (%.3f)\n%.0f",
                                beta_ATE_coef_CV, beta_ATE_p_neglog_CV, duration_M)), size = 3) +
  scale_fill_gradient(limits = c(0.040, 0.550), low = "#1976D2", high = "white", name = "CV of -log(P)") + 
  labs(title = "Grid of beta_ATE", x = "# of trees", y = "# of seeds") +
  theme_minimal()

## > III-3. GRID (2): beta_{ITE} ####
fig2c.ITE <- ggplot(total.df, aes(x = tree, y = seed, fill = beta_ITE_p_neglog_CV)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f (%.3f)\n%.0f",
                                beta_ITE_coef_CV, beta_ITE_p_neglog_CV, duration_M)), size = 3) +
  scale_fill_gradient(limits = c(0.040, 0.550), low = "#1976D2", high = "white", name = "CV of -log(P)") + 
  labs(title = "Grid of beta_ITE", x = "# of trees", y = "# of seeds") +
  theme_minimal()
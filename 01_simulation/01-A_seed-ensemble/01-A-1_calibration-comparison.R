### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 01-A | simulation > seed_ensemble
### (1)  | calibration_comparison.R
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### [NOTE]  
### This is a data analysis script for Fig 2a and 2b.
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### written by
### ---
### Jinwoo Lee
### SNU Connectome Lab
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(grf)
library(dplyr)
library(tidyr)
library(ggplot2)

# assume we have simulation data using ../train_data_simulation.R.
data <- read.csv("../data/train_nonlin_weak.csv") 

### _______________________________________________________________________ ####
### @@@ PART I. DEFAULT SETTING @@@ ####
## > I-1. Random Seeds ####
# for reproduction using seeds below
set.seed(2025)
seed.list <- sample(c(1000:9999), size = 20, replace = FALSE)

ensemble.list <- split(seed.list, ceiling(seq_along(seed.list) / 5))
ensemble.list <- unname(ensemble.list)

single.seed.list <- sapply(ensemble.list, function(x) x[1])

## > I-2. Variable Definition ####
Y <- data[, "y"]
W <- data[, "w"]
X <- data[, c(5:154)]
X <- subset(X, select = -site_id_l) # excluding site variables 
C <- data$site_id_l                 # site would be used as a clustering variable)

## > I-3. Result Template ####
result.df <- data.frame(
  condition = character(),
  seed = character(),
  beta_ATE_est = numeric(),
  beta_ATE_se = numeric(),
  beta_ATE_p = numeric(),
  beta_ITE_est = numeric(),
  beta_ITE_se = numeric(),
  beta_ITE_p = numeric()
)

### _______________________________________________________________________ ####
### @@@ PART II. SINGLE SEED ANALYSIS @@@ ####
single.n.tree <- 1000

for (seed_idx in 1:length(single.seed.list)) {
  # model fitting without seed ensemble
  tmp.seed <- single.seed.list[seed_idx]
  
  # model fitting
  tmp.cf <- causal_forest(X, Y, W, 
                          Y.hat = NULL, W.hat = NULL,
                          clusters = C,
                          tune.parameters = "all",
                          compute.oob.predictions = FALSE,
                          num.trees = single.n.tree,
                          num.threads = 12,
                          seed = tmp.seed)
  
  # calibration test
  tmp.fit <- test_calibration(tmp.cf)
  tmp.beta_ATE.est <- tmp.fit[1]
  tmp.beta_ATE.se <- tmp.fit[3]
  tmp.beta_ATE.p <- tmp.fit[7]
  tmp.beta_ITE.est <- tmp.fit[2]
  tmp.beta_ITE.se <- tmp.fit[4]
  tmp.beta_ITE.p <- tmp.fit[8]
  
  # result archiving
  result.df <- rbind(result.df, c("not ensembled", tmp.seed,
                                  tmp.beta_ATE.est, tmp.beta_ATE.se, tmp.beta_ATE.p,
                                  tmp.beta_ITE.est, tmp.beta_ITE.se, tmp.beta_ITE.p))
  
  print(paste0("A model under seed ", tmp.seed, " was assessed!"))
}

colnames(result.df) <- c("condition", "seeds", "beta_ATE_est", "beta_ATE_se", "beta_ATE_p", "beta_ITE_est", "beta_ITE_se", "beta_ITE_p")
rm(list = ls(pattern = "^tmp"))

### _______________________________________________________________________ ####
### @@@ PART III. SEED ENSEMBLED CONDITION @@@ ####
for (list.idx in 1:length(ensemble.list)) {
  tmp.seed.list <- ensemble.list[[list.idx]]
  tmp.mini.cf.list <- list()
  
  for (seed.idx in 1:length(tmp.seed.list)) {
    tmp.seed <- tmp.seed.list[seed.idx]
    
    # iteratively fitting 'mini' models
    tmp.mini.cf <- causal_forest(X, Y, W, 
                                 Y.hat = NULL, W.hat = NULL,
                                 clusters = C,
                                 tune.parameters = "all",
                                 compute.oob.predictions = FALSE,
                                 num.trees = single.n.tree/length(ensemble.list),
                                 num.threads = 12,
                                 seed = tmp.seed)
    
    # saving 'mini' models
    tmp.mini.cf.list[[paste0("mini.cf.", seed.idx)]] <- tmp.mini.cf
  }
  
  # performing seed ensembles and making 'big' model
  tmp.big.cf <- merge_forests(tmp.mini.cf.list, compute.oob.predictions = FALSE)
  
  # calibration test with 'big' model
  tmp.fit <- test_calibration(tmp.big.cf)
  tmp.beta_ATE.est <- tmp.fit[1]
  tmp.beta_ATE.se <- tmp.fit[3]
  tmp.beta_ATE.p <- tmp.fit[7]
  tmp.beta_ITE.est <- tmp.fit[2]
  tmp.beta_ITE.se <- tmp.fit[4]
  tmp.beta_ITE.p <- tmp.fit[8]
  
  # result archiving
  result.df <- rbind(result.df, c("ensembled", c(paste(tmp.seed.list, collapse = ", "),
                                  tmp.beta_ATE.est, tmp.beta_ATE.se, tmp.beta_ATE.p,
                                  tmp.beta_ITE.est, tmp.beta_ITE.se, tmp.beta_ITE.p)))
                     
  print("A model ensembled across five seeds was assessed!")
}

rm(list = ls(pattern = "^tmp"))
result.df[, c(1:2)] <- lapply(result.df[, c(1:2)], as.factor)
result.df[, c(3:8)] <- lapply(result.df[, c(3:8)], as.numeric)

result.df[result.df$condition == 'not ensembled', 'beta_ATE_p_fdr'] <- p.adjust(result.df[result.df$condition == 'not ensembled', 'beta_ATE_p'], method = 'fdr')
result.df[result.df$condition == 'not ensembled', 'beta_ITE_p_fdr'] <- p.adjust(result.df[result.df$condition == 'not ensembled', 'beta_ITE_p'], method = 'fdr')
result.df[result.df$condition == 'ensembled', 'beta_ATE_p_fdr'] <- p.adjust(result.df[result.df$condition == 'ensembled', 'beta_ATE_p'], method = 'fdr')
result.df[result.df$condition == 'ensembled', 'beta_ITE_p_fdr'] <- p.adjust(result.df[result.df$condition == 'ensembled', 'beta_ITE_p'], method = 'fdr')

write.csv(result.df, file = "results/calibration-comparison-results.csv", row.names = FALSE)


### _____________________________________________________________________ ####
### @@@ PART IV. VISUALIZATION - FIG 2A/B @@@ #### 
## IV-1. Fig 2a (single-seed) ####
fig2a.df <- result.df[result.df$condition == 'not ensembled', ] %>%
  mutate(
    beta_ATE_lower = beta_ATE_est - 1.96 * beta_ATE_se,
    beta_ATE_upper = beta_ATE_est + 1.96 * beta_ATE_se,
    beta_ITE_lower = beta_ITE_est - 1.96 * beta_ITE_se,
    beta_ITE_upper = beta_ITE_est + 1.96 * beta_ITE_se
  ) %>%
  select(seeds,
         beta_ATE_est, beta_ATE_lower, beta_ATE_upper,
         beta_ITE_est, beta_ITE_lower, beta_ITE_upper) %>%
  pivot_longer(
    cols = -seeds,
    names_to = c("type", ".value"),
    names_pattern = "beta_(.*)_(.*)"
  )

fig2a.df$trials <- rep(c("trial 1", "trial 2", "trial 3", "trial 4"), each = 2)

fig2a <- ggplot(fig2a.df, aes(x = est, y = trials, color = type)) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(width = 0.6),
                width = 0.0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("ATE" = "lightgray", "ITE" = "#1976D2")) +
  labs(
    x = "Calibration Estimates (95% CI)",
    y = "Trials",
    color = "Type"
  ) +
  xlim(-1.0, 3.0) + 
  theme_classic()

## IV-2. Fig 2b (seed-ensembled) ####
fig2b.df <- result.df[result.df$condition == 'ensembled', ] %>%
  mutate(
    beta_ATE_lower = beta_ATE_est - 1.96 * beta_ATE_se,
    beta_ATE_upper = beta_ATE_est + 1.96 * beta_ATE_se,
    beta_ITE_lower = beta_ITE_est - 1.96 * beta_ITE_se,
    beta_ITE_upper = beta_ITE_est + 1.96 * beta_ITE_se
  ) %>%
  select(seeds,
         beta_ATE_est, beta_ATE_lower, beta_ATE_upper,
         beta_ITE_est, beta_ITE_lower, beta_ITE_upper) %>%
  pivot_longer(
    cols = -seeds,
    names_to = c("type", ".value"),
    names_pattern = "beta_(.*)_(.*)"
  )

fig2b.df$trials <- rep(c("trial 1", "trial 2", "trial 3", "trial 4"), each = 2)

fig2b <- ggplot(fig2b.df, aes(x = est, y = trials, color = type)) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(width = 0.6),
                width = 0.0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("ATE" = "lightgray", "ITE" = "#1976D2")) +
  labs(
    x = "Calibration Estimates (95% CI)",
    y = "Trials",
    color = "Type"
  ) +
  xlim(-1.0, 3.0) + 
  theme_classic()
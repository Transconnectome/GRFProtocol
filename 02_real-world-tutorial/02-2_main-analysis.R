### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 02  | real-world-tutorial
### (2) | main-analysis.R
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### [NOTE]  
### This is a script for all analyses in Fig 4. 
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### written by
### ---
### Jinwoo Lee
### SNU Connectome Lab
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(grf)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(purrr)

### _______________________________________________________________________ ####
### @@@ PART I. DEFAULT SETTING @@@ ####
## > I-1. Preparing Preprocessed Dataset ####
data <- read.csv("GRFProtocol_real-world-dataset_pped-by-JW.csv",  header = TRUE)

outcome_var <- "cbcl_dsm_depression_t"
Y <- data[, outcome_var]

treatment_var <- "Bullying"
W <- data[, treatment_var]

X <- data[, c(12:150)]
X <- X[, !colnames(X) %in% "site_id_l"] # excluding the site variable from covariates

site_id_l <- as.numeric(factor(data$site_id_l))
site.aov.model <- aov(cbcl_dsm_depression_t ~ site_id_l, data = data)
summary(site.aov.model)


### _______________________________________________________________________ ####
### @@@ PART II. BACKWARD ELIMINATION-BASED MODEL SELECTION @@@ ####
## > II-1. Set the Results Template ####
result.models <- data.frame(var_num = numeric(),
                            beta_ATE_est = numeric(), beta_ATE_p = numeric(),
                            beta_ITE_est = numeric(), beta_ITE_p = numeric(),
                            min_e_hat = numeric(), max_e_hat = numeric())
varimp.list <- colnames(X)
worst.var.list <- list()

n <- nrow(data)
full.var.num <- length(X)

## > II-2. Set-up for Seed Ensemble ####
n.seeds <- 5
n.trees <- 2000

seeds <- c(6406, 1553, 2976, 4302, 5960) # for fixing the seeds across analyses

## > II-3. Main Iterations ####
for (threshold in c(0:(full.var.num-1))) {
  # basic setting
  options(warn = 1)
  current.var.num <- full.var.num - threshold
  chosen.feature <- varimp.list
  
  assign(paste0("var_list_", current.var.num), chosen.feature)
  
  cf.thr.seed.list <- list() # for archiving mini forests
  
  X.tmp <- subset(X, select = chosen.feature)
  
  # fitting the multiple 'mini forests' under each random seed in 'seeds'
  for (seed_idx in c(1:length(seeds))) {
    current.cf.thr.seed <- causal_forest(X.tmp, Y, W, 
                                         Y.hat = NULL, W.hat = NULL,
                                         clusters = site_id_l,
                                         tune.parameters = "all",
                                         compute.oob.predictions = FALSE,
                                         num.trees = n.trees,
                                         num.threads = 12,
                                         seed = seeds[seed_idx])
    
    cf.thr.seed.list[[paste0("cf.thr.seed", seed_idx)]] = current.cf.thr.seed
  }
  
  # merging mini forests to the single 'big forest'
  cf.thr.big <- merge_forests(cf.thr.seed.list, compute.oob.predictions = FALSE)
  print(paste0("[", threshold+1, "/", full.var.num, "] STEP 1 - The 'big forest' with the ", current.var.num, " variables was created."))
  
  # cleaning the environment
  rm(cf.thr.seed.list)
  
  # assessing the model fit with the big forest
  fit.thr <- test_calibration(cf.thr.big)
  
  # archiving the results with the big forest
  beta_ATE.est <- fit.thr[1]
  beta_ATE.p <- fit.thr[7]
  beta_ITE.est <- fit.thr[2]
  beta_ITE.p <- fit.thr[8]
  min.e.hat <- min(cf.thr.big$W.hat)
  max.e.hat <- max(cf.thr.big$W.hat)
  
  result.models <- rbind(result.models, c(current.var.num, 
                                          beta_ATE.est, beta_ATE.p,
                                          beta_ITE.est, beta_ITE.p, 
                                          min.e.hat, max.e.hat))
  
  print(paste0("[", threshold+1, "/", full.var.num, "] STEP 2 - The results from the big forest with the ", current.var.num, " variables were saved."))
  
  # updating the variable importance for the next iteration
  current.varimp <- c(variable_importance(cf.thr.big))
  names(current.varimp) <- colnames(X.tmp)
  current.varimp <- sort(current.varimp, decreasing = TRUE)
  worst.var <- names(current.varimp)[current.var.num]
  varimp.list <- varimp.list[varimp.list != worst.var] # exclude the feature with the lowest importance
  worst.var.list <- append(worst.var.list, worst.var)  # save the order of the excluded variables
  
  print(paste0("[", threshold+1, "/", full.var.num, "] STEP 3 - The covariates for next iteration were successfully updated."))
  print(paste0("[", threshold+1, "/", full.var.num, "] The excluded variable is ", worst.var, " and the number of left covariates is ", length(varimp.list), "."))
}

colnames(result.models) <- c("var_num", "beta_ATE_est", "beta_ATE_p", "beta_ITE_est", "beta_ITE_p", "min_e_hat", "max_e_hat")

## > II-4. Model Selection based on Calibration Test and Overlap Assumptions ####
result.models$beta_ATE_p_fdr <- p.adjust(result.models$beta_ATE_p, "fdr")
result.models$beta_ITE_p_fdr <- p.adjust(result.models$beta_ITE_p, "fdr")
result.models$ate_p_fdr <- p.adjust(result.models$ATE_p, "fdr")
result.models$fit_index <- abs(1 - result.models$beta_ATE_est) + abs(1 - result.models$beta_ITE_est)

write.csv(result.models, paste0("GRF Protocol Results_", treatment_var, "_", outcome_var, ".csv"), row.names = FALSE) # Save iterated models as CSV file

result.models.pass <- result.models[result.models$beta_ATE_p_fdr < .05 & 
                                      result.models$beta_ITE_p_fdr < .05 & # models passing the calibration test
                                      result.models$min_e_hat >= .05 & 
                                      result.models$max_e_hat < .95, ]     # and fulfilling the overlap assumption

best.model.var_num <- result.models.pass[which.min(result.models.pass$fit_index), 'var_num'] 
best.model.var <- get(paste0('var_list_', as.character(best.model.var_num)))

selected <- best.model.var
selected.raw <- result.models[result.models$var_num == length(selected), ]


### _______________________________________________________________________ ####
### @@@ PART III. GATE TEST FOR HETEROGENEITY ASSESSMENT @@@ ####
## > III-1. Fitting the Best Model Again #### 
X.best <- subset(X, select = selected)

cf.best.seed.list <- list() # for archiving mini forests

for (seed_idx in c(1:length(seeds))) {
  current.cf.best <- causal_forest(X.best, Y, W, 
                                   Y.hat = NULL, W.hat = NULL,
                                   clusters = site_id_l,
                                   tune.parameters = "all",
                                   compute.oob.predictions = FALSE,
                                   num.trees = n.trees,
                                   num.threads = 12,
                                   seed = seeds[seed_idx])
  
  cf.best.seed.list[[paste0("cf.best.seed", seed_idx)]] <- current.cf.best
}

# merging mini forests into the big forest
cf.best.big <- merge_forests(cf.best.seed.list, compute.oob.predictions = TRUE)

# check whether the re-fitting yields the same results
test_calibration(cf.best.big)

## > III-2. ATE Estimation #### 
ATE_est <- average_treatment_effect(cf.best.big)[[1]]
ATE_se <- average_treatment_effect(cf.best.big)[[2]]
ATE_tstat <- ATE_est / ATE_se
ATE_p <- 1.96 * (1 - pnorm(abs(ATE_tstat)))

## > III-3. Predicting ITE in an Out-of-Bag Manner #### 
ite.pred <- predict(cf.best.big)$predictions
ite.pred.tertile <- quantile(ite.pred, seq(0, 1, by = 1/3))
hist(ite.pred, breaks = 100)

## > III-4. Estimating GATEs for Three Groups ####
result.GATE.tertile <- data.frame(matrix(ncol = 4, nrow = 0))
ite.pred.tertile

for (group_idx in c(1:3)) {
  # define the subjects' boolean condition for each quantile (TRUE if subject in Q_n)
  group.conditions <- c()
  
  for (sub_idx in c(1:n)) {
    sub.group.condition <- ifelse(group_idx == 1, ite.pred[sub_idx] < ite.pred.tertile[[2]],
                                  ifelse(group_idx == 3, ite.pred[sub_idx] >= ite.pred.tertile[[3]],
                                         ite.pred[sub_idx] >= ite.pred.tertile[[group_idx]] & 
                                           ite.pred[sub_idx] < ite.pred.tertile[[group_idx + 1]]))
    
    group.conditions <- append(group.conditions, sub.group.condition)
  }
  
  # assign values to results matrix
  current_group <- paste0("Q", group_idx)
  current_GATE_est <- average_treatment_effect(cf.best.big, subset = group.conditions)[[1]]
  current_GATE_se <- average_treatment_effect(cf.best.big, subset = group.conditions)[[2]]
  current_GATE_p <- 1.96 * (1 - pnorm(abs(current_GATE_est / current_GATE_se)))
  
  result.GATE.tertile <- rbind(result.GATE.tertile, c(current_group, current_GATE_est, current_GATE_se, current_GATE_p))
}

colnames(result.GATE.tertile) <- c("group", "GATE_est", "GATE_se", "GATE_p")
result.GATE.tertile$group <- as.factor(result.GATE.tertile$group)
result.GATE.tertile$GATE_est <- as.numeric(result.GATE.tertile$GATE_est)
result.GATE.tertile$GATE_se <- as.numeric(result.GATE.tertile$GATE_se)
result.GATE.tertile$GATE_p <- as.numeric(result.GATE.tertile$GATE_p)

## > III-5. Statistical Testing with Estimated Three GATEs ####
group.pairs <- list(c("Q2", "Q1"),
                    c("Q3", "Q1"),
                    c("Q3", "Q2"))

result.GATE.posthoc <- data.frame(
  contrast = character(),
  diff = numeric(),
  se_diff = numeric(),
  t = numeric(),
  p = numeric(),
  stringsAsFactors = FALSE
)

for (pair in group.pairs) {
  g1 <- pair[1]
  g2 <- pair[2]
  
  est1 <- result.GATE.tertile$GATE_est[result.GATE.tertile$group == g1]
  se1  <- result.GATE.tertile$GATE_se[result.GATE.tertile$group == g1]
  
  est2 <- result.GATE.tertile$GATE_est[result.GATE.tertile$group == g2]
  se2  <- result.GATE.tertile$GATE_se[result.GATE.tertile$group == g2]
  
  diff  <- est1 - est2
  se_diff <- sqrt(se1^2 + se2^2)
  t_stat <- diff / se_diff
  p_val  <- 1.96 * (1 - pnorm(abs(t_stat))) 
  
  result.GATE.posthoc <- rbind(result.GATE.posthoc, data.frame(
    contrast = paste0(g1, " - ", g2),
    diff = diff,
    se_diff = se_diff,
    t = t_stat,
    p = p_val
  ))
}

result.GATE.posthoc$p_fdr <- p.adjust(result.GATE.posthoc$p, method = 'fdr')


### _______________________________________________________________________ ####
### @@@ PART IV. IDENTIFICATION OF KEY MODERATORS @@@ ####
## > IV-1. Group Comparison ####
X.best$pred_ITE <- ite.pred
X.best <- X.best %>%
  mutate(group = case_when(
    pred_ITE <= ite.pred.tertile[[2]] ~ "Q1",
    pred_ITE <= ite.pred.tertile[[3]] ~ "Q2",
    pred_ITE <= ite.pred.tertile[[4]] ~ "Q3",
  )) %>%
  mutate(group = factor(group, levels = c("Q1", "Q2", "Q3")))

result.group.comp <- list()

for (var in selected) {
  formula <- as.formula(paste(var, "~ group"))
  model <- aov(formula, data = X.best)
  
  summary_model <- summary(model)[[1]]
  
  result <- data.frame(
    variable = var,
    df_between = summary_model[1, "Df"],
    df_within = summary_model[2, "Df"],
    F_value = summary_model[1, "F value"],
    p_value = summary_model[1, "Pr(>F)"]
  )
  
  result.group.comp[[var]] <- result
}

result.group.comp <- do.call(rbind, result.group.comp)
result.group.comp$p_value_fdr <- p.adjust(result.group.comp$p_value, method = 'fdr')
rownames(result.group.comp) <- NULL

# summary statistics
summary_table <- X.best %>%
  select(group, all_of(selected)) %>%
  pivot_longer(cols = -group, names_to = "variable", values_to = "value") %>%
  group_by(group, variable) %>%
  summarise(
    mean = round(mean(value, na.rm = TRUE), 3),
    sd   = round(sd(value, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = paste0(mean, " (", sd, ")")) %>%
  select(group, variable, mean_sd) %>%
  pivot_wider(names_from = group, values_from = mean_sd)

## > IV-2. Best Linear Projection ####
result.blp <- best_linear_projection(cf.best.big, X.best[c(1:length(selected))])
result.blp <- as.data.frame(result.blp[, ])
result.blp <- tibble::rownames_to_column(result.blp, var = "covariate")
colnames(result.blp) <- c("covariate", "blp_coef", "blp_se", "blp_t", "blp_p")

result.blp$ci_lower <- result.blp$blp_coef - 1.96*result.blp$blp_se
result.blp$ci_upper <- result.blp$blp_coef + 1.96*result.blp$blp_se
result.blp$blp_coef <- round(result.blp$blp_coef, 3)
result.blp$ci_lower <- round(result.blp$ci_lower, 3)
result.blp$ci_upper <- round(result.blp$ci_upper, 3)
result.blp$blp_p <- round(result.blp$blp_p, 3)

## > IV-3. Partial Dependence Simulation ####
result.part.dep <- data.frame(matrix(ncol = 0, nrow = 100))
part.dep.model <- formula(paste0("~ 0 +", paste0(selected, collapse = "+")))

for (cov_idx in c(1:length(selected))) {
  current_cov <- selected[cov_idx]
  other_covs <- selected[which(selected != current_cov)]
  grid.size <- 100
  grid <- seq(min(X.best[, current_cov]), max(X.best[, current_cov]), length.out = grid.size)
  
  medians <- apply(X.best[, other_covs, F], 2, median) 
  grid.frame <- data.frame(sapply(medians, function(x) rep(x, grid.size)), grid)
  colnames(grid.frame) <- c(other_covs, current_cov)
  
  grid.X <- model.matrix(part.dep.model, grid.frame)
  simulated_pred_ITE <- predict(cf.best.big, newdata = grid.X)$predictions
  
  result.part.dep <- cbind(result.part.dep, simulated_pred_ITE)
}

colnames(result.part.dep) <- selected


### _______________________________________________________________________ ####
### @@@ PART V. VISUALIZATION FOR FIG 4 @@@ ####
## > Fig 4a. Causal Model ####
# This panel was manually visualized by Adobe Illustrator.

## > Fig 4b. Calibration Tests Across Iterations ####
fig4b <- ggplot(data = result.models, aes(x = var_num)) +
  geom_line(aes(y = beta_ATE_est, group = 1), color = "coral1", size = 0.5) + 
  geom_line(aes(y = beta_ITE_est, group = 1), color = "#2196F3", size = 0.5) + 
  geom_vline(xintercept = length(selected), linetype = "dashed") + 
  geom_point(data = subset(result.models, var_num %in% result.models[result.models$beta_ATE_p < .05 & result.models$beta_ITE_p <.05, 'var_num']),
             aes(x = var_num, y = rep(-2, length(var_num))), color = "lightgray", size = 2) + 
  geom_point(data = subset(result.models, var_num %in% result.models.pass$var_num), 
             aes(x = var_num, y = rep(-2, length(var_num))), color = "#9E9E9E", size = 2) + 
  geom_hline(yintercept = 1, linetype = "dotted") + 
  scale_y_continuous(limits = c(-2.5, 2.0)) +
  labs(x = "Number of Covariates", y = "Calibration Test Coefficients", title = "Calibration Test") + 
  scale_x_reverse() +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

## > Fig 4c. Key Moderators ####
# This panel was manually visualized by BrainNetViewer.

## > Fig 4d. Distribution of Predicted ITEs ####
fig4d.df <- data.frame(ite = ite.pred)

fig4d <- ggplot(fig4d.df, aes(x = ite)) +
  geom_histogram(
    bins = 100,
    fill = "lightgray",
    color = NA
  ) +
  geom_vline(xintercept = c(ite.pred.tertile[[2]], ite.pred.tertile[[3]]),
             color = c("black", "black"),
             linetype = "dotted",
             size = 0.75) +
  theme_classic()

## > Fig 4e. GATE Statistics ####
fig4e <- ggplot(result.GATE.tertile) +
  aes(x = group, y = GATE_est) + 
  geom_point(position=position_dodge(0.4)) +
  geom_errorbar(aes(ymin = GATE_est - 1.96*GATE_se, ymax = GATE_est + 1.96*GATE_se), 
                width = .2, 
                position = position_dodge(0.4)) +
  ylab("Group ATE estimates") + xlab("Group") +
  ggtitle("GATE test") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_hline(yintercept = 0, linetype = 'dotted')

## > Fig 4f. Group Comparison ####
fig4f.df <- X.best %>%
  select(1:length(selected), group) %>%
  pivot_longer(cols = 1:length(selected), names_to = "variable", values_to = "value")

fig4f.df.stat <- fig4f.df %>%
  group_by(group, variable) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

fig4f <- ggplot(fig4f.df.stat, aes(x = group, y = mean)) +
  geom_point(size = 0.50, color = "black") +  
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se),
                width = 0.0, color = "black") +  
  facet_wrap(~ variable, nrow = 4, ncol = 2) +  
  scale_y_continuous(limits = c(-0.5, 0.5)) +  
  labs(x = "Group", y = "95% CI") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 9),    
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

vul.factor <- c("history_ratio", "smri_vol_cdk_ehinallh_norm", "smri_vol_cdk_iftmlh_norm", "smri_vol_cdk_parsopclh_norm")
res.factor <- c("bmi", "smri_vol_cdk_precnlh_norm", "smri_vol_cdk_postcnlh_norm", "smri_vol_scs_amygdalalh_norm")

## > Fig 4g. Partial Dependence Simulation ####
fig4g.df <- result.part.dep %>%
  pivot_longer(cols = starts_with(selected), names_to = "covariate", values_to = "simulated_pred_ITE")

fig4g.df$covariate <- as.factor(fig4g.df$covariate)
fig4g.df$index <- rep(seq(0.01, 1, by = 0.01), each = length(selected))

fig4g.vul <- ggplot(fig4g.df[fig4g.df$covariate %in% vul.factor, ], 
                    aes(x = index, y = simulated_pred_ITE, color = covariate)) +
  geom_line() +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "Percentage", y = "Simulated ITEs", title = "Partial dependence")

fig4g.res <- ggplot(fig4g.df[fig4g.df$covariate %in% res.factor, ], 
                    aes(x = index, y = simulated_pred_ITE, color = covariate)) +
  geom_line() +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "Percentage", y = "Simulated ITEs", title = "Partial dependence")

write.csv(result.part.dep, "result-partial-dependency.csv")
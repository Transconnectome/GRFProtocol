### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### 01-A | simulation > seed_ensemble
### (3)  | ITE-prediction-comparison.R
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### [NOTE]  
### This is a data analysis script for Fig 2d.
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### written by
### ---
### Jinwoo Lee
### SNU Connectome Lab
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(ggplot2)
library(dplyr)
library(car)

### _____________________________________________________________________ ####
### @@@ PART I. INITIAL SETUP @@@ ####
## > I-1. Data Preparation
df.nonlin_weak <- read.csv("results/grid-analysis-results_data_nonlin_weak.csv")

target_comb <- data.frame(
  n_seed = c(1, 2, 5, 10),
  n_tree = c(10000, 5000, 2000, 1000)
)

df.chosen <- df.nonlin_weak %>%
  inner_join(target_comb, by = c("n_seed", "n_tree"))

df.chosen <- df.chosen %>%
  mutate(pair = case_when(
    n_seed == 1  & n_tree == 10000 ~ "1/10000",
    n_seed == 2  & n_tree == 5000  ~ "2/5000",
    n_seed == 5  & n_tree == 2000  ~ "5/2000",
    n_seed == 10 & n_tree == 1000  ~ "10/1000",
    TRUE ~ NA_character_  
  )) %>%
  filter(!is.na(pair))

df.chosen$pair <- factor(df.chosen$pair, levels = c("1/10000", "2/5000", "5/2000", "10/1000"))


### _____________________________________________________________________ ####
### @@@ PART II. BIAS COMPARSION @@@ ####
bias.test <- summary(aov(recov_coef ~ pair, data = df.chosen))
bias.test

### _____________________________________________________________________ ####
### @@@ PART III. VARIANCE COMPARSION @@@ ####
## > III-1. Levene's Testing ####
var.test <- leveneTest(recov_coef ~ pair, data = df.chosen)
var.test

## > III-2. Post-Hoc Testing ####
pair.levels <- levels(df.chosen$pair)
pair.comb <- combn(pair.levels, 2, simplify = FALSE)

levene.f <- c()
levene.df1 <- c()
levene.df2 <- c()
levene.p <- c()

for (pair in pair.comb) {
  subset_data <- df.chosen[df.chosen$pair %in% pair, ]
  test_result <- leveneTest(recov_coef ~ pair, data = subset_data)
  
  levene.f <- c(levene.f, test_result$`F value`[1])
  levene.df1 <- c(levene.df1, test_result$Df[1])
  levene.df2 <- c(levene.df2, test_result$Df[2])
  levene.p <- c(levene.p, test_result$`Pr(>F)`[1])
}

posthoc.levene.df <- data.frame(
  comparison = sapply(pair.comb, function(p) paste(p, collapse = " vs ")),
  f_stats <- levene.f,
  df1 <- levene.df1,
  df2 <- levene.df2,
  p_uncorrected = levene.p,
  p_fdr = p.adjust(levene.p, method = "fdr")
)

df.chosen$pair2 <- rep(c("low", "high"), each = 100)
df.chosen$pair2 <- as.factor(df.chosen$pair2)

var.test.lv2 <- leveneTest(recov_coef ~ pair2, data = df.chosen) # (1/10000, 2/5000) vs. (5/2000, 10/1000)
var.test.lv2

### _____________________________________________________________________ ####
### @@@ PART III. VISUALIZATION - FIG 2D @@@ ####
df.stats <- df.chosen %>%
  group_by(pair) %>%
  summarise(
    mean_coef = mean(recov_coef),
    se_coef = sd(recov_coef) / sqrt(n()),  # standard error
    .groups = "drop"
  )

fig2d <- ggplot() +  
  geom_jitter(data = df.chosen, 
              aes(x = pair, y = recov_coef), 
              width = 0.1, alpha = 0.1, size = 1.0) +
  geom_point(data = df.stats, 
             aes(x = pair, y = mean_coef), 
             color = "#1976D2", size = 2.0) +
  geom_errorbar(data = df.stats, 
                aes(x = pair, 
                    ymin = mean_coef - 1.96*se_coef, 
                    ymax = mean_coef + 1.96*se_coef),
                width = 0.0, color = "#1976D2") +
  labs(x = "Grid pair", y = "Pearson's R", 
       title = "") +
  theme_classic()
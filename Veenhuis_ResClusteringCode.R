# ============================================================
# Code for:
# Identification of distinct HIV reservoir phenotypes and their associated immune landscapes
#
# Manuscript authorship:
# - Main author: Rebecca Veenhuis
# - Code and analysis: Raha Dastgheyb, Aparna Bhattacharyya
#
# Purpose:
# - Reproducible script to run PCA on reservoir variables,
#   fit latent profile models (tidyLPA), and export cluster assignments.
#
# What this script assumes:
# 1) The raw data Excel file exists at: Data/Raw/MergedReservoirData.xlsx
# 2) The codebook exists at: Data/Clean/<codebook_filename>
# 3) VarTypes includes columns:
#    - Variable (variable name in df_Raw)
#    - Reservoir_subset (1 indicates reservoir variables used for clustering)
# 4) df_Raw contains an ID column (e.g., study_id) used to link outputs back to people.
# 5) df_Revalued has already been cleaned enough that your reservoir vars are numeric and without missingness.
#

# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(SciDataReportR)
  # To install SciDataReportR run:
  # remotes::install_github("RDastgh1/SciDataReportR")
  library(tidyLPA)
  library(mclust)
  library(ggplot2)
  library(here)
  library(readxl)
  library(kableExtra)
  library(paletteer)
  
})

# Reproducibility seed ----
SEED <- 93421
set.seed(SEED)

# Paths / filenames ----
raw_data_path <- here("Data", "Raw", "MergedReservoirData.xlsx")
codebook_filename <- "BrainCodeBook_2025_09_28_AB.xlsx"
codebook_path <- here("Data", "Clean", codebook_filename)

# Load data ----
df_Raw <- read_xlsx(raw_data_path)         # Raw merged data
VarTypes <- read_xlsx(codebook_path)       # Codebook to relabel and filter variables

# Revalue / relabel data using SciDataReportR ----
df_Revalued <- SciDataReportR::RevalueData(df_Raw, VarTypes)

# Define reservoir variables to cluster  ----
# These are the reservoir variables to cluster, as flagged in the codebook.
vars_Reservoir <- VarTypes %>%
  filter(Reservoir_subset == 1) %>%
  pull(Variable) %>%
  unique()

# Define analysis dataset ----
# Assumption: df_Revalued is already the intended analysis set (filtered, de-duplicated, etc).
# If you need a specific inclusion rule (complete cases, visit selection, etc),
# implement it here.
df_Revalued_complete <- df_Revalued

# Log transform reservoir variables ----
# Note: log1p(x) computes log(1 + x)
# Assumption: selected columns are numeric and non-negative as appropriate.
df_Reservoir_log <- df_Revalued_complete %>%
  select(all_of(vars_Reservoir)) %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(across(everything(), log1p)) %>%
  as.data.frame()

# PCA ----

PCAObj <- SciDataReportR::CreatePCATable(
  Data = df_Reservoir_log,
  VarsToReduce = vars_Reservoir,
  minThresh = 0.8
)

Fig_S1A <- PCAObj$p_scree
Fig_S1B <- PCAObj$Lollipop

# Clustering (Latent Profile Analysis via tidyLPA) ----
X <- PCAObj$Scores

# Loop through different clusters and paramaters
model_mclust_PCA <- estimate_profiles(
  X,
  n_profiles = 2:10,
  models = c(1, 2, 3)
)

d <- compare_solutions(
  model_mclust_PCA,
  statistics = c("AIC", "BIC", "Entropy", "LogLik")
)

# Label model forms for readability in plots/tables
d$fit$Model <- factor(
  d$fit$Model,
  levels = c(1, 2, 3),
  labels = c(
    "Equal Variance, Covariance = 0",
    "Varying Variance, Covariance = 0",
    "Equal Variance, Equal Covariance"
  )
)

model_data <- d$fits %>%
  as.data.frame() %>%
  pivot_longer(cols = c("AIC", "BIC", "Entropy", "LogLik"))

model_data$Model <- factor(
  model_data$Model,
  levels = c(1, 2, 3),
  labels = c(
    "1:Equal Variance, Covariance = 0",
    "2:Varying Variance, Covariance = 0",
    "3:Equal Variance, Equal Covariance"
  )
)

Fig_S1C <- model_data %>%
  ggplot(aes(x = Classes, y = value, color = Model)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free", ncol = 1) +
  theme(legend.position = "top") +
  labs(x = "Number of Clusters") +
  theme_bw() +
  scale_color_manual(values = c("#1553f4", "#ac2ca1", "#5ff733"))

# Choose final model parameters ---
# Note: Based on Figure S1C, Model 2 with 5 clusters was chosen because it showed the best balance of minimizing AIC/BIC while maintaining high Entropy.
model_ideal <- 2
k_5 <- 5

# Final model fit ----
set.seed(12345) 
Model_2_5 <- estimate_profiles(X, n_profiles = k_5, models = model_ideal)

# Print fit table for record
get_fit(Model_2_5) %>% kable()

# Extract cluster assignments ----
ClusterAssignments <- get_data(Model_2_5) %>% pluck("Class")

### Probability Distributions of Cluster Assignment (Fig_S1D)

# 1. Prepare probability data from Model_2_5
df_probabilities5_all <- get_data(Model_2_5) %>%
  # Extract posterior probabilities for all clusters
  mutate(
    maxProb = apply(select(., starts_with("CPROB")), 1, max),
    row_id  = as.character(rownames(.)),
    ID      = df_Revalued_complete$ID,
    Study   = df_Revalued_complete$Study,
    cluster = as.factor(Class)
  )

# Posterior Probability

Fig_S1D <- df_probabilities5_all %>%
  ggplot(aes(x = maxProb, fill = cluster)) +
  # Density plot to show distribution of assignment certainty
  geom_density(alpha = 0.6) +
  # Threshold line at 0.75 for manuscript reporting
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "black") +
  # Individual data points (jittered at bottom) to show N per cluster
  geom_jitter(aes(y = 0),
              height = 0.01,
              color = "black",
              size = 1.5,
              alpha = 0.7) +
  # Facet by cluster 
  facet_wrap(~paste("Cluster", Class), 
             scales = "free_y", 
             nrow = 1) +
  labs(
    title = "Distribution of Maximum Assignment Probabilities by Cluster",
    subtitle = "Dashed line indicates 75% probability threshold used for robust assignment",
    x = "Maximum Assignment Probability",
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_paletteer_d("PrettyCols::Bold")

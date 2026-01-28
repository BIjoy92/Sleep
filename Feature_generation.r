# ==============================================================================
# Title: Wearable-Derived Sleep Signatures as Digital Biomarkers of Cognitive Decline in Alzheimer’s Disease
# Author: Hyeonseul Park
# Description: Processing wearable-derived hypnogram and RMSSD data 
#              to identify digital biomarkers for MCI and Dementia.
# ==============================================================================

# 0. Global Settings & Libraries -----------------------------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(zoo)
library(stats)
library(broom.mixed)
library(lme4)
library(lmerTest)

options(scipen = 999)

# 1. Utility Functions ---------------------------------------------------------

#' Calculate Z-score with Epsilon to prevent division by zero
calc_zscore <- function(obs, mu, sigma, eps = 0.1) {
  sigma_adj <- if_else(sigma == 0, eps, sigma)
  return((obs - mu) / sigma_adj)
}

#' Extract Run-Length Encoding (RLE) Statistics
get_rle_stats <- function(data_vector) {
  r <- rle(data_vector)
  list(
    mean_run = mean(r$lengths),
    max_run  = max(r$lengths),
    deep_stage = r$values[which.max(r$lengths)]
  )
}

#' Apply Moving Window Statistics (Median, SD, Mean, CV)
apply_rolling_stats <- function(df, target_vars, width = 7, min_obs = 4) {
  df %>%
    group_by(email) %>%
    arrange(date) %>%
    mutate(across(all_of(target_vars), list(
      median_win = ~rollapplyr(.x, width, median, fill = NA, na.rm = TRUE),
      sd_win     = ~rollapplyr(.x, width, sd,     fill = NA, na.rm = TRUE),
      mean_win   = ~rollapplyr(.x, width, mean,   fill = NA, na.rm = TRUE),
      nobs       = ~rollapplyr(!is.na(.x), width, sum, fill = NA)
    ), .names = "{.col}_{.fn}")) %>%
    filter(if_all(ends_with("_nobs"), ~ .x >= min_obs)) %>%
    select(-ends_with("_nobs"))
}

# 2. Data Loading & Initial Cleaning -------------------------------------------
# setwd("D:/Lifelog/") # Set your path accordingly

sleep_cycle    <- read.csv("Contain_RMSSD_zero_Cycle_merged.csv")
sleep_variable <- read.csv("sleep_night_to_analysis.csv")

# Standardize Columns
sleep_cycle <- sleep_cycle %>%
  select(email, date, epoch_time, cycle_id, sleep_hypnogram_5min, sleep_rmssd_5min, 
         Group = 8, Group_code = 9) %>%
  mutate(date = as.Date(date))

pheno <- sleep_cycle %>% distinct(email, Group, Group_code)

# 3. Baseline Definition (Healthy Controls) ------------------------------------
good_subjects_filtered <- sleep_cycle %>%
  filter(Group == "Stable_CN") %>%
  mutate(epoch_time = as.POSIXct(epoch_time, format = "%Y-%m-%d %H:%M:%S", tz = "Asia/Seoul")) %>%
  arrange(email, date, epoch_time) %>%
  group_by(email, date) %>%
  summarise(
    n_cycles = n_distinct(cycle_id),
    allow = {
      s <- sleep_hypnogram_5min
      idx <- which(s == 4)
      length(idx) == 0 || all(idx %in% c(1, length(s)))
    }, .groups = "drop"
  ) %>%
  filter(allow, between(n_cycles, 3, 7)) %>%
  select(email, date)

# 4. Insomnia Scoring (PCA Based) ----------------------------------------------
insomnia_raw <- sleep_cycle %>%
  filter(cycle_id == 1) %>%
  group_by(email, date) %>%
  summarise(
    awake_min = {
      r <- rle(sleep_hypnogram_5min)
      (if (r$values[1] == 4) r$lengths[1] else 0) * 5
    }, .groups = "drop"
  ) %>%
  left_join(sleep_variable %>% select(email, date, SOL = sleep_onset_latency), by = c("email", "date")) %>%
  mutate(SOL_min = SOL / 60,
         insomnia_score = awake_min + SOL_min)

# PCA for Insomnia PC1
pca_mat <- insomnia_raw %>% 
  select(awake_min, SOL_min) %>% 
  drop_na() %>% 
  scale()
pc_res  <- prcomp(pca_mat)
insomnia_raw$insomnia_pc1 <- NA
insomnia_raw[!is.na(insomnia_raw$awake_min) & !is.na(insomnia_raw$SOL_min), "insomnia_pc1"] <- pc_res$x[,1]

# 5. Hypnogram Feature Engineering ---------------------------------------------

# A. RLE & Stability Features
cycle_feats <- sleep_cycle %>%
  group_by(email, date, cycle_id) %>%
  summarise(
    stats = list(get_rle_stats(sleep_hypnogram_5min)),
    transitions = sum(diff(sleep_hypnogram_5min) != 0),
    .groups = "drop"
  ) %>%
  unnest_wider(stats)

# B. Weighted Quality Score (Run-Frag)
baseline_stats <- cycle_feats %>%
  semi_join(good_subjects_filtered, by = c("email", "date")) %>%
  summarise(across(c(mean_run, transitions), list(mu = mean, sd = sd), .names = "{.fn}_{.col}"))

cycle_feats_final <- cycle_feats %>%
  mutate(
    z_stable = calc_zscore(mean_run, baseline_stats$mu_mean_run, baseline_stats$sd_mean_run),
    z_frag   = calc_zscore(transitions, baseline_stats$mu_transitions, baseline_stats$sd_transitions),
    run_frag_score = z_stable - z_frag # Simplified for example
  ) %>%
  group_by(email, date) %>%
  summarise(mean_run_frag = mean(run_frag_score, na.rm = TRUE),
            cv_run_frag   = sd(run_frag_score, na.rm = TRUE) / mean(run_frag_score, na.rm = TRUE),
            n_cycles      = n(), .groups = "drop") %>%
  filter(n_cycles > 1)

# 6. Autonomic Function (RMSSD) ------------------------------------------------
rmssd_baseline <- sleep_cycle %>%
  semi_join(good_subjects_filtered, by = c("email", "date")) %>%
  filter(sleep_hypnogram_5min %in% 1:4) %>%
  group_by(sleep_hypnogram_5min) %>%
  summarise(mu_b = mean(sleep_rmssd_5min, na.rm = TRUE),
            sd_b = sd(sleep_rmssd_5min, na.rm = TRUE), .groups = "drop")

rmssd_daily <- sleep_cycle %>%
  filter(sleep_hypnogram_5min %in% 1:4) %>%
  left_join(rmssd_baseline, by = "sleep_hypnogram_5min") %>%
  group_by(email, date, sleep_hypnogram_5min) %>%
  summarise(
    obs_mean = mean(sleep_rmssd_5min, na.rm = TRUE),
    z_score  = calc_zscore(obs_mean, first(mu_b), first(sd_b)),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = sleep_hypnogram_5min, values_from = z_score, names_prefix = "z_RMSSD_")

# 7. Final Integration & Moving Window -----------------------------------------
final_df <- insomnia_raw %>%
  left_join(cycle_feats_final, by = c("email", "date")) %>%
  left_join(rmssd_daily, by = c("email", "date")) %>%
  left_join(pheno, by = "email") %>%
  filter(Group %in% c("Stable_CN", "Stable_MCI", "Dem"))

# Feature list for Windowing
target_vars <- names(final_df)[!names(final_df) %in% c("email", "date", "Group", "Group_code")]

# Apply Windowing
df_win_all <- apply_rolling_stats(final_df, target_vars, width = 7, min_obs = 4)

# Calculate CV for all windowed variables
for (v in target_vars) {
  df_win_all[[paste0(v, "_cv_win")]] <- df_win_all[[paste0(v, "_sd_win")]] / df_win_all[[paste0(v, "_mean_win")]]
}

# 8. Export --------------------------------------------------------------------
write.csv(df_win_all, "Final_Sleep_Features_Windowed.csv", row.names = FALSE)

message("Processing complete. Final dataset contains ", nrow(df_win_all), " rows.")

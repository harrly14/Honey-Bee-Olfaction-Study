# ==============================================================================
# Fits GLMMs for:
#   - first arm choice (binary)
#   - proportion of time in treatment arm
#   - visit count to treatment and control arms
# using a set of pre-specified candidate models, then compares them by AICc.
# All the pre-specified candidate models were intentionally chosen by
# considering what combinations are biologically meaningful and feasible.
# Requires 01_clean_data.R to have been run.
# Outputs results/model_results.rds, results/model_comparison_tables.txt
# ==============================================================================

library(lme4)
library(glmmTMB)
library(MuMIn)
library(dplyr)

trial_data <- readRDS(here::here("data", "processed", "trial_data_clean.rds"))

# model comparisons and selection here!! old code pasted below


choice_glmm_full <- glmer(
  chose_trt ~ trt_arm + time_of_day_hrs_z +
    (1 | batch_id) + (1 | location),
  data = trial_data,
  family = binomial,
  na.action = "na.fail"
)
choice_glmm_best <- glmer(
  chose_trt ~ 1 +
    (1 | batch_id) + (1 | location),
  data = trial_data,
  family = binomial,
  na.action = "na.fail"
)

time_glmm_full <- glmmTMB(
  adj_prop_trt_time_secs ~ trt_arm + time_of_day_hrs_z +
    (1 | batch_id) + (1 | location),
  data = trial_data,
  family = beta_family(),
  na.action = "na.fail"
)

time_glmm_best <- glmmTMB(
  adj_prop_trt_time_secs ~ 1 +
    (1 | batch_id) + (1 | location),
  data = trial_data, 
  family = beta_family(),
  na.action = "na.fail"
)

visits_glmm_full <- glmer(cbind(trt_visits, ctrl_visits) ~ trt_arm +
                            time_of_day_hrs_z +
                            (1 | batch_id) + (1 | location),
                          data = trial_data,
                          family = binomial,
                          na.action = "na.fail"
)
visits_glmm_best <- glmer(cbind(trt_visits, ctrl_visits) ~ 1 +
                            (1 | batch_id) + (1 | location),
                          data = trial_data,
                          family = binomial,
                          na.action = "na.fail"
)

# ==============================================================================
dir.create(here::here("results"), showWarnings = FALSE, recursive = TRUE)

saveRDS(
    list(choice = choice_results, time = time_results, visits = visits_results),
    here::here("results", "model_results.rds")
)

message("Model results saved and written to results/")
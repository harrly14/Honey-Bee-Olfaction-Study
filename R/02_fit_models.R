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

library(glmmTMB)
library(MuMIn)
library(dplyr)

trial_data <- readRDS(here::here("data", "processed", "trial_data_clean.rds"))

# Each response variable gets the same set of 8 candidate models:
# - null: no fixed effects
# - only starvation_duration_mins_z
# - only trt_arm
# - only sodium_level
# - starv+arm
# - starv+salt
# - arm+salt
# - full
# glmmTMB must be used for time preference, because glmer does not support
# beta_family, so it has been used for all the models over glmer for consistency
fit_candidates <- function(response_var, data, family) {
  create_formula <- function(fixed_effects) {
    as.formula(paste(response_var, "~", fixed_effects, "+ (1|parent_batch_id)"))
  }
  models <- list(
    null = glmmTMB(formula = create_formula("1"), data = data, family = family, na.action = "na.fail"),
    starv = glmmTMB(formula = create_formula("starvation_duration_mins_z"), data = data, family = family, na.action = "na.fail"),
    arm = glmmTMB(formula = create_formula("trt_arm"), data = data, family = family, na.action = "na.fail"),
    salt = glmmTMB(formula = create_formula("sodium_level"), data = data, family = family, na.action = "na.fail"),
    starv_arm  = glmmTMB(formula = create_formula("starvation_duration_mins_z + trt_arm"), data = data, family = family, na.action = "na.fail"),
    starv_salt = glmmTMB(formula = create_formula("starvation_duration_mins_z + sodium_level"), data = data, family = family, na.action = "na.fail"),
    arm_salt = glmmTMB(formula = create_formula("trt_arm + sodium_level"), data = data, family = family, na.action = "na.fail"),
    full = glmmTMB(formula = create_formula("starvation_duration_mins_z + trt_arm + sodium_level"), data = data, family = family, na.action = "na.fail")
  )

  selection_table <- model.sel(models)
  
  # best model should be the one with the best AICc
  # OR the simplest model with a delta AIC <= 2 to the best model
  model_options <- filter(as.data.frame(selection_table), delta <= 2)
  # 'df' represents the parameter count. because the table is already sorted by
  # AICc, the first match picks the simplest model with the best AICc
  best_model_name <- rownames(model_options)[which.min(model_options$df)]

  message(sprintf("Best model selected for %s: %s", response_var, best_model_name))
  
  # append selection table and best model
  c(models, list(selection_table = selection_table, best = models[[best_model_name]]))
}

choice_results <- fit_candidates(
  "chose_trt", data = trial_data, family = binomial
)
time_results <- fit_candidates(
  "adj_prop_trt_time_secs", data = trial_data, family = beta_family()
)
visits_results <- fit_candidates(
  "cbind(trt_visits, ctrl_visits)", data = trial_data, family = binomial
)

dir.create(here::here("results"), showWarnings = FALSE, recursive = TRUE)

saveRDS(
  list(choice = choice_results, time = time_results, visits = visits_results),
  here::here("results", "model_results.rds")
)

message("Model results saved and written to results/")
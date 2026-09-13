# ==============================================================================
# Summarizes the fitted candidate models and produces the numerical results
# used in the manuscript, including:
# - proportions for choice, time, and visits
# - AICc model-selection tables for each response variable
# - selected-model estimates and 95% confidence intervals
# Requires 01_clean_data.R and 02_fit_models.R to have been run.
# Uses data/processed/trial_data_clean.rds and results/model_results.rds.
# ==============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(knitr)

# descriptive stats
trial_data <- readRDS(here::here("data", "processed", "trial_data_clean.rds"))
choice_prop <- mean(trial_data$chose_trt, na.rm = TRUE)
time_prop <- mean(trial_data$prop_trt_time_secs, na.rm = TRUE)
visits_prop <- sum(trial_data$trt_visits, na.rm = TRUE) / sum(
  trial_data$trt_visits + trial_data$ctrl_visits,
  na.rm = TRUE
)


# make model selection table
model_results <- readRDS(here::here("results", "model_results.rds"))

make_model_table <- function(results, response_name) {

  results$selection_table |>
    as.data.frame() |>
    rownames_to_column("Model name") |>

    select(
      `Model name`,
      delta,
      logLik,
      df
    ) |>

    # Best-supported model first
    arrange(delta) |>
    mutate(
      Response = response_name
    ) |>
    select(
      Response,
      `Model name`,
      delta,
      logLik,
      df
    )
}

results_table <- bind_rows(
  make_model_table(
    model_results$choice,
    "Choice"
  ),
  make_model_table(
    model_results$time,
    "Time"
  ),
  make_model_table(
    model_results$visits,
    "Visits"
  )
)

results_table <- results_table |>
  mutate(
    delta = round(delta, 2),
    logLik = round(logLik, 2)
  )

# Only show response name on first row of each section
results_table <- results_table |>
  group_by(Response) |>
  mutate(
    Response = if_else(
      row_number() == 1,
      Response,
      ""
    )
  ) |>
  ungroup()

results_table <- results_table |>
  rename(
    `ΔAICc` = delta,
    `Log likelihood` = logLik
  )
# Markdown + the GDocifyMd extension was the easiest way to get
# this formatted quickly in Google Docs (at least that i found)
markdown_table <- knitr::kable(
  results_table,
  format = "pipe",
  align = c("l", "l", "r", "r", "r")
)
markdown_table


# selected model stats
choice_model <- model_results$choice$best
time_model <- model_results$time$best
visits_model <- model_results$visits$best

choice_est <- plogis(summary(choice_model)$coefficients$cond["(Intercept)", "Estimate"])
time_est <- plogis(summary(time_model)$coefficients$cond["(Intercept)", "Estimate"])
visits_est <- plogis(summary(visits_model)$coefficients$cond["(Intercept)", "Estimate"])

choice_confint <- plogis(confint(choice_model, parm = "beta_", method = "wald"))
time_confint <- plogis(confint(time_model,   parm = "beta_", method = "wald"))
visits_confint <- plogis(confint(visits_model, parm = "beta_", method = "wald"))

message("Choice proportion: ", round(choice_prop, 3))
message("Time proportion: ", round(time_prop, 3))
message("Visits proportion: ", round(visits_prop, 3))
message(
  "Choice: estimate = ", round(choice_est, 3),
  ", 95% CI [", round(choice_confint[1, 1], 3),
  ", ", round(choice_confint[1, 2], 3), "]"
)

message(
  "Time: estimate = ", round(time_est, 3),
  ", 95% CI [", round(time_confint[1, 1], 3),
  ", ", round(time_confint[1, 2], 3), "]"
)

message(
  "Visits: estimate = ", round(visits_est, 3),
  ", 95% CI [", round(visits_confint[1, 1], 3),
  ", ", round(visits_confint[1, 2], 3), "]"
)
library(dplyr)
library(tidyr)
library(tibble)
library(knitr)

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


choice_model <- model_results$choice$best
time_model <- model_results$time$best
visits_model <- model_results$visits$best

# get confints (back transformed)
choice_confint <- plogis(confint(choice_model, parm = "beta_", method = "wald"))
time_confint <- plogis(confint(time_model,   parm = "beta_", method = "wald"))
visits_confint <- plogis(confint(visits_model, parm = "beta_", method = "wald"))


# get p values
choice_pvalues <- summary(choice_model)$coefficients$cond[, "Pr(>|z|)"]
time_pvalues <- summary(time_model)$coefficients$cond[, "Pr(>|z|)"]
visits_pvalues <- summary(visits_model)$coefficients$cond[, "Pr(>|z|)"]
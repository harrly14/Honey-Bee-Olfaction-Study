# ==============================================================================
# Reads the raw trial, collection, and batch CSVs from data/, joins them, and
# creates the variables used downstream by the models and plots.
# Outputs data/processed/trial_data_clean.{rds,csv}
# ==============================================================================

library(lubridate)
library(dplyr)
library(tidyr)

trial_data <- read.csv(here::here("data", "trial_data.csv"))
collection_data <- read.csv(here::here("data", "collection_data.csv"))
batch_data <- read.csv(here::here("data", "batch_data.csv"))

exclusion_rate <- mean(trial_data$is_excluded)
trial_data <- filter(trial_data, !is_excluded)


# helper function used in the following mutates
parse_datetime <- function(date, time) {
  as.POSIXct(paste(date, time), format = "%m/%d/%Y %I:%M:%S %p")
}

trial_data <- trial_data |>
  mutate(
    across(c(species, cage_id, batch_id, trt_arm, first_choice), as.factor),

    starving_start = parse_datetime(start_date, starving_since),
    trial_start = parse_datetime(start_date, start_time),

    # kept as a factor for possible use as a fixed/random effect later
    start_date = as.factor(start_date),

    across(c(latency, left_time, right_time), \(x) period_to_seconds(ms(x))),

    starvation_duration_mins = as.numeric(difftime(trial_start,
                                                   starving_start,
                                                   units = "mins")),

    chose_trt = as.numeric(as.character(first_choice) == as.character(trt_arm)),

    # replace L/R columns with trt/ctrl
    trt_visits = ifelse(trt_arm == "L", left_visits, right_visits),
    ctrl_visits = ifelse(trt_arm == "R", left_visits, right_visits),
    trt_time_secs = ifelse(trt_arm == "L", left_time, right_time),
    ctrl_time_secs = ifelse(trt_arm == "R", left_time, right_time),

    # batch num needs to be gotten from batch_id. e.g. "4a" -> "4"
    parent_batch_id = sub("^([0-9]+).*", "\\1", as.character(batch_id)),
    cage_id_char = as.character(cage_id),
    batch_id_char = as.character(batch_id)
  ) |>
  left_join(
    collection_data |>
      mutate(
        cage_id_char = as.character(cage_id),
        parent_batch_id = as.character(parent_batch_id)
      ) |>
      select(location, flower, cage_id_char, parent_batch_id),
    by = c("cage_id_char", "parent_batch_id")
  ) |>
  left_join(
    batch_data |>
      mutate(
        batch_id_char = as.character(batch_id),
        pct_sodium_error = (target_sodium - actual_sodium) / target_sodium,
        sodium_add_datetime = parse_datetime(sodium_add_date, sodium_add_time),
      ) |>
      select(pct_sodium_error, batch_id_char, sodium_add_datetime),
    by = "batch_id_char"
  ) |>
  mutate(
    time_since_sodium_added = as.numeric(difftime(trial_start,
                                                  sodium_add_datetime,
                                                  units = "mins")),
    time_of_day_hrs = hour(trial_start) + minute(trial_start) / 60,
    across(c(starvation_duration_mins,
             temp, pct_sodium_error,
             time_since_sodium_added,
             time_of_day_hrs
             ),
           \(x) as.numeric(scale(x)),
          .names = "{.col}_z"
          ),
    prop_trt_time_secs = trt_time_secs / (trt_time_secs + ctrl_time_secs),
    # adjusted for beta regression
    adj_prop_trt_time_secs = (prop_trt_time_secs * (n() - 1) + 0.5) / n()
  ) |>
  arrange(trial_start) |>
  tibble::rowid_to_column("trial_ID") |>
  select(
    -starving_since,
    -start_time,
    -left_visits,
    -right_visits,
    -left_time,
    -right_time
    )


dir.create(here::here("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(trial_data, here::here("data", "processed", "trial_data_clean.rds"))
write.csv(trial_data, here::here("data", "processed", "trial_data_clean.csv"), row.names = FALSE)

message("Data has been cleaned and written to data/processed/")
library(lubridate)
library(dplyr)
library(ggplot2)
library(lme4)
library(performance)

trial_data <- read.csv("trial_data.csv")
collection_data <- read.csv("collection_data.csv")
batch_data <- read.csv("batch_data.csv")

exclusion_rate <- mean(trial_data$is_excluded)
trial_data <- filter(trial_data, !is_excluded)

nrow(trial_data)

trial_data <- trial_data %>%
  filter(!is.na(start_date)) %>%
  mutate(
    across(c('species', 'cage_id', 'batch_id', 'trt_arm', 'first_choice'), as.factor),
    
    # convert dates times to date times
    starving_start = as.POSIXct(
      paste(start_date, starving_since), 
      format="%m/%d/%Y %I:%M:%S %p"
      ),
    trial_start = as.POSIXct(
      paste(start_date, start_time), 
      format="%m/%d/%Y %I:%M:%S %p"
      ),
    
    # convert durations/times to number of seconds
    latency = period_to_seconds(ms(latency)),
    left_time = period_to_seconds(ms(left_time)),
    right_time = period_to_seconds(ms(right_time)),
    
    # add new columns for analysis
    starvation_duration_mins = as.numeric(difftime(trial_start, starving_start, units = "mins")),
    chose_trt = first_choice == trt_arm,
    
    # replace L/R columns with trt_ctrl
    trt_visits = ifelse(trt_arm == "L", left_visits, right_visits),
    ctrl_visits = ifelse(trt_arm == "R", left_visits, right_visits),
    trt_time_secs = ifelse(trt_arm == "L", left_time, right_time),
    ctrl_time_secs = ifelse(trt_arm == "R", left_time, right_time)
  ) %>%
  # remove replaced columns
  select(-start_date, -starving_since, -start_time, -(left_visits:right_time))

prop_chose_trt = mean(trial_data$chose_trt)
preference_binomial_test <- binom.test(x = sum(trial_data$chose_trt == TRUE), 
                                       n = nrow(trial_data), 
                                       p = 0.5,
                                       alternative = "greater", 
                                       conf.level = 0.95)

# do the joins with collection and batch data
trial_data <- trial_data %>% 
  mutate(
    parent_batch_id = sub("^([0-9]+).*", "\\1", as.character(batch_id)), # "4a" -> "4"
    cage_id_char = as.character(cage_id),
    batch_id_char = as.character(batch_id)
  ) %>% 
  left_join(
    collection_data %>% 
      mutate(
        cage_id_char = as.character(cage_id),
        parent_batch_id = as.character(parent_batch_id)
      ) %>% 
      select(location, flower, cage_id_char, parent_batch_id), 
    by = c("cage_id_char", "parent_batch_id")
  ) %>% 
  left_join(
    batch_data %>% 
      mutate(
        batch_id_char = as.character(batch_id),
        pct_sodium_error = (target_sodium - actual_sodium)/target_sodium,
        sodium_add_datetime = as.POSIXct(
          paste(sodium_add_date, sodium_add_time), 
          format="%m/%d/%Y %I:%M:%S %p"
        )
      ) %>% 
      select(pct_sodium_error, batch_id_char, sodium_add_datetime), 
    by = "batch_id_char"
  ) %>% 
  mutate(
    time_since_sodium_added = as.numeric(difftime(trial_start, sodium_add_datetime, units = "mins")),
    time_of_day_hrs = hour(trial_start) + minute(trial_start)/60,
    ) %>%
  arrange(trial_start) %>%
  tibble::rowid_to_column("trial_ID")

nrow(trial_data) # should match pre-join row count

# standardize scaling of continuous predictors
trial_data <- trial_data %>% 
  mutate(
    z_starvation_duration_mins = as.numeric(scale(starvation_duration_mins)),
    z_temp = as.numeric(scale(temp)),
    z_pct_sodium_error = as.numeric(scale(pct_sodium_error)),
    z_time_since_sodium_added = as.numeric(scale(time_since_sodium_added)),
    z_time_of_day_hrs = as.numeric(scale(time_of_day_hrs))
  )

preference_model <- glmer(
  chose_trt ~ z_starvation_duration_mins + z_temp + z_pct_sodium_error + 
    z_time_since_sodium_added + z_time_of_day_hrs +
    (1 | batch_id) + (1 | location) + (1 | flower),  
  data = trial_data, family = binomial, na.action = "na.fail"
)
summary(preference_model)
check_collinearity(preference_model)

# next: dredge using the dredge function then refit the model using the dredge results using: 
  # 1. dredge results <- dredge(...)
  # 2. best model <- get.models(dredge results, 1)[[1]]
  # use saveRDS() on the good model
# things i may want to store from the model: 
  # fixed effect coefficients, in log-odds: fixef(sodium_model)
  # same, converted to odds ratios (poster-friendly): exp(fixef(sodium_model))
  # confidence intervals on those odds ratios: exp(confint(sodium_model, method = "Wald"))
  # full coefficient table: estimate, SE, z, p - this is your reportable stats table: coefs <- summary(sodium_model)$coefficients
  # random effect variances - tells you how much batch/location/flower actually mattered: VarCorr(sodium_model)
  # the specific "your headline result" number: is baseline preference different from 50/50? intercept row of `coefs` - its p-value tests exactly this: coefs["(Intercept)", ]
  # model fit quality, for comparing candidate models later (this is what dredge ranks by): AIC(sodium_model)
  # individual predicted probabilities per bee, for building your predicted-probability plot: predict(sodium_model, type = "response")
# plot model results with a simple bar/point chart with CI showing overall proportion choosing treatment vs. 0.5 chance line to start 
# predict(model, type = "response") or the ggeffects/sjPlot package's ggpredict() 
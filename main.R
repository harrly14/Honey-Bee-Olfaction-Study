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
    
    # i'm also keeping start date as a factor for use in the model later
    start_date = as.factor(start_date),
    
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
  select(-starving_since, -start_time, -(left_visits:right_time))

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
    z_time_of_day_hrs = as.numeric(scale(time_of_day_hrs)),
    prop_trt_time_secs = trt_time_secs / (15*60),
    prop_ctrl_time_secs = ctrl_time_secs / (15*60)
  )

preference_model <- glmer(
  chose_trt ~ trt_arm + start_date + z_time_of_day_hrs + 
    (1 | batch_id) + (1 | location),  
  data = trial_data, family = binomial, na.action = "na.fail"
)
summary(preference_model)
check_collinearity(preference_model)
check_model(preference_model)

plot_data <- data.frame(
  proportion = preference_binomial_test$estimate,
  lower = preference_binomial_test$conf.int[1],
  upper = preference_binomial_test$conf.int[2]
)

ggplot(plot_data, aes(x = "Sodium arm", y = proportion)) +
  geom_col() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    x = NULL,
    y = "Proportion choosing sodium arm",
    title = "First-choice preference for sodium-enriched arm"
  )

trial_data %>%
  ggplot(aes(x = "", y = chose_trt)) +
  geom_jitter(height = 0) +
  stat_summary(fun.data = "mean_cl_boot")

hist(trial_data$prop_trt_time_secs)

library(glmmTMB)

trt_time_model <- glmmTMB(
  prop_trt_time_secs ~ trt_arm + start_date + z_time_of_day_hrs + 
    (1 | batch_id) + (1 | location),  
  data = trial_data %>% filter(!(prop_trt_time_secs == 0)), family = beta_family(), na.action = "na.fail"
)

summary(preference_model)
check_collinearity(preference_model)

wilcox.test(trial_data$trt_time_secs, trial_data$ctrl_time_secs, paired = TRUE)

# next: dredge using the dredge function then refit the model using the dredge results using: 
  # 1. dredge results <- dredge(...)
  # 2. best model <- get.models(dredge results, 1)[[1]]
  # use saveRDS() on the good model

# predict(model, type = "response") or the ggeffects/sjPlot package's ggpredict()

# also: investigate if theres a relationship between time spent in trt arm and control arm

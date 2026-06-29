library(lubridate)
library(dplyr)

trial_data <- read.csv("y_maze_trial_data.csv") %>%
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
    
    latency = period_to_seconds(ms(latency)),
    left_time = period_to_seconds(ms(left_time)),
    right_time = period_to_seconds(ms(right_time)),
    
    # add starvation duration col
    # chose_treatment col (first choice == trt_arm)
    # add total arm visits col?
    # add total time in arms col?
    # add trt_time and comtrol_time cols? then remove L/R cols?
    # add preference time col? trt_time - control_time
    # preference ratio col? trt_time / control_time
    # visitation ratio col? trt visits / total visits
    # add unique ID col?
  ) %>%
  select(-start_date, -starving_since, -start_time) %>% # no longer needed because of above datetime conversions

exclusion_rate <- mean(trial_data$is_excluded)
trial_data <- trial_data %>% filter(!is_excluded)

# ==================== next step ====================
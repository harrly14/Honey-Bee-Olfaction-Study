library(lubridate)
library(dplyr)
library(ggplot2)

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
    
    # convert durations/times to number of seconds
    latency = period_to_seconds(ms(latency)),
    left_time = period_to_seconds(ms(left_time)),
    right_time = period_to_seconds(ms(right_time)),
    
    # add new columns for analysis
    starvation_duration_hrs = as.numeric(trial_start - starving_start),
    chose_trt = first_choice == trt_arm,
    
    # replace L/R columns with trt_ctrl
    trt_visits = ifelse(trt_arm == "L", left_visits, right_visits),
    ctrl_visits = ifelse(trt_arm == "R", left_visits, right_visits),
    trt_time_secs = ifelse(trt_arm == "L", left_time, right_time),
    ctrl_time_secs = ifelse(trt_arm == "R", left_time, right_time)
  ) %>%
  # remove replaced columns
  select(-start_date, -starving_since, -start_time, -(left_visits:right_time)) 

trial_data <- trial_data[order(trial_data$trial_start),] %>% tibble::rowid_to_column("trial_ID")

exclusion_rate <- mean(trial_data$is_excluded)
trial_data <- trial_data %>% filter(!is_excluded)

prop_chose_trt = mean(trial_data$chose_trt)



arm_times <- data.frame(
  time_secs = c(trial_data$trt_time_secs, trial_data$ctrl_time_secs),
  group = c(rep("Treatment", nrow(trial_data)),
            rep("Control", nrow(trial_data)))
)

arm_times <- (lm(time_secs ~ group, arm_times))

arm_visits <- data.frame(
  visit_num = c(trial_data$trt_visits, trial_data$ctrl_visits),
  group = c(rep("Treatment", nrow(trial_data)),
            rep("Control", nrow(trial_data)))
)

arm_visit_model <- (lm(visit_num ~ group, arm_visits))


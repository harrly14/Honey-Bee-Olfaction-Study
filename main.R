library(lubridate)
library(dplyr)
library(ggplot2)

trial_data <- read.csv("trial_data.csv") 

exclusion_rate <- mean(trial_data$is_excluded)

trial_data <- trial_data[order(trial_data$trial_start),] %>% 
  tibble::rowid_to_column("trial_ID") %>%
  filter(!is_excluded)

nrow(trial_data)

trial_data <-trial_data %>%
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

nrow(trial_data)

prop_chose_trt = mean(trial_data$chose_trt)

preference_binomial_test <- binom.test(x = sum(trial_data$chose_trt == TRUE), 
           n = nrow(trial_data), 
           p = 0.5,
           alternative = "greater", 
           conf.level = 0.95)

# 1. wrangle collection data and batch data. i need location, flower from collection. i need % Na error and time since Na added from batch data
# 2. add to mutate: parent_batch_id (use regex to get just the digit part)
# note for 3 and 4, relevant vars to join on must be chars not factors (convert to chars just before joins then back to factors)
# 3. add to pipeline: left join (trial data, collection data, by = c("cage_id", "parent_batch_id"))
# 4. add to pipeline: left join (trial data, collection data, by = "parent_batch_id") or use by = "x" = "y" if the cols are called different things in the two datasets
# 5. heres my model
preference_model <- glmer(
  chose_trt ~ starvation_duration_hrs + temp + pct_sodium_error + 
    hrs_since_sodium_added + time_of_day_hrs + cage_id + num_bees +
    (1 | batch_id) + (1 | location) + (1 | flower),  # or flower as fixed, pending your answer above
  data = trial_data, family = binomial, na.action = "na.fail"
)
summary(preference_model) # look for convergence warnings
check_collinearity() #idk how to use this, i just know i should

# 6. dredge (SHOILD ONLY DONE AFTER I GET ALL DATA) using the dredge function then refit the model using the dredge results using: 
  # 1. dredge results <- dredge(...)
  # 2. best model <- get.models(dredge results, 1)[[1]]
  # use saveRDS() on the good model
# 7. things i may want to store from the model: 
  # fixed effect coefficients, in log-odds: fixef(sodium_model)
  # same, converted to odds ratios (poster-friendly): exp(fixef(sodium_model))
  # confidence intervals on those odds ratios: exp(confint(sodium_model, method = "Wald"))
  # full coefficient table: estimate, SE, z, p - this is your reportable stats table: coefs <- summary(sodium_model)$coefficients
  # random effect variances - tells you how much batch/location/flower actually mattered: VarCorr(sodium_model)
  # the specific "your headline result" number: is baseline preference different from 50/50? intercept row of `coefs` - its p-value tests exactly this: coefs["(Intercept)", ]
  # model fit quality, for comparing candidate models later (this is what dredge ranks by): AIC(sodium_model)
  # individual predicted probabilities per bee, for building your predicted-probability plot: predict(sodium_model, type = "response")
# 8. plot model results with a simple bar/point chart with CI showing overall proportion choosing treatment vs. 0.5 chance line to start 
# 9. predict(model, type = "response") or the ggeffects/sjPlot package's ggpredict() 
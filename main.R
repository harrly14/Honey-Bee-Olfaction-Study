library(lubridate)
library(dplyr)
library(ggplot2)
library(lme4)
library(performance)
library(MuMIn)
library(glmmTMB)
library(DHARMa)
library(tidyr)
library(showtext)

trial_data <- read.csv("trial_data.csv")
collection_data <- read.csv("collection_data.csv")
batch_data <- read.csv("batch_data.csv")

exclusion_rate <- mean(trial_data$is_excluded)
trial_data <- filter(trial_data, !is_excluded)

# ========================== data processing ================================

nrow(trial_data)

trial_data <- trial_data |>
  filter(!is.na(start_date)) |>
  mutate(
    across(c(species, cage_id, batch_id, trt_arm, first_choice), as.factor),

    starving_start = as.POSIXct(
      paste(start_date, starving_since),
      format = "%m/%d/%Y %I:%M:%S %p"
      ),
    trial_start = as.POSIXct(
      paste(start_date, start_time),
      format = "%m/%d/%Y %I:%M:%S %p"
      ),

    # i'm keeping start date as a factor for possible use in the model later
    start_date = as.factor(start_date),

    latency = period_to_seconds(ms(latency)),
    left_time = period_to_seconds(ms(left_time)),
    right_time = period_to_seconds(ms(right_time)),

    starvation_duration_mins = as.numeric(difftime(trial_start,
                                                   starving_start,
                                                   units = "mins")),
    chose_trt = as.numeric(first_choice == trt_arm),

  )
trial_data <- trial_data |>
  mutate(
    # replace L/R columns with trt/ctrl
    trt_visits = ifelse(trt_arm == "L", left_visits, right_visits),
    ctrl_visits = ifelse(trt_arm == "R", left_visits, right_visits),
    trt_time_secs = ifelse(trt_arm == "L", left_time, right_time),
    ctrl_time_secs = ifelse(trt_arm == "R", left_time, right_time)
  ) |>
  # remove replaced columns
  select(-starving_since, -start_time, -(left_visits:right_time))
trial_data <- trial_data |>
  mutate(
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
        sodium_add_datetime = as.POSIXct(
          paste(sodium_add_date, sodium_add_time),
          format = "%m/%d/%Y %I:%M:%S %p"
        )
      ) |>
      select(pct_sodium_error, batch_id_char, sodium_add_datetime),
    by = "batch_id_char"
  )
trial_data <- trial_data |>
  mutate(
    time_since_sodium_added = as.numeric(difftime(trial_start,
                                                  sodium_add_datetime,
                                                  units = "mins")),
    time_of_day_hrs = hour(trial_start) + minute(trial_start) / 60,
    ) |>
  arrange(trial_start) |>
  tibble::rowid_to_column("trial_ID") |>
  mutate(
    starvation_duration_mins_z = as.numeric(scale(starvation_duration_mins)),
    temp_z = as.numeric(scale(temp)),
    pct_sodium_error_z = as.numeric(scale(pct_sodium_error)),
    time_since_sodium_added_z = as.numeric(scale(time_since_sodium_added)),
    time_of_day_hrs_z = as.numeric(scale(time_of_day_hrs)),
    prop_trt_time_secs = trt_time_secs / (trt_time_secs + ctrl_time_secs),
    # adjusted for beta regression
    adj_prop_trt_time_secs = (prop_trt_time_secs * (nrow(trial_data) - 1) + 0.5) / nrow(trial_data)
  )

nrow(trial_data) # should match pre-join row count

# ========================== stat tests ====================================

prop_chose_trt <- mean(trial_data$chose_trt)

preference_binomial_test <- binom.test(x = sum(trial_data$chose_trt == TRUE),
                                       n = nrow(trial_data),
                                       p = 0.5,
                                       alternative = "greater",
                                       conf.level = 0.95
                                       )

time_in_arm_wilcox_test <- wilcox.test(trial_data$trt_time_secs,
                                       trial_data$ctrl_time_secs,
                                       paired = TRUE
                                       )

# ============================== models =====================================

pref_glmm_full <- glmer(
  chose_trt ~ trt_arm + start_date + time_of_day_hrs_z +
    (1 | batch_id) + (1 | location),
  data = trial_data,
  family = binomial,
  na.action = "na.fail"
)
pref_glmm_best <- glmer(
  chose_trt ~ 1 +
    (1 | batch_id) + (1 | location),
  data = trial_data,
  family = binomial,
  na.action = "na.fail"
)

sim_pref <- simulateResiduals(pref_glmm_best)
plot(sim_pref)
testUniformity(sim_pref)
testDispersion(sim_pref)
testOutliers(sim_pref)
plotResiduals(sim_pref, form = trial_data$batch_id)
plotResiduals(sim_pref, form = as.factor(trial_data$location))


time_glmm_full <- glmmTMB(
  adj_prop_trt_time_secs ~ trt_arm + start_date + time_of_day_hrs_z +
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

sim_time <- simulateResiduals(time_glmm_best)
plot(sim_time)
testUniformity(sim_time)
testDispersion(sim_time)
testOutliers(sim_time)
plotResiduals(sim_time, form = trial_data$batch_id)
plotResiduals(sim_time, form = as.factor(trial_data$location))


visits_glmm_full <- glmer(cbind(trt_visits, ctrl_visits) ~ trt_arm +
                            start_date + time_of_day_hrs_z +
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

sim_visits <- simulateResiduals(visits_glmm_best, n = 1000)
plot(sim_visits)
testUniformity(sim_visits)
testDispersion(sim_visits)
testOutliers(sim_visits)
testZeroInflation(sim_visits)
plotResiduals(sim_visits, form = trial_data$batch_id)
plotResiduals(sim_visits, form = as.factor(trial_data$location))

# =============================== plots =====================================

font_add_google("Open Sans", "opensans")
showtext_auto()
showtext_opts(dpi = 600)

choice_plot <- ggplot(
  trial_data,
  aes(x = "", y = chose_trt, color = factor(chose_trt))
) +
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "red") +
  geom_jitter(
    aes(shape = factor(chose_trt)),
    width = 0.3,
    height = 0.0,
    alpha = 0.6,
    size = 4
  ) +
  stat_summary(
    fun.data = "mean_cl_boot",
    color = "#2C5F8A",
    linewidth = 0.8
  ) +
  scale_color_manual(values = c("0" = "#D59C55", "1" = "#558ed5")) +
  labs(
    x = NULL,
    y = "Proportion Choosing Treatment"
  ) +
  theme_bw(base_size = 20, base_family = "opensans") +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )

plot_data <- trial_data |>
  select(trial_ID, trt_time_secs, ctrl_time_secs, trt_visits, ctrl_visits) |>
  pivot_longer(
    cols = -trial_ID,
    names_to = c("arm", "metric"),
    names_pattern = "^(trt|ctrl)_(.*)$",
    values_to = "value"
  )

base_plot <- ggplot(
                    plot_data,
                    aes(x = arm, y = value, color = arm, fill = arm)
                    ) +
  geom_boxplot(alpha = 0.5, width = 0.7, linewidth = 0.5) +
  geom_point(
             aes(shape = arm), 
             position = position_jitter(width = 0.05, height = 0), 
             alpha = 0.5
             ) +
  geom_line(aes(group = trial_ID), alpha = 0.2, color = "grey") +
  labs(x = NULL, color = "Arm", fill = "Arm") +
  scale_color_manual(
    name = "Arm",
    values = c("ctrl" = "#D59C55", "trt" = "#558ed5"),
    labels = c("Control", "Treatment")
  ) +
  scale_fill_manual(
    name = "Arm",
    values = c("ctrl" = "#D59C55", "trt" = "#558ed5"),
    labels = c("Control", "Treatment")
  ) +
  scale_shape_manual(
    name = "Arm",
    values = c("ctrl" = 16, "trt" = 17),
    labels = c("Control", "Treatment")
  ) +
  scale_x_discrete(labels = c("ctrl" = "Control", "trt" = "Treatment")) +
  theme_bw(base_size = 20, base_family = "opensans") +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )

time_plot <-
  base_plot +
  filter(plot_data, metric == "time_secs") +
  labs(y = "Time Spent (s)")

visit_plot <-
  base_plot +
  filter(plot_data, metric == "visits") +
  labs(y = "Number of Visits")

ggsave(
  "choice_plot.png",
  choice_plot,
  width = 9,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "transparent"
)
ggsave(
  "time_plot.png",
  time_plot,
  width = 6.75,
  height = 5.25,
  units = "in",
  dpi = 600,
  bg = "transparent"
)
ggsave(
  "visit_plot.png",
  visit_plot,
  width = 6.75,
  height = 5.25,
  units = "in",
  dpi = 600,
  bg = "transparent"
)

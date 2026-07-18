# ==============================================================================
# Produces one figure for each of the three key metrics:
#   - first arm choice (binary)
#   - proportion of time in treatment arm
#   - visit count to treatment and control arms
# Requires 01_clean_data.R to have been run.
# Outputs plots/choice_plot.png, plots/time_plot.png, plots/visit_plot.png
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(showtext)

trial_data <- readRDS(here::here("data", "processed", "trial_data_clean.rds"))

# setting seed for reproducability of stat(summary(fun.data = "mean_cl_boot"))
set.seed(20260718)

arm_colors <- c(ctrl = "#D59C55", trt = "#558ed5")
arm_labels <- c(ctrl = "Control", trt = "Salty")

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
    height = 0,
    alpha = 0.6,
    size = 4
  ) +
  stat_summary(
    fun.data = "mean_cl_boot",
    color = "#2C5F8A",
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = c("0" = arm_colors[["ctrl"]],
               "1" = arm_colors[["trt"]]
               )
    ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 35, base_family = "opensans") +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA),
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
  scale_color_manual(
    name = "Arm",
    values = arm_colors,
    labels = arm_labels
  ) +
  scale_fill_manual(
    name = "Arm",
    values = arm_colors,
    labels = arm_labels
  ) +
  scale_shape_manual(
    name = "Arm",
    values = c("ctrl" = 16, "trt" = 17),
    labels = c("Control", "Salty")
  ) +
  scale_x_discrete(labels = arm_labels) +
  labs(x = NULL) +
  theme_classic(base_size = 35, base_family = "opensans") +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none"
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
  here::here("plots", "choice_plot.png"),
  choice_plot,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "transparent"
)
ggsave(
  here::here("plots", "time_plot.png"),
  time_plot,
  height = 8,
  width = 7,
  units = "in",
  dpi = 600,
  bg = "transparent"
)
ggsave(
  here::here("plots", "visits_plot.png"),
  visit_plot,
  height = 8,
  width = 7,
  units = "in",
  dpi = 600,
  bg = "transparent"
)

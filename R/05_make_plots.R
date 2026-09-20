# ==============================================================================
# Produces Figure 2 (Behavioral Metrics) combining three panels:
#   A: first arm choice (binary)
#   B: proportion of time in treatment arm
#   C: visit count to treatment and control arms
# Requires 01_clean_data.R and 02_fit_models.R to have been run.
# Outputs plots/Figure2_Behavioral_Metrics.png
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(showtext)
library(patchwork)
library(DHARMa)

model_results <- readRDS(here::here("results", "model_results.rds"))

sim_choice <- simulateResiduals(model_results$choice$best)
sim_time   <- simulateResiduals(model_results$time$best)
sim_visits <- simulateResiduals(model_results$visits$best)

# Output Figure 1
png(here::here("plots", "Figure1_DHARMa_Checks.png"), 
    width = 10, height = 12, units = "in", res = 600)

par(mfrow = c(3, 2), mar = c(4, 4, 6, 2))

plotQQunif(sim_choice)
mtext(
  "A: First Choice Preference",
  side = 3,
  line = 4.5,
  adj = 0,
  cex = 1.2,
  font = 2
)
plotResiduals(sim_choice)

plotQQunif(sim_time)
mtext(
  "B: Proportion of Time",
  side = 3,
  line = 4.5,
  adj = 0,
  cex = 1.2,
  font = 2
)
plotResiduals(sim_time)

plotQQunif(sim_visits)
mtext(
  "C: Number of Visits",
  side = 3,
  line = 4.5,
  adj = 0,
  cex = 1.2,
  font = 2
)
plotResiduals(sim_visits)

dev.off()

trial_data <- readRDS(here::here("data", "processed", "trial_data_clean.rds")) |>
  # Ensure trial_ID is a factor so geom_line connects paired observations correctly
  mutate(trial_ID = as.factor(trial_ID)) 

# setting seed for reproducibility of stat(summary(fun.data = "mean_cl_boot"))
set.seed(20260718)

arm_colors <- c(ctrl = "gray50", trt = "black")
arm_labels <- c(ctrl = "Control", trt = "Salty")

font_add_google("Open Sans", "opensans")
showtext_auto()
showtext_opts(dpi = 600)

# panel A: Choice Plot
choice_plot <- ggplot(
  trial_data,
  aes(x = "", y = chose_trt, color = factor(chose_trt))
) +
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "black") +
  geom_jitter(
    aes(shape = factor(chose_trt)),
    width = 0.3,
    height = 0,
    alpha = 0.6,
    size = 4
  ) +
  stat_summary(
    fun.data = "mean_cl_boot",
    color = "black",
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = c("0" = arm_colors[["ctrl"]],
               "1" = arm_colors[["trt"]])
  ) +
  scale_shape_manual(
    values = c("0" = 16, "1" = 17) # 16 = circle (control), 17 = triangle (salty)
  ) +
  labs(x = NULL, y = "Probability of Choosing Salty First") +
  theme_classic(base_size = 35, base_family = "opensans") +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA),
  )

# Reshape data for continuous metrics
plot_data <- trial_data |>
  select(trial_ID, trt_time_secs, ctrl_time_secs, trt_visits, ctrl_visits) |>
  pivot_longer(
    cols = -trial_ID,
    names_to = c("arm", "metric"),
    names_pattern = "^(trt|ctrl)_(.*)$",
    values_to = "value"
  )

# Base template for panels B and C
base_plot <- ggplot(
  plot_data,
  aes(x = arm, y = value, color = arm, fill = arm)
) +
  geom_boxplot(
    alpha = 0.5,
    width = 0.7,
    linewidth = 0.5,
    outlier.shape = NA) +
  geom_point(
    aes(shape = arm), 
    position = position_jitter(width = 0.05, height = 0),
    alpha = 0.6,
    size = 3
  ) +
  geom_line(aes(group = trial_ID), 
            alpha = 0.2, 
            color = "gray30",
            linewidth = 0.5
  ) +
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

# panel B: Time Plot
time_plot <- base_plot %+% 
  filter(plot_data, metric == "time_secs") +
  labs(y = "Time Spent (s)")

# panel C: Visit Plot
visit_plot <- base_plot %+% 
  filter(plot_data, metric == "visits") +
  labs(y = "Number of Visits")

figure_2 <- choice_plot + time_plot + visit_plot + 
  plot_annotation(tag_levels = 'A') + 
  plot_layout(ncol = 3, widths = c(1, 1.2, 1.2))

ggsave(
  here::here("plots", "Figure2_Behavioral_Metrics.png"),
  figure_2,
  height = 8,
  width = 20,
  units = "in",
  dpi = 600,
  bg = "white"
)
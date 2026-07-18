# ==============================================================================
# DHARMa residual diagnostics for the best model for the three metrics.
# Meant to be run interactively while reviewing plots, so it is not sourced
# automatically in main.R
# ==============================================================================

library(DHARMa)

model_results <- readRDS(here::here("results", "model_results.rds"))
trial_data <- readRDS(here::here("data", "processed", "trial_data_clean.rds"))

# ------------------------------------------------------------------------------
sim_choice <- simulateResiduals(model_results$choice$best_model)
plot(sim_choice)
testUniformity(sim_choice)
testDispersion(sim_choice)
testOutliers(sim_choice)
plotResiduals(sim_choice, form = trial_data$batch_id)
plotResiduals(sim_choice, form = as.factor(trial_data$location))

# ------------------------------------------------------------------------------
sim_time <- simulateResiduals(model_results$time$best_model)
plot(sim_time)
testUniformity(sim_time)
testDispersion(sim_time)
testOutliers(sim_time)
plotResiduals(sim_time, form = trial_data$batch_id)
plotResiduals(sim_time, form = as.factor(trial_data$location))

# ------------------------------------------------------------------------------
sim_visits <- simulateResiduals(model_results$visits$best_model, n = 1000)
plot(sim_visits)
testUniformity(sim_visits)
testDispersion(sim_visits)
testOutliers(sim_visits)
testZeroInflation(sim_visits)
plotResiduals(sim_visits, form = trial_data$batch_id)
plotResiduals(sim_visits, form = as.factor(trial_data$location))
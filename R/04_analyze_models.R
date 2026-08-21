library(dplyr)
library(tidyr)

model_results <- readRDS(here::here("results", "model_results.rds"))

# create summary table here

model_results$choice$selection_table |> pivot_longer()

# get confints (back transformed)
plogis(confint(model_results$choice$best, parm = "beta_", method = "wald"))
plogis(confint(model_results$time$best,   parm = "beta_", method = "wald"))
plogis(confint(model_results$visits$best, parm = "beta_", method = "wald"))

# get p values (?)
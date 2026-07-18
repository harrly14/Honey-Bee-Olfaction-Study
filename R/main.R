# ==============================================================================
# Entry point for the full pipeline
# Produces:
# - clean processed data (csv and rds) located in data/processed
# - model results located in results/
# - all plots specified in "03_make_plots.R" located in plots/

# 03_check_diagnostics.R is not run here intentionally, as those diagnostics
# should be manually reviewed for model fit when necessary
# ==============================================================================

if(!requireNamespace("here", quietly = TRUE)) install.packages("here")

source(here::here("R", "01_clean_data.R"))
source(here::here("R", "02_fit_models.R"))
source(here::here("R", "04_make_plots.R"))

message("Pipeline complete. See data/processed/, results/, and plots/.")
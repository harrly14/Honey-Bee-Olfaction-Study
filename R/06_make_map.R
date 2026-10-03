# ==============================================================================
# Produces Figure S1 showing collection locations of honey bees.
#
# Requires 01_clean_data.R to have been run.
# Uses data/processed/trial_data_clean.rds.
#
# Outputs plots/FigureS1_Collection_Locations.png
# ==============================================================================

library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(ggspatial)

trial_data <- readRDS(
  here::here("data", "processed", "trial_data_clean.rds")
)

# Locations are stored as "latitude, longitude" in the location column.
collection_sites <- trial_data |>
  filter(!is.na(location)) |>
  distinct(location) |>
  separate(
    location,
    into = c("latitude", "longitude"),
    sep = ",",
    convert = TRUE
  ) |>
  mutate(
    latitude = as.numeric(trimws(latitude)),
    longitude = as.numeric(trimws(longitude))
  ) |>
  filter(
    !is.na(latitude),
    !is.na(longitude)
  )

# Convert collection locations to a spatial object.
collection_sf <- st_as_sf(
  collection_sites,
  coords = c("longitude", "latitude"),
  crs = 4326
)

fig_s1 <- ggplot() +
  annotation_map_tile(
    type = "osm",
    zoomin = 0,
    alpha = 0.8,
    progress = "none"
  ) +
  geom_sf(
    data = collection_sf,
    shape = 21,
    size = 3.2,
    stroke = 0.7,
    color = "white",
    fill = "black"
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    unit_category = "metric",
    text_cex = 0.7,
    pad_x = unit(0.3, "cm"),
    pad_y = unit(0.3, "cm")
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    pad_x = unit(0.3, "cm"),
    pad_y = unit(0.3, "cm"),
    style = north_arrow_minimal(
      text_size = 9
    )
  ) +
  coord_sf(
    expand = TRUE
  ) +
  labs(
    x = NULL,
    y = NULL,
    caption = "\u00a9 OpenStreetMap contributors"
  ) +
  theme_void() +
  theme(
    plot.caption = element_text(
      size = 7,
      color = "grey30",
      hjust = 1,
      margin = margin(t = 3)
    ),
    plot.margin = margin(5, 5, 5, 5)
  )

fig_s1

ggsave(
  here::here("plots", "FigureS1_Collection_Locations.png"),
  plot = fig_s1,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600
)
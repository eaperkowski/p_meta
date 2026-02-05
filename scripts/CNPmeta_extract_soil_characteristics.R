## Libraries
library(tidyverse)
library(raster)
library(soilDB)

## Read in compiled meta-analysis file
df <- read.csv("../data/CNP_data_compiled.csv")

# Create data frame with only field experiments (since field experiments
# only have latitude and longitude)
full_df_field <- df %>%
  filter(experiment_type == "field")

# Create file that includes the latitude and longitude of all unique
# sites in meta-analysis
experiment_summary_field <- distinct(full_df_field, exp, .keep_all = TRUE) %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  dplyr::select(exp, latitude, longitude)

unique(experiment_summary_field$exp)

# Get soilGrids data for nitrogen
soilGrids_nitrogen <- fetchSoilGrids(experiment_summary_field,
                                     loc.names = c("exp", "latitude", "longitude"),
                                     depth_intervals = c("0-5", "5-15"),
                                     variables = "nitrogen")

# Calculate sum of 0-5 and 0-15 horizon
soilGrids_nitrogen_sum <- aqp::horizons(soilGrids_nitrogen) %>%
  group_by(id) %>%
  summarize(background_n_mgkg = sum(nitrogenmean, na.rm = TRUE) * 10 / 2) %>%
  dplyr::select(exp = id, background_n_mgkg)

# Merge soil nitrogen data to compiled dataset
df_nitrogen <- df %>%
  full_join(soilGrids_nitrogen_sum, by = "exp")
write.csv(df_nitrogen, "../data/CNP_data_compiled.csv", row.names = F)

# Read .tif with Olsen-P
olsen_map <- raster::brick("../gridded_products/soil/OlsenP_kgha1_World_Aug2022_ver_COG.tif")
plot(olsen_map)

# Convert lat and long to spatial points
experiment_summary_field_sp <- sp::SpatialPoints(
  cbind(experiment_summary_field$longitude, experiment_summary_field$latitude),
  proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs")  # EPSG:4326
)

# Transform spatial points to align with crs of olsen map
experiment_summary_field_sptrans <- spTransform(experiment_summary_field_sp,
                                                crs(olsen_map))

# Extract Olsen-P point data from each field experiment coordinate
olsenP_site <- raster::extract(olsen_map, experiment_summary_field_sptrans)

olsenP_clean <- experiment_summary_field %>%
  cbind(olsenP_site) %>%
  mutate(olsenP_mgkg = Band_1) %>%
  dplyr::select(exp:longitude, olsenP_mgkg)

# Write final .csv
df_nitrogen %>%
  full_join(olsenP_clean, by = c("exp", "latitude", "longitude")) %>%
  write.csv("../data/CNP_data_compiled.csv", row.names = F)


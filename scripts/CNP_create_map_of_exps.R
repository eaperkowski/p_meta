# Script for visualizing sites included in the CNP meta-analysis. Also
# extracts 1901-2024 climate data from CRU at 0.5deg grid cell and merges
# with file that includes all other data for meta-analysis

#####################################################################
# Load libraries and datasets
#####################################################################

# Load libraries
library(tidyverse)
library(tmap)
library(rnaturalearth)
library(sf)
library(plotbiomes)


# Read raw data file
df <- read.csv("../data/CNP_data_compiled.csv")
df_field <- filter(df, experiment_type == "field")


# Create file that includes the latitude and longitude of all unique
# sites in meta-analysis
experiment_summary_field <- distinct(df_field, exp, .keep_all = TRUE) %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  dplyr::select(exp, latitude, longitude)

exp_plot_prep <- st_as_sf(experiment_summary_field,
                          coords = c("longitude", "latitude"),
                          crs = 4326)

#####################################################################
# Create map of all experiments included in meta-analysis
#####################################################################

# Load world map as sf
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define Robinson projection (EPSG:54030 or proj4string)
robinson_crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Transform both datasets to Robinson
world_robin <- st_transform(world, crs = robinson_crs)
coords_robin <- st_transform(exp_plot_prep, crs = robinson_crs)


ggplot(data = world_robin) +
  geom_sf(fill = "antiquewhite", color = "black", size = 0.2) +
  geom_sf(data = coords_robin, size = 1) +
  theme_minimal(base_size = 12) +
  coord_sf(expand = F,
           lims_method = "box") +
  scale_y_continuous(breaks = seq(from = -90, to = 90, by = 30)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.3),
    #panel.background = element_rect(fill = "lightblue1"),
    legend.position = "right"
  )


CNP_meta_experiment_map <- ggplot() +
  borders(database = "world", colour = "black", fill = "antiquewhite") +
  geom_point(data = experiment_summary_field,
             aes(x = longitude, y = latitude), 
             color = "red", size = 0.5) +
  scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 180, 90)) +
  scale_y_continuous(limits = c(-80, 90), breaks = seq(-60, 90, 30)) +
  labs(x = expression("Longitude ("*degree*")"),
       y = expression("Latitutde ("*degree*")")) +
  theme_classic(base_size = 18) +
  theme(panel.grid = element_blank())

# png("../plots/CNPmeta_site_map.png", width = 3600, height = 2400,
#     res = 600)
CNP_meta_experiment_map 
# dev.off()

#####################################################################
# WorldClim v2.1: Temperature
#####################################################################

# Let's start with temperature. Folder contains multiple .tif files
# that represent mean grid-cell temperature for each month. I reckon
# a rasterStack might be best way to go about this
tavg_files <- list.files("../gridded_products/worldclim/tavg/", 
                         pattern = "\\.tif$",
                         full.names = TRUE)
wc_tavg <- raster::stack(tavg_files)

# Extract average monthly temperature point data from rasterStack
tavg_extracted <-  data.frame(
  raster::extract(
    wc_tavg, experiment_summary_field[, c("longitude", "latitude")]))
names(tavg_extracted) <- gsub("wc2.1_30s_", "", names(tavg_extracted))

# Add experiment information
tavg_extracted2 <- cbind(experiment_summary_field, tavg_extracted)

# Pivot longer to make for easy growing season estimation
tavg_long <- tavg_extracted2 %>%
  pivot_longer(cols = tavg_01:tavg_12,
               names_to = "var", values_to = "value") %>%
  separate(var, into = c("var", "month"), sep = "_")

# Growing season temperature (where degC > 0)
mat <- tavg_long %>%
  group_by(exp, latitude, longitude) %>%
  summarize(mat = mean(value, na.rm = TRUE))

#####################################################################
# WorldClim v2.1: Precipitation
#####################################################################


# Create precipitation rasterBrick
prcp_files <- list.files("../gridded_products/worldclim/prcp/", 
                         pattern = "\\.tif$",
                         full.names = TRUE)
wc_prcp <- raster::stack(prcp_files)

# Extract average monthly precipitation point data from rasterStack
prcp_extracted <-  data.frame(
  raster::extract(
    wc_prcp, experiment_summary_field[, c("longitude", "latitude")]))
names(prcp_extracted) <- gsub("wc2.1_30s_", "", names(prcp_extracted))

# Add experiment information
prcp_extracted2 <- cbind(experiment_summary_field, prcp_extracted)

# Pivot longer to make for easy growing season estimation
prcp_long <- prcp_extracted2 %>%
  pivot_longer(cols = prcp_01:prcp_12,
               names_to = "var", values_to = "value") %>%
  separate(var, into = c("var", "month"), sep = "_")

# Growing season precipitation (where degC > 0)
map <- prcp_long %>%
  group_by(exp, latitude, longitude) %>%
  summarize(map = mean(value, na.rm = TRUE))


#####################################################################
# Combine extracted climate data into single data frame
#####################################################################
climate_norms <- experiment_summary_field %>%
  full_join(mat) %>% full_join(map)

#####################################################################
# Some plots
#####################################################################
# Whittaker plot
png("../plots/CNP_fig1_whittaker_plot.png",
    width = 10, height = 6, units = "in", res = 600)
whittaker_base_plot() +
  geom_point(data = climate_norms,
             aes(x = mat, y = map), alpha = 0.7, size = 2) +
  scale_y_continuous(limits = c(0, 500),
                     breaks = seq(0, 500, 100),
                     labels = seq(0, 5000, 1000)) +
  scale_x_continuous(limits = c(-17, 30),
                     breaks = seq(-15, 30, 15)) +
  labs(x = expression(bold("MAT ("*degree*"C)")),
       y = "MAP (mm)") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.box.background = element_blank())
dev.off()
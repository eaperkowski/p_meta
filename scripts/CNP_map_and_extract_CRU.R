# Script for visualizing sites included in the CNP meta-analysis. Also
# extracts 1901-2024 climate data from CRU at 0.5deg grid cell and merges
# with file that includes all other data for meta-analysis

#####################################################################
# Load libraries and datasets
#####################################################################
library(R.utils)
library(raster)
library(ncdf4)
library(plotbiomes)
library(maps)
library(tidyverse)
library(rnaturalearth)
library(naniar)
library(lubridate)
library(sf)
library(rsofun)

# Read PAR dataset
par <- read.csv("../cru/cru_par_climExtract_growingseason_globe.csv") %>%
  dplyr::select(-X)

# Read data sources (MESI, NutNet, EAP manual compilation)
mesi <- read.csv("../data/mesi_main_manual.csv")
nutnet <- read.csv("../data/nutnet_main.csv")
eap <- read.csv("../data/eap_main.csv")

# Merge data sources into single data frame
full_df <- mesi %>% full_join(nutnet) %>% full_join(eap) %>%
  replace_with_na_all(~.x == "<NA>") %>%
  dplyr::select(-doi, -x_units, -se_c, -se_t) %>%
  mutate(fert = ifelse(npk == "_100", 
                       "n", ifelse(npk == "_010",
                                   "p", ifelse(npk == "_110", "np", NA))))

# Create data frame with only field experiments (since field experiments
# only have latitude and longitude)
full_df_field <- full_df %>%
  filter(experiment_type == "field")

# Create file that includes the latitude and longitude of all unique
# sites in meta-analysis
experiment_summary_field <- distinct(full_df_field, exp, .keep_all = TRUE) %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  dplyr::select(exp, latitude, longitude)

unique(experiment_summary_field$exp)

#####################################################################
# Create map of all experiments included in meta-analysis
#####################################################################
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
# Extract PAR values and merge with climate_long
#####################################################################
# Convert experiment coordinates and CRU PAR coordinates into sf objects
site_coords_sf <- st_as_sf(experiment_summary_field, 
                           coords = c("longitude", "latitude"),
                           crs = 4326)
par_coords_sf <- st_as_sf(par, coords = c("lon", "lat"), crs = 4326)

# Find nearest gridded point for each site
site_coords_sf$nearest_point <- st_nearest_feature(site_coords_sf, par_coords_sf)

# Extract the closest points, merge into experiment_summary_field
experiment_summary_field <- site_coords_sf %>%
  mutate(nearest_PAR = par_coords_sf$par[nearest_point]) %>%
  data.frame() %>%
  full_join(experiment_summary_field) %>%
  dplyr::select(exp, latitude, longitude, par = nearest_PAR)

#####################################################################
# Convert CRU data to RasterBrick, then extract data from coordinates
# of each site
#####################################################################

# Precipitation (mm/month)
pre <- brick("../cru/cru_ts4.09.1901.2024.pre.dat.nc", varname = "pre")
pre_extracted <- data.frame(
  raster::extract(pre, experiment_summary_field[, c("longitude", "latitude")]))
names(pre_extracted)
names(pre_extracted) <- str_c("pre", names(pre_extracted))

# Temperature (degC/month)
temp <- brick("../cru/cru_ts4.09.1901.2024.tmp.dat.nc", varname = "tmp")
temp_extracted <- data.frame(
  raster::extract(temp, experiment_summary_field[, c("longitude", "latitude")]))
names(temp_extracted)
names(temp_extracted) <- str_c("temp", names(temp_extracted))

# PET (mm/day)
pet <- brick("../cru/cru_ts4.09.1901.2024.pet.dat.nc", varname = "pet")
pet_extracted <- data.frame(
  raster::extract(pet, experiment_summary_field[, c("longitude", "latitude")]))
names(pet_extracted)
names(pet_extracted) <- str_c("pet", names(pet_extracted))

#####################################################################
# Combine extracted climate data into single data frame
#####################################################################
climate_all <- cbind(experiment_summary_field,
                     pre_extracted, temp_extracted, 
                     pet_extracted)
names(climate_all)

# Reshape dataframe into long format (i.e. one row per month for each trait
# for each site)
climate_long <- climate_all %>%
  pivot_longer(cols = preX1901.01.16:petX2024.12.16,
               names_to = "variable", values_to = "value") %>%
  separate(variable, into = c("var", "date"), 
           sep = "X", extra = "merge") %>%
  mutate(month = month(ymd(date)),
         day = day(ymd(date)),
         year = year(ymd(date)))

#####################################################################
# Let's calculate summary statistics! First, define the growing
# season, then calculate mean annual growing season temperature, 
# precipitation, pet
#####################################################################

# Define growing season (where monthly temp > 0)
growing_szn <- climate_long %>%
  filter(var == "temp" & year %in% c(1901:2024) & value > 0) %>%
  dplyr::select(exp, latitude, longitude, date)

# Filter climate_long to only include months of the growing 
# season for each experiment
climate_gs_long <- climate_long %>%
  semi_join(growing_szn, by = c("exp", "date")) %>%
  pivot_wider(names_from = var, values_from = value) %>%
  mutate(pet = pet * 30, # scale to mm/month
         ai = pre / pet) %>%
  dplyr::select(exp:date, month, year, pre:ai, par)

# Mean precipitation, potential evapotranspiration, and AI
climate_means <- climate_gs_long %>%
  group_by(exp, latitude, longitude, year, par) %>%
  summarize(annual_gs_pre = sum(pre, na.rm = TRUE),
            annual_gs_pet = sum(pet, na.rm = TRUE),
            annual_gs_temp = mean(temp, na.rm = TRUE),
            annual_gs_ai = annual_gs_pre / annual_gs_pet) %>%
  ungroup() %>%
  group_by(exp, latitude, longitude, par) %>%
  summarize(gs_map = mean(annual_gs_pre, na.rm = TRUE),
            gs_mapet = mean(annual_gs_pet, na.rm = TRUE),
            gs_mat = mean(annual_gs_temp, na.rm = TRUE),
            gs_ai = mean(annual_gs_ai, na.rm = TRUE))

# Merge climate summary with compiled dataset
compiled_df <- full_df %>%
  full_join(climate_means, by = c("exp", "latitude", "longitude")) %>%
  dplyr::select(source:elevation, gs_map:gs_ai, gs_par = par, ecosystem_type:npk, 
                fert, n_c:rep_t)
write.csv(compiled_df, "../data/CNP_data_compiled.csv", row.names = F)

#####################################################################
# Some plots
#####################################################################
# Whittaker plot
png("../plots/CNPmeta_whittaker_plot.png",
    width = 5400, height = 3000, res = 600)
biome_type_plot <- whittaker_base_plot() +
  geom_point(data = climate_means,
             aes(x = gs_mat, y = gs_map / 10)) +
  scale_y_continuous(limits = c(0, 500),
                     breaks = seq(0, 500, 100),
                     labels = seq(0, 5000, 1000)) +
  scale_x_continuous(limits = c(-17, 30),
                     breaks = seq(-15, 30, 15)) +
  labs(x = expression("MAT ("*degree*"C)"),
       y = "MAP (mm)") +
  theme_classic(base_size = 18) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        legend.box.background = element_blank(),
        legend.position = "bottom")
#dev.off()
biome_type_plot

png("../drafts/plots/CNPmeta_sites.png", height = 12, width = 8,
    units = "in", res = 600)
CNP_meta_experiment_map / biome_type_plot
dev.off()


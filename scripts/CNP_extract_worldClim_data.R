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
# WorldClim v2.1: Temperature
#####################################################################

# Let's start with temperature. Folder contains multiple .tif files
# that represent mean grid-cell temperature for each month. I reckon
# a rasterStack might be best way to go about this
tavg_files <- list.files("../cru/world_clim/wc2.1_30s_tavg/", 
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
gs_temp <- tavg_long %>%
  filter(value > 0) %>%
  group_by(exp, latitude, longitude) %>%
  summarize(gs_mat = mean(value, na.rm = TRUE))

#####################################################################
# WorldClim v2.1: Precipitation
#####################################################################

# Create growing season filter (average monthly temperature > 0)
growing_szn <- tavg_long %>% filter(value > 0)

# Create precipitation rasterBrick
prcp_files <- list.files("../cru/world_clim/wc2.1_30s_prec/", 
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
  pivot_longer(cols = prec_01:prec_12,
               names_to = "var", values_to = "value") %>%
  separate(var, into = c("var", "month"), sep = "_") %>%
  mutate(var = "prcp")

# Growing season precipitation (where degC > 0)
gs_prcp <- prcp_long %>%
  semi_join(growing_szn, by = c("exp", "month")) %>%
  group_by(exp, latitude, longitude) %>%
  summarize(gs_map = mean(value, na.rm = TRUE))

#####################################################################
# WorldClim v2.1: Solar radiation (then convert to PAR)
#####################################################################
# Create solar radiation rasterBrick
srad_files <- list.files("../cru/world_clim/wc2.1_30s_srad/", 
                         pattern = "\\.tif$",
                         full.names = TRUE)
wc_srad <- raster::stack(srad_files)

# Extract average monthly precipitation point data from rasterStack
srad_extracted <-  data.frame(
  raster::extract(
    wc_srad, experiment_summary_field[, c("longitude", "latitude")]))
names(srad_extracted) <- gsub("wc2.1_30s_", "", names(srad_extracted))

# Add experiment information
srad_extracted2 <- cbind(experiment_summary_field, srad_extracted)

# Pivot longer to make for easy growing season estimation
srad_long <- srad_extracted2 %>%
  pivot_longer(cols = srad_01:srad_12,
               names_to = "var", values_to = "value") %>%
  separate(var, into = c("var", "month"), sep = "_") %>%
  mutate(var = "par",
         value = value * 1000 / 86400 * 2.1, # convert to umol/m2/s
         value = value * 2) # crude estimator for daytime par (assuming daylength = 0.5 days)

# Growing season precipitation (where degC > 0)
gs_srad <- srad_long %>%
  semi_join(growing_szn, by = c("exp", "month")) %>%
  group_by(exp, latitude, longitude) %>%
  summarize(gs_par = mean(value, na.rm = TRUE))
  
#####################################################################
# WorldClim v2.1: Aridity and elevation
#####################################################################

#################
# Aridity index
#################
# Create raster
ai_index <- raster("../cru/world_clim/ai_v3_yr.tif")

# Extract point data
ai_extracted <- data.frame(
  raster::extract(
    ai_index, experiment_summary_field[, c("longitude", "latitude")]))

# Add experiment information
ai_extracted2 <- cbind(experiment_summary_field, ai_extracted)
names(ai_extracted2)[4] <- "ai"

# Correct AI units
ai_extracted2$ai <- ai_extracted2$ai / 10000

#################
# Elevation
#################
# Create raster
z_index <- raster("../cru/world_clim/wc2.1_30s_elev.tif")

# Extract point data
z_extracted <- raster::extract(
  z_index, experiment_summary_field[, c("longitude", "latitude")])

# Add experiment information
z_extracted2 <- cbind(experiment_summary_field, z_extracted)
names(z_extracted2)[4] <- "z"


#####################################################################
# MERGE WORLDCLIM DATA TOGETHER & INTO COMPILED DATASET
#####################################################################
worldClim_sites <- z_extracted2 %>%
  full_join(gs_temp) %>%
  full_join(gs_prcp) %>%
  full_join(gs_srad) %>%
  full_join(ai_extracted2) %>%
  filter(exp != "hasselt_fnp" & exp != "hammersmith_fnp")

# Merge climate summary with compiled dataset
compiled_df <- full_df %>%
  full_join(worldClim_sites, by = c("exp", "latitude", "longitude")) %>%
  dplyr::select(source:elevation, z:ai, ecosystem_type:npk, 
                fert, n_c:rep_t)
write.csv(compiled_df, "../data/CNP_data_compiled.csv", row.names = F)


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
# WorldClim v2.1: Solar radiation
#####################################################################


# Create precipitation rasterBrick
srad_files <- list.files("../cru/world_clim/wc2.1_30s_srad/", 
                         pattern = "\\.tif$",
                         full.names = TRUE)
wc_srad <- raster::stack(srad_files)

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
  












ai_index <- brick("../cru/world_clim/ai_v3_yr.tif")
ai_extracted <- data.frame(
  raster::extract(ai_index, experiment_summary_field[, c("longitude", "latitude")]))


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


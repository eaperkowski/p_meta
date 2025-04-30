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

# Vapor pressure (HPA/month)
vp <- brick("../cru/cru_ts4.09.1901.2024.vap.dat.nc", varname = "vap")
vp_extracted <- data.frame(
  raster::extract(vp, experiment_summary_field[, c("longitude", "latitude")]))
names(vp_extracted)
names(vp_extracted) <- str_c("vp", names(vp_extracted))

#####################################################################
# Combine extracted climate data into single data frame
#####################################################################
climate_all <- cbind(experiment_summary_field,
                     pre_extracted, temp_extracted, 
                     pet_extracted, vp_extracted)
names(climate_all)

# Reshape dataframe into long format (i.e. one row per month for each trait
# for each site)
climate_long <- climate_all %>%
  pivot_longer(cols = preX1901.01.16:vpX2024.12.16,
               names_to = "variable", values_to = "value") %>%
  separate(variable, into = c("var", "date"), 
           sep = "X", extra = "merge") %>%
  mutate(year = year(ymd(date)))

#####################################################################
# Let's calculate summary statistics! First, calculate total
# precipitation and PET, mean temperature, mean vapor pressure 
# every year, then calculate mean annual values of each variable
#####################################################################

# Mean annual precipitation
map <- climate_long %>%
  filter(var == "pre" & year %in% c(1901:2024)) %>%
  group_by(exp, latitude, longitude, year) %>%
  summarize(annual_precip = sum(value)) %>%
  ungroup(year) %>%
  summarize(map = mean(annual_precip))

# Mean annual potential evapotranspiration
mapet <- climate_long %>%
  filter(var == "pet" & year %in% c(1901:2024)) %>%
  group_by(exp, latitude, longitude, year) %>%
  summarize(annual_pet = sum(value * 30)) %>% # *30 to scale to mm/month
  ungroup(year) %>%
  summarize(mapet = mean(annual_pet))

# Mean annual temperature
mat <- climate_long %>%
  filter(var == "temp" & year %in% c(1901:2024)) %>%
  group_by(exp, latitude, longitude, year) %>%
  summarize(annual_temp = mean(value)) %>%
  ungroup(year) %>%
  summarize(mat = mean(annual_temp))

# Mean annual vapor pressure
mavp <- climate_long %>%
  filter(var == "vp" & year %in% c(1901:2024)) %>%
  group_by(exp, latitude, longitude, year) %>%
  summarize(annual_vaporPressure = mean(value)) %>%
  ungroup(year) %>%
  summarize(mavp = mean(annual_vaporPressure))

# Merge summary statistics and calculate aridity index (P/PET)
complete_climate_summary <- map %>%
  full_join(mat) %>%
  full_join(mapet) %>%
  full_join(mavp) %>%
  mutate(ai = map / mapet)

# Merge climate summary with compiled dataset
compiled_df <- full_df %>%
  full_join(complete_climate_summary, by = c("exp", "latitude", "longitude")) %>%
  dplyr::select(source:elevation, map:ai, ecosystem_type:npk, 
                fert, n_c:rep_t)
# write.csv(compiled_df, "../data/CNP_data_compiled.csv", row.names = F)

#####################################################################
# Some plots
#####################################################################
# Whittaker plot
# png("../plots/CNPmeta_whittaker_plot.png",
#     width = 5400, height = 3000, res = 600)
biome_type_plot <- whittaker_base_plot() +
  geom_point(data = complete_climate_summary,
             aes(x = mat, y = map / 10)) +
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


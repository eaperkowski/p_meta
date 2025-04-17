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

# Create experiment metadata summary
experiment_summary <- full_df %>%
  dplyr::select(citation, exp:experiment_type) %>%
  distinct(citation, exp, .keep_all = TRUE)

# Reproject points to land on Robinson projection
experiment_summary_map <- experiment_summary %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)


map(database = "world")
points(x = experiment_summary$longitude, 
       y = experiment_summary$latitude,
       col = "red")

# Create map of all experiments included in meta-analysis
ggplot() +
  geom_sf(data = world_sf, fill = "antiquewhite", inherit.aes = FALSE) +
  geom_sf(data = experiment_summary_map, 
          color = "red", size = 0.5, inherit.aes = FALSE) +
  coord_sf(crs = "+proj=robin", expand = F) +
  scale_x_continuous(limits = c(-180, 180)) +
  theme_bw(base_size = 18)

  
  
  borders(database = "world", colour = "black") +
  geom_point(data = experiment_summary,
             aes(x = longitude, y = latitude), color = "red", size = 0.5) +
  scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 180, 90)) +
  scale_y_continuous(limits = c(-80, 90), breaks = seq(-60, 90, 30)) +
  labs(x = expression("Longitude ("*degree*")"),
       y = expression("Latitutde ("*degree*")")) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank())

png("../plots/CNPmeta_site_map.png", width = 3600, height = 2400,
    res = 600)
CNP_meta_experiment_map
dev.off()

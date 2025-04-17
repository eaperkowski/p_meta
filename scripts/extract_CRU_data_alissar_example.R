############## packages to extract climate data from climate CRU files
library(R.utils)
library(raster)
library(ncdf4)


######################### Extract climate variables: precipitations, potential evapotranspiration, temperature, VPD
########################### CRU version 4.09 from 1901 to 2024: https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.09/cruts.2503051245.v4.09/

######## load the file containing the coordinates (lat and lon) of all NutNet sites ##
sites <- read.csv("../data/sites_coord.csv")
names(sites)

############ Open precipitation file downloaded from CRU: 1901-2024#########
getwd()
setwd("C:/Users/Alissar/NutNet/CWM_NutNet/data")

nc.pre <- nc_open("cru_ts4.09.1901.2024.pre.dat.nc")
print(nc.pre)
############ Open tmp file downloaded from CRU: 1901-2024#########
nc.tmp <- nc_open("cru_ts4.09.1901.2024.tmp.dat.nc")
print(nc.tmp)

############ Open pet file downloaded from CRU: 1901-2024#########
nc.pet <- nc_open("cru_ts4.09.1901.2024.pet.dat.nc")
print(nc.pet)

############ Open pet file downloaded from CRU: 1901-2024#########
nc.vpd <- nc_open("cru_ts4.09.1901.2024.vap.dat.nc")
print(nc.vpd)

## load the data as rasterBrick object ###############

############# precipitation mm/month ##############
pre <- brick("C:/Users/Alissar/NutNet/CWM_NutNet/data/cru_ts4.09.1901.2024.pre.dat.nc", varname="pre")
pre.sites <- data.frame(extract(pre, sites[, c("longitude", "latitude")])) ## extract precip for each site
names(pre.sites)
colnames(pre.sites) <- paste0("pre", colnames(pre.sites)) ## add precip 
#####################################################################################

############# temperature celsius##############
tmp <- brick("C:/Users/Alissar/NutNet/CWM_NutNet/data/cru_ts4.09.1901.2024.tmp.dat.nc", varname="tmp")
tmp.sites <- data.frame(extract(tmp, sites[, c("longitude", "latitude")])) ## extract tmp for each site
names(tmp.sites)
colnames(tmp.sites) <- paste0("tmp", colnames(tmp.sites))
######################################################################################

######################### potential evapotranspiration mm/day #####################
pet <- brick("C:/Users/Alissar/NutNet/CWM_NutNet/data/cru_ts4.09.1901.2024.pet.dat.nc", varname="pet")
pet.sites <- data.frame(extract(pet, sites[, c("longitude", "latitude")])) ## extract pet for each site
names(pet.sites)
colnames(pet.sites) = paste0("pet", colnames(pet.sites))
########################################################################

######################## vapor pressure in hpa ###########################################
vpd <- brick("C:/Users/Alissar/NutNet/CWM_NutNet/data/cru_ts4.09.1901.2024.vap.dat.nc", varname="vap")
vpd.sites <- data.frame(extract(vpd, sites[, c("longitude", "latitude")])) ## extract vpd for each site
names(vpd.sites)
colnames (vpd.sites)= paste0("vpd", colnames(vpd.sites))
########################################################################################

################ combine all variables ############################################
climate_sites <- cbind(sites, pre.sites, tmp.sites, pet.sites, vpd.sites ) ## combine wiht the original to conserve the pther information
names(climate_sites)
write.csv(climate_sites, "../output/climate_sites.csv")

############## select all years #######################
selected_years <- 1901:2024
## extract the years for each climatic variable 
pattern <- paste0("(", paste(selected_years, collapse = "|"), ")\\.")

# Use grep() to find matching column names
selected_cols <- grep(pattern, names(climate_sites), value = TRUE)

# Subset the data frame with the selected columns
select_climate_cols <- climate_sites[selected_cols]
select_climate_sites_all <- cbind(sites,select_climate_cols)
names(select_climate_sites_all)
write.csv(select_climate_sites_all, "../output/select_climate_sites_all.csv")
view(select_climate_sites_all)

############ calculate the sum of precip and pet per year (we have to multiply pet by the number of days per years if we have to use it), calculate the average tmp and vpd ##########
# Convert the data frame from wide to long format (each row for a single site, year, and variable)
df1 = select_climate_sites_all 
id_vars <- c("site_name", "site_code", "continent","country",
             "habitat", "first_year", "final_year","active", "number_experiment_years", "latitude", "longitude", "experiment_type")
# Reshaping the dataframe and extracting the year
df1_long <- df1 %>%
  pivot_longer(cols = -all_of(id_vars), names_to = "Variable", values_to = "Value") %>%
  separate(Variable, into = c("Var", "Date"), sep = "X", extra = "merge") %>%
  mutate(Year = as.character(str_sub(Date, 1, 4))) %>%  # Extract year as a character
  dplyr::select(all_of(id_vars), Year, Var, Value)  # Keep all categorical variables

# View the reshaped data
print(df1_long)
names(df1_long)
df1_long$Var
write.csv(df1_long, "../output/df1_long.csv")

############### calculate temperature and vpd mean, precip and pet sum per year, and average across years ##########

################# annual mean tmp + average across years #############################
df1_tmp_summary <- df1_long %>%  # anuual mean
  filter(Var == "tmp") %>%
  group_by(across(all_of(id_vars)), Year) %>%
  summarise(Mean_tmp = mean(Value, na.rm = TRUE), .groups = "drop")
view(df1_tmp_summary)
write.csv(df1_tmp_summary, "../output/df1_tmp_summary.csv")

df1_tmp_allyears <- df1_tmp_summary %>%  # Mean 1901-2024
  group_by(across(all_of(id_vars))) %>%
  summarise(MAT = mean(Mean_tmp, na.rm = TRUE), .groups = "drop")
view(df1_tmp_allyears)

# Compute yearly mean for VPD
df1_vpd_summary <- df1_long %>%
  filter(Var == "vpd") %>%
  group_by(across(all_of(id_vars)), Year) %>%
  summarise(Mean_VPD = mean(Value, na.rm = TRUE), .groups = "drop")
view(df1_vpd_summary)

df1_vpd_allyears <- df1_vpd_summary %>%  # Mean 1901-2024
  group_by(across(all_of(id_vars))) %>%
  summarise(MVPD = mean(Mean_VPD, na.rm = TRUE), .groups = "drop")
view(df1_vpd_allyears)

# Compute sum for Precipitation
df1_precip_summary <- df1_long %>%
  filter(Var == "pre") %>%
  group_by(across(all_of(id_vars)), Year) %>%
  summarise(Sum_Precip = sum(Value, na.rm = TRUE), .groups = "drop")
view(df1_precip_summary)

df1_precip_allyears <- df1_precip_summary %>%  # Mean 1901-2024
  group_by(across(all_of(id_vars))) %>%
  summarise(MAP = mean(Sum_Precip, na.rm = TRUE), .groups = "drop")
view(df1_precip_allyears)

# Compute sum for PET
df1_pet_summary <- df1_long %>%
  filter(Var == "pet") %>%
  group_by(across(all_of(id_vars)), Year) %>%
  summarise(Sum_PET = sum(Value*30, na.rm = TRUE), .groups = "drop")
view(df1_pet_summary)

df1_pet_allyears <- df1_pet_summary %>%  # Mean 1901-2024
  group_by(across(all_of(id_vars))) %>%
  summarise(MPET = mean(Sum_PET, na.rm = TRUE), .groups = "drop")
view(df1_pet_allyears)


########## merge annual pet and precip to calculate airidty index ################
precip_pet <- merge(df1_precip_summary, df1_pet_summary, 
                    by= c("site_name", "site_code", "continent","country",
                          "habitat", "first_year", "final_year","active", "number_experiment_years", 
                          "latitude", "longitude", "experiment_type", "Year"),
                    all = TRUE)
view(precip_pet)
names(precip_pet)
precip_pet$AI=precip_pet$Sum_Precip/precip_pet$Sum_PET

################## calculate AI average across years for each site ###############
df1_AI_allyears <- precip_pet %>%  # Mean 1901-2024
  group_by(across(all_of(id_vars))) %>%
  summarise(MAI = mean(AI, na.rm = TRUE), .groups = "drop")
view(df1_AI_allyears)


########################## merge means across years ###################################


vpd_tmp_allyears <- merge(df1_tmp_allyears, df1_vpd_allyears, 
                          by= c("site_name", "site_code", "continent","country",
                                "habitat", "first_year", "final_year","active", "number_experiment_years", 
                                "latitude", "longitude", "experiment_type"),
                          all = TRUE)
view(vpd_tmp_allyears)


vpd_tmp_pre_allyears <- merge(vpd_tmp_allyears, df1_precip_allyears, 
                              by= c("site_name", "site_code", "continent","country",
                                    "habitat", "first_year", "final_year","active", "number_experiment_years", 
                                    "latitude", "longitude", "experiment_type"),
                              all = TRUE)
view(vpd_tmp_pre_allyears)

vpd_tmp_pre_pet_allyears <- merge(vpd_tmp_pre_allyears, df1_pet_allyears, 
                                  by= c("site_name", "site_code", "continent","country",
                                        "habitat", "first_year", "final_year","active", "number_experiment_years", 
                                        "latitude", "longitude", "experiment_type"),
                                  all = TRUE)
view(vpd_tmp_pre_pet_allyears)


climate_summary_allyears <- merge(vpd_tmp_pre_pet_allyears, df1_AI_allyears, 
                                  by= c("site_name", "site_code", "continent","country",
                                        "habitat", "first_year", "final_year","active", "number_experiment_years", 
                                        "latitude", "longitude", "experiment_type"),
                                  all = TRUE)
climate_summary_allyears
view(climate_summary_allyears)
write.csv(climate_summary_allyears, "../output/climate_summary_allyearscsv", row.names = FALSE)

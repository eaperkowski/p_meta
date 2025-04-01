# determine Nutnet treatment mean +/- SE

# Libraries
library(tidyverse)

# Load NutNet data from Alissar paper and Firn et al., 2019, filter to
# include only control, P, N, and NP data
firn_data <- read.csv("../../raw_data/Firn_2019_nature.csv") %>%
  filter(trt %in% c("P", "Control", "N", "NP") & Exclose == 0 & year_trt > 1) %>%
  mutate(trt = factor(trt, levels = c("Control", "N", "P", "NP")),
         lma = (1/SLA_v2) * 1000000,
         narea = (leaf_pct_N/100) * lma,
         parea = (leaf_pct_P/100) * lma,
         leaf_np = leaf_pct_N / leaf_pct_P)

head(firn_data)

# How many sites?
unique(firn_data$site_code)


# Calculate site-level treatment means +/- SE
firn_data_summary <- firn_data %>%
  group_by(site_code,  trt, Taxon) %>%
  summarize(
    
    exp_duration = mean(year_trt),
    
    lma_n = sum(!is.na(lma)),
    lma_mean = mean(lma, na.rm = TRUE),
    lma_sd = sd(lma, na.rm = TRUE),
    lma_se = lma_sd / sqrt(lma_n),
    
    nmass_n = sum(!is.na(leaf_pct_N)),
    nmass_mean = mean(leaf_pct_N, na.rm = TRUE),
    nmass_sd = sd(leaf_pct_N, na.rm = TRUE),
    nmass_se = nmass_sd / sqrt(nmass_n),
    
    narea_n = sum(!is.na(narea)),
    narea_mean = mean(narea, na.rm = TRUE),
    narea_sd = sd(narea, na.rm = TRUE),
    narea_se = narea_sd / sqrt(narea_n),
    
    pmass_n = sum(!is.na(leaf_pct_P)),
    pmass_mean = mean(leaf_pct_P, na.rm = TRUE),
    pmass_sd = sd(leaf_pct_P, na.rm = TRUE),
    pmass_se = pmass_sd / sqrt(pmass_n),
    
    parea_n = sum(!is.na(parea)),
    parea_mean = mean(parea, na.rm = TRUE),
    parea_sd = sd(parea, na.rm = TRUE),
    parea_se = parea_sd / sqrt(parea_n),
    
    leafnp_n = sum(!is.na(leaf_np)),
    leafnp_mean = mean(leaf_np, na.rm = TRUE),
    leafnp_sd = sd(leaf_np, na.rm = TRUE),
    leafnp_se = leafnp_sd / sqrt(leafnp_n)
    
  )

# Prep for easy merge into compiled datasheet
firn_data_summary_control <- firn_data_summary %>%
  filter(trt == "Control") %>%
  ungroup(trt) %>%
  select(-trt)
names(firn_data_summary_control)[4:27] <- str_c(names(firn_data_summary_control)[4:27], "_control")

firn_data_summary_treatment <- firn_data_summary %>%
  filter(trt != "Control")
names(firn_data_summary_treatment)[5:28] <- str_c(names(firn_data_summary_treatment)[5:28], "_trt")

# Format into easy merge into compiled datasheet, write to .csv
firn_data_summary_control %>%
  full_join(firn_data_summary_treatment, by = c("site_code", "Taxon", "exp_duration")) %>%
  dplyr::select(site_code:exp_duration, trt,
                
                lma_mean_control, lma_mean_trt, 
                lma_sd_control, lma_sd_trt, 
                lma_se_control, lma_se_trt, 
                lma_n_control, lma_n_trt,
                
                nmass_mean_control, nmass_mean_trt, 
                nmass_sd_control, nmass_sd_trt, 
                nmass_se_control, nmass_se_trt, 
                nmass_n_control, nmass_n_trt,
                
                narea_mean_control, narea_mean_trt, 
                narea_sd_control, narea_sd_trt, 
                narea_se_control, narea_se_trt, 
                narea_n_control, narea_n_trt,
                
                pmass_mean_control, pmass_mean_trt, 
                pmass_sd_control, pmass_sd_trt, 
                pmass_se_control, pmass_se_trt, 
                pmass_n_control, pmass_n_trt,
                
                parea_mean_control, parea_mean_trt, 
                parea_sd_control, parea_sd_trt, 
                parea_se_control, parea_se_trt, 
                parea_n_control, parea_n_trt,
                
                leafnp_mean_control, leafnp_mean_trt, 
                leafnp_sd_control, leafnp_sd_trt, 
                leafnp_se_control, leafnp_se_trt, 
                leafnp_n_control, leafnp_n_trt) %>%
  pivot_longer(cols = lma_mean_control:leafnp_n_trt,
               names_to = c("trait", "stat", "treatment"), 
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = stat:treatment,
              values_from = value) %>%
  mutate(trait = factor(trait, levels = c("lma", "nmass", "narea",
                                          "pmass", "parea", "leafnp")),
         Taxon = gsub(tolower(Taxon), pattern = " ", replacement = "_")) %>%
  group_by(site_code, Taxon, trait) %>%
  filter(n_control > 1 & n_trt > 1) %>%
  filter(n_distinct(trt) == 3) %>%
  arrange(trait, site_code, Taxon, trt) %>%
  write.csv("../../data_summaries/Firn2019_nature_summary_spp.csv", row.names = F)


# Climate data summaries
firn_data %>%
  group_by(site_code) %>%
  summarize(latitude = mean(latitude), 
            longitude = mean(longitude), 
            MAT = mean(MAT), 
            MAP = mean(MAP)) %>%
  write.csv("../data_summaries/Firn2019_nature_climate_summary.csv", row.names = F)


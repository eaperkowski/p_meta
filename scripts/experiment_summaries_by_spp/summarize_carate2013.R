# Reformat Carate et al. (2013) summary table for easy load
# into compiled meta-analysis datasheet

# Library
library(tidyverse)

# Read data
carate_df <- read.csv("../../raw_data/Carate_Tandalla_2013.csv") %>%
  dplyr::select(Site_name, PlotID, trt = PlotTreat, Tree_species,
                sla = specific_leaf_area, biomass_root:leaf_P)

# Convert traits to numeric
carate_df[, 5:10] <- lapply(carate_df[, 5:10], as.numeric)

# Add a few calculations
carate_df <- carate_df %>%
  mutate(sla_m2g = sla / 10000,
         agb = biomass_shoot + biomass_leaf,
         root_shoot = biomass_root / agb,
         rmf = biomass_root / (biomass_root + biomass_shoot + biomass_leaf),
         narea = leaf_N / 1000 / sla_m2g,
         parea = leaf_P / 1000 / sla_m2g,
         leaf_np = leaf_N / leaf_P)
head(carate_df)

# Summarize dataset
carate_data_summary_plot <- carate_df %>%
  group_by(Site_name, PlotID, trt, Tree_species) %>%
  summarize(agb_plot = mean(agb, na.rm = TRUE),
            bgb_plot = mean(biomass_root, na.rm = TRUE),
            rootshoot_plot = mean(root_shoot, na.rm = TRUE),
            rmf_plot = mean(rmf, na.rm = TRUE),
            sla_plot = mean(sla, na.rm = TRUE),
            nmass_plot = mean(leaf_N, na.rm = TRUE),
            narea_plot = mean(narea, na.rm = TRUE),
            pmass_plot = mean(leaf_P, na.rm = TRUE),
            parea_plot = mean(parea, na.rm = TRUE),
            leafnp_plot = mean(leaf_np, na.rm = TRUE))
carate_data_summary_plot

# Data summary by plot mean
carate2013_data_summary <- carate_data_summary_plot %>%
  group_by(Site_name, trt, Tree_species) %>%
  summarize(
    
    agb_n = sum(!is.na(agb_plot)),
    agb_mean = mean(agb_plot, na.rm = TRUE),
    agb_sd = sd(agb_plot, na.rm = TRUE),
    agb_se = agb_sd / sqrt(agb_n),
    
    bgb_n = sum(!is.na(bgb_plot)),
    bgb_mean = mean(bgb_plot, na.rm = TRUE),
    bgb_sd = sd(bgb_plot, na.rm = TRUE),
    bgb_se = bgb_sd / sqrt(bgb_n),
    
    rootshoot_n = sum(!is.na(rootshoot_plot)),
    rootshoot_mean = mean(rootshoot_plot, na.rm = TRUE),
    rootshoot_sd = sd(rootshoot_plot, na.rm = TRUE),
    rootshoot_se = rootshoot_sd / sqrt(rootshoot_n),
    
    rmf_n = sum(!is.na(rmf_plot)),
    rmf_mean = mean(rmf_plot, na.rm = TRUE),
    rmf_sd = sd(rmf_plot, na.rm = TRUE),
    rmf_se = rmf_sd / sqrt(rmf_n),
    
    sla_n = sum(!is.na(nmass_plot)),
    sla_mean = mean(sla_plot, na.rm = TRUE),
    sla_sd = sd(sla_plot, na.rm = TRUE),
    sla_se = sla_sd / sqrt(sla_n),
    
    nmass_n = sum(!is.na(nmass_plot)),
    nmass_mean = mean(nmass_plot, na.rm = TRUE),
    nmass_sd = sd(nmass_plot, na.rm = TRUE),
    nmass_se = nmass_sd / sqrt(nmass_n),
    
    narea_n = sum(!is.na(narea_plot)),
    narea_mean = mean(narea_plot, na.rm = TRUE),
    narea_sd = sd(narea_plot, na.rm = TRUE),
    narea_se = narea_sd / sqrt(narea_n),
    
    pmass_n = sum(!is.na(pmass_plot)),
    pmass_mean = mean(pmass_plot, na.rm = TRUE),
    pmass_sd = sd(pmass_plot, na.rm = TRUE),
    pmass_se = pmass_sd / sqrt(pmass_n),
    
    parea_n = sum(!is.na(parea_plot)),
    parea_mean = mean(parea_plot, na.rm = TRUE),
    parea_sd = sd(parea_plot, na.rm = TRUE),
    parea_se = parea_sd / sqrt(parea_n),
    
    leafnp_n = sum(!is.na(leafnp_plot)),
    leafnp_mean = mean(leafnp_plot, na.rm = TRUE),
    leafnp_sd = sd(leafnp_plot, na.rm = TRUE),
    leafnp_se = leafnp_sd / sqrt(leafnp_n)
    )


# Prep for easy merge into compiled datasheet
carate_summary_control <- carate2013_data_summary %>%
  filter(trt == "con") %>%
  ungroup(trt) %>%
  select(-trt)
names(carate_summary_control)[3:42] <- str_c(names(carate_summary_control)[3:42], "_control")

carate_summary_treatment <- carate2013_data_summary %>%
  filter(trt != "con")
names(carate_summary_treatment)[4:43] <- str_c(names(carate_summary_treatment)[4:43], "_trt")

# Format into easy merge into compiled datasheet, write to .csv
carate_summary_control %>%
  full_join(carate_summary_treatment, by = c("Site_name", "Tree_species")) %>%
  dplyr::select(Site_name, trt, Tree_species,
                
                agb_mean_control, agb_mean_trt, 
                agb_sd_control, agb_sd_trt, 
                agb_se_control, agb_se_trt, 
                agb_n_control, agb_n_trt,
    
                bgb_mean_control, bgb_mean_trt, 
                bgb_sd_control, bgb_sd_trt, 
                bgb_se_control, bgb_se_trt, 
                bgb_n_control, bgb_n_trt,
                
                rootshoot_mean_control, rootshoot_mean_trt, 
                rootshoot_sd_control, rootshoot_sd_trt, 
                rootshoot_se_control, rootshoot_se_trt, 
                rootshoot_n_control, rootshoot_n_trt,
                
                rmf_mean_control, rmf_mean_trt, 
                rmf_sd_control, rmf_sd_trt, 
                rmf_se_control, rmf_se_trt, 
                rmf_n_control, rmf_n_trt,
                
                sla_mean_control, sla_mean_trt, 
                sla_sd_control, sla_sd_trt, 
                sla_se_control, sla_se_trt, 
                sla_n_control, sla_n_trt,
                
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
  pivot_longer(cols = agb_mean_control:leafnp_n_trt,
               names_to = c("trait", "stat", "treatment"), 
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = stat:treatment,
              values_from = value) %>%
  mutate(trait = factor(trait, levels = c("agb", "bgb", "rootshoot", "rmf", 
                                          "sla", "nmass", "narea", "pmass",
                                          "parea", "leafnp")),
         trt = factor(trt, levels = c("N", "P", "NP")),
         Tree_species = gsub(tolower(Tree_species), pattern = " ", replacement = "_")) %>%
  group_by(Site_name, Tree_species, trait) %>%
  filter(n_control > 1 & n_trt > 1) %>%
  filter(n_distinct(trt) == 3) %>%
  arrange(Site_name, trait, Tree_species, trt) %>%
  write.csv("../../data_summaries/species_level/Carate_2013_summarized.csv", row.names = F)


















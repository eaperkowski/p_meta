# Script that explores the effect of N, P, and N+P addition on leaf- 
# and whole-plant functional traits. Meta-analysis uses data compiled
# from full-factorial N*P manipulation experiments. Whole-plant
# measurements are aggregated at the experiment level and leaf-level
# measurements are aggregated at the species level.

# Libraries
library(tidyverse)
library(metafor)
library(ggpubr)
library(forcats)
library(patchwork)
library(naniar) # to resolve NA/<NA> issue
library(orchaRd)

# Read compiled dataset
full_df <- read.csv("../data/CNP_data_compiled.csv") %>%
  replace_with_na_all(~.x == "<NA>") %>%
  filter(treatment == "fnp")

# Create experiment metadata summary
experiment_summary <- full_df %>%
  dplyr::select(citation, exp:experiment_type) %>%
  distinct(citation, exp, .keep_all = TRUE)

# Species identity summary
species_summary <- full_df %>% 
  distinct(species, family, growth_form, growth_duration,
           photo_path, n_fixer, myc_assoc, .keep_all = TRUE) %>%
  dplyr::select(species:myc_assoc) %>%
  slice(-1)

# How many experiments?
unique(experiment_summary$citation)

# How many sites?
unique(filter(experiment_summary, experiment_type == "field")$exp)

# How many field experiments?
unique(filter(experiment_summary, experiment_type == "field")$citation)
unique(filter(experiment_summary, experiment_type == "greenhouse")$citation)
unique(filter(experiment_summary, experiment_type == "chamber")$citation)

# How many species?
distinct(full_df, species)
distinct(full_df, family)

# Add helper function for aggregating lnRRs across experiments
# (determines meta-analytic mean)
source("../helper_fxns/analyse_meta.R")

# Add helper function for calculating interaction effect sizes
# (uses Hedge's d)
source("../helper_fxns/calc_intxn_effSize_meta.R")

#####################################################################
# Set up and carry out N addition meta-analysis
#####################################################################

# Subset `full_df` to include only N addition treatments
nfert_only <- full_df %>% filter(npk == "_100")
head(nfert_only)

# What traits are included?
unique(nfert_only$response)

# Which experiments are included?
unique(nfert_only$citation)

unique(nfert_only$response)

# Select variables
use_response_n <- c("anpp_n", "anpp_p", "anpp", "leaf_n_mass", "leaf_p_mass", "anet",
                    "tbio_gm2", "bnpp", "rootshoot", "bgb", "lma",
                    "lai", "rmf", "tla", "gsw", "agb", "leaf_np", "total_biomass",
                    "asat", "leaf_wue", "leaf_n_area", "leaf_p_area", "vcmax", 
                    "jmax", "jmax_vcmax", "rd", "leaf_pnue", "leaf_ppue",
                    "leaf_structural_p", "leaf_metabolic_p", "leaf_nucleic_p",
                    "leaf_residual_p", "amax", "leaf_pi", "leaf_nre", "leaf_pre", 
                    "fine_root_biomass")

nfert_responses <- nfert_only %>%
  filter(response %in% use_response_n) %>%
  mutate(myvar = response) %>%
  mutate(myvar = ifelse(myvar %in% c("amax", "anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("bgb", "fine_root_biomass"),
                        "bgb", myvar),
         myvar = ifelse(myvar %in% "lai", "tla", myvar))

use_vars_n <- unique(nfert_responses$myvar)

# Calculate log-response ratios using `escalc` in `metafor`. Also, convert SLA 
# to Marea (logr = opposite sign, stdev/error/variance remains the same)

nfert_lnRR <- nfert_responses %>% 
  escalc(measure = "ROM",
         m1i = x_t, sd1i = sd_t, n1i = rep_t, 
         m2i = x_c, sd2i = sd_c, n2i = rep_c, 
         data = ., 
         append = TRUE, 
         var.names = c("logr", "logr_var")) %>%
  mutate(logr_se = sqrt(logr_var) / sqrt(rep_t))  %>%
  mutate(logr = ifelse(myvar == "sla", logr * -1, logr), 
         myvar = ifelse(myvar == "sla", "lma", myvar))

# Determine lnRR across multiple experiments. Takes into account
# experiment variance and uses experiment identity as grouping
# factor for random intercepts
out_n <- purrr::map(as.list(use_vars_n),
                    ~analyse_meta(nfert_lnRR %>%
                                    rename(var = myvar), 
                                  nam_target = .))
names(out_n) <- use_vars_n

# Extract summary statistics for each model
out_n_modl <- lapply(out_n, `[[`, "modl")

out_n_modl_df <- data.frame(
  var = names(out_n_modl),
  nut_add = "n",
  k = sapply(out_n_modl, \(x) x$k),
  estimate = sapply(out_n_modl, \(x) as.numeric(coef(x))[1]),
  SE = sapply(out_n_modl, \(x) x$se[1]),
  zval = sapply(out_n_modl, \(x) x$zval[1]),
  pval = sapply(out_n_modl, \(x) x$pval[1]),
  ci.lb = sapply(out_n_modl, \(x) x$ci.lb[1]),
  ci.ub = sapply(out_n_modl, \(x) x$ci.ub[1]),
  row.names = NULL) |>
  mutate(across(estimate:zval, \(x) round(x, 3)),
         across(ci.lb:ci.ub, \(x) round(x, 3)),
         pval = ifelse(pval < 0.001, "<0.001",
                       round(pval, digits = 3)))



#####################################################################
# Set up and carry out P fertilization meta-analysis
#####################################################################

# Subset `full_df` to include only P addition treatments
pfert_only <- full_df %>% filter(npk == "_010")
head(pfert_only)

# What traits are included?
unique(pfert_only$response)

# Select variables
use_response_p <- c("anpp_n", "anpp_p", "anpp", "leaf_n_mass", "leaf_p_mass", "anet",
                    "tbio_gm2", "bnpp", "rootshoot", "bgb", "lma",
                    "lai", "rmf", "tla", "gsw", "agb", "leaf_np", "total_biomass",
                    "asat", "leaf_wue", "leaf_n_area", "leaf_p_area", "vcmax", 
                    "jmax", "jmax_vcmax", "rd", "leaf_pnue", "leaf_ppue",
                    "leaf_structural_p", "leaf_metabolic_p", "leaf_nucleic_p",
                    "leaf_residual_p", "amax", "leaf_pi", "leaf_nre", "leaf_pre", 
                    "fine_root_biomass")

pfert_responses <- pfert_only %>%
  filter(response %in% use_response_p) %>%
  mutate(myvar = response) %>%
  mutate(myvar = ifelse(myvar %in% c("amax", "anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("bgb", "fine_root_biomass"),
                        "bgb", myvar),
         myvar = ifelse(myvar %in% "lai", "tla", myvar))

use_vars_p <- unique(pfert_responses$myvar)

# Calculate log-response ratios using `escalc` in `metafor`
pfert_lnRR <- pfert_responses %>% 
  escalc(measure = "ROM",
         m1i = x_t, sd1i = sd_t, n1i = rep_t, 
         m2i = x_c, sd2i = sd_c, n2i = rep_c, 
         data = ., 
         append = TRUE, 
         var.names = c("logr", "logr_var")) %>%
  mutate(logr_se = sqrt(logr_var) / sqrt(rep_t))

# Determine lnRR across multiple experiments. Takes into account
# experiment variance and uses experiment identity as grouping
# factor for random intercepts
out_p <- purrr::map(as.list(use_vars_p),
                    ~analyse_meta(pfert_lnRR %>%
                                    rename(var = myvar), nam_target = .))
names(out_p) <- use_vars_p

# Extract summary statistics for each model
out_p_modl <- lapply(out_p, `[[`, "modl")

out_p_modl_df <- data.frame(
  var = names(out_p_modl),
  nut_add = "p",
  k = sapply(out_p_modl, \(x) x$k),
  estimate = sapply(out_p_modl, \(x) as.numeric(coef(x))[1]),
  SE = sapply(out_p_modl, \(x) x$se[1]),
  zval = sapply(out_p_modl, \(x) x$zval[1]),
  pval = sapply(out_p_modl, \(x) x$pval[1]),
  ci.lb = sapply(out_p_modl, \(x) x$ci.lb[1]),
  ci.ub = sapply(out_p_modl, \(x) x$ci.ub[1]),
  row.names = NULL) |>
  mutate(across(estimate:zval, \(x) round(x, 3)),
         across(ci.lb:ci.ub, \(x) round(x, 3)),
         pval = ifelse(pval < 0.001, "<0.001",
                       round(pval, digits = 3)))

#####################################################################
# Set up and carry out N+P fertilization meta-analysis
#####################################################################

# Subset `full_df` to include only N+P addition treatments
npfert_only <- full_df %>% filter(npk == "_110")
head(npfert_only)

# What traits are included?
unique(npfert_only$response)

# Select variables
use_response_np <- c("anpp_n", "anpp_p", "anpp", "leaf_n_mass", "leaf_p_mass", "anet",
                    "tbio_gm2", "bnpp", "rootshoot", "bgb", "lma",
                    "lai", "rmf", "tla", "gsw", "agb", "leaf_np", "total_biomass",
                    "asat", "leaf_wue", "leaf_n_area", "leaf_p_area", "vcmax", 
                    "jmax", "jmax_vcmax", "rd", "leaf_pnue", "leaf_ppue",
                    "leaf_structural_p", "leaf_metabolic_p", "leaf_nucleic_p",
                    "leaf_residual_p", "amax", "leaf_pi", "leaf_nre", "leaf_pre", 
                    "fine_root_biomass")

npfert_responses <- npfert_only %>%
  filter(response %in% use_response_np) %>%
  mutate(myvar = response) %>%
  mutate(myvar = ifelse(myvar %in% c("amax", "anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("bgb", "fine_root_biomass"),
                        "bgb", myvar),
         myvar = ifelse(myvar %in% "lai", "tla", myvar))


use_vars_np <- unique(npfert_responses$myvar)

# Calculate log-response ratios using `escalc` in `metafor`
npfert_lnRR <- npfert_responses %>% 
  escalc(measure = "ROM",
         m1i = x_t, sd1i = sd_t, n1i = rep_t, 
         m2i = x_c, sd2i = sd_c, n2i = rep_c, 
         data = ., 
         append = TRUE, 
         var.names = c("logr", "logr_var")) %>%
  mutate(logr_se = sqrt(logr_var) / sqrt(rep_t))

# Determine lnRR across multiple experiments. Takes into account
# experiment variance and uses experiment identity as grouping
# factor for random intercepts
out_np <- purrr::map(as.list(use_vars_np),
                    ~analyse_meta(npfert_lnRR %>%
                                    rename(var = myvar), nam_target = .))
names(out_np) <- use_vars_np

# Extract summary statistics for each model
out_np_modl <- lapply(out_np, `[[`, "modl")

out_np_modl_df <- data.frame(
  var = names(out_np_modl),
  nut_add = "np",
  k = sapply(out_np_modl, \(x) x$k),
  estimate = sapply(out_np_modl, \(x) as.numeric(coef(x))[1]),
  SE = sapply(out_np_modl, \(x) x$se[1]),
  zval = sapply(out_np_modl, \(x) x$zval[1]),
  pval = sapply(out_np_modl, \(x) x$pval[1]),
  ci.lb = sapply(out_np_modl, \(x) x$ci.lb[1]),
  ci.ub = sapply(out_np_modl, \(x) x$ci.ub[1]),
  row.names = NULL) |>
  mutate(across(estimate:zval, \(x) round(x, 3)),
         across(ci.lb:ci.ub, \(x) round(x, 3)),
         pval = ifelse(pval < 0.001, "<0.001",
                       round(pval, digits = 3)))

##############################################################################
# Some prep work for calculating interaction effect sizes (need
# to turn summary statistics into long format)
##############################################################################

# Create df that only includes control trt, prep for merge to long format
full_df_control <- full_df %>%
  mutate(unique_id = ifelse(
    is.na(species),
    str_c(citation, "XXX", exp, "XXX", response, "XXX", "none"),
    str_c(citation, "XXX", exp, "XXX", response, "XXX", species))) %>%
  group_by(unique_id) %>%
  summarize(x_c = mean(x_c),
            sd_c = mean(sd_c),
            rep_c = mean(rep_c))

# Convert df into long format to make for easy effect size calculation
full_df_trt_long <- full_df %>% 
  mutate(unique_id = ifelse(is.na(species), 
                            str_c(citation, "XXX", exp, "XXX", response, "XXX", "none"),
                            str_c(citation, "XXX", exp, "XXX", response, "XXX", species))) %>%
  dplyr::select(unique_id, fert, x_t, sd_t, rep_t) %>%
  pivot_wider(names_from = fert, values_from = c(x_t, sd_t, rep_t),
              values_fn = mean) %>%
  full_join(full_df_control) %>%
  dplyr::select(unique_id,
                x_c, sd_c, rep_c,
                x_n = x_t_n, sd_n = sd_t_n, rep_n = rep_t_n,
                x_p = x_t_p, sd_p = sd_t_p, rep_p = rep_t_p,
                x_np = x_t_np, sd_np = sd_t_np, rep_np = rep_t_np)


###############################################################################
# Calculate individual, main, and interaction effect sizes
###############################################################################
CNP_effect_sizes <- data.frame(
  full_df_trt_long,
  calc_intxn_effSize_meta(x_a = full_df_trt_long$x_n, 
                          s_a = full_df_trt_long$sd_n, 
                          n_a = full_df_trt_long$rep_n,
                          
                          x_b = full_df_trt_long$x_p, 
                          s_b = full_df_trt_long$sd_p, 
                          n_b = full_df_trt_long$rep_p,
                          
                          x_c = full_df_trt_long$x_c, 
                          s_c = full_df_trt_long$sd_c, 
                          n_c = full_df_trt_long$rep_c,
                          
                          x_ab = full_df_trt_long$x_np, 
                          s_ab = full_df_trt_long$sd_np, 
                          n_ab = full_df_trt_long$rep_np)) %>%
  separate(unique_id, 
           into = c("citation", "exp", "response", "species"), 
           sep = "XXX", remove = TRUE) %>%
  full_join(experiment_summary, by = c("citation", "exp")) %>%
  replace_with_na_all(~.x == "none") %>%
  full_join(species_summary, by = "species") %>%
  mutate(pft = str_c(photo_path, n_fixer, myc_assoc, sep = "_"),
         dAB_se = sqrt(v_ab_int) / sqrt(rep_np)) %>%
  dplyr::select(citation, exp, latitude:experiment_type, 
                species, family:myc_assoc, pft, response,
                
                # Individual effects
                gNi = g_a, vNi = v_a, wNi = w_a,
                gPi = g_b, vPi = v_b, wPi = w_b,
                gNPi = g_ab, vNPi = v_ab, wNPi = w_ab,
                
                # Main effects
                dNi = dA, dvNi = v_a_main, dwNi = w_a_main,
                dPi = dB, dvPi = v_b_main, dwPi = w_b_main,
                
                # Interaction effects
                dNPi = dAB, dvNPi = v_ab_int, dwNPi = w_ab_int,
                dNPi_se = dAB_se)

###############################################################################
# Clean CNP_effect_sizes to only include relevant traits. Also tidy up to
# group similar traits without losing integrity of original trait value in
# `response` column
###############################################################################
CNP_effect_sizes_reduced <- CNP_effect_sizes %>%
  filter(response %in% c("anpp_n", "anpp_p", "anpp", "leaf_n_mass", "leaf_p_mass", "anet",
         "tbio_gm2", "bnpp", "rootshoot", "bgb", "lma",
         "lai", "rmf", "tla", "gsw", "agb", "leaf_np", "total_biomass",
         "asat", "leaf_wue", "leaf_n_area", "leaf_p_area", "vcmax", 
         "jmax", "jmax_vcmax", "rd", "leaf_pnue", "leaf_ppue",
         "leaf_structural_p", "leaf_metabolic_p", "leaf_nucleic_p",
         "leaf_residual_p", "amax", "leaf_pi", "leaf_nre", "leaf_pre", 
         "fine_root_biomass")) %>%
  mutate(response = ifelse(response %in% c("amax", "anet", "asat"),
                        "asat", response),
         response = ifelse(response %in% c("bgb", "fine_root_biomass"),
                        "bgb", response),
         response = ifelse(response %in% "lai", "tla", response))

write.csv(CNP_effect_sizes_reduced, "../data/CNPmeta_logr_results_int.csv", row.names = F)

###############################################################################
# Use helper fxn from `analyse_meta.R` to iterate through traits and determine
# meta-analytic mean of the interaction effect size and associated 95% CI
###############################################################################

# Determine interaction effect size across multiple experiments. 
# Takes into account experiment variance and uses experiment 
# identity as grouping factor for random intercepts
use_vars_int <- unique(CNP_effect_sizes_reduced$response)

out_int <- purrr::map(as.list(use_vars_int),
                      ~analyse_meta_int(CNP_effect_sizes_reduced %>%
                                          rename(var = response), nam_target = .))
names(out_int) <- use_vars_int

# Extract summary statistics for each model
out_int_modl <- lapply(out_int, `[[`, "modl")

out_int_modl_df <- data.frame(
  var = names(out_int_modl),
  nut_add = "n",
  k = sapply(out_int_modl, \(x) x$k),
  estimate = sapply(out_int_modl, \(x) as.numeric(coef(x))[1]),
  SE = sapply(out_int_modl, \(x) x$se[1]),
  zval = sapply(out_int_modl, \(x) x$zval[1]),
  pval = sapply(out_int_modl, \(x) x$pval[1]),
  ci.lb = sapply(out_int_modl, \(x) x$ci.lb[1]),
  ci.ub = sapply(out_int_modl, \(x) x$ci.ub[1]),
  row.names = NULL) |>
  mutate(across(estimate:zval, \(x) round(x, 3)),
         across(ci.lb:ci.ub, \(x) round(x, 3)),
         pval = ifelse(pval < 0.001, "<0.001",
                       round(pval, digits = 3)))


##############################################################################
# Plot prep
##############################################################################
# Add exp type to all data frames to merge together
npfert_lnRR$nut_add <- "np"
pfert_lnRR$nut_add <- "p"
nfert_lnRR$nut_add <- "n"

# Merge N, P, and NP meta results
nfert_lnRR %>%
  full_join(pfert_lnRR) %>%
  full_join(npfert_lnRR) %>%
  write.csv("../data/CNPmeta_logr_results.csv", row.names = F)

# Merge N, P, and NP meta results
out_n_modl_df %>%
  full_join(out_p_modl_df) %>%
  full_join(out_np_modl_df) %>%
  mutate(nut_add = factor(nut_add, levels = c("n", "p", "np")),
         var = factor(var, 
                      levels = c("rootshoot", "rmf", "bgb", "bnpp", "anpp_p", 
                                 "anpp_n", "agb", "anpp", "total_biomass", "tbio_gm2", 
                                 "tla", "leaf_wue", "leaf_ppue", "leaf_pnue", "jmax_vcmax",
                                 "jmax", "vcmax", "gsw", "rd", "asat", "leaf_residual_p", 
                                 "leaf_structural_p", "leaf_nucleic_p", 
                                 "leaf_metabolic_p", "leaf_pi", "leaf_pre",
                                 "leaf_nre", "leaf_np", "leaf_p_area", "leaf_p_mass", 
                                 "leaf_n_area", "leaf_n_mass", "lma")),
         estimate_se = str_c(estimate, "±", SE),
         ci_range = str_c("[", ci.lb, ", ", ci.ub, "]")) %>%
  arrange(var, nut_add) %>%
  dplyr::select(var:SE, estimate_se, zval:pval, ci.lb, ci.ub, ci_range) #%>%
  #write.csv("../data/CNPmeta_ci.csv", row.names = F)

# Factor interaction effect size variables in a certain order
out_int_modl_df %>%
  mutate(var = factor(var, 
                      levels = c("rootshoot", "rmf", "bgb", "bnpp", "anpp_p", 
                                 "anpp_n", "agb", "anpp", "total_biomass", "tbio_gm2", 
                                 "tla", "leaf_wue", "leaf_ppue", "leaf_pnue", "jmax_vcmax",
                                 "jmax", "vcmax", "gsw", "rd", "asat", "leaf_residual_p", 
                                 "leaf_structural_p", "leaf_nucleic_p", 
                                 "leaf_metabolic_p", "leaf_pi", "leaf_pre",
                                 "leaf_nre", "leaf_np", "leaf_p_area", "leaf_p_mass", 
                                 "leaf_n_area", "leaf_n_mass", "lma")),
         int_type = ifelse(var %in% c("leaf_np", "anpp", "leaf_p_mass"),
                           "synergistic", "additive"),
         estimate_se = str_c(estimate, "±", SE),
         ci_range = str_c("[", ci.lb, ", ", ci.ub, "]")) %>%
  arrange(var, nut_add) %>%
  dplyr::select(var:SE, estimate_se, zval:pval, ci.lb, ci.ub, ci_range) # %>%
  #write.csv("../data/CNPmeta_ci_int.csv", row.names = F)

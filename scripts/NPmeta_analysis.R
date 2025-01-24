# Script that explores the effect of P addition on leaf and whole-plant 
# functional traits using P fertilization and N+P fertilization 
# experiments. The meta-analysis includes data from the MESI database as of 
# January 07, 2025 and additional experiments compiled by Evan Perkowski.
# 
# Script summarizes the number of observations per trait and then conducts a 
# meta-analysis to summarize plant responses to P addition, then conducts a 
# second meta-analysis that summarizes responses when P is added in 
# concert with N.
# 
# NOTE: NEED TO CHECK THAT ADDED SAMPLE SIZE OF ADDED EXPERIMENTS IS 
# APPROPRIATE FOR MESI (SHOULD BE LISTED AS 'PLOT' REPS, NOT INDIVIDUAL REP 
# NUMBER)

# Libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readr)
library(metafor)
library(MAd)
library(patchwork)
library(ggpubr)

# MESI data
df_mesi <- read_csv("../data/mesi_main.csv") %>%
  group_by(exp) %>%
  filter(treatment == "f") %>%
  filter(any(npk == "_100")) %>%
  filter(any(npk == "_010")) %>%
  filter(any(npk == "_110"))
# Note: filter is setup to only include MESI data from experiments that have
# full-factorial combinations of N and P addition

# Manual data compilation
df_manual <- read_csv("../data/CNP_compiled_data.csv") %>%
  mutate(sampling_year = as.character(sampling_year),
         treatment = ifelse(treatment == FALSE, "f", treatment),
         sampling_date = as.character(sampling_date))

# Merge MESI database with manual data compilation
df_total <- df_mesi %>%
  full_join(df_manual)

##############################################################################
# NITROGEN FERTILIZATION RESPONSES
##############################################################################
# Subset to only include N addition treatments
explore_nfert_exps <- df_total %>%
  
  # field experiments only
  #filter(experiment_type == "field") %>%
  
  # fertilisation experiments only
  filter(treatment == "f") %>%
  
  # P-fertilisation only (without N or K addition)
  filter(npk == "_100")

head(explore_nfert_exps)

## How many experiments?
length(unique(explore_nfert_exps$exp))

## What traits are available?
unique(explore_nfert_exps$response)

# Select variables
use_response_n <- c("agb", "bgb", "total_biomass", "fine_root_biomass",
                    "rmf", "rootshoot", "lma", "leaf_n_mass", "leaf_n_area",
                    "leaf_p_mass", "leaf_p_area", "leaf_np", "sla", "vcmax",                 
                    "jmax", "tpu", "leaf_nue", "leaf_pue", "asat", "rd",                    
                    "gsw", "gs", "spad", "anet", "leaf_structure_p", 
                    "leaf_metabolic_p", "leaf_nucleic_p", "leaf_residual_p",
                    "leaf_thickness", "leaf_pi", "leaf_sugar_p", "E", 
                    "leaf_wue")              

nfert_responses <- explore_nfert_exps %>% 
  filter(response %in% use_response_n) %>% 
  mutate(myvar = response) %>%
  mutate(myvar = ifelse(myvar %in% c("anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("gs", "gsw"),
                        "gs", myvar))

use_vars_n <- unique(nfert_responses$myvar)

## Analysis

# Calculate `"ROM"` - the log transformed ratio of means (Hedges et al., 1999; 
# Lajeunesse, 2011) for each # observation pair (ambient and N added).


nfert_responses2 <- nfert_responses %>%
  
  ## keep only essential variables and drop rows containing missing values for 
  ## essential variables
  select(id, duplicate_id, exp, myvar, treatment, sampling_year, 
         x_t, x_c, sd_t, sd_c, rep_t, rep_c) %>%
  
  ## Get logarithm of response ratio and its variance
  metafor::escalc( 
    measure = "ROM", 
    m1i = x_t, sd1i = sd_t, n1i = rep_t, 
    m2i = x_c, sd2i = sd_c, n2i = rep_c, 
    data = ., 
    append = TRUE, 
    var.names = c("logr", "logr_var") 
  ) %>% 
  
  ## to keep the output readable from the console output
  as_tibble() %>% 
  
  ## get standard error
  mutate( logr_se = sqrt(logr_var) / sqrt(rep_t) )

head(nfert_responses2)

# Aggregate all measurements (multiple years, sampling dates and plots) by 
# experiment (and response variable - although here only one) for meta analysis.
nfert_responses3 <- nfert_responses2 %>% 
  
  # suggested addition by Kevin, email 02.10.2023 10:03
  dplyr::distinct(duplicate_id, x_t, x_c, .keep_all = TRUE) |> 
  
  filter(!is.na(logr_var) & !is.na(logr)) %>% 
  
  # re-create ID (common ID per experiment and response variable)
  select(-id) %>%
  mutate( id = paste(exp, myvar, sep = "_XXX_")) %>% 
  
  MAd::agg( 
    id = id, 
    es = logr, 
    var = logr_var,
    cor = 1.0, 
    method = "BHHR", 
    data = . 
  ) %>% 
  
  ## to keep the output readable from the console output
  as_tibble() %>% 
  
  # separate ID again for ease of data use
  mutate( id = str_split(id, "_XXX_") ) %>%
  mutate( exp = purrr::map_chr(id, 1),
          myvar = purrr::map_chr(id, 2) ) %>%
  
  ## rename again
  select(exp, myvar, logr = es, logr_var = var) %>%
  
  ## add number of observations (sum of plots and repeated samplings)
  left_join(
    nfert_responses2 %>%
      group_by(exp, myvar, treatment) %>%
      summarise(n_c = sum(rep_c), n_t = sum(rep_t)),
    by = c("exp", "myvar")
  ) %>% 
  
  ## get standard error. Verify if number available observations are identical
  ## for ambient and elevated. Use N from control here (n_c).
  mutate( logr_se = sqrt(logr_var) / sqrt(n_c) ,
          
          # merge SLA and LMA measurements by taking inverse of logr (keep SE)
          logr = ifelse(myvar == "sla", -logr, logr),
          myvar = ifelse(myvar == "sla", "lma", myvar))

head(nfert_responses3)


# Meta-analysis

# Aggregate log-ratios across multiple experiments, taking into account 
# their respective variance and using the experiment identity as a 
# grouping factor for random intercepts.

source("../helper_fxns/analyse_meta.R")

out_n  <- purrr::map(as.list(use_vars_n), 
                     ~analyse_meta(nfert_responses3 %>% 
                                     rename(var = myvar), nam_target = .))
names(out_n) <- use_vars_n
df_box_n <- purrr::map_dfr(out_n, "df_box") |> 
  left_join(
    nfert_responses3 |> 
      group_by(myvar) |> 
      summarise(logr_min = min(logr), logr_max = max(logr)) |> 
      rename(var = myvar),
    by = "var")

## Final data size

# Number of data points (plot-level measurements) per variable:
nfert_responses3 %>% 
  group_by(myvar) %>% 
  summarise(n_plots = sum(n_c, na.rm = TRUE), n_exp = n()) %>% 
  rename("Variable"="myvar", "N observations"="n_plots", "N experiments"="n_exp") %>% 
  knitr::kable()


# Some quick plots:
nfert_responses3 <- nfert_responses3 %>%
mutate(myvar = factor(myvar,
                      levels = c("rootshoot", "rmf", "bgb", "agb", 
                                 "total_biomass", "leaf_wue", 
                                 "leaf_pue", "leaf_nue", "rd", 
                                 "tpu", "jmax", "vcmax", 
                                 "E", "asat", "leaf_residual_p",
                                 "leaf_structure_p",
                                 "leaf_nucleic_p", "leaf_metabolic_p",
                                 "leaf_sugar_p", "leaf_pi",
                                 "leaf_np", "leaf_p_area", 
                                 "leaf_p_mass", "leaf_n_area",
                                 "leaf_n_mass", "lma")))
df_box_n$var <- factor(df_box_n$var, levels = c("rootshoot", "rmf", "bgb", "agb", 
                                                "total_biomass", "leaf_wue", 
                                                "leaf_pue", "leaf_nue", "rd", 
                                                "tpu", "jmax", "vcmax", 
                                                "E", "asat", "leaf_residual_p",
                                                "leaf_structure_p",
                                                "leaf_nucleic_p", "leaf_metabolic_p",
                                                "leaf_sugar_p", "leaf_pi",
                                                "leaf_np", "leaf_p_area", 
                                                "leaf_p_mass", "leaf_n_area",
                                                "leaf_n_mass", "lma"))

meta_plot_n <- ggplot(data = nfert_responses3, 
                      aes(x = myvar,
                          y = logr)) +
  geom_jitter(color = rgb(0,0,0,0.3), 
              aes(x = myvar,y = logr,
                  size = 1/logr_se), 
              position = position_jitter(w = 0.2, h = 0),
              show.legend = FALSE) +
  geom_crossbar(data = df_box_n %>% drop_na(var),
                 aes(x = fct_inorder(var), y = middle, ymin = ymin, ymax = ymax),
                 fill = "royalblue", 
                 color = "royalblue4", 
                 alpha = 0.6, 
                 width = 0.5) +
  geom_hline( yintercept = 0.0, linewidth = 0.5, linetype = "dotted" ) +
  scale_x_discrete(limits = c("rootshoot", "rmf", "bgb", "agb", 
                                "total_biomass", "leaf_wue", 
                                "leaf_pue", "leaf_nue", "rd", 
                                "tpu", "jmax", "vcmax", 
                                "E", "asat", "leaf_residual_p",
                                "leaf_structure_p",
                                "leaf_nucleic_p", "leaf_metabolic_p",
                                "leaf_sugar_p", "leaf_pi",
                                "leaf_np", "leaf_p_area", 
                                "leaf_p_mass", "leaf_n_area",
                                "leaf_n_mass", "lma"),
                   labels = c("Root:shoot", "Root mass fraction", 
                              "BG biomass", 
                              "AG biomass", 
                              "Total biomass", 
                              "WUEi", 
                              "PPUE", 
                              "PNUE", 
                              expression("R"["d"]), 
                              "TPU", 
                              expression("J"["max"]), 
                              expression("V"["cmax"]),
                              "E",
                              expression("A"["sat"]),  
                              "Leaf residual P",
                              "Leaf structural P",
                              "Leaf nucleic P", 
                              "Leaf metabolic P",
                              "Leaf soluble P", 
                              "Leaf inorganic P",
                              "Leaf N:P",
                              expression("P"["area"]), 
                              expression("P"["mass"]), 
                              expression("N"["area"]),
                              expression("N"["mass"]),
                              expression("M"["area"]))) +
  labs(x = "", 
       y = "Log response to N addition") +
  coord_flip() +
  theme_classic(base_size = 18)
meta_plot_n

##############################################################################
# PHOSPHORUS FERTILIZATION RESPONSES
##############################################################################
# Subset to only include N addition treatments
explore_pfert_exps <- df_total %>%
  
  # field experiments only
  #filter(experiment_type == "field") %>%
  
  # fertilisation experiments only
  filter(treatment == "f") %>%
  
  # P-fertilisation only (without N or K addition)
  filter(npk == "_010")

head(explore_pfert_exps)

## How many experiments?
length(unique(explore_pfert_exps$exp))

## What traits are available?
unique(explore_pfert_exps$response)

# Select variables
use_response_p <- c("agb", "bgb", "total_biomass", "fine_root_biomass",
                    "rmf", "rootshoot", "lma", "leaf_n_mass", "leaf_n_area",
                    "leaf_p_mass", "leaf_p_area", "leaf_np", "sla", "vcmax",
                    "jmax", "tpu", "leaf_nue", "leaf_pue", "asat", "rd",
                    "gsw", "gs", "spad", "anet", "leaf_structure_p",
                    "leaf_metabolic_p", "leaf_nucleic_p", "leaf_residual_p",
                    "leaf_thickness", "leaf_pi", "leaf_sugar_p", "E",
                    "leaf_wue")              

pfert_responses <- explore_pfert_exps %>% 
  filter(response %in% use_response_p) %>% 
  mutate(myvar = response) %>%
  mutate(myvar = ifelse(myvar %in% c("anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("gs", "gsw"),
                        "gs", myvar))

use_vars_p <- unique(pfert_responses$myvar)

## Analysis

# Calculate `"ROM"` - the log transformed ratio of means (Hedges et al., 1999; 
# Lajeunesse, 2011) for each # observation pair (ambient and N added).


pfert_responses2 <- pfert_responses %>%
  
  ## keep only essential variables and drop rows containing missing values for 
  ## essential variables
  select(id, duplicate_id, exp, myvar, treatment, sampling_year, 
         x_t, x_c, sd_t, sd_c, rep_t, rep_c) %>%
  
  ## Get logarithm of response ratio and its variance
  metafor::escalc( 
    measure = "ROM", 
    m1i = x_t, sd1i = sd_t, n1i = rep_t, 
    m2i = x_c, sd2i = sd_c, n2i = rep_c, 
    data = ., 
    append = TRUE, 
    var.names = c("logr", "logr_var") 
  ) %>% 
  
  ## to keep the output readable from the console output
  as_tibble() %>% 
  
  ## get standard error
  mutate( logr_se = sqrt(logr_var) / sqrt(rep_t) )

head(pfert_responses2)

# Aggregate all measurements (multiple years, sampling dates and plots) by 
# experiment (and response variable - although here only one) for meta analysis.
pfert_responses3 <- pfert_responses2 %>% 
  
  # suggested addition by Kevin, email 02.10.2023 10:03
  dplyr::distinct(duplicate_id, x_t, x_c, .keep_all = TRUE) |> 
  
  filter(!is.na(logr_var) & !is.na(logr)) %>% 
  
  # re-create ID (common ID per experiment and response variable)
  select(-id) %>%
  mutate( id = paste(exp, myvar, sep = "_XXX_")) %>% 
  
  MAd::agg( 
    id = id, 
    es = logr, 
    var = logr_var,
    cor = 1.0, 
    method = "BHHR", 
    data = . 
  ) %>% 
  
  ## to keep the output readable from the console output
  as_tibble() %>% 
  
  # separate ID again for ease of data use
  mutate( id = str_split(id, "_XXX_") ) %>%
  mutate( exp = purrr::map_chr(id, 1),
          myvar = purrr::map_chr(id, 2) ) %>%
  
  ## rename again
  select(exp, myvar, logr = es, logr_var = var) %>%
  
  ## add number of observations (sum of plots and repeated samplings)
  left_join(
    pfert_responses2 %>%
      group_by(exp, myvar, treatment) %>%
      summarise(n_c = sum(rep_c), n_t = sum(rep_t)),
    by = c("exp", "myvar")
  ) %>% 
  
  ## get standard error. Verify if number available observations are identical
  ## for ambient and elevated. Use N from control here (n_c).
  mutate( logr_se = sqrt(logr_var) / sqrt(n_c) ,
          
          # merge SLA and LMA measurements by taking inverse of logr (keep SE)
          logr = ifelse(myvar == "sla", -logr, logr),
          myvar = ifelse(myvar == "sla", "lma", myvar))

head(pfert_responses3)


# Meta-analysis

# Aggregate log-ratios across multiple experiments, taking into account 
# their respective variance and using the experiment identity as a 
# grouping factor for random intercepts.
out_p  <- purrr::map(as.list(use_vars_p), 
                     ~analyse_meta(pfert_responses3 %>% 
                                     rename(var = myvar), nam_target = .))
names(out_p) <- use_vars_p
df_box_p <- purrr::map_dfr(out_p, "df_box") |> 
  left_join(
    pfert_responses3 |> 
      group_by(myvar) |> 
      summarise(logr_min = min(logr), logr_max = max(logr)) |> 
      rename(var = myvar),
    by = "var")

## Final data size

# Number of data points (plot-level measurements) per variable:
pfert_responses3 %>% 
  group_by(myvar) %>% 
  summarise(n_plots = sum(n_c, na.rm = TRUE), n_exp = n()) %>% 
  rename("Variable"="myvar", "N observations"="n_plots", "N experiments"="n_exp") %>% 
  knitr::kable()


# Some quick plots:
pfert_responses3 <- pfert_responses3 %>%
  mutate(myvar = factor(myvar,
                        levels = c("rootshoot", "rmf", "bgb", "agb", 
                                   "total_biomass", "leaf_wue", 
                                   "leaf_pue", "leaf_nue", "rd", 
                                   "tpu", "jmax", "vcmax", 
                                   "E", "asat", "leaf_residual_p",
                                   "leaf_structure_p",
                                   "leaf_nucleic_p", "leaf_metabolic_p",
                                   "leaf_sugar_p", "leaf_pi",
                                   "leaf_np", "leaf_p_area", 
                                   "leaf_p_mass", "leaf_n_area",
                                   "leaf_n_mass", "lma")))
df_box_p$var <- factor(df_box_p$var, levels = c("rootshoot", "rmf", "bgb", "agb", 
                                                "total_biomass", "leaf_wue", 
                                                "leaf_pue", "leaf_nue", "rd", 
                                                "tpu", "jmax", "vcmax", 
                                                "E", "asat", "leaf_residual_p",
                                                "leaf_structure_p",
                                                "leaf_nucleic_p", "leaf_metabolic_p",
                                                "leaf_sugar_p", "leaf_pi",
                                                "leaf_np", "leaf_p_area", 
                                                "leaf_p_mass", "leaf_n_area",
                                                "leaf_n_mass", "lma"))

meta_plot_p <- ggplot(data = pfert_responses3, 
                      aes(x = myvar,
                          y = logr)) +
  geom_jitter(color = rgb(0,0,0,0.3), 
              aes(x = myvar,y = logr,
                  size = 1/logr_se), 
              position = position_jitter(w = 0.2, h = 0),
              show.legend = FALSE) +
  geom_crossbar(data = df_box_p %>% drop_na(var),
                aes(x = fct_inorder(var), y = middle, ymin = ymin, ymax = ymax),
                fill = "royalblue", 
                color = "royalblue4", 
                alpha = 0.6, 
                width = 0.5) +
  geom_hline( yintercept = 0.0, linewidth = 0.5, linetype = "dotted" ) +
  scale_x_discrete(limits = c("rootshoot", "rmf", "bgb", "agb", 
                              "total_biomass", "leaf_wue", 
                              "leaf_pue", "leaf_nue", "rd", 
                              "tpu", "jmax", "vcmax", 
                              "E", "asat", "leaf_residual_p",
                              "leaf_structure_p",
                              "leaf_nucleic_p", "leaf_metabolic_p",
                              "leaf_sugar_p", "leaf_pi",
                              "leaf_np", "leaf_p_area", 
                              "leaf_p_mass", "leaf_n_area",
                              "leaf_n_mass", "lma"),
                   labels = c("Root:shoot", "Root mass fraction", 
                              "BG biomass", 
                              "AG biomass", 
                              "Total biomass", 
                              "WUEi", 
                              "PPUE", 
                              "PNUE", 
                              expression("R"["d"]), 
                              "TPU", 
                              expression("J"["max"]), 
                              expression("V"["cmax"]),
                              "E",
                              expression("A"["sat"]),  
                              "Leaf residual P",
                              "Leaf structural P",
                              "Leaf nucleic P", 
                              "Leaf metabolic P",
                              "Leaf soluble P", 
                              "Leaf inorganic P",
                              "Leaf N:P",
                              expression("P"["area"]), 
                              expression("P"["mass"]), 
                              expression("N"["area"]),
                              expression("N"["mass"]),
                              expression("M"["area"]))) +
  labs(x = "", 
       y = "Log response to P addition") +
  coord_flip() +
  theme_classic(base_size = 18)
meta_plot_p

##############################################################################
# NITROGEN + PHOSPHORUS FERTILIZATION RESPONSES
##############################################################################
# Subset to only include N+P addition treatments
explore_npfert_exps <- df_total %>%
  
  # field experiments only
  #filter(experiment_type == "field") %>%
  
  # fertilisation experiments only
  filter(treatment == "f") %>%
  
  # P-fertilisation only (without N or K addition)
  filter(npk == "_110")

head(explore_npfert_exps)

## How many experiments?
length(unique(explore_npfert_exps$exp))

## What traits are available?
unique(explore_npfert_exps$response)

# Select variables
use_response_np <- c("agb", "bgb", "total_biomass", "fine_root_biomass",
                    "rmf", "rootshoot", "lma", "leaf_n_mass", "leaf_n_area",
                    "leaf_p_mass", "leaf_p_area", "leaf_np", "sla", "vcmax",
                    "jmax", "tpu", "leaf_nue", "leaf_pue", "asat", "rd",
                    "gsw", "gs", "spad", "anet", "leaf_structure_p",
                    "leaf_metabolic_p", "leaf_nucleic_p", "leaf_residual_p",
                    "leaf_thickness", "leaf_pi", "leaf_sugar_p", "E",
                    "leaf_wue")              

npfert_responses <- explore_npfert_exps %>% 
  filter(response %in% use_response_np) %>% 
  mutate(myvar = response) %>%
  mutate(myvar = ifelse(myvar %in% c("anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("gs", "gsw"),
                        "gs", myvar))

use_vars_np <- unique(npfert_responses$myvar)

## Analysis

# Calculate `"ROM"` - the log transformed ratio of means (Hedges et al., 1999; 
# Lajeunesse, 2011) for each # observation pair (ambient and N added).


npfert_responses2 <- npfert_responses %>%
  
  ## keep only essential variables and drop rows containing missing values for 
  ## essential variables
  select(id, duplicate_id, exp, myvar, treatment, sampling_year, 
         x_t, x_c, sd_t, sd_c, rep_t, rep_c) %>%
  
  ## Get logarithm of response ratio and its variance
  metafor::escalc( 
    measure = "ROM", 
    m1i = x_t, sd1i = sd_t, n1i = rep_t, 
    m2i = x_c, sd2i = sd_c, n2i = rep_c, 
    data = ., 
    append = TRUE, 
    var.names = c("logr", "logr_var") 
  ) %>% 
  
  ## to keep the output readable from the console output
  as_tibble() %>% 
  
  ## get standard error
  mutate( logr_se = sqrt(logr_var) / sqrt(rep_t) )

head(npfert_responses2)

# Aggregate all measurements (multiple years, sampling dates and plots) by 
# experiment (and response variable - although here only one) for meta analysis.
npfert_responses3 <- npfert_responses2 %>% 
  
  # suggested addition by Kevin, email 02.10.2023 10:03
  dplyr::distinct(duplicate_id, x_t, x_c, .keep_all = TRUE) |> 
  
  filter(!is.na(logr_var) & !is.na(logr)) %>% 
  
  # re-create ID (common ID per experiment and response variable)
  select(-id) %>%
  mutate( id = paste(exp, myvar, sep = "_XXX_")) %>% 
  
  MAd::agg( 
    id = id, 
    es = logr, 
    var = logr_var,
    cor = 1.0, 
    method = "BHHR", 
    data = . 
  ) %>% 
  
  ## to keep the output readable from the console output
  as_tibble() %>% 
  
  # separate ID again for ease of data use
  mutate( id = str_split(id, "_XXX_") ) %>%
  mutate( exp = purrr::map_chr(id, 1),
          myvar = purrr::map_chr(id, 2) ) %>%
  
  ## rename again
  select(exp, myvar, logr = es, logr_var = var) %>%
  
  ## add number of observations (sum of plots and repeated samplings)
  left_join(
    pfert_responses2 %>%
      group_by(exp, myvar, treatment) %>%
      summarise(n_c = sum(rep_c), n_t = sum(rep_t)),
    by = c("exp", "myvar")
  ) %>% 
  
  ## get standard error. Verify if number available observations are identical
  ## for ambient and elevated. Use N from control here (n_c).
  mutate( logr_se = sqrt(logr_var) / sqrt(n_c) ,
          
          # merge SLA and LMA measurements by taking inverse of logr (keep SE)
          logr = ifelse(myvar == "sla", -logr, logr),
          myvar = ifelse(myvar == "sla", "lma", myvar))

head(npfert_responses3)


# Meta-analysis

# Aggregate log-ratios across multiple experiments, taking into account 
# their respective variance and using the experiment identity as a 
# grouping factor for random intercepts.
out_np  <- purrr::map(as.list(use_vars_np), 
                     ~analyse_meta(npfert_responses3 %>% 
                                     rename(var = myvar), nam_target = .))
names(out_np) <- use_vars_np
df_box_np <- purrr::map_dfr(out_np, "df_box") |> 
  left_join(
    npfert_responses3 |> 
      group_by(myvar) |> 
      summarise(logr_min = min(logr), logr_max = max(logr)) |> 
      rename(var = myvar),
    by = "var")

## Final data size

# Number of data points (plot-level measurements) per variable:
npfert_responses3 %>% 
  group_by(myvar) %>% 
  summarise(n_plots = sum(n_c, na.rm = TRUE), n_exp = n()) %>% 
  rename("Variable"="myvar", "N observations"="n_plots", "N experiments"="n_exp") %>% 
  knitr::kable()


# Some quick plots:
npfert_responses3 <- npfert_responses3 %>%
  mutate(myvar = factor(myvar,
                        levels = c("rootshoot", "rmf", "bgb", "agb", 
                                   "total_biomass", "leaf_wue", 
                                   "leaf_pue", "leaf_nue", "rd", 
                                   "tpu", "jmax", "vcmax", 
                                   "E", "asat", "leaf_residual_p",
                                   "leaf_structure_p",
                                   "leaf_nucleic_p", "leaf_metabolic_p",
                                   "leaf_sugar_p", "leaf_pi",
                                   "leaf_np", "leaf_p_area", 
                                   "leaf_p_mass", "leaf_n_area",
                                   "leaf_n_mass", "lma")))
df_box_np$var <- factor(df_box_np$var, levels = c("rootshoot", "rmf", "bgb", "agb", 
                                                "total_biomass", "leaf_wue", 
                                                "leaf_pue", "leaf_nue", "rd", 
                                                "tpu", "jmax", "vcmax", 
                                                "E", "asat", "leaf_residual_p",
                                                "leaf_structure_p",
                                                "leaf_nucleic_p", "leaf_metabolic_p",
                                                "leaf_sugar_p", "leaf_pi",
                                                "leaf_np", "leaf_p_area", 
                                                "leaf_p_mass", "leaf_n_area",
                                                "leaf_n_mass", "lma"))

meta_plot_np <- ggplot(data = npfert_responses3, 
                      aes(x = myvar,
                          y = logr)) +
  geom_jitter(color = rgb(0,0,0,0.3), 
              aes(x = myvar,y = logr,
                  size = 1/logr_se), 
              position = position_jitter(w = 0.2, h = 0),
              show.legend = FALSE) +
  geom_crossbar(data = df_box_np %>% drop_na(var),
                aes(x = fct_inorder(var), y = middle, ymin = ymin, ymax = ymax),
                fill = "royalblue", 
                color = "royalblue4", 
                alpha = 0.6, 
                width = 0.5) +
  geom_hline( yintercept = 0.0, linewidth = 0.5, linetype = "dotted" ) +
  scale_x_discrete(limits = c("rootshoot", "rmf", "bgb", "agb", 
                              "total_biomass", "leaf_wue", 
                              "leaf_pue", "leaf_nue", "rd", 
                              "tpu", "jmax", "vcmax", 
                              "E", "asat", "leaf_residual_p",
                              "leaf_structure_p",
                              "leaf_nucleic_p", "leaf_metabolic_p",
                              "leaf_sugar_p", "leaf_pi",
                              "leaf_np", "leaf_p_area", 
                              "leaf_p_mass", "leaf_n_area",
                              "leaf_n_mass", "lma"),
                   labels = c("Root:shoot", "Root mass fraction", 
                              "BG biomass", 
                              "AG biomass", 
                              "Total biomass", 
                              "WUEi", 
                              "PPUE", 
                              "PNUE", 
                              expression("R"["d"]), 
                              "TPU", 
                              expression("J"["max"]), 
                              expression("V"["cmax"]),
                              "E",
                              expression("A"["sat"]),  
                              "Leaf residual P",
                              "Leaf structural P",
                              "Leaf nucleic P", 
                              "Leaf metabolic P",
                              "Leaf soluble P", 
                              "Leaf inorganic P",
                              "Leaf N:P",
                              expression("P"["area"]), 
                              expression("P"["mass"]), 
                              expression("N"["area"]),
                              expression("N"["mass"]),
                              expression("M"["area"]))) +
  labs(x = "", 
       y = "Log response to N+P addition") +
  coord_flip() +
  theme_classic(base_size = 18)
meta_plot_np

##############################################################################
# Let's put together some plots
##############################################################################

# Add exp type to all data frames to merge together
df_box_np$manip_type <- "np"
df_box_p$manip_type <- "p"
df_box_n$manip_type <- "n"

npfert_responses3$manip_type <- "np"
pfert_responses3$manip_type <- "p"
nfert_responses3$manip_type <- "n"

# Merge P and NP meta results
df_box_all <- df_box_n %>%
  full_join(df_box_p) %>%
  full_join(df_box_np) %>%
  mutate(manip_type = factor(manip_type, levels = c("np", "p", "n")))

fert_exp_responses_all <- nfert_responses3 %>%
  full_join(pfert_responses3) %>%
  full_join(npfert_responses3) %>%
  mutate(manip_type = factor(manip_type, levels = c("np", "p", "n")))
head(fert_exp_responses_all)

# Plot nutrients. Separating by trait type to avoid plot overwhelm
meta_plot_all_leaf_nutrients <- ggplot(
  data = subset(fert_exp_responses_all, 
                myvar %in% c("lma", "leaf_n_mass", "leaf_n_area", "leaf_p_mass",
                             "leaf_p_area", "leaf_np", "spad")),
  aes(x = myvar, y = logr, fill = manip_type)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/logr_se)) +
  geom_crossbar(data = df_box_all %>% drop_na(var) %>%
                  filter(var %in% c("lma", "leaf_n_mass", "leaf_n_area", 
                                    "leaf_p_mass", "leaf_p_area", "leaf_np", 
                                    "spad")),
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.5,
                position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("SPAD", 
                              "Leaf N:P",
                              expression("P"["area"]), 
                              expression("P"["mass"]), 
                              expression("N"["area"]),
                              expression("N"["mass"]), 
                              expression("M"["area"]))) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_fill_manual(limits = c("n", "p", "np"),
                    values = c("red", "blue", "magenta")) +
  scale_size(range = c(0.25, 4)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition",
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_plot_all_leaf_nutrients

# Plot photosynthetic traits. Separating by trait type to avoid plot overwhelm
meta_plot_all_photo <- ggplot(
  data = subset(fert_exp_responses_all, 
                myvar %in% c("asat", "vcmax", "jmax",
                             "rd", "leaf_nue", "leaf_pue")),
  aes(x = myvar, y = logr, fill = manip_type)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/logr_se)) +
  geom_crossbar(data = df_box_all %>% drop_na(var) %>%
                  filter(var %in% c("asat", "vcmax", "jmax",
                                    "rd", "leaf_nue", "leaf_pue")),
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.5,
                position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("PPUE",
                              "PNUE",
                              expression("Rd"), 
                              expression("J"["max"]),
                              expression("V"["cmax"]),
                              expression("g"["s"]),
                              expression("A"["net"]))) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_fill_manual(limits = c("n", "p", "np"),
                    values = c("red", "blue", "magenta")) +
  scale_size(range = c(0.25, 4)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition",
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_plot_all_photo

# Plot biomass traits. Separating by trait type to avoid plot overwhelm
meta_plot_all_biomass <- ggplot(
  data = subset(fert_exp_responses_all, 
                myvar %in% c("rootshoot", "rmf", "bgb", "agb", "total_biomass")),
  aes(x = myvar, y = logr, fill = manip_type)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/logr_se)) +
  geom_crossbar(data = df_box_all %>% drop_na(var) %>%
                  filter(var %in% c("rootshoot", "rmf", "bgb", "agb", "total_biomass")),
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.5,
                position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Root:shoot",
                              "RMF",
                              "Belowground biomass",
                              "Aboveground biomass",
                              "Total biomass")) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_fill_manual(limits = c("n", "p", "np"),
                    values = c("red", "blue", "magenta")) +
  scale_size(range = c(0.25, 4)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition",
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
meta_plot_all_biomass


meta_plot_all_leaf_nutrients / meta_plot_all_photo / meta_plot_all_biomass +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 12, face = "bold"))






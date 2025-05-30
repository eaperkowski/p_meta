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
  mutate(pft = str_c(photo_path, n_fixer, myc_assoc, sep = "_"))

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
unique(experiment_summary$exp)

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
# Investigate N fertilization responses
#####################################################################

# Subset `full_df` to include only N addition treatments
nfert_only <- full_df %>% filter(npk == "_100")
head(nfert_only)

# What traits are included?
unique(nfert_only$response)

# Which experiments are included?
unique(nfert_only$citation)

# Select variables
use_response_n <- c("agb_n", "agb_p", "agb", "anpp", "agb_n_mass",
                    "agb_p_mass", "leaf_n_mass", "leaf_p_mass", "anet",
                    "total_biomass", "bgb", "rootshoot", "bnpp", "lma", "sla", "rmf",
                    "gsw", "amax", "leaf_wue", "fine_root_biomass", "leaf_n_area",
                    "leaf_p_area", "asat", "vcmax", "jmax", "leaf_pnue", "leaf_np",
                    "leaf_ppue", "leaf_structural_p", "leaf_metabolic_p",
                    "leaf_nucleic_p", "leaf_residual_p", "leaf_pi")

nfert_responses <- nfert_only %>%
  filter(response %in% use_response_n) %>%
  mutate(myvar = response) %>%
  mutate(myvar = ifelse(myvar %in% c("amax", "anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("agb", "anpp"),
                        "agb", myvar),
         myvar = ifelse(myvar %in% c("bgb", "fine_root_biomass", "bnpp"),
                        "bgb", myvar),
         myvar = ifelse(myvar %in% c("total_leaf_area", "tla"),
                        "tla", myvar))

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
                                    rename(var = myvar), nam_target = .))
names(out_n) <- use_vars_n

df_box_n_k <- data.frame(var = use_vars_n,
                          k = c("(18)", "(12)", "(125)", "(30)", "(23)", "(139)",
                                "(133)", "(85)", "(42)", "(63)", "(40)",
                                NA, "(37)", "(43)", "(115)", "(113)", "(37)", 
                                "(84)", "(79)", "(42)", "(40)", "(58)", "(59)", 
                                "(20)", "(21)", "(21)", "(20)", "(5)"),
                         sig.level = c("***", "***", "***", "***", "***", "***",
                                       "*. ", "ns ", "ns ", "ns ", "***", 
                                       NA, "***", "ns ", "** ", "*. ", "ns ", 
                                       "** ", "ns ", "ns ", "ns ", "ns ", "ns ", 
                                       "ns ", "ns ", "ns ", "ns ", "ns "))

df_box_n <- purrr::map_dfr(out_n, "df_box") |> 
  full_join(df_box_n_k) |>
  left_join(
    nfert_lnRR |> 
      group_by(myvar) |> 
      summarise(logr_min = min(logr), logr_max = max(logr)) |> 
      rename(var = myvar),
    by = "var")

#####################################################################
# Investigate P fertilization responses
#####################################################################

# Subset `full_df` to include only P addition treatments
pfert_only <- full_df %>% filter(npk == "_010")
head(pfert_only)

# What traits are included?
unique(pfert_only$response)

# Select variables
use_response_p <- c("agb_n", "agb_p", "agb", "anpp", "agb_n_mass",
                    "agb_p_mass", "leaf_n_mass", "leaf_p_mass", "anet",
                    "total_biomass", "bgb", "rootshoot", "bnpp", "lma", "sla", "rmf",
                    "gsw", "amax", "leaf_wue", "fine_root_biomass", "leaf_n_area",
                    "leaf_p_area", "asat", "vcmax", "jmax", "leaf_pnue", "leaf_np",
                    "leaf_ppue", "leaf_structural_p", "leaf_metabolic_p",
                    "leaf_nucleic_p", "leaf_residual_p", "leaf_pi")

pfert_responses <- pfert_only %>%
  filter(response %in% use_response_p) %>%
  mutate(myvar = response) %>%
  mutate(myvar = ifelse(myvar %in% c("amax", "anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("agb", "anpp"),
                        "agb", myvar),
         myvar = ifelse(myvar %in% c("bgb", "fine_root_biomass", "bnpp"),
                        "bgb", myvar),
         myvar = ifelse(myvar %in% c("total_leaf_area", "tla"),
                        "tla", myvar))

use_vars_p <- unique(pfert_responses$myvar)

# Calculate log-response ratios using `escalc` in `metafor`
pfert_lnRR <- pfert_responses %>% 
  escalc(measure = "ROM",
         m1i = x_t, sd1i = sd_t, n1i = rep_t, 
         m2i = x_c, sd2i = sd_c, n2i = rep_c, 
         data = ., 
         append = TRUE, 
         var.names = c("logr", "logr_var")) %>%
  mutate(logr_se = sqrt(logr_var) / sqrt(rep_t)) %>%
  mutate(logr = ifelse(myvar == "sla", logr * -1, logr), 
         myvar = ifelse(myvar == "sla", "lma", myvar))

# Determine lnRR across multiple experiments. Takes into account
# experiment variance and uses experiment identity as grouping
# factor for random intercepts
out_p <- purrr::map(as.list(use_vars_p),
                    ~analyse_meta(pfert_lnRR %>%
                                    rename(var = myvar), nam_target = .))
names(out_p) <- use_vars_p

df_box_p_k <- data.frame(var = use_vars_p,
                         k = c("(18)", "(12)", "(125)", "(30)", "(23)", "(139)",
                               "(133)", "(85)", "(42)", "(63)", "(40)", NA,
                               "(37)", "(43)", "(115)", "(113)", "(37)", "(84)", 
                               "(79)", "(42)", "(40)", "(58)", "(59)", "(20)",
                               "(21)", "(21)", "(20)", "(5)"),
                         sig.level = c(".  ", "***", "***", "***", "***", "ns ",
                                       "***", "ns ", "* ", "ns ", "*. ", NA,
                                       ".  ", ".  ", "***", "ns ", "ns ", "ns ", 
                                       "***", "ns ", "*. ", "ns ", "ns ", "ns ",
                                       "*. ", "ns ", "ns ", "**"))


df_box_p <- purrr::map_dfr(out_p, "df_box") |> 
  full_join(df_box_p_k) |>
  left_join(
    pfert_lnRR |> 
      group_by(myvar) |> 
      summarise(logr_min = min(logr), logr_max = max(logr)) |> 
      rename(var = myvar),
    by = "var")

#####################################################################
# Investigate N+P fertilization responses
#####################################################################

# Subset `full_df` to include only N addition treatments
npfert_only <- full_df %>% filter(npk == "_110")
head(npfert_only)

# What traits are included?
unique(npfert_only$response)

# Select variables
use_response_np <- c("agb_n", "agb_p", "agb", "anpp", "agb_n_mass",
                    "agb_p_mass", "leaf_n_mass", "leaf_p_mass", "anet",
                    "total_biomass", "bgb", "rootshoot", "bnpp", "lma", "sla", "rmf",
                    "gsw", "amax", "leaf_wue", "fine_root_biomass", "leaf_n_area",
                    "leaf_p_area", "asat", "vcmax", "jmax", "leaf_pnue", "leaf_np",
                    "leaf_ppue", "leaf_structural_p", "leaf_metabolic_p",
                    "leaf_nucleic_p", "leaf_residual_p", "leaf_pi")

npfert_responses <- npfert_only %>%
  filter(response %in% use_response_np) %>%
  mutate(myvar = response) %>%
  mutate(myvar = ifelse(myvar %in% c("amax", "anet", "asat"),
                        "asat", myvar),
         myvar = ifelse(myvar %in% c("agb", "anpp"),
                        "agb", myvar),
         myvar = ifelse(myvar %in% c("bgb", "fine_root_biomass", "bnpp"),
                        "bgb", myvar),
         myvar = ifelse(myvar %in% c("total_leaf_area", "tla"),
                        "tla", myvar)) %>%
  filter(myvar != "r_eco" & myvar != "nee" & 
           myvar != "spad" & myvar != "tpu" &
           myvar != "rd" & myvar != "stom_lim" &
           myvar != "anet_mass")

use_vars_np <- unique(npfert_responses$myvar)

# Calculate log-response ratios using `escalc` in `metafor`
npfert_lnRR <- npfert_responses %>% 
  escalc(measure = "ROM",
         m1i = x_t, sd1i = sd_t, n1i = rep_t, 
         m2i = x_c, sd2i = sd_c, n2i = rep_c, 
         data = ., 
         append = TRUE, 
         var.names = c("logr", "logr_var")) %>%
  mutate(logr_se = sqrt(logr_var) / sqrt(rep_t)) %>%
  mutate(logr = ifelse(myvar == "sla", logr * -1, logr), 
         myvar = ifelse(myvar == "sla", "lma", myvar))

# Determine lnRR across multiple experiments. Takes into account
# experiment variance and uses experiment identity as grouping
# factor for random intercepts
out_np <- purrr::map(as.list(use_vars_np),
                    ~analyse_meta(npfert_lnRR %>%
                                    rename(var = myvar), nam_target = .))
names(out_np) <- use_vars_np

df_box_np_k <- data.frame(var = use_vars_np,
                          k = c("(18)", "(12)", "(125)", "(30)", "(23)", "(139)",
                                "(133)", "(85)", "(42)", "(63)", "(40)", NA,
                                "(37)", "(43)", "(115)", "(113)", "(37)", "(84)", 
                                "(79)", "(42)", "(40)", "(58)", "(59)", "(20)",
                                "(21)", "(21)", "(20)", "(5)"),
                          sig.level = c("***", "***", "***", "** ", "*  ", "***",
                                        "***", "*  ", "***", "ns ", "***", NA,
                                        "** ", "*  ", "***", "** ", ".  ", "***", 
                                        "***", ".  ", "***", "ns ", "ns ", "ns ",
                                        "*  ", "ns ", "ns ", "** "))

df_box_np <- purrr::map_dfr(out_np, "df_box") |> 
  full_join(df_box_np_k) %>%
  left_join(
    npfert_lnRR |> 
      group_by(myvar) |> 
      summarise(logr_min = min(logr), logr_max = max(logr)) |> 
      rename(var = myvar),
    by = "var")

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
  pivot_wider(names_from = fert, values_from = c(x_t, sd_t, rep_t)) %>%
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
  mutate(myvar = response,
         myvar = ifelse(myvar == "amax", "asat", myvar),
         myvar = ifelse(myvar == "anet", "asat", myvar),
         myvar = ifelse(myvar == "fine_root_biomass" | myvar == "bnpp", 
                        "bgb", myvar),
         myvar = ifelse(myvar == "anpp", "agb", myvar),
         myvar = ifelse(myvar == "total_leaf_area", "tla", myvar)) %>%
  filter(myvar != "r_eco" & myvar != "nee" & 
           myvar != "spad" & myvar != "tpu" &
           myvar != "rd" & myvar != "stom_lim" &
           myvar != "anet_mass" & myvar!= "lai" & myvar != "agb_np" &
           myvar != "leaf_sugar_p" & myvar != "gpp")

###############################################################################
# Use helper fxn from `analyse_meta.R` to iterate through traits and determine
# meta-analytic mean of the interaction effect size and associated 95% CI
###############################################################################

# Determine interaction effect size across multiple experiments. 
# Takes into account experiment variance and uses experiment 
# identity as grouping factor for random intercepts
use_vars_int <- unique(CNP_effect_sizes_reduced$myvar)

out_int <- purrr::map(as.list(use_vars_int),
                      ~analyse_meta_int(CNP_effect_sizes_reduced %>%
                                          rename(var = myvar), nam_target = .))
names(out_int) <- use_vars_int

df_box_int <- purrr::map_dfr(out_int, "df_box") |> 
  left_join(
    CNP_effect_sizes_reduced |> 
      group_by(myvar) |> 
      summarise(intES_min = min(dNPi), intES_max = max(dNPi)) |> 
      rename(var = myvar),
    by = "var")


##############################################################################
# Plot prep
##############################################################################
# Add exp type to all data frames to merge together
df_box_np$manip_type <- "np"
df_box_p$manip_type <- "p"
df_box_n$manip_type <- "n"

npfert_lnRR$manip_type <- "np"
pfert_lnRR$manip_type <- "p"
nfert_lnRR$manip_type <- "n"

# Merge N, P, and NP meta results
df_box_all <- df_box_n %>%
  full_join(df_box_p) %>%
  full_join(df_box_np) %>%
  filter(!is.na(middle)) %>%
  mutate(manip_type = factor(manip_type, levels = c("np", "p", "n")),
         var = factor(var, 
                        levels = c("rootshoot", "rmf", "bgb", "agb_p", 
                                   "agb_p_mass", "agb_n", "agb_n_mass",
                                   "agb", "total_biomass", 
                                   "leaf_wue", "leaf_ppue", "leaf_pnue", 
                                   "jmax", "vcmax", "gsw", "asat", "leaf_residual_p", 
                                   "leaf_structural_p", "leaf_nucleic_p", 
                                   "leaf_metabolic_p", "leaf_pi", 
                                   "leaf_np", "leaf_p_area", "leaf_p_mass", 
                                   "leaf_n_area", "leaf_n_mass", "lma"))) %>%
  mutate(middle_fixed = ifelse(middle >= 0,
                               sprintf(" %.3f", middle),
                               sprintf("%.3f", middle)),
         plot_label = paste0("paste('", middle_fixed, "'^'", sig.level, "')"),
         bold_yn = ifelse(sig.level %in% c("*", "**", "***"),
                          "bold", "plain"))

# Merge N, P, and NP lnRR values
fert_exp_responses_all <- nfert_lnRR %>%
  full_join(pfert_lnRR) %>%
  full_join(npfert_lnRR) %>%
  mutate(manip_type = factor(manip_type, levels = c("np", "p", "n")),
         myvar = factor(myvar, 
                        levels = c("rootshoot", "rmf", "bgb", "agb_p", 
                                   "agb_p_mass", "agb_n", "agb_n_mass",
                                   "agb", "total_biomass", 
                                   "leaf_wue", "leaf_ppue", "leaf_pnue", 
                                   "jmax", "vcmax", "gsw", "asat", "leaf_residual_p", 
                                   "leaf_structural_p", "leaf_nucleic_p", 
                                   "leaf_metabolic_p", "leaf_pi", 
                                   "leaf_np", "leaf_p_area", "leaf_p_mass", 
                                   "leaf_n_area", "leaf_n_mass", "lma")))

# Factor interaction effect size variables to appear in a certain order
CNP_effect_sizes_reduced <- CNP_effect_sizes_reduced %>%
  mutate(myvar = factor(myvar, 
                        levels = c("rootshoot", "rmf", "bgb", "agb_p", 
                                   "agb_n", "agb", "total_biomass", 
                                   "leaf_wue", "leaf_ppue", "leaf_pnue", 
                                   "jmax", "vcmax", "asat", "leaf_residual_p", 
                                   "leaf_structural_p", "leaf_nucleic_p", 
                                   "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi", 
                                   "leaf_np", "leaf_p_area", "leaf_p_mass", 
                                   "leaf_n_area", "leaf_n_mass", "lma", "sla")))

# Factor interaction effect size variables in a certain order
df_box_int <- df_box_int %>%
  mutate(var = factor(var, 
                      levels = c("rootshoot", "rmf", "bgb", "agb_p", 
                                 "agb_n", "agb", "total_biomass", 
                                 "leaf_wue", "leaf_ppue", "leaf_pnue", 
                                 "jmax", "vcmax", "asat", "leaf_residual_p", 
                                 "leaf_structural_p", "leaf_nucleic_p", 
                                 "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi", 
                                 "leaf_np", "leaf_p_area", "leaf_p_mass", 
                                 "leaf_n_area", "leaf_n_mass", "lma", "sla")),
         int_type = ifelse(var == "leaf_np", 
                           "synergistic",
                           ifelse(var == "agb", 
                                  "synergistic",
                                  "additive")))

##############################################################################
# Leaf nutrient plots
##############################################################################
# Individual leaf nutrient responses
meta_plot_all_leaf_nutrients <- ggplot(
  data = subset(fert_exp_responses_all, 
                myvar %in% c("lma", "leaf_n_mass", "leaf_n_area", "leaf_p_mass",
                             "leaf_p_area", "leaf_np")),
  aes(x = myvar, y = logr, fill = manip_type)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/logr_se), alpha = 0.3) +
  geom_crossbar(data = df_box_all %>% drop_na(var) %>%
                  filter(var %in% c("lma", "leaf_n_mass", "leaf_n_area", 
                                    "leaf_p_mass", "leaf_p_area", "leaf_np")),
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.6,
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf N:P",
                              expression("P"["area"]), 
                              expression("P"["mass"]), 
                              expression("N"["area"]),
                              expression("N"["mass"]), 
                              expression("M"["area"]),
                              "SLA")) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_fill_manual(limits = c("n", "p", "np"),
                    values = c("red", "blue", "magenta")) +
  scale_size_continuous(limits = c(0, 224), range = c(0.2, 4)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition",
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))+
  guides(fill = guide_legend(order = 1),
         size = guide_legend(override.aes = list(alpha = 1),
                             order = 2))
meta_plot_all_leaf_nutrients

# Interaction effect sizes for leaf nutrients
meta_intPlot_leaf_nutrients <- ggplot(
  data = CNP_effect_sizes_reduced %>%
    filter(myvar %in% c("lma", "leaf_n_mass", "leaf_n_area", 
                        "leaf_p_mass", "leaf_p_area", "leaf_np")),
  aes(x = myvar, y = dNPi)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/dNPi_se), alpha = 0.3) +
  geom_crossbar(data = df_box_int %>% drop_na(var) %>%
                  filter(var %in% c("lma", "leaf_n_mass", "leaf_n_area", 
                                    "leaf_p_mass", "leaf_p_area", "leaf_np")),
                aes(x = var, y = middle, ymin = ymin, 
                    ymax = ymax, fill = int_type),
                alpha = 0.6, width = 0.4,
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf N:P",
                              expression("P"["area"]), 
                              expression("P"["mass"]), 
                              expression("N"["area"]),
                              expression("N"["mass"]), 
                              expression("M"["area"]),
                              "SLA")) +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1)) +
  scale_size_continuous(limits = c(0, 15), range = c(0.2, 4)) +
  labs(x = "", 
       y = "Interaction effect size (Hedge's d)",
       size = expression(bold("Error"^"-1")),
       fill = "Interaction type") +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold")) +
  guides(fill = guide_legend(order = 1),
         size = guide_legend(override.aes = list(alpha = 1),
                             order = 2))
meta_intPlot_leaf_nutrients

png("../plots/CNPmeta_figXX_leaf_nutrients.png", width = 10, height = 12,
    units = "in", res = 600)
ggarrange(meta_plot_all_leaf_nutrients, meta_intPlot_leaf_nutrients,
          nrow = 2, labels = c("(a)", "(b)"), font.label = list(size = 18))
dev.off()

##############################################################################
# Leaf phosphorus fractionation plots (removing residual P pool)
##############################################################################
# Individual leaf nutrient responses
meta_plot_all_phosFract <- ggplot(
  data = subset(fert_exp_responses_all, 
                myvar %in% c("leaf_structural_p", "leaf_nucleic_p", 
                             "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi")),
  aes(x = myvar, y = logr, fill = manip_type)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/logr_se), alpha = 0.3) +
  geom_crossbar(data = df_box_all %>% drop_na(var) %>%
                  filter(var %in% c("leaf_structural_p", "leaf_nucleic_p", 
                                    "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi")),
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.6,
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf structural P",
                              "Leaf nucleic acid P",
                              "Leaf metabolic P",
                              "Leaf Pi")) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_fill_manual(limits = c("n", "p", "np"),
                    values = c("red", "blue", "magenta")) +
  scale_size_continuous(limits = c(0, 224), range = c(0.2, 4)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition",
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold")) +
  guides(fill = guide_legend(order = 1),
         size = guide_legend(override.aes = list(alpha = 1),
                             order = 2))
meta_plot_all_phosFract

# Interaction effect sizes for leaf phosphorus fractionation
meta_intPlot_phosFract <- ggplot(
  data = CNP_effect_sizes_reduced %>%
    filter(myvar %in%  c("leaf_structural_p", "leaf_nucleic_p", 
                         "leaf_metabolic_p", "leaf_pi")),
  aes(x = myvar, y = dNPi)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1 / dNPi_se), alpha = 0.3) +
  geom_crossbar(data = df_box_int %>% drop_na(var) %>%
                  filter(var %in% c("leaf_structural_p", "leaf_nucleic_p", 
                                    "leaf_metabolic_p", "leaf_sugar_p", "leaf_pi")),
                aes(x = var, y = middle, ymin = ymin, 
                    ymax = ymax, fill = int_type),
                alpha = 0.6, width = 0.4,
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf structural P",
                              "Leaf nucleic acid P",
                              "Leaf metabolic P",
                              "Leaf Pi")) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_size_continuous(limits = c(0, 8), breaks = seq(0, 9, 3), range = c(0.2, 4)) +
  labs(x = "", 
       y = "Interaction effect size (Hedge's d)",
       size = expression(bold("Error"^"-1")),
       fill = "Interaction type") +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold")) +
  guides(fill = guide_legend(order = 1),
         size = guide_legend(override.aes = list(alpha = 1),
                             order = 2))
meta_intPlot_phosFract

png("../plots/CNPmeta_figXX_phos_fract.png", width = 10, height = 12,
    units = "in", res = 600)
ggarrange(meta_plot_all_phosFract, meta_intPlot_phosFract,
          nrow = 2, labels = c("(a)", "(b)"), font.label = list(size = 18))
dev.off()

##############################################################################
# Leaf photosynthesis plots
##############################################################################
# Individual leaf nutrient responses
meta_plot_all_photo <- ggplot(
  data = subset(fert_exp_responses_all, 
                myvar %in% c("asat", "vcmax", "jmax",
                             "rd", "leaf_pnue", "leaf_ppue")),
  aes(x = myvar, y = logr, fill = manip_type)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/logr_se), alpha = 0.3) +
  geom_crossbar(data = df_box_all %>% drop_na(var) %>%
                  filter(var %in% c("asat", "vcmax", "jmax",
                                    "rd", "leaf_pnue", "leaf_ppue")),
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.6,
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("PPUE",
                              "PNUE",
                              expression("J"["max"]),
                              expression("V"["cmax"]),
                              expression("A"["sat"]))) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_fill_manual(limits = c("n", "p", "np"),
                    values = c("red", "blue", "magenta")) +
  scale_size_continuous(limits = c(0, 224), range = c(0.2, 4)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition",
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold")) +
  guides(fill = guide_legend(order = 1),
         size = guide_legend(override.aes = list(alpha = 1),
                             order = 2))
meta_plot_all_photo

# Interaction effect sizes for photosynthetic traits
meta_intPlot_photo <- ggplot(
  data = CNP_effect_sizes_reduced %>%
    filter(myvar %in% c("asat", "vcmax", "jmax",
                        "rd", "leaf_pnue", "leaf_ppue")),
  aes(x = myvar, y = dNPi)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/dNPi_se), alpha = 0.3) +
  geom_crossbar(data = df_box_int %>% drop_na(var) %>%
                  filter(var %in% c("asat", "vcmax", "jmax",
                                    "rd", "leaf_pnue", "leaf_ppue")),
                aes(x = var, y = middle, ymin = ymin, 
                    ymax = ymax, fill = int_type),
                alpha = 0.6, width = 0.4, 
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("PPUE",
                              "PNUE",
                              expression("J"["max"]),
                              expression("V"["cmax"]),
                              expression("A"["sat"]))) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_size_continuous(limits = c(0, 12), breaks = seq(0, 12, 4), range = c(0.2, 4)) +
  labs(x = NULL, 
       y = "Interaction effect size (Hedge's d)",
       fill = "Interaction type",
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold")) +
    guides(fill = guide_legend(order = 1),
           size = guide_legend(override.aes = list(alpha = 1),
                               order = 2))
meta_intPlot_photo

# Write plot
png("../plots/CNPmeta_figXX_photo.png", width = 10, height = 12,
    units = "in", res = 600)
ggarrange(meta_plot_all_photo, meta_intPlot_photo,
          nrow = 2, labels = c("(a)", "(b)"), font.label = list(size = 18))
dev.off()


##############################################################################
# Whole-plant trait plots
##############################################################################
# Individual leaf nutrient responses
meta_plot_all_biomass <- ggplot(
  data = subset(fert_exp_responses_all, 
                myvar %in% c("rootshoot", "rmf", "bgb", "agb", "agb_n", 
                             "agb_p", "total_biomass")),
  aes(x = myvar, y = logr, fill = manip_type)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1/logr_var), alpha = 0.3) +
  geom_crossbar(data = df_box_all %>% drop_na(var) %>%
                  filter(var %in% c("rootshoot", "rmf", "bgb", "agb", "agb_n", 
                                    "agb_p", "total_biomass")),
                aes(x = var, y = middle, ymin = ymin, ymax = ymax),
                alpha = 0.6, width = 0.6,
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Root:shoot",
                              "RMF",
                              "Belowground biomass",
                              "Aboveground P standing stock",
                              "Aboveground N standing stock",
                              "Aboveground biomass",
                              "Total biomass")) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_fill_manual(limits = c("n", "p", "np"),
                    values = c("red", "blue", "magenta")) +
  scale_size_continuous(limits = c(0, 224), range = c(0.2, 4)) +
  labs(x = "", 
       y = "Log response to N, P, or N+P addition",
       fill = "Nutrient addition",
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold")) +
  guides(fill = guide_legend(order = 1),
         size = guide_legend(override.aes = list(alpha = 1),
                             order = 2))
meta_plot_all_biomass

# Interaction effect sizes for whole-plant traits
meta_intPlot_biomass <- ggplot(
  data = subset(CNP_effect_sizes_reduced, 
                myvar %in% c("rootshoot", "rmf", "bgb", "agb", "agb_n", 
                             "agb_p", "total_biomass")),
  aes(x = myvar, y = dNPi)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75),
              shape = 21, aes(size = 1 / dNPi_se), alpha = 0.3) +
  geom_crossbar(data = df_box_int %>% drop_na(var) %>%
                  filter(var %in% c("rootshoot", "rmf", "bgb", "agb", "agb_n", 
                                    "agb_p", "total_biomass")),
                aes(x = var, y = middle, ymin = ymin, 
                    ymax = ymax, fill = int_type),
                alpha = 0.6, width = 0.4,
                position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Root:shoot",
                              "RMF",
                              "Belowground biomass",
                              "Aboveground P standing stock",
                              "Aboveground N standing stock",
                              "Aboveground biomass",
                              "Total biomass")) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_size_continuous(limits = c(0, 60), breaks = seq(0, 60, 15), range = c(0.2, 4)) +
  labs(x = "", 
       y = "Interaction effect size (Hedge's d)",
       fill = "Interaction type",
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold")) +
  guides(fill = guide_legend(order = 1),
         size = guide_legend(override.aes = list(alpha = 1),
                             order = 2))
meta_intPlot_biomass

# Write plot
png("../plots/CNPmeta_figXX_biomass.png", width = 12, height = 12,
    units = "in", res = 600)
ggarrange(meta_plot_all_biomass, meta_intPlot_biomass,
          nrow = 2, labels = c("(a)", "(b)"), font.label = list(size = 18))
dev.off()

###############################################################################
# lnRR trait summaries
###############################################################################
tableXX_lnRR_summary <- data.frame(
  trait = df_box_all$var,
  response = df_box_all$manip_type,
  lnRR = round(df_box_all$middle, digits = 3),
  lowerCI = round(df_box_all$ymin, digits = 3),
  upperCI = round(df_box_all$ymax, digits = 3),
  ci_range = str_c("(", 
                   round(df_box_all$ymin, digits = 3), ", ", 
                   round(df_box_all$ymax, digits = 3), ")")) %>%
  filter(trait %in% c("lma", "leaf_n_mass", "leaf_n_area",
                      "leaf_p_mass", "leaf_p_area", "leaf_np",
                      "leaf_pi", "leaf_sugar_p", "leaf_metabolic_p",
                      "leaf_nucleic_p", "leaf_structural_p", "asat", 
                      "vcmax", "jmax", "leaf_pnue", "leaf_ppue", 
                      "leaf_wue", "total_biomass", "agb", "agb_n", 
                      "agb_p", "bgb", "rmf", "rootshoot")) %>%
  mutate(trait = factor(trait, 
                        levels = c("lma", "leaf_n_mass", "leaf_n_area",
                                   "leaf_p_mass", "leaf_p_area", "leaf_np",
                                   "leaf_pi", "leaf_sugar_p", "leaf_metabolic_p",
                                   "leaf_nucleic_p", "leaf_structural_p", "asat", 
                                   "vcmax", "jmax", "leaf_pnue", "leaf_ppue", 
                                   "leaf_wue", "total_biomass", "agb", "agb_n", 
                                   "agb_p", "bgb", "rmf", "rootshoot")),
         sig = ifelse((lowerCI > 0 & upperCI > 0) |
                        (lowerCI < 0 & upperCI < 0), "***",
                      "NS")) %>%
  arrange(trait)

###############################################################################
# Re-run Nmass models, but include moderator variables
###############################################################################

# Nitrogen addition effect on Nmass - full model
nfert_nmass_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = nfert_lnRR %>% 
                                  filter(myvar == "leaf_n_mass" & !is.na(gs_mat)))
summary(nfert_nmass_fullModel)

# Nfert: Nmass - MAT plot
nfert_nmass_mat_plot <- mod_results(nfert_nmass_fullModel, 
                              mod = "gs_mat",
                              group = "exp", subset = TRUE
                              )$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_n_mass" & 
                             !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of N"["mass"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18)

# Nfert: Nmass - AI plot
nfert_nmass_ai_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_n_mass" & 
                             !is.na(gs_ai)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = expression(bold("ln RR of N"["mass"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Nfert: Nmass - PAR plot
nfert_nmass_par_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_n_mass" & !is.na(gs_par)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of N"["mass"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))


mod_results(nfert_nmass_fullModel, 
            mod = "photo_path",
            group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = name, y = estimate)) +
  geom_point(data = subset(nfert_lnRR, 
                           myvar == "leaf_n_mass" & !is.na(gs_par) & 
                             !is.na(photo_path)),
             aes(x = toupper(photo_path), y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL), width = 0.1) +
  geom_point(aes(fill = name), size = 7, shape = 21) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "Photosynthetic pathway", 
       y = expression(bold("ln RR of N"["mass"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  coord_flip() +
  theme_classic(base_size = 18)

# Phosphorus addition effect on Nmass - full model
pfert_nmass_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = pfert_lnRR %>% 
                                  filter(myvar == "leaf_n_mass" & 
                                           !is.na(gs_mat)))
summary(pfert_nmass_fullModel)

# Pfert: Nmass - MAT plot
pfert_nmass_mat_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_n_mass" & 
                             !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of N"["mass"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18)

# Pfert: Nmass - AI plot
pfert_nmass_ai_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_n_mass" & !is.na(gs_ai)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = expression(bold("ln RR of N"["mass"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Pfert: Nmass - PAR plot
pfert_nmass_par_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_n_mass" & !is.na(gs_par)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of N"["mass"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+P addition effect on Nmass - full model
npfert_nmass_fullModel <- rma.mv(logr, 
                                 logr_var,
                                 method = "REML", 
                                 random = ~ 1 | exp, 
                                 mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                 slab = exp, control = list(stepadj = 0.3), 
                                 data = npfert_lnRR %>% 
                                   filter(myvar == "leaf_n_mass" & !is.na(gs_mat)))
summary(npfert_nmass_fullModel)

# N+P fert: Nmass - MAT plot
npfert_nmass_mat_plot <- mod_results(npfert_nmass_fullModel, 
                                     mod = "gs_mat",
                                     group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_n_mass" & 
                             !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "magenta") +
  geom_smooth(method = "loess", linewidth = 2, color = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of N"["mass"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18)

# N+P fert: Nmass - AI plot
npfert_nmass_ai_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_n_mass" & 
                             !is.na(gs_ai)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = expression(bold("ln RR of N"["mass"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+P fert: Nmass - PAR plot
npfert_nmass_par_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_n_mass" & 
                             !is.na(gs_par)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of N"["mass"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

png("../plots/CNPmeta_Nmass_clim_moderators.png", units = "in",
    height = 16, width = 16, res = 600)
ggarrange(nfert_nmass_mat_plot, nfert_nmass_ai_plot, nfert_nmass_par_plot,
          pfert_nmass_mat_plot, pfert_nmass_ai_plot, pfert_nmass_par_plot,
          npfert_nmass_mat_plot, npfert_nmass_ai_plot, npfert_nmass_par_plot,
          nrow = 3, ncol = 3, common.legend = TRUE, hjust = 0,
          legend = "bottom",
          labels = c("(a)", "(b)", "(c)", 
                     "(d)", "(e)", "(f)",
                     "(g)", "(h)", "(i)"),
          font.label = list(size = 18))
dev.off()

###############################################################################
# Re-run Pmass models, but include moderator variables
###############################################################################

# Nitrogen addition effect on Pmass - full model
nfert_pmass_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = nfert_lnRR %>% 
                                  filter(myvar == "leaf_p_mass" & !is.na(gs_mat)))
summary(nfert_pmass_fullModel)

# Nfert: Pmass - MAT plot
nfert_pmass_mat_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_p_mass" &
                             !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-1, 0.6), breaks = seq(-1, 0.6, 0.4)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of P"["mass"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18)

# Nfert: Pmass - AI plot
nfert_pmass_ai_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-1, 0.6), breaks = seq(-1, 0.6, 0.4)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = expression(bold("ln RR of P"["mass"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Nfert: Pmass - PAR plot
nfert_pmass_par_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-1, 0.6), breaks = seq(-1, 0.6, 0.4)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["mass"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Phosphorus addition effect on Pmass - full model
pfert_pmass_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = pfert_lnRR %>% 
                                  filter(myvar == "leaf_p_mass" & !is.na(gs_mat)))

summary(pfert_pmass_fullModel)

# Pfert: Pmass - MAT plot
pfert_pmass_mat_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of P"["mass"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18)

# Pfert: Pmass - AI plot
pfert_pmass_ai_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = expression(bold("ln RR of P"["mass"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Pfert: Pmass - PAR plot
pfert_pmass_par_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["mass"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))


# N+P addition effect on Pmass - full model
npfert_pmass_fullModel <- rma.mv(logr, 
                                 logr_var,
                                 method = "REML", 
                                 random = ~ 1 | exp, 
                                 mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                 slab = exp, control = list(stepadj = 0.3), 
                                 data = npfert_lnRR %>% 
                                   filter(myvar == "leaf_p_mass" & !is.na(gs_mat)))
summary(npfert_pmass_fullModel)

# N+P fert: Pmass - MAT plot
npfert_pmass_mat_plot <- mod_results(npfert_pmass_fullModel, 
                                     mod = "gs_mat",
                                     group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_mass" & 
                             !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "magenta") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "magenta", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of P"["mass"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18)

# N+P fert: Pmass - AI plot
npfert_pmass_ai_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = expression(bold("ln RR of P"["mass"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+P fert: Pmass - PAR plot
npfert_pmass_par_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_mass"),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["mass"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Make plot
png("../plots/CNPmeta_Pmass_clim_moderators.png", units = "in",
    height = 14, width = 16, res = 600)
ggarrange(nfert_pmass_mat_plot, nfert_pmass_ai_plot, nfert_pmass_par_plot,
          pfert_pmass_mat_plot, pfert_pmass_ai_plot, pfert_pmass_par_plot,
          npfert_pmass_mat_plot, npfert_pmass_ai_plot, npfert_pmass_par_plot,
          nrow = 3, ncol = 3, common.legend = TRUE, legend = "bottom", hjust = 0,
          labels = c("(a)", "(b)", "(c)", 
                     "(d)", "(e)", "(f)",
                     "(g)", "(h)", "(i)"),
          font.label = list(size = 18))
dev.off()

###############################################################################
# Re-run leaf_np models, but include moderator variables
###############################################################################

# Nitrogen addition effect on Pmass - full model
nfert_leafnp_fullModel <- rma.mv(logr, 
                                 logr_var,
                                 method = "REML", 
                                 random = ~ 1 | exp, 
                                 mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                 slab = exp, control = list(stepadj = 0.3), 
                                 data = nfert_lnRR %>% 
                                   filter(myvar == "leaf_np" & !is.na(gs_mat)))
summary(nfert_leafnp_fullModel)

# Nfert: Leaf N:P - MAT plot
nfert_leafnp_mat_plot <- mod_results(nfert_leafnp_fullModel,
                                     mod = "gs_mat", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_np" & 
                             !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = "ln RR of leaf N:P to N addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Nfert: Leaf N:P - AI plot
nfert_leafnp_ai_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_np" & 
                             !is.na(gs_mat)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = "ln RR of leaf N:P to N addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Nfert: Leaf N:P - PAR plot
nfert_leafnp_par_plot <- mod_results(nfert_leafnp_fullModel,
                                     mod = "gs_par", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_np" & 
                             !is.na(gs_mat)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       y = "ln RR of leaf N:P to N addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Phosphorus addition effect on leaf N:P - full model
pfert_leafnp_fullModel <- rma.mv(logr, 
                                 logr_var,
                                 method = "REML", 
                                 random = ~ 1 | exp, 
                                 mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                 slab = exp, control = list(stepadj = 0.3), 
                                 data = pfert_lnRR %>% 
                                   filter(myvar == "leaf_np" & !is.na(gs_mat)))
summary(pfert_leafnp_fullModel)

# Pfert: Leaf N:P - MAT plot
pfert_leafnp_mat_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_np"),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 25), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-2, 0.5), breaks = seq(-2, 0.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = "ln RR of leaf N:P to P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Pfert: Leaf N:P - AI plot
pfert_leafnp_ai_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_np"),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-2, 0.5), breaks = seq(-2, 0.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = "ln RR of leaf N:P to P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Pfert: Leaf N:P - PAR plot
pfert_leafnp_par_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_np"),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-2, 0.5), breaks = seq(-2, 0.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       y = expression(bold("ln RR of leaf N:P to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))


# N+P addition effect on leaf N:P - full model
npfert_leafnp_fullModel <- rma.mv(logr, 
                                  logr_var,
                                  method = "REML", 
                                  random = ~ 1 | exp, 
                                  mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                  slab = exp, control = list(stepadj = 0.3), 
                                  data = npfert_lnRR %>% 
                                    filter(myvar == "leaf_np" & !is.na(gs_mat)))
summary(npfert_leafnp_fullModel)

# N+P fert: Leaf N:P - MAT plot
npfert_leafnp_mat_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_np"),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(-1.5, 0.6), breaks = seq(-1.5, 0.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = "ln RR of leaf N:P to N+P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# N+P fert: Leaf N:P - AI plot
npfert_leafnp_ai_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_np"),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-1.5, 0.6), breaks = seq(-1.5, 0.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = "ln RR of leaf N:P to N+P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+P fert: Leaf N:P - PAR plot
npfert_leafnp_par_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_np" & 
                             !is.na(gs_mat)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-1.5, 0.6), breaks = seq(-1.5, 0.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = "ln RR of leaf N:P to N+P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Interaction effect on Leaf N:P - full model
int_leafnp_fullModel <- rma.mv(yi = dNPi,
                               V = dvNPi,
                               W = dwNPi,
                               method = "REML", 
                               random = ~ 1 | exp, 
                               mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                               slab = exp, control = list(stepadj = 0.3), 
                               data = CNP_effect_sizes_reduced %>% 
                                 filter(myvar == "leaf_np" & !is.na(gs_mat)))
summary(int_leafnp_fullModel)

# Interaction: Pmass - MAT plot
int_leafnp_mat_plot <- mod_results(int_leafnp_fullModel,
                                   mod = "gs_mat",
                                   group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(CNP_effect_sizes_reduced, myvar == "leaf_np"),
             aes(x = gs_mat, y = dNPi, size = 1/dvNPi), 
             alpha = 0.30) +
  geom_ribbon(aes(ymin = lowerCL, ymax = upperCL),
              alpha = 0.3, fill = "black") +
  geom_smooth(method = "loess", color = "black", linewidth = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = "Interaction effect size of leaf N:P",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Interaction: Pmass - AI plot
int_leafnp_ai_plot <- ggplot() +
  geom_point(data = subset(CNP_effect_sizes_reduced, myvar == "leaf_np"),
             aes(x = gs_ai, y = dNPi, size = 1/dvNPi), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +  
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = "Interaction effect size of leaf N:P",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Interaction: Pmass - PAR plot
int_leafnp_par_plot <- mod_results(int_leafnp_fullModel, 
                                   mod = "gs_par", 
                                   group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(CNP_effect_sizes_reduced, myvar == "leaf_np" & 
                             !is.na(gs_mat)),
             aes(x = gs_par, y = dNPi, size = 1/dvNPi), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "black") +
  geom_smooth(method = "loess", linewidth = 2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       y = "Interaction effect size of leaf N:P",       
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Make plot
png("../plots/CNPmeta_leafNP_clim_moderators.png", units = "in",
    height = 20, width = 20, res = 600)
ggarrange(nfert_leafnp_mat_plot, nfert_leafnp_ai_plot, nfert_leafnp_par_plot,
          pfert_leafnp_mat_plot, pfert_leafnp_ai_plot, pfert_leafnp_par_plot,
          npfert_leafnp_mat_plot, npfert_leafnp_ai_plot, npfert_leafnp_par_plot,
          int_leafnp_mat_plot, int_leafnp_ai_plot, int_leafnp_par_plot,
          nrow = 4, ncol = 3, common.legend = TRUE, legend = "bottom", hjust = 0,
          labels = c("(a)", "(b)", "(c)", 
                     "(d)", "(e)", "(f)",
                     "(g)", "(h)", "(i)",
                     "(j)", "(k)", "(l)"),
          font.label = list(size = 18))
dev.off()

###############################################################################
# Re-run Narea models, but include moderator variables
###############################################################################

# Nitrogen addition effect on Narea - full model
nfert_narea_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = nfert_lnRR %>% 
                                  filter(myvar == "leaf_n_area" & !is.na(gs_mat)))
summary(nfert_narea_fullModel)

# Nfert: Narea - MAT plot
nfert_narea_mat_plot <- mod_results(nfert_narea_fullModel,
                                    mod = "gs_mat", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_n_area" & 
                             !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of N"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Nfert: Narea - AI plot
nfert_narea_ai_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_n_area" & 
                             !is.na(gs_mat)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = expression(bold("ln RR of N"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Nfert: Narea - PAR plot
nfert_narea_par_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_n_area" & !is.na(gs_mat)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       y = expression(bold("ln RR of N"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Phosphorus addition effect on Narea - full model
pfert_narea_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = pfert_lnRR %>% 
                                  filter(myvar == "leaf_n_area" & !is.na(gs_mat)))
summary(pfert_narea_fullModel)

# Pfert: Narea - MAT plot
pfert_narea_mat_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_n_area" & 
                             !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of N"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Pfert: Narea - AI plot
pfert_narea_ai_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_n_area" & 
                             !is.na(gs_mat)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = expression(bold("ln RR of N"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Pfert: Narea - PAR plot
pfert_narea_par_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_n_area" & 
                             !is.na(gs_mat)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       y = expression(bold("ln RR of N"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+P addition effect on Narea - full model
npfert_narea_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = npfert_lnRR %>% 
                                  filter(myvar == "leaf_n_area" & !is.na(gs_mat)))
summary(npfert_narea_fullModel)

# N+Pfert: Narea - MAT plot
npfert_narea_mat_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_n_area" & 
                             !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of N"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# N+Pfert: Narea - AI plot
npfert_narea_ai_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_n_area" & 
                             !is.na(gs_mat)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = expression(bold("ln RR of N"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+Pfert: Narea - PAR plot
npfert_narea_par_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_n_area" & !is.na(gs_mat)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       y = expression(bold("ln RR of N"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Make plot
png("../plots/CNPmeta_Narea_clim_moderators.png", units = "in",
    height = 14, width = 16, res = 600)
ggarrange(nfert_narea_mat_plot, nfert_narea_ai_plot, nfert_narea_par_plot,
          pfert_narea_mat_plot, pfert_narea_ai_plot, pfert_narea_par_plot,
          npfert_narea_mat_plot, npfert_narea_ai_plot, npfert_narea_par_plot,
          nrow = 3, ncol = 3, common.legend = TRUE, legend = "bottom", hjust = 0,
          labels = c("(a)", "(b)", "(c)", 
                     "(d)", "(e)", "(f)",
                     "(g)", "(h)", "(i)"),
          font.label = list(size = 18))
dev.off()

###############################################################################
# Re-run Parea models, but include moderator variables
###############################################################################

# Nitrogen addition effect on Parea - full model
nfert_parea_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = nfert_lnRR %>% 
                                  filter(myvar == "leaf_p_area" & !is.na(gs_mat)))
summary(nfert_parea_fullModel)

# Nfert: Parea - MAT plot
nfert_parea_mat_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_p_area" & !
                             is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of P"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# NPfert: Parea - AI plot
nfert_parea_ai_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_p_area" & 
                             !is.na(gs_mat)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless",
       y = expression(bold("ln RR of P"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Nfert: Parea - PAR plot
nfert_parea_par_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "leaf_p_area" & 
                             !is.na(gs_mat)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Phosphorus addition effect on Parea - full model
pfert_parea_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = pfert_lnRR %>% 
                                  filter(myvar == "leaf_p_area" & !is.na(gs_mat) & logr > -0.5))
summary(pfert_parea_fullModel)

# Pfert: Parea - MAT plot
pfert_parea_mat_plot <- mod_results(pfert_parea_fullModel,
                                    mod = "gs_mat", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_area" 
                           & logr > -0.5 & !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of P"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Pfert: Parea - AI plot
pfert_parea_ai_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_area" & 
                             !is.na(gs_mat)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = expression(bold("ln RR of P"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Pfert: Parea - PAR plot
pfert_parea_par_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "leaf_p_area" & 
                             logr > -0.5 & !is.na(gs_mat)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.6, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+P addition effect on Parea - full model
npfert_parea_fullModel <- rma.mv(logr, 
                                 logr_var,
                                 method = "REML", 
                                 random = ~ 1 | exp, 
                                 mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                 slab = exp, control = list(stepadj = 0.3), 
                                 data = npfert_lnRR %>% 
                                   filter(myvar == "leaf_p_area"  & !is.na(gs_mat)))
summary(npfert_parea_fullModel)

# N+Pfert: Parea - MAT plot
npfert_parea_mat_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_area" & 
                             !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of P"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# N+Pfert: Parea - AI plot
npfert_parea_ai_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_area" & 
                             !is.na(gs_mat)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless",
       y = expression(bold("ln RR of P"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+Pfert: Parea - PAR plot
npfert_parea_par_plot <- mod_results(npfert_parea_fullModel,
                                    mod = "gs_par", group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(npfert_lnRR, myvar == "leaf_p_area" & 
                             !is.na(gs_mat)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "magenta") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of P"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Make plot
png("../plots/CNPmeta_Parea_clim_moderators.png", units = "in",
    height = 16, width = 16, res = 600)
ggarrange(nfert_parea_mat_plot, nfert_parea_ai_plot, nfert_parea_par_plot,
          pfert_parea_mat_plot, pfert_parea_ai_plot, pfert_parea_par_plot,
          npfert_parea_mat_plot, npfert_parea_ai_plot, npfert_parea_par_plot,
          nrow = 3, ncol = 3, common.legend = TRUE, legend = "bottom", hjust = 0,
          labels = c("(a)", "(b)", "(c)", 
                     "(d)", "(e)", "(f)",
                     "(g)", "(h)", "(i)"),
          font.label = list(size = 18))
dev.off()

###############################################################################
# Re-run Marea models, but include moderator variables
###############################################################################

# Nitrogen addition effect on Marea - full model
nfert_marea_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = nfert_lnRR %>% 
                                  filter(myvar == "lma" & !is.na(gs_mat)))
summary(nfert_marea_fullModel)

# Nfert: Marea - MAT plot
nfert_marea_mat_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "lma" & 
                             !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.8, 0.8), breaks = seq(-0.8, 0.8, 0.4)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of M"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# NPfert: Marea - AI plot
nfert_marea_ai_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "lma" & 
                             !is.na(gs_mat)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.8, 0.8), breaks = seq(-0.8, 0.8, 0.4)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = expression(bold("ln RR of M"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Nfert: Marea - PAR plot
nfert_marea_par_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "lma" & 
                             !is.na(gs_mat)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.8, 0.8), breaks = seq(-0.8, 0.8, 0.4)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of M"["area"]*" to N addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Phosphorus addition effect on Marea - full model
pfert_marea_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = pfert_lnRR %>% 
                                  filter(myvar == "lma" & !is.na(gs_mat) & logr > -0.5))
summary(pfert_marea_fullModel)

# Pfert: Marea - MAT plot
pfert_marea_mat_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "lma" & 
                             !is.na(gs_mat) & logr > -0.5),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.4, 0.5), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of M"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Pfert: Parea - AI plot
pfert_marea_ai_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "lma" & 
                             !is.na(gs_mat)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.4, 0.5), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = expression(bold("ln RR of M"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Pfert: Parea - PAR plot
pfert_marea_par_plot <-  ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "lma" & 
                             !is.na(gs_mat)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.4, 0.5), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of M"["area"]*" to P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+P addition effect on Parea - full model
npfert_marea_fullModel <- rma.mv(logr, 
                                 logr_var,
                                 method = "REML", 
                                 random = ~ 1 | exp, 
                                 mods = ~ gs_mat + gs_ai + gs_par + n_fixer + photo_path + myc_assoc,
                                 slab = exp, control = list(stepadj = 0.3), 
                                 data = npfert_lnRR %>% 
                                   filter(myvar == "lma" & logr > -0.5 & !is.na(gs_mat)))
summary(npfert_marea_fullModel)

# N+Pfert: Parea - MAT plot
npfert_marea_mat_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "lma" & 
                             !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, 0.25)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("ln RR of M"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# N+Pfert: Parea - AI plot
npfert_marea_ai_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "lma" & 
                             !is.na(gs_mat)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, 0.25)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = expression(bold("ln RR of M"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+Pfert: Parea - PAR plot
npfert_marea_par_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "lma" & 
                             !is.na(gs_mat)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1010), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, 0.25)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = expression(bold("ln RR of M"["area"]*" to N+P addition")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Make plot
png("../plots/CNPmeta_Marea_clim_moderators.png", units = "in",
    height = 14, width = 16, res = 600)
ggarrange(nfert_marea_mat_plot, nfert_marea_ai_plot, nfert_marea_par_plot,
          pfert_marea_mat_plot, pfert_marea_ai_plot, pfert_marea_par_plot,
          npfert_marea_mat_plot, npfert_marea_ai_plot, npfert_marea_par_plot,
          nrow = 3, ncol = 3, common.legend = TRUE, hjust = 0, vjust = 0,
          labels = c("(a)", "(b)", "(c)", 
                     "(d)", "(e)", "(f)",
                     "(g)", "(h)", "(i)"),
          font.label = list(size = 18))
dev.off()

###############################################################################
# Re-run AGB models, but include moderator variables
###############################################################################

# Nitrogen addition effect on Marea - full model
nfert_agb_fullModel <- rma.mv(logr, 
                              logr_var,
                              method = "REML", 
                              random = ~ 1 | exp, 
                              mods = ~ gs_mat + gs_ai + gs_par,
                              slab = exp, control = list(stepadj = 0.3), 
                              data = nfert_lnRR %>% 
                                filter(myvar == "agb" & !is.na(gs_mat)))
summary(nfert_agb_fullModel)

# Nfert: AGB - MAT plot
nfert_agb_mat_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "agb"),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = "ln RR of AGB to N addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Nfert: AGB - AI plot
nfert_agb_ai_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "agb"& 
                             !is.na(gs_mat)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = "ln RR of AGB to N addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Nfert: AGB - PAR plot
nfert_agb_par_plot <- ggplot() +
  geom_point(data = subset(nfert_lnRR, myvar == "agb"& 
                             !is.na(gs_mat)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1210), breaks = seq(500, 1200, 100)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = "ln RR of AGB to N addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Phosphorus addition effect on aboveground biomass - full model
pfert_agb_fullModel <- rma.mv(logr, 
                              logr_var,
                              method = "REML", 
                              random = ~ 1 | exp, 
                              mods = ~ gs_mat + gs_ai + gs_par,
                              slab = exp, control = list(stepadj = 0.3), 
                              data = pfert_lnRR %>% 
                                filter(myvar == "agb" & !is.na(gs_mat) & logr < 2))
summary(pfert_agb_fullModel)

# Pfert: AGB - MAT plot
pfert_agb_mat_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "agb"& 
                             !is.na(gs_mat) & logr < 2),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = "ln RR of AGB to P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Pfert: AGB - AI plot
pfert_agb_ai_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "agb"& 
                             !is.na(gs_mat)),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = "ln RR of AGB to P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Pfert: AGB - PAR plot
pfert_agb_par_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "agb"& 
                             !is.na(gs_mat)),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1210), breaks = seq(500, 1200, 100)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = "ln RR of AGB to P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# N+P addition effect on aboveground biomass - full model
npfert_agb_fullModel <- rma.mv(logr, 
                               logr_var,
                               method = "REML", 
                               random = ~ 1 | exp, 
                               mods = ~ gs_mat + gs_ai + gs_par,
                               slab = exp, control = list(stepadj = 0.3), 
                               data = npfert_lnRR %>% 
                                 filter(myvar == "agb" & !is.na(gs_mat) & logr < 2))
summary(npfert_agb_fullModel)

# Pfert: AGB - MAT plot
npfert_agb_mat_plot <- ggplot() +
  geom_point(data = subset(npfert_lnRR, myvar == "agb"& 
                             !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = "ln RR of AGB to N+P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Pfert: AGB - AI plot
npfert_agb_ai_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "agb"& 
                             !is.na(gs_mat) & logr < 2),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = "ln RR of AGB to N+P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))

# Pfert: AGB - PAR plot
npfert_agb_par_plot <- ggplot() +
  geom_point(data = subset(pfert_lnRR, myvar == "agb"& 
                             !is.na(gs_mat) &  logr < 2),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1210), breaks = seq(500, 1200, 100)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = "ln RR of AGB to N+P addition",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Interaction effect on aboveground biomass - full model
int_agb_fullModel <- rma.mv(yi = dNPi,
                            V = dvNPi,
                            W = dwNPi,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ gs_mat + gs_ai + gs_par,
                            slab = exp, control = list(stepadj = 0.3), 
                            data = CNP_effect_sizes_reduced %>% 
                              filter(myvar == "agb" & !is.na(gs_mat) & dNPi < 2.5 & dNPi > -1.5))

summary(int_agb_fullModel)


# Interaction: AGB - MAT plot
int_agb_mat_plot <- mod_results(int_agb_fullModel, 
                                           mod = "gs_mat", 
                                           group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_point(data = subset(CNP_effect_sizes_reduced, 
                           myvar == "agb" & !is.na(gs_mat) & 
                             dNPi < 2.5 & dNPi > -2.5),
             aes(x = gs_mat, y = dNPi, size = 1/dvNPi), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "black") +
  geom_smooth(method = "loess", linewidth = 2, color = "black", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-1.5, 2), breaks = seq(-1.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = "Interaction effect size of AGB",       
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Interaction: AGB - AI plot
int_agb_ai_plot <- ggplot() +
  geom_point(data = subset(CNP_effect_sizes_reduced, 
                           myvar == "agb" & !is.na(gs_mat) & 
                             dNPi < 2.5 & dNPi > -2.5),
             aes(x = gs_ai, y = dNPi, size = 1/dvNPi), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-1.5, 2), breaks = seq(-1.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = "AI (unitless)",
       y = "Interaction effect size of AGB",       
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Interaction: AGB - PAR plot
int_agb_par_plot <- ggplot() +
  geom_point(data = subset(CNP_effect_sizes_reduced, 
                           myvar == "agb" & !is.na(gs_mat) & 
                             dNPi < 2.5 & dNPi > -2.5),
             aes(x = gs_par, y = dNPi, size = 1/dvNPi), 
             alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(490, 1210), breaks = seq(500, 1200, 100)) +
  scale_y_continuous(limits = c(-1.5, 2), breaks = seq(-1.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("I"["L"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = "Interaction effect size of AGB",       
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 18) + 
  theme(axis.title = element_text(face = "bold"))

# Make plot
png("../plots/CNPmeta_AGB_clim_moderators.png", units = "in",
    height = 20, width = 20, res = 600)
ggarrange(nfert_agb_mat_plot, nfert_agb_ai_plot, nfert_agb_par_plot,
          pfert_agb_mat_plot, pfert_agb_ai_plot, pfert_agb_par_plot,
          npfert_agb_mat_plot, npfert_agb_ai_plot, npfert_agb_par_plot,
          int_agb_mat_plot, int_agb_ai_plot, int_agb_par_plot,
          nrow = 4, ncol = 3, common.legend = TRUE, hjust = 0, vjust = 0,
          labels = c("(a)", "(b)", "(c)", 
                     "(d)", "(e)", "(f)",
                     "(g)", "(h)", "(i)"),
          font.label = list(size = 18))
dev.off()

###############################################################################
# Some additional exploratory models
###############################################################################

# No effect of climate demand on total biomass

# BGB response to P fertilization increases with AI (wetter = stronger response)
# and PAR:
pfert_bgb_fullModel <- rma.mv(logr, 
                               logr_var,
                               method = "REML", 
                               random = ~ 1 | exp, 
                               mods = ~ mat + ai + par + n_fixer + photo_path + myc_assoc,
                               slab = exp, control = list(stepadj = 0.3), 
                               data = pfert_lnRR %>% 
                                 filter(myvar == "bgb" & !is.na(mat)))
summary(pfert_bgb_fullModel)

# BGB response to N+P fertilization decreases with MAT and strongly increases
# with increasing AI (wetter = stronger response) and with increasing PAR:
npfert_bgb_fullModel <- rma.mv(logr, 
                                logr_var,
                                method = "REML", 
                                random = ~ 1 | exp, 
                                mods = ~ mat + ai + par + n_fixer + photo_path + myc_assoc,
                                slab = exp, control = list(stepadj = 0.3), 
                                data = npfert_lnRR %>% 
                                  filter(myvar == "bgb" & !is.na(mat)))
summary(nfert_bgb_fullModel)

# RMF response to N fertilization decreases with increasing MAT and increases
# with increasing PAR (all marginal)
nfert_rmf_fullModel <- rma.mv(logr, 
                               logr_var,
                               method = "REML", 
                               random = ~ 1 | exp, 
                               mods = ~ mat + ai + par + n_fixer + photo_path + myc_assoc,
                               slab = exp, control = list(stepadj = 0.3), 
                               data = nfert_lnRR %>% 
                                 filter(myvar == "rmf" & !is.na(mat)))
summary(nfert_rmf_fullModel)



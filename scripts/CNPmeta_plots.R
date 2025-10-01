# Libraries
library(tidyverse)
library(metafor)
library(ggpubr)
library(naniar) # to resolve NA/<NA> issue
library(orchaRd)

# Load meta-analysis results
meta_results <- read.csv("../data/CNPmeta_logr_results.csv")
meta_results_int <- read.csv("../data/CNPmeta_logr_results_int.csv")

# Load meta-analysis confidence intervals
meta_ci <- read.csv("../data/CNPmeta_ci.csv") %>%
  filter(var %in% c("leaf_np", "leaf_p_area", "leaf_p_mass", "leaf_n_area", "leaf_n_mass",
                    "lma", "leaf_ppue", "leaf_pnue", "jmax_vcmax", "jmax", "vcmax", "asat",
                    "rootshoot", "rmf", "bgb", "agb", "total_biomass")) %>%
  mutate(var = factor(var, levels = c("leaf_np", "leaf_p_area", "leaf_p_mass", 
                                      "leaf_n_area", "leaf_n_mass",
                                      "lma", "leaf_ppue", "leaf_pnue", 
                                      "jmax_vcmax", "jmax", "vcmax", "asat", "rootshoot",
                                      "rmf", "bgb", "agb", "total_biomass")))

meta_ci_int <- read.csv("../data/CNPmeta_ci_int.csv") %>%
  filter(var %in% c("leaf_np", "leaf_p_area", "leaf_p_mass", "leaf_n_area", "leaf_n_mass",
                    "lma", "leaf_ppue", "leaf_pnue", "jmax", "vcmax", "asat",
                    "rootshoot", "rmf", "bgb", "agb", "total_biomass")) %>%
  mutate(var = factor(var, levels = c("leaf_np", "leaf_p_area", "leaf_p_mass", 
                                      "leaf_n_area", "leaf_n_mass",
                                      "lma", "leaf_ppue", "leaf_pnue", 
                                      "jmax", "vcmax", "asat", "rootshoot",
                                      "rmf", "bgb", "agb", "total_biomass")))

# Load species moderator results
meta_photo_results <- read.csv("../data/CNPmeta_photo_moderators.csv")
meta_nfix_results <- read.csv("../data/CNPmeta_nfix_moderators.csv")
meta_myc_results <- read.csv("../data/CNPmeta_myc_moderators.csv")

# Check file structure
head(meta_results)
head(meta_results_int)
head(meta_ci)
head(meta_ci_int)
head(meta_photo_results)
head(meta_nfix_results)
head(meta_myc_results)

#####################################################################
# N addition plots
#####################################################################
nadd_chemistry_plot <- ggplot(data = meta_ci %>% 
                                drop_na(var) %>% 
                                filter(manip_type == "n") %>%
                                filter(var %in% c("lma", 
                                                  "leaf_n_mass", 
                                                  "leaf_n_area", 
                                                  "leaf_p_mass", 
                                                  "leaf_p_area",
                                                  "leaf_np")),
                              aes(x = var, y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), size = 1, width = 0.25) +
  geom_point(size = 4, fill = "red", shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf N:P",
                              expression("P"["area"]), 
                              expression("P"["mass"]), 
                              expression("N"["area"]),
                              expression("N"["mass"]),
                              expression("M"["area"]))) +
  geom_text(aes(label = k_sig), y = 135, fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(-50, 150), breaks = seq(-50, 150, 50)) +
  labs(x = "", 
       y = NULL) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 18, color = "black"))
nadd_chemistry_plot

nadd_photo_plot <- ggplot(data = meta_ci %>% 
                            drop_na(var) %>% 
                            filter(manip_type == "n") %>%
                            filter(var %in% c("asat", 
                                              "vcmax", 
                                              "jmax",
                                              "leaf_pnue", 
                                              "leaf_ppue")),
                          aes(x = var, y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), size = 1, width = 0.25) +
  geom_point(size = 4, fill = "red", shape = 21) +
  geom_text(aes(label = k_sig), y = 85, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("PPUE",
                              "PNUE",
                              expression("J"["max"]),
                              expression("V"["cmax"]),
                              expression("A"["sat"]))) +
  scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, 50)) +
  scale_fill_manual(values = "red") +
  labs(x = "", y = "Response to N addition (%)") +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 18, color = "black"))
nadd_photo_plot

nadd_bio_plot <- ggplot(data = meta_ci %>% 
                          drop_na(var) %>% 
                          filter(manip_type == "n") %>%
                          filter(var %in% c("rootshoot", 
                                            "rmf", 
                                            "bgb", 
                                            "agb", 
                                            "total_biomass")),
                        aes(x = var, y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), size = 1, width = 0.25) +
  geom_point(size = 4, fill = "red", shape = 21) +
  geom_text(aes(label = k_sig), y = 135, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Root:shoot",
                              "RMF",
                              "BGB",
                              "AGB",
                              "Total biomass")) +
  scale_y_continuous(limits = c(-50, 150), breaks = seq(-50, 150, 50)) +
  labs(x = "", 
       y = NULL)  +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 18, color = "black"))
nadd_bio_plot

#####################################################################
# P addition plots
#####################################################################
padd_chemistry_plot <- ggplot(data = meta_ci %>% 
                                drop_na(var) %>% 
                                filter(manip_type == "p") %>%
                                filter(var %in% c("lma", 
                                                  "leaf_n_mass", 
                                                  "leaf_n_area", 
                                                  "leaf_p_mass", 
                                                  "leaf_p_area",
                                                  "leaf_np")),
                              aes(x = var, y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), size = 1, width = 0.25) +
  geom_point(size = 4, fill = "blue", shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf N:P",
                              expression("P"["area"]),
                              expression("P"["mass"]),
                              expression("N"["area"]),
                              expression("N"["mass"]),
                              expression("M"["area"]))) +
  geom_text(aes(label = k_sig), y = 135, fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(-50, 150), breaks = seq(-50, 150, 50)) +
  labs(x = "", 
       y = NULL) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 18, color = "black"))
padd_chemistry_plot

padd_photo_plot <- ggplot(data = meta_ci %>% 
                            drop_na(var) %>% 
                            filter(manip_type == "p") %>%
                            filter(var %in% c("asat", 
                                              "vcmax", 
                                              "jmax",
                                              "leaf_pnue", 
                                              "leaf_ppue")),
                          aes(x = var, y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), size = 1, width = 0.25) +
  geom_point(size = 4, fill = "blue", shape = 21) +
  geom_text(aes(label = k_sig), y = 85, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("PPUE",
                              "PNUE",
                              expression("J"["max"]),
                              expression("V"["cmax"]),
                              expression("A"["sat"]))) +
  scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, 50)) +
  scale_fill_manual(values = "blue") +
  labs(x = "", y = "Response to P addition (%)") +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 18, color = "black"))
padd_photo_plot

padd_bio_plot <- ggplot(data = meta_ci %>% 
                          drop_na(var) %>% 
                          filter(manip_type == "p") %>%
                          filter(var %in% c("rootshoot", 
                                            "rmf", 
                                            "bgb", 
                                            "agb", 
                                            "total_biomass")),
                        aes(x = var, y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), 
                size = 1, width = 0.25) +
  geom_point(size = 4, fill = "blue", shape = 21) +
  geom_text(aes(label = k_sig), y = 135, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Root:shoot",
                              "RMF",
                              "BGB",
                              "AGB",
                              "Total biomass")) +
  scale_y_continuous(limits = c(-50, 150), breaks = seq(-50, 150, 50)) +
  labs(x = "", 
       y = NULL)  +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 18, color = "black"))
padd_bio_plot

#####################################################################
# N+P addition plots
#####################################################################
npadd_chemistry_plot <- ggplot(data = meta_ci %>% 
                                 drop_na(var) %>% 
                                 filter(manip_type == "np") %>%
                                 filter(var %in% c("lma", 
                                                   "leaf_n_mass", 
                                                   "leaf_n_area", 
                                                   "leaf_p_mass", 
                                                   "leaf_p_area",
                                                   "leaf_np")),
                               aes(x = var, y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), 
                size = 1, width = 0.25) +
  geom_point(size = 4, fill = "magenta", shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf N:P",
                              expression("P"["area"]), 
                              expression("P"["mass"]), 
                              expression("N"["area"]),
                              expression("N"["mass"]),
                              expression("M"["area"]))) +
  geom_text(aes(label = k_sig), y = 135, fontface = "bold", size = 5) +
  scale_y_continuous(limits = c(-50, 150), breaks = seq(-50, 150, 50)) +
  labs(x = "", 
       y = NULL) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 18, color = "black"))
npadd_chemistry_plot

npadd_photo_plot <- ggplot(data = meta_ci %>% 
                             drop_na(var) %>% 
                             filter(manip_type == "np") %>%
                             filter(var %in% c("asat", 
                                               "vcmax", 
                                               "jmax",
                                               "leaf_pnue", 
                                               "leaf_ppue")),
                           aes(x = var, y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), size = 1, width = 0.25) +
  geom_point(size = 4, fill = "magenta", shape = 21) +
  geom_text(aes(label = k_sig), y = 85, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("PPUE",
                              "PNUE",
                              expression("J"["max"]),
                              expression("V"["cmax"]),
                              expression("A"["sat"]))) +
  scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, 50)) +
  labs(x = "", y = "Response to N+P addition (%)") +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 18, color = "black"))
npadd_photo_plot

npadd_bio_plot <- ggplot(data = meta_ci %>% 
                           drop_na(var) %>% 
                           filter(manip_type == "np") %>%
                           filter(var %in% c("rootshoot", 
                                             "rmf", 
                                             "bgb", 
                                             "agb", 
                                             "total_biomass")),
                         aes(x = var, y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), size = 1, width = 0.25) +
  geom_point(size = 4, fill = "magenta", shape = 21) +
  geom_text(aes(label = k_sig), y = 135, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Root:shoot",
                              "RMF",
                              "BGB",
                              "AGB",
                              "Total biomass")) +
  scale_y_continuous(limits = c(-50, 150), breaks = seq(-50, 150, 50)) +
  labs(x = "", 
       y = NULL)  +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 18),
        axis.text.y = element_text(size = 18, color = "black"))
npadd_bio_plot

#####################################################################
# N+P interaction plots
#####################################################################
npint_chemistry_plot <- ggplot(data = meta_ci_int %>% 
                                 drop_na(var) %>%
                                 filter(var %in% c("lma",
                                                   "leaf_n_mass", 
                                                   "leaf_n_area", 
                                                   "leaf_p_mass", 
                                                   "leaf_p_area",
                                                   "leaf_np")),
                               aes(x = var, 
                                   y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), 
                size = 1, width = 0.25) +
  geom_point(aes(fill = int_type),size = 4, shape = 21) +
  geom_text(aes(label = k_sig), y = 55, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Leaf N:P",
                              expression("P"["area"]), 
                              expression("P"["mass"]), 
                              expression("N"["area"]),
                              expression("N"["mass"]),
                              expression("M"["area"]))) +
  scale_y_continuous(limits = c(-60, 60), breaks = seq(-60, 60, 30)) +
  scale_fill_manual(values = c("black", "pink", "red")) +
  labs(x = "", 
       y = NULL) +
  coord_flip() +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(color = "black", size = 18))
npint_chemistry_plot

npint_photo_plot <- ggplot(data = meta_ci_int %>% 
                             drop_na(var) %>%
                             filter(var %in% c("asat", 
                                               "vcmax", 
                                               "jmax",
                                               "leaf_pnue", 
                                               "leaf_ppue")),
                           aes(x = factor(var,
                                          levels = c("asat", 
                                                     "vcmax", 
                                                     "jmax",
                                                     "leaf_pnue", 
                                                     "leaf_ppue")), 
                               y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), 
                size = 1, width = 0.25) +
  geom_point(size = 4, fill = "black", shape = 21) +
  geom_text(aes(label = k_sig), y = 185, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("PPUE",
                              "PNUE",
                              expression("J"["max"]),
                              expression("V"["cmax"]),
                              expression("A"["sat"]))) +
  scale_y_continuous(limits = c(-200, 200), breaks = seq(-200, 200, 100)) +
  labs(x = "", y = "") +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(color = "black", size = 18))
npint_photo_plot

npint_bio_plot <- ggplot(data = meta_ci_int %>% 
                           drop_na(var) %>%
                           filter(var %in% c("rootshoot", 
                                             "rmf", 
                                             "bgb", 
                                             "agb", 
                                             "total_biomass")),
                         aes(x = factor(var,
                                        levels = c("rootshoot", 
                                        "rmf", 
                                        "bgb", 
                                        "agb", 
                                        "total_biomass")), 
                             y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), size = 1, width = 0.25) +
  geom_point(aes(fill = int_type), size = 4, shape = 21) +
  geom_text(aes(label = k_sig), y = 75, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("Root:shoot",
                              "RMF",
                              "BGB",
                              "AGB",
                              "Total biomass")) +
  scale_y_continuous(limits = c(-80, 80), breaks = seq(-80, 80, 40)) +
  scale_fill_manual(values = c("black", "red")) +
  labs(y = "Interaction effect (%)", 
       x = "")  +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(color = "black", size = 18)) +
  guides(fill = "none")
npint_bio_plot

#####################################################################
# Nmass -- climate and species identity moderators
#####################################################################

# Nmass climate moderator model
nadd_nmass_clim <- rma.mv(logr, logr_var, method = "REML", 
                          random = ~ 1 | exp, mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, control = list(stepadj = 0.3), 
                          data = meta_results %>% filter(manip_type == "n" & 
                                     myvar == "leaf_n_mass" & 
                                     !is.na(gs_mat) & logr > -0.2))
summary(nadd_nmass_clim)

# Nmass - temperature plot
nmass_tg_plot <- mod_results(nadd_nmass_clim, mod = "gs_mat",
                             group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = (exp(estimate) - 1) * 100)) +
  geom_point(data = subset(meta_results, myvar == "leaf_n_mass" & 
                             manip_type == "n" & !is.na(gs_mat) & 
                             logr > -0.2),
             aes(x = gs_mat, y = (exp(logr) - 1) * 100, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = (exp(upperCL) - 1) * 100, ymin = (exp(lowerCL) - 1) * 100),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-25, 100), breaks = seq(-25, 100, 25)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("N"["mass"]*" response to N addition (%)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
nmass_tg_plot 

# Nmass - aridity plot
nmass_ai_plot <- mod_results(nadd_nmass_clim, mod = "gs_ai",
                             group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = (exp(estimate) - 1) * 100)) +
  geom_point(data = subset(meta_results, myvar == "leaf_n_mass" & 
                             manip_type == "n" & !is.na(gs_ai) & 
                             logr > -0.2),
             aes(x = gs_ai, y = (exp(logr) - 1) * 100, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = (exp(upperCL) - 1) * 100, ymin = (exp(lowerCL) - 1) * 100),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3.4), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-25, 100), breaks = seq(-25, 100, 25)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("AI"["g"]*" (unitless)")),
       y = expression(bold("N"["mass"]*" response to N addition (%)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
nmass_ai_plot 

# Nmass - photosynthetic pathway
nmass_photo_plot <- ggplot(data = meta_photo_results %>% 
                             filter(trait == "nmass" &
                                      nut_add == "n"),
                           aes(x = photo, y = (exp(estimate) - 1) * 100)) +
  geom_errorbar(aes(ymin = (exp(lowerCL) - 1) * 100, 
                    ymax = (exp(upperCL) - 1) * 100),
                width = 0.3, linewidth = 1) +
  geom_point(aes(fill = photo), shape = 21, size = 5) +
  geom_bracket(xmin = "C3", xmax = "C4", y.position = 98, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(x = 1.5, y = 100, label = "*", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("pink", "darkred")) +
  scale_y_continuous(limits = c(-25, 100), breaks = seq(-25, 100, 25)) +
  labs(x = "Photosynthetic pathway",
       y = expression(bold("N"["mass"]*" response to N addition"))) +
  guides(fill = "none") +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
nmass_photo_plot

# Nmass - N fixation
nmass_nfix_plot <- ggplot(data = meta_nfix_results %>% 
                             filter(trait == "nmass" &
                                      nut_add == "n"),
                           aes(x = nfix, y = (exp(estimate) - 1) * 100)) +
  geom_errorbar(aes(ymin = (exp(lowerCL) - 1) * 100, ymax = (exp(upperCL) - 1) * 100),
                width = 0.3, linewidth = 1) +
  geom_point(aes(fill = nfix), shape = 21, size = 5) +
  geom_bracket(xmin = "No", xmax = "Yes", y.position = 98, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(x = 1.5, y = 100, label = "***", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("pink", "darkred")) +
  scale_x_discrete(labels = c("non-fixer", expression("N"["2"]*"-fixer"))) +
  scale_y_continuous(limits = c(-25, 100), breaks = seq(-25, 100, 25)) +
  labs(x = expression(bold("N"["2"]*" fixation ability")),
       y = expression(bold("N"["mass"]*" response to N addition (%)"))) +
  guides(fill = "none") +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
nmass_nfix_plot

# Nmass - mycorrhizal acquisition strategy
nmass_myc_plot <- ggplot(data = meta_myc_results %>% 
                            filter(trait == "nmass" &
                                     nut_add == "n"),
                          aes(x = myc, y = (exp(estimate) - 1) * 100)) +
  geom_errorbar(aes(ymin = (exp(lowerCL) - 1) * 100, ymax = (exp(upperCL) - 1) * 100),
                width = 0.3, linewidth = 1) +
  geom_point(aes(fill = myc), shape = 21, size = 5) +
  geom_bracket(xmin = "Mining", xmax = "Scavenging", y.position = 98, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(x = 1.5, y = 102, label = "NS", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("pink", "darkred")) +
  scale_x_discrete(labels = c("mining", "scavenging")) +
  scale_y_continuous(limits = c(-25, 100), breaks = seq(-25, 100, 25)) +
  labs(x = expression(bold("Mycorrhizal-NAS")),
       y = expression(bold("N"["mass"]*" response to N addition (%)"))) +
  guides(fill = "none") +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
nmass_myc_plot

#####################################################################
# Pmass -- climate and species identity moderators
#####################################################################

# Nmass climate moderator model
padd_pmass_clim <- rma.mv(logr, logr_var, method = "REML", 
                          random = ~ 1 | exp, mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_p_mass" & 
                                     !is.na(gs_mat)))
summary(padd_pmass_clim)

# Pmass - temperature plot
pmass_tg_plot <- mod_results(padd_pmass_clim, mod = "gs_mat",
                             group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = (exp(estimate) - 1) * 100)) +
  geom_point(data = subset(meta_results, myvar == "leaf_p_mass" & 
                             manip_type == "p" & !is.na(gs_mat)),
             aes(x = gs_mat, y = (exp(logr) - 1) * 100, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = (exp(upperCL) - 1) * 100, ymin = (exp(lowerCL) - 1) * 100),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-50, 400), breaks = seq(0, 400, 100)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("P"["mass"]*" response to P addition (%)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
pmass_tg_plot 

# Pmass - aridity plot
pmass_ai_plot <- mod_results(padd_pmass_clim, mod = "gs_ai",
                             group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = (exp(estimate) - 1) * 100)) +
  geom_point(data = subset(meta_results, myvar == "leaf_p_mass" & 
                             manip_type == "p" & !is.na(gs_ai)),
             aes(x = gs_ai, y = (exp(logr) - 1) * 100, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = (exp(upperCL) - 1) * 100, ymin = (exp(lowerCL) - 1) * 100),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, 
              color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 3.4), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-50, 400), breaks = seq(0, 400, 100)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("AI"["g"]*" (unitless)")),
       y = expression(bold("P"["mass"]*" response to P addition (%)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
pmass_ai_plot 

# Pmass - photosynthetic pathway
pmass_photo_plot <- ggplot(data = meta_photo_results %>% 
                             filter(trait == "pmass" &
                                      nut_add == "p"),
                           aes(x = photo, y = (exp(estimate) - 1) * 100)) +
  geom_errorbar(aes(ymin = (exp(lowerCL) - 1) * 100, ymax = (exp(upperCL) - 1) * 100),
                width = 0.3, linewidth = 1) +
  geom_point(aes(fill = photo), shape = 21, size = 5) +
  geom_bracket(xmin = "C3", xmax = "C4", y.position = 390, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(x = 1.5, y = 400, label = "***", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("lightblue", "darkblue")) +
  scale_y_continuous(limits = c(-50, 400), breaks = seq(0, 400, 100)) +
  labs(x = "Photosynthetic pathway",
       y = expression(bold("P"["mass"]*" response to P addition (%)"))) +
  guides(fill = "none") +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
pmass_photo_plot

# Pmass - N fixation
pmass_nfix_plot <- ggplot(data = meta_nfix_results %>% 
                            filter(trait == "pmass" &
                                     nut_add == "p"),
                          aes(x = nfix, y = (exp(estimate) - 1) * 100)) +
  geom_errorbar(aes(ymin = (exp(lowerCL) - 1) * 100,
                    ymax = (exp(upperCL) - 1) * 100),
                width = 0.3, linewidth = 1) +
  geom_point(aes(fill = nfix), shape = 21, size = 5) +
  geom_bracket(xmin = "No", xmax = "Yes", y.position = 390, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(x = 1.5, y = 400, label = "***", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("lightblue", "darkblue")) +
  scale_x_discrete(labels = c("non-fixer", expression("N"["2"]*"-fixer"))) +
  scale_y_continuous(limits = c(-50, 400), breaks = seq(0, 400, 100)) +
  labs(x = expression(bold("N"["2"]*" fixation ability")),
       y = expression(bold("P"["mass"]*" response to P addition (%)"))) +
  guides(fill = "none") +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
pmass_nfix_plot

# Pmass - mycorrhizal acquisition strategy
pmass_myc_plot <- ggplot(data = meta_myc_results %>% 
                           filter(trait == "pmass" &
                                    nut_add == "p"),
                         aes(x = myc, y = (exp(estimate) - 1) * 100)) +
  geom_errorbar(aes(ymin = (exp(lowerCL) - 1) * 100, ymax = (exp(upperCL) - 1) * 100),
                width = 0.3, linewidth = 1) +
  geom_point(aes(fill = myc), shape = 21, size = 5) +
  geom_bracket(xmin = "Mining", xmax = "Scavenging", y.position = 390, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(x = 1.5, y = 400, label = "***", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(label = c("mining", "scavenging")) +
  scale_fill_manual(values = c("lightblue", "darkblue")) +
  scale_y_continuous(limits = c(-50, 400), breaks = seq(0, 400, 100)) +
  labs(x = expression(bold("Mycorrhizal-NAS")),
       y = expression(bold("P"["mass"]*" response to P addition (%)"))) +
  guides(fill = "none") +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
pmass_myc_plot

#####################################################################
# Write plots
#####################################################################

# Figure 2 - individual effects
png("../plots/CNP_fig2_indEffects.png", height = 16, width = 16, 
    units = "in", res = 600)
ggarrange(nadd_chemistry_plot, nadd_photo_plot, nadd_bio_plot,
          padd_chemistry_plot, padd_photo_plot, padd_bio_plot,
          npadd_chemistry_plot, npadd_photo_plot, npadd_bio_plot,
          nrow = 3, ncol = 3, labels = c("(a)", "(b)", "(c)",
                                         "(d)", "(e)", "(f)",
                                         "(g)", "(h)", "(i)"), 
          font.label = list(size = 18),
          align = "hv")
dev.off()

# Figure 3 - interaction effects
png("../plots/CNP_fig3_intEffects.png", height = 12, width = 5.5, 
    units = "in", res = 600)
ggarrange(npint_chemistry_plot, npint_photo_plot, npint_bio_plot, 
          nrow = 3, ncol = 1, labels = c("(a)", "(b)", "(c)"), 
          font.label = list(size = 18), align = "hv")
dev.off()

# Figure 4 - moderator plots
png("../plots/CNP_fig4_moderators.png", height = 12.5, width = 24,
    units = "in", res = 600)
ggarrange(nmass_tg_plot, nmass_ai_plot, nmass_nfix_plot, nmass_myc_plot, nmass_photo_plot,
          pmass_tg_plot, pmass_ai_plot, pmass_nfix_plot, pmass_myc_plot, pmass_photo_plot,
          nrow = 2, ncol = 5, legend = "bottom", common.legend = TRUE,
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)",
                     "(g)", "(h)", "(i)", "(j)"),
          font.label = list(size = 20), align = "hv")
dev.off()


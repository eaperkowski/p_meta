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
meta_ci <- read.csv("../data/CNPmeta_ci.csv")
meta_ci_int <- read.csv("../data/CNPmeta_ci_int.csv")

# Check file structure
head(meta_results)
head(meta_results_int)
head(meta_ci)
head(meta_ci_int)

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
                                                  "leaf_p_area")),
                              aes(x = var, y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), size = 1, width = 0.25) +
  geom_point(size = 4, fill = "red", shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c(expression("P"["area"]), 
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
                                                  "leaf_p_area")),
                              aes(x = var, y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), size = 1, width = 0.25) +
  geom_point(size = 4, fill = "blue", shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c(expression("P"["area"]), 
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
                                                   "leaf_p_area")),
                               aes(x = var, y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), 
                size = 1, width = 0.25) +
  geom_point(size = 4, fill = "magenta", shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c(expression("P"["area"]), 
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
        axis.text = element_text(size = 18, color = "black"))
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
                                                   "leaf_p_area")),
                               aes(x = var, y = middle_perc)) +
  geom_errorbar(aes(ymin = lower_perc, ymax = upper_perc), 
                size = 1, width = 0.25) +
  geom_point(size = 4, fill = "black", shape = 21) +
  geom_text(aes(label = k_sig), y = 55, fontface = "bold", size = 5) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c(expression("P"["area"]), 
                              expression("P"["mass"]), 
                              expression("N"["area"]),
                              expression("N"["mass"]),
                              expression("M"["area"]))) +
  scale_y_continuous(limits = c(-60, 60), breaks = seq(-60, 60, 30)) +
  labs(x = "", 
       y = NULL) +
  coord_flip() +
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
                           aes(x = var, y = middle_perc)) +
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
  labs(x = "", y = NULL) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(color = "black", size = 18))
npint_photo_plot

npint_bio_plot <- ggplot(data = meta_ci_int %>% 
                           drop_na(var) %>%
                           filter(var %in% c("rootshoot", 
                                             "rmf", 
                                             "bgb", 
                                             "agb", 
                                             "total_biomass")),
                         aes(x = var, y = middle_perc)) +
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
# Write plots
#####################################################################

# Figure 2 - individual effects
#png("../plots/CNP_fig2_indEffects.png", height = 12, width = 16, 
#    units = "in", res = 600)
ggarrange(nadd_chemistry_plot, nadd_photo_plot, nadd_bio_plot,
          padd_chemistry_plot, padd_photo_plot, padd_bio_plot,
          npadd_chemistry_plot, npadd_photo_plot, npadd_bio_plot,
          nrow = 3, ncol = 3, labels = c("(a)", "(b)", "(c)",
                                         "(d)", "(e)", "(f)",
                                         "(g)", "(h)", "(i)"), 
          font.label = list(size = 18),
          align = "hv")
#dev.off()

# Figure 3 - interaction effects
#png("../plots/CNP_fig3_intEffects.png", height = 12, width = 6, 
#    units = "in", res = 600)
ggarrange(npint_chemistry_plot, npint_photo_plot, npint_bio_plot, 
          nrow = 3, ncol = 1, labels = c("(a)", "(b)", "(c)"), 
          font.label = list(size = 18), align = "hv")
#dev.off()


ci_summary <- df_box_all %>%
  mutate(middle_perc = round(middle_perc, 1),
         upper_perc = round(upper_perc, 1),
         lower_perc = round(lower_perc, 1)) %>%
  dplyr::select(var, manip_type, middle_perc, upper_perc, lower_perc) %>%
  arrange(var)





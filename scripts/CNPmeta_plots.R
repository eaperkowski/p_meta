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
                    "lma", "leaf_ppue", "leaf_pnue", "jmax_vcmax", "jmax", "vcmax", "gsw", "rd", "asat",
                    "rootshoot", "rmf", "bgb", "bnpp", "agb", "anpp_p", "anpp_n", "anpp", "total_biomass",
                    "tbio_gm2", "tla")) %>%
  mutate(var = factor(var, levels = c("rootshoot", "rmf", "bgb", "bnpp", "agb", "anpp_p", "anpp_n",
                                      "anpp", "total_biomass", "tbio_gm2", "leaf_ppue", "leaf_pnue", 
                                      "jmax_vcmax", "jmax", "vcmax", "gsw", "rd", "asat",
                                      "leaf_np", "leaf_p_area", "leaf_p_mass", 
                                      "leaf_n_area", "leaf_n_mass", "lma", "tla")),
         nut_add = factor(nut_add, levels = c("np", "p", "n")),
         ci_range_plot = str_c("[", sprintf("%.3f", ci.lb), ", ", sprintf("%.3f", ci.ub), "]"),
         k_plot = str_c("(", k, ")"))

meta_ci_int <- read.csv("../data/CNPmeta_ci_int.csv") %>%
  filter(var %in% c("leaf_np", "leaf_p_area", "leaf_p_mass", "leaf_n_area", "leaf_n_mass",
                    "lma", "leaf_ppue", "leaf_pnue", "jmax_vcmax", "jmax", "vcmax", "gsw", "rd", "asat",
                    "rootshoot", "rmf", "bgb", "bnpp", "agb", "anpp", "total_biomass", "tla", "tbio_gm2", "tla")) %>%
  mutate(var = factor(var, levels = c("leaf_np", "leaf_p_area", "leaf_p_mass", 
                                      "leaf_n_area", "leaf_n_mass",
                                      "lma", "leaf_ppue", "leaf_pnue", 
                                      "jmax_vcmax", "jmax", "vcmax", "gsw", "rd", "asat", "rootshoot",
                                      "rmf", "bgb", "bnpp", "agb", "anpp", "total_biomass", "tbio_gm2", "tla")),
         ci_range_plot = str_c("[", sprintf("%.3f", ci.lb), ", ", sprintf("%.3f", ci.ub), "]"),
         k_plot = str_c("(", k, ")"),
         int_type = ifelse(var == "anpp" | var == "leaf_np" | var == "jmax",
                           "synergistic", 
                           ifelse(var == "leaf_p_area",
                                  "antagonistic", "additive")))

# Load species moderator results for individual and interaction effects
ind_pft_results <- read.csv("../data/CNPmeta_pft_moderators.csv")
int_pft_results <- read.csv("../data/CNPmeta_clim_moderators_int.csv")

# Check file structure
head(meta_results)
head(meta_results_int)
head(meta_ci)
head(meta_ci_int)

#####################################################################
# Leaf chemistry individual plot
#####################################################################
chemistry_ind_plot <- ggplot(data = meta_ci %>% 
                               drop_na(var) %>% 
                               filter(var %in% c("lma", 
                                                 "leaf_n_mass", 
                                                 "leaf_n_area", 
                                                 "leaf_p_mass", 
                                                 "leaf_p_area",
                                                 "leaf_np")),
                             aes(x = var, y = estimate, 
                                 group = factor(nut_add, 
                                                levels = c("np", "p", "n")))) +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                position = position_dodge(0.75), size = 1, width = 0.5) +
  geom_point(aes(fill = nut_add),
             position = position_dodge(0.75), size = 4, shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_text(data = subset(meta_ci_int, var %in% c("lma", 
                                                  "leaf_n_mass", 
                                                  "leaf_n_area", 
                                                  "leaf_p_mass", 
                                                  "leaf_p_area",
                                                  "leaf_np")),
            aes(label = k_plot, y = 1.1), size = 6) +
  scale_fill_manual(values = c("red", "blue", "magenta"),
                    breaks = c("n", "p", "np"),
                    labels = c("N", "P", "N+P")) +
  scale_x_discrete(labels = c("Leaf N:P",
                              expression(italic("P")["area"]), 
                              expression(italic("P")["mass"]), 
                              expression(italic("N")["area"]),
                              expression(italic("N")["mass"]),
                              expression(italic("M")["area"]))) +
  scale_y_continuous(limits = c(-0.6, 1.2), breaks = seq(-0.6, 1.2, 0.6)) +
  labs(title = "Leaf responses to N, P, N+P addition",
       x = "", 
       y = "",
       fill = "Nutrient addition") +
  coord_flip() +
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
chemistry_ind_plot

#####################################################################
# Leaf chemistry interaction plot
#####################################################################
chemistry_int_plot <- ggplot(data = meta_ci_int %>% 
                               drop_na(var) %>% 
                               filter(var %in% c("lma", 
                                                 "leaf_n_mass", 
                                                 "leaf_n_area", 
                                                 "leaf_p_mass", 
                                                 "leaf_p_area",
                                                 "leaf_np")),
                             aes(x = var, y = estimate)) +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 1, width = 0.25) +
  geom_point(aes(fill = int_type), size = 4, shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_text(aes(label = k_plot, y = 0.95), size = 6) +
  scale_fill_manual(values = c("black", "red", "blue")) +
  scale_x_discrete(labels = c("Leaf N:P",
                              expression(italic("P")["area"]), 
                              expression(italic("P")["mass"]), 
                              expression(italic("N")["area"]),
                              expression(italic("N")["mass"]),
                              expression(italic("M")["area"]))) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(title = "Leaf interaction responses",
       x = "", 
       y = "",
       fill = "Interaction type") +
  coord_flip() +
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
chemistry_int_plot

#####################################################################
# Photosynthesis individual plot
#####################################################################
photo_ind_plot <- ggplot(data = meta_ci %>% 
                                    drop_na(var) %>% 
                                    filter(var %in% c("asat",
                                                      "rd",
                                                      "gsw",
                                                      "vcmax", 
                                                      "jmax",
                                                      "leaf_pnue", 
                                                      "leaf_ppue")),
                                  aes(x = var, y = estimate, 
                                      group = factor(nut_add, 
                                                     levels = c("np", "p", "n")))) +
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                position = position_dodge(0.75), size = 1, width = 0.5) +
  geom_point(aes(fill = nut_add),
             position = position_dodge(0.75), size = 4, shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_text(data = subset(meta_ci_int, var %in% c("asat",
                                                  "rd",
                                                  "gsw",
                                                  "vcmax", 
                                                  "jmax",
                                                  "leaf_pnue", 
                                                  "leaf_ppue")),
            aes(label = k_plot, y = 1.1), size = 6) +
  #geom_text(aes(label = ci_range, y = 1.2, fontface = bold_yn), size = 5,
  #          position = position_dodge(0.75), hjust = 1) +
  scale_fill_manual(values = c("red", "blue", "magenta"),
                    breaks = c("n", "p", "np"),
                    labels = c("N", "P", "N+P")) +
  scale_x_discrete(labels = c("PPUE",
                              "PNUE",
                              expression(italic("J")["max"]),
                              expression(italic("V")["cmax"]),
                              expression(italic("g")["sw"]),
                              expression(italic("R")["d"]),
                              expression(italic("A")["sat"]))) +
  scale_y_continuous(limits = c(-0.6, 1.2), breaks = seq(-0.6, 1.2, 0.6)) +
  labs(title = "Photo. responses to N, P, N+P addition",
       x = NULL, 
       y = "",
       fill = "Nutrient addition") +
  coord_flip() +
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
photo_ind_plot

#####################################################################
# Photosynthesis interaction plot
#####################################################################
photo_int_plot <- ggplot(data = meta_ci_int %>% 
                               drop_na(var) %>% 
                               filter(var %in% c("asat",
                                                 "rd",
                                                 "gsw",
                                                 "vcmax", 
                                                 "jmax",
                                                 "leaf_pnue", 
                                                 "leaf_ppue")),
                             aes(x = var, y = estimate)) +
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 1, width = 0.25) +
  geom_point(aes(fill = int_type), size = 4, shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_text(aes(label = k_plot, y = 0.95), size = 6) +
  scale_fill_manual(values = c("black", "blue")) +
  scale_x_discrete(labels = c("PPUE",
                              "PNUE",
                              expression(italic("J")["max"]),
                              expression(italic("V")["cmax"]),
                              expression(italic("g")["sw"]),
                              expression(italic("R")["d"]),
                              expression(italic("A")["sat"]))) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(title = "Photo. interaction responses",
       x = "", 
       y = "",
       fill = "Interaction type") +
  coord_flip() +
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
photo_int_plot

#####################################################################
# Biomass individual plot
#####################################################################
bio_ind_plot <- ggplot(data = meta_ci %>% 
                         drop_na(var) %>% 
                         filter(var %in% c("rootshoot", 
                                           "rmf", 
                                           "bnpp",
                                           "anpp",
                                           "tbio_gm2",
                                           "tla")),
                       aes(x = var, y = estimate, 
                           group = factor(nut_add, 
                                          levels = c("np", "p", "n")))) +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                position = position_dodge(0.75), size = 1, width = 0.5) +
  geom_point(aes(fill = nut_add),
             position = position_dodge(0.75), size = 4, shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_text(data = subset(meta_ci_int, var %in% c("rootshoot", 
                                                  "rmf", 
                                                  "bnpp",
                                                  "anpp",
                                                  "tbio_gm2",
                                                  "tla")),
            aes(label = k_plot, y = 1.1), size = 6) +
  #geom_text(aes(label = ci_range, y = 1.5, fontface = bold_yn), size = 5,
  #          position = position_dodge(0.75), hjust = 1) +
  scale_fill_manual(values = c("red", "blue", "magenta"),
                    breaks = c("n", "p", "np"),
                    labels = c("N", "P", "N+P")) +
  scale_x_discrete(labels = c("Root:shoot",
                              "Root mass fraction",
                              "Belowground prod.",
                              "Aboveground prod.",
                              "Biomass prod.",
                              "Leaf prod.")) +
  scale_y_continuous(limits = c(-0.6, 1.2), breaks = seq(-0.6, 1.2, 0.6)) +
  labs(title = "Biomass responses to N, P, N+P addition",
       x = NULL, 
       y = "Log-response ratio",
       fill = "Nutrient addition ") +
  coord_flip() +
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
bio_ind_plot

#####################################################################
# Biomass interaction plot
#####################################################################
bio_int_plot <- ggplot(data = meta_ci_int %>% 
                           drop_na(var) %>% 
                           filter(var %in% c("rootshoot", 
                                             "rmf", 
                                             "bnpp",
                                             "anpp",
                                             "tbio_gm2",
                                             "tla")),
                         aes(x = var, y = estimate)) +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 1, width = 0.25) +
  geom_point(aes(fill = int_type), size = 4, shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_text(aes(label = k_plot, y = 0.95), size = 6) +
  scale_fill_manual(values = c("black", "blue")) +
  scale_x_discrete(labels = c("Root:shoot",
                              "Root mass fraction",
                              "Belowground prod.",
                              "Aboveground prod.",
                              "Biomass prod.",
                              "Leaf prod.")) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(title = "Biomass interaction responses",
       x = "", 
       y = expression(bold("Interaction effect size (")*bar(bolditalic("d")[bold("NP")])*bold(")")),
       fill = "Interaction type") +
  coord_flip() +
  theme_classic(base_size = 20) +
  theme(title = element_text(face = "bold", hjust = 0),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22, hjust = 0.5),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 18, color = "black"),
        panel.grid = element_blank())
bio_int_plot


#####################################################################
# Figure 2: Individual effects 
#####################################################################

png("../plots/CNP_fig2_ind_responses.png", height = 16, width = 10,
    units = "in", res = 600)
ggarrange(chemistry_ind_plot, photo_ind_plot, bio_ind_plot,
          ncol = 1, nrow = 3, align = "hv",
          common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)", "(c)"),
          font.label = list(size = 20, face = "bold"))
dev.off()


#####################################################################
# Figure 3: Interaction effects 
#####################################################################

png("../plots/CNP_fig3_int_responses.png", height = 16, width = 8,
    units = "in", res = 600)
ggarrange(chemistry_int_plot, photo_int_plot, bio_int_plot,
          ncol = 1, nrow = 3, align = "hv",
          common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)", "(c)"),
          font.label = list(size = 20, face = "bold"))
dev.off()

#####################################################################
# Narea -- climate moderators
#####################################################################
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "leaf_n_area" & 
                                        !is.na(gs_mat))) +
  geom_point(aes(x = gs_mat, y = logr))


# Narea climate moderator model
nadd_narea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "n" & 
                                     myvar == "leaf_n_area" & 
                                     !is.na(gs_mat) & gs_ai < 3 & logr > -1 & logr < 0.95))

summary(nadd_narea_clim)

# Narea - temperature plot
narea_tg_plot <- mod_results(nadd_narea_clim, mod = "gs_mat",
                             group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_area" & 
                             nut_add == "n" & !is.na(gs_mat) & gs_ai < 3 & 
                             logr > -1 & logr < 0.95),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("area")]*bold(" N response to ")*bolditalic("T")[bold("g")]),
       x = "",
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
narea_tg_plot 

# Narea - aridity plot
narea_ai_plot <- mod_results(nadd_narea_clim, mod = "gs_ai",
                             group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_area" & 
                             nut_add == "n" & !is.na(gs_mat) & gs_ai < 3 & 
                             logr > -1 & logr < 0.95),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("area")]*bold(" N response to ")*bolditalic("MI")[bold("g")]),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
narea_ai_plot 

# Narea - PAR plot
narea_par_plot <- ggplot(data = subset(meta_results, myvar == "leaf_n_area" & 
                                         nut_add == "n" & !is.na(gs_par) & gs_ai < 3 & 
                                         logr > -1 & logr < 0.95),
                         aes(x = gs_par, y = logr, size = 1/logr_se)) +
  geom_point(alpha = 0.30, fill = "red", shape = 21) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("area")]*bold(" N response to ")*bolditalic("PAR")[bold("g")]),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
narea_par_plot 

#####################################################################
# Parea -- climate moderators
#####################################################################
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "leaf_p_area" & 
                                        !is.na(gs_mat)& logr >-1 & gs_ai < 3)) +
  geom_point(aes(x = gs_mat, y = logr))

# Narea climate moderator model
padd_parea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "p" & 
                                     myvar == "leaf_p_area" & 
                                     !is.na(gs_mat) & gs_ai < 3 & logr > -1))
summary(padd_parea_clim)

# Parea - temperature plot
parea_tg_plot <- mod_results(padd_parea_clim, mod = "gs_mat",
                             group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             nut_add == "p" & !is.na(gs_mat)  & logr >-1 & gs_ai < 3),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, color = "blue") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("area")]*bold(" P response to ")*bolditalic("T")[bold("g")]),
       x = "",
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
parea_tg_plot 

# Parea - aridity plot
parea_ai_plot <- mod_results(padd_parea_clim, mod = "gs_ai",
                             group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             nut_add == "p" & !is.na(gs_mat)  & logr >-1 & gs_ai < 3),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, color = "blue") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("area")]*bold(" P response to ")*bolditalic("MI")[bold("g")]),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
parea_ai_plot 

# Parea - PAR plot
parea_par_plot <- mod_results(padd_parea_clim, mod = "gs_par",
                              group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             nut_add == "p" & !is.na(gs_mat)  & logr >-1 & gs_ai < 3),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, color = "blue") +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("area")]*bold(" P response to ")*bolditalic("PAR")[bold("g")]),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
parea_par_plot 

#####################################################################
# Leaf N:P interaction -- climate moderators
#####################################################################
ggplot(data = meta_results_int %>% filter(response == "leaf_np" & 
                                        !is.na(gs_mat) & gs_ai < 3)) +
  geom_point(aes(x = gs_mat, y = dNPi))

# Model
int_leafnp_clim <- rma.mv(yi = dNPi,
                          V = vNPi,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results_int %>% 
                            filter(response == "leaf_np" & 
                                     !is.na(gs_mat) & gs_ai < 3 & dNPi > -2))

# Leaf N:P interaction - temperature plot
leafnp_int_tg_plot <- mod_results(int_leafnp_clim, mod = "gs_mat",
                             group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results_int, response == "leaf_np" & 
                             !is.na(gs_mat) & dNPi > -2),
             aes(x = gs_mat, y = dNPi, size = 1/dNPi_se), 
             alpha = 0.30, fill = "black", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "black") +
  geom_smooth(method = "loess", linewidth = 2, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-1.4, 1.4), breaks = seq(-1.4, 1.4, 0.7)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bold("Leaf N:P int. response to ")*bolditalic("T")[bold("g")]),
       x = expression(bolditalic("T")[bold("g")]*bold(" ("*degree*"C)")),
       y = expression(bold("Hedge's ")*bolditalic("d")*" ("*bar(bolditalic("d")[bold("NP")])*bold(")")),
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
leafnp_int_tg_plot 

# Leaf N:P interaction - aridity plot
leafnp_int_ai_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results_int, response == "leaf_np" & 
                             !is.na(gs_mat)  & dNPi >-2 & gs_ai < 3),
             aes(x = gs_ai, y = dNPi, size = 1/dNPi_se), 
             alpha = 0.30, fill = "black", shape = 21) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-1.4, 1.4), breaks = seq(-1.4, 1.4, 0.7)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bold("Leaf N:P int. response to ")*bolditalic("MI")[bold("g")]),       x = expression(bolditalic("MI")[bold("g")]*bold(" (unitless)")),
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
leafnp_int_ai_plot 

# Leaf N:P - PAR plot
leafnp_int_par_plot <- mod_results(int_leafnp_clim, mod = "gs_par",
                                   group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results_int, response == "leaf_np" & 
                             !is.na(gs_mat)  & dNPi >-2 & gs_ai < 3),
             aes(x = gs_par, y = dNPi, size = 1/dNPi_se), 
             alpha = 0.30, fill = "black", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "black") +
  geom_smooth(method = "loess", linewidth = 2, color = "black") +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-1.4, 1.4), breaks = seq(-1.4, 1.4, 0.7)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bold("Leaf N:P int. response to ")*bolditalic("PAR")[bold("g")]),       x = expression(bolditalic("MI")[bold("g")]*bold(" (unitless)")),
       x = expression(bolditalic("PAR")[bold("g")]*bold(" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
leafnp_int_par_plot 


#####################################################################
# Figure 4: climate responses
#####################################################################

png("../plots/CNP_fig4_climate_responses.png", height = 14, width = 14,
    units = "in", res = 600)
ggarrange(narea_tg_plot, narea_ai_plot, narea_par_plot,
          parea_tg_plot, parea_ai_plot, parea_par_plot,
          leafnp_int_tg_plot, leafnp_int_ai_plot, leafnp_int_par_plot,
          ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)", "(c)", 
                     "(d)", "(e)", "(f)", 
                     "(g)", "(h)", "(i)"),
          font.label = list(size = 22))
dev.off()


#####################################################################
# N addition: Photosynthetic pathway
#####################################################################
# Plot
nadd_photo_plot <- ggplot(data = ind_pft_results %>% 
                            filter(trait %in% c("nmass", "narea",
                                                "asat", "rd",
                                                "vcmax", "jmax") & nut_add == "n" & mod == "photo"),
                          aes(x = factor(trait, 
                                         levels = c("jmax", "vcmax", "rd", 
                                                    "asat", 
                                                    "narea", "nmass")), 
                              y = estimate, 
                              shape = factor(comp, levels = c("C3", "C4")),
                              fill = factor(comp, levels = c("C3", "C4")))) +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL),
                position = position_dodge(width = 0.75), width = 0.2, color = "black", size = 1) +
  geom_point(position = position_dodge(width = 0.75), size = 5) +
  geom_text(aes(x = 1, y = 0.95, label = " NS "), size = 7, fontface = "bold") +
  geom_text(aes(x = 2, y = 0.95, label = " NS "), size = 7, fontface = "bold") +
  geom_text(aes(x = 3, y = 0.95, label = " ** "), size = 7, fontface = "bold") +
  geom_text(aes(x = 4, y = 0.95, label = " NS "), size = 7, fontface = "bold") +
  geom_text(aes(x = 5, y = 0.95, label = " * "), size = 7, fontface = "bold") +
  geom_text(aes(x = 6, y = 0.95, label = " * "), size = 7, fontface = "bold") +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"["3"]),
                                expression("C"["4"]))) +
  scale_fill_manual(values = c("ivory", "red2"),
                    labels = c(expression("C"["3"]),
                               expression("C"["4"]))) +
  scale_x_discrete(labels = c(expression(italic("J")["max"]),
                              expression(italic("V")["cmax"]),
                              expression(italic("R")["d"]),
                              expression(italic("A")["sat"]),
                              expression(italic("N")["area"]),
                              expression(italic("N")["mass"]))) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  coord_flip() +
  labs(title = "N addition responses to photo. pathway",
       x = "",
       y = "",
       shape = expression(bold("Photo.\npathway")),
       fill = expression(bold("Photo.\npathway"))) +   
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
nadd_photo_plot

#####################################################################
# N addition: N2 fixation
#####################################################################
# Plot
nadd_nfix_plot <- ggplot(data = ind_pft_results %>% 
                            filter(trait %in% c("nmass", "narea",
                                                "asat", "rd",
                                                "vcmax", "jmax") & nut_add == "n" & mod == "nfix"),
                          aes(x = factor(trait, 
                                         levels = c("jmax", "vcmax", "rd", "asat", 
                                                    "narea", "nmass")), 
                              y = estimate, 
                              shape = factor(comp, levels = c("Yes", "No")),
                              fill = factor(comp, levels = c("Yes", "No")))) +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_text(aes(x = 1, y = 0.95, label = " * "), size = 7, fontface = "bold") +
  geom_text(aes(x = 2, y = 0.95, label = " ** "), size = 7, fontface = "bold") +
  geom_text(aes(x = 3, y = 0.95, label = " * "), size = 7, fontface = "bold") +
  geom_text(aes(x = 4, y = 0.95, label = "***"), size = 7, fontface = "bold") +
  geom_text(aes(x = 5, y = 0.95, label = " NS "), size = 7, fontface = "bold") +
  geom_text(aes(x = 6, y = 0.95, label = "***"), size = 7, fontface = "bold") +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL),
                position = position_dodge(width = 0.75),
                width = 0.2, color = "black", size = 1) +
  geom_point(position = position_dodge(width = 0.75), size = 5) +
  scale_shape_manual(values = c(23, 21),
                     labels = c(expression("N"["2"]*"-fixer"),
                                "non-fixer")) +
  scale_fill_manual(values = c("red2", "ivory"),
                     labels = c(expression("N"["2"]*"-fixer"),
                                "non-fixer")) +
  scale_x_discrete(labels = c(expression(italic("J")["max"]),
                              expression(italic("V")["cmax"]),
                              expression(italic("R")["d"]),
                              expression(italic("A")["sat"]),
                              expression(italic("N")["area"]),
                              expression(italic("N")["mass"]))) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  coord_flip() +
  labs(title = "N addition responses to N-fixation ability",
       x = "",
       y = "",
       shape = expression(bold("N-fixation\nability")),
       fill = expression(bold("N-fixation\nability"))) +   
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
nadd_nfix_plot

#####################################################################
# N addition: Mycorrhizal acquisition strategy
#####################################################################
# Plot
nadd_myc_plot <- ggplot(data = ind_pft_results %>% 
                          filter(trait %in% c("nmass", "narea",
                                              "asat", "rd",
                                              "vcmax", "jmax") & nut_add == "n" & mod == "myc_nas"),
                        aes(x = factor(trait, levels = c("jmax", "vcmax", "rd", "asat", 
                                                         "narea", "nmass")), 
                            y = estimate, 
                            shape = factor(comp, levels = c("Scavenging", "Mining")),
                            fill = factor(comp, levels = c("Scavenging", "Mining")))) +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL),
                position = position_dodge(width = 0.75), 
                width = 0.2, color = "black", size = 1) +
  geom_point(position = position_dodge(width = 0.75), size = 5) +
  geom_text(aes(x = 1, y = 0.95, label = "***"), size = 7, fontface = "bold") +
  geom_text(aes(x = 2, y = 0.95, label = "***"), size = 7, fontface = "bold") +
  geom_text(aes(x = 3, y = 0.95, label = "***"), size = 7, fontface = "bold") +
  geom_text(aes(x = 4, y = 0.95, label = "***"), size = 7, fontface = "bold") +
  geom_text(aes(x = 5, y = 0.95, label = " ** "), size = 7, fontface = "bold") +
  geom_text(aes(x = 6, y = 0.95, label = " NS "), size = 7, fontface = "bold") +
  scale_shape_manual(values = c(21, 24),
                     labels = c("scavenging", "mining")) +
  scale_fill_manual(values = c("ivory", "red2"),
                    labels = c("scavenging",
                               "mining")) +
  scale_x_discrete(labels = c(expression(italic("J")["max"]),
                              expression(italic("V")["cmax"]),
                              expression(italic("R")["d"]),
                              expression(italic("A")["sat"]),
                              expression(italic("N")["area"]),
                              expression(italic("N")["mass"]))) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  coord_flip() +
  labs(title = "N addition responses to mycorrhizal strategy",
       x = "",
       y = "",
       shape = expression(bold("Mycorrhizal\nacquisition\nstrategy")),
       fill = expression(bold("Mycorrhizal\nacquisition\nstrategy"))) +   
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
nadd_myc_plot

#####################################################################
# Figure 5: N addition PFT moderators
#####################################################################
png("../plots/CNP_fig5_Nadd_pft_responses.png", height = 16, width = 8.5,
    units = "in", res = 600)
ggarrange(nadd_nfix_plot, nadd_myc_plot, nadd_photo_plot,
          nrow = 3, labels = c("(a)", "(b)", "(c)"), align = "hv",
          font.label = list(size = 20))
dev.off()



#####################################################################
# N addition: Photosynthetic pathway
#####################################################################
# Plot
padd_photo_plot <- ggplot(data = ind_pft_results %>% 
                            filter(trait %in% c("pmass", "parea",
                                                "asat", "rd",
                                                "vcmax", "jmax") & nut_add == "p" & mod == "photo"),
                          aes(x = factor(trait, 
                                         levels = c("jmax", "vcmax", "rd", 
                                                    "asat", "parea", "pmass")), 
                              y = estimate, 
                              shape = factor(comp, levels = c("C3", "C4")),
                              fill = factor(comp, levels = c("C3", "C4")))) +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL),
                position = position_dodge(width = 0.75), width = 0.2, color = "black", size = 1) +
  geom_point(position = position_dodge(width = 0.75), size = 5) +
  geom_text(aes(x = 1, y = 1.45, label = "NS"), size = 7, fontface = "bold") +
  geom_text(aes(x = 2, y = 1.45, label = "NS"), size = 7, fontface = "bold") +
  geom_text(aes(x = 3, y = 1.45, label = "NS"), size = 7, fontface = "bold") +
  geom_text(aes(x = 4, y = 1.45, label = "NS"), size = 7, fontface = "bold") +
  geom_text(aes(x = 5, y = 1.45, label = "NS"), size = 7, fontface = "bold") +
  geom_text(aes(x = 6, y = 1.45, label = "NS"), size = 7, fontface = "bold") +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"["3"]),
                                expression("C"["4"]))) +
  scale_fill_manual(values = c("ivory", "blue2"),
                    labels = c(expression("C"["3"]),
                               expression("C"["4"]))) +
  scale_x_discrete(labels = c(expression(italic("J")["max"]),
                              expression(italic("V")["cmax"]),
                              expression(italic("R")["d"]),
                              expression(italic("A")["sat"]),
                              expression(italic("P")["area"]),
                              expression(italic("P")["mass"]))) +
  scale_y_continuous(limits = c(-0.5, 1.5), breaks = seq(-0.5, 1.5, 0.5)) +
  coord_flip() +
  labs(title = "P addition responses to photo. pathway",
       x = "",
       y = "",
       shape = expression(bold("Photo.\npathway")),
       fill = expression(bold("Photo.\npathway"))) +   
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
padd_photo_plot

#####################################################################
# P addition: N2 fixation
#####################################################################
# Plot
padd_nfix_plot <- ggplot(data = ind_pft_results %>% 
                           filter(trait %in% c("pmass", "parea",
                                               "asat", "rd",
                                               "vcmax", "jmax") & 
                                    nut_add == "p" & mod == "nfix"),
                         aes(x = factor(trait, 
                                        levels = c("jmax", "vcmax", "rd", "asat", 
                                                   "parea", "pmass")), 
                             y = estimate, 
                             shape = factor(comp, levels = c("Yes", "No")),
                             fill = factor(comp, levels = c("Yes", "No")))) +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_text(aes(x = 1, y = 1.45, label = " NS "), size = 7, fontface = "bold") +
  geom_text(aes(x = 2, y = 1.45, label = " * "), size = 7, fontface = "bold") +
  geom_text(aes(x = 3, y = 1.45, label = " NS "), size = 7, fontface = "bold") +
  geom_text(aes(x = 4, y = 1.45, label = "***"), size = 7, fontface = "bold") +
  geom_text(aes(x = 5, y = 1.45, label = " NS "), size = 7, fontface = "bold") +
  geom_text(aes(x = 6, y = 1.45, label = "***"), size = 7, fontface = "bold") +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL),
                position = position_dodge(width = 0.75),
                width = 0.2, color = "black", size = 1) +
  geom_point(position = position_dodge(width = 0.75), size = 5) +
  scale_shape_manual(values = c(23, 21),
                     labels = c(expression("N"["2"]*"-fixer"),
                                "non-fixer")) +
  scale_fill_manual(values = c("blue2", "ivory"),
                    labels = c(expression("N"["2"]*"-fixer"),
                               "non-fixer")) +
  scale_x_discrete(labels = c(expression(italic("J")["max"]),
                              expression(italic("V")["cmax"]),
                              expression(italic("R")["d"]),
                              expression(italic("A")["sat"]),
                              expression(italic("P")["area"]),
                              expression(italic("P")["mass"]))) +
  scale_y_continuous(limits = c(-0.5, 1.5), breaks = seq(-0.5, 1.5, 0.5)) +
  coord_flip() +
  labs(title = "P addition responses to N-fixation ability",
       x = "",
       y = "",
       shape = expression(bold("N-fixation\nability")),
       fill = expression(bold("N-fixation\nability"))) +   
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
padd_nfix_plot

#####################################################################
# N addition: Mycorrhizal acquisition strategy
#####################################################################
# Plot
padd_myc_plot <- ggplot(data = ind_pft_results %>% 
                          filter(trait %in% c("pmass", "parea",
                                              "asat", "rd",
                                              "vcmax", "jmax") & nut_add == "p" & mod == "myc_nas"),
                        aes(x = factor(trait, levels = c("jmax", "vcmax", "rd", "asat", 
                                                         "parea", "pmass")), 
                            y = estimate, 
                            shape = factor(comp, levels = c("Scavenging", "Mining")),
                            fill = factor(comp, levels = c("Scavenging", "Mining")))) +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL),
                position = position_dodge(width = 0.75), 
                width = 0.2, color = "black", size = 1) +
  geom_point(position = position_dodge(width = 0.75), size = 5) +
  geom_text(aes(x = 1, y = 1.45, label = " NS "), size = 7, fontface = "bold") +
  geom_text(aes(x = 2, y = 1.45, label = " NS "), size = 7, fontface = "bold") +
  geom_text(aes(x = 3, y = 1.45, label = " NS "), size = 7, fontface = "bold") +
  geom_text(aes(x = 4, y = 1.45, label = " NS "), size = 7, fontface = "bold") +
  geom_text(aes(x = 5, y = 1.45, label = "***"), size = 7, fontface = "bold") +
  geom_text(aes(x = 6, y = 1.45, label = "***"), size = 7, fontface = "bold") +
  scale_shape_manual(values = c(21, 24),
                     labels = c("scavenging", "mining")) +
  scale_fill_manual(values = c("ivory", "blue2"),
                    labels = c("scavenging",
                               "mining")) +
  scale_x_discrete(labels = c(expression(italic("J")["max"]),
                              expression(italic("V")["cmax"]),
                              expression(italic("R")["d"]),
                              expression(italic("A")["sat"]),
                              expression(italic("P")["area"]),
                              expression(italic("P")["mass"]))) +
  scale_y_continuous(limits = c(-0.5, 1.5), breaks = seq(-0.5, 1.5, 0.5)) +
  coord_flip() +
  labs(title = "P addition responses to mycorrhizal strategy",
       x = "",
       y = "",
       shape = expression(bold("Mycorrhizal\nacquisition\nstrategy")),
       fill = expression(bold("Mycorrhizal\nacquisition\nstrategy"))) +   
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
padd_myc_plot

#####################################################################
# Figure 5: P addition PFT moderators
#####################################################################
png("../plots/CNP_fig6_Padd_pft_responses.png", height = 16, width = 8.5,
    units = "in", res = 600)
ggarrange(padd_nfix_plot, padd_myc_plot, padd_photo_plot,
          nrow = 3, labels = c("(a)", "(b)", "(c)"), align = "hv",
          font.label = list(size = 20))
dev.off()






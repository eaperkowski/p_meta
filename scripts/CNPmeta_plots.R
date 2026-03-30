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
         int_type = ifelse(var == "leaf_p_mass" | var == "anpp" | var == "leaf_np",
                           "synergistic", "additive"))

# Load species moderator results for individual and interaction effects
ind_pft_results <- read.csv("../data/CNPmeta_pft_moderators.csv")
int_pft_results <- read.csv("../data/CNPmeta_clim_moderators_int.csv")

# Check file structure
head(meta_results)
head(meta_results_int)
head(meta_ci)
head(meta_ci_int)

#####################################################################
# Leaf hemistry individual plot
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
                position = position_dodge(0.75), size = 1, width = 0.6) +
  geom_point(aes(fill = nut_add),
             position = position_dodge(0.75), size = 4, shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_text(aes(label = k_plot, y = 1.1), size = 6) +
  scale_fill_manual(values = c("red", "blue", "magenta"),
                    breaks = c("n", "p", "np")) +
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
  theme_classic(base_size = 18) +
  theme(title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "plain"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 18, color = "black"),
        panel.grid = element_blank())
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
  scale_fill_manual(values = c("black", "red")) +
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
  theme_classic(base_size = 18) +
  theme(title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "plain"),
        #axis.text.y = element_text(size = 20, color = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 18, color = "black"),
        panel.grid = element_blank())

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
                position = position_dodge(0.75), size = 1, width = 0.6) +
  geom_point(aes(fill = nut_add),
             position = position_dodge(0.75), size = 4, shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_text(aes(label = k_plot, y = 1.1), size = 6) +
  #geom_text(aes(label = ci_range, y = 1.2, fontface = bold_yn), size = 5,
  #          position = position_dodge(0.75), hjust = 1) +
  scale_fill_manual(values = c("red", "blue", "magenta"),
                    breaks = c("n", "p", "np")) +
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
  theme_classic(base_size = 18) +
  theme(title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "plain"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 18, color = "black"),
        panel.grid = element_blank())
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
  scale_fill_manual(values = c("black", "red")) +
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
  theme_classic(base_size = 18) +
  theme(title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "plain"),
        #axis.text.y = element_text(size = 20, color = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 18, color = "black"),
        panel.grid = element_blank())
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
                position = position_dodge(0.75), size = 1, width = 0.6) +
  geom_point(aes(fill = nut_add),
             position = position_dodge(0.75), size = 4, shape = 21) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_text(aes(label = k_plot, y = 1.1), size = 6) +
  #geom_text(aes(label = ci_range, y = 1.5, fontface = bold_yn), size = 5,
  #          position = position_dodge(0.75), hjust = 1) +
  scale_fill_manual(values = c("red", "blue", "magenta"),
                    breaks = c("n", "p", "np")) +
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
  theme_classic(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "plain", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 18, color = "black"),
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
  scale_fill_manual(values = c("black", "red")) +
  scale_x_discrete(labels = c("Root:shoot",
                              "Root mass fraction",
                              "Belowground prod.",
                              "Aboveground prod.",
                              "Biomass prod.",
                              "Leaf prod.")) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(title = "Biomass interaction responses",
       x = "", 
       y = expression("Interaction effect size ("*bar(italic("d")[NP])*")"),
       fill = "Interaction type") +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(title = element_text(face = "bold", hjust = 0),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "plain", size = 22),
        #axis.text.y = element_text(size = 20, color = "black"),
        axis.text.y = element_blank(),
        panel.grid = element_blank())
bio_int_plot


#####################################################################
# Write plot that summarizes individual and interactive effects of
# N and P addition on leaf chemistry, photosynthesis, and biomass
# target variables
#####################################################################

png("../plots/fig2_full_responses.png", height = 15, width = 15,
    units = "in", res = 600)
ggarrange(ggarrange(chemistry_ind_plot, photo_ind_plot, bio_ind_plot,
                    ncol = 1, nrow = 3, align = "hv",
                    common.legend = TRUE, legend = "bottom",
                    labels = c("(a)", "(c)", "(e)"),
                    font.label = list(size = 20, face = "bold")),
          ggarrange(chemistry_int_plot, photo_int_plot, bio_int_plot,
                    ncol = 1, nrow = 3, align = "hv",
                    common.legend = TRUE, legend = "bottom",
                    labels = c("(b)", "(d)", "(f)"),
                    font.label = list(size = 20, face = "bold"),
                    hjust = 0.25),
          ncol = 2, nrow = 1, align = "hv", widths = c(1, 0.6))
dev.off()


#####################################################################
# Nmass -- climate moderators
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
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_mass" & 
                             manip_type == "n" & !is.na(gs_mat) & 
                             logr > -0.2),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.25, 0.5), breaks = seq(-0.25, 0.5, 0.25)) +
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
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_mass" & 
                             manip_type == "n" & !is.na(gs_ai) & 
                             logr > -0.2),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red") +
  scale_x_continuous(limits = c(0, 3.4), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.25, 0.5), breaks = seq(-0.25, 0.5, 0.25)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("MI"["g"]*" (unitless)")),
       y = expression(bold("N"["mass"]*" response to N addition (%)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
nmass_ai_plot 

# Nmass - PAR plot
nmass_par_plot <- ggplot(data = subset(meta_results, myvar == "leaf_n_mass" & 
                                         manip_type == "n" & !is.na(gs_par) & 
                                         logr > -0.2),
                         aes(x = gs_par, y = logr, size = 1/logr_se)) +
  geom_point(alpha = 0.30) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.25, 0.5), breaks = seq(-0.25, 0.5, 0.25)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("PAR"["g"]*" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       y = expression(bold("N"["mass"]*" response to N addition (lnRR)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
nmass_par_plot 

#####################################################################
# Narea -- climate moderators
#####################################################################
ggplot(data = meta_results %>% filter(manip_type == "n" & 
                                        myvar == "leaf_n_area" & 
                                        !is.na(gs_mat))) +
  geom_point(aes(x = gs_mat, y = logr))


# Narea climate moderator model
nadd_narea_clim <- rma.mv(logr, logr_var, method = "REML", 
                          random = ~ 1 | exp, mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, control = list(stepadj = 0.3), 
                          data = meta_results %>% filter(manip_type == "n" & 
                                                           myvar == "leaf_n_area" & 
                                                           !is.na(gs_mat)))
summary(nadd_narea_clim)

# Narea - temperature plot
narea_tg_plot <- mod_results(nadd_narea_clim, mod = "gs_mat",
                             group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_area" & 
                             manip_type == "n" & !is.na(gs_mat)),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("N"["area"]*" response to N addition (lnRR)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
narea_tg_plot 

# Narea - aridity plot
narea_ai_plot <- ggplot() +
  geom_point(data = subset(meta_results, myvar == "leaf_n_area" & 
                             manip_type == "n" & !is.na(gs_ai) & 
                             logr > -0.2),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  scale_x_continuous(limits = c(0, 3.4), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("MI"["g"]*" (unitless)")),
       y = expression(bold("N"["area"]*" response to N addition (%)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
narea_ai_plot 

# Narea - PAR plot
narea_par_plot <- ggplot(data = subset(meta_results, myvar == "leaf_n_area" & 
                                         manip_type == "n" & !is.na(gs_par)),
                         aes(x = gs_par, y = logr, size = 1/logr_se)) +
  geom_point(alpha = 0.30) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("PAR"["g"]*" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       y = expression(bold("N"["mass"]*" response to N addition (lnRR)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
narea_par_plot 


#####################################################################
# Pmass -- climate moderators
#####################################################################
ggplot(data = meta_results %>% filter(manip_type == "p" & 
                                        myvar == "leaf_p_mass" & 
                                        !is.na(gs_mat))) +
  geom_point(aes(x = gs_mat, y = logr))

# Pmass climate moderator model
padd_pmass_clim <- rma.mv(logr, logr_var, method = "REML", 
                          random = ~ 1 | exp, mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, control = list(stepadj = 0.3), 
                          data = meta_results %>% filter(manip_type == "p" & 
                                                           myvar == "leaf_p_mass" & 
                                                           !is.na(gs_mat) & gs_ai < 3))
summary(padd_pmass_clim)

# Pmass - temperature plot
pmass_tg_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_mass" & 
                             manip_type == "p" & !is.na(gs_mat) & 
                             gs_ai < 3),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold("P"["mass"]*" response to P addition (lnRR)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
pmass_tg_plot 

# Pmass - aridity plot
pmass_ai_plot <- mod_results(padd_pmass_clim, mod = "gs_ai",
                             group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_mass" & 
                             manip_type == "p" & !is.na(gs_ai) & gs_ai < 3),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("MI"["g"]*" (unitless)")),
       y = expression(bold("P"["mass"]*" response to P addition (%)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
pmass_ai_plot 

# Pmass - PAR plot
pmass_par_plot <- ggplot(data = subset(meta_results, myvar == "leaf_p_mass" & 
                                         manip_type == "p" & !is.na(gs_par) & 
                                         gs_ai < 3),
                         aes(x = gs_par, y = logr, size = 1/logr_se)) +
  geom_point(alpha = 0.30) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("PAR"["g"]*" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       y = expression(bold("P"["mass"]*" response to P addition (lnRR)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
pmass_par_plot 

#####################################################################
# Parea -- climate moderators
#####################################################################
ggplot(data = meta_results %>% filter(manip_type == "p" & 
                                        myvar == "leaf_p_area" & 
                                        !is.na(gs_mat)& logr >-1 & gs_ai < 3)) +
  geom_point(aes(x = gs_mat, y = logr))


# Narea climate moderator model
padd_parea_clim <- rma.mv(logr, logr_var, method = "REML", 
                          random = ~ 1 | exp, mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, control = list(stepadj = 0.3), 
                          data = meta_results %>% filter(manip_type == "p" & 
                                                           myvar == "leaf_p_area" & 
                                                           !is.na(gs_mat) & logr >-1 & gs_ai < 3))
summary(padd_parea_clim)

# Narea - temperature plot
parea_tg_plot <- mod_results(padd_parea_clim, mod = "gs_mat",
                             group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             manip_type == "p" & !is.na(gs_mat)  & logr >-1 & gs_ai < 3),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, color = "blue") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = expression(bold(bolditalic("P")["area"]*" response to P addition (lnRR)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
parea_tg_plot 

# Narea - aridity plot
parea_ai_plot <- mod_results(padd_parea_clim, mod = "gs_ai",
                             group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             manip_type == "p" & !is.na(gs_mat)  & logr >-1 & gs_ai < 3),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             manip_type == "p" & !is.na(gs_ai) & 
                             logr > -1 & gs_ai < 3),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("MI"["g"]*" (unitless)")),
       y = expression(bold("P"["area"]*" response to P addition (%)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
parea_ai_plot 

# Parea - PAR plot
parea_par_plot <- mod_results(padd_parea_clim, mod = "gs_par",
                              group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             manip_type == "p" & !is.na(gs_mat)  & logr >-1 & gs_ai < 3),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             manip_type == "p" & !is.na(gs_ai) & 
                             logr > -1 & gs_ai < 3),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(x = expression(bold("PAR"["g"]*" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       y = expression(bold("P"["area"]*" response to P addition (lnRR)")),
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
parea_par_plot 
































# Nmass - photosynthetic pathway
nmass_photo_plot <- ggplot(data = meta_photo_results %>% 
                             filter(trait == "nmass" &
                                      nut_add == "n"),
                           aes(x = photo, y = estimate)) +
  geom_errorbar(aes(ymin = lowerCL, 
                    ymax = upperCL),
                width = 0.3, linewidth = 1) +
  geom_point(aes(fill = photo), shape = 21, size = 5) +
  geom_bracket(xmin = "C3", xmax = "C4", y.position = 0.5, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(x = 1.5, y = 0.52, label = "**", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("pink", "darkred")) +
  scale_y_continuous(limits = c(-0.25, 0.5), breaks = seq(-0.25, 0.5, 0.25)) +
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
                           aes(x = nfix, y = estimate)) +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL),
                width = 0.3, linewidth = 1) +
  geom_point(aes(fill = nfix), shape = 21, size = 5) +
  geom_bracket(xmin = "No", xmax = "Yes", y.position = 0.5, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(x = 1.5, y = 0.52, label = "***", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("pink", "darkred")) +
  scale_x_discrete(labels = c("non-fixer", expression("N"["2"]*"-fixer"))) +
  scale_y_continuous(limits = c(-0.25, 0.5), breaks = seq(-0.25, 0.5, 0.25)) +
  labs(x = expression(bold("N"["2"]*" fixation ability")),
       y = expression(bold("N"["mass"]*" response to N addition"))) +
  guides(fill = "none") +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
nmass_nfix_plot

# Nmass - mycorrhizal acquisition strategy
nmass_myc_plot <- ggplot(data = meta_myc_results %>% 
                            filter(trait == "nmass" &
                                     nut_add == "n"),
                          aes(x = myc, y =  estimate)) +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL),
                width = 0.3, linewidth = 1) +
  geom_point(aes(fill = myc), shape = 21, size = 5) +
  geom_bracket(xmin = "Mining", xmax = "Scavenging", y.position = 0.5, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(x = 1.5, y = 0.52, label = "NS", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("pink", "darkred")) +
  scale_x_discrete(labels = c("mining", "scavenging")) +
  scale_y_continuous(limits = c(-0.25, 0.5), breaks = seq(-0.25, 0.5, 0.25)) +
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
  labs(x = expression(bold("MI"["g"]*" (unitless)")),
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
# Leaf N:P interaction -- climate and species identity moderators
#####################################################################
# Climate Model
int_leafnp_clim <- rma.mv(yi = dNPi, V = vNPi, W = wNPi, method = "REML", 
                          random = ~ 1 | exp, mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, control = list(stepadj = 0.3), 
                          data = meta_results_int %>% 
                            filter(myvar == "leaf_np" & !is.na(gs_mat) & dNPi < 1))
summary(int_leafnp_clim)

# Temperature plot
leafnp_int_tg_plot <- mod_results(int_leafnp_clim, 
                                  mod = "gs_mat", 
                                  group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = (exp(estimate) - 1) * 100)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = subset(meta_results_int, myvar == "leaf_np" & dNPi < 1), 
             aes(x = gs_mat, y = (exp(dNPi) - 1) * 100,
             size = 1 / dNPi_se)) +
  geom_ribbon(aes(ymax = (exp(upperCL) - 1) * 100, ymin = (exp(lowerCL) - 1) * 100),
              alpha = 0.3) +
  geom_smooth(linewidth = 1.5, color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-100, 200)) +
  scale_size_continuous(limits = c(0, 120), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = "Leaf N:P interaction effect (%)",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
leafnp_int_tg_plot

# Moisture index plot
leafnp_int_par_plot <- mod_results(int_leafnp_clim, 
                                  mod = "gs_par", 
                                  group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = (exp(estimate) - 1) * 100)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = subset(meta_results_int, myvar == "leaf_np" & dNPi < 1), 
             aes(x = gs_par, y = (exp(dNPi) - 1) * 100,
                 size = 1 / dNPi_se)) +
  geom_ribbon(aes(ymax = (exp(upperCL) - 1) * 100, ymin = (exp(lowerCL) - 1) * 100),
              alpha = 0.3) +
  geom_smooth(linewidth = 1.5, color = "black") +
  scale_y_continuous(limits = c(-100, 200)) +
  scale_x_continuous(limits = c(600, 1000), breaks = seq(600, 1000, 100)) +
  scale_size_continuous(limits = c(0, 120), range = c(1, 7)) +
  labs(x = expression(bold("PAR"["g"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       y = "Leaf N:P interaction effect (%)",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
leafnp_int_par_plot

#####################################################################
# AGB interaction -- climate moderators
#####################################################################

# Climate Model
int_agb_clim <- rma.mv(yi = dNPi, V = vNPi, W = wNPi, method = "REML", 
                       random = ~ 1 | exp, mods = ~ gs_mat + gs_ai + gs_par,
                       slab = exp, control = list(stepadj = 0.3), 
                       data = meta_results_int %>% 
                         filter(myvar == "agb" & !is.na(gs_mat) & dNPi < 2))
summary(int_agb_clim)


# Temperature plot
agb_int_tg_plot <- mod_results(int_agb_clim, 
                               mod = "gs_mat", 
                               group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = (exp(estimate) - 1) * 100)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = subset(meta_results_int, myvar == "agb" & dNPi < 2), 
             aes(x = gs_mat, y = (exp(dNPi) - 1) * 100,
                 size = 1 / dNPi_se)) +
  geom_ribbon(aes(ymax = (exp(upperCL) - 1) * 100, ymin = (exp(lowerCL) - 1) * 100),
              alpha = 0.3) +
  geom_smooth(linewidth = 1.5, color = "black") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-100, 400)) +
  scale_size_continuous(limits = c(0, 120), range = c(1, 7)) +
  labs(x = expression(bold("T"["g"]*" ("*degree*"C)")),
       y = "AGB interaction effect (%)",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))

# Moisture index plot
agb_int_ai_plot <- mod_results(int_agb_clim, 
                                  mod = "gs_ai", 
                                  group = "exp")$mod_table %>%
  ggplot(aes(x = moderator, y = (exp(estimate) - 1) * 100)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = subset(meta_results_int, myvar == "agb" & dNPi < 2), 
             aes(x = gs_ai, y = (exp(dNPi) - 1) * 100,
                 size = 1 / dNPi_se)) +
  geom_ribbon(aes(ymax = (exp(upperCL) - 1) * 100, ymin = (exp(lowerCL) - 1) * 100),
              alpha = 0.3) +
  geom_smooth(linewidth = 1.5, color = "black") +
  scale_y_continuous(limits = c(-100, 400)) +
  scale_x_continuous(limits = c(0, 3.4), breaks = seq(0, 3, 1)) +
  scale_size_continuous(limits = c(0, 120), range = c(1, 7)) +
  labs(x = expression(bold("MI"["g"]*" (unitless)")),
       y = "AGB interaction effect (%)",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))


#####################################################################
# Write plots
#####################################################################

# Figure 1 - individual effects
png("../plots/CNP_fig1_indEffects2.png", height = 12, width = 16, 
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

# Figure 2 - moderator plots
png("../plots/CNP_fig2_moderators_ind.png", height = 12.5, width = 24,
    units = "in", res = 600)
ggarrange(nmass_tg_plot, nmass_ai_plot, nmass_nfix_plot, nmass_myc_plot, nmass_photo_plot,
          pmass_tg_plot, pmass_ai_plot, pmass_nfix_plot, pmass_myc_plot, pmass_photo_plot,
          nrow = 2, ncol = 5, legend = "bottom", common.legend = TRUE,
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", 
                     "(f)", "(g)", "(h)", "(i)", "(j)"),
          font.label = list(size = 20), align = "hv")
dev.off()

# Figure 3 - interaction effects
png("../plots/CNP_fig3_intEffects.png", height = 12, width = 5.5, 
    units = "in", res = 600)
ggarrange(npint_chemistry_plot, npint_photo_plot, npint_bio_plot, 
          nrow = 3, ncol = 1, labels = c("(a)", "(b)", "(c)"), 
          font.label = list(size = 18), align = "hv")
dev.off()

# Figure 4 - climate interactions plots
png("../plots/CNP_fig4_climate_mod_int.png", height = 12, width = 12,
    units = "in", res = 600)
ggarrange(leafnp_int_tg_plot, leafnp_int_par_plot, 
          agb_int_tg_plot, agb_int_ai_plot, 
          nrow = 2, ncol = 2, legend = "bottom", common.legend = TRUE,
          labels = c("(a)", "(b)", "(c)", "(d)"),
          font.label = list(size = 20), align = "hv")
dev.off()

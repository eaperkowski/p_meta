# Libraries
library(tidyverse)
library(metafor)
library(ggpubr)
library(naniar) # to resolve NA/<NA> issue
library(orchaRd)

# Load meta-analysis results
meta_results <- read.csv("../data/CNPmeta_logr_results.csv")
meta_results_int <- read.csv("../data/CNPmeta_logr_results_int.csv")

# Check file structure
head(meta_results)
head(meta_results_int)

##############################################################################
# Marea climate moderators
##############################################################################
###############
# N addition
###############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "lma" & 
                                        !is.na(gs_mat) & logr > -0.5 & 
                                        logr < 0.4 & gs_ai < 3)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
nadd_marea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "n" & 
                                     myvar == "lma" & 
                                     !is.na(gs_mat) & logr > -0.5 & 
                                     logr < 0.4 & gs_ai < 3))

# Climate summary
nadd_marea_clim_summary <- data.frame(trait = "marea",
                                      nut_add = "n",
                                      k = 78,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_marea_clim)),
                                      row.names = NULL)

# Temperature plot
nadd_marea_tg_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "lma" & 
                             nut_add == "n" & !is.na(gs_mat) & 
                             logr > -0.5 & logr < 0.4 & gs_ai < 3),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("M")[bold("area")]*bold(" N response to "*bolditalic("T")[bold("g")])),
       x = "",
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_marea_tg_plot 

# Aridity plot
nadd_marea_ai_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "lma" & 
                             nut_add == "n" & !is.na(gs_mat) & 
                             logr > -0.5 & logr < 0.4 & gs_ai < 3),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("M")[bold("area")]*bold(" N response to "*bolditalic("MI")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_marea_ai_plot 

# PAR plot
nadd_marea_par_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "lma" & 
                             nut_add == "n" & !is.na(gs_mat) & 
                             logr > -0.5 & logr < 0.4 & gs_ai < 3),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("M")[bold("area")]*bold(" N response to "*bolditalic("PAR")[bold("g")])),,
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_marea_par_plot 

###############
# P addition
###############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "lma" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr > -0.5 & logr < 0.4)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
padd_marea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "p" & 
                                     myvar == "lma" & 
                                     !is.na(gs_mat) & gs_ai < 3 & 
                                     logr > -0.5 & logr < 0.4))

# Climate summary
padd_marea_clim_summary <- data.frame(trait = "marea",
                                      nut_add = "p",
                                      k = 75,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_marea_clim)),
                                      row.names = NULL)

# Temperature plot
padd_marea_tg_plot <- mod_results(padd_marea_clim, mod = "gs_mat",
                                  group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "lma" & 
                             nut_add == "p" & !is.na(gs_mat)  & gs_ai < 3 & 
                             logr > -0.5 & logr < 0.4),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, color = "blue", linetype = "dashed") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("M")[bold("area")]*bold(" P response to "*bolditalic("T")[bold("g")])),,
       x = "",
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_marea_tg_plot 

# Aridity plot
padd_marea_ai_plot <- mod_results(padd_marea_clim, mod = "gs_ai",
                                  group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "lma" & 
                             nut_add == "p" & !is.na(gs_mat)  & gs_ai < 3 & 
                             logr > -0.5 & logr < 0.4),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, color = "blue") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("M")[bold("area")]*bold(" P response to "*bolditalic("MI")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 20))
padd_marea_ai_plot 

# PAR plot
padd_marea_par_plot <- mod_results(padd_marea_clim, mod = "gs_par",
                                 group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "lma" & 
                             nut_add == "p" & !is.na(gs_mat)  & gs_ai < 3 & 
                             logr > -0.5 & logr < 0.4),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, color = "blue") +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("M")[bold("area")]*bold(" P response to "*bolditalic("PAR")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_marea_par_plot 

###############
# N+P addition
###############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "lma" & 
                                        !is.na(gs_mat) & gs_ai < 3 & 
                                        logr > -0.5 & logr < 0.4)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
npadd_marea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "np" & 
                                     myvar == "lma" & 
                                     !is.na(gs_mat) & gs_ai < 3 & logr > -0.5 & 
                                     logr < 0.4))

# Climate summary
npadd_marea_clim_summary <- data.frame(trait = "marea",
                                       nut_add = "np",
                                       k = 75,
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_marea_clim)),
                                       row.names = NULL)

# Temperature plot
npadd_marea_tg_plot <- mod_results(npadd_marea_clim, mod = "gs_mat",
                                  group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "lma" & 
                             nut_add == "np" & !is.na(gs_mat) & gs_ai < 3 & 
                             logr > -0.5 & logr < 0.4),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "magenta") +
  geom_smooth(method = "loess", linewidth = 2, color = "magenta", 
              linetype = "dashed") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("M")[bold("area")]*bold(" N+P response to "*bolditalic("T")[bold("g")])),
       x = expression(bolditalic("T")[bold("g")]*bold(" ("*degree*"C)")),
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_marea_tg_plot 

# Aridity plot
npadd_marea_ai_plot <- mod_results(npadd_marea_clim, mod = "gs_ai",
                                  group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "lma" & 
                             nut_add == "np" & !is.na(gs_mat)  & gs_ai < 3 & 
                             logr > -0.5 & logr < 0.4),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, shape = 21, fill = "magenta") +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "magenta") +
  geom_smooth(method = "loess", linewidth = 2, color = "magenta", 
              linetype = "dashed") +
  scale_fill_manual(values = "magenta") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("M")[bold("area")]*bold(" N+P response to "*bolditalic("MI")[bold("g")])),,
       x = expression(bolditalic("MI")[bold("g")]*bold(" (unitless)")),
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_marea_ai_plot 

# PAR plot
npadd_marea_par_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "lma" & 
                             nut_add == "np" & !is.na(gs_mat)  & gs_ai < 3 & 
                             logr > -0.5 & logr < 0.4),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("M")[bold("area")]*bold(" N+P response to "*bolditalic("PAR")[bold("g")])),
       x = expression(bolditalic("PAR")[bold("g")]*bold(" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       y = "",
       size = expression(bold("Error"^"-1"))) +
  theme_classic(base_size = 20) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_marea_par_plot

###################
# Merge Marea results, plots
###################
# Merge Marea moderator results, with some light cleaning
marea_clim_summary <- rbind(nadd_marea_clim_summary, 
                            padd_marea_clim_summary, 
                            npadd_marea_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

# Create plot
png("../plots/supp/CNP_figSX_marea_climate.png", height = 14, width = 14,
    units = "in", res = 600)
ggarrange(nadd_marea_tg_plot, nadd_marea_ai_plot, nadd_marea_par_plot,
          padd_marea_tg_plot, padd_marea_ai_plot, padd_marea_par_plot,
          npadd_marea_tg_plot, npadd_marea_ai_plot, npadd_marea_par_plot,
          ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"),
          font.label = list(size = 22))
dev.off()

##############################################################################
# Nmass climate moderators
##############################################################################
###############
# N addition
###############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "leaf_n_mass" & 
                                        !is.na(gs_mat) & gs_ai < 3)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
nadd_nmass_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "n" & 
                                     myvar == "leaf_n_mass" & 
                                     !is.na(gs_mat) & gs_ai < 3))
# Climate summary
nadd_nmass_clim_summary <- data.frame(trait = "nmass",
                                      nut_add = "n",
                                      k = 103,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_nmass_clim)),
                                      row.names = NULL)

# Temperature plot
nadd_nmass_tg_plot <- mod_results(nadd_nmass_clim, mod = "gs_mat",
                                  group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_mass" & 
                             nut_add == "n" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.3, 0.6), breaks = seq(-0.3, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("mass")]*bold(" N response to "*bolditalic("T")[bold("g")])),,
       x = "",
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_nmass_tg_plot 

# Aridity plot
nadd_nmass_ai_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_mass" & 
                             nut_add == "n" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.3, 0.6), breaks = seq(-0.3, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("mass")]*bold(" N response to "*bolditalic("MI")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_nmass_ai_plot 

# PAR plot
nadd_nmass_par_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_mass" & 
                             nut_add == "n" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.3, 0.6), breaks = seq(-0.3, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("mass")]*bold(" N response to "*bolditalic("PAR")[bold("g")])),,
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_nmass_par_plot 

###############
# P addition
###############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "leaf_n_mass" & 
                                        !is.na(gs_mat) & gs_ai < 3 & 
                                        logr < 0.5 & logr > -0.45)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
padd_nmass_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "p" & 
                                     myvar == "leaf_n_mass" & 
                                     !is.na(gs_mat) & gs_ai < 3 & 
                                     logr < 0.5 & logr > -0.45))

# Create climate summary
padd_nmass_clim_summary <- data.frame(trait = "nmass",
                                      nut_add = "p",
                                      k = 101,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_nmass_clim)),
                                      row.names = NULL)

# Temperature plot
padd_nmass_tg_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_mass" & 
                             nut_add == "p" & !is.na(gs_mat) & 
                             logr > -0.5 & logr < 0.4 & gs_ai < 3),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("mass")]*bold(" P response to "*bolditalic("T")[bold("g")])),
       x = "",
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_nmass_tg_plot 

# Aridity plot
padd_nmass_ai_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_mass" & 
                             nut_add == "p" & !is.na(gs_mat) & 
                             logr > -0.5 & logr < 0.4 & gs_ai < 3),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("mass")]*bold(" P response to "*bolditalic("MI")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_nmass_ai_plot 

# PAR plot
padd_nmass_par_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_mass" & 
                             nut_add == "p" & !is.na(gs_mat) & 
                             logr > -0.5 & logr < 0.4 & gs_ai < 3),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("mass")]*bold(" P response to "*bolditalic("PAR")[bold("g")])),,
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_nmass_par_plot 

###############
# N+P addition
###############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "leaf_n_mass" & 
                                        !is.na(gs_mat) & gs_ai < 3 & 
                                        logr < 0.7 & logr > -0.3)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
npadd_nmass_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(nut_add == "np" & 
                                      myvar == "leaf_n_mass" & 
                                      !is.na(gs_mat) & gs_ai < 3 & 
                                      logr < 0.7 & logr > -0.3))

# Create climate summary
npadd_nmass_clim_summary <- data.frame(trait = "nmass",
                                       nut_add = "np",
                                       k = 101,
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_nmass_clim)),
                                       row.names = NULL)

# Temperature plot
npadd_nmass_tg_plot <- mod_results(npadd_nmass_clim, mod = "gs_mat",
                                  group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_mass" & 
                             nut_add == "np" & !is.na(gs_mat) & gs_ai < 3 & 
                             logr < 0.7 & logr > -0.3),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "magenta") +
  geom_smooth(method = "loess", linewidth = 2, color = "magenta") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.3, 0.6), breaks = seq(-0.3, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("mass")]*bold(" N+P response to "*bolditalic("T")[bold("g")])),
       x = expression(bolditalic("T")[bold("g")]*bold(" ("*degree*"C)")),
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_nmass_tg_plot

# Aridity plot
npadd_nmass_ai_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_mass" & 
                             nut_add == "np" & !is.na(gs_mat) & gs_ai < 3 & 
                             logr < 0.7 & logr > -0.3),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.3, 0.6), breaks = seq(-0.3, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("mass")]*bold(" N+P response to "*bolditalic("MI")[bold("g")])),
       x = expression(bolditalic("MI")[bold("g")]*bold(" (unitless)")),
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_nmass_ai_plot 

# PAR plot
npadd_nmass_par_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_mass" & 
                             nut_add == "np" & !is.na(gs_mat) & gs_ai < 3 & 
                             logr < 0.7 & logr > -0.3),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.3, 0.6), breaks = seq(-0.3, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("mass")]*bold(" N+P response to "*bolditalic("PAR")[bold("g")])),,
       x = expression(bolditalic("PAR")[bold("g")]*bold(" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_nmass_par_plot

###################
# Merge Nmass results, plots
###################
# Merge Nmass moderator results, with some light cleaning
nmass_clim_summary <- rbind(nadd_nmass_clim_summary, 
                            padd_nmass_clim_summary, 
                            npadd_nmass_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

# Create plot
png("../plots/supp/CNP_figSX_nmass_climate.png", height = 14, width = 14,
    units = "in", res = 600)
ggarrange(nadd_nmass_tg_plot, nadd_nmass_ai_plot, nadd_nmass_par_plot,
          padd_nmass_tg_plot, padd_nmass_ai_plot, padd_nmass_par_plot,
          npadd_nmass_tg_plot, npadd_nmass_ai_plot, npadd_nmass_par_plot,
          ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"),
          font.label = list(size = 22))
dev.off()

##############################################################################
# Narea climate moderators
##############################################################################
###############
# N addition
###############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "leaf_n_area" & 
                                        !is.na(gs_mat) & gs_ai < 3 & 
                                        logr > -1 & logr < 0.95)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
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

# Climate summary
nadd_narea_clim_summary <- data.frame(trait = "narea",
                                      nut_add = "n",
                                      k = 49,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_narea_clim)),
                                      row.names = NULL)

# Temperature plot
nadd_narea_tg_plot <- mod_results(nadd_narea_clim, mod = "gs_mat",
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
  labs(title = expression(bolditalic("N")[bold("area")]*bold(" N response to "*bolditalic("T")[bold("g")])),
       x = "",
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_narea_tg_plot 

# Aridity plot
nadd_narea_ai_plot <- mod_results(nadd_narea_clim, mod = "gs_ai",
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
  labs(title = expression(bolditalic("N")[bold("area")]*bold(" N response to "*bolditalic("MI")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_narea_ai_plot 

# PAR plot
nadd_narea_par_plot <- ggplot(data = subset(meta_results, myvar == "leaf_n_area" & 
                                         nut_add == "n" & !is.na(gs_par) & gs_ai < 3 & 
                                         logr > -1 & logr < 0.95),
                         aes(x = gs_par, y = logr, size = 1/logr_se)) +
  geom_point(alpha = 0.30, fill = "red", shape = 21) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("area")]*bold(" N response to "*bolditalic("PAR")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_narea_par_plot 

###############
# P addition
###############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "leaf_n_area" & 
                                        !is.na(gs_mat) & gs_ai < 3 & 
                                        logr > -1 & logr < 0.6)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
padd_narea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "p" & 
                                     myvar == "leaf_n_area" & 
                                     !is.na(gs_mat) & gs_ai < 3 & logr > -1))

# Climate summary
padd_narea_clim_summary <- data.frame(trait = "narea",
                                      nut_add = "p",
                                      k = 49,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_narea_clim)),
                                      row.names = NULL)

# Temperature plot
padd_narea_tg_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_area" & 
                             nut_add == "p" & !is.na(gs_mat) &
                             gs_ai < 3 & logr > -1),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("area")]*bold(" P response to "*bolditalic("T")[bold("g")])),
       x = "",
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_narea_tg_plot 

# Aridity plot
padd_narea_ai_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_area" & 
                             nut_add == "p" & !is.na(gs_mat) & 
                             gs_ai < 3 & logr > -1),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("area")]*bold(" P response to "*bolditalic("MI")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_narea_ai_plot 

# PAR plot
padd_narea_par_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_area" & 
                             nut_add == "p" & !is.na(gs_mat) & 
                             gs_ai < 3 & logr > -1),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("area")]*bold(" P response to "*bolditalic("PAR")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_narea_par_plot 

###############
# N+P addition
###############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "leaf_n_area" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr > -0.2)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
npadd_narea_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(nut_add == "np" & 
                                      myvar == "leaf_n_area" & 
                                      !is.na(gs_mat) & gs_ai < 3 & logr > -0.2))

# Climate summary
npadd_narea_clim_summary <- data.frame(trait = "narea",
                                       nut_add = "np",
                                       k = 46,
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_narea_clim)),
                                       row.names = NULL)

# Temperature plot
npadd_narea_tg_plot <- mod_results(npadd_narea_clim, mod = "gs_mat",
                                   group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_area" & 
                             nut_add == "np" & !is.na(gs_mat) & 
                             gs_ai < 3 & logr > -0.2),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "magenta") +
  geom_smooth(method = "loess", linewidth = 2, color = "magenta", linetype = "dashed") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.3, 0.6), breaks = seq(-0.3, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("area")]*bold(" N+P response to "*bolditalic("T")[bold("g")])),
       x = expression(bolditalic("T")[bold("g")]*bold(" ("*degree*"C)")),
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_narea_tg_plot 

# Aridity plot
npadd_narea_ai_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_area" & 
                             nut_add == "np" & !is.na(gs_mat) & 
                             gs_ai < 3 & logr > -0.2),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.3, 0.6), breaks = seq(-0.3, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("area")]*bold(" N+P response to "*bolditalic("MI")[bold("g")])),
       x = expression(bolditalic("MI")[bold("g")]*bold(" (unitless)")),
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_narea_ai_plot 

# PAR plot
npadd_narea_par_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_n_area" & 
                             nut_add == "np" & !is.na(gs_mat) & 
                             gs_ai < 3 & logr > -0.2),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.3, 0.6), breaks = seq(-0.3, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("N")[bold("area")]*bold(" N+P response to "*bolditalic("PAR")[bold("g")])),
       x = expression(bolditalic("PAR")[bold("g")]*bold(" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_narea_par_plot 

###################
# Merge Narea results, plots
###################
# Merge Narea moderator results, with some light cleaning
narea_clim_summary <- rbind(nadd_narea_clim_summary, 
                            padd_narea_clim_summary, 
                            npadd_narea_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

# Create plot
png("../plots/supp/CNP_figSX_narea_climate.png", height = 14, width = 14,
    units = "in", res = 600)
ggarrange(nadd_narea_tg_plot, nadd_narea_ai_plot, nadd_narea_par_plot,
          padd_narea_tg_plot, padd_narea_ai_plot, padd_narea_par_plot,
          npadd_narea_tg_plot, npadd_narea_ai_plot, npadd_narea_par_plot,
          ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"),
          font.label = list(size = 22))
dev.off()

##############################################################################
# Pmass climate moderators
##############################################################################
###############
# N addition
###############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "leaf_p_mass" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr > -0.7)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
nadd_pmass_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "n" & 
                                     myvar == "leaf_p_mass" & 
                                     !is.na(gs_mat) & gs_ai < 3 & logr > -0.7))

# Climate summary
nadd_pmass_clim_summary <- data.frame(trait = "pmass",
                                      nut_add = "n",
                                      k = 95,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_pmass_clim)),
                                      row.names = NULL)

# Temperature plot
nadd_pmass_tg_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_mass" & 
                             nut_add == "n" & !is.na(gs_mat) & gs_ai < 3 & logr > -0.7),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("mass")]*bold(" N response to "*bolditalic("T")[bold("g")])),
       x = "",
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_pmass_tg_plot 

# Aridity plot
nadd_pmass_ai_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_mass" & 
                             nut_add == "n" & !is.na(gs_mat) & gs_ai < 3 & logr > -0.7),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("mass")]*bold(" N response to "*bolditalic("MI")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_pmass_ai_plot 

# PAR plot
nadd_pmass_par_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_mass" & 
                             nut_add == "n" & !is.na(gs_mat) & gs_ai < 3 & logr > -0.7),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, 0.3)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("mass")]*bold(" N response to "*bolditalic("PAR")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_pmass_par_plot 

###############
# P addition
###############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "leaf_p_mass" & 
                                        !is.na(gs_mat) & gs_ai < 3)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
padd_pmass_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "p" & 
                                     myvar == "leaf_p_mass" & 
                                     !is.na(gs_mat) & gs_ai < 3))

# Climate summary
padd_pmass_clim_summary <- data.frame(trait = "pmass",
                                      nut_add = "p",
                                      k = 97,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_pmass_clim)),
                                      row.names = NULL)

# Temperature plot
padd_pmass_tg_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_mass" & 
                             nut_add == "p" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("mass")]*bold(" P response to "*bolditalic("T")[bold("g")])),
       x = "",
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_pmass_tg_plot 

# Aridity plot
padd_pmass_ai_plot <- mod_results(padd_pmass_clim, mod = "gs_ai",
                                  group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_mass" & 
                             nut_add == "p" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, color = "blue", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("mass")]*bold(" P response to "*bolditalic("MI")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_pmass_ai_plot 

# PAR plot
padd_pmass_par_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_mass" & 
                             nut_add == "p" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("mass")]*bold(" P response to "*bolditalic("PAR")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_pmass_par_plot 

###############
# N+P addition
###############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "leaf_p_mass" & 
                                        !is.na(gs_mat) & gs_ai < 3)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
npadd_pmass_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(nut_add == "np" & 
                                      myvar == "leaf_p_mass" & 
                                      !is.na(gs_mat) & gs_ai < 3))

# Climate summary
npadd_pmass_clim_summary <- data.frame(trait = "pmass",
                                       nut_add = "np",
                                       k = 97,
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_pmass_clim)),
                                       row.names = NULL)

# Temperature plot
npadd_pmass_tg_plot <- mod_results(npadd_parea_clim, mod = "gs_mat",
                                  group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_mass" & 
                             nut_add == "np" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "magenta") +
  geom_smooth(method = "loess", linewidth = 2, color = "magenta", linetype = "dashed") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("mass")]*bold(" N+P response to "*bolditalic("T")[bold("g")])),
       x = expression(bolditalic("T")[bold("g")]*bold(" ("*degree*"C)")),
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_pmass_tg_plot 

# Aridity plot
npadd_pmass_ai_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_mass" & 
                             nut_add == "np" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("mass")]*bold(" N+P response to "*bolditalic("MI")[bold("g")])),
       x = expression(bolditalic("MI")[bold("g")]*bold(" (unitless)")),
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_pmass_ai_plot 

# PAR plot
npadd_pmass_par_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_mass" & 
                             nut_add == "np" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("mass")]*bold(" N+P response to "*bolditalic("PAR")[bold("g")])),
       x = expression(bolditalic("PAR")[bold("g")]*bold(" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_pmass_par_plot 

###################
# Merge Pmass results, plots
###################
# Merge Pmass moderator results, with some light cleaning
pmass_clim_summary <- rbind(nadd_pmass_clim_summary, 
                            padd_pmass_clim_summary, 
                            npadd_pmass_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

# Create plot
png("../plots/supp/CNP_figSX_pmass_climate.png", height = 14, width = 14,
    units = "in", res = 600)
ggarrange(nadd_pmass_tg_plot, nadd_pmass_ai_plot, nadd_pmass_par_plot,
          padd_pmass_tg_plot, padd_pmass_ai_plot, padd_pmass_par_plot,
          npadd_pmass_tg_plot, npadd_pmass_ai_plot, npadd_pmass_par_plot,
          ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"),
          font.label = list(size = 22))
dev.off()

##############################################################################
# Parea climate moderators
##############################################################################
###############
# N addition
###############

# Data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "leaf_p_area" & 
                                        !is.na(gs_mat) & gs_ai < 3)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
nadd_parea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "n" & 
                                     myvar == "leaf_p_area" & 
                                     !is.na(gs_mat) & gs_ai < 3))

# Climate summary
nadd_parea_clim_summary <- data.frame(trait = "parea",
                                      nut_add = "n",
                                      k = 45,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_parea_clim)),
                                      row.names = NULL)

# Temperature plot
nadd_parea_tg_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             nut_add == "n" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("area")]*bold(" N response to "*bolditalic("T")[bold("g")])),
       x = "",
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_parea_tg_plot 

# Aridity plot
nadd_parea_ai_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             nut_add == "n" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("area")]*bold(" N response to "*bolditalic("MI")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_parea_ai_plot 

# PAR plot
nadd_parea_par_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             nut_add == "n" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("area")]*bold(" N response to "*bolditalic("PAR")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_parea_par_plot 

###############
# P addition
###############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "leaf_p_area" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr > -1)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
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

# Climate summary
padd_parea_clim_summary <- data.frame(trait = "parea",
                                      nut_add = "p",
                                      k = 44,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_parea_clim)),
                                      row.names = NULL)

# Temperature plot
padd_parea_tg_plot <- mod_results(padd_parea_clim, mod = "gs_mat",
                             group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             nut_add == "p" & !is.na(gs_mat) & logr >-1 & gs_ai < 3),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, color = "blue") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0.5, 2, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("area")]*bold(" P response to "*bolditalic("T")[bold("g")])),
       x = "",
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_parea_tg_plot 

# Aridity plot
padd_parea_ai_plot <- mod_results(padd_parea_clim, mod = "gs_ai",
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
  labs(title = expression(bolditalic("P")[bold("area")]*bold(" P response to "*bolditalic("MI")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_parea_ai_plot 

# PAR plot
padd_parea_par_plot <- mod_results(padd_parea_clim, mod = "gs_par",
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
  labs(title = expression(bolditalic("P")[bold("area")]*bold(" P response to "*bolditalic("PAR")[bold("g")])),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_parea_par_plot 

###############
# N+P addition
###############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "leaf_p_area" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr > -1 & logr < 1.25)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
npadd_parea_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(nut_add == "np" & 
                                      myvar == "leaf_p_area" & 
                                      !is.na(gs_mat) & gs_ai < 3 & logr > -1 & logr < 1.25))

# Climate summary
npadd_parea_clim_summary <- data.frame(trait = "parea",
                                       nut_add = "np",
                                       k = 43,
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_parea_clim)),
                                       row.names = NULL)

# Temperature plot
npadd_parea_tg_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             nut_add == "np" & !is.na(gs_mat) & 
                             gs_ai < 3 & logr > -1 & logr < 1.25),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 1.5), breaks = seq(-0.5, 1.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("area")]*bold(" N+P response to "*bolditalic("T")[bold("g")])),
       x = expression(bolditalic("T")[bold("g")]*bold(" ("*degree*"C)")),
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_parea_tg_plot 

# Aridity plot
npadd_parea_ai_plot <- mod_results(npadd_parea_clim, mod = "gs_ai",
                                  group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             nut_add == "np" & !is.na(gs_mat) & gs_ai < 3 & 
                             logr > -1 & logr < 1.25),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "magenta") +
  geom_smooth(method = "loess", linewidth = 2, color = "magenta") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 1.5), breaks = seq(-0.5, 1.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("area")]*bold(" N+P response to "*bolditalic("MI")[bold("g")])),
       x = expression(bolditalic("MI")[bold("g")]*bold(" (unitless)")),
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_parea_ai_plot 

# PAR plot
npadd_parea_par_plot <- mod_results(padd_parea_clim, mod = "gs_par",
                                   group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_p_area" & 
                             nut_add == "np" & !is.na(gs_mat) & 
                             gs_ai < 3 & logr > -1 & logr < 1.25),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "magenta") +
  geom_smooth(method = "loess", linewidth = 2, color = "magenta") +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 1.5), breaks = seq(-0.5, 1.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bolditalic("P")[bold("area")]*bold(" N+P response to "*bolditalic("PAR")[bold("g")])),
       x = expression(bolditalic("PAR")[bold("g")]*bold(" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_parea_par_plot 


###################
# Merge Parea results, plots
###################
# Merge Parea moderator results, with some light cleaning
parea_clim_summary <- rbind(nadd_parea_clim_summary, 
                            padd_parea_clim_summary, 
                            npadd_parea_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

# Create plot
png("../plots/supp/CNP_figSX_parea_climate.png", height = 14, width = 14,
    units = "in", res = 600)
ggarrange(nadd_parea_tg_plot, nadd_parea_ai_plot, nadd_parea_par_plot,
          padd_parea_tg_plot, padd_parea_ai_plot, padd_parea_par_plot,
          npadd_parea_tg_plot, npadd_parea_ai_plot, npadd_parea_par_plot,
          ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"),
          font.label = list(size = 22))
dev.off()

##############################################################################
# Leaf N:P climate moderators
##############################################################################
###############
# N addition
###############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "leaf_np" & 
                                        !is.na(gs_mat) & gs_ai < 3)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
nadd_leafnp_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "n" & 
                                     myvar == "leaf_np" & 
                                     !is.na(gs_mat) & gs_ai < 3))

# Climate summary
nadd_leafnp_clim_summary <- data.frame(trait = "leaf_np",
                                      nut_add = "n",
                                      k = 84,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_leafnp_clim)),
                                      row.names = NULL)

# Temperature plot
nadd_leafnp_tg_plot <- mod_results(nadd_leafnp_clim, mod = "gs_mat",
                                   group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_np" & 
                             nut_add == "n" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bold("Leaf N:P N response to ")*bolditalic("T")[bold("g")]),
       x = "",
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_leafnp_tg_plot 

# Aridity plot
nadd_leafnp_ai_plot <- mod_results(nadd_leafnp_clim, mod = "gs_ai",
                                  group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_np" & 
                             nut_add == "n" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "red") +
  geom_smooth(method = "loess", linewidth = 2, color = "red", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bold("Leaf N:P N response to ")*bolditalic("MI")[bold("g")]),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_leafnp_ai_plot 

# PAR plot
nadd_leafnp_par_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_np" & 
                             nut_add == "n" & !is.na(gs_mat) & gs_ai < 3),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "red", shape = 21) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bold("Leaf N:P N response to ")*bolditalic("PAR")[bold("g")]),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
nadd_leafnp_par_plot 

###############
# P addition
###############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "leaf_np" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr > -1.6)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
padd_leafnp_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "p" & 
                                     myvar == "leaf_np" & 
                                     !is.na(gs_mat) & gs_ai < 3 & logr > -1.6))

# Climate summary
padd_leafnp_clim_summary <- data.frame(trait = "leaf_np",
                                      nut_add = "p",
                                      k = 83,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_leafnp_clim)),
                                      row.names = NULL)

# Temperature plot
padd_leafnp_tg_plot <- mod_results(padd_leafnp_clim, mod = "gs_mat",
                                   group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_np" & 
                             nut_add == "p" & !is.na(gs_mat) & gs_ai < 3 & logr > -1.6),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, color = "blue") +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-2, 0.5), breaks = seq(-2, 0.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bold("Leaf N:P P response to ")*bolditalic("T")[bold("g")]),
       x = "",
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_leafnp_tg_plot 

# Aridity plot
padd_leafnp_ai_plot <- mod_results(padd_leafnp_clim, mod = "gs_ai",
                                   group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_np" & 
                             nut_add == "p" & !is.na(gs_mat) & gs_ai < 3 & logr > -1.6),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  geom_ribbon(aes(ymax = upperCL, ymin = lowerCL),
              alpha = 0.3, fill = "blue") +
  geom_smooth(method = "loess", linewidth = 2, color = "blue") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-2, 0.5), breaks = seq(-2, 0.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bold("Leaf N:P P response to ")*bolditalic("MI")[bold("g")]),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_leafnp_ai_plot 

# PAR plot
padd_leafnp_par_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_np" & 
                             nut_add == "p" & !is.na(gs_mat) & gs_ai < 3 & logr > -1.6),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "blue", shape = 21) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-2, 0.5), breaks = seq(-2, 0.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bold("Leaf N:P P response to ")*bolditalic("PAR")[bold("g")]),
       x = "",
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
padd_leafnp_par_plot 

###############
# N+P addition
###############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "leaf_np" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr > -1.45)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
npadd_leafnp_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(nut_add == "np" & 
                                      myvar == "leaf_np" & 
                                      !is.na(gs_mat) & gs_ai < 3 & logr > -1.45))

# Climate summary
npadd_leafnp_clim_summary <- data.frame(trait = "leaf_np",
                                       nut_add = "np",
                                       k = 83,
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_leafnp_clim)),
                                       row.names = NULL)

# Temperature plot
npadd_leafnp_tg_plot <-  ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_np" & 
                             nut_add == "np" & !is.na(gs_mat) & gs_ai < 3 & logr > -1.45),
             aes(x = gs_mat, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  scale_x_continuous(limits = c(5, 27), breaks = seq(5, 25, 5)) +
  scale_y_continuous(limits = c(-1.5, 0.5), breaks = seq(-1.5, 0.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bold("Leaf N:P N+P response to ")*bolditalic("T")[bold("g")]),
       x = expression(bolditalic("T")[bold("g")]*bold(" ("*degree*"C)")),
       y = "Log-response ratio",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_leafnp_tg_plot 

# Aridity plot
npadd_leafnp_ai_plot <- mod_results(padd_leafnp_clim, mod = "gs_ai",
                                   group = "exp", subset = TRUE)$mod_table %>%
  ggplot(aes(x = moderator, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_np" & 
                             nut_add == "np" & !is.na(gs_mat) & gs_ai < 3 & logr > -1.45),
             aes(x = gs_ai, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(-1.5, 0.5), breaks = seq(-1.5, 0.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bold("Leaf N:P N+P response to ")*bolditalic("MI")[bold("g")]),
       x = expression(bolditalic("MI")[bold("g")]*bold(" (unitless)")),
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 15.5),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_leafnp_ai_plot 

# PAR plot
npadd_leafnp_par_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(meta_results, myvar == "leaf_np" & 
                             nut_add == "np" & !is.na(gs_mat) & gs_ai < 3 & logr > -1.45),
             aes(x = gs_par, y = logr, size = 1/logr_se), 
             alpha = 0.30, fill = "magenta", shape = 21) +
  scale_x_continuous(limits = c(500, 1000), breaks = seq(500, 1000, 100)) +
  scale_y_continuous(limits = c(-1.5, 0.5), breaks = seq(-1.5, 0.5, 0.5)) +
  scale_size_continuous(limits = c(0, 224), range = c(1, 7)) +
  labs(title = expression(bold("Leaf N:P N+P response to ")*bolditalic("PAR")[bold("g")]),
       x = expression(bolditalic("PAR")[bold("g")]*bold(" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       y = "",
       size = expression(bold("Error"^"-1"))) +
  guides(size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_classic(base_size = 20) +
  theme(title = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(color = "black", size = 20))
npadd_leafnp_par_plot 

###################
# Merge leaf N:P results, plots
###################
# Merge leaf N:P moderator results, with some light cleaning
leafnp_clim_summary <- rbind(nadd_leafnp_clim_summary, 
                            padd_leafnp_clim_summary, 
                            npadd_leafnp_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

# Write plots
png("../plots/supp/CNP_figSX_leafnp_climate.png", height = 14, width = 14.5,
    units = "in", res = 600)
ggarrange(nadd_leafnp_tg_plot, nadd_leafnp_ai_plot, nadd_leafnp_par_plot,
          padd_leafnp_tg_plot, padd_leafnp_ai_plot, padd_leafnp_par_plot,
          npadd_leafnp_tg_plot, npadd_leafnp_ai_plot, npadd_leafnp_par_plot,
          ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"),
          font.label = list(size = 22))
dev.off()

##############################################################################
# Total biomass climate moderators
##############################################################################
###############
# N addition
###############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "tbio_gm2" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr > -0.75)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
nadd_tbio_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "n" & 
                                     myvar == "tbio_gm2" & 
                                     !is.na(gs_mat) & gs_ai < 3 & logr > -0.75))

# Climate summary
nadd_tbio_clim_summary <- data.frame(trait = "tbio_gm2",
                                     nut_add = "n",
                                     k = 29,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(nadd_tbio_clim)),
                                     row.names = NULL)

###############
# P addition
###############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "tbio_gm2" & 
                                        !is.na(gs_mat) & gs_ai < 3)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
padd_tbio_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "p" & 
                                    myvar == "tbio_gm2" & 
                                    !is.na(gs_mat) & gs_ai < 3))

# Climate summary
padd_tbio_clim_summary <- data.frame(trait = "tbio",
                                     nut_add = "p",
                                     k = 30,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(padd_tbio_clim)),
                                     row.names = NULL)

###############
# N+P addition
###############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "tbio_gm2" & 
                                        !is.na(gs_mat) & gs_ai < 3)) +
  geom_point(aes(x = gs_mat, y = logr))

# Model
npadd_tbio_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(nut_add == "np" & 
                                      myvar == "tbio_gm2" & 
                                      !is.na(gs_mat) & gs_ai < 3))

# Climate summary
npadd_tbio_clim_summary <- data.frame(trait = "tbio",
                                      nut_add = "np",
                                      k = 30,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(npadd_tbio_clim)),
                                      row.names = NULL)

###################
# Merge total biomass results, plots
###################
# Merge total biomass moderator results, with some light cleaning
tbio_clim_summary <- rbind(nadd_tbio_clim_summary, 
                            padd_tbio_clim_summary, 
                            npadd_tbio_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

# Create plot


##############################################################################
# Aboveground biomass (g/m2) climate moderators
##############################################################################

# N addition
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "anpp" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr < 1.5)) +
  geom_point(aes(x = gs_mat, y = logr))

nadd_anpp_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "n" & 
                                    myvar == "anpp" & 
                                    !is.na(gs_mat) & gs_ai < 3 & logr < 1.5))

nadd_agb_clim_summary <- data.frame(trait = "anpp",
                                    nut_add = "n",
                                    k = 111,
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(nadd_anpp_clim)),
                                    row.names = NULL)

# P addition
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "anpp" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr < 0.9 & logr > -0.8)) +
  geom_point(aes(x = gs_mat, y = logr))

padd_anpp_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "p" & 
                                    myvar == "anpp" & 
                                    !is.na(gs_mat) & gs_ai < 3 & logr < 0.9 & logr > -0.8))

padd_agb_clim_summary <- data.frame(trait = "anpp",
                                    nut_add = "p",
                                    k = 104,
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(padd_anpp_clim)),
                                    row.names = NULL)

# N+P addition
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "anpp" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr < 1.5 & logr > -0.5)) +
  geom_point(aes(x = gs_mat, y = logr))

npadd_anpp_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "np" & 
                                     myvar == "anpp" & 
                                     !is.na(gs_mat) & gs_ai < 3 & logr < 1.5 & logr > -0.5))

npadd_agb_clim_summary <- data.frame(trait = "anpp",
                                     nut_add = "np",
                                     k = 102,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(npadd_anpp_clim)),
                                     row.names = NULL)

## Merge Parea moderator results, with some light cleaning
agb_clim_summary <- rbind(nadd_agb_clim_summary, 
                           padd_agb_clim_summary, 
                           npadd_agb_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

##############################################################################
# Belowground biomass climate moderators
##############################################################################

# N addition
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "bnpp" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr > -0.9 & logr < 0.75)) +
  geom_point(aes(x = gs_mat, y = logr))

nadd_bnpp_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "n" & 
                                   myvar == "bnpp" & 
                                   !is.na(gs_mat) & gs_ai < 3 & logr > -0.9 & logr < 0.75))

nadd_bnpp_clim_summary <- data.frame(trait = "bnpp",
                                     nut_add = "n",
                                     k = 45,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(nadd_bnpp_clim)),
                                     row.names = NULL)

# P addition
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "bnpp" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr < 1)) +
  geom_point(aes(x = gs_mat, y = logr))

padd_bnpp_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "p" & 
                                   myvar == "bnpp" & 
                                   !is.na(gs_mat) & gs_ai < 3 & logr < 1))

padd_bnpp_clim_summary <- data.frame(trait = "bnpp",
                                     nut_add = "p",
                                     k = 51,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(padd_bnpp_clim)),
                                     row.names = NULL)

# N+P addition
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "bnpp" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr > -2 & logr < 1)) +
  geom_point(aes(x = gs_mat, y = logr))

npadd_bnpp_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "np" & 
                                    myvar == "bnpp" & 
                                    !is.na(gs_mat) & gs_ai < 3 & logr > -2 & logr < 1))

npadd_bnpp_clim_summary <- data.frame(trait = "bnpp",
                                      nut_add = "np",
                                      k = 48,
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(npadd_bnpp_clim)),
                                      row.names = NULL)

## Merge Parea moderator results, with some light cleaning
bgb_clim_summary <- rbind(nadd_bnpp_clim_summary, 
                          padd_bnpp_clim_summary, 
                          npadd_bnpp_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

##############################################################################
# Root mass fraction climate moderators
##############################################################################

# N addition
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "rmf" & 
                                        !is.na(gs_mat) & gs_ai < 3)) +
  geom_point(aes(x = gs_mat, y = logr))

nadd_rmf_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "n" & 
                                   myvar == "rmf" & 
                                   !is.na(gs_mat) & gs_ai < 3))

nadd_rmf_clim_summary <- data.frame(trait = "rmf",
                                    nut_add = "n",
                                    k = 34,
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(nadd_rmf_clim)),
                                    row.names = NULL)

# P addition
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "rmf" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr > -0.6)) +
  geom_point(aes(x = gs_mat, y = logr))

padd_rmf_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "p" & 
                                   myvar == "rmf" & 
                                   !is.na(gs_mat) & gs_ai < 3 & logr > -0.6))

padd_rmf_clim_summary <- data.frame(trait = "rmf",
                                    nut_add = "p",
                                    k = 32,
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(padd_rmf_clim)),
                                    row.names = NULL)

# N+P addition
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "rmf" & 
                                        !is.na(gs_mat) & gs_ai < 3)) +
  geom_point(aes(x = gs_mat, y = logr))

npadd_rmf_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "np" & 
                                    myvar == "rmf" & 
                                    !is.na(gs_mat) & gs_ai < 3))

npadd_rmf_clim_summary <- data.frame(trait = "rmf",
                                     nut_add = "np",
                                     k = 34,
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(npadd_rmf_clim)),
                                     row.names = NULL)

## Merge Parea moderator results, with some light cleaning
rmf_clim_summary <- rbind(nadd_rmf_clim_summary, 
                          padd_rmf_clim_summary, 
                          npadd_rmf_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

##############################################################################
# Root:shoot climate moderators
##############################################################################

# N addition
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "rootshoot" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr > -1.5)) +
  geom_point(aes(x = gs_mat, y = logr))

nadd_rootshoot_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "n" & 
                                   myvar == "rootshoot" & 
                                   !is.na(gs_mat) & gs_ai < 3 & logr > -1.5))

nadd_rootshoot_clim_summary <- data.frame(trait = "rootshoot",
                                          nut_add = "n",
                                          k = 34,
                                          mod = c("intrcpt", "gs_mat",
                                                  "gs_ai", "gs_par"),
                                          coef(summary(nadd_rootshoot_clim)),
                                          row.names = NULL)

# P addition
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "rootshoot" & 
                                        !is.na(gs_mat) & gs_ai < 3 & logr > -1.5 & logr < 1)) +
  geom_point(aes(x = gs_mat, y = logr))

padd_rootshoot_clim <- rma.mv(logr, 
                              logr_var,
                              method = "REML", 
                              random = ~ 1 | exp, 
                              mods = ~ gs_mat + gs_ai + gs_par,
                              slab = exp, 
                              control = list(stepadj = 0.3), 
                              data = meta_results %>% 
                                filter(nut_add == "p" & 
                                         myvar == "rootshoot" & 
                                         !is.na(gs_mat) & gs_ai < 3 & logr > -1.5 & logr < 1))

padd_rootshoot_clim_summary <- data.frame(trait = "rmf",
                                          nut_add = "p",
                                          k = 33,
                                          mod = c("intrcpt", "gs_mat",
                                                  "gs_ai", "gs_par"),
                                          coef(summary(padd_rootshoot_clim)),
                                          row.names = NULL)

# N+P addition
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "rootshoot" & 
                                        !is.na(gs_mat) & gs_ai < 3)) +
  geom_point(aes(x = gs_mat, y = logr))

npadd_rootshoot_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "np" & 
                                    myvar == "rootshoot" & 
                                    !is.na(gs_mat) & gs_ai < 3 & logr < 1))

npadd_rootshoot_clim_summary <- data.frame(trait = "rootshoot",
                                           nut_add = "np",
                                           k = 35,
                                           mod = c("intrcpt", "gs_mat",
                                                   "gs_ai", "gs_par"),
                                           coef(summary(npadd_rootshoot_clim)),
                                           row.names = NULL)

## Merge Parea moderator results, with some light cleaning
rootshoot_clim_summary <- rbind(nadd_rootshoot_clim_summary, 
                          padd_rootshoot_clim_summary, 
                          npadd_rootshoot_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 4)))

##############################################################################
# Merge climate moderators and write to .csv
##############################################################################
marea_clim_summary %>%
  full_join(nmass_clim_summary) %>%
  full_join(narea_clim_summary) %>%
  full_join(pmass_clim_summary) %>%
  full_join(parea_clim_summary) %>%
  full_join(leafnp_clim_summary) %>%
  full_join(tbio_clim_summary) %>%
  full_join(agb_clim_summary) %>%
  full_join(bgb_clim_summary) %>%
  full_join(rmf_clim_summary) %>%
  full_join(rootshoot_clim_summary) %>%
  mutate(across(estimate:ci.ub, \(x) round(x, digits = 3)),
         estimate_se = str_c(sprintf("%.3f", estimate), "±", sprintf("%.3f", se)),
         ci_range = str_c("[", sprintf("%.3f", ci.lb), ", ", sprintf("%.3f", ci.ub), "]")) %>%
  dplyr::select(trait:mod, estimate, se, estimate_se, zval, pval, ci.lb, ci.ub, ci_range) %>%
  write_excel_csv("../data/CNPmeta_clim_moderators.csv")


##############################################################################
# Marea plant functional type moderators
##############################################################################

#############
# N addition
#############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "lma" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model 
nadd_marea_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "n" & 
                                    myvar == "lma" & !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_marea_photo <- data.frame(trait = "marea", 
                               nut_add = "n",
                               mod = "photo",
                               mod_results(nadd_marea_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_marea_pft))[2,3],
                               p = coef(summary(nadd_marea_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_marea_myc <- data.frame(trait = "marea", 
                             nut_add = "n",
                             mod = "myc_nas",
                             mod_results(nadd_marea_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_marea_pft))[3,3],
                             p = coef(summary(nadd_marea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_marea_nfix <- data.frame(trait = "marea", 
                              nut_add = "n",
                              mod = "nfix",
                              mod_results(nadd_marea_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_marea_pft))[4,3],
                              p = coef(summary(nadd_marea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_marea_pft_results <- nadd_marea_photo %>% 
  rbind(nadd_marea_myc) %>% 
  rbind(nadd_marea_nfix) %>%
  mutate(k = 110)

#############
# P addition
#############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "lma" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
padd_marea_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "p" & 
                                     myvar == "lma" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_marea_photo <- data.frame(trait = "marea", 
                               nut_add = "p",
                               mod = "photo",
                               mod_results(padd_marea_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_marea_pft))[2,3],
                               p = coef(summary(padd_marea_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_marea_myc <- data.frame(trait = "marea", 
                             nut_add = "p",
                             mod = "myc_nas",
                             mod_results(padd_marea_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_marea_pft))[3,3],
                             p = coef(summary(padd_marea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_marea_nfix <- data.frame(trait = "marea", 
                              nut_add = "p",
                              mod = "nfix",
                              mod_results(padd_marea_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_marea_pft))[4,3],
                              p = coef(summary(padd_marea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_marea_pft_results <- padd_marea_photo %>% 
  rbind(padd_marea_myc) %>% 
  rbind(padd_marea_nfix) %>%
  mutate(k = 111)

#############
# N+P addition
#############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "lma" & 
                                        !is.na(photo_path) & logr > -1 & logr < 1)) +
  geom_point(aes(x = photo_path, y = logr))

# Model
npadd_marea_pft <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path + myc_nas + n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(nut_add == "np" & 
                                      myvar == "lma" & 
                                      !is.na(photo_path) & logr > -1 & logr < 1))

# Extract photosynthetic pathway summary statistics
npadd_marea_photo <- data.frame(trait = "marea", 
                               nut_add = "np",
                               mod = "photo",
                               mod_results(npadd_marea_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_marea_pft))[2,3],
                               p = coef(summary(npadd_marea_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_marea_myc <- data.frame(trait = "marea", 
                             nut_add = "np",
                             mod = "myc_nas",
                             mod_results(npadd_marea_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_marea_pft))[3,3],
                             p = coef(summary(npadd_marea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_marea_nfix <- data.frame(trait = "marea", 
                              nut_add = "np",
                              mod = "nfix",
                              mod_results(npadd_marea_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_marea_pft))[4,3],
                              p = coef(summary(npadd_marea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_marea_pft_results <- npadd_marea_photo %>% 
  rbind(npadd_marea_myc) %>% 
  rbind(npadd_marea_nfix) %>%
  mutate(k = 111)

#############
# Merge Marea moderator results, with some light cleaning
#############
marea_pft_summary <- rbind(nadd_marea_pft_results, 
                             padd_marea_pft_results, 
                             npadd_marea_pft_results) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)),
         nut_add = factor(nut_add, levels = c("n", "p", "np"))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

#############
# Make Marea moderator plots
#############

# Photosynthetic pathway plot
marea_photo_plot <- ggplot(data = subset(marea_pft_summary, mod == "photo"), 
                           aes(x = nut_add, y = estimate, fill = nut_add, shape = comp)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL),
                position = position_dodge(width = 0.75),
                width = 0.2, size = 1) +
  geom_point(size = 5, position = position_dodge(width = 0.75)) +
  scale_shape_manual(values = c(21, 24),
                     labels = c(expression("C"["3"]),
                                expression("C"["4"]))) +
  scale_fill_manual(values = c("red", "blue", "magenta")) +
  scale_x_discrete(labels = c("N", "P", "N+P")) +
  scale_y_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, 0.15)) +
  labs(title = expression(bolditalic("M")[bold("area")]*bold(" by photo. pathway")),
       x = "",
       y = "",
       shape = expression(bold("Photosynthetic\npathway"))) +
  guides(fill = "none") +
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
marea_photo_plot

# N-fixer plot
marea_nfix_plot <- ggplot(data = subset(marea_pft_summary, mod == "nfix"), 
                           aes(x = nut_add, y = estimate, fill = nut_add, shape = comp)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL),
                position = position_dodge(width = 0.75),
                width = 0.2, size = 1) +
  geom_point(size = 5, position = position_dodge(width = 0.75)) +
  scale_shape_manual(values = c(21, 24),
                     labels = c(expression("N"["2"]*"-fixer"),
                                "non-fixer")) +
  scale_fill_manual(values = c("red", "blue", "magenta")) +
  scale_x_discrete(labels = c("N", "P", "N+P")) +
  scale_y_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, 0.15)) +
  labs(title = expression(bolditalic("M")[bold("area")]*bold(" by N"["2"]*" fixation ability")),
       x = "",
       y = "Log-response ratio",
       shape = expression(bold("N"["2"]*" fixation ability"))) +
  guides(fill = "none") +
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
marea_nfix_plot

# Mycorrhizal acquisition strategy plot
marea_myc_plot <- ggplot(data = subset(marea_pft_summary, mod == "myc_nas"), 
                          aes(x = nut_add, y = estimate, fill = nut_add, shape = comp)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL),
                position = position_dodge(width = 0.75),
                width = 0.2, size = 1) +
  geom_point(size = 5, position = position_dodge(width = 0.75)) +
  scale_shape_manual(values = c(21, 24),
                     labels = c("mining",
                                "scavenging")) +
  scale_fill_manual(values = c("red", "blue", "magenta")) +
  scale_x_discrete(labels = c("N", "P", "N+P")) +
  scale_y_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, 0.15)) +
  labs(title = expression(bolditalic("M")[bold("area")]*bold(" by mycorrhizal strategy")),
       x = "Nutrient addition",
       y = "",
       shape = expression(bold("Mycorrhizal\nacquisition\nstrategy"))) +
  guides(fill = "none") +
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
marea_myc_plot

# Make plot
png("../plots/supp/CNP_figSX_marea_pft.png", height = 14, width = 8, units = "in",
    res = 600)
ggarrange(marea_photo_plot, marea_nfix_plot, marea_myc_plot,
          nrow = 3, ncol = 1, align = "hv",
          labels = c("(a)", "(b)", "(c)"), font.label = list(size = 20))
dev.off()

##############################################################################
# Nmass
##############################################################################

#############
# N addition
#############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "leaf_n_mass" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
nadd_nmass_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "n" & 
                                    myvar == "leaf_n_mass" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_nmass_photo <- data.frame(trait = "nmass", 
                               nut_add = "n",
                               mod = "photo",
                               mod_results(nadd_nmass_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_nmass_pft))[2,3],
                               p = coef(summary(nadd_nmass_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_nmass_myc <- data.frame(trait = "nmass", 
                              nut_add = "n",
                              mod = "myc_nas",
                              mod_results(nadd_nmass_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_nmass_pft))[3,3],
                              p = coef(summary(nadd_nmass_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_nmass_nfix <- data.frame(trait = "nmass", 
                               nut_add = "n",
                               mod = "nfix",
                               mod_results(nadd_nmass_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_nmass_pft))[4,3],
                               p = coef(summary(nadd_nmass_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_nmass_pft_results <- nadd_nmass_photo %>% 
  rbind(nadd_nmass_myc) %>% 
  rbind(nadd_nmass_nfix) %>%
  mutate(k = 136)

##############
# P addition
##############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "leaf_n_mass" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
padd_nmass_pft <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path + myc_nas + n_fixer,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(nut_add == "p" & 
                                      myvar == "leaf_n_mass" & 
                                      !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_nmass_photo <- data.frame(trait = "nmass", 
                               nut_add = "p",
                               mod = "photo",
                               mod_results(padd_nmass_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_nmass_pft))[2,3],
                               p = coef(summary(padd_nmass_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_nmass_myc <- data.frame(trait = "nmass", 
                             nut_add = "p",
                             mod = "myc_nas",
                             mod_results(padd_nmass_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_nmass_pft))[3,3],
                             p = coef(summary(padd_nmass_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_nmass_nfix <- data.frame(trait = "nmass", 
                              nut_add = "p",
                              mod = "nfix",
                              mod_results(padd_nmass_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_nmass_pft))[4,3],
                              p = coef(summary(padd_nmass_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_nmass_pft_results <- padd_nmass_photo %>% 
  rbind(padd_nmass_myc) %>% 
  rbind(padd_nmass_nfix) %>%
  mutate(k = 136)

################
# N+P addition
################
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "leaf_n_mass" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
npadd_nmass_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "np" & 
                                     myvar == "leaf_n_mass" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_nmass_photo <- data.frame(trait = "nmass", 
                               nut_add = "np",
                               mod = "photo",
                               mod_results(npadd_nmass_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_nmass_pft))[2,3],
                               p = coef(summary(npadd_nmass_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_nmass_myc <- data.frame(trait = "nmass", 
                             nut_add = "np",
                             mod = "myc_nas",
                             mod_results(npadd_nmass_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_nmass_pft))[3,3],
                             p = coef(summary(npadd_nmass_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_nmass_nfix <- data.frame(trait = "nmass", 
                              nut_add = "np",
                              mod = "nfix",
                              mod_results(npadd_nmass_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_nmass_pft))[4,3],
                              p = coef(summary(npadd_nmass_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_nmass_pft_results <- npadd_nmass_photo %>% 
  rbind(npadd_nmass_myc) %>% 
  rbind(npadd_nmass_nfix) %>%
  mutate(k = 136)

#############
# Merge Nmass moderator results, with some light cleaning
#############
nmass_pft_summary <- rbind(nadd_nmass_pft_results, 
                             padd_nmass_pft_results, 
                             npadd_nmass_pft_results) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)),
         nut_add = factor(nut_add, levels = c("n", "p", "np"))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

#############
# Make Nmass moderator plots
#############

# Photosynthetic pathway plot
nmass_photo_plot <- ggplot(data = subset(nmass_pft_summary, mod == "photo"), 
                           aes(x = nut_add, y = estimate, fill = nut_add, shape = comp)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL),
                position = position_dodge(width = 0.75),
                width = 0.2, size = 1) +
  geom_point(size = 5, position = position_dodge(width = 0.75)) +
  scale_shape_manual(values = c(21, 24),
                     labels = c(expression("C"["3"]),
                                expression("C"["4"]))) +
  scale_fill_manual(values = c("red", "blue", "magenta")) +
  scale_x_discrete(labels = c("N", "P", "N+P")) +
  scale_y_continuous(limits = c(-0.2, 0.4), breaks = seq(-0.2, 0.4, 0.2)) +
  labs(title = expression(bolditalic("N")[bold("mass")]*bold(" by photo. pathway")),
       x = "",
       y = "",
       shape = expression(bold("Photosynthetic\npathway"))) +
  guides(fill = "none") +
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
nmass_photo_plot

# N-fixer plot
nmass_nfix_plot <- ggplot(data = subset(nmass_pft_summary, mod == "nfix"), 
                          aes(x = nut_add, y = estimate, fill = nut_add, shape = comp)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL),
                position = position_dodge(width = 0.75),
                width = 0.2, size = 1) +
  geom_point(size = 5, position = position_dodge(width = 0.75)) +
  scale_shape_manual(values = c(21, 24),
                     labels = c(expression("N"["2"]*"-fixer"),
                                "non-fixer")) +
  scale_fill_manual(values = c("red", "blue", "magenta")) +
  scale_x_discrete(labels = c("N", "P", "N+P")) +
  scale_y_continuous(limits = c(-0.2, 0.4), breaks = seq(-0.2, 0.4, 0.2)) +
  labs(title = expression(bolditalic("N")[bold("mass")]*bold(" by N"["2"]*" fixation ability")),
       x = "",
       y = "Log-response ratio",
       shape = expression(bold("N"["2"]*" fixation ability"))) +
  guides(fill = "none") +
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
nmass_nfix_plot

# Mycorrhizal acquisition strategy plot
nmass_myc_plot <- ggplot(data = subset(nmass_pft_summary, mod == "myc_nas"), 
                         aes(x = nut_add, y = estimate, fill = nut_add, shape = comp)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_errorbar(aes(ymin = lowerCL, ymax = upperCL),
                position = position_dodge(width = 0.75),
                width = 0.2, size = 1) +
  geom_point(size = 5, position = position_dodge(width = 0.75)) +
  scale_shape_manual(values = c(21, 24),
                     labels = c("mining",
                                "scavenging")) +
  scale_fill_manual(values = c("red", "blue", "magenta")) +
  scale_x_discrete(labels = c("N", "P", "N+P")) +
  scale_y_continuous(limits = c(-0.2, 0.4), breaks = seq(-0.2, 0.4, 0.2)) +
  labs(title = expression(bolditalic("N")[bold("mass")]*bold(" by mycorrhizal strategy")),
       x = "Nutrient addition",
       y = "",
       shape = expression(bold("Mycorrhizal\nacquisition\nstrategy"))) +
  guides(fill = "none") +
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        title = element_text(face = "bold"))
nmass_myc_plot

# Make plot
png("../plots/supp/CNP_figSX_nmass_pft.png", height = 14, width = 14, units = "in",
    res = 600)
ggarrange(nmass_photo_plot, nmass_nfix_plot, nmass_myc_plot,
          nrow = 3, ncol = 1, align = "hv",
          labels = c("(a)", "(b)", "(c)"), font.label = list(size = 20))
dev.off()

##############################################################################
# Narea - pft
##############################################################################

##############
# N addition
##############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "leaf_n_area" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
nadd_narea_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "n" & 
                                    myvar == "leaf_n_area" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_narea_photo <- data.frame(trait = "narea", 
                             nut_add = "n",
                             mod = "photo",
                             mod_results(nadd_narea_pft, 
                                         mod = "photo_path", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_narea_pft))[2,3],
                             p = coef(summary(nadd_narea_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_narea_myc <- data.frame(trait = "narea", 
                              nut_add = "n",
                              mod = "myc_nas",
                              mod_results(nadd_narea_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_narea_pft))[3,3],
                              p = coef(summary(nadd_narea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_narea_nfix <- data.frame(trait = "narea", 
                               nut_add = "n",
                               mod = "nfix",
                               mod_results(nadd_narea_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_narea_pft))[4,3],
                               p = coef(summary(nadd_narea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_narea_pft_results <- nadd_narea_photo %>% 
  rbind(nadd_narea_myc) %>% 
  rbind(nadd_narea_nfix) %>%
  mutate(k = 87)

##############
# P addition
##############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "leaf_n_area" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
padd_narea_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "p" & 
                                    myvar == "leaf_n_area" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_narea_photo <- data.frame(trait = "narea", 
                               nut_add = "p",
                               mod = "photo",
                               mod_results(padd_narea_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_narea_pft))[2,3],
                               p = coef(summary(padd_narea_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_narea_myc <- data.frame(trait = "narea", 
                              nut_add = "p",
                              mod = "myc_nas",
                              mod_results(padd_narea_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_narea_pft))[3,3],
                              p = coef(summary(padd_narea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_narea_nfix <- data.frame(trait = "narea", 
                               nut_add = "p",
                               mod = "nfix",
                               mod_results(padd_narea_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_narea_pft))[4,3],
                               p = coef(summary(padd_narea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_narea_pft_results <- padd_narea_photo %>% 
  rbind(padd_narea_myc) %>% 
  rbind(padd_narea_nfix) %>%
  mutate(k = 86)

################
# N+P addition
################

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "leaf_n_area" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
npadd_narea_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "np" & 
                                     myvar == "leaf_n_area" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_narea_photo <- data.frame(trait = "narea", 
                                nut_add = "np",
                                mod = "photo",
                                mod_results(npadd_narea_pft, 
                                            mod = "photo_path", 
                                            group = "exp")$mod_table,
                                z = coef(summary(npadd_narea_pft))[2,3],
                                p = coef(summary(npadd_narea_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_narea_myc <- data.frame(trait = "narea", 
                              nut_add = "np",
                              mod = "myc_nas",
                              mod_results(npadd_narea_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_narea_pft))[3,3],
                              p = coef(summary(npadd_narea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_narea_nfix <- data.frame(trait = "narea", 
                               nut_add = "np",
                               mod = "nfix",
                               mod_results(npadd_narea_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_narea_pft))[4,3],
                               p = coef(summary(npadd_narea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_narea_pft_results <- npadd_narea_photo %>% 
  rbind(npadd_narea_myc) %>% 
  rbind(npadd_narea_nfix) %>%
  mutate(k = 87)

#############
# Merge Narea moderator results, with some light cleaning
#############
narea_pft_summary <- rbind(nadd_narea_pft_results, 
                           padd_narea_pft_results, 
                           npadd_narea_pft_results) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

#############
# Make Narea moderator plots
#############

##############################################################################
# Pmass - plant functional type
##############################################################################

##############
# N addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "leaf_p_mass" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
nadd_pmass_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "n" & 
                                    myvar == "leaf_p_mass" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_pmass_photo <- data.frame(trait = "pmass", 
                               nut_add = "n",
                               mod = "photo",
                               mod_results(nadd_pmass_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_pmass_pft))[2,3],
                               p = coef(summary(nadd_pmass_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_pmass_myc <- data.frame(trait = "pmass", 
                             nut_add = "n",
                             mod = "myc_nas",
                             mod_results(nadd_pmass_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_pmass_pft))[3,3],
                             p = coef(summary(nadd_pmass_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_pmass_nfix <- data.frame(trait = "pmass", 
                              nut_add = "n",
                              mod = "nfix",
                              mod_results(nadd_pmass_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_pmass_pft))[4,3],
                              p = coef(summary(nadd_pmass_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_pmass_pft_results <- nadd_pmass_photo %>% 
  rbind(nadd_pmass_myc) %>% 
  rbind(nadd_pmass_nfix) %>%
  mutate(k = 129)

##############
# P addition
##############
# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "leaf_p_mass" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
padd_pmass_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "p" & 
                                    myvar == "leaf_p_mass" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_pmass_photo <- data.frame(trait = "pmass", 
                               nut_add = "p",
                               mod = "photo",
                               mod_results(padd_pmass_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_pmass_pft))[2,3],
                               p = coef(summary(padd_pmass_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_pmass_myc <- data.frame(trait = "pmass", 
                             nut_add = "p",
                             mod = "myc_nas",
                             mod_results(padd_pmass_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_pmass_pft))[3,3],
                             p = coef(summary(padd_pmass_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_pmass_nfix <- data.frame(trait = "pmass", 
                              nut_add = "p",
                              mod = "nfix",
                              mod_results(padd_pmass_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_pmass_pft))[4,3],
                              p = coef(summary(padd_pmass_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_pmass_pft_results <- padd_pmass_photo %>% 
  rbind(padd_pmass_myc) %>% 
  rbind(padd_pmass_nfix) %>%
  mutate(k = 128)

################
# N+P addition
################

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "leaf_p_mass" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
npadd_pmass_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "np" & 
                                     myvar == "leaf_p_mass" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_pmass_photo <- data.frame(trait = "pmass", 
                                nut_add = "np",
                                mod = "photo",
                                mod_results(npadd_pmass_pft, 
                                            mod = "photo_path", 
                                            group = "exp")$mod_table,
                                z = coef(summary(npadd_pmass_pft))[2,3],
                                p = coef(summary(npadd_pmass_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_pmass_myc <- data.frame(trait = "pmass", 
                              nut_add = "np",
                              mod = "myc_nas",
                              mod_results(npadd_pmass_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_pmass_pft))[3,3],
                              p = coef(summary(npadd_pmass_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_pmass_nfix <- data.frame(trait = "pmass", 
                               nut_add = "np",
                               mod = "nfix",
                               mod_results(npadd_pmass_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_pmass_pft))[4,3],
                               p = coef(summary(npadd_pmass_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_pmass_pft_results <- npadd_pmass_photo %>% 
  rbind(npadd_pmass_myc) %>% 
  rbind(npadd_pmass_nfix) %>%
  mutate(k = 130)

#############
# Merge Pmass moderator results, with some light cleaning
#############
pmass_pft_summary <- rbind(nadd_pmass_pft_results, 
                           padd_pmass_pft_results, 
                           npadd_pmass_pft_results) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Parea - plant functional type
##############################################################################

##############
# N addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "leaf_p_area" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
nadd_parea_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "n" & 
                                    myvar == "leaf_p_area" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_parea_photo <- data.frame(trait = "parea", 
                               nut_add = "n",
                               mod = "photo",
                               mod_results(nadd_parea_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_parea_pft))[2,3],
                               p = coef(summary(nadd_parea_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_parea_myc <- data.frame(trait = "parea", 
                             nut_add = "n",
                             mod = "myc_nas",
                             mod_results(nadd_parea_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_parea_pft))[3,3],
                             p = coef(summary(nadd_parea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_parea_nfix <- data.frame(trait = "parea", 
                              nut_add = "n",
                              mod = "nfix",
                              mod_results(nadd_parea_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_parea_pft))[4,3],
                              p = coef(summary(nadd_parea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_parea_pft_results <- nadd_parea_photo %>% 
  rbind(nadd_parea_myc) %>% 
  rbind(nadd_parea_nfix) %>%
  mutate(k = 82)

##############
# P addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "leaf_p_area" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
padd_parea_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "p" & 
                                    myvar == "leaf_p_area" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_parea_photo <- data.frame(trait = "parea", 
                               nut_add = "p",
                               mod = "photo",
                               mod_results(padd_parea_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_parea_pft))[2,3],
                               p = coef(summary(padd_parea_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_parea_myc <- data.frame(trait = "parea", 
                             nut_add = "p",
                             mod = "myc_nas",
                             mod_results(padd_parea_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_parea_pft))[3,3],
                             p = coef(summary(padd_parea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_parea_nfix <- data.frame(trait = "parea", 
                              nut_add = "p",
                              mod = "nfix",
                              mod_results(padd_parea_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_parea_pft))[4,3],
                              p = coef(summary(padd_parea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_parea_pft_results <- padd_parea_photo %>% 
  rbind(padd_parea_myc) %>% 
  rbind(padd_parea_nfix) %>%
  mutate(k = 82)

################
# N+P addition
################

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "leaf_p_area" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
npadd_parea_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "np" & 
                                     myvar == "leaf_p_area" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_parea_photo <- data.frame(trait = "parea", 
                                nut_add = "np",
                                mod = "photo",
                                mod_results(npadd_parea_pft, 
                                            mod = "photo_path", 
                                            group = "exp")$mod_table,
                                z = coef(summary(npadd_parea_pft))[2,3],
                                p = coef(summary(npadd_parea_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_parea_myc <- data.frame(trait = "parea", 
                              nut_add = "np",
                              mod = "myc_nas",
                              mod_results(npadd_parea_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_parea_pft))[3,3],
                              p = coef(summary(npadd_parea_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_parea_nfix <- data.frame(trait = "parea", 
                               nut_add = "np",
                               mod = "nfix",
                               mod_results(npadd_parea_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_parea_pft))[4,3],
                               p = coef(summary(npadd_parea_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_parea_pft_results <- npadd_parea_photo %>% 
  rbind(npadd_parea_myc) %>% 
  rbind(npadd_parea_nfix) %>%
  mutate(k = 82)

#############
# Merge Parea moderator results, with some light cleaning
#############
parea_pft_summary <- rbind(nadd_parea_pft_results, 
                           padd_parea_pft_results, 
                           npadd_parea_pft_results) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Leaf N:P - plant functional type
##############################################################################

##############
# N addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "leaf_np" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
nadd_leafnp_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "n" & 
                                    myvar == "leaf_np" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_leafnp_photo <- data.frame(trait = "leaf_np", 
                               nut_add = "n",
                               mod = "photo",
                               mod_results(nadd_leafnp_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_leafnp_pft))[2,3],
                               p = coef(summary(nadd_leafnp_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_leafnp_myc <- data.frame(trait = "leaf_np", 
                             nut_add = "n",
                             mod = "myc_nas",
                             mod_results(nadd_leafnp_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_leafnp_pft))[3,3],
                             p = coef(summary(nadd_leafnp_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_leafnp_nfix <- data.frame(trait = "leaf_np", 
                              nut_add = "n",
                              mod = "nfix",
                              mod_results(nadd_leafnp_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_leafnp_pft))[4,3],
                              p = coef(summary(nadd_leafnp_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_leafnp_pft_results <- nadd_leafnp_photo %>% 
  rbind(nadd_leafnp_myc) %>% 
  rbind(nadd_leafnp_nfix) %>%
  mutate(k = 115)

##############
# P addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "leaf_np" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
padd_leafnp_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "p" & 
                                    myvar == "leaf_np" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_leafnp_photo <- data.frame(trait = "leaf_np", 
                               nut_add = "p",
                               mod = "photo",
                               mod_results(padd_leafnp_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_leafnp_pft))[2,3],
                               p = coef(summary(padd_leafnp_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_leafnp_myc <- data.frame(trait = "leaf_np", 
                             nut_add = "p",
                             mod = "myc_nas",
                             mod_results(padd_leafnp_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_leafnp_pft))[3,3],
                             p = coef(summary(padd_leafnp_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_leafnp_nfix <- data.frame(trait = "leaf_np", 
                              nut_add = "p",
                              mod = "nfix",
                              mod_results(padd_leafnp_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_leafnp_pft))[4,3],
                              p = coef(summary(padd_leafnp_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_leafnp_pft_results <- padd_leafnp_photo %>% 
  rbind(padd_leafnp_myc) %>% 
  rbind(padd_leafnp_nfix) %>%
  mutate(k = 115)

################
# N+P addition
################

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "leaf_np" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
npadd_leafnp_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "np" & 
                                     myvar == "leaf_np" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_leafnp_photo <- data.frame(trait = "leaf_np", 
                                 nut_add = "np",
                                 mod = "photo",
                                 mod_results(npadd_leafnp_pft, 
                                             mod = "photo_path", 
                                             group = "exp")$mod_table,
                                 z = coef(summary(npadd_leafnp_pft))[2,3],
                                 p = coef(summary(npadd_leafnp_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_leafnp_myc <- data.frame(trait = "leaf_np", 
                               nut_add = "np",
                               mod = "myc_nas",
                               mod_results(npadd_leafnp_pft, 
                                           mod = "myc_nas", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_leafnp_pft))[3,3],
                               p = coef(summary(npadd_leafnp_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_leafnp_nfix <- data.frame(trait = "leaf_np", 
                                nut_add = "np",
                                mod = "nfix",
                                mod_results(npadd_leafnp_pft, 
                                            mod = "n_fixer", 
                                            group = "exp")$mod_table,
                                z = coef(summary(npadd_leafnp_pft))[4,3],
                                p = coef(summary(npadd_leafnp_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_leafnp_pft_results <- npadd_leafnp_photo %>% 
  rbind(npadd_leafnp_myc) %>% 
  rbind(npadd_leafnp_nfix) %>%
  mutate(k = 115)

#############
# Merge leaf N:P moderator results, with some light cleaning
#############
leafnp_pft_summary <- rbind(nadd_leafnp_pft_results, 
                            padd_leafnp_pft_results, 
                            npadd_leafnp_pft_results) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Asat - plant functional type
##############################################################################

##############
# N addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "asat" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
nadd_asat_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "n" & 
                                   myvar == "asat" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_asat_photo <- data.frame(trait = "asat", 
                              nut_add = "n",
                              mod = "photo",
                              mod_results(nadd_asat_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_asat_pft))[2,3],
                              p = coef(summary(nadd_asat_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_asat_myc <- data.frame(trait = "asat", 
                            nut_add = "n",
                            mod = "myc_nas",
                            mod_results(nadd_asat_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(nadd_asat_pft))[3,3],
                            p = coef(summary(nadd_asat_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_asat_nfix <- data.frame(trait = "asat", 
                             nut_add = "n",
                             mod = "nfix",
                             mod_results(nadd_asat_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_asat_pft))[4,3],
                             p = coef(summary(nadd_asat_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_asat_pft_results <- nadd_asat_photo %>% 
  rbind(nadd_asat_myc) %>% 
  rbind(nadd_asat_nfix) %>%
  mutate(k = 93)

##############
# P addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "asat" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
padd_asat_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "p" & 
                                   myvar == "asat" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_asat_photo <- data.frame(trait = "asat", 
                              nut_add = "p",
                              mod = "photo",
                              mod_results(padd_asat_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_asat_pft))[2,3],
                              p = coef(summary(padd_asat_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_asat_myc <- data.frame(trait = "asat", 
                            nut_add = "p",
                            mod = "myc_nas",
                            mod_results(padd_asat_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(padd_asat_pft))[3,3],
                            p = coef(summary(padd_asat_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_asat_nfix <- data.frame(trait = "asat", 
                             nut_add = "p",
                             mod = "nfix",
                             mod_results(padd_asat_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_asat_pft))[4,3],
                             p = coef(summary(padd_asat_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_asat_pft_results <- padd_asat_photo %>% 
  rbind(padd_asat_myc) %>% 
  rbind(padd_asat_nfix) %>%
  mutate(k = 93)

################
# N+P addition
################

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "asat" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
npadd_asat_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "np" & 
                                    myvar == "asat" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_asat_photo <- data.frame(trait = "asat", 
                               nut_add = "np",
                               mod = "photo",
                               mod_results(npadd_asat_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_asat_pft))[2,3],
                               p = coef(summary(npadd_asat_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_asat_myc <- data.frame(trait = "asat", 
                             nut_add = "np",
                             mod = "myc_nas",
                             mod_results(npadd_asat_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_asat_pft))[3,3],
                             p = coef(summary(npadd_asat_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_asat_nfix <- data.frame(trait = "asat", 
                              nut_add = "np",
                              mod = "nfix",
                              mod_results(npadd_asat_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_asat_pft))[4,3],
                              p = coef(summary(npadd_asat_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_asat_pft_results <- npadd_asat_photo %>% 
  rbind(npadd_asat_myc) %>% 
  rbind(npadd_asat_nfix) %>%
  mutate(k = 92)

#############
# Merge Asat moderator results, with some light cleaning
#############
asat_pft_summary <- rbind(nadd_asat_pft_results, 
                             padd_asat_pft_results, 
                             npadd_asat_pft_results) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# gsw - plant functional type
##############################################################################

##############
# N addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "gsw" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
nadd_gsw_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "n" & 
                                    myvar == "gsw" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_gsw_photo <- data.frame(trait = "gsw", 
                               nut_add = "n",
                               mod = "photo",
                               mod_results(nadd_gsw_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_gsw_pft))[2,3],
                               p = coef(summary(nadd_gsw_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_gsw_myc <- data.frame(trait = "gsw", 
                             nut_add = "n",
                             mod = "myc_nas",
                             mod_results(nadd_gsw_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_gsw_pft))[3,3],
                             p = coef(summary(nadd_gsw_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_gsw_nfix <- data.frame(trait = "gsw", 
                              nut_add = "n",
                              mod = "nfix",
                              mod_results(nadd_gsw_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_gsw_pft))[4,3],
                              p = coef(summary(nadd_gsw_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_gsw_pft_results <- nadd_gsw_photo %>% 
  rbind(nadd_gsw_myc) %>% 
  rbind(nadd_gsw_nfix) %>%
  mutate(k = 47)

##############
# P addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "gsw" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
padd_gsw_pft <- rma.mv(logr, 
                       logr_var,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ photo_path + myc_nas + n_fixer,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results %>% 
                         filter(nut_add == "p" & 
                                  myvar == "gsw" & 
                                  !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_gsw_photo <- data.frame(trait = "gsw", 
                               nut_add = "p",
                               mod = "photo",
                               mod_results(padd_gsw_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_gsw_pft))[2,3],
                               p = coef(summary(padd_gsw_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_gsw_myc <- data.frame(trait = "gsw", 
                             nut_add = "p",
                             mod = "myc_nas",
                             mod_results(padd_gsw_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_gsw_pft))[3,3],
                             p = coef(summary(padd_gsw_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_gsw_nfix <- data.frame(trait = "gsw", 
                              nut_add = "p",
                              mod = "nfix",
                              mod_results(padd_gsw_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_gsw_pft))[4,3],
                              p = coef(summary(padd_gsw_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_gsw_pft_results <- padd_gsw_photo %>% 
  rbind(padd_gsw_myc) %>% 
  rbind(padd_gsw_nfix) %>%
  mutate(k = 47)

################
# N+P addition
################

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "gsw" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
npadd_gsw_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "np" & 
                                   myvar == "gsw" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_gsw_photo <- data.frame(trait = "gsw", 
                                nut_add = "np",
                                mod = "photo",
                                mod_results(npadd_gsw_pft, 
                                            mod = "photo_path", 
                                            group = "exp")$mod_table,
                                z = coef(summary(npadd_gsw_pft))[2,3],
                                p = coef(summary(npadd_gsw_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_gsw_myc <- data.frame(trait = "gsw", 
                              nut_add = "np",
                              mod = "myc_nas",
                              mod_results(npadd_gsw_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_gsw_pft))[3,3],
                              p = coef(summary(npadd_gsw_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_gsw_nfix <- data.frame(trait = "gsw", 
                               nut_add = "np",
                               mod = "nfix",
                               mod_results(npadd_gsw_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_gsw_pft))[4,3],
                               p = coef(summary(npadd_gsw_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_gsw_pft_results <- npadd_gsw_photo %>% 
  rbind(npadd_gsw_myc) %>% 
  rbind(npadd_gsw_nfix) %>%
  mutate(k = 47)

#############
# Merge gsw moderator results, with some light cleaning
#############
gsw_pft_summary <- rbind(nadd_gsw_pft_results, 
                           padd_gsw_pft_results, 
                           npadd_gsw_pft_results) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Rd - plant functional type
##############################################################################

##############
# N addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "rd" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
nadd_rd_pft <- rma.mv(logr, 
                       logr_var,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ photo_path + myc_nas + n_fixer,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results %>% 
                         filter(nut_add == "n" & 
                                  myvar == "rd" & 
                                  !is.na(photo_path) & logr < 1))

# Extract photosynthetic pathway summary statistics
nadd_rd_photo <- data.frame(trait = "rd", 
                             nut_add = "n",
                             mod = "photo",
                             mod_results(nadd_rd_pft, 
                                         mod = "photo_path", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_rd_pft))[2,3],
                             p = coef(summary(nadd_rd_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_rd_myc <- data.frame(trait = "rd", 
                           nut_add = "n",
                           mod = "myc_nas",
                           mod_results(nadd_rd_pft, 
                                       mod = "myc_nas", 
                                       group = "exp")$mod_table,
                           z = coef(summary(nadd_rd_pft))[3,3],
                           p = coef(summary(nadd_rd_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_rd_nfix <- data.frame(trait = "rd", 
                            nut_add = "n",
                            mod = "nfix",
                            mod_results(nadd_rd_pft, 
                                        mod = "n_fixer", 
                                        group = "exp")$mod_table,
                            z = coef(summary(nadd_rd_pft))[4,3],
                            p = coef(summary(nadd_rd_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_rd_pft_results <- nadd_rd_photo %>% 
  rbind(nadd_rd_myc) %>% 
  rbind(nadd_rd_nfix) %>%
  mutate(k = 31)
  

##############
# P addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "rd" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
padd_rd_pft <- rma.mv(logr, 
                       logr_var,
                       method = "REML", 
                       random = ~ 1 | exp, 
                       mods = ~ photo_path + myc_nas + n_fixer,
                       slab = exp, 
                       control = list(stepadj = 0.3), 
                       data = meta_results %>% 
                         filter(nut_add == "p" & 
                                  myvar == "rd" & 
                                  !is.na(photo_path) & logr < 1 & logr > -1))

# Extract photosynthetic pathway summary statistics
padd_rd_photo <- data.frame(trait = "rd", 
                             nut_add = "p",
                             mod = "photo",
                             mod_results(padd_rd_pft, 
                                         mod = "photo_path", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_rd_pft))[2,3],
                             p = coef(summary(padd_rd_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_rd_myc <- data.frame(trait = "rd", 
                           nut_add = "p",
                           mod = "myc_nas",
                           mod_results(padd_rd_pft, 
                                       mod = "myc_nas", 
                                       group = "exp")$mod_table,
                           z = coef(summary(padd_rd_pft))[3,3],
                           p = coef(summary(padd_rd_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_rd_nfix <- data.frame(trait = "rd", 
                            nut_add = "p",
                            mod = "nfix",
                            mod_results(padd_rd_pft, 
                                        mod = "n_fixer", 
                                        group = "exp")$mod_table,
                            z = coef(summary(padd_rd_pft))[4,3],
                            p = coef(summary(padd_rd_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_rd_pft_results <- padd_rd_photo %>% 
  rbind(padd_rd_myc) %>% 
  rbind(padd_rd_nfix) %>%
  mutate(k = 30)

################
# N+P addition
################

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "rd" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
npadd_rd_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "np" & 
                                   myvar == "rd" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_rd_photo <- data.frame(trait = "rd", 
                              nut_add = "np",
                              mod = "photo",
                              mod_results(npadd_rd_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_rd_pft))[2,3],
                              p = coef(summary(npadd_rd_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_rd_myc <- data.frame(trait = "rd", 
                            nut_add = "np",
                            mod = "myc_nas",
                            mod_results(npadd_rd_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(npadd_rd_pft))[3,3],
                            p = coef(summary(npadd_rd_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_rd_nfix <- data.frame(trait = "rd", 
                             nut_add = "np",
                             mod = "nfix",
                             mod_results(npadd_rd_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_rd_pft))[4,3],
                             p = coef(summary(npadd_rd_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_rd_pft_results <- npadd_rd_photo %>% 
  rbind(npadd_rd_myc) %>% 
  rbind(npadd_rd_nfix) %>%
  mutate(k = 32)

#############
# Merge Rd moderator results, with some light cleaning
#############
rd_pft_summary <- rbind(nadd_rd_pft_results, 
                        padd_rd_pft_results, 
                        npadd_rd_pft_results) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Vcmax - plant functional type
##############################################################################

##############
# N addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "vcmax" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
nadd_vcmax_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "n" & 
                                   myvar == "vcmax" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_vcmax_photo <- data.frame(trait = "vcmax", 
                              nut_add = "n",
                              mod = "photo",
                              mod_results(nadd_vcmax_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_vcmax_pft))[2,3],
                              p = coef(summary(nadd_vcmax_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_vcmax_myc <- data.frame(trait = "vcmax", 
                            nut_add = "n",
                            mod = "myc_nas",
                            mod_results(nadd_vcmax_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(nadd_vcmax_pft))[3,3],
                            p = coef(summary(nadd_vcmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_vcmax_nfix <- data.frame(trait = "vcmax", 
                             nut_add = "n",
                             mod = "nfix",
                             mod_results(nadd_vcmax_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_vcmax_pft))[4,3],
                             p = coef(summary(nadd_vcmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_vcmax_pft_results <- nadd_vcmax_photo %>% 
  rbind(nadd_vcmax_myc) %>% 
  rbind(nadd_vcmax_nfix) %>%
  mutate(k = 41)

##############
# P addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "vcmax" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
padd_vcmax_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "p" & 
                                   myvar == "vcmax" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_vcmax_photo <- data.frame(trait = "vcmax", 
                              nut_add = "p",
                              mod = "photo",
                              mod_results(padd_vcmax_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_vcmax_pft))[2,3],
                              p = coef(summary(padd_vcmax_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_vcmax_myc <- data.frame(trait = "vcmax", 
                            nut_add = "p",
                            mod = "myc_nas",
                            mod_results(padd_vcmax_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(padd_vcmax_pft))[3,3],
                            p = coef(summary(padd_vcmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_vcmax_nfix <- data.frame(trait = "vcmax", 
                             nut_add = "p",
                             mod = "nfix",
                             mod_results(padd_vcmax_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_vcmax_pft))[4,3],
                             p = coef(summary(padd_vcmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_vcmax_pft_results <- padd_vcmax_photo %>% 
  rbind(padd_vcmax_myc) %>% 
  rbind(padd_vcmax_nfix) %>%
  mutate(k = 42)

################
# N+P addition
################

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "vcmax" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
npadd_vcmax_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "np" & 
                                    myvar == "vcmax" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_vcmax_photo <- data.frame(trait = "vcmax", 
                               nut_add = "np",
                               mod = "photo",
                               mod_results(npadd_vcmax_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_vcmax_pft))[2,3],
                               p = coef(summary(npadd_vcmax_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_vcmax_myc <- data.frame(trait = "vcmax", 
                             nut_add = "np",
                             mod = "myc_nas",
                             mod_results(npadd_vcmax_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_vcmax_pft))[3,3],
                             p = coef(summary(npadd_vcmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_vcmax_nfix <- data.frame(trait = "vcmax", 
                              nut_add = "np",
                              mod = "nfix",
                              mod_results(npadd_vcmax_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_vcmax_pft))[4,3],
                              p = coef(summary(npadd_vcmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_vcmax_pft_results <- npadd_vcmax_photo %>% 
  rbind(npadd_vcmax_myc) %>% 
  rbind(npadd_vcmax_nfix) %>%
  mutate(k = 41)

#############
# Merge Vcmax moderator results, with some light cleaning
#############
vcmax_pft_summary <- rbind(nadd_vcmax_pft_results, 
                             padd_vcmax_pft_results, 
                             npadd_vcmax_pft_results) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Jmax - plant functional type
##############################################################################

##############
# N addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "jmax" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
nadd_jmax_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "n" & 
                                    myvar == "jmax" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_jmax_photo <- data.frame(trait = "jmax", 
                               nut_add = "n",
                               mod = "photo",
                               mod_results(nadd_jmax_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(nadd_jmax_pft))[2,3],
                               p = coef(summary(nadd_jmax_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_jmax_myc <- data.frame(trait = "jmax", 
                             nut_add = "n",
                             mod = "myc_nas",
                             mod_results(nadd_jmax_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_jmax_pft))[3,3],
                             p = coef(summary(nadd_jmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_jmax_nfix <- data.frame(trait = "jmax", 
                              nut_add = "n",
                              mod = "nfix",
                              mod_results(nadd_jmax_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_jmax_pft))[4,3],
                              p = coef(summary(nadd_jmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_jmax_pft_results <- nadd_jmax_photo %>% 
  rbind(nadd_jmax_myc) %>% 
  rbind(nadd_jmax_nfix) %>%
  mutate(k = 39)

##############
# P addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "jmax" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
padd_jmax_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "p" & 
                                    myvar == "jmax" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_jmax_photo <- data.frame(trait = "jmax", 
                               nut_add = "p",
                               mod = "photo",
                               mod_results(padd_jmax_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(padd_jmax_pft))[2,3],
                               p = coef(summary(padd_jmax_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_jmax_myc <- data.frame(trait = "jmax", 
                             nut_add = "p",
                             mod = "myc_nas",
                             mod_results(padd_jmax_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_jmax_pft))[3,3],
                             p = coef(summary(padd_jmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_jmax_nfix <- data.frame(trait = "jmax", 
                              nut_add = "p",
                              mod = "nfix",
                              mod_results(padd_jmax_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_jmax_pft))[4,3],
                              p = coef(summary(padd_jmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_jmax_pft_results <- padd_jmax_photo %>% 
  rbind(padd_jmax_myc) %>% 
  rbind(padd_jmax_nfix) %>%
  mutate(k = 40)

################
# N+P addition
################

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "jmax" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
npadd_jmax_pft <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path + myc_nas + n_fixer,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(nut_add == "np" & 
                                     myvar == "jmax" & 
                                     !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_jmax_photo <- data.frame(trait = "jmax", 
                                nut_add = "np",
                                mod = "photo",
                                mod_results(npadd_jmax_pft, 
                                            mod = "photo_path", 
                                            group = "exp")$mod_table,
                                z = coef(summary(npadd_jmax_pft))[2,3],
                                p = coef(summary(npadd_jmax_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_jmax_myc <- data.frame(trait = "jmax", 
                              nut_add = "np",
                              mod = "myc_nas",
                              mod_results(npadd_jmax_pft, 
                                          mod = "myc_nas", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_jmax_pft))[3,3],
                              p = coef(summary(npadd_jmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_jmax_nfix <- data.frame(trait = "jmax", 
                               nut_add = "np",
                               mod = "nfix",
                               mod_results(npadd_jmax_pft, 
                                           mod = "n_fixer", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_jmax_pft))[4,3],
                               p = coef(summary(npadd_jmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_jmax_pft_results <- npadd_jmax_photo %>% 
  rbind(npadd_jmax_myc) %>% 
  rbind(npadd_jmax_nfix) %>%
  mutate(k = 39)

#############
# Merge Jmax moderator results, with some light cleaning
#############
jmax_pft_summary <- rbind(nadd_jmax_pft_results, 
                              padd_jmax_pft_results, 
                              npadd_jmax_pft_results) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Jmax:Vcmax - plant functional type
##############################################################################

##############
# N addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "jmax_vcmax" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
nadd_jmaxvcmax_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "n" & 
                                   myvar == "jmax_vcmax" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_jmaxvcmax_photo <- data.frame(trait = "jmax_vcmax", 
                              nut_add = "n",
                              mod = "photo",
                              mod_results(nadd_jmaxvcmax_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_jmaxvcmax_pft))[2,3],
                              p = coef(summary(nadd_jmaxvcmax_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_jmaxvcmax_myc <- data.frame(trait = "jmax_vcmax", 
                            nut_add = "n",
                            mod = "myc_nas",
                            mod_results(nadd_jmaxvcmax_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(nadd_jmaxvcmax_pft))[3,3],
                            p = coef(summary(nadd_jmaxvcmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_jmaxvcmax_nfix <- data.frame(trait = "jmax_vcmax", 
                             nut_add = "n",
                             mod = "nfix",
                             mod_results(nadd_jmaxvcmax_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_jmaxvcmax_pft))[4,3],
                             p = coef(summary(nadd_jmaxvcmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_jmaxvcmax_pft_results <- nadd_jmaxvcmax_photo %>% 
  rbind(nadd_jmaxvcmax_myc) %>% 
  rbind(nadd_jmaxvcmax_nfix) %>%
  mutate(k = 31)

##############
# P addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "jmax_vcmax" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
padd_jmaxvcmax_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "p" & 
                                   myvar == "jmax_vcmax" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_jmaxvcmax_photo <- data.frame(trait = "jmax_vcmax", 
                              nut_add = "p",
                              mod = "photo",
                              mod_results(padd_jmaxvcmax_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_jmaxvcmax_pft))[2,3],
                              p = coef(summary(padd_jmaxvcmax_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_jmaxvcmax_myc <- data.frame(trait = "jmax_vcmax", 
                            nut_add = "p",
                            mod = "myc_nas",
                            mod_results(padd_jmaxvcmax_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(padd_jmaxvcmax_pft))[3,3],
                            p = coef(summary(padd_jmaxvcmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_jmaxvcmax_nfix <- data.frame(trait = "jmax_vcmax", 
                             nut_add = "p",
                             mod = "nfix",
                             mod_results(padd_jmaxvcmax_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_jmaxvcmax_pft))[4,3],
                             p = coef(summary(padd_jmaxvcmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_jmaxvcmax_pft_results <- padd_jmaxvcmax_photo %>% 
  rbind(padd_jmaxvcmax_myc) %>% 
  rbind(padd_jmaxvcmax_nfix) %>%
  mutate(k = 32)

################
# N+P addition
################

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "jmax_vcmax" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
npadd_jmaxvcmax_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "np" & 
                                    myvar == "jmax_vcmax" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_jmaxvcmax_photo <- data.frame(trait = "jmax_vcmax", 
                               nut_add = "np",
                               mod = "photo",
                               mod_results(npadd_jmaxvcmax_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_jmaxvcmax_pft))[2,3],
                               p = coef(summary(npadd_jmaxvcmax_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_jmaxvcmax_myc <- data.frame(trait = "jmax_vcmax", 
                             nut_add = "np",
                             mod = "myc_nas",
                             mod_results(npadd_jmaxvcmax_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_jmaxvcmax_pft))[3,3],
                             p = coef(summary(npadd_jmaxvcmax_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_jmaxvcmax_nfix <- data.frame(trait = "jmax_vcmax", 
                              nut_add = "np",
                              mod = "nfix",
                              mod_results(npadd_jmaxvcmax_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_jmaxvcmax_pft))[4,3],
                              p = coef(summary(npadd_jmaxvcmax_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_jmaxvcmax_pft_results <- npadd_jmaxvcmax_photo %>% 
  rbind(npadd_jmaxvcmax_myc) %>% 
  rbind(npadd_jmaxvcmax_nfix) %>%
  mutate(k = 30)

#############
# Merge Jmax:Vcmax moderator results, with some light cleaning
#############
jmaxvcmax_pft_summary <- rbind(nadd_jmaxvcmax_pft_results, 
                                  padd_jmaxvcmax_pft_results, 
                                  npadd_jmaxvcmax_pft_results) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# PNUE - plant functional type
##############################################################################

##############
# N addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "leaf_pnue" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
nadd_pnue_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "n" & 
                                   myvar == "leaf_pnue" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_pnue_photo <- data.frame(trait = "pnue", 
                                   nut_add = "n",
                                   mod = "photo",
                                   mod_results(nadd_pnue_pft, 
                                               mod = "photo_path", 
                                               group = "exp")$mod_table,
                                   z = coef(summary(nadd_pnue_pft))[2,3],
                                   p = coef(summary(nadd_pnue_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_pnue_myc <- data.frame(trait = "pnue", 
                                 nut_add = "n",
                                 mod = "myc_nas",
                                 mod_results(nadd_pnue_pft, 
                                             mod = "myc_nas", 
                                             group = "exp")$mod_table,
                                 z = coef(summary(nadd_pnue_pft))[3,3],
                                 p = coef(summary(nadd_pnue_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_pnue_nfix <- data.frame(trait = "pnue", 
                                  nut_add = "n",
                                  mod = "nfix",
                                  mod_results(nadd_pnue_pft, 
                                              mod = "n_fixer", 
                                              group = "exp")$mod_table,
                                  z = coef(summary(nadd_pnue_pft))[4,3],
                                  p = coef(summary(nadd_pnue_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_pnue_pft_results <- nadd_pnue_photo %>% 
  rbind(nadd_pnue_myc) %>% 
  rbind(nadd_pnue_nfix) %>%
  mutate(k = 64)

##############
# P addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "leaf_pnue" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
padd_pnue_pft <- rma.mv(logr, 
                             logr_var,
                             method = "REML", 
                             random = ~ 1 | exp, 
                             mods = ~ photo_path + myc_nas + n_fixer,
                             slab = exp, 
                             control = list(stepadj = 0.3), 
                             data = meta_results %>% 
                               filter(nut_add == "p" & 
                                        myvar == "leaf_pnue" & 
                                        !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_pnue_photo <- data.frame(trait = "pnue", 
                                   nut_add = "p",
                                   mod = "photo",
                                   mod_results(padd_pnue_pft, 
                                               mod = "photo_path", 
                                               group = "exp")$mod_table,
                                   z = coef(summary(padd_pnue_pft))[2,3],
                                   p = coef(summary(padd_pnue_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_pnue_myc <- data.frame(trait = "pnue", 
                                 nut_add = "p",
                                 mod = "myc_nas",
                                 mod_results(padd_pnue_pft, 
                                             mod = "myc_nas", 
                                             group = "exp")$mod_table,
                                 z = coef(summary(padd_pnue_pft))[3,3],
                                 p = coef(summary(padd_pnue_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_pnue_nfix <- data.frame(trait = "pnue", 
                                  nut_add = "p",
                                  mod = "nfix",
                                  mod_results(padd_pnue_pft, 
                                              mod = "n_fixer", 
                                              group = "exp")$mod_table,
                                  z = coef(summary(padd_pnue_pft))[4,3],
                                  p = coef(summary(padd_pnue_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_pnue_pft_results <- padd_pnue_photo %>% 
  rbind(padd_pnue_myc) %>% 
  rbind(padd_pnue_nfix) %>%
  mutate(k = 64)

################
# N+P addition
################

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "leaf_pnue" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
npadd_pnue_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "np" & 
                                    myvar == "leaf_pnue" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_pnue_photo <- data.frame(trait = "pnue", 
                                    nut_add = "np",
                                    mod = "photo",
                                    mod_results(npadd_pnue_pft, 
                                                mod = "photo_path", 
                                                group = "exp")$mod_table,
                                    z = coef(summary(npadd_pnue_pft))[2,3],
                                    p = coef(summary(npadd_pnue_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_pnue_myc <- data.frame(trait = "pnue", 
                                  nut_add = "np",
                                  mod = "myc_nas",
                                  mod_results(npadd_pnue_pft, 
                                              mod = "myc_nas", 
                                              group = "exp")$mod_table,
                                  z = coef(summary(npadd_pnue_pft))[3,3],
                                  p = coef(summary(npadd_pnue_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_pnue_nfix <- data.frame(trait = "pnue", 
                                   nut_add = "np",
                                   mod = "nfix",
                                   mod_results(npadd_pnue_pft, 
                                               mod = "n_fixer", 
                                               group = "exp")$mod_table,
                                   z = coef(summary(npadd_pnue_pft))[4,3],
                                   p = coef(summary(npadd_pnue_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_pnue_pft_results <- npadd_pnue_photo %>% 
  rbind(npadd_pnue_myc) %>% 
  rbind(npadd_pnue_nfix) %>%
  mutate(k = 64)

#############
# Merge PNUE moderator results, with some light cleaning
#############
pnue_pft_summary <- rbind(nadd_pnue_pft_results, 
                          padd_pnue_pft_results, 
                          npadd_pnue_pft_results) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# PPUE - plant functional type
##############################################################################

##############
# N addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "n" & 
                                        myvar == "leaf_ppue" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
nadd_ppue_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "n" & 
                                   myvar == "leaf_ppue" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
nadd_ppue_photo <- data.frame(trait = "ppue", 
                              nut_add = "n",
                              mod = "photo",
                              mod_results(nadd_ppue_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(nadd_ppue_pft))[2,3],
                              p = coef(summary(nadd_ppue_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
nadd_ppue_myc <- data.frame(trait = "ppue", 
                            nut_add = "n",
                            mod = "myc_nas",
                            mod_results(nadd_ppue_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(nadd_ppue_pft))[3,3],
                            p = coef(summary(nadd_ppue_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
nadd_ppue_nfix <- data.frame(trait = "ppue", 
                             nut_add = "n",
                             mod = "nfix",
                             mod_results(nadd_ppue_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(nadd_ppue_pft))[4,3],
                             p = coef(summary(nadd_ppue_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
nadd_ppue_pft_results <- nadd_ppue_photo %>% 
  rbind(nadd_ppue_myc) %>% 
  rbind(nadd_ppue_nfix) %>%
  mutate(k = 66)

##############
# P addition
##############

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "p" & 
                                        myvar == "leaf_ppue" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
padd_ppue_pft <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ photo_path + myc_nas + n_fixer,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(nut_add == "p" & 
                                   myvar == "leaf_ppue" & 
                                   !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
padd_ppue_photo <- data.frame(trait = "ppue", 
                              nut_add = "p",
                              mod = "photo",
                              mod_results(padd_pnue_pft, 
                                          mod = "photo_path", 
                                          group = "exp")$mod_table,
                              z = coef(summary(padd_ppue_pft))[2,3],
                              p = coef(summary(padd_ppue_pft))[2,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
padd_ppue_myc <- data.frame(trait = "ppue", 
                            nut_add = "p",
                            mod = "myc_nas",
                            mod_results(padd_ppue_pft, 
                                        mod = "myc_nas", 
                                        group = "exp")$mod_table,
                            z = coef(summary(padd_ppue_pft))[3,3],
                            p = coef(summary(padd_ppue_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
padd_ppue_nfix <- data.frame(trait = "ppue", 
                             nut_add = "p",
                             mod = "nfix",
                             mod_results(padd_ppue_pft, 
                                         mod = "n_fixer", 
                                         group = "exp")$mod_table,
                             z = coef(summary(padd_ppue_pft))[4,3],
                             p = coef(summary(padd_ppue_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
padd_ppue_pft_results <- padd_ppue_photo %>% 
  rbind(padd_ppue_myc) %>% 
  rbind(padd_ppue_nfix) %>%
  mutate(k = 65)

################
# N+P addition
################

# Visualize data distribution
ggplot(data = meta_results %>% filter(nut_add == "np" & 
                                        myvar == "leaf_ppue" & 
                                        !is.na(photo_path))) +
  geom_point(aes(x = photo_path, y = logr))

# Model
npadd_ppue_pft <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path + myc_nas + n_fixer,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(nut_add == "np" & 
                                    myvar == "leaf_ppue" & 
                                    !is.na(photo_path)))

# Extract photosynthetic pathway summary statistics
npadd_ppue_photo <- data.frame(trait = "ppue", 
                               nut_add = "np",
                               mod = "photo",
                               mod_results(npadd_ppue_pft, 
                                           mod = "photo_path", 
                                           group = "exp")$mod_table,
                               z = coef(summary(npadd_ppue_pft))[2,3],
                               p = coef(summary(npadd_ppue_pft))[2, 4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract mycorrhizal acquisition strategy summary statistics
npadd_ppue_myc <- data.frame(trait = "ppue", 
                             nut_add = "np",
                             mod = "myc_nas",
                             mod_results(npadd_ppue_pft, 
                                         mod = "myc_nas", 
                                         group = "exp")$mod_table,
                             z = coef(summary(npadd_ppue_pft))[3,3],
                             p = coef(summary(npadd_ppue_pft))[3,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Extract N-fixation ability summary statistics
npadd_ppue_nfix <- data.frame(trait = "ppue", 
                              nut_add = "np",
                              mod = "nfix",
                              mod_results(npadd_ppue_pft, 
                                          mod = "n_fixer", 
                                          group = "exp")$mod_table,
                              z = coef(summary(npadd_ppue_pft))[4,3],
                              p = coef(summary(npadd_ppue_pft))[4,4]) %>%
  dplyr::select(trait, nut_add, mod, comp = name, estimate, z, p, lowerCL, upperCL)

# Merge summary statistics into single data frame
npadd_ppue_pft_results <- npadd_ppue_photo %>% 
  rbind(npadd_ppue_myc) %>% 
  rbind(npadd_ppue_nfix) %>%
  mutate(k = 66)

#############
# Merge PPUE moderator results, with some light cleaning
#############
ppue_pft_summary <- rbind(nadd_ppue_pft_results, 
                          padd_ppue_pft_results, 
                          npadd_ppue_pft_results) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p))) %>%
  dplyr::select(trait:comp, k, estimate:upperCL)

##############################################################################
# Merge PFT summary data frames, write to .csv
##############################################################################
marea_pft_summary %>%
  rbind(nmass_pft_summary) %>%
  rbind(narea_pft_summary) %>%
  rbind(pmass_pft_summary) %>%
  rbind(parea_pft_summary) %>%
  rbind(leafnp_pft_summary) %>%
  rbind(asat_pft_summary) %>%
  rbind(gsw_pft_summary) %>%
  rbind(rd_pft_summary) %>%
  rbind(vcmax_pft_summary) %>%
  rbind(jmax_pft_summary) %>%
  rbind(jmaxvcmax_pft_summary) %>%
  rbind(pnue_pft_summary) %>%
  rbind(ppue_pft_summary) %>%
  mutate(ci_range = str_c("[", sprintf("%.3f", lowerCL), ", ", sprintf("%.3f", upperCL), "]")) %>%
  write.csv("../data/CNPmeta_pft_moderators.csv", 
            row.names = F)








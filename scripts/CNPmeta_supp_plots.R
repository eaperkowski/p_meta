###############################################################################
# Script prep
###############################################################################

# Libraries
library(tidyverse)
library(metafor)
library(ggpubr)
library(forcats)
library(patchwork)
library(naniar)
library(orchaRd)
library(ggh4x)

# Read in .csv files from CNPmeta_analysis.R
spp_mod_df <- read.csv("../data/CNP_spp_moderators.csv") %>%
  mutate(moderator = factor(moderator, levels = c("c4", "c3",
                                                  "n_fixer", "non_fixer",
                                                  "scavenging", "mining")),
         nut_add = factor(nut_add, levels = c("n", "p", "np")))
meta_results <- read.csv("../data/CNP_metaanalysis_results.csv")

# Add code for facet labels
facet.labs <- c("N add.", "P add.", "N+P add.")
names(facet.labs) <- c("n", "p", "np")

###############################################################################
# LMA - species identity
###############################################################################

# Photosynthetic pathway
marea_photo_plot <- ggplot(data = spp_mod_df %>% filter(trait == "marea" & 
                                      moderator_group == "photo"),
       aes(group = nut_add)) +
  geom_errorbar(aes(x = moderator, y = estimate, 
                    ymin = lower_ci, ymax = upper_ci),
                width = 0.5, linewidth = 1) +
  geom_point(aes(x = moderator, y = estimate, fill = moderator),
             shape = 21, size = 5) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "marea" & 
                                              moderator_group == "photo" &
                                              nut_add == "n"),
               xmin = "c4", xmax = "c3", y.position = 0.3, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "marea" & 
                                           moderator_group == "photo" &
                                           nut_add == "n"),
            x = 1.5, y = 0.35,
            label = "**", size = 6) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "marea" & 
                                              moderator_group == "photo" &
                                              nut_add == "np"),
               xmin = "c4", xmax = "c3", y.position = 0.3, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "marea" & 
                                           moderator_group == "photo" &
                                           nut_add == "np"),
            x = 1.5, y = 0.35,
            label = "***", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  scale_x_discrete(labels = c("C4", "C3")) +
  labs(x = "Photo. pathway",
       y = NULL) +
  coord_flip() +
  facet_grid(nut_add~., labeller = labeller(nut_add = facet.labs)) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(strip.background = element_rect(linewidth = 1.5),
        strip.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"))
marea_photo_plot


# N fixation
marea_nfix_plot <- ggplot(data = spp_mod_df %>% filter(trait == "marea" & 
                                      moderator_group == "nfix"),
       aes(group = nut_add)) +
  geom_errorbar(aes(x = moderator, y = estimate, 
                    ymin = lower_ci, ymax = upper_ci),
                width = 0.5, linewidth = 1) +
  geom_point(aes(x = moderator, y = estimate, fill = moderator),
             shape = 21, size = 5) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "marea" & 
                                              moderator_group == "nfix" &
                                              nut_add == "n"),
               xmin = "n_fixer", xmax = "non_fixer", y.position = 0.3, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "marea" & 
                                           moderator_group == "nfix" &
                                           nut_add == "n"),
            x = 1.5, y = 0.35,
            label = "***", size = 6) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "marea" & 
                                              moderator_group == "nfix" &
                                              nut_add == "p"),
               xmin = "n_fixer", xmax = "non_fixer", y.position = 0.3, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "marea" & 
                                           moderator_group == "nfix" &
                                           nut_add == "p"),
            x = 1.5, y = 0.35,
            label = "***", size = 6) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "marea" & 
                                              moderator_group == "nfix" &
                                              nut_add == "np"),
               xmin = "n_fixer", xmax = "non_fixer", y.position = 0.3, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "marea" & 
                                           moderator_group == "nfix" &
                                           nut_add == "np"),
            x = 1.5, y = 0.35,
            label = "*", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  scale_x_discrete(labels = c("N-fixer", "non-fixer")) +
  labs(x = "N-fixation ability",
       y = NULL) +
  coord_flip() +
  facet_grid(nut_add~., labeller = labeller(nut_add = facet.labs)) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(strip.background = element_rect(linewidth = 1.5),
        strip.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"))

# Mycorrhizal acquisition strategy
marea_myc_plot <- ggplot(data = spp_mod_df %>% filter(trait == "marea" & 
                                      moderator_group == "myc"),
       aes(group = nut_add)) +
  geom_errorbar(aes(x = moderator, y = estimate, 
                    ymin = lower_ci, ymax = upper_ci),
                width = 0.5, linewidth = 1) +
  geom_point(aes(x = moderator, y = estimate, fill = moderator),
             shape = 21, size = 5) +
  
  geom_bracket(data = spp_mod_df %>% filter(trait == "marea" & 
                                              moderator_group == "myc" &
                                              nut_add == "n"),
               xmin = "mining", xmax = "scavenging", y.position = 0.3, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "marea" & 
                                           moderator_group == "myc" &
                                           nut_add == "n"),
            x = 1.5, y = 0.35,
            label = "*", size = 6) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "marea" & 
                                              moderator_group == "myc" &
                                              nut_add == "np"),
               xmin = "mining", xmax = "scavenging", y.position = 0.3, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "marea" & 
                                           moderator_group == "myc" &
                                           nut_add == "np"),
            x = 1.5, y = 0.35,
            label = "***", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  labs(x = "Myc. acquisition strategy",
       y = "Log-response ratio") +
  coord_flip() +
  facet_grid(nut_add~., labeller = labeller(nut_add = facet.labs)) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(strip.background = element_rect(linewidth = 1.5),
        strip.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"))

# png("../plots/CNPmeta_spp_moderator/CNPmeta_marea_sppMod.png",
#     height = 12.5, width = 8, units = "in", res = 600)
ggarrange(marea_photo_plot, marea_nfix_plot, marea_myc_plot,
          ncol = 1, nrow = 3, align = "hv", labels = c("(a)", "(b)", "(c)"),
          font.label = list(size = 18), hjust = -0.5)
# dev.off()

###############################################################################
# Nmass - species identity
###############################################################################

# Photosynthetic pathway
nmass_photo_plot <- ggplot(data = spp_mod_df %>% 
                             filter(trait == "nmass" &
                                      moderator_group == "photo"),
                           aes(group = nut_add)) +
  geom_errorbar(aes(x = moderator, y = estimate, 
                    ymin = lower_ci, ymax = upper_ci),
                width = 0.5, linewidth = 1) +
  geom_point(aes(x = moderator, y = estimate, fill = moderator),
             shape = 21, size = 5) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "nmass" & 
                                              moderator_group == "photo" &
                                              nut_add == "n"),
               xmin = "c4", xmax = "c3", y.position = 0.5, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "nmass" & 
                                           moderator_group == "photo" &
                                           nut_add == "n"),
            x = 1.5, y = 0.55,
            label = "*", size = 6) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "nmass" & 
                                              moderator_group == "photo" &
                                              nut_add == "np"),
               xmin = "c4", xmax = "c3", y.position = 0.5, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "nmass" & 
                                           moderator_group == "photo" &
                                           nut_add == "np"),
            x = 1.5, y = 0.55,
            label = "**", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.2, 0.6)) +
  scale_x_discrete(labels = c("C4", "C3")) +
  labs(x = "Photo. pathway",
       y = NULL) +
  coord_flip() +
  facet_grid(nut_add~., labeller = labeller(nut_add = facet.labs)) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(strip.background = element_rect(linewidth = 1.5),
        strip.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"))
nmass_photo_plot


# N fixation
nmass_nfix_plot <- ggplot(data = spp_mod_df %>%
                            filter(trait == "nmass" &
                                     moderator_group == "nfix"),
                          aes(group = nut_add)) +
  geom_errorbar(aes(x = moderator, y = estimate, 
                    ymin = lower_ci, ymax = upper_ci),
                width = 0.5, linewidth = 1) +
  geom_point(aes(x = moderator, y = estimate, fill = moderator),
             shape = 21, size = 5) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "nmass" & 
                                              moderator_group == "nfix" &
                                              nut_add == "n"),
               xmin = "n_fixer", xmax = "non_fixer", y.position = 0.5, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "nmass" & 
                                           moderator_group == "nfix" &
                                           nut_add == "n"),
            x = 1.5, y = 0.55,
            label = "***", size = 6) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "nmass" & 
                                              moderator_group == "nfix" &
                                              nut_add == "p"),
               xmin = "n_fixer", xmax = "non_fixer", y.position = 0.5, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "nmass" & 
                                           moderator_group == "nfix" &
                                           nut_add == "p"),
            x = 1.5, y = 0.55,
            label = "***", size = 6) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "nmass" & 
                                              moderator_group == "nfix" &
                                              nut_add == "np"),
               xmin = "n_fixer", xmax = "non_fixer", y.position = 0.5, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "nmass" & 
                                           moderator_group == "nfix" &
                                           nut_add == "np"),
            x = 1.5, y = 0.55,
            label = "*", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.2, 0.6)) +
  scale_x_discrete(labels = c("N-fixer", "non-fixer")) +
  labs(x = "N-fixation ability",
       y = NULL) +
  coord_flip() +
  facet_grid(nut_add~., labeller = labeller(nut_add = facet.labs)) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(strip.background = element_rect(linewidth = 1.5),
        strip.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"))

# Mycorrhizal acquisition strategy
nmass_myc_plot <- ggplot(data = spp_mod_df %>%
                           filter(trait == "nmass" &
                                    moderator_group == "myc"),
                         aes(group = nut_add)) +
  geom_errorbar(aes(x = moderator, y = estimate, 
                    ymin = lower_ci, ymax = upper_ci),
                width = 0.5, linewidth = 1) +
  geom_point(aes(x = moderator, y = estimate, fill = moderator),
             shape = 21, size = 5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.2, 0.6)) +
  labs(x = "Myc. acquisition strategy",
       y = "Log-response ratio") +
  coord_flip() +
  facet_grid(nut_add~., labeller = labeller(nut_add = facet.labs)) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(strip.background = element_rect(linewidth = 1.5),
        strip.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"))

png("../plots/CNPmeta_spp_moderator/CNPmeta_nmass_sppMod.png",
    height = 12.5, width = 8, units = "in", res = 600)
ggarrange(nmass_photo_plot, nmass_nfix_plot, nmass_myc_plot,
          ncol = 1, nrow = 3, align = "hv", labels = c("(a)", "(b)", "(c)"),
          font.label = list(size = 18), hjust = -0.5)
dev.off()


###############################################################################
# Nmass - species identity
###############################################################################

# Photosynthetic pathway
nmass_photo_plot <- ggplot(data = spp_mod_df %>% 
                             filter(trait == "narea" &
                                      moderator_group == "photo"),
                           aes(group = nut_add)) +
  geom_errorbar(aes(x = moderator, y = estimate, 
                    ymin = lower_ci, ymax = upper_ci),
                width = 0.5, linewidth = 1) +
  geom_point(aes(x = moderator, y = estimate, fill = moderator),
             shape = 21, size = 5) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "narea" & 
                                              moderator_group == "photo" &
                                              nut_add == "n"),
               xmin = "c4", xmax = "c3", y.position = 0.5, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "narea" & 
                                           moderator_group == "photo" &
                                           nut_add == "n"),
            x = 1.5, y = 0.55,
            label = "**", size = 6) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "narea" & 
                                              moderator_group == "photo" &
                                              nut_add == "np"),
               xmin = "c4", xmax = "c3", y.position = 0.5, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "narea" & 
                                           moderator_group == "photo" &
                                           nut_add == "np"),
            x = 1.5, y = 0.55,
            label = "***", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.2, 0.6)) +
  scale_x_discrete(labels = c("C4", "C3")) +
  labs(x = "Photo. pathway",
       y = NULL) +
  coord_flip() +
  facet_grid(nut_add~., labeller = labeller(nut_add = facet.labs)) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(strip.background = element_rect(linewidth = 1.5),
        strip.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"))
nmass_photo_plot


# N fixation
nmass_nfix_plot <- ggplot(data = spp_mod_df %>%
                            filter(trait == "nmass" &
                                     moderator_group == "nfix"),
                          aes(group = nut_add)) +
  geom_errorbar(aes(x = moderator, y = estimate, 
                    ymin = lower_ci, ymax = upper_ci),
                width = 0.5, linewidth = 1) +
  geom_point(aes(x = moderator, y = estimate, fill = moderator),
             shape = 21, size = 5) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "nmass" & 
                                              moderator_group == "nfix" &
                                              nut_add == "n"),
               xmin = "n_fixer", xmax = "non_fixer", y.position = 0.5, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "nmass" & 
                                           moderator_group == "nfix" &
                                           nut_add == "n"),
            x = 1.5, y = 0.55,
            label = "***", size = 6) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "nmass" & 
                                              moderator_group == "nfix" &
                                              nut_add == "p"),
               xmin = "n_fixer", xmax = "non_fixer", y.position = 0.5, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "nmass" & 
                                           moderator_group == "nfix" &
                                           nut_add == "p"),
            x = 1.5, y = 0.55,
            label = "***", size = 6) +
  geom_bracket(data = spp_mod_df %>% filter(trait == "narea" & 
                                              moderator_group == "nfix" &
                                              nut_add == "np"),
               xmin = "n_fixer", xmax = "non_fixer", y.position = 0.5, size = 1,
               label = "", tip.length = 0.02) +
  geom_text(data = spp_mod_df %>% filter(trait == "nmass" & 
                                           moderator_group == "nfix" &
                                           nut_add == "np"),
            x = 1.5, y = 0.55,
            label = "*", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.2, 0.6)) +
  scale_x_discrete(labels = c("N-fixer", "non-fixer")) +
  labs(x = "N-fixation ability",
       y = NULL) +
  coord_flip() +
  facet_grid(nut_add~., labeller = labeller(nut_add = facet.labs)) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(strip.background = element_rect(linewidth = 1.5),
        strip.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"))

# Mycorrhizal acquisition strategy
nmass_myc_plot <- ggplot(data = spp_mod_df %>%
                           filter(trait == "narea" &
                                    moderator_group == "myc"),
                         aes(group = nut_add)) +
  geom_errorbar(aes(x = moderator, y = estimate, 
                    ymin = lower_ci, ymax = upper_ci),
                width = 0.5, linewidth = 1) +
  geom_point(aes(x = moderator, y = estimate, fill = moderator),
             shape = 21, size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.2, 0.6)) +
  labs(x = "Myc. acquisition strategy",
       y = "Log-response ratio") +
  coord_flip() +
  facet_grid(nut_add~., labeller = labeller(nut_add = facet.labs)) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(strip.background = element_rect(linewidth = 1.5),
        strip.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"))

png("../plots/CNPmeta_spp_moderator/CNPmeta_narea_sppMod.png",
    height = 12.5, width = 8, units = "in", res = 600)
ggarrange(narea_photo_plot, narea_nfix_plot, narea_myc_plot,
          ncol = 1, nrow = 3, align = "hv", labels = c("(a)", "(b)", "(c)"),
          font.label = list(size = 18), hjust = -0.5)
dev.off()




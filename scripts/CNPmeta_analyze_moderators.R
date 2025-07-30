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

# N addition
nadd_marea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "lma" & 
                                     !is.na(gs_mat)))

nadd_marea_clim_summary <- data.frame(trait = "marea",
                                      nut_add = "n",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_marea_clim)),
                                      row.names = NULL)


# P addition
padd_marea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "lma" & 
                                     !is.na(gs_mat)))

padd_marea_clim_summary <- data.frame(trait = "marea",
                                      nut_add = "p",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_marea_clim)),
                                      row.names = NULL)

# N+P addition
npadd_marea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "lma" & 
                                     !is.na(gs_mat)))

npadd_marea_clim_summary <- data.frame(trait = "marea",
                                       nut_add = "np",
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_marea_clim)),
                                       row.names = NULL)

## Merge Marea moderator results, with some light cleaning
marea_clim_summary <- rbind(nadd_marea_clim_summary, 
      padd_marea_clim_summary, 
      npadd_marea_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 3)))

##############################################################################
# Nmass climate moderators
##############################################################################

# N addition
nadd_nmass_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_n_mass" & 
                                     !is.na(gs_mat) & logr > -0.2))

nadd_nmass_clim_summary <- data.frame(trait = "nmass",
                                      nut_add = "n",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_nmass_clim)),
                                      row.names = NULL)

# P addition
padd_nmass_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_n_mass" & 
                                     !is.na(gs_mat)))

padd_nmass_clim_summary <- data.frame(trait = "nmass",
                                      nut_add = "p",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_nmass_clim)),
                                      row.names = NULL)

# N+P addition
npadd_nmass_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_n_mass" & 
                                      !is.na(gs_mat)))

npadd_nmass_clim_summary <- data.frame(trait = "nmass",
                                       nut_add = "np",
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_nmass_clim)),
                                       row.names = NULL)

## Merge Nmass moderator results, with some light cleaning
nmass_clim_summary <- rbind(nadd_nmass_clim_summary, 
                            padd_nmass_clim_summary, 
                            npadd_nmass_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 3)))

##############################################################################
# Narea climate moderators
##############################################################################

# N addition
nadd_narea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_n_area" & 
                                     !is.na(gs_mat) & 
                                     gs_ai < 2.5 & 
                                     logr < 0.75))

nadd_narea_clim_summary <- data.frame(trait = "narea",
                                      nut_add = "n",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_narea_clim)),
                                      row.names = NULL)

# P addition
padd_narea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_n_area" & 
                                     !is.na(gs_mat)))

padd_narea_clim_summary <- data.frame(trait = "narea",
                                      nut_add = "p",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_narea_clim)),
                                      row.names = NULL)

# N+P addition
npadd_narea_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_n_area" & 
                                      !is.na(gs_mat)))

npadd_narea_clim_summary <- data.frame(trait = "narea",
                                       nut_add = "np",
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_narea_clim)),
                                       row.names = NULL)

## Merge Narea moderator results, with some light cleaning
narea_clim_summary <- rbind(nadd_narea_clim_summary, 
                            padd_narea_clim_summary, 
                            npadd_narea_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 3)))

##############################################################################
# Pmass climate moderators
##############################################################################

# N addition
nadd_pmass_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_p_mass" & 
                                     !is.na(gs_mat)))

nadd_pmass_clim_summary <- data.frame(trait = "pmass",
                                      nut_add = "n",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_pmass_clim)),
                                      row.names = NULL)

# P addition
padd_pmass_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_p_mass" & 
                                     !is.na(gs_mat)))

padd_pmass_clim_summary <- data.frame(trait = "pmass",
                                      nut_add = "p",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_pmass_clim)),
                                      row.names = NULL)

# N+P addition
npadd_pmass_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_p_mass" & 
                                      !is.na(gs_mat)))

npadd_pmass_clim_summary <- data.frame(trait = "pmass",
                                       nut_add = "np",
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_pmass_clim)),
                                       row.names = NULL)

## Merge Pmass moderator results, with some light cleaning
pmass_clim_summary <- rbind(nadd_pmass_clim_summary, 
                            padd_pmass_clim_summary, 
                            npadd_pmass_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 3)))

##############################################################################
# Parea climate moderators
##############################################################################

# N addition
nadd_parea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_p_area" & 
                                     !is.na(gs_mat)))

nadd_parea_clim_summary <- data.frame(trait = "parea",
                                      nut_add = "n",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_parea_clim)),
                                      row.names = NULL)

# P addition
padd_parea_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_p_area" & 
                                     !is.na(gs_mat) & logr > -0.5))

padd_parea_clim_summary <- data.frame(trait = "parea",
                                      nut_add = "p",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_parea_clim)),
                                      row.names = NULL)

# N+P addition
npadd_parea_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_p_area" & 
                                      !is.na(gs_mat)))

npadd_parea_clim_summary <- data.frame(trait = "parea",
                                       nut_add = "np",
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_parea_clim)),
                                       row.names = NULL)

## Merge Parea moderator results, with some light cleaning
parea_clim_summary <- rbind(nadd_parea_clim_summary, 
                            padd_parea_clim_summary, 
                            npadd_parea_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 3)))

##############################################################################
# Total biomass climate moderators
##############################################################################

# N addition
nadd_tbio_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "total_biomass" & 
                                     !is.na(gs_mat)))

nadd_tbio_clim_summary <- data.frame(trait = "tbio",
                                      nut_add = "n",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(nadd_tbio_clim)),
                                      row.names = NULL)

# P addition
padd_tbio_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "total_biomass" & 
                                     !is.na(gs_mat)))

padd_tbio_clim_summary <- data.frame(trait = "tbio",
                                      nut_add = "p",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(padd_tbio_clim)),
                                      row.names = NULL)

# N+P addition
npadd_tbio_clim <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ gs_mat + gs_ai + gs_par,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "total_biomass" & 
                                      !is.na(gs_mat)))

npadd_tbio_clim_summary <- data.frame(trait = "tbio",
                                       nut_add = "np",
                                       mod = c("intrcpt", "gs_mat",
                                               "gs_ai", "gs_par"),
                                       coef(summary(npadd_tbio_clim)),
                                       row.names = NULL)

## Merge Parea moderator results, with some light cleaning
tbio_clim_summary <- rbind(nadd_tbio_clim_summary, 
                            padd_tbio_clim_summary, 
                            npadd_tbio_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 3)))

##############################################################################
# Aboveground biomass climate moderators
##############################################################################

# N addition
nadd_agb_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "n" & 
                                    myvar == "agb" & 
                                    !is.na(gs_mat)))

nadd_agb_clim_summary <- data.frame(trait = "agb",
                                     nut_add = "n",
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(nadd_agb_clim)),
                                     row.names = NULL)

# P addition
padd_agb_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "p" & 
                                    myvar == "agb" & 
                                    !is.na(gs_mat)))

padd_agb_clim_summary <- data.frame(trait = "agb",
                                     nut_add = "p",
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(padd_agb_clim)),
                                     row.names = NULL)

# N+P addition
npadd_agb_clim <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ gs_mat + gs_ai + gs_par,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "agb" & 
                                     !is.na(gs_mat)))

npadd_agb_clim_summary <- data.frame(trait = "agb",
                                      nut_add = "np",
                                      mod = c("intrcpt", "gs_mat",
                                              "gs_ai", "gs_par"),
                                      coef(summary(npadd_agb_clim)),
                                      row.names = NULL)

## Merge Parea moderator results, with some light cleaning
agb_clim_summary <- rbind(nadd_agb_clim_summary, 
                           padd_agb_clim_summary, 
                           npadd_agb_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 3)))

##############################################################################
# Belowground biomass climate moderators
##############################################################################

# N addition
nadd_bgb_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "n" & 
                                   myvar == "bgb" & 
                                   !is.na(gs_mat)))

nadd_bgb_clim_summary <- data.frame(trait = "bgb",
                                    nut_add = "n",
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(nadd_bgb_clim)),
                                    row.names = NULL)

# P addition
padd_bgb_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "p" & 
                                   myvar == "bgb" & 
                                   !is.na(gs_mat)))

padd_bgb_clim_summary <- data.frame(trait = "bgb",
                                    nut_add = "p",
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(padd_bgb_clim)),
                                    row.names = NULL)

# N+P addition
npadd_bgb_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "np" & 
                                    myvar == "bgb" & 
                                    !is.na(gs_mat)))

npadd_bgb_clim_summary <- data.frame(trait = "bgb",
                                     nut_add = "np",
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(npadd_bgb_clim)),
                                     row.names = NULL)

## Merge Parea moderator results, with some light cleaning
bgb_clim_summary <- rbind(nadd_bgb_clim_summary, 
                          padd_bgb_clim_summary, 
                          npadd_bgb_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 3)))

##############################################################################
# Root mass fraction climate moderators
##############################################################################

# N addition
nadd_rmf_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "n" & 
                                   myvar == "rmf" & 
                                   !is.na(gs_mat)))

nadd_rmf_clim_summary <- data.frame(trait = "rmf",
                                    nut_add = "n",
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(nadd_rmf_clim)),
                                    row.names = NULL)

# P addition
padd_rmf_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "p" & 
                                   myvar == "rmf" & 
                                   !is.na(gs_mat)))

padd_rmf_clim_summary <- data.frame(trait = "rmf",
                                    nut_add = "p",
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(padd_rmf_clim)),
                                    row.names = NULL)

# N+P addition
npadd_rmf_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "np" & 
                                    myvar == "rmf" & 
                                    !is.na(gs_mat)))

npadd_rmf_clim_summary <- data.frame(trait = "rmf",
                                     nut_add = "np",
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(npadd_rmf_clim)),
                                     row.names = NULL)

## Merge Parea moderator results, with some light cleaning
rmf_clim_summary <- rbind(nadd_rmf_clim_summary, 
                          padd_rmf_clim_summary, 
                          npadd_rmf_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 3)))

##############################################################################
# Root:shoot climate moderators
##############################################################################

# N addition
nadd_rootshoot_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "n" & 
                                   myvar == "rootshoot" & 
                                   !is.na(gs_mat)))

nadd_rootshoot_clim_summary <- data.frame(trait = "rootshoot",
                                    nut_add = "n",
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(nadd_rootshoot_clim)),
                                    row.names = NULL)

# P addition
padd_rootshoot_clim <- rma.mv(logr, 
                        logr_var,
                        method = "REML", 
                        random = ~ 1 | exp, 
                        mods = ~ gs_mat + gs_ai + gs_par,
                        slab = exp, 
                        control = list(stepadj = 0.3), 
                        data = meta_results %>% 
                          filter(manip_type == "p" & 
                                   myvar == "rootshoot" & 
                                   !is.na(gs_mat)))

padd_rootshoot_clim_summary <- data.frame(trait = "rmf",
                                    nut_add = "p",
                                    mod = c("intrcpt", "gs_mat",
                                            "gs_ai", "gs_par"),
                                    coef(summary(padd_rootshoot_clim)),
                                    row.names = NULL)

# N+P addition
npadd_rootshoot_clim <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ gs_mat + gs_ai + gs_par,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "np" & 
                                    myvar == "rootshoot" & 
                                    !is.na(gs_mat)))

npadd_rootshoot_clim_summary <- data.frame(trait = "rootshoot",
                                     nut_add = "np",
                                     mod = c("intrcpt", "gs_mat",
                                             "gs_ai", "gs_par"),
                                     coef(summary(npadd_rootshoot_clim)),
                                     row.names = NULL)

## Merge Parea moderator results, with some light cleaning
rootshoot_clim_summary <- rbind(nadd_rootshoot_clim_summary, 
                          padd_rootshoot_clim_summary, 
                          npadd_rootshoot_clim_summary) %>%
  mutate(across(estimate:se, ~ round(.x, digits = 4)),
         across(zval:ci.ub, ~ round(.x, digits = 3)))

##############################################################################
# Merge climate moderators and write to .csv
##############################################################################
marea_clim_summary %>%
  full_join(nmass_clim_summary) %>%
  full_join(narea_clim_summary) %>%
  full_join(pmass_clim_summary) %>%
  full_join(parea_clim_summary) %>%
  full_join(tbio_clim_summary) %>%
  full_join(agb_clim_summary) %>%
  full_join(bgb_clim_summary) %>%
  full_join(rmf_clim_summary) %>%
  full_join(rootshoot_clim_summary) %>%
  write.csv("../data/CNPmeta_clim_moderators.csv", row.names = F)

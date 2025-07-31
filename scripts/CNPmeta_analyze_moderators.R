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

##############################################################################
# Marea photosynthetic pathway
##############################################################################

# N addition
nadd_marea_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "lma" & 
                                     !is.na(photo_path)))

nadd_marea_photo_summary <- data.frame(trait = "marea", 
                                       nut_add = "n",
                                       mod_results(nadd_marea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_marea_photo))[2,3],
                                       p = coef(summary(nadd_marea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_marea_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "lma" & 
                                     !is.na(photo_path)))

padd_marea_photo_summary <- data.frame(trait = "marea", 
                                       nut_add = "p",
                                       mod_results(padd_marea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_marea_photo))[2,3],
                                       p = coef(summary(padd_marea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_marea_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "lma" & 
                                      !is.na(photo_path)))

npadd_marea_photo_summary <- data.frame(trait = "marea", 
                                       nut_add = "np",
                                       mod_results(npadd_marea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_marea_photo))[2,3],
                                       p = coef(summary(npadd_marea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Marea moderator results, with some light cleaning
marea_photo_summary <- rbind(nadd_marea_photo_summary, 
                            padd_marea_photo_summary, 
                            npadd_marea_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Nmass photosynthetic pathway
##############################################################################

# N addition
nadd_nmass_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "leaf_n_mass" & 
                                      !is.na(photo_path)))

nadd_nmass_photo_summary <- data.frame(trait = "nmass", 
                                       nut_add = "n",
                                       mod_results(nadd_nmass_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_nmass_photo))[2,3],
                                       p = coef(summary(nadd_nmass_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_nmass_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "leaf_n_mass" & 
                                      !is.na(photo_path)))

padd_nmass_photo_summary <- data.frame(trait = "nmass", 
                                       nut_add = "p",
                                       mod_results(padd_nmass_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_nmass_photo))[2,3],
                                       p = coef(summary(padd_nmass_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_nmass_photo <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ photo_path,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "leaf_n_mass" & 
                                       !is.na(photo_path)))

npadd_nmass_photo_summary <- data.frame(trait = "nmass", 
                                        nut_add = "np",
                                        mod_results(npadd_nmass_photo, 
                                                    mod = "photo_path", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_nmass_photo))[2,3],
                                        p = coef(summary(npadd_nmass_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Nmass moderator results, with some light cleaning
nmass_photo_summary <- rbind(nadd_nmass_photo_summary, 
                             padd_nmass_photo_summary, 
                             npadd_nmass_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Narea photosynthetic pathway
##############################################################################

# N addition
nadd_narea_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "leaf_n_area" & 
                                      !is.na(photo_path)))

nadd_narea_photo_summary <- data.frame(trait = "narea", 
                                       nut_add = "n",
                                       mod_results(nadd_narea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_narea_photo))[2,3],
                                       p = coef(summary(nadd_narea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_narea_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "leaf_n_area" & 
                                      !is.na(photo_path)))

padd_narea_photo_summary <- data.frame(trait = "narea", 
                                       nut_add = "p",
                                       mod_results(padd_narea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_narea_photo))[2,3],
                                       p = coef(summary(padd_narea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_narea_photo <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ photo_path,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "leaf_n_area" & 
                                       !is.na(photo_path)))

npadd_narea_photo_summary <- data.frame(trait = "narea", 
                                        nut_add = "np",
                                        mod_results(npadd_narea_photo, 
                                                    mod = "photo_path", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_narea_photo))[2,3],
                                        p = coef(summary(npadd_narea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Nmass moderator results, with some light cleaning
narea_photo_summary <- rbind(nadd_narea_photo_summary, 
                             padd_narea_photo_summary, 
                             npadd_narea_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Pmass photosynthetic pathway
##############################################################################

# N addition
nadd_pmass_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "leaf_p_mass" & 
                                      !is.na(photo_path)))

nadd_pmass_photo_summary <- data.frame(trait = "pmass", 
                                       nut_add = "n",
                                       mod_results(nadd_pmass_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_pmass_photo))[2,3],
                                       p = coef(summary(nadd_pmass_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_pmass_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "leaf_p_mass" & 
                                      !is.na(photo_path)))

padd_pmass_photo_summary <- data.frame(trait = "pmass", 
                                       nut_add = "p",
                                       mod_results(padd_pmass_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_pmass_photo))[2,3],
                                       p = coef(summary(padd_pmass_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_pmass_photo <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ photo_path,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "leaf_p_mass" & 
                                       !is.na(photo_path)))

npadd_pmass_photo_summary <- data.frame(trait = "pmass", 
                                        nut_add = "np",
                                        mod_results(npadd_pmass_photo, 
                                                    mod = "photo_path", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_pmass_photo))[2,3],
                                        p = coef(summary(npadd_pmass_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Pmass moderator results, with some light cleaning
pmass_photo_summary <- rbind(nadd_pmass_photo_summary, 
                             padd_pmass_photo_summary, 
                             npadd_pmass_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Parea photosynthetic pathway
##############################################################################

# N addition
nadd_parea_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "leaf_p_area" & 
                                      !is.na(photo_path)))

nadd_parea_photo_summary <- data.frame(trait = "parea", 
                                       nut_add = "n",
                                       mod_results(nadd_parea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_parea_photo))[2,3],
                                       p = coef(summary(nadd_parea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_parea_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "leaf_p_area" & 
                                      !is.na(photo_path)))

padd_parea_photo_summary <- data.frame(trait = "parea", 
                                       nut_add = "p",
                                       mod_results(padd_parea_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_parea_photo))[2,3],
                                       p = coef(summary(padd_parea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_parea_photo <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ photo_path,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "leaf_p_area" & 
                                       !is.na(photo_path)))

npadd_parea_photo_summary <- data.frame(trait = "parea", 
                                        nut_add = "np",
                                        mod_results(npadd_parea_photo, 
                                                    mod = "photo_path", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_parea_photo))[2,3],
                                        p = coef(summary(npadd_parea_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Parea moderator results, with some light cleaning
parea_photo_summary <- rbind(nadd_parea_photo_summary, 
                             padd_parea_photo_summary, 
                             npadd_parea_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Asat photosynthetic pathway
##############################################################################

# N addition
nadd_asat_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "asat" & 
                                      !is.na(photo_path)))

nadd_asat_photo_summary <- data.frame(trait = "asat", 
                                       nut_add = "n",
                                       mod_results(nadd_asat_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_asat_photo))[2,3],
                                       p = coef(summary(nadd_asat_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_asat_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "asat" & 
                                      !is.na(photo_path)))

padd_asat_photo_summary <- data.frame(trait = "asat", 
                                       nut_add = "p",
                                       mod_results(padd_asat_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_asat_photo))[2,3],
                                       p = coef(summary(padd_asat_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_asat_photo <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ photo_path,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "asat" & 
                                       !is.na(photo_path)))

npadd_asat_photo_summary <- data.frame(trait = "asat", 
                                        nut_add = "np",
                                        mod_results(npadd_asat_photo, 
                                                    mod = "photo_path", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_asat_photo))[2,3],
                                        p = coef(summary(npadd_asat_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Asat moderator results, with some light cleaning
asat_photo_summary <- rbind(nadd_asat_photo_summary, 
                             padd_asat_photo_summary, 
                             npadd_asat_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Vcmax photosynthetic pathway
##############################################################################

# N addition
nadd_vcmax_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "vcmax" & 
                                     !is.na(photo_path)))

nadd_vcmax_photo_summary <- data.frame(trait = "vcmax", 
                                      nut_add = "n",
                                      mod_results(nadd_vcmax_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_vcmax_photo))[2,3],
                                      p = coef(summary(nadd_vcmax_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_vcmax_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "vcmax" & 
                                     !is.na(photo_path)))

padd_vcmax_photo_summary <- data.frame(trait = "vcmax", 
                                      nut_add = "p",
                                      mod_results(padd_vcmax_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_vcmax_photo))[2,3],
                                      p = coef(summary(padd_vcmax_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_vcmax_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "vcmax" & 
                                      !is.na(photo_path)))

npadd_vcmax_photo_summary <- data.frame(trait = "vcmax", 
                                       nut_add = "np",
                                       mod_results(npadd_vcmax_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_vcmax_photo))[2,3],
                                       p = coef(summary(npadd_vcmax_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Vcmax moderator results, with some light cleaning
vcmax_photo_summary <- rbind(nadd_vcmax_photo_summary, 
                            padd_vcmax_photo_summary, 
                            npadd_vcmax_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Jmax photosynthetic pathway
##############################################################################

# N addition
nadd_jmax_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "n" & 
                                      myvar == "jmax" & 
                                      !is.na(photo_path)))

nadd_jmax_photo_summary <- data.frame(trait = "jmax", 
                                       nut_add = "n",
                                       mod_results(nadd_jmax_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(nadd_jmax_photo))[2,3],
                                       p = coef(summary(nadd_jmax_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_jmax_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "p" & 
                                      myvar == "jmax" & 
                                      !is.na(photo_path)))

padd_jmax_photo_summary <- data.frame(trait = "jmax", 
                                       nut_add = "p",
                                       mod_results(padd_jmax_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(padd_jmax_photo))[2,3],
                                       p = coef(summary(padd_jmax_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_jmax_photo <- rma.mv(logr, 
                            logr_var,
                            method = "REML", 
                            random = ~ 1 | exp, 
                            mods = ~ photo_path,
                            slab = exp, 
                            control = list(stepadj = 0.3), 
                            data = meta_results %>% 
                              filter(manip_type == "np" & 
                                       myvar == "jmax" & 
                                       !is.na(photo_path)))

npadd_jmax_photo_summary <- data.frame(trait = "jmax", 
                                        nut_add = "np",
                                        mod_results(npadd_jmax_photo, 
                                                    mod = "photo_path", 
                                                    group = "exp")$mod_table,
                                        z = coef(summary(npadd_jmax_photo))[2,3],
                                        p = coef(summary(npadd_jmax_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Jmax moderator results, with some light cleaning
jmax_photo_summary <- rbind(nadd_jmax_photo_summary, 
                             padd_jmax_photo_summary, 
                             npadd_jmax_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# PNUE photosynthetic pathway
##############################################################################

# N addition
nadd_pnue_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_pnue" & 
                                     !is.na(photo_path)))

nadd_pnue_photo_summary <- data.frame(trait = "pnue", 
                                      nut_add = "n",
                                      mod_results(nadd_pnue_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_pnue_photo))[2,3],
                                      p = coef(summary(nadd_pnue_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_pnue_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_pnue" & 
                                     !is.na(photo_path)))

padd_pnue_photo_summary <- data.frame(trait = "pnue", 
                                      nut_add = "p",
                                      mod_results(padd_pnue_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_pnue_photo))[2,3],
                                      p = coef(summary(padd_pnue_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_pnue_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_pnue" & 
                                      !is.na(photo_path)))

npadd_pnue_photo_summary <- data.frame(trait = "pnue", 
                                       nut_add = "np",
                                       mod_results(npadd_pnue_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_pnue_photo))[2,3],
                                       p = coef(summary(npadd_pnue_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Jmax moderator results, with some light cleaning
pnue_photo_summary <- rbind(nadd_pnue_photo_summary, 
                            padd_pnue_photo_summary, 
                            npadd_pnue_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# PPUE photosynthetic pathway
##############################################################################

# N addition
nadd_ppue_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "leaf_ppue" & 
                                     !is.na(photo_path)))

nadd_ppue_photo_summary <- data.frame(trait = "ppue", 
                                      nut_add = "n",
                                      mod_results(nadd_ppue_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_ppue_photo))[2,3],
                                      p = coef(summary(nadd_ppue_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_ppue_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "leaf_ppue" & 
                                     !is.na(photo_path)))

padd_ppue_photo_summary <- data.frame(trait = "ppue", 
                                      nut_add = "p",
                                      mod_results(padd_ppue_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_ppue_photo))[2,3],
                                      p = coef(summary(padd_ppue_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_ppue_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "leaf_ppue" & 
                                      !is.na(photo_path)))

npadd_ppue_photo_summary <- data.frame(trait = "ppue", 
                                       nut_add = "np",
                                       mod_results(npadd_ppue_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_ppue_photo))[2,3],
                                       p = coef(summary(npadd_ppue_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge Jmax moderator results, with some light cleaning
ppue_photo_summary <- rbind(nadd_ppue_photo_summary, 
                            padd_ppue_photo_summary, 
                            npadd_ppue_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Total biomass photosynthetic pathway
##############################################################################

# N addition
nadd_tbio_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "total_biomass" & 
                                     !is.na(photo_path)))

nadd_tbio_photo_summary <- data.frame(trait = "tbio", 
                                      nut_add = "n",
                                      mod_results(nadd_tbio_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_tbio_photo))[2,3],
                                      p = coef(summary(nadd_tbio_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_tbio_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "total_biomass" & 
                                     !is.na(photo_path)))

padd_tbio_photo_summary <- data.frame(trait = "tbio", 
                                      nut_add = "p",
                                      mod_results(padd_tbio_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_tbio_photo))[2,3],
                                      p = coef(summary(padd_tbio_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_tbio_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "total_biomass" & 
                                      !is.na(photo_path)))

npadd_tbio_photo_summary <- data.frame(trait = "tbio", 
                                       nut_add = "np",
                                       mod_results(npadd_tbio_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_tbio_photo))[2,3],
                                       p = coef(summary(npadd_tbio_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge total biomass moderator results, with some light cleaning
tbio_photo_summary <- rbind(nadd_tbio_photo_summary, 
                            padd_tbio_photo_summary, 
                            npadd_tbio_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Aboveground biomass photosynthetic pathway
##############################################################################

# N addition
nadd_agb_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "n" & 
                                     myvar == "agb" & 
                                     !is.na(photo_path)))

nadd_agb_photo_summary <- data.frame(trait = "agb", 
                                      nut_add = "n",
                                      mod_results(nadd_agb_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(nadd_agb_photo))[2,3],
                                      p = coef(summary(nadd_agb_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_agb_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "p" & 
                                     myvar == "agb" & 
                                     !is.na(photo_path)))

padd_agb_photo_summary <- data.frame(trait = "agb", 
                                      nut_add = "p",
                                      mod_results(padd_agb_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(padd_agb_photo))[2,3],
                                      p = coef(summary(padd_agb_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_agb_photo <- rma.mv(logr, 
                           logr_var,
                           method = "REML", 
                           random = ~ 1 | exp, 
                           mods = ~ photo_path,
                           slab = exp, 
                           control = list(stepadj = 0.3), 
                           data = meta_results %>% 
                             filter(manip_type == "np" & 
                                      myvar == "agb" & 
                                      !is.na(photo_path)))

npadd_agb_photo_summary <- data.frame(trait = "agb", 
                                       nut_add = "np",
                                       mod_results(npadd_agb_photo, 
                                                   mod = "photo_path", 
                                                   group = "exp")$mod_table,
                                       z = coef(summary(npadd_agb_photo))[2,3],
                                       p = coef(summary(npadd_agb_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge aboveground biomass moderator results, with some light cleaning
agb_photo_summary <- rbind(nadd_agb_photo_summary, 
                            padd_agb_photo_summary, 
                            npadd_agb_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))

##############################################################################
# Belowground biomass photosynthetic pathway
##############################################################################

# N addition
nadd_bgb_photo <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "n" & 
                                    myvar == "bgb" & 
                                    !is.na(photo_path)))

nadd_bgb_photo_summary <- data.frame(trait = "bgb", 
                                     nut_add = "n",
                                     mod_results(nadd_bgb_photo, 
                                                 mod = "photo_path", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(nadd_bgb_photo))[2,3],
                                     p = coef(summary(nadd_bgb_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# P addition
padd_bgb_photo <- rma.mv(logr, 
                         logr_var,
                         method = "REML", 
                         random = ~ 1 | exp, 
                         mods = ~ photo_path,
                         slab = exp, 
                         control = list(stepadj = 0.3), 
                         data = meta_results %>% 
                           filter(manip_type == "p" & 
                                    myvar == "bgb" & 
                                    !is.na(photo_path)))

padd_bgb_photo_summary <- data.frame(trait = "bgb", 
                                     nut_add = "p",
                                     mod_results(padd_bgb_photo, 
                                                 mod = "photo_path", 
                                                 group = "exp")$mod_table,
                                     z = coef(summary(padd_bgb_photo))[2,3],
                                     p = coef(summary(padd_agb_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

# N+P addition
npadd_bgb_photo <- rma.mv(logr, 
                          logr_var,
                          method = "REML", 
                          random = ~ 1 | exp, 
                          mods = ~ photo_path,
                          slab = exp, 
                          control = list(stepadj = 0.3), 
                          data = meta_results %>% 
                            filter(manip_type == "np" & 
                                     myvar == "bgb" & 
                                     !is.na(photo_path)))

npadd_bgb_photo_summary <- data.frame(trait = "bgb", 
                                      nut_add = "np",
                                      mod_results(npadd_bgb_photo, 
                                                  mod = "photo_path", 
                                                  group = "exp")$mod_table,
                                      z = coef(summary(npadd_bgb_photo))[2,3],
                                      p = coef(summary(npadd_bgb_photo))[2, 4]) %>%
  dplyr::select(trait, nut_add, photo = name, estimate, z, p, lowerCL, upperCL)

## Merge aboveground biomass moderator results, with some light cleaning
bgb_photo_summary <- rbind(nadd_bgb_photo_summary, 
                           padd_bgb_photo_summary, 
                           npadd_bgb_photo_summary) %>%
  mutate(estimate = round(estimate, digits = 3),
         across(z:upperCL, ~ round(.x, digits = 3)),
         p = as.character(ifelse(p < 0.001, "<0.001", p)))


##############################################################################
# Merge photo moderator results and write to .csv
##############################################################################
marea_photo_summary %>%
  full_join(nmass_photo_summary) %>%
  full_join(narea_photo_summary) %>%
  full_join(pmass_photo_summary) %>%
  full_join(parea_photo_summary) %>%
  full_join(asat_photo_summary) %>%
  full_join(vcmax_photo_summary) %>%
  full_join(jmax_photo_summary) %>%
  full_join(pnue_photo_summary) %>%
  full_join(ppue_photo_summary) %>%
  full_join(tbio_photo_summary) %>%
  full_join(agb_photo_summary) %>%
  full_join(bgb_photo_summary) %>%
  write.csv("../data/CNPmeta_photo_moderators.csv", row.names = F)

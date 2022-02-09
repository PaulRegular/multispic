
source("analysis/001_NL_case_study_helpers.R")


list2env(nl_inputs_and_priors(region = "2J3K", species = NULL), envir = globalenv())

inputs$landings$spp_group <- ifelse(inputs$landings$species %in% c("Atlantic Cod-2J3K"),
                                    "cod", "not-cod")

inputs$index$spp_group <- ifelse(inputs$index$species %in% c("Atlantic Cod-2K3K"),
                                 "cod", "not-cod")

fit <- multispic(inputs, species_cor = "none", temporal_cor = "AR1",
                 log_K_option = par_option(option = "normal_prior",
                                           mean = mean_log_K, sd = sd_log_K),
                 log_B0_option = par_option(option = "normal_prior",
                                            mean = mean_log_B0, sd = sd_log_B0),
                 log_r_option = par_option(option = "normal_prior",
                                           mean = mean_log_r, sd = sd_log_r),
                 log_sd_B_option = par_option(option = "normal_prior",
                                              mean = mean_log_sd_B, sd = sd_log_sd_B),
                 log_q_option = par_option(option = "normal_prior",
                                           mean = mean_log_q, sd = sd_log_q),
                 log_sd_I_option = par_option(option = "normal_prior",
                                              mean = mean_log_sd_I, sd = sd_log_sd_I),
                 logit_rho_option = par_option(option = "normal_prior",
                                               mean = mean_logit_rho, sd = sd_logit_rho),
                 logit_phi_option = par_option(option = "normal_prior",
                                               mean = mean_logit_phi, sd = sd_logit_phi),
                 n_forecast = 1, K_groups = ~spp_group, survey_groups = ~gear + species,
                 pe_covariates = ~0)

vis_multispic(fit)




source("analysis/001_NL_case_study_helpers.R")


list2env(nl_inputs_and_priors(region = "3LNO",
                              species = c("Atlantic Cod", "American Plaice", "Redfish spp.")), envir = globalenv())


fit <- multispic(inputs, species_cor = "all", temporal_cor = "ar1",
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
                     n_forecast = 1, K_groups = ~1, survey_groups = ~gear + season * species,
                     pe_covariates = ~0)

vis_multispic(fit)



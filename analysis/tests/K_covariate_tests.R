
source("analysis/001_NL_case_study_helpers.R")

list2env(nl_inputs_and_priors(region = "2J3K", species = NULL), envir = globalenv())

inputs$landings$shift <- ifelse(inputs$landings$year < 1991, "pre-1991", "post-1991")
inputs$landings$is_cod <- ifelse(inputs$landings$species == "Atlantic Cod-2J3K", "cod", "not-cod")
inputs$index$is_cod <- ifelse(inputs$index$species == "Atlantic Cod-2J3K", "cod", "not-cod")

fit <- multispic(inputs, species_cor = "one", temporal_cor = "ar1",
                 log_K_option = par_option(option = "prior",
                                           mean = c(mean_log_K, mean_log_K + 2), sd = rep(sd_log_K, 2)),
                 log_B0_option = par_option(option = "prior",
                                            mean = mean_log_B0, sd = sd_log_B0),
                 log_r_option = par_option(option = "prior",
                                           mean = mean_log_r, sd = sd_log_r),
                 log_sd_B_option = par_option(option = "prior",
                                              mean = mean_log_sd_B, sd = sd_log_sd_B),
                 log_q_option = par_option(option = "prior",
                                           mean = mean_log_q, sd = sd_log_q),
                 log_sd_I_option = par_option(option = "prior",
                                              mean = mean_log_sd_I, sd = sd_log_sd_I),
                 logit_rho_option = par_option(option = "prior",
                                               mean = mean_logit_rho, sd = sd_logit_rho),
                 logit_phi_option = par_option(option = "prior",
                                               mean = mean_logit_phi, sd = sd_logit_phi),
                 K_betas_option = par_option(option = "prior",
                                             mean = mean_K_betas, sd = sd_K_betas),
                 n_forecast = 1, K_groups = ~is_cod, survey_groups = ~species_survey,
                 pe_covariates = ~0, K_covariates = ~shift)

vis_multispic(fit)




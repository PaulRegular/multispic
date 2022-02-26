

source("analysis/001_NL_case_study_helpers.R")

tot_landings <- multispic::landings %>%
    filter(region == "3LNO") %>%
    group_by(species) %>%
    summarise(total_landings = sum(landings)) %>%
    arrange(-total_landings)

spp <- head(tot_landings$species, 7)

## 3LNO does not work with Hake in the mix

list2env(nl_inputs_and_priors(region = "3LNO", species = spp,
                              K_groups = ~species), envir = globalenv())

sp_fit <- multispic(inputs, species_cor = "none", temporal_cor = "none",
                    log_K_option = par_option(option = "prior",
                                              mean = mean_log_K, sd = sd_log_K),
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
                    pe_betas_option = par_option(option = "prior",
                                                 mean = mean_pe_betas, sd = sd_pe_betas),
                    n_forecast = 1, K_groups = ~species, survey_groups = ~species_survey,
                    pe_covariates = ~0, K_covariates = ~0)

vis_multispic(sp_fit)



tot_landings <- multispic::landings %>%
    filter(region == "3Ps") %>%
    group_by(species) %>%
    summarise(total_landings = sum(landings)) %>%
    arrange(-total_landings)

spp <- head(tot_landings$species, 3)

spp <- c("American Plaice", "Greenland Halibut",
         "Skate spp.", "Witch Flounder", "Haddock",
         "Monkfish", "White Hake", "Atlantic Halibut",
         "Yellowtail Flounder", "Silver Hake")

spp <- c("American Plaice", "Haddock", "Witch Flounder", "White Hake",
         "Skate spp.", "Yellowtail Flounder", "Atlantic Halibut",
         "Wolffish spp.")

list2env(nl_inputs_and_priors(region = "3Ps", species = spp,
                              K_groups = ~species), envir = globalenv())

sp_fit <- multispic(inputs, species_cor = "none", temporal_cor = "none",
                    log_K_option = par_option(option = "prior",
                                              mean = mean_log_K, sd = sd_log_K),
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
                    pe_betas_option = par_option(option = "prior",
                                                 mean = mean_pe_betas, sd = sd_pe_betas),
                    n_forecast = 1, K_groups = ~species, survey_groups = ~species_survey,
                    pe_covariates = ~0, K_covariates = ~0)
sp_fit$sd_rep

vis_multispic(sp_fit)


list2env(nl_inputs_and_priors(region = "3Ps", species = "Atlantic Cod",
                              K_groups = ~species), envir = globalenv())

sp_fit <- multispic(inputs, species_cor = "none", temporal_cor = "none",
                    log_K_option = par_option(option = "prior",
                                              mean = mean_log_K, sd = sd_log_K),
                    log_B0_option = par_option(option = "prior",
                                               mean = mean_log_B0, sd = sd_log_B0),
                    log_r_option = par_option(option = "prior",
                                              mean = mean_log_r, sd = sd_log_r),
                    log_sd_B_option = par_option(option = "prior",
                                                 mean = mean_log_sd_B, sd = sd_log_sd_B),
                    log_q_option = par_option(option = "prior",
                                              mean = mean_log_q, sd = sd_log_q * 5),
                    log_sd_I_option = par_option(option = "prior",
                                                 mean = mean_log_sd_I, sd = sd_log_sd_I),
                    logit_rho_option = par_option(option = "prior",
                                                  mean = mean_logit_rho, sd = sd_logit_rho),
                    logit_phi_option = par_option(option = "prior",
                                                  mean = mean_logit_phi, sd = sd_logit_phi),
                    K_betas_option = par_option(option = "prior",
                                                mean = mean_K_betas, sd = sd_K_betas),
                    pe_betas_option = par_option(option = "prior",
                                                 mean = mean_pe_betas, sd = sd_pe_betas),
                    n_forecast = 1, K_groups = ~1, survey_groups = ~species_survey,
                    pe_covariates = ~0, K_covariates = ~0)


sp_fit$sd_rep



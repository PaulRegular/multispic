
source("analysis/001_NL_case_study_helpers.R")


list2env(nl_inputs_and_priors(region = "3LNO", species = NULL), envir = globalenv())

inputs$landings <- inputs$landings %>%
    mutate(three_groups = recode(species,
                                 "Atlantic Cod-3LNO" = "gadid",
                                 "Haddock-3LNO" = "gadid",
                                 "American Plaice-3LNO" = "flat",
                                 "Atlantic Halibut-3LNO" = "flat",
                                 "Greenland Halibut-3LNO" = "flat",
                                 "Yellowtail Flounder-3LNO" = "flat",
                                 "Witch Flounder-3LNO" = "flat",
                                 "Wolffish spp.-3LNO" = "other",
                                 "Redfish spp.-3LNO" = "other",
                                 "Skate spp.-3LNO" = "other"))
inputs$index <- inputs$index %>%
    mutate(three_groups = recode(species,
                                 "Atlantic Cod-3LNO" = "gadid",
                                 "Haddock-3LNO" = "gadid",
                                 "American Plaice-3LNO" = "flat",
                                 "Atlantic Halibut-3LNO" = "flat",
                                 "Greenland Halibut-3LNO" = "flat",
                                 "Yellowtail Flounder-3LNO" = "flat",
                                 "Witch Flounder-3LNO" = "flat",
                                 "Wolffish spp.-3LNO" = "other",
                                 "Redfish spp.-3LNO" = "other",
                                 "Skate spp.-3LNO" = "other"))


fit <- multispic(inputs, species_cor = "one", temporal_cor = "ar1",
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
                 n_forecast = 1, K_groups = ~three_groups, survey_groups = ~species_survey,
                 pe_covariates = ~0)

vis_multispic(fit)



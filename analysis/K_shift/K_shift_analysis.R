
source("analysis/001_NL_case_study_helpers.R")

for (r in c("2J3K", "3LNO", "3Ps")) {

    list2env(nl_inputs_and_priors(region = r, species = NULL), envir = globalenv())

    inputs$landings$shift <- ifelse(inputs$landings$year < 1991, "pre-1991", "post-1991")
    inputs$landings$shift <- factor(inputs$landings$shift, levels = c("pre-1991", "post-1991"))

    one_K <- multispic(inputs, species_cor = "one", temporal_cor = "ar1",
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
                       n_forecast = 1, K_groups = ~1, survey_groups = ~species_survey,
                       pe_covariates = ~0, K_covariates = ~0)

    two_K <- update(one_K, K_covariates = ~shift)

    vis_multispic(one_K, output_file = paste0("analysis/K_shift/one_K_", r, ".html"))

    vis_multispic(two_K, output_file = paste0("analysis/K_shift/two_K_", r, ".html"))

}



library(multispic)
library(plotly)
library(dplyr)

## Capelin data ------------------------------------------------------------------------------------

index <- read.csv("data-raw/covariates/capelin/capelin_biomass.csv") %>%
    filter(survey != "Trinity Bay (3L)",  # too small of an area
           !is.na(biomass)) %>%
    rename(index = biomass)

index %>%
    plot_ly(x = ~year, y = ~index, color = ~survey,
            colors = viridis::viridis(100)) %>%
    add_lines()

landings <- read.csv("data-raw/landings/STATLANT21A_Extraction.csv")
names(landings) <- c("year", "country", "division", "species", "landings")

landings <- landings %>%
    filter(species == "CAPELIN - CAP",
           division %in% c("2J", "3K", "3L", "3N", "3O")) %>%
    group_by(year) %>%
    summarise(landings = sum(landings) / 1000) %>%
    ungroup() %>%
    as.data.frame()

landings %>%
    plot_ly(x = ~year, y = ~landings) %>%
    add_lines()


## Set-up priors -----------------------------------------------------------------------------------

tot_landings <- landings %>%
    group_by(year) %>%
    summarise(tot_landings = sum(landings))

max_tot_landings <- max(tot_landings$tot_landings)

lower_log_r <- log(0.01)
upper_log_r <- log(1)
mean_log_r <- (lower_log_r + upper_log_r) / 2
sd_log_r <- (upper_log_r - lower_log_r) / 2

lower_log_K <- log(max_tot_landings) - upper_log_r
upper_log_K <- log(max_tot_landings * 100) - lower_log_r
mean_log_K <- (lower_log_K + upper_log_K) / 2
sd_log_K <- (upper_log_K - lower_log_K) / 2

L0 <- landings$landings[landings$year == min(index$year)]
lower_log_B0 <- log(L0) - upper_log_r
upper_log_B0 <- log(L0 * 100) - lower_log_r
mean_log_B0 <- (lower_log_B0 + upper_log_B0) / 2
sd_log_B0 <- c(upper_log_B0 - lower_log_B0) / 2

lower_log_sd_B <- log(0.01)
upper_log_sd_B <- log(1)
mean_log_sd_B <- (lower_log_sd_B + upper_log_sd_B) / 2
sd_log_sd_B <- (upper_log_sd_B - lower_log_sd_B) / 2

lower_log_sd_I <- log(0.01)
upper_log_sd_I <- log(1)
mean_log_sd_I <- (lower_log_sd_I + upper_log_sd_I) / 2
sd_log_sd_I <- (upper_log_sd_I - lower_log_sd_I) / 2

lower_log_q <- log(0.5)
upper_log_q <- log(1)
mean_log_q <- (lower_log_q + upper_log_q) / 2
sd_log_q <- (upper_log_q - lower_log_q) / 2

lower_logit_rho <- logit(-0.9, shift = TRUE)
upper_logit_rho <- logit(0.9, shift = TRUE)
mean_logit_rho <- (lower_logit_rho + upper_logit_rho) / 2
sd_logit_rho <- (upper_logit_rho - lower_logit_rho) / 2

lower_logit_phi <- logit(0.1)
upper_logit_phi <- logit(0.9)
mean_logit_phi <- (lower_logit_phi + upper_logit_phi) / 2
sd_logit_phi <- (upper_logit_phi - lower_logit_phi) / 2

mean_pe_betas <- 0
sd_pe_betas <- 10

mean_K_betas <- 0
sd_K_betas <- 10

## Use full landings time series to inform prior for K but limit analysis to
## index time series as there is nothing to inform the process errors prior to the start
## of the surveys
landings <- landings[landings$year >= min(index$year) &
                         landings$year <= max(index$year), ]

index$species <- landings$species <- "Capelin"

inputs <- list(landings = landings, index = index)


## Run model ---------------------------------------------------------------------------------------

fit <- multispic(inputs, species_cor = "none", temporal_cor = "ar1",
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
                 n_forecast = 1, K_groups = ~1, survey_groups = ~survey,
                 pe_covariates = ~0, K_covariates = ~0)
fit$sd_rep

vis_multispic(fit)



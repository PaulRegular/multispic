
## TODO:
## - Calculate one-step ahead residuals

## - Add data from 2J3K and 3Ps eco-regions <-- done
## - Fit to each region (single-species, 5 species, 10 species, all species)
## - Try to combine regions (K_group = ~region)
## - Calculate species-specific leave one out scores to assess predictive ability of each model
## - Hypothesis: multispecies >> single-species inference


library(units)
library(plotly)
library(TMB)
library(multispic)
library(dplyr)
library(zoo)


## All regions ------------------------------------------------------------------------------

index <- multispic::index %>% filter(region == "3LNO")
landings <- multispic::landings %>% filter(region == "3LNO")
covariates <- multispic::covariates
landings <- merge(landings, covariates, by = "year", all.x = TRUE)

## Limit analysis to start of survey series and to species with indices
sp_region <- table(index$species, index$region) > 0 # present-absent
sp_ind <- rowSums(sp_region) == max(rowSums(sp_region))
sub_sp <- rownames(sp_region)[sp_ind]

sub_sp <- sub_sp[sub_sp != "Silver Hake"]
# sub_sp <- c("American Plaice", "Atlantic Cod", "Greenland Halibut", "Redfish spp.")

index <- index[index$species %in% sub_sp, ]
landings <- landings[landings$year >= min(index$year) &
                         landings$year <= max(index$year) &
                         landings$species %in% sub_sp, ]

## Survey = species-region-gear-season (i.e. groups for catchability estimates)
index$survey <- paste0(index$species, "-", index$region, "-", index$season, "-", index$gear)

## Species (stock) = species-region
index$species <- paste0(index$species, "-", index$region)
landings$species <- paste0(landings$species, "-", landings$region)

p <- index %>%
    plot_ly() %>%
    add_trace(x = ~year, y = ~index, color = ~survey,
              colors = viridis::viridis(100), mode = "markers+lines",
              type = "scatter")
p
p %>% layout(yaxis = list(type = "log"))


p <- landings %>%
    plot_ly() %>%
    add_trace(x = ~year, y = ~landings, color = ~species,
              colors = viridis::viridis(100), mode = "markers+lines",
              type = "scatter")
p
p %>% layout(yaxis = list(type = "log"))

landings %>%
    group_by(year) %>%
    summarise(total_landings = sum(landings)) %>%
    ungroup() %>%
    plot_ly() %>%
    add_lines(x = ~year, y = ~total_landings)


landings %>%
    plot_ly(x = ~year) %>%
    add_lines(y = ~winter_nao, name = "Winter NAO") %>%
    add_lines(y = ~spring_nao, name = "Spring NAO") %>%
    add_lines(y = ~summer_nao, name = "Summer NAO") %>%
    add_lines(y = ~fall_nao, name = "Fall NAO")


inputs <- list(landings = landings, index = index)


## Run model -------------------------------------------------------------------

## Set-up prior settings
## Note: during testing the fixed, random, and uniform_prior options rarely converged


## Find the smallest max aggregate landings among the groups to inform lower range for K

tot_landings <- landings %>%
    group_by(year, region) %>%
    summarise(tot_landings = sum(landings))

max_tot_landings <- tot_landings %>%
    group_by(region) %>%
    summarise(max = max(tot_landings)) %>%
    with(min(max))

lower_log_r <- log(0.01)
upper_log_r <- log(1)
mean_log_r <- (lower_log_r + upper_log_r) / 2
sd_log_r <- (upper_log_r - lower_log_r) / 2

lower_log_K <- log(max_tot_landings) - upper_log_r
upper_log_K <- log(max_tot_landings * 100) - lower_log_r
mean_log_K <- (lower_log_K + upper_log_K) / 2
sd_log_K <- (upper_log_K - lower_log_K) / 2

lower_log_B0 <- log(0.01 * exp(lower_log_K) / length(unique(landings$species)))
upper_log_B0 <- upper_log_K
mean_log_B0 <- (lower_log_B0 + upper_log_B0) / 2
sd_log_B0 <- c(upper_log_B0 - lower_log_B0) / 2

lower_log_sd <- log(0.01)
upper_log_sd <- log(1)
mean_log_sd <- (lower_log_sd + upper_log_sd) / 2
sd_log_sd <- (upper_log_sd - lower_log_sd) / 2

lower_log_q <- log(0.1)
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


## Multivariate AR1 process now working
## Forcing the RW structure results in unusual process errors for some species
fit <- multispic(inputs, species_cor = "all", temporal_cor = "ar1",
                 log_K_option = par_option(option = "normal_prior",
                                           mean = mean_log_K, sd = sd_log_K),
                 log_B0_option = par_option(option = "normal_prior",
                                            mean = mean_log_B0, sd = sd_log_B0),
                 log_r_option = par_option(option = "normal_prior",
                                           mean = mean_log_r, sd = sd_log_r),
                 log_sd_B_option = par_option(option = "normal_prior",
                                              mean = mean_log_sd, sd = sd_log_sd),
                 log_q_option = par_option(option = "normal_prior",
                                           mean = mean_log_q, sd = sd_log_q),
                 log_sd_I_option = par_option(option = "normal_prior",
                                              mean = mean_log_sd, sd = sd_log_sd),
                 logit_rho_option = par_option(option = "normal_prior",
                                               mean = mean_logit_rho, sd = sd_logit_rho),
                 logit_phi_option = par_option(option = "normal_prior",
                                               mean = mean_logit_phi, sd = sd_logit_phi),
                 n_forecast = 1, K_groups = NULL, pe_covariates = NULL)

fit$opt$message
fit$sd_rep
fit$opt$objective
fit$mAIC

# loo_fit <- run_loo(fit)
# loo_fit$mse

if (is.null(fit$call$K_groups)) {
    K_groups <- NULL
    K_label <- "All species"
} else {
    K_groups <- as.formula(fit$call$K_groups)
    K_label <- unique(fit$landings[, all.vars(K_groups)])
}


## Raw par
par <- as.list(fit$sd_rep, "Est")
hist(unlist(par), breaks = 30)

## Prior and posterior
post_mean <- as.list(fit$sd_rep, "Est")
post_sd <- as.list(fit$sd_rep, "Std. Error")
plot_prior_post(prior_mean = mean_log_K, prior_sd = sd_log_K,
                post_mean = post_mean$log_K,
                post_sd = post_sd$log_K,
                post_names = K_label,
                xlab = "log(K)")
plot_prior_post(prior_mean = mean_log_r, prior_sd = sd_log_r,
                post_mean = post_mean$log_r,
                post_sd = post_sd$log_r,
                post_names = levels(fit$landings$species),
                xlab = "log(r)")
plot_prior_post(prior_mean = mean_log_B0, prior_sd = sd_log_B0,
                post_mean = post_mean$log_B0,
                post_sd = post_sd$log_B0,
                post_names = levels(fit$landings$species),
                xlab = "log(B0)")
plot_prior_post(prior_mean = mean_log_sd, prior_sd = sd_log_sd,
                post_mean = post_mean$log_sd_B,
                post_sd = post_sd$log_sd_B,
                post_names = levels(fit$landings$species),
                xlab = "log(SD<sub>B</sub>)")
plot_prior_post(prior_mean = mean_log_q, prior_sd = sd_log_q,
                post_mean = post_mean$log_q,
                post_sd = post_sd$log_q,
                post_names = levels(fit$index$survey),
                xlab = "log(q)")
plot_prior_post(prior_mean = mean_log_sd, prior_sd = sd_log_sd,
                post_mean = post_mean$log_sd_I,
                post_sd = post_sd$log_sd_I,
                post_names = levels(fit$index$survey),
                xlab = "log(SD<sub>I</sub>)")

sp_rho <- sp_nm_mat <- matrix(NA, nrow = nlevels(fit$pop$species), ncol = nlevels(fit$pop$species))
rownames(sp_rho) <- colnames(sp_rho) <- levels(fit$pop$species)
sp_rho[lower.tri(sp_rho)] <- inv_logit(post_mean$logit_rho, shift = TRUE)
sp_rho <- t(sp_rho)
sp_rho[lower.tri(sp_rho)] <- inv_logit(post_mean$logit_rho, shift = TRUE)
diag(sp_rho) <- 1
round(sp_rho, 2)
for (i in seq(nrow(sp_rho))) {
    for (j in seq(ncol(sp_rho))) {
        sp_nm_mat[i, j] <- paste(levels(fit$pop$species)[i], "-", levels(fit$pop$species)[j])
    }
}

plot_prior_post(prior_mean = mean_logit_rho, prior_sd = sd_logit_rho,
                post_mean = post_mean$logit_rho,
                post_sd = post_sd$logit_rho,
                post_names = sp_nm_mat[lower.tri(sp_nm_mat)],
                xlab = "logit(rho)")# , trans_fun = function(x) inv_logit(x, shift = TRUE))

plot_prior_post(prior_mean = mean_logit_phi, prior_sd = sd_logit_phi,
                post_mean = post_mean$logit_phi,
                post_sd = post_sd$logit_phi,
                post_names = "logit(phi)",
                xlab = "logit(phi)")# , trans_fun = inv_logit)


## Visually assess par
par <- fit$par
q <- exp(par$log_q)
names(q) <- levels(fit$index$survey)
round(q, 2)
sd_I <- exp(par$log_sd_I)
names(sd_I) <- levels(fit$index$survey)
round(sd_I, 2)
K <- exp(par$log_K)
names(K) <- K_label
signif(K, 2)
r <- exp(par$log_r)
names(r) <- levels(fit$index$species)
round(r, 2)
sd_B <- exp(par$log_sd_B)
names(sd_B) <- levels(fit$index$species)
round(sd_B, 2)
B0 <- exp(par$log_B0)
round(B0)
rho <- inv_logit(par$logit_rho, shift = TRUE)
round(rho, 2)
round(sp_rho, 2)
phi <- inv_logit(par$logit_phi)
round(phi, 2)
fit$par$pe_betas; fit$par_lwr$pe_betas; fit$par_upr$pe_betas

## Explore parameter correlations
sd_rep <- fit$sd_rep
dsd <- sqrt(diag(sd_rep$cov.fixed))
cor_mat <- diag(1 / dsd) %*% sd_rep$cov.fixed %*% diag(1 / dsd)
rownames(cor_mat) <- paste0(seq(nrow(cor_mat)), ":", names(dsd))
colnames(cor_mat) <- paste0(seq(nrow(cor_mat)), ":", names(dsd))
cor_tab <- as.data.frame.table(cor_mat)
names(cor_tab) <- c("x", "y", "z")
# corrplot::corrplot.mixed(cor_mat, diag = "n", lower = "ellipse", upper = "number")
cor_tab %>%
    plot_ly(x = ~x, y = ~y, z = ~z, text = ~text,
            colors = c("#B2182B", "white", "#2166AC")) %>%
    add_heatmap() %>% colorbar(limits = c(-1, 1)) %>%
    layout(xaxis = list(title = "", showticklabels = TRUE),
           yaxis = list(title = "", showticklabels = TRUE))

## Index residuals
p <- fit$index %>%
    plot_ly(color = ~species, colors = viridis::viridis(100))
p %>% add_markers(x = ~year, y = ~std_res)
p %>% add_markers(x = ~log(pred), y = ~std_res)
p %>% add_markers(x = ~survey, y = ~std_res)
p %>% add_markers(x = ~species, y = ~std_res)

fit$index %>%
    plot_ly(x = ~year, y = ~std_res, color = ~survey,
            colors = viridis::viridis(100)) %>%
    add_lines()


## Process error residuals
p <- fit$pop %>%
    plot_ly(color = ~species, colors = viridis::viridis(100))
p %>% add_lines(x = ~year, y = ~pe)

## pe vs covariates
fit$pop %>%
    plot_ly(color = ~species, colors = viridis::viridis(100)) %>%
    add_markers(x = ~winter_nao, y = ~pe, text = ~year)

## Correlation in pe - post hoc
pe_wide <- tidyr::spread(fit$pop[, c("year", "species", "pe")], species, pe)
pe_wide$year <- NULL
plot(pe_wide)
pe_wide <- as.matrix(pe_wide)
pe_wide[is.infinite(pe_wide)] <- NA
cor_mat <- cor(pe_wide, use = "na.or.complete")
plot_ly(x = rownames(cor_mat), y = rownames(cor_mat), z = ~cor_mat,
        colors = c("#B2182B", "white", "#2166AC")) %>%
    add_heatmap() %>% colorbar(limits = c(-1, 1))
corrplot::corrplot.mixed(cor_mat, diag = "n", lower = "ellipse", upper = "number")


## Correlation in pe - parameter estimates
plot_ly(x = rownames(sp_rho), y = rownames(sp_rho), z = ~sp_rho,
        colors = c("#B2182B", "white", "#2166AC")) %>%
    add_heatmap() %>% colorbar(limits = c(-1, 1), title = "ρ")
corrplot::corrplot.mixed(sp_rho, diag = "n", lower = "ellipse", upper = "number")

## Fits to the index
p <- fit$index %>%
    plot_ly(x = ~year, color = ~survey, colors = viridis::viridis(100),
            legendgroup = ~survey) %>%
    add_ribbons(ymin = ~pred_lwr, ymax = ~pred_upr, line = list(width = 0),
                alpha = 0.2, showlegend = FALSE) %>%
    add_lines(y = ~pred) %>%
    add_markers(y = ~index, showlegend = FALSE)
p
p %>% layout(yaxis = list(type = "log"))


## Biomass
p <- fit$pop %>%
    plot_ly(x = ~year, color = ~species, colors = viridis::viridis(100),
            legendgroup = ~species) %>%
    add_ribbons(ymin = ~B_lwr, ymax = ~B_upr, line = list(width = 0),
                alpha = 0.2, showlegend = FALSE) %>%
    add_lines(y = ~B)
# if (!is.null(K_groups)) {
#     p <- p %>% add_lines(x = ~year, y = ~K, linetype = I(3), name = "K")
# } else {
#     p <- p %>% add_lines(x = ~year, y = ~K, linetype = I(3), name = "K",
#                          color = I("black"), legendgroup = NULL)
# }

p
p %>% layout(yaxis = list(type = "log"))


## Total biomass

p <- fit$tot_pop %>%
    plot_ly(x = ~year, frame = K_groups) %>%
    add_ribbons(ymin = ~B_lwr, ymax = ~B_upr, line = list(width = 0),
                alpha = 0.2, showlegend = FALSE, legendgroup = "B",
                color = I("steelblue")) %>%
    add_lines(y = ~B, name = "B", color = I("steelblue"),
              legendgroup = "B") %>%
    add_lines(y = ~K, name = "K", legendgroup = "K",
              linetype = I(1), color = I("black")) %>%
    add_lines(y = ~K_lwr, legendgroup = "K", showlegend = FALSE,
              linetype = I(3), color = I("black"), size = I(1)) %>%
    add_lines(y = ~K_upr, legendgroup = "K", showlegend = FALSE,
              linetype = I(3), color = I("black"), size = I(1)) %>%
    layout(title = "Total biomass")
p
p %>% layout(yaxis = list(type = "log"))


## F
p <- fit$pop %>%
    # filter(is.finite(F)) %>%
    plot_ly(x = ~year, color = ~species, colors = viridis::viridis(100),
            legendgroup = ~species) %>%
    add_ribbons(ymin = ~F_lwr, ymax = ~F_upr, line = list(width = 0),
                alpha = 0.2, showlegend = FALSE) %>%
    add_lines(y = ~F)
p

## K
p <- fit$pop %>%
    plot_ly(x = ~year, color = ~species, legendgroup = ~species) %>%
    add_ribbons(ymin = ~K_lwr, ymax = ~K_upr, line = list(width = 0),
                alpha = 0.2, showlegend = FALSE, name = "95% CI") %>%
    add_lines(y = ~K)
p


## Compare to accepted assessment model results --------------------------------

## Note: there is a degree of spatial missmatch here as some of these assessments are 3NO

assess <- read.csv("analysis/stock_assessment_estimates.csv")
names(assess) <- c("species_div", "year", "B", "B_type")
x <- data.table::tstrsplit(assess$species, split = " ")
assess$species <- x[[1]]
unique(assess$species)
assess$species[assess$species == "Cod"] <- "Atlantic Cod-3LNO"
assess$species[assess$species == "Plaice"] <- "American Plaice-3LNO"
assess$species[assess$species == "Redfish"] <- "Redfish spp.-3LNO"
assess$species[assess$species == "Yellowtail"] <- "Yellowtail Flounder-3LNO"
assess$species[assess$species == "Witch"] <- "Witch Flounder-3LNO"
assess$division <- x[[2]]
assess$source <- "assessment"
assess$B_lwr <- NA
assess$B_upr <- NA

spm <- fit$pop[fit$pop$species %in% unique(assess$species), ]
spm$source <- "multispic"

keep <- c("year", "species", "B", "B_lwr", "B_upr", "source")
comp <- rbind(assess[, keep], spm[, keep])
comp <- comp[order(comp$species, comp$source), ]

comp %>%
    group_by(species, source) %>%
    mutate(scaled_B = c(scale(B)),
           center = attr(scale(B), "scaled:center"),
           scale = attr(scale(B), "scaled:scale")) %>%
    mutate(lwr = (B_lwr - center) / scale,
           upr = (B_upr - center) / scale) %>%
    filter(year > min(index$year)) %>%
    ungroup() %>%
    plot_ly(x = ~year, color = ~species, colors = viridis::viridis(100),
            linetype = ~source, legendgroup = ~species,
            text = ~B) %>%
    add_ribbons(ymin = ~lwr, ymax = ~upr, line = list(width = 0),
                alpha = 0.2, showlegend = FALSE) %>%
    add_lines(y = ~scaled_B)


## Model comparison ------------------------------------------------------------

full <- multispic(inputs, species_cor = "all", temporal_cor = "ar1",
                  log_K_option = par_option(option = "normal_prior",
                                            mean = mean_log_K, sd = sd_log_K),
                  log_B0_option = par_option(option = "normal_prior",
                                             mean = mean_log_B0, sd = sd_log_B0),
                  log_r_option = par_option(option = "normal_prior",
                                            mean = mean_log_r, sd = sd_log_r),
                  log_sd_B_option = par_option(option = "normal_prior",
                                               mean = mean_log_sd, sd = sd_log_sd),
                  log_q_option = par_option(option = "normal_prior",
                                            mean = mean_log_q, sd = sd_log_q),
                  log_sd_I_option = par_option(option = "normal_prior",
                                               mean = mean_log_sd, sd = sd_log_sd),
                  logit_rho_option = par_option(option = "normal_prior",
                                                mean = mean_logit_rho, sd = sd_logit_rho),
                  logit_phi_option = par_option(option = "normal_prior",
                                                mean = mean_logit_phi, sd = sd_logit_phi),
                  n_forecast = 1, K_groups = ~region, pe_covariates = ~winter_nao)

no_nao <- update(full, pe_covariates = NULL)
one_species_cor <- update(no_nao, species_cor = "one")
no_species_cor <- update(one_species_cor, species_cor = "none")
no_temporal_cor <- update(no_species_cor, temporal_cor = "none")

loo_full <- run_loo(full)
loo_no_nao <- run_loo(no_nao)
loo_one_species_cor <- run_loo(one_species_cor)
loo_no_species_cor <- run_loo(no_species_cor)
loo_no_temporal_cor <- run_loo(no_temporal_cor)

full$mAIC
no_nao$mAIC
one_species_cor$mAIC
no_species_cor$mAIC
no_temporal_cor$mAIC


## make a tidy_par function and name par using factor levels
## consider imposing a mean change in the collapse era (current fit is using K to cause the collapse)
## consider using full time series of catch to inform lower bound for K
## check factor levels and make sure indexing isn't messed up for K when K_groups is used




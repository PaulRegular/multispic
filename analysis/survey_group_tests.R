
source("analysis/001_NL_case_study_helpers.R")


list2env(nl_inputs_and_priors(region = "3LNO", species = NULL), envir = globalenv())

## Dev notes:
## - Forcing the RW structure results in unusual process errors for some species
## - During testing the fixed, random, and uniform_prior options rarely converged
int_eff <- multispic(inputs, species_cor = "all", temporal_cor = "ar1",
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
                  n_forecast = 1, K_groups = ~1, survey_groups = ~species * gear * season,
                  pe_covariates = ~0)


mix_eff <- update(int_eff, survey_groups = ~gear + species * season)

main_eff <- update(int_eff, survey_groups = ~gear + species + season)


int_eff$mAIC
mix_eff$mAIC
main_eff$mAIC

int_eff$loo <- run_loo(int_eff)
mix_eff$loo <- run_loo(mix_eff)
main_eff$loo <- run_loo(main_eff)

int_eff$loo$mse
mix_eff$loo$mse
main_eff$loo$mse

int_eff$retro <- run_retro(int_eff, folds = 10)
mix_eff$retro <- run_retro(int_eff, folds = 10)
main_eff$retro <- run_retro(main_eff, folds = 10)

int_eff$retro$mse
mix_eff$retro$mse
main_eff$retro$mse

fit <- main_eff

if (fit$call$K_groups == ~1) {
    K_groups <- NULL
    if (length(unique(fit$landings$species)) == 1) {
        K_label <- "K"
    } else {
        K_label <- "All species"
    }
} else {
    K_groups <- as.formula(fit$call$K_groups)
    K_label <- unique(fit$landings[, all.vars(K_groups)])
}

p <- fit$index %>%
    plot_ly() %>%
    add_trace(x = ~year, y = ~index, color = ~survey,
              colors = viridis::viridis(100), mode = "markers+lines",
              type = "scatter")
p
p %>% layout(yaxis = list(type = "log"))


p <- fit$landings %>%
    plot_ly() %>%
    add_trace(x = ~year, y = ~landings, color = ~species,
              colors = viridis::viridis(100), mode = "markers+lines",
              type = "scatter")
p
p %>% layout(yaxis = list(type = "log"))

fit$tot_pop %>%
    plot_ly(frame = K_groups) %>%
    add_lines(x = ~year, y = ~landings) %>%
    layout(title = "Total landings")

fit$landings %>%
    plot_ly(x = ~year) %>%
    add_lines(y = ~winter_nao, name = "Winter NAO") %>%
    add_lines(y = ~spring_nao, name = "Spring NAO") %>%
    add_lines(y = ~summer_nao, name = "Summer NAO") %>%
    add_lines(y = ~fall_nao, name = "Fall NAO")



## Raw par
par <- as.list(fit$sd_rep, "Est")
hist(unlist(par), breaks = 30)

## Prior and posterior
plot_prior_post(prior_mean = fit$par$mean_log_K,
                prior_sd = exp(fit$par$log_sd_log_K),
                post_mean = fit$par$log_K,
                post_sd = fit$se$log_K,
                post_names = K_label,
                xlab = "log(K)")
plot_prior_post(prior_mean = fit$par$mean_log_r,
                prior_sd = exp(fit$par$log_sd_log_r),
                post_mean = fit$par$log_r,
                post_sd = fit$se$log_r,
                post_names = levels(fit$landings$species),
                xlab = "log(r)")
plot_prior_post(prior_mean = fit$par$mean_log_B0,
                prior_sd = exp(fit$par$log_sd_log_B0),
                post_mean = fit$par$log_B0,
                post_sd = fit$se$log_B0,
                post_names = levels(fit$landings$species),
                xlab = "log(B0)")
plot_prior_post(prior_mean = fit$par$mean_log_sd_B,
                prior_sd = exp(fit$par$log_sd_log_sd_B),
                post_mean = fit$par$log_sd_B,
                post_sd = fit$se$log_sd_B,
                post_names = levels(fit$landings$species),
                xlab = "log(SD<sub>B</sub>)")
plot_prior_post(prior_mean = fit$par$mean_log_q,
                prior_sd = exp(fit$par$log_sd_log_q),
                post_mean = fit$par$log_q,
                post_sd = fit$se$log_q,
                post_names = levels(fit$index$survey),
                xlab = "log(q)")
plot_prior_post(prior_mean = fit$par$mean_log_sd_I,
                prior_sd = exp(fit$par$log_sd_log_sd_I),
                post_mean = fit$par$log_sd_I,
                post_sd = fit$se$log_sd_I,
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

plot_prior_post(prior_mean = fit$par$mean_logit_rho,
                prior_sd = exp(fit$par$log_sd_logit_rho),
                post_mean = fit$par$logit_rho,
                post_sd = fit$se$logit_rho,
                post_names = sp_nm_mat[lower.tri(sp_nm_mat)],
                xlab = "logit(rho)")# , trans_fun = function(x) inv_logit(x, shift = TRUE))

plot_prior_post(prior_mean = fit$obj$env$data$mean_logit_phi,
                prior_sd = fit$obj$env$data$sd_logit_phi,
                post_mean = fit$par$logit_phi,
                post_sd = fit$se$logit_phi,
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
p %>% add_markers(x = ~survey_id, y = ~std_res)
p %>% add_markers(x = ~species, y = ~std_res)
p %>% add_markers(x = ~season, y = ~std_res)
p %>% add_markers(x = ~gear, y = ~std_res)
p %>% add_markers(x = ~paste(gear, season), y = ~std_res)
p %>% add_markers(x = ~paste(species, season), y = ~std_res)
p %>% add_markers(x = ~paste(species, gear), y = ~std_res)


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
    add_lines(y = ~landings, name = "L", color = I("red")) %>%
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






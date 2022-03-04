
library(plotly)
library(dplyr)
library(ggplot2)

normalize <- function(x) {(x - min(x)) / (max(x) - min(x)) }

## Leave-one-out results ---------------------------------------------------------------------------

loo_dat <- lapply(c("2J3K", "3LNO", "3Ps"), function(r) {

    spp_fits <- readRDS(paste0("analysis/exports/spp_fits_", r, ".rds"))
    models <- names(spp_fits)
    names(models) <- c("Full", "Just covariates", "Just shift", "Just climate",
                       "Just correlation", "Shared correlation", "Just species correlation",
                       "Just temporal correlation", "No correlation")

    spp_loo_dat <- lapply(models, function (m) {
        spp_fit <- spp_fits[[m]]
        data.frame(model = names(models)[models == m],
                   n_index = nrow(spp_fit$index),
                   n_random = length(spp_fit$sd_rep$par.random),
                   n_fixed = length(spp_fit$sd_rep$par.fixed),
                   spp_fit$loo$preds)
    })
    spp_loo_dat <- do.call(rbind, spp_loo_dat)

    sp_fit <- readRDS(paste0("analysis/exports/sp_fits_", r, ".rds"))
    sp_loo_dat <- data.frame(model = "Single-species",
                             n_index = nrow(sp_fit$index),
                             n_random = length(sp_fit$sd_rep$par.random),
                             n_fixed = length(sp_fit$sd_rep$par.fixed),
                             sp_fit$loo$preds)

    r_loo_dat <- rbind(spp_loo_dat, sp_loo_dat)
    r_loo_dat$model <- factor(r_loo_dat$model, levels = c(names(models), "Single-species"))
    r_loo_dat

})
loo_dat <- do.call(rbind, loo_dat)
rownames(loo_dat) <- NULL

sr <- do.call(rbind, strsplit(as.character(loo_dat$species), "-"))
loo_dat$species <- sr[, 1]
loo_dat$region <- sr[, 2]

loo_scores <- loo_dat %>%
    group_by(model, species, region) %>%
    summarise(n = n(), n_species = length(unique(species)),
              rmse = sqrt(mean((log_index - log_pred_index) ^ 2, na.rm = TRUE))) %>%
    group_by(species, region) %>%
    mutate(scaled_rmse = normalize(rmse),
           delta_rmse = rmse - min(rmse)) %>%
    ungroup()

overall_loo_scores <- loo_dat %>%
    group_by(model, region) %>%
    summarise(n = n(), n_species = length(unique(species)),
              n_index = unique(n_index), n_random = unique(n_random),
              n_fixed = unique(n_fixed),
              rmse = sqrt(mean((log_index - log_pred_index) ^ 2, na.rm = TRUE))) %>%
    group_by(region) %>%
    mutate(scaled_rmse = normalize(rmse),
           delta_rmse = rmse - min(rmse)) %>%
    ungroup()


loo_scores %>%
    filter(region == "2J3K") %>%
    mutate(model = factor(model)) %>%
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) %>%
    add_lines()

loo_scores %>%
    filter(region == "3LNO") %>%
    mutate(model = factor(model)) %>%
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) %>%
    add_lines()

loo_scores %>%
    filter(region == "3Ps") %>%
    mutate(model = factor(model)) %>%
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) %>%
    add_lines()

overall_loo_scores %>%
    mutate(model = factor(model)) %>%
    plot_ly(x = ~model, y = ~rmse, color = ~region,
            colors = viridis::viridis(100)) %>%
    add_lines() %>% add_text(text = ~n_fixed)


## Hindcast results ---------------------------------------------------------------------------

hind_dat <- lapply(c("2J3K", "3LNO", "3Ps"), function(r) {

    spp_fits <- readRDS(paste0("analysis/exports/spp_fits_", r, ".rds"))
    models <- names(spp_fits)
    names(models) <- c("Full", "Just covariates", "Just shift", "Just climate",
                       "Just correlation", "Shared correlation", "Just species correlation",
                       "Just temporal correlation", "No correlation")

    spp_hind_dat <- lapply(models, function (m) {
        data.frame(model = names(models)[models == m], spp_fits[[m]]$retro$hindcasts)
    })
    spp_hind_dat <- do.call(rbind, spp_hind_dat)

    sp_fit <- readRDS(paste0("analysis/exports/sp_fits_", r, ".rds"))
    sp_hind_dat <- data.frame(model = "Single-species", sp_fit$retro$hindcasts)

    r_hind_dat <- rbind(spp_hind_dat, sp_hind_dat)
    r_hind_dat$model <- factor(r_hind_dat$model, levels = c(names(models), "Single-species"))
    r_hind_dat

})
hind_dat <- do.call(rbind, hind_dat)
rownames(hind_dat) <- NULL

sr <- do.call(rbind, strsplit(as.character(hind_dat$species), "-"))
hind_dat$species <- sr[, 1]
hind_dat$region <- sr[, 2]

hind_scores <- hind_dat %>%
    group_by(model, species, region) %>%
    summarise(n = n(), n_species = length(unique(species)),
              rmse = sqrt(mean((log_index - log_pred_index) ^ 2, na.rm = TRUE))) %>%
    group_by(species, region) %>%
    mutate(scaled_rmse = normalize(rmse),
           delta_rmse = rmse - min(rmse)) %>%
    ungroup()

overall_hind_scores <- hind_dat %>%
    group_by(model, region) %>%
    summarise(rmse = sqrt(mean((log_index - log_pred_index) ^ 2, na.rm = TRUE))) %>%
    group_by(region) %>%
    mutate(scaled_rmse = normalize(rmse),
           delta_rmse = rmse - min(rmse)) %>%
    ungroup()


hind_scores %>%
    filter(region == "2J3K") %>%
    mutate(model = factor(model)) %>%
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) %>%
    add_lines()

hind_scores %>%
    filter(region == "3LNO") %>%
    mutate(model = factor(model)) %>%
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) %>%
    add_lines()

hind_scores %>%
    filter(region == "3Ps") %>%
    mutate(model = factor(model)) %>%
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) %>%
    add_lines()

overall_hind_scores %>%
    plot_ly(x = ~model, y = ~rmse, color = ~region,
            colors = viridis::viridis(100)) %>%
    add_lines()



## Combined plot -----------------------------------------------------------------------------------


loo_p <- overall_loo_scores %>%
    mutate(model = factor(model)) %>%
    plot_ly(x = ~model, y = ~rmse, color = ~region,
            colors = viridis::viridis(100),
            legendgroup = ~region) %>%
    add_lines() %>%
    layout(yaxis = list(title = "LOO-CV Score"))

hind_p <- overall_hind_scores %>%
    plot_ly(x = ~model, y = ~rmse, color = ~region,
            colors = viridis::viridis(100),
            legendgroup = ~region, showlegend = FALSE) %>%
    add_lines() %>%
    layout(yaxis = list(title = "Hind-CV Score"))

subplot(loo_p, hind_p, nrows = 2, shareX = TRUE, titleY = TRUE, titleX = FALSE)



scores <- merge(overall_loo_scores, overall_hind_scores, by = c("model", "region"),
      suffixes = c("_loo", "_hind"))


scores %>%
  plot_ly(y = ~model, frame = ~region) %>%
  add_markers(x = ~rmse_loo, name = "LOO") %>%
  add_markers(x = ~rmse_hind, name = "hindcast")

scores <- scores %>%
    group_by(model, region) %>%
    mutate(mean_score = mean(c(rmse_loo, rmse_hind))) %>%
    group_by(model) %>%
    mutate(overall_mean_score = mean(mean_score)) %>%
    arrange(-overall_mean_score) %>%
    ungroup() %>%
    mutate(model = factor(model, levels = unique(model))) %>%
    as.data.frame()


scores %>%
  filter(region == "3Ps") %>%
  plot_ly(y = ~model, text = ~round(rmse_loo, 1),
          size = ~n_fixed, sizes = c(50, 500)) %>%
  add_paths(x = ~rmse_loo, name = "LOO-CV") %>%
  add_paths(x = ~rmse_hind, name = "Hind-CV") %>%
  # add_markers(x = ~rmse_hind, name = "Hind-CV") %>%
  # add_trace(type = "scatter", mode = "markers+text",
  #           textfont = list(size = 10, color = "white"),
  #           showlegend = FALSE) %>%
  layout(yaxis = list(title = "",
                      categoryorder = "array",
                      categoryarray = levels(scores$model)))


## Figure idea: http://www.dataplusscience.com/images/PewAlt2.PNG
## Think better / worse predictive performance | forecasting skill over single-species approach

scores <- merge(overall_loo_scores, overall_hind_scores, by = c("model", "region"),
                suffixes = c("_loo", "_hind"))

## Scores relative to single-species model
rel_scores <- scores %>%
  group_by(region) %>%
  mutate(ss_rmse_loo = rmse_loo[model == "Single-species"],
         rel_rmse_loo = ss_rmse_loo - rmse_loo,
         ss_rmse_hind = rmse_hind[model == "Single-species"],
         rel_rmse_hind = ss_rmse_hind - rmse_hind) %>%
  filter(model != "Single-species") %>%
  arrange(region, n_fixed) %>%
  ungroup() %>%
  mutate(model = factor(model, levels = unique(model))) %>%
  as.data.frame()

rel_scores %>%
  plot_ly(y = ~model, x = ~rel_rmse_hind, color = ~region) %>%
  add_markers(size = I(100))



## Export select dashboards ------------------------------------------------------------------------

# spp_fits <- readRDS("analysis/exports/spp_fits_2J3K.rds")
# fit <- spp_fits$one_species_cor
# vis_multispic(fit, output_file = "analysis/exports/dashboards/spp_fit_2J3KL_one_species_cor.html")
#
# spp_fits <- readRDS("analysis/exports/spp_fits_3LNO.rds")
# fit <- spp_fits$no_nao
# vis_multispic(fit, output_file = "analysis/exports/dashboards/spp_fit_3LNO_no_nao.html")
# fit <- spp_fits$one_species_cor
# vis_multispic(fit, output_file = "analysis/exports/dashboards/spp_fit_3LNO_one_species_cor.html")
#
# spp_fits <- readRDS("analysis/exports/spp_fits_3Ps.rds")
# fit <- spp_fits$full
# vis_multispic(fit, output_file = "analysis/exports/dashboards/spp_fit_3Ps_full.html")
# fit <- spp_fits$full
# vis_multispic(fit, output_file = "analysis/exports/dashboards/spp_fit_3Ps_just_nao.html")
# fit <- spp_fits$no_nao
# vis_multispic(fit, output_file = "analysis/exports/dashboards/spp_fit_3Ps_no_nao.html")



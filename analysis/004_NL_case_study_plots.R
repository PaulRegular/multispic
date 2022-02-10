
library(plotly)
library(dplyr)
library(ggplot2)

normalize <- function(x) {(x - min(x)) / (max(x) - min(x)) }


## Export select dashboards ------------------------------------------------------------------------

spp_fits <- readRDS("analysis/exports/spp_fits_2J3K.rds")
fit <- spp_fits$one_species_cor
vis_multispic(fit, output_file = "analysis/exports/dashboards/spp_fit_2J3KL_one_species_cor.html")


## Leave-one-out results ---------------------------------------------------------------------------

loo_dat <- lapply(c("2J3K", "3LNO", "3Ps"), function(r) {

    spp_fits <- readRDS(paste0("analysis/exports/spp_fits_", r, ".rds"))
    sp_fits <- readRDS(paste0("analysis/exports/sp_fits_", r, ".rds"))
    models <- names(spp_fits)
    names(models) <- c("Full", "No NAO", "Just NAO", "Shared correlation",
                       "Just temporal correlation", "No correlation")

    spp_loo_dat <- lapply(models, function (m) {
        spp_fit <- spp_fits[[m]]
        data.frame(model = names(models)[models == m],
                   year = spp_fit$index$year, species = spp_fit$index$species,
                   obs = spp_fit$loo$obs, pred = spp_fit$loo$pred)
    })
    spp_loo_dat <- do.call(rbind, spp_loo_dat)

    species <- names(sp_fits)
    sp_loo_dat <- lapply(species, function (sp) {
        sp_fit <- sp_fits[[sp]]
        if (!(length(sp_fit) == 1 && sp_fit == "Did not converge")) {
            if (sp_fit$sd_rep$pdHess) {
                data.frame(model = "Single-species",
                           year = sp_fit$index$year, species = sp_fit$index$species,
                           obs = sp_fit$loo$obs, pred = sp_fit$loo$pred)
            }
        }
    })
    sp_loo_dat <- do.call(rbind, sp_loo_dat)

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
              rmse = sqrt(mean((obs - pred) ^ 2, na.rm = TRUE))) %>%
    group_by(species, region) %>%
    mutate(scaled_rmse = normalize(rmse),
           delta_rmse = rmse - min(rmse)) %>%
    ungroup()

sp_loo_groups <- loo_dat %>%
    filter(model == "Single-species") %>%
    mutate(ysr = paste0(year, "-", species, "-", region))
sp_loo_groups <- unique(sp_loo_groups$ysr)

overall_loo_scores <- loo_dat %>%
    mutate(ysr = paste0(year, "-", species, "-", region)) %>%
    filter(ysr %in% sp_loo_groups) %>%  # limit to species that converged in the single-species analysis
    group_by(model, region) %>%
    summarise(n = n(), n_species = length(unique(species)),
              rmse = sqrt(mean((obs - pred) ^ 2, na.rm = TRUE))) %>%
    group_by(region) %>%
    mutate(scaled_rmse = normalize(rmse),
           delta_rmse = rmse - min(rmse)) %>%
    ungroup()

loo_scores %>%
    filter(region == "2J3K") %>%
    filter(!model %in% c("No NAO", "Just NAO")) %>%
    mutate(model = factor(model)) %>%
    plot_ly(x = ~model, y = ~delta_rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) %>%
    add_lines()

loo_scores %>%
    filter(region == "3LNO") %>%
    filter(!model %in% c("No NAO", "Just NAO")) %>%
    mutate(model = factor(model)) %>%
    plot_ly(x = ~model, y = ~delta_rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) %>%
    add_lines()

loo_scores %>%
    filter(region == "3Ps") %>%
    filter(!model %in% c("No NAO", "Just NAO")) %>%
    mutate(model = factor(model)) %>%
    plot_ly(x = ~model, y = ~delta_rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) %>%
    add_lines()


overall_loo_scores %>%
    filter(!model %in% c("No NAO", "Just NAO")) %>%
    mutate(model = factor(model)) %>%
    plot_ly(x = ~model, y = ~delta_rmse, color = ~region,
            colors = viridis::viridis(100)) %>%
    add_lines()


## Hindcast results ---------------------------------------------------------------------------

hind_dat <- lapply(c("2J3K", "3LNO", "3Ps"), function(r) {

    spp_fits <- readRDS(paste0("analysis/exports/spp_fits_", r, ".rds"))
    sp_fits <- readRDS(paste0("analysis/exports/sp_fits_", r, ".rds"))
    models <- names(spp_fits)
    names(models) <- c("Full", "No NAO", "Just NAO", "Shared correlation",
                       "Just temporal correlation", "No correlation")

    spp_hind_dat <- lapply(models, function (m) {
        data.frame(model = names(models)[models == m], spp_fits[[m]]$retro$hindcasts)
    })
    spp_hind_dat <- do.call(rbind, spp_hind_dat)

    species <- names(sp_fits)
    sp_hind_dat <- lapply(species, function (sp) {
        sp_fit <- sp_fits[[sp]]
        if (!(length(sp_fit) == 1 && sp_fit == "Did not converge")) {
            data.frame(model = "Single-species", sp_fits[[sp]]$retro$hindcasts)
        }
    })
    sp_hind_dat <- do.call(rbind, sp_hind_dat)

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

sp_hind_groups <- hind_dat %>%
    filter(model == "Single-species") %>%
    mutate(ysr = paste0(year, "-", species, "-", region))
sp_hind_groups <- unique(sp_hind_groups$ysr)

overall_hind_scores <- hind_dat %>%
    mutate(ysr = paste0(year, "-", species, "-", region)) %>%
    filter(ysr %in% sp_hind_groups) %>%  # limit to species that converged in the single-species analysis
    group_by(model, region) %>%
    summarise(n = n(), n_species = length(unique(species)),
              rmse = sqrt(mean((log_index - log_pred_index) ^ 2, na.rm = TRUE))) %>%
    group_by(region) %>%
    mutate(scaled_rmse = normalize(rmse),
           delta_rmse = rmse - min(rmse)) %>%
    ungroup()

## May not be plotting correctly given different species in different regions
hind_scores %>%
    plot_ly(x = ~model, y = ~delta_rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) %>%
    add_lines()

overall_hind_scores %>%
    plot_ly(x = ~model, y = ~delta_rmse, color = ~region,
            colors = viridis::viridis(100)) %>%
    add_lines()


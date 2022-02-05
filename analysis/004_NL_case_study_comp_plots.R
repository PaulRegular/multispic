
library(plotly)
library(dplyr)
library(ggplot2)

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
            data.frame(model = "Single-species",
                       year = sp_fit$index$year, species = sp_fit$index$species,
                       obs = sp_fit$loo$obs, pred = sp_fit$loo$pred)
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
    mutate(scaled_rmse = scale(rmse)[, 1]) %>%
    ungroup()

overall_loo_scores <- loo_dat %>%
    group_by(model, region) %>%
    summarise(n = n(), n_species = length(unique(species)),
              rmse = sqrt(mean((obs - pred) ^ 2, na.rm = TRUE))) %>%
    group_by(region) %>%
    mutate(scaled_rmse = scale(rmse)[, 1]) %>%
    ungroup()

## Not a perfect comparison because some species are missing
overall_loo_scores %>%
    filter(region == "3Ps")

spp_names <- sort(unique(loo_dat$species))
loo_scores$species <- factor(loo_scores$species, levels = c("Overall", spp_names))

loo_scores %>%
    plot_ly(x = ~model, y = ~scaled_rmse, color = ~species, frame = ~region) %>%
    add_lines()


overall_loo_scores %>%
    plot_ly(x = ~model, y = ~scaled_rmse, color = ~region) %>%
    add_lines()



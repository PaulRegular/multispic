
library(dplyr)

normalize <- function(x) {(x - min(x)) / (max(x) - min(x)) }

## Leave-one-out results ---------------------------------------------------------------------------

loo_dat <- lapply(c("2J3K", "3LNO", "3Ps"), function(r) {

    fits <- readRDS(paste0("analysis/NL_case_study/exports/fits_", r, ".rds"))
    models <- names(fits)
    names(models) <- c("Full", "Just shift", "Just correlation", "Shared correlation",
                       "Just species correlation", "Just temporal correlation", "No correlation",
                       "Single-species")

    loo_dat <- lapply(models, function (m) {
        fit <- fits[[m]]
        data.frame(model = names(models)[models == m],
                   n_index = nrow(fit$index),
                   n_random = length(fit$sd_rep$par.random),
                   n_fixed = length(fit$sd_rep$par.fixed),
                   rec_id = seq(nrow(fit$loo$preds)),
                   fit$loo$preds)
    })
    loo_dat <- do.call(rbind, loo_dat)

    loo_dat$model <- factor(loo_dat$model, levels = names(models))
    loo_dat

})
loo_dat <- do.call(rbind, loo_dat)
rownames(loo_dat) <- NULL

sr <- do.call(rbind, strsplit(as.character(loo_dat$species), "-"))
loo_dat$species <- sr[, 1]
loo_dat$region <- sr[, 2]

## Drop cases with convergence issues to ensure metrics are comparable across all cases
ind <- is.na(loo_dat$log_pred_index)
sum(ind)
drop_recs <- unique(loo_dat[ind, c("region", "rec_id")])
ind <- paste0(loo_dat$region, "-", loo_dat$rec_id) %in%
    paste0(drop_recs$region, "-", drop_recs$rec_id)
loo_dat <- loo_dat[!ind, ]

loo_scores <- loo_dat |>
    group_by(model, species, region) |>
    summarise(n = n(), n_species = length(unique(species)),
              rmse = sqrt(mean((log_index - log_pred_index) ^ 2, na.rm = TRUE))) |>
    group_by(species, region) |>
    mutate(scaled_rmse = normalize(rmse),
           delta_rmse = rmse - min(rmse),
           ranked_rmse = rank(rmse)) |>
    ungroup()

overall_loo_scores <- loo_dat |>
    group_by(model, region) |>
    summarise(n = n(), n_species = length(unique(species)),
              n_index = unique(n_index), n_random = unique(n_random),
              n_fixed = unique(n_fixed),
              rmse = sqrt(mean((log_index - log_pred_index) ^ 2, na.rm = TRUE))) |>
    group_by(region) |>
    mutate(scaled_rmse = normalize(rmse),
           delta_rmse = rmse - min(rmse),
           ranked_rmse = rank(rmse)) |>
    ungroup()


## Hindcast results ---------------------------------------------------------------------------

drop_years <- vector(mode = "character")
hind_dat <- lapply(c("2J3K", "3LNO", "3Ps"), function(r) {

    fits <- readRDS(paste0("analysis/NL_case_study/exports/fits_", r, ".rds"))
    models <- names(fits)
    names(models) <- c("Full", "Just shift", "Just correlation", "Shared correlation",
                       "Just species correlation", "Just temporal correlation", "No correlation",
                       "Single-species")

    hind_dat <- lapply(models, function (m) {
        missing_years <<- is.na(fits[[m]]$retro$retro_fits) |> which() |> names()
        if (length(missing_years) > 0) missing_years <<- paste0(r, "-", missing_years)
        drop_years <<- c(drop_years, missing_years)
        # cat(r, "|", m, "| missing years:", missing_years, "\n\n")
        data.frame(model = names(models)[models == m], fits[[m]]$retro$hindcasts)
    })
    hind_dat <- do.call(rbind, hind_dat)

    hind_dat$model <- factor(hind_dat$model, levels = names(models))
    hind_dat

})
hind_dat <- do.call(rbind, hind_dat)
rownames(hind_dat) <- NULL

sr <- do.call(rbind, strsplit(as.character(hind_dat$species), "-"))
hind_dat$species <- sr[, 1]
hind_dat$region <- sr[, 2]

## Drop cases with convergence issues to ensure metrics are comparable across all cases
drop_years <- unique(drop_years)
drop_years
ind <- paste0(hind_dat$region, "-", hind_dat$retro_year) %in% drop_years
hind_dat <- hind_dat[!ind, ]

hind_scores <- hind_dat |>
    group_by(model, species, region) |>
    summarise(n = n(), n_species = length(unique(species)),
              rmse = sqrt(mean((log_index - log_pred_index) ^ 2, na.rm = TRUE))) |>
    group_by(species, region) |>
    mutate(scaled_rmse = normalize(rmse),
           delta_rmse = rmse - min(rmse),
           ranked_rmse = rank(rmse)) |>
    ungroup()

overall_hind_scores <- hind_dat |>
    group_by(model, region) |>
    summarise(rmse = sqrt(mean((log_index - log_pred_index) ^ 2, na.rm = TRUE))) |>
    group_by(region) |>
    mutate(scaled_rmse = normalize(rmse),
           delta_rmse = rmse - min(rmse),
           ranked_rmse = rank(rmse)) |>
    ungroup()



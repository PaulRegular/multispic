
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
                   rec_id = seq(nrow(spp_fit$loo$preds)),
                   spp_fit$loo$preds)
    })
    spp_loo_dat <- do.call(rbind, spp_loo_dat)

    sp_fit <- readRDS(paste0("analysis/exports/sp_fits_", r, ".rds"))
    sp_loo_dat <- data.frame(model = "Single-species",
                             n_index = nrow(sp_fit$index),
                             n_random = length(sp_fit$sd_rep$par.random),
                             n_fixed = length(sp_fit$sd_rep$par.fixed),
                             rec_id = seq(nrow(sp_fit$loo$preds)),
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


loo_scores |>
    filter(region == "2J3K") |>
    mutate(model = factor(model)) |>
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) |>
    add_lines()

loo_scores |>
    filter(region == "3LNO") |>
    mutate(model = factor(model)) |>
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) |>
    add_lines()

loo_scores |>
    filter(region == "3Ps") |>
    mutate(model = factor(model)) |>
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) |>
    add_lines()

overall_loo_scores |>
    mutate(model = factor(model)) |>
    plot_ly(x = ~model, y = ~rmse, color = ~region,
            colors = viridis::viridis(100)) |>
    add_lines() |> add_text(text = ~n_fixed)


## Hindcast results ---------------------------------------------------------------------------

drop_years <- vector(mode = "character")
hind_dat <- lapply(c("2J3K", "3LNO", "3Ps"), function(r) {

    spp_fits <- readRDS(paste0("analysis/exports/spp_fits_", r, ".rds"))
    models <- names(spp_fits)
    names(models) <- c("Full", "Just covariates", "Just shift", "Just climate",
                       "Just correlation", "Shared correlation", "Just species correlation",
                       "Just temporal correlation", "No correlation")

    spp_hind_dat <- lapply(models, function (m) {
        missing_years <<- is.na(spp_fits[[m]]$retro$retro_fits) |> which() |> names()
        if (length(missing_years) > 0) missing_years <<- paste0(r, "-", missing_years)
        drop_years <<- c(drop_years, missing_years)
        # cat(r, "|", m, "| missing years:", missing_years, "\n\n")
        data.frame(model = names(models)[models == m], spp_fits[[m]]$retro$hindcasts)
    })
    spp_hind_dat <- do.call(rbind, spp_hind_dat)

    sp_fit <- readRDS(paste0("analysis/exports/sp_fits_", r, ".rds"))
    sp_hind_dat <- data.frame(model = "Single-species", sp_fit$retro$hindcasts)

    missing_years <<- is.na(sp_fit$retro$retro_fits) |> which() |> names()
    if (length(missing_years) > 0) missing_years <<- paste0(r, "-", missing_years)
    drop_years <<- c(drop_years, missing_years)

    r_hind_dat <- rbind(spp_hind_dat, sp_hind_dat)
    r_hind_dat$model <- factor(r_hind_dat$model, levels = c(names(models), "Single-species"))
    r_hind_dat

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


hind_scores |>
    filter(region == "2J3K") |>
    mutate(model = factor(model)) |>
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) |>
    add_lines()

hind_scores |>
    filter(region == "3LNO") |>
    mutate(model = factor(model)) |>
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) |>
    add_lines()

hind_scores |>
    filter(region == "3Ps") |>
    mutate(model = factor(model)) |>
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) |>
    add_lines()

overall_hind_scores |>
    plot_ly(x = ~model, y = ~rmse, color = ~region,
            colors = viridis::viridis(100)) |>
    add_lines()



## Combined plot -----------------------------------------------------------------------------------


scores <- merge(overall_loo_scores, overall_hind_scores,
                by = c("model", "region"),
                suffixes = c("_loo", "_hind"))

a <- scores |>
    arrange(model) |>
    plot_ly(y = ~model, color = ~region, legendgroup = ~region,
            colors = viridis::viridis(3)) |>
    add_paths(x = ~rmse_loo) |>
    add_trace(x = ~rmse_loo, text = ~ranked_rmse_loo,
              type = "scatter", mode = "lines+markers+text",
              marker = list(color = "white", size = 20,
                            line = list(width = 2)),
              textfont = list(size = 10),
              showlegend = FALSE) |>
    layout(yaxis = list(title = "", autorange = "reversed"),
           xaxis = list(title = "LOO-CV Score", side ="top"))


b <- scores |>
    arrange(model) |>
    plot_ly(y = ~model, color = ~region, legendgroup = ~region,
            colors = viridis::viridis(3), showlegend = FALSE) |>
    add_paths(x = ~rmse_hind) |>
    add_trace(x = ~rmse_hind, text = ~ranked_rmse_hind,
              type = "scatter", mode = "lines+markers+text",
              marker = list(color = "white", size = 20,
                            line = list(width = 2)),
              textfont = list(size = 10)) |>
    layout(yaxis = list(title = "", autorange = "reversed"),
           xaxis = list(title = "Hindcast-CV Score", side ="top"))

subplot(a, b, nrows = 1, shareY = TRUE, titleX = TRUE)



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



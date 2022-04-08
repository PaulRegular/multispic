
library(multispic)
library(plotly)

save_html <- function(p, file = file, selfcontained = TRUE, ...) {
    ext <- tools::file_ext(file)
    tmp_file <- paste0("tmp.", ext)
    htmlwidgets::saveWidget(partial_bundle(p), file = tmp_file,
                            selfcontained = selfcontained, ...)
    invisible(file.rename(tmp_file, file))
}

## Plot trends estimated by the just correlation model

fits_2J3K <- readRDS("analysis/NL_case_study/exports/spp_fits_2J3K.rds")
fit_2J3K <- fits_2J3K$just_cor

fits_3LNO <- readRDS("analysis/NL_case_study/exports/spp_fits_3LNO.rds")
fit_3LNO <- fits_3LNO$just_cor

fits_3Ps <- readRDS("analysis/NL_case_study/exports/spp_fits_3Ps.rds")
fit_3Ps <- fits_3Ps$just_cor

spp_2J3K <- gsub("-2J3K", "", levels(fit_2J3K$landings$species))
spp_3LNO <- gsub("-3LNO", "", levels(fit_3LNO$landings$species))
spp_3Ps <- gsub("-3Ps", "", levels(fit_3Ps$landings$species))

all_spp <- c(spp_2J3K, spp_3LNO, spp_3Ps) |>
    unique() |>
    sort()
all_spp

## Rough phylogenetic order + limit species to improve visual contrast

all_spp <- c("Redfish spp.", "Wolffish spp.",
             "Yellowtail Flounder", "Witch Flounder", "American Plaice", "Greenland Halibut",
             "White Hake", "Haddock", "Atlantic Cod",
             "Skate spp.")
core_spp <- c("Redfish spp.",
             "Yellowtail Flounder", "American Plaice", "Greenland Halibut",
             "Atlantic Cod",
             "Skate spp.")
spp_cols <- viridis::viridis(length(core_spp))
# see::material_colors()
# spp_cols <- see::material_colors() |> head(length(core_spp))
# spp_cols <- see::material_colors()[c("red", "pink",
#                                      "yellow", "brown", "deep purple", "green",
#                                      "purple", "teal", "blue",
#                                      "amber")]
names(spp_cols) <- core_spp
spp_cols

## Identify species included across all regions
spp_tab <- table(c(spp_2J3K, spp_3LNO, spp_3Ps))[core_spp]
common_spp <- names(which(spp_tab == 3))


## Function for extracting species correlations
get_spp_rho <- function(fit, tag) {
    n_spp_cor <- length(fit$par$logit_rho)

    sp_rho <- sp_nm_mat <- matrix(NA, nrow = nlevels(fit$pop$species), ncol = nlevels(fit$pop$species))
    rownames(sp_rho) <- colnames(sp_rho) <- gsub(tag, "", levels(fit$pop$species))
    sp_rho[lower.tri(sp_rho)] <- inv_logit(fit$par$logit_rho, shift = TRUE)
    sp_rho <- t(sp_rho)
    sp_rho[lower.tri(sp_rho)] <- inv_logit(fit$par$logit_rho, shift = TRUE)
    diag(sp_rho) <- 1
    round(sp_rho, 2)
    for (i in seq(nrow(sp_rho))) {
        for (j in seq(ncol(sp_rho))) {
            sp_nm_mat[i, j] <- paste(levels(fit$pop$species)[i], "-", levels(fit$pop$species)[j])
        }
    }
    ind <- match(all_spp, rownames(sp_rho)) |> na.omit() # consistent order
    sp_rho[ind, ind]
}

fit_2J3K$spp_rho <- get_spp_rho(fit_2J3K, "-2J3K")
fit_3LNO$spp_rho <- get_spp_rho(fit_3LNO, "-3LNO")
fit_3Ps$spp_rho <- get_spp_rho(fit_3Ps, "-3Ps")

## Function for dropping region tag from species name
simplify_spp_name <- function(fit, tag) {
    for (df in c("pop", "index")) {
        fit[[df]]$species <- gsub(tag, "", fit[[df]]$species)
        fit[[df]]$species <- factor(fit[[df]]$species, levels = core_spp)
    }
    fit
}

fit_2J3K <- simplify_spp_name(fit_2J3K, "-2J3K")
fit_3LNO <- simplify_spp_name(fit_3LNO, "-3LNO")
fit_3Ps <- simplify_spp_name(fit_3Ps, "-3Ps")

## create dummy data for legend
dummy_data <- expand.grid(year = min(fit_2J3K$tot_pop$year), x = 0,
                          species = factor(core_spp, levels = core_spp))






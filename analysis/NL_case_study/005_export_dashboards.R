
## Export select dashboards ------------------------------------------------------------------------

## Single species
for (r in c("2J3K", "3LNO", "3Ps")) {
    sp_fits <- readRDS(paste0("analysis/NL_case_study/exports/sp_fits_", r, ".rds"))
    vis_multispic(sp_fits,
                  output_file = paste0("analysis/NL_case_study/exports/dashboards/sp_fits_", r, ".html"))
}

## Just correlation model as it is generally ranked high
for (r in c("2J3K", "3LNO", "3Ps")) {
    spp_fits <- readRDS(paste0("analysis/NL_case_study/exports/spp_fits_", r, ".rds"))
    spp_fits_just_cor <- spp_fits$just_cor
    vis_multispic(spp_fits_just_cor,
                  output_file = paste0("analysis/NL_case_study/exports/dashboards/spp_fits_just_cor_", r, ".html"))
}



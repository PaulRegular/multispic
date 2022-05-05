
## Export select dashboards ------------------------------------------------------------------------

## Just correlation model as it is generally ranked high
for (r in c("2J3K", "3LNO", "3Ps")) {
    fits <- readRDS(paste0("analysis/NL_case_study/exports/fits_", r, ".rds"))
    vis_multispic(fits$null,
                  output_file = paste0("analysis/NL_case_study/exports/dashboards/null_", r, ".html"))
    vis_multispic(fits$just_cor,
                  output_file = paste0("analysis/NL_case_study/exports/dashboards/just_cor_", r, ".html"))
}



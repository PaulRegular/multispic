
## Export select dashboards ------------------------------------------------------------------------

region_name <- c("2J3K" = "Northeast NL Shelf",
                 "3LNO" = "Grand Bank",
                 "3Ps" = "Southern NL")

## Just correlation model as it is generally ranked high
for (r in c("2J3K", "3LNO", "3Ps")) {

    fits <- readRDS(paste0("analysis/NL_case_study/exports/fits_", r, ".rds"))

    vis_multispic(fits$null,
                  output_file = paste0("analysis/NL_case_study/exports/dashboards/Single-species - ", region_name[r], ".html"))
    vis_multispic(fits$just_cor,
                  output_file = paste0("analysis/NL_case_study/exports/dashboards/Just correlation - ", region_name[r], ".html"))

}

zip("analysis/paper/supplement/model_dashboards.zip",
    list.files("analysis/NL_case_study/exports/dashboards/", full.names = TRUE),
    extras = "-j")

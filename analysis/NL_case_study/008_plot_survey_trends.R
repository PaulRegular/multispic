
library(plotly)

source("analysis/NL_case_study/006_compile_trends.R")

region_trends <- function(fit, title, show_spp = TRUE) {

    spp_plots <- vector("list", length(common_spp))
    names(spp_plots) <- common_spp

    for (s in common_spp) {

        sub_index <- fit$index |>
            filter(species == s)
        sub_pop <- fit$pop |>
            filter(species == s)
        sub_col <- spp_cols[s]

        spp_plots[[s]] <- plot_ly(x = ~year) |>
            add_ribbons(data = sub_pop, ymin = ~B_lwr, ymax = ~B_upr,
                        line = list(width = 0), color = I(sub_col),
                        alpha = 0.2, showlegend = FALSE, legendgroup = s) %>%
            add_lines(data = sub_pop, y = ~B,
                      color = I(sub_col), name = s,
                      line = list(width = 2), legendgroup = s,
                      alpha = 0.5, showlegend = FALSE) |>
            add_lines(data = sub_index, y = ~pred, split = ~survey,
                      color = I(sub_col), name = s, #
                      line = list(width = 1),
                      showlegend = FALSE, legendgroup = s) |>
            add_markers(data = sub_index, y = ~index,
                        color = I(sub_col), name = s,
                        showlegend = FALSE, legendgroup = s,
                        marker = list(size = 2)) |>
            add_annotations(x = 0.5, y = 1, xref = "paper", yref = "paper", text = title,
                            font = list(size = 16), showarrow = FALSE,
                            xanchor = "center",  yanchor = "bottom",
                            visible = s == "Redfish spp.") |>
            layout(xaxis = list(title = "Year"),
                   yaxis = list(title = "Biomass (kt)"))


    }

    subplot(spp_plots, nrows = 4, shareX = TRUE)

}

p_2J3K <- region_trends(fit_2J3K, "Northeast NL Shelf")
p_3LNO <- region_trends(fit_3LNO, "Grand Bank")
p_3Ps <- region_trends(fit_3Ps, "Southern NL")

## Vertical annotations produced using the function did not align properly; hence manual approach below
p <- subplot(p_2J3K, p_3LNO, p_3Ps) |>
    add_annotations(x = rep(1, 4), y = c(0.96, 0.74, 0.31, 0.06), text = common_spp,
                    xref = "paper", yref = "paper",
                    font = list(size = 16), showarrow = FALSE, textangle = 90,
                    xanchor = "left", yanchor = "center") |>
    add_annotations(x = -0.04, y = 0.5, text = "Biomass (kt)",
                    xref = "paper", yref = "paper",
                    font = list(size = 16), showarrow = FALSE, textangle = 270,
                    xanchor = "right", yanchor = "center") |>
    layout(margin = list(r = 20))
p

reticulate::py_run_string("import sys")
save_image(p, file = "analysis/NL_case_study/exports/plots/survey_trends.svg",
           width = 1000, height = 700)
file.copy("analysis/NL_case_study/exports/plots/survey_trends.svg", "analysis/paper/figures/survey_trends.svg",
          overwrite = TRUE)
save_html(p, file = "analysis/NL_case_study/exports/plots/survey_trends.html")
saveRDS(p, file = "analysis/NL_case_study/exports/plots/survey_trends.rds")


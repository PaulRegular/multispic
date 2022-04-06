
library(plotly)

source("analysis/NL_case_study/006_compile_trends.R")


fit <- fit_2J3K

spp_plots <- vector("list", length(all_spp))
names(spp_plots) <- all_spp

for (s in all_spp) {

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
                  line = list(width = 1.5), legendgroup = s) |>
        add_lines(data = sub_index, y = ~pred,
                  color = I(sub_col), name = s,
                  line = list(width = 1.5, dash = "dot"),
                  showlegend = FALSE, legendgroup = s) |>
        add_markers(data = sub_index, y = ~index,
                    color = I(sub_col), name = s,
                    showlegend = FALSE, legendgroup = s,
                    marker = list(size = 3, color = "white", line = list(width = 1.2))) |>
        layout(xaxis = list(title = "Year"),
               yaxis = list(title = "Biomass (kt)"))


}

subplot(spp_plots[c(1, 3, 4)], nrows = 3, shareX = TRUE)


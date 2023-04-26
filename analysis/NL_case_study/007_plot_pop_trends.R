
library(plotly)

source("analysis/NL_case_study/006_compile_trends.R")

region_trends <- function(fit, title, show_legend = TRUE, show_axis_titles = TRUE) {

    pop <- fit$pop
    pop$B_pe <- pop$B - (pop$B / pop$pe)
    pop$B_net <- pop$B_growth - pop$landings

    extra_tots <- pop |>
        group_by(year) |>
        summarise(total_growth = sum(B_growth),
                  total_pe = sum(B_pe),
                  total_net = sum(B_net))

    tot_pop <- fit$tot_pop |>
        merge(extra_tots, by = "year")


    a <- plot_ly(colors = spp_cols) |>
        add_trace(data = tot_pop, x = ~year, y = ~B, color = I("lightgrey"),
                  type = 'scatter', mode = 'lines', fill = 'tozeroy',
                  name = "Total", legendgroup = "Total", line = list(width = 0),
                  showlegend = show_legend) |>
        add_lines(data = tot_pop, x = ~year, y = ~K, color = I("black"),
                  name = "K", legendgroup = "K",
                  showlegend = show_legend) |>
        add_lines(data = pop, x = ~year, y = ~B,
                  color = ~species, legendgroup = ~species, showlegend = FALSE,
                  line = list(width = 2)) |>
        add_lines(data = dummy_data, x = ~year, y = ~x, color = ~species,
                  legendgroup = ~species, line = list(width = 2),
                  showlegend = show_legend) |>
        add_annotations(x = 0.5, y = 1, xref = "paper", yref = "paper", text = title,
                        xanchor = "center",  yanchor = "bottom",
                        font = list(size = 16), showarrow = FALSE) |>
        layout(xaxis = list(title = ifelse(show_axis_titles, "Year", "")),
               yaxis = list(title = ifelse(show_axis_titles, "Biomass (kt)", "")))

    b <- plot_ly(colors = spp_cols) |>
        add_trace(data = tot_pop, x = ~year, y = ~total_growth, color = I("lightgrey"),
                  type = 'scatter', mode = 'lines', fill = 'tozeroy',
                  name = "Total", legendgroup = "Total", line = list(width = 0),
                  showlegend = FALSE) |>
        add_lines(data = pop, x = ~year, y = ~B_growth,
                  color = ~species, legendgroup = ~species,
                  showlegend = FALSE, line = list(width = 2)) |>
        layout(xaxis = list(title = ifelse(show_axis_titles, "Year", "")),
               yaxis = list(title = ifelse(show_axis_titles, "Expected production", "")))

    c <- plot_ly(colors = spp_cols) |>
        add_trace(data = tot_pop, x = ~year, y = ~-landings, color = I("lightgrey"),
                  type = 'scatter', mode = 'lines', fill = 'tozeroy',
                  name = "Total", legendgroup = "Total", line = list(width = 0),
                  showlegend = FALSE) |>
        add_lines(data = pop, x = ~year, y = ~-landings,
                  color = ~species, legendgroup = ~species,
                  showlegend = FALSE, line = list(width = 2)) |>
        layout(xaxis = list(title = ifelse(show_axis_titles, "Year", "")),
               yaxis = list(title = ifelse(show_axis_titles, "Landings", "")))

    d <- plot_ly(colors = spp_cols) |>
        add_trace(data = tot_pop, x = ~year, y = ~total_pe, color = I("lightgrey"),
                  type = 'scatter', mode = 'lines', fill = 'tozeroy',
                  name = "Total", legendgroup = "Total", line = list(width = 0),
                  showlegend = FALSE) |>
        add_lines(data = pop, x = ~year, y = ~B_pe,
                  color = ~species, legendgroup = ~species,
                  showlegend = FALSE, line = list(width = 2)) |>
        layout(xaxis = list(title = ifelse(show_axis_titles, "Year", "")),
               yaxis = list(title = ifelse(show_axis_titles, "Residual change", "")))

    subplot(a, b, c, d, nrows = 4, shareX = TRUE, titleY = TRUE,
            heights = c(0.4, 0.2, 0.2, 0.2))


}


p_2J3K <- region_trends(fit_2J3K, "Northeast NL Shelf", TRUE, TRUE)
p_3LNO <- region_trends(fit_3LNO, "Grand Bank", FALSE, FALSE)
p_3Ps <- region_trends(fit_3Ps, "Southern NL", FALSE, FALSE)

p <- subplot(p_2J3K, p_3LNO, p_3Ps, titleY = TRUE)
p

reticulate::use_miniconda('r-reticulate')
reticulate::py_run_string("import sys")
save_image(p, file = "analysis/NL_case_study/exports/plots/pop_trends.svg",
           width = 1000, height = 700)
file.copy("analysis/NL_case_study/exports/plots/pop_trends.svg", "analysis/paper/figures/pop_trends.svg",
          overwrite = TRUE)
save_html(p, file = "analysis/NL_case_study/exports/plots/pop_trends.html")
saveRDS(p, file = "analysis/NL_case_study/exports/plots/pop_trends.rds")





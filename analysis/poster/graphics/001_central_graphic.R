
library(plotly)
library(bezier)

source("analysis/NL_case_study/006_compile_trends.R")

bezier_smooth <- function(x, y, nout = 1000) {
    t <- seq(0, 1, length = nout)
    p <- cbind(x, y)
    s <- bezier::bezier(t = t, p = p)
}

fit <- fit_2J3K

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

smooth_tot_pop <- bezier_smooth(tot_pop$year, tot_pop$B) |>
    as.data.frame()
names(smooth_tot_pop) <- c("year", "B")

split_pop <- split(pop, pop$species)
for (nm in names(split_pop)) {
    if (nrow(split_pop[[nm]]) > 0) {
        s <- bezier_smooth(split_pop[[nm]]$year, split_pop[[nm]]$B)
        split_pop[[nm]] <- data.frame(species = nm, year = s[, 1], B = s[, 2])
    } else {
        split_pop[[nm]] <- NULL
    }
}
smooth_pop <- do.call(rbind, split_pop)


p <- plot_ly(colors = viridis::viridis(100)) |>
    add_trace(data = smooth_tot_pop, x = ~year, y = ~B, color = I("lightgrey"),
              type = 'scatter', mode = 'lines', fill = 'tozeroy',
              name = "Total", legendgroup = "Total",
              line = list(width = 0),
              showlegend = FALSE) |>
    add_lines(data = tot_pop, x = ~year, y = ~K, color = I("darkgrey"),
              name = "Carrying capacity", legendgroup = "Carrying capacity",
              showlegend = FALSE, line = list(dash = "dot", width = 5)) |>
    add_lines(data = smooth_pop, x = ~year, y = ~B,
              color = ~species, legendgroup = ~species, showlegend = FALSE,
              line = list(width = 5)) |>
    layout(
        xaxis = list(
            title = "",
            range = c(1996, 2020),
            linewidth = 5,
            showgrid = FALSE,
            zeroline = TRUE,
            showline = TRUE,
            showticklabels = FALSE
        ),
        yaxis = list(
            title = "",
            range = c(0, 1500),
            linewidth = 5,
            showgrid = FALSE,
            zeroline = TRUE,
            showline = TRUE,
            showticklabels = FALSE
        )
    )

save_image(p, file = "analysis/poster/graphics/pop_trends.svg",
           width = 700, height = 700)


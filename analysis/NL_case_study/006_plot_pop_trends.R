
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

all_spp <- c(spp_2J3K, spp_3LNO, spp_3Ps)|>
    unique() |>
    sort()
all_spp

## Rough phylogenetic order + limit species to improve visual contrast

# all_spp <- c("Redfish spp.", "Wolffish spp.",
#              "Yellowtail Flounder", "Witch Flounder", "American Plaice", "Greenland Halibut",
#              "White Hake", "Haddock", "Atlantic Cod",
#              "Skate spp.")
all_spp <- c("Redfish spp.",
             "Yellowtail Flounder", "American Plaice", "Greenland Halibut",
             "Atlantic Cod",
             "Skate spp.")
spp_cols <- viridis::viridis(length(all_spp))
# see::material_colors()
# spp_cols <- see::material_colors() |> head(length(all_spp))
# spp_cols <- see::material_colors()[c("red", "pink",
#                                      "yellow", "brown", "deep purple", "green",
#                                      "purple", "teal", "blue",
#                                      "amber")]
names(spp_cols) <- all_spp
spp_cols

fit_2J3K$pop$species <- gsub("-2J3K", "", fit_2J3K$pop$species)
fit_2J3K$pop$species <- factor(fit_2J3K$pop$species, levels = all_spp)

fit_3LNO$pop$species <- gsub("-3LNO", "", fit_3LNO$pop$species)
fit_3LNO$pop$species <- factor(fit_3LNO$pop$species, levels = all_spp)

fit_3Ps$pop$species <- gsub("-3Ps", "", fit_3Ps$pop$species)
fit_3Ps$pop$species <- factor(fit_3Ps$pop$species, levels = all_spp)

## create dummy data for legend
dummy_data <- expand.grid(year = min(fit_2J3K$tot_pop$year), x = 0,
                          species = factor(all_spp, levels = all_spp))


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
                  line = list(width = 1.2)) |>
        add_lines(data = dummy_data, x = ~year, y = ~x, color = ~species,
                  legendgroup = ~species, line = list(width = 1.2),
                  showlegend = show_legend) |>
        add_annotations(x = 0.5, y = 1.07, xref = "paper", yref = "paper", text = title,
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
                  showlegend = FALSE, line = list(width = 1.2)) |>
        layout(xaxis = list(title = ifelse(show_axis_titles, "Year", "")),
               yaxis = list(title = ifelse(show_axis_titles, "Expected production", "")))

    c <- plot_ly(colors = spp_cols) |>
        add_trace(data = tot_pop, x = ~year, y = ~-landings, color = I("lightgrey"),
                  type = 'scatter', mode = 'lines', fill = 'tozeroy',
                  name = "Total", legendgroup = "Total", line = list(width = 0),
                  showlegend = FALSE) |>
        add_lines(data = pop, x = ~year, y = ~-landings,
                  color = ~species, legendgroup = ~species,
                  showlegend = FALSE, line = list(width = 1.2)) |>
        layout(xaxis = list(title = ifelse(show_axis_titles, "Year", "")),
               yaxis = list(title = ifelse(show_axis_titles, "Landings", "")))

    d <- plot_ly(colors = spp_cols) |>
        add_trace(data = tot_pop, x = ~year, y = ~total_pe, color = I("lightgrey"),
                  type = 'scatter', mode = 'lines', fill = 'tozeroy',
                  name = "Total", legendgroup = "Total", line = list(width = 0),
                  showlegend = FALSE) |>
        add_lines(data = pop, x = ~year, y = ~B_pe,
                  color = ~species, legendgroup = ~species,
                  showlegend = FALSE, line = list(width = 1.2)) |>
        layout(xaxis = list(title = ifelse(show_axis_titles, "Year", "")),
               yaxis = list(title = ifelse(show_axis_titles, "Residual change", "")))

    subplot(a, b, c, d, nrows = 4, shareX = TRUE, titleY = TRUE,
            heights = c(0.4, 0.2, 0.2, 0.2))


}


p_2J3K <- region_trends(fit_2J3K, "2J3K", TRUE, TRUE)
p_3LNO <- region_trends(fit_3LNO, "3LNO", FALSE, FALSE)
p_3Ps <- region_trends(fit_3Ps, "3Ps", FALSE, FALSE)

p <- subplot(p_2J3K, p_3LNO, p_3Ps, titleY = TRUE)
p

save_image(p, file = "analysis/NL_case_study/exports/plots/pop_trends.svg",
           width = 1000, height = 700)
save_html(p, file = "analysis/NL_case_study/exports/plots/pop_trends.html")
saveRDS(p, file = "analysis/NL_case_study/exports/plots/pop_trends.rds")


## Consider exporting and displaying density dependent subtractions
## along with additions and subtractions from process error
## And subtractions from landings

## Figure out how to display uncertainty without making the plot too busy






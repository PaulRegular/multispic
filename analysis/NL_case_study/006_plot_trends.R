
library(plotly)

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

p_2J3K <- plot_ly(colors = spp_cols) |>
    add_trace(data = fit_2J3K$tot_pop, x = ~year, y = ~B, color = I("grey"),
              type = 'scatter', mode = 'lines', fill = 'tozeroy',
              name = "Total", legendgroup = "Total", line = list(width = 0)) |>
    add_lines(data = fit_2J3K$tot_pop, x = ~year, y = ~K, color = I("black"),
              name = "K", legendgroup = "K") |>
    add_lines(data = fit_2J3K$pop, x = ~year, y = ~B,
              color = ~species, legendgroup = ~species, showlegend = FALSE) |>
    add_lines(data = dummy_data, x = ~year, y = ~x, color = ~species,
              legendgroup = ~species) |>
    add_annotations(x = 0.5, y = 1.07, xref = "paper", yref = "paper", text = "2J3K",
                    font = list(size = 16), showarrow = FALSE) |>
    layout(xaxis = list(title = "Year"),
           yaxis = list(title = "Biomass (kt)"))

p_3LNO <- plot_ly(colors = spp_cols) |>
    add_trace(data = fit_3LNO$tot_pop, x = ~year, y = ~B, color = I("grey"),
              type = 'scatter', mode = 'lines', fill = 'tozeroy',
              name = "Total", legendgroup = "Total", line = list(width = 0),
              showlegend = FALSE) |>
    add_lines(data = fit_3LNO$tot_pop, x = ~year, y = ~K, color = I("black"),
              name = "K", legendgroup = "K", showlegend = FALSE) |>
    add_lines(data = fit_3LNO$pop, x = ~year, y = ~B,
              color = ~species, legendgroup = ~species, showlegend = FALSE) |>
    add_annotations(x = 0.5, y = 1.07, xref = "paper", yref = "paper", text = "3LNO",
                    font = list(size = 16), showarrow = FALSE) |>
    layout(xaxis = list(title = "Year"),
           yaxis = list(title = "Biomass (kt)"))

p_3Ps <- plot_ly(colors = spp_cols) |>
    add_trace(data = fit_3Ps$tot_pop, x = ~year, y = ~B, color = I("grey"),
              type = 'scatter', mode = 'lines', fill = 'tozeroy',
              name = "Total", legendgroup = "Total", line = list(width = 0),
              showlegend = FALSE) |>
    add_lines(data = fit_3Ps$tot_pop, x = ~year, y = ~K, color = I("black"),
              name = "K", legendgroup = "K", showlegend = FALSE) |>
    add_lines(data = fit_3Ps$pop, x = ~year, y = ~B,
              color = ~species, legendgroup = ~species, showlegend = FALSE) |>
    add_annotations(x = 0.5, y = 1.07, xref = "paper", yref = "paper", text = "3Ps",
                    font = list(size = 16), showarrow = FALSE) |>
    layout(xaxis = list(title = "Year"),
           yaxis = list(title = "Biomass (kt)"))


subplot(p_2J3K, p_3LNO, p_3Ps, nrows = 1, shareY = FALSE)



## Consider exporting and displaying density dependent subtractions
## along with additions and subtractions from process error
## And subtractions from landings

## Figure out how to display uncertainty without making the plot too busy






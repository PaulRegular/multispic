
library(plotly)
library(dplyr)
library(ggplot2)

normalize <- function(x) {(x - min(x)) / (max(x) - min(x)) }

save_html <- function(p, file = file, selfcontained = TRUE, ...) {
    ext <- tools::file_ext(file)
    tmp_file <- paste0("tmp.", ext)
    htmlwidgets::saveWidget(partial_bundle(p), file = tmp_file,
                            selfcontained = selfcontained, ...)
    invisible(file.rename(tmp_file, file))
}

source("analysis/NL_case_study/003_compile_scores.R")

## Visualize some of the detailed results ----------------------------------------------------------

loo_dat |>
    plot_ly(x = ~year, y = ~log_pred_index - log_index,
            color = ~survey, colors = viridis::viridis(100),
            frame = ~model) |>
    add_markers()

loo_dat |>
    filter(region == "3LNO") |>
    plot_ly(x = ~year, y = ~log_pred_index - log_index,
            color = ~factor(survey), colors = viridis::viridis(100),
            frame = ~model) |>
    add_markers()

loo_dat |>
    filter(region == "3LNO") |>
    plot_ly(x = ~year, y = ~log_pred_index - log_index,
            color = ~model, colors = viridis::viridis(100),
            frame = ~factor(survey)) |>
    add_markers()

loo_scores |>
    filter(region == "2J3K") |>
    mutate(model = factor(model)) |>
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) |>
    add_lines()

loo_scores |>
    filter(region == "3LNO") |>
    mutate(model = factor(model)) |>
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) |>
    add_lines()

loo_scores |>
    filter(region == "3Ps") |>
    mutate(model = factor(model)) |>
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) |>
    add_lines()

overall_loo_scores |>
    mutate(model = factor(model)) |>
    plot_ly(x = ~model, y = ~rmse, color = ~region,
            colors = viridis::viridis(100)) |>
    add_lines() |> add_text(text = ~n_fixed)

hind_scores |>
    filter(region == "2J3K") |>
    mutate(model = factor(model)) |>
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) |>
    add_lines()

hind_scores |>
    filter(region == "3LNO") |>
    mutate(model = factor(model)) |>
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) |>
    add_lines()

hind_scores |>
    filter(region == "3Ps") |>
    mutate(model = factor(model)) |>
    plot_ly(x = ~model, y = ~rmse, color = ~species, frame = ~region,
            colors = viridis::viridis(100)) |>
    add_lines()

overall_hind_scores |>
    plot_ly(x = ~model, y = ~rmse, color = ~region,
            colors = viridis::viridis(100)) |>
    add_lines()


## Combined plot -----------------------------------------------------------------------------------

scores <- merge(overall_loo_scores, overall_hind_scores,
                by = c("model", "region"),
                suffixes = c("_loo", "_hind"))

mean_score <- bind_rows(overall_loo_scores, overall_hind_scores) |>
    group_by(model) |>
    summarise(mean_score = mean(rmse)) |>
    mutate(ranked_score = rank(mean_score)) |>
    ungroup()

a <- scores |>
    arrange(model) |>
    plot_ly(y = ~model, color = ~region, legendgroup = ~region,
            colors = viridis::viridis(3)) |>
    add_paths(x = ~rmse_loo) |>
    add_trace(x = ~rmse_loo, text = ~ranked_rmse_loo,
              type = "scatter", mode = "lines+markers+text",
              line = list(width = 2),
              marker = list(color = "white", size = 20,
                            line = list(width = 2)),
              textfont = list(size = 10),
              showlegend = FALSE) |>
    layout(yaxis = list(title = "", autorange = "reversed"),
           xaxis = list(title = "LOO-CV score", side ="top"))


b <- scores |>
    arrange(model) |>
    plot_ly(y = ~model, color = ~region, legendgroup = ~region,
            colors = viridis::viridis(3), showlegend = FALSE) |>
    add_paths(x = ~rmse_hind) |>
    add_trace(x = ~rmse_hind, text = ~ranked_rmse_hind,
              type = "scatter", mode = "lines+markers+text",
              line = list(width = 2),
              marker = list(color = "white", size = 20,
                            line = list(width = 2)),
              textfont = list(size = 10)) |>
    layout(yaxis = list(title = "", autorange = "reversed"),
           xaxis = list(title = "Hindcast-CV score", side ="top"))

c <- mean_score |>
    arrange(model) |>
    plot_ly(y = ~model, color = I("grey30"), showlegend = FALSE) |>
    add_trace(x = ~mean_score, text = ~ranked_score,
              type = "scatter", mode = "lines+markers+text",
              line = list(width = 2),
              marker = list(color = "white", size = 20,
                            line = list(width = 2)),
              textfont = list(size = 10)) |>
    layout(yaxis = list(title = "", autorange = "reversed"),
           xaxis = list(title = "Mean score", side ="top"))

p <- subplot(a, b, c, nrows = 1, shareY = TRUE, titleX = TRUE,
             widths = c(0.4, 0.4, 0.2)) |>
    layout(legend = list(orientation = "h", x = 0.3, y = 0.01))
p

save_image(p, file = "analysis/NL_case_study/exports/plots/scores.svg", width = 750, height = 350)
save_html(p, file = "analysis/NL_case_study/exports/plots/scores.html")
saveRDS(p, file = "analysis/NL_case_study/exports/plots/scores.rds")


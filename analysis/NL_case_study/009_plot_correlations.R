
library(plotly)

source("analysis/NL_case_study/006_compile_trends.R")

cor_plot <- function(fit, title, show_scale = TRUE) {

    spp_rho <- fit$spp_rho

    spp_rho_df <- as.data.frame.table(spp_rho)
    names(spp_rho_df) <- c("spp1", "spp2", "rho")
    spp_rho_df$size <- abs(spp_rho_df$rho)
    spp_rho_df$lower_tri <- c(lower.tri(spp_rho, diag = TRUE))
    spp_rho_df$upper_tri <- c(upper.tri(spp_rho, diag = TRUE))
    spp_rho_df$bubble_rho <- spp_rho_df$text_rho <- spp_rho_df$rho
    spp_rho_df$bubble_rho[spp_rho_df$lower_tri] <- NA
    spp_rho_df$text_rho[spp_rho_df$upper_tri] <- NA

    n_spp <- nlevels(spp_rho_df$spp1)

    p <- plot_ly(data = spp_rho_df, x = ~as.numeric(spp1), y = ~as.numeric(spp2)) |>
        add_markers(text = ~round(rho, 2), color = ~bubble_rho,
                    colors = c("#B2182B", "white", "#2166AC"),
                    size = ~abs(bubble_rho), sizes = c(10, 200),
                    marker = list(opacity = 1, sizeref = 0.5, line = list(width = 0)),
                    name = "ρ", showlegend = FALSE) |>
        add_text(text = ~round(text_rho, 2), showlegend = FALSE) |>
        add_segments(x = rep(-0.5, n_spp + 1), xend = rep(7.5, n_spp + 1),
                     y = seq(0.5, n_spp + 1), yend = seq(0.5, n_spp + 1),
                     line = list(width = 1, color = "lightgrey"),
                     showlegend = FALSE, hoverinfo = "none") |>
        add_segments(y = rep(-0.5, n_spp + 1), yend = rep(7.5, n_spp + 1),
                     x = seq(0.5, n_spp + 1), xend = seq(0.5, n_spp + 1),
                     line = list(width = 1, color = "lightgrey"),
                     showlegend = FALSE, hoverinfo = "none") |>
        add_annotations(x = 0.5, y = 1, yshift = 110,
                        xref = "paper", yref = "paper", text = title,
                        font = list(size = 16), showarrow = FALSE,
                        xanchor = "center",  yanchor = "bottom") |>
        colorbar(limits = c(-1.001, 1.001), y = 0.42, title = "ρ") |>
        layout(yaxis = list(title = "", showgrid = FALSE, zerolinecolor = "white",
                            tickvals = seq(nlevels(spp_rho_df$spp1)),
                            ticktext = levels(spp_rho_df$spp1),
                            range = c(n_spp + 0.45, 0.55),
                            fixedrange = TRUE),
               xaxis = list(title = "", side = "top", tickangle = 90,
                            showgrid = FALSE, zerolinecolor = "white",
                            tickvals = seq(nlevels(spp_rho_df$spp1)),
                            ticktext = levels(spp_rho_df$spp1),
                            range = c(0.55, n_spp + 0.45),
                            fixedrange = TRUE))
    if (!show_scale) p <- p |> hide_colorbar()
    p
}


cor_2J3K <- cor_plot(fit_2J3K, "", TRUE)
cor_3LNO <- cor_plot(fit_3LNO, "", FALSE)
cor_3Ps <- cor_plot(fit_3Ps, "", FALSE)

cor_p <- subplot(cor_2J3K, cor_3LNO, cor_3Ps,
                 margin = c(0.12, 0.01, 0.12, 0.01),
                 widths = c(0.27, 0.38, 0.35))


pe_plot <- function(fit, title, show_legend = TRUE) {
    fit$pop |>
        plot_ly(color = ~species, colors = spp_cols,
                legendgroup = ~species) |>
        add_lines(x = ~year, y = ~log_std_res_pe, showlegend = show_legend) |>
        add_annotations(x = 0.5, y = 1,
                        xref = "paper", yref = "paper", text = title,
                        font = list(size = 16), showarrow = FALSE,
                        xanchor = "center",  yanchor = "bottom") |>
        layout(xaxis = list(title = "Year"),
               yaxis = list(title = "Standardized process error"))
}


pe_2J3K <- pe_plot(fit_2J3K, "Northeast NL Shelf", TRUE)
pe_3LNO <- pe_plot(fit_3LNO, "Grand Bank", FALSE)
pe_3Ps <- pe_plot(fit_3Ps, "Southern NL", FALSE)

pe_p <- subplot(pe_2J3K, pe_3LNO, pe_3Ps,
                margin = c(0.12, 0.01, 0.12, 0.01),
                widths = c(0.27, 0.38, 0.35))

p <- subplot(pe_p, cor_p, nrows = 2, heights = c(0.4, 0.6),
             margin = c(0.01, 0.01, 0.22, 0.01)) |>
    add_annotations(x = 0, y = 0.65, xshift = -25,
                    xref = "paper", yref = "paper", text = "Standardized process error",
                    font = list(size = 13), showarrow = FALSE, textangle = 270,
                    xanchor = "right",  yanchor = "bottom") |>
    layout(legend = list(y = 1))
p

save_image(p, file = "analysis/NL_case_study/exports/plots/pe_cor.svg",
           width = 1300, height = 700)
file.copy("analysis/NL_case_study/exports/plots/pe_cor.svg", "analysis/paper/figures/pe_cor.svg",
          overwrite = TRUE)
save_html(p, file = "analysis/NL_case_study/exports/plots/pe_cor.html")
saveRDS(p, file = "analysis/NL_case_study/exports/plots/pe_cor.rds")


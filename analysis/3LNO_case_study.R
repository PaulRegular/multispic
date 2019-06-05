
## TODO:
## - Think more about the cor option; add "couple" and "assume" options to par_option??
## - Calculate one-step ahead residuals

library(units)
library(plotly)
library(TMB)
library(multispic)

unique(multispic::landings$species)

index <- multispic::index
landings <- multispic::landings
covariates <- multispic::covariates


## Setup the data --------------------------------------------------------------

## Subset the data
sub_sp <- unique(multispic::landings$species)
# sub_sp <- sub_sp[sub_sp != "Haddock"]
# sub_sp <- c("Yellowtail", "Witch", "Cod", "Plaice", "Redfish", "Skate")
# sub_sp <- c("Cod", "Plaice", "Yellowtail", "Redfish", "Witch")
# sub_sp <- c("Yellowtail", "Plaice", "Skate", "Cod", "Witch", "Redfish")
# sub_sp <- c("Cod", "Yellowtail", "Plaice")
start_year <- 1977
end_year <- 2018
index <- index[index$year >= start_year & index$year <= end_year &
                   index$species %in% sub_sp, ]
landings <- landings[landings$year >= start_year & landings$year <= end_year &
                         landings$species %in% sub_sp, ]
covariates <- covariates[covariates$year >= start_year & covariates$year <= end_year, ]

## Assume Yankee Q = Engel Q
# index$gear[index$gear == "Yankee"] <- "Engel"


## Set-up indices for TMB
landings$species <- factor(landings$species)
landings$y <- factor(landings$year)
landings$sy <- factor(paste0(landings$species, "-", landings$year))
landings <- landings[order(landings$sy), ]
index$sy <- factor(paste0(index$species, "-", index$year), levels = levels(landings$sy))
index$survey <- factor(paste0(index$species, "-", index$season, "-", index$gear))
index$gear_season <- factor(paste0(index$gear, "-", index$season))
index$gear_species <- factor(paste0(index$gear, "-", index$species))
index$species <- factor(index$species)
index$null <- factor(rep("null", nrow(index)))

p <- index %>%
    group_by(survey) %>%
    plot_ly() %>%
    add_lines(x = ~year, y = ~index, color = ~species,
              colors = viridis::viridis(100))
p
p %>% layout(yaxis = list(type = "log"))

## Exploratory plots
index %>%
    group_by(survey) %>%
    mutate(scaled_index = scale(index)) %>%
    plot_ly() %>%
    add_lines(x = ~year, y = ~scaled_index, color = ~species,
              colors = viridis::viridis(100)) %>%
    layout(yaxis2 = list(overlaying = "y"))

p <- landings %>%
    group_by(stock) %>%
    plot_ly() %>%
    add_lines(x = ~year, y = ~landings, color = ~species,
              colors = viridis::viridis(100))
p
p %>% layout(yaxis = list(type = "log"))

landings %>%
    group_by(year) %>%
    summarise(total_landings = sum(landings)) %>%
    plot_ly() %>%
    add_lines(x = ~year, y = ~total_landings)


## Run model -------------------------------------------------------------------

inputs <- list(landings = landings, index = index, covariates = covariates)
fit <- fit_model(inputs, survey_group = "survey", cor_str = "all",
                 logit_cor_option = par_option(option = "fixed", mean = -1, sd = 1),
                 log_B0_option = par_option(option = "fixed", mean = -1, sd = 1),
                 log_r_option = par_option(option = "prior", mean = -1, sd = 1),
                 log_sd_B_option = par_option(option = "prior", mean = -1, sd = 1),
                 log_q_option = par_option(option = "prior", mean = -1, sd = 1),
                 log_sd_I_option = par_option(option = "prior", mean = -1, sd = 1),
                 formula = NULL)
fit$opt$message
fit$sd_rep
fit$opt$objective
fit$mAIC

## Raw par
par <- as.list(fit$sd_rep, "Est")
hist(unlist(par), breaks = 30)

## Mean and 95% limits of general prior (mean = -1, sd = 1)
mean <- exp(-1)
lwr <- exp(-1 - (qnorm(0.975)))
upr <- exp(-1 + (qnorm(0.975)))
c(lwr, mean, upr)

## Prior and posterior
post_mean <- as.list(fit$sd_rep, "Est")
post_sd <- as.list(fit$sd_rep, "Std. Error")
plot_prior_post(prior_mean = -1, prior_sd = 1,
                post_mean = post_mean$log_r,
                post_sd = post_sd$log_r,
                post_names = levels(landings$species),
                xlab = "log(r)")
plot_prior_post(prior_mean = -1, prior_sd = 1,
                post_mean = post_mean$log_sd_B,
                post_sd = post_sd$log_sd_B,
                post_names = levels(landings$species),
                xlab = "log(SD<sub>B</sub>)")
plot_prior_post(prior_mean = -1, prior_sd = 1,
                post_mean = post_mean$log_q,
                post_sd = post_sd$log_q,
                post_names = levels(index$survey),
                xlab = "log(q)")
plot_prior_post(prior_mean = -1, prior_sd = 1,
                post_mean = post_mean$log_sd_I,
                post_sd = post_sd$log_sd_I,
                post_names = levels(index$survey),
                xlab = "log(SD<sub>I</sub>)")


## Visually assess par
par <- fit$par
q <- exp(par$log_q)
names(q) <- levels(index$survey)
round(q, 2)
sd_I <- exp(par$log_sd_I)
names(sd_I) <- levels(index$survey)
round(sd_I, 2)
K <- exp(par$log_K)
signif(K, 2)
r <- exp(par$log_r)
names(r) <- levels(index$species)
round(r, 2)
sd_B <- exp(par$log_sd_B)
names(sd_B) <- levels(index$species)
round(sd_B, 2)
B0 <- exp(par$log_B0)
round(B0)
cor <- 2.0 / (1.0 + exp(-par$logit_cor)) - 1.0
round(cor, 2)


## Explore parameter correlations
sd_rep <- fit$sd_rep
dsd <- sqrt(diag(sd_rep$cov.fixed))
cor_mat <- diag(1 / dsd) %*% sd_rep$cov.fixed %*% diag(1 / dsd)
rownames(cor_mat) <- paste0(seq(nrow(cor_mat)), ":", names(dsd))
colnames(cor_mat) <- paste0(seq(nrow(cor_mat)), ":", names(dsd))
cor_tab <- as.data.frame.table(cor_mat)
names(cor_tab) <- c("x", "y", "z")
# corrplot::corrplot.mixed(cor_mat, diag = "n", lower = "ellipse", upper = "number")
cor_tab %>%
    plot_ly(x = ~x, y = ~y, z = ~z, text = ~text,
            colors = c("#B2182B", "white", "#2166AC")) %>%
    add_heatmap() %>% colorbar(limits = c(-1, 1)) %>%
    layout(xaxis = list(title = "", showticklabels = TRUE),
           yaxis = list(title = "", showticklabels = TRUE))

## Index residuals
p <- fit$index %>%
    plot_ly(color = ~species, colors = viridis::viridis(100))
p %>% add_markers(x = ~year, y = ~std_res)
p %>% add_markers(x = ~log(pred), y = ~std_res)
p %>% add_markers(x = ~survey, y = ~std_res)

## Process error residuals
p <- fit$pop %>%
    plot_ly(color = ~species, colors = viridis::viridis(100))
p %>% add_lines(x = ~year, y = ~pe)

## pe vs covariates
fit$pop %>%
    dplyr::left_join(inputs$covariates, by = c("year")) %>%
    plot_ly(color = ~species, colors = viridis::viridis(100)) %>%
    add_markers(x = ~core_cil, y = ~pe, text = ~year)
fit$pop %>%
    dplyr::left_join(inputs$covariates, by = c("year")) %>%
    plot_ly(color = ~species, colors = viridis::viridis(100)) %>%
    add_markers(x = ~nao, y = ~pe, text = ~year)


## Correlation in pe
pe_wide <- tidyr::spread(fit$pop[, c("year", "species", "pe")], species, pe)
pe_wide$year <- NULL
plot(pe_wide)
pe_wide <- as.matrix(pe_wide)
pe_wide[is.infinite(pe_wide)] <- NA
cor_mat <- cor(pe_wide, use = "na.or.complete")
plot_ly(x = rownames(cor_mat), y = rownames(cor_mat), z = ~cor_mat,
        colors = c("#B2182B", "white", "#2166AC")) %>%
    add_heatmap() %>% colorbar(limits = c(-1, 1))
corrplot::corrplot.mixed(cor_mat, diag = "n", lower = "ellipse", upper = "number")


## Fits to the index
p <- fit$index %>%
    group_by(survey) %>%
    plot_ly(x = ~year, color = ~species, colors = viridis::viridis(100),
            legendgroup = ~species) %>%
    add_ribbons(ymin = ~pred_lwr, ymax = ~pred_upr, line = list(width = 0),
                alpha = 0.2, showlegend = FALSE) %>%
    add_lines(y = ~pred) %>%
    add_markers(y = ~index, showlegend = FALSE)
p
p %>% layout(yaxis = list(type = "log"))


## Biomass
p <- fit$pop %>%
    plot_ly(x = ~year, color = ~species, colors = viridis::viridis(100),
            legendgroup = ~species) %>%
    add_ribbons(ymin = ~B_lwr, ymax = ~B_upr, line = list(width = 0),
                alpha = 0.2, showlegend = FALSE) %>%
    add_lines(y = ~B)
p
p %>% layout(yaxis = list(type = "log"))


## F
p <- fit$pop %>%
    # filter(is.finite(F)) %>%
    plot_ly(x = ~year, color = ~species, colors = viridis::viridis(100),
            legendgroup = ~species) %>%
    add_ribbons(ymin = ~F_lwr, ymax = ~F_upr, line = list(width = 0),
                alpha = 0.2, showlegend = FALSE) %>%
    add_lines(y = ~F)
p


## Compare to accepted assessment model results --------------------------------


assess <- read.csv("analysis/stock_assessment_estimates.csv")
names(assess) <- c("species_div", "year", "B", "B_type")
x <- data.table::tstrsplit(assess$species, split = " ")
assess$species <- x[[1]]
assess$division <- x[[2]]
assess$source <- "assessment"
assess$B_lwr <- NA
assess$B_upr <- NA

spm <- fit$pop
spm$source <- "multispic"

keep <- c("year", "species", "B", "B_lwr", "B_upr", "source")
comp <- rbind(assess[, keep], spm[, keep])

comp %>%
    group_by(species, source) %>%
    mutate(scaled_B = scale(B),
           center = attr(scale(B), "scaled:center"),
           scale = attr(scale(B), "scaled:scale")) %>%
    mutate(lwr = (B_lwr - center) / scale,
           upr = (B_upr - center) / scale) %>%
    plot_ly(x = ~year, color = ~species, colors = viridis::viridis(100),
            linetype = ~source, legendgroup = ~species) %>%
    add_ribbons(ymin = ~lwr, ymax = ~upr, line = list(width = 0),
                alpha = 0.2, showlegend = FALSE) %>%
    add_lines(y = ~scaled_B)



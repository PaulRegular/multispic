
## TODO:
## - Revert to one K?
## - Use priors from Miller and Myers?
## - Consider estimating one sd for the process error
## - Think more about B0 now that it works back to the 60s
## - Add a rho parameter of sorts to the K for the relative contribution
##   of each species to total K
## - Return to converted data for species you can
## - Restrict range of q and r in nlminb?
## - Fix q of Campelen-Fall to 1? For cod...maybe?
## - Start comparing parameter estimates to those in the assessments
## - Get latest catch numbers for 3O redfish and 3LNO hake (2016 and 2017 are guesses)
## - Calculate one-step ahead residuals
## - Use covariates to estimate time varrying K or r?
## - Get hake landings for 2017

library(units)
library(plotly)
library(TMB)
library(multispic)

unique(multispic::landings$species)

index <- multispic::index
landings <- multispic::landings


## Setup the data --------------------------------------------------------------

## Subset the data
sub_sp <- unique(multispic::landings$species)
# sub_sp <- c("Yellowtail", "Witch", "Cod", "Plaice", "Redfish", "Skate")
# sub_sp <- c("Cod", "Plaice", "Yellowtail", "Redfish", "Witch")
# sub_sp <- c("Yellowtail", "Plaice", "Skate", "Cod", "Witch", "Redfish")
# sub_sp <- c("Cod", "Yellowtail", "Plaice")
start_year <- 1984
end_year <- 2017
index <- index[index$year >= start_year & index$year <= end_year &
                   index$species %in% sub_sp, ]
landings <- landings[landings$year >= start_year & landings$year <= end_year &
                         landings$species %in% sub_sp, ]

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
              colors = viridis::viridis(100))

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

## Prior visuals
curve(dlnorm(x, meanlog = 5, sdlog = 0.5), 0, 1000)
curve(dlnorm(x, meanlog = -2, sdlog = 0.2), 0, 1.5)
curve(dlnorm(x, meanlog = 0, sdlog = 1), 0, 5)

inputs <- list(landings = landings, index = index)
fit <- fit_model(inputs, survey_group = "survey", cor_str = "all",
                 log_K_option = par_option(option = "prior", mean = 5, sd = 0.5),
                 log_r_option = par_option(option = "prior", mean = -1, sd = 0.5),
                 log_sd_B_option = par_option(option = "prior", mean = -1, sd = 0.5),
                 log_q_option = par_option(option = "prior", mean = -1, sd = 0.5),
                 log_sd_I_option = par_option(option = "prior", mean = -1, sd = 0.5))
fit$opt$message
fit$sd_rep

## Visually assess par
par <- fit$par
hist(unlist(par), breaks = 30)
q <- exp(par$log_q)
names(q) <- levels(index$survey)
round(q, 2)
sd_I <- exp(par$log_sd_I)
names(sd_I) <- levels(index$survey)
round(sd_I, 2)
K <- exp(par$log_K)
names(K) <- levels(index$species)
signif(K, 2)
r <- exp(par$log_r)
names(r) <- levels(index$species)
round(r, 2)
sd_B <- exp(par$log_sd_B)
names(sd_B) <- levels(index$species)
round(sd_B, 2)


## Explore parameter correlations
sd_rep <- fit$sd_rep
dsd <- sqrt(diag(sd_rep$cov.fixed))
cor_mat <- diag(1 / dsd) %*% sd_rep$cov.fixed %*% diag(1 / dsd)
rownames(cor_mat) <- paste0(seq(nrow(cor_mat)), ":", names(dsd))
colnames(cor_mat) <- paste0(seq(nrow(cor_mat)), ":", names(dsd))
cor_tab <- as.data.frame.table(cor_mat)
names(cor_tab) <- c("x", "y", "z")
corrplot::corrplot.mixed(cor_mat, diag = "n", lower = "ellipse", upper = "number")
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
p <- fit$pe %>%
    plot_ly(color = ~species, colors = viridis::viridis(100))
p %>% add_lines(x = ~year, y = ~pe)


## Correlation in pe
pe_wide <- tidyr::spread(fit$pe[, c("year", "species", "pe")], species, pe)
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
p <- fit$biomass %>%
    plot_ly(x = ~year, color = ~species, colors = viridis::viridis(100),
            legendgroup = ~species) %>%
    add_ribbons(ymin = ~B_lwr, ymax = ~B_upr, line = list(width = 0),
                alpha = 0.2, showlegend = FALSE) %>%
    add_lines(y = ~B)
p
p %>% layout(yaxis = list(type = "log"))



## TODO:
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

## Scale index and landings
index$index <- index$index / 100
landings$landings <- landings$landings / 100

## Subset the data
sub_sp <- unique(multispic::landings$species)
# sub_sp <- c("Yellowtail", "Witch", "Cod", "Plaice", "Redfish", "Skate")
# sub_sp <- c("Cod", "Plaice", "Yellowtail", "Redfish", "Witch")
sub_sp <- c("Yellowtail", "Plaice", "Skate", "Cod", "Witch", "Redfish")
start_year <- 1975
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
index$species <- factor(index$species)

p <- index %>%
    group_by(survey) %>%
    plot_ly() %>%
    add_lines(x = ~year, y = ~index, color = ~species,
              colors = viridis::viridis(100))
p
p %>% layout(yaxis = list(type = "log"))


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


n_cor <- sum(lower.tri(matrix(NA, nrow = nlevels(landings$species),
                              ncol = nlevels(landings$species))))
dat <- list(L = as.numeric(landings$landings),
            L_species = as.numeric(landings$species) - 1,
            L_year = as.numeric(landings$y) - 1,
            I = as.numeric(index$index),
            I_species = as.numeric(index$species) - 1,
            I_survey = as.numeric(index$gear) - 1,
            I_sy = as.numeric(index$sy) - 1,
            min_B = 0.001,
            nY = max(as.numeric(landings$y)),
            nS = max(as.numeric(landings$species)))
par <- list(log_B = matrix(0, nrow = dat$nY, ncol = dat$nS),
            log_sd_B = rep(-1, nlevels(landings$species)),
            logit_cor = rep(0, n_cor),
            log_K = 0,
            log_r = rep(-1, nlevels(landings$species)),
            log_m = rep(log(2), nlevels(landings$species)),
            log_q = rep(-1, nlevels(index$gear)),
            log_sd_I = rep(-1, nlevels(index$gear)))
map <- list(log_m = factor(rep(NA, nlevels(landings$species))),
            logit_cor = factor(rep(1, length(par$logit_cor))))

obj <- MakeADFun(dat, par, map = map, random = "log_B", DLL = "multispic")
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(eval.max = 1000, iter.max = 1000))
sd_rep <- sdreport(obj)
obj$report()
exp(opt$par)
sd_rep

## Parameter estimates
lapply(as.list(sd_rep, "Est"), exp)

## Extract some estimates
est <- split(unname(sd_rep$value), names(sd_rep$value))
sd <- split(sd_rep$sd, names(sd_rep$value))
lwr <- split(unname(sd_rep$value) - 1.96 * sd_rep$sd, names(sd_rep$val))
upr <- split(unname(sd_rep$value) + 1.96 * sd_rep$sd, names(sd_rep$val))
index$pred <- exp(est$log_pred_I)
index$pred_lwr <- exp(lwr$log_pred_I)
index$pred_upr <- exp(upr$log_pred_I)
index$std_res <- est$std_res_I

## Explore parameter correlations
dsd <- sqrt(diag(sd_rep$cov.fixed))
cor_mat <- diag(1 / dsd) %*% sd_rep$cov.fixed %*% diag(1 / dsd)
rownames(cor_mat) <- paste0(seq(nrow(cor_mat)), ":", names(dsd))
colnames(cor_mat) <- paste0(seq(nrow(cor_mat)), ":", names(dsd))
cor_tab <- as.data.frame.table(cor_mat)
names(cor_tab) <- c("x", "y", "z")
cor_tab %>%
    plot_ly(x = ~x, y = ~y, z = ~z, text = ~text,
            colors = c("#B2182B", "white", "#2166AC")) %>%
    add_heatmap() %>%
    layout(xaxis = list(title = "", showticklabels = FALSE),
           yaxis = list(title = "", showticklabels = FALSE))

## Index residuals
p <- index %>%
    plot_ly(color = ~species, colors = viridis::viridis(100))
p %>% add_markers(x = ~year, y = ~std_res)
p %>% add_markers(x = ~log(pred), y = ~std_res)
p %>% add_markers(x = ~survey, y = ~std_res)

## Process error residuals
pe <- data.frame(year = landings$year,
                 species = landings$species,
                 stock = landings$stock,
                 pe = est$log_B_res,
                 pe_lwr = lwr$log_B_res,
                 pe_upr = upr$log_B_res)
p <- plot_ly()
for (nm in unique(pe$species)) {
    p <- p %>% add_fit(data = pe[pe$species == nm, ],
                       x = ~year, color = ~species, colors = viridis::viridis(100),
                       yline = ~pe,
                       legendgroup = ~species)
}
p

## Correlation in pe
pe_wide <- tidyr::spread(pe[, c("year", "species", "pe")], species, pe)
pe_wide$year <- NULL
plot(pe_wide)
pe_wide <- as.matrix(pe_wide)
pe_wide[is.infinite(pe_wide)] <- NA
cor_mat <- cor(pe_wide, use = "na.or.complete")
plot_ly(x = rownames(cor_mat), y = rownames(cor_mat), z = ~cor_mat) %>%
    add_heatmap()
corrplot::corrplot.mixed(cor_mat, diag = "n", lower = "ellipse", upper = "number")

p <- plot_ly()
for (nm in levels(index$survey)) {
    p <- p %>% add_fit(data = index[index$survey == nm, ],
                       x = ~year, color = ~survey, colors = viridis::viridis(100),
                       ymin = ~pred_lwr, ymax = ~pred_upr,
                       ymarker = ~index, yline = ~pred,
                       legendgroup = ~species) %>%
        layout(yaxis = list(type = "log", title = "index"))
}
p


biomass <- data.frame(year = landings$year,
                      species = landings$species,
                      stock = landings$stock,
                      B = exp(est$log_B_vec),
                      B_lwr = exp(lwr$log_B_vec),
                      B_upr = exp(upr$log_B_vec))

p <- plot_ly()
for (nm in unique(biomass$species)) {
    p <- p %>% add_fit(data = biomass[biomass$species == nm, ],
                       x = ~year, color = ~species, colors = viridis::viridis(100),
                       yline = ~B, ymin = ~B_lwr, ymax = ~B_upr,
                       legendgroup = ~species)
}
p
p %>% layout(yaxis = list(type = "log"))





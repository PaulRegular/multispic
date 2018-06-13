
library(units)
library(plotly)
library(TMB)

unique(MSP::landings$species)

index <- MSP::index
landings <- MSP::landings

## Subset the data
# landings <- landings[landings$year >= min(index$year), ]
start_year <- 1985
index <- index[index$year >= start_year, ]
landings <- landings[landings$year >= start_year, ]
sub_sp <- c("Cod", "Yellowtail", "Witch", "Plaice")
index <- index[index$species %in% sub_sp, ]
landings <- landings[landings$species %in% sub_sp, ]

## Set-up indices for TMB
landings$species <- factor(landings$species)
landings$y <- factor(landings$year)
landings$sy <- factor(paste0(landings$species, "-", landings$year))
landings <- landings[order(landings$sy), ]
index$sy <- factor(paste0(index$species, "-", index$year), levels = levels(landings$sy))
index$ss <- factor(paste0(index$species, "-", index$survey))
index$species <- factor(index$species)

index %>%
    plot_ly() %>%
    add_lines(x = ~year, y = ~index, color = ~ss,
              colors = viridis::viridis(100)) %>%
    layout(yaxis = list(type = "log", title = "index"))

landings %>%
    plot_ly() %>%
    add_lines(x = ~year, y = ~landings, color = ~species,
              colors = viridis::viridis(100)) %>%
    layout(yaxis = list(title = "Landings"))


dat <- list(L = as.numeric(landings$landings),
            L_species = as.numeric(landings$species) - 1,
            L_year = as.numeric(landings$y) - 1,
            I = as.numeric(index$index),
            I_species = as.numeric(index$species) - 1,
            I_survey = as.numeric(index$ss) - 1,
            I_sy = as.numeric(index$sy) - 1,
            min_P = 0.0001)
par <- list(log_P = rep(0, nrow(landings)),
            log_sd_P = rep(0, nlevels(landings$species)),
            log_K = rep(5, nlevels(landings$species)),
            log_mu_r = rep(log(0.5), nlevels(landings$species)),
            log_sd_r = 0,
            log_res_r = rep(0, nrow(landings)),
            log_m = rep(log(2), nlevels(landings$species)),
            log_q = rep(0, nlevels(index$ss)),
            log_sd_I = rep(0, nlevels(index$ss)))
map <- list(log_m = factor(rep(NA, nlevels(landings$species))),
            log_res_r = landings$y)

obj <- MakeADFun(dat, par, map = map, random = c("log_P", "log_res_r"), DLL = "MSP")
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(eval.max = 1000, iter.max = 1000))
sd_rep <- sdreport(obj)
obj$report()
exp(opt$par)
sd_rep

est <- split(unname(sd_rep$value), names(sd_rep$value))
sd <- split(sd_rep$sd, names(sd_rep$value))
lwr <- split(unname(sd_rep$value) - 1.96 * sd_rep$sd, names(sd_rep$val))
upr <- split(unname(sd_rep$value) + 1.96 * sd_rep$sd, names(sd_rep$val))

index$pred <- exp(est$log_pred_I)
index$pred_lwr <- exp(lwr$log_pred_I)
index$pred_upr <- exp(upr$log_pred_I)

p <- plot_ly()
for (nm in levels(index$ss)) {
    p <- p %>% add_fit(data = index[index$ss == nm, ],
                       x = ~year, color = ~ss, colors = viridis::viridis(100),
                       ymin = ~pred_lwr, ymax = ~pred_upr,
                       ymarker = ~index, yline = ~pred,
                       legendgroup = ~ss) %>%
        layout(yaxis = list(type = "log", title = "index"))
}
p


biomass <- data.frame(year = landings$year,
                      species = landings$species,
                      B = exp(est$log_B),
                      B_lwr = exp(lwr$log_B),
                      B_upr = exp(upr$log_B))

p <- plot_ly()
for (nm in unique(biomass$species)) {
    p <- p %>% add_fit(data = biomass[biomass$species == nm, ],
                       x = ~year, color = ~species, colors = viridis::viridis(100),
                       yline = ~B, ymin = ~B_lwr, ymax = ~B_upr,
                       legendgroup = ~species)
}
p


pe <- data.frame(year = landings$year,
                 species = landings$species,
                 pe = exp(est$log_res_P),
                 pe_lwr = exp(lwr$log_res_P),
                 pe_upr = exp(upr$log_res_P))

p <- plot_ly()
for (nm in unique(pe$species)) {
    p <- p %>% add_fit(data = pe[pe$species == nm, ],
                       x = ~year, color = ~species, colors = viridis::viridis(100),
                       yline = ~pe, ymin = ~pe_lwr, ymax = ~pe_upr,
                       legendgroup = ~species)
}
p



r <- data.frame(year = landings$year,
                 species = landings$species,
                 r = exp(est$log_r),
                 r_lwr = exp(lwr$log_r),
                 r_upr = exp(upr$log_r))
p <- plot_ly()
for (nm in unique(r$species)) {
    p <- p %>% add_fit(data = r[r$species == nm, ],
                       x = ~year, color = ~species, colors = viridis::viridis(100),
                       yline = ~r, ymin = ~r_lwr, ymax = ~r_upr,
                       legendgroup = ~species)
}
p


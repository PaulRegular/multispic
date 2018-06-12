
library(units)
library(plotly)
library(TMB)

unique(MSP::landings$species)

index <- MSP::index[MSP::index$species == "Witch", ]
landings <- MSP::landings[MSP::landings$species == "Witch" &
                              MSP::landings$year >= min(index$year), ]

plot_ly() %>%
    add_lines(data = landings, x = ~year, y = ~landings, name = "Landings") %>%
    add_lines(data = index, x = ~year, y = ~index, color = ~survey) %>%
    layout(yaxis = list(title = "inputs"))

dat <- list(L = as.numeric(landings$landings),
            I = as.numeric(index$index),
            I_year = index$year - min(landings$year),
            I_survey = as.numeric(factor(index$survey)) - 1,
            min_P = 0.0001)
par <- list(log_P = rep(0, nrow(landings)),
            log_sd_P = 0,
            log_K = 5,
            log_r = 0,
            log_m = log(2),
            log_q = rep(0, length(unique(index$survey))),
            log_sd_I = rep(0, length(unique(index$survey))))
map <- list(log_P = factor(c(0, seq(nrow(landings) - 1))),
            log_m = factor(NA))

obj <- MakeADFun(dat, par, map = map, random = "log_P", DLL = "MSP")
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

index %>%
    plot_ly(x = ~year, color = ~survey) %>%
    add_ribbons(ymin = ~pred_lwr, ymax = ~pred_upr) %>%
    add_markers(y = ~index) %>%
    add_lines(y = ~pred) %>%
    layout(yaxis = list(type = "log"))

biomass <- data.frame(year = landings$year, B = exp(est$log_B),
                      B_lwr = exp(lwr$log_B),
                      B_upr = exp(upr$log_B))
biomass %>%
    plot_ly(x = ~year) %>%
    add_ribbons(ymin = ~B_lwr, ymax = ~B_upr, color = I("grey"), name = "95% CI") %>%
    add_lines(y = ~B, color = I("darkgrey"), name = "Est")


pe <- data.frame(year = landings$year,
                 pe = exp(est$log_res_P),
                 pe_lwr = exp(lwr$log_res_P),
                 pe_upr = exp(upr$log_res_P))
pe %>%
    plot_ly(x = ~year) %>%
    add_ribbons(ymin = ~pe_lwr, ymax = ~pe_upr, color = I("grey"), name = "95% CI") %>%
    add_lines(y = ~pe, color = I("darkgrey"), name = "Est")



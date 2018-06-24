
## TODO:
## - Get Rstrap calls from DE's and run yourself in data-raw
## - Make year 0 species specific, and start model when landings start
## - Look into calculating one-step ahead residuals
## - Continue to think of ways to share information across species

library(units)
library(plotly)
library(TMB)

unique(multispic::landings$species)

index <- multispic::index
landings <- multispic::landings

## Subset the data
sub_sp <- unique(multispic::landings$species)
start_year <- 1985 # restricted by hake and skate landings
index <- index[index$year >= start_year & index$species %in% sub_sp, ]
landings <- landings[landings$year >= start_year & landings$species %in% sub_sp, ]

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
    layout(yaxis = list(type = "log", title = "Landings"))


## Set up correlations to estimate
cor_nms <- t(combn(levels(landings$species), 2))
cor_ind <- t(combn(seq(nlevels(landings$species)) - 1, 2))

dat <- list(L = as.numeric(landings$landings),
            L_species = as.numeric(landings$species) - 1,
            L_year = as.numeric(landings$y) - 1,
            I = as.numeric(index$index),
            I_species = as.numeric(index$species) - 1,
            I_survey = as.numeric(index$ss) - 1,
            I_sy = as.numeric(index$sy) - 1,
            min_P = 0.0001,
            cor_ind = cor_ind)
par <- list(log_P = rep(0, nrow(landings)),
            log_sd_P = rep(0, nlevels(landings$species)),
            logit_cor = rnorm(nrow(cor_ind)), # rep(0, nrow(cor_ind)),
            log_K = rep(5, nlevels(landings$species)),
            log_r = rep(0, nlevels(landings$species)),
            log_m = rep(log(2), nlevels(landings$species)),
            log_q = rep(0, nlevels(index$ss)),
            log_sd_I = rep(0, nlevels(index$ss)))
map <- list(log_m = factor(rep(NA, nlevels(landings$species))),
            logit_cor = factor(rep(NA, length(par$logit_cor))))

obj <- MakeADFun(dat, par, map = map, random = "log_P", DLL = "multispic")
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(eval.max = 1000, iter.max = 1000))
sd_rep <- sdreport(obj)
obj$report()
exp(opt$par)
sd_rep

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
head(index)
p <- index %>% plot_ly(color = ~ss, colors = viridis::viridis(100))
p %>% add_markers(x = ~year, y = ~std_res)
p %>% add_markers(x = ~log(pred), y = ~std_res)

## Process error residuals
pe <- data.frame(year = landings$year,
                 species = landings$species,
                 pe = est$log_res_P,
                 pe_lwr = lwr$log_res_P,
                 pe_upr = upr$log_res_P)
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
p %>% layout(yaxis = list(type = "log"))


prop <- data.frame(year = landings$year,
                   species = landings$species,
                   P = exp(est$log_P),
                   P_lwr = exp(lwr$log_P),
                   P_upr = exp(upr$log_P))

p <- plot_ly()
for (nm in unique(prop$species)) {
    p <- p %>% add_fit(data = prop[prop$species == nm, ],
                       x = ~year, color = ~species, colors = viridis::viridis(100),
                       yline = ~P, ymin = ~P_lwr, ymax = ~P_upr,
                       legendgroup = ~species)
}
p
p %>% layout(yaxis = list(type = "log"))


par_est <- split(opt$par, names(opt$par))
lapply(names(par_est), function(nm) hist(exp(par_est[[nm]]), xlab = nm, main = nm))
## consider random effect on process error, observation error
## and q (especially when you replace mwpt with total estimates)



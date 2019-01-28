
library(Rstrap)

load("data-raw/Rstrap_runs/converted_set_details_2018-09-19.Rdata")

## Note: index converted to kt

## 3LNO Plaice
out <- strat.fun(setdet = con.setdet, program = "strat2",
                 data.series = c("Engel","Campelen"), species = 889,
                 survey.year = c(1985:2005, 2007:2014, 2016),
                 season = "spring", NAFOdiv = c("3L", "3N", "3O"), export = NULL)
spring <- out$strat2$biomass$summary[, c("survey.year", "total")]
spring <- data.frame(species = "Plaice", stock = "3LNO", survey = "Spring",
                     year = spring$survey.year, index = spring$total / 1000000)

out <- strat.fun(setdet = con.setdet, program = "strat2",
                 data.series = c("Engel","Campelen"), species = 889,
                 survey.year = c(1990:2003, 2005:2013, 2015:2017),
                 season = "fall", NAFOdiv = c("3L", "3N", "3O"), export = NULL)
fall <- out$strat2$biomass$summary[, c("survey.year", "total")]
fall <- data.frame(species = "Plaice", stock = "3LNO", survey = "Fall",
                   year = fall$survey.year, index = fall$total / 1000000)

plaice <- rbind(spring, fall)

write.csv(plaice, file = "data-raw/Rstrap_runs/plaice.csv", row.names = FALSE)




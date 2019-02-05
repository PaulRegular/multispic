
library(Rstrap)

load("data-raw/Rstrap_runs/set_details_2018-11-26.Rdata")

## Note: index converted to kt

some_strat <- function(years, season, series, series_name) {
    out <- strat.fun(setdet = setdet, program = "strat2",
                     data.series = series, species = 447,
                     survey.year = years,
                     season = season, NAFOdiv = c("3L", "3N", "3O"), export = NULL)
    tab <- out$strat2$biomass$summary[, c("survey.year", "total")]
    tab <- data.frame(species = "Hake", stock = "3LNO", survey = series_name,
                         year = tab$survey.year, index = tab$total / 1000000)
    tab
}

## 3LNO Hake
spring_yankee <- some_strat(c(1975:1982), "spring", "Yankee", "Spring Yankee")
spring_engel <- some_strat(c(1984:1995), "spring", "Engel", "Spring Engel")
spring_campelen <- some_strat(c(1996:2016), "spring", "Campelen", "Spring Campelen")
fall_engel <- some_strat(c(1990:1994), "fall", "Engel", "Fall Engel")
fall_campelen <- some_strat(c(1995:2013, 2015:2017), "fall", "Campelen", "Fall Campelen")

skate <- rbind(spring_yankee, spring_engel, spring_campelen, fall_engel, fall_campelen)

write.csv(skate, file = "data-raw/Rstrap_runs/hake.csv", row.names = FALSE)




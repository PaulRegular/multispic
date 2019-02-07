
library(Rstrap)

load("data-raw/Rstrap_runs/set_details_2018-11-26.Rdata")

## Note: index converted to kt

one_strat <- function(years, season, series, NAFOdiv, species, species_name) {

    out <- strat.fun(setdet = setdet, program = "strat2",
                     data.series = series, species = species,
                     survey.year = years,
                     season = season, NAFOdiv = NAFOdiv,
                     export = NULL, plot.results = FALSE)
    tab <- out$strat2$biomass$summary[, c("survey.year", "total")]

    stock <- paste0("3", paste0(gsub("3", "", NAFOdiv), collapse = ""))

    tab <- data.frame(species = species_name, stock = stock, gear = series,
                      season = Hmisc::capitalize(season),
                      year = tab$survey.year, index = tab$total / 1000000)
    tab
}

stack_strat <- function(NAFOdiv, species, species_name) {

    spring_yankee <- one_strat(c(1975:1982), "spring", "Yankee", NAFOdiv, species, species_name)
    spring_engel <- one_strat(c(1984:1995), "spring", "Engel", NAFOdiv, species, species_name)
    spring_campelen <- one_strat(c(1996:2016), "spring", "Campelen", NAFOdiv, species, species_name)
    fall_engel <- one_strat(c(1990:1994), "fall", "Engel", NAFOdiv, species, species_name)
    fall_campelen <- one_strat(c(1995:2013, 2015:2017), "fall", "Campelen", NAFOdiv, species, species_name)

    rbind(spring_yankee, spring_engel, spring_campelen, fall_engel, fall_campelen)

}


yellowtail <- stack_strat(c("3L", "3N", "3O"), 891, "Yellowtail")
witch <- stack_strat(c("3N", "3O"), 890, "Witch")
cod <- stack_strat(c("3N", "3O"), 438, "Cod")
plaice <- stack_strat(c("3L", "3N", "3O"), 889, "Plaice")
redfish <- stack_strat(c("3L", "3N", "3O"), 794, "Redfish")
skate <- stack_strat(c("3L", "3N", "3O"), 90, "Skate")
hake <- stack_strat(c("3L", "3N", "3O"), 447, "Hake")

index <- rbind(yellowtail, witch, cod, plaice, redfish, skate, hake)

write.csv(index, file = "data-raw/index.csv", row.names = FALSE)



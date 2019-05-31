
library(Rstrap)

load("data-raw/Rstrap_runs/set_details_2019-02-15.Rdata")

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

    ## Keep years were > 60 strat were covered
    spring_yankee <- one_strat(c(1979, 1980, 1982), "spring", "Yankee", NAFOdiv, species, species_name)
    spring_engel <- one_strat(c(1985:1995), "spring", "Engel", NAFOdiv, species, species_name)
    spring_campelen <- one_strat(c(1996:2005, 2007:2016, 2018),
                                 "spring", "Campelen", NAFOdiv, species, species_name)
    fall_engel <- one_strat(c(1990:1994), "fall", "Engel", NAFOdiv, species, species_name)
    fall_campelen <- one_strat(c(1995:2013, 2015:2018), "fall", "Campelen", NAFOdiv, species, species_name)

    rbind(spring_yankee, spring_engel, spring_campelen, fall_engel, fall_campelen)

}

divs <- c("3L", "3N", "3O")
yellowtail <- stack_strat(divs, 891, "Yellowtail")
witch <- stack_strat(divs, 890, "Witch")
cod <- stack_strat(divs, 438, "Cod")
plaice <- stack_strat(divs, 889, "Plaice")
redfish <- stack_strat(divs, 794, "Redfish")
skate <- stack_strat(divs, 90, "Skate")
hake <- stack_strat(divs, 447, "Hake")
haddock <- stack_strat(divs, 441, "Haddock")

index <- rbind(yellowtail, witch, cod, plaice, redfish, skate, hake, haddock)

write.csv(index, file = "data-raw/index.csv", row.names = FALSE)



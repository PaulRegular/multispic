
library(Rstrap)
library(dplyr)
library(plotly)

load("data-raw/Rstrap_runs/set_details_2019-02-15.Rdata")

## Note: index converted to kt

## Keep years where 90% of the core area is covered
## Core area = strat that have been sampled for more than 20 years
setrec <- setdet %>%
    filter(rec == 5, NAFOdiv %in% c("3L", "3N", "3O"),
           which.survey == "multispecies")

strat_freq <- setrec %>%
    group_by(survey.year, season) %>%
    distinct(NAFOdiv, strat) %>%
    group_by(season, NAFOdiv, strat) %>%
    summarise(n_years_sampled = n()) %>%
    ungroup()

strat_freq %>%
    plot_ly(x = ~factor(strat), y = ~n_years_sampled, color = ~season) %>%
    add_bars() %>%
    layout(xaxis = list(title = "Strat"),
           yaxis = list(title = "Years sampled"))

core_strat <- strat_freq %>%
    filter(n_years_sampled > 20) %>%
    .$strat %>% sort() %>% unique()

coverage <- setrec %>%
    filter(strat %in% core_strat) %>%
    group_by(season, survey.year) %>%
    distinct(NAFOdiv, strat, strat.area) %>%
    summarise(n_div = length(unique(NAFOdiv)),
              n_strat = length(unique(strat)),
              area_covered = sum(strat.area)) %>%
    ungroup() %>%
    mutate(percent_covered = (area_covered / max(area_covered)) * 100) %>%
    as.data.frame()

coverage %>%
    plot_ly(x = ~survey.year, y = ~percent_covered, color = ~season) %>%
    add_bars() %>%
    layout(xaxis = list(title = "Year"),
           yaxis = list(title = "Percent coverage"))

fall_years <- coverage %>%
    filter(percent_covered >= 90, season == "fall") %>%
    .$survey.year

spring_years <- coverage %>%
    filter(percent_covered >= 90, season == "spring") %>%
    .$survey.year


## Functions

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
    spring_yankee <- one_strat(spring_years[spring_years <= 1982],
                               "spring", "Yankee", NAFOdiv, species, species_name)
    spring_engel <- one_strat(spring_years[spring_years >= 1982 & spring_years <= 1995],
                              "spring", "Engel", NAFOdiv, species, species_name)
    spring_campelen <- one_strat(spring_years[spring_years >= 1996],
                                 "spring", "Campelen", NAFOdiv, species, species_name)
    fall_engel <- one_strat(fall_years[fall_years >= 1982 & fall_years <= 1995],
                            "fall", "Engel", NAFOdiv, species, species_name)
    fall_campelen <- one_strat(fall_years[fall_years >= 1996],
                               "fall", "Campelen", NAFOdiv, species, species_name)

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



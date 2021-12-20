
library(Rstrap)
library(dplyr)
library(plotly)

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


## Common species = those that have been consistently sampled in
## core strata for > 40 years

totals <- setdet %>%
    filter(!is.na(spec) & strat %in% core_strat,
           survey.year %in% unique(c(fall_years, spring_years))) %>%
    group_by(survey.year, spec, common.name) %>%
    summarise(n_sets = n(),
              total_n = sum(number, na.rm = TRUE),
              total_weight = sum(weight, na.rm = TRUE)) %>%
    group_by(survey.year) %>%
    mutate(ranking = dense_rank(desc(total_weight))) %>%
    arrange(ranking)

# length(unique(setdet$survey.year))
common_spp <- totals %>%
    filter(n_sets > 1) %>%
    group_by(spec, common.name) %>%
    summarise(n_cases = n()) %>%
    arrange(-n_cases) %>%
    filter(n_cases > 40)
data.frame(common_spp)

## Tidy up common name
x <- common_spp$common.name
x <- gsub(",", ", ", x)
x <- sub("(^.*),\\s(.*$)","\\2 \\1", x)
x <- tolower(x)
x <- tools::toTitleCase(x)
x <- gsub(" \\(Ns\\)| \\(Common)| \\(Marlin|\\(Monkfish\\)", "", x)
x[x == "Turbot"] <- "Greenland Halibut"
x[x == "Halibut (Atlantic)"] <- "Atlantic Halibut"
x[x == " Deep Water Redfish"] <- "Redfish"
x[x == "Offshore Sand Launce"] <- "Sand Lance"
x[x == "Lanternfishes"] <- "Lanternfish"
names(x) <- common_spp$spec
common_spp <- x

totals$common.name <- common_spp[as.character(totals$spec)]

totals %>%
    filter(!is.na(common.name)) %>%
    plot_ly(x = ~survey.year, y = ~total_weight, color = ~common.name) %>%
    add_bars() %>%
    layout(barmode = "stack")


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



## Generate index across common species
index <- lapply(seq_along(common_spp), function(i) {
    stack_strat(c("3L", "3N", "3O"),
                as.numeric(names(common_spp[i])),
                unname(common_spp[i]))
})
index <- do.call(rbind, index)

write.csv(index, file = "data-raw/index.csv", row.names = FALSE)


# ## Export some yellowtail and witch data
# yt <- strat.fun(setdet = setdet, program = "strat2",
#                 data.series = "Campelen", species = 891,
#                 survey.year = 1996:2018,
#                 season = "spring", NAFOdiv = c("3L", "3N", "3O"),
#                 export = NULL, plot.results = FALSE)
# yt_setdet <- yt$raw.data$set.details
# write.csv(yt_setdet, file = "data-raw/Rstrap_runs/exports/yellowtail_3LNO.csv",
#           row.names = FALSE)
#
# wi <- strat.fun(setdet = setdet, program = "strat2",
#                 data.series = "Campelen", species = 890,
#                 survey.year = 1996:2018,
#                 season = "spring", NAFOdiv = c("3N", "3O"),
#                 export = NULL, plot.results = FALSE)
# wi_setdet <- wi$raw.data$set.details
# write.csv(wi_setdet, file = "data-raw/Rstrap_runs/exports/witch_3NO.csv",
#           row.names = FALSE)
#



library(Rstrap)
library(dplyr)
library(plotly)

## Prep data

## Keep years where 80% of the core area is covered
## Core area = strat that have been sampled for more than 20 years

## Also identify common species
## Common species = those that have been consistently sampled in
## core strata for > 30 years

region_data <- function(divisions) {

    setrec <- Rstrap::setdet %>%
        filter(rec == 5, NAFOdiv %in% divisions,
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
        filter(percent_covered >= 80, season == "fall") %>%
        .$survey.year

    spring_years <- coverage %>%
        filter(percent_covered >= 80, season == "spring") %>%
        .$survey.year


    totals <- Rstrap::setdet %>%
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
        filter(n_cases > 30)
    data.frame(common_spp)

    ## Tidy up common name
    x <- common_spp$common.name
    x <- gsub(",", ", ", x)
    x <- sub("(^.*),\\s(.*$)","\\2 \\1", x)
    x <- tolower(x)
    x <- tools::toTitleCase(x)
    x <- gsub(" \\(Ns\\)| \\(Common)| \\(Marlin|\\(Monkfish\\)|\\(Marinus\\)|\\(Round\\) ", "", x)
    x[x == "Turbot"] <- "Greenland Halibut"
    x[x == "Halibut (Atlantic)"] <- "Atlantic Halibut"
    x[x == " Deep Water Redfish"] <- "Redfish spp."
    x[x == "Offshore Sand Launce"] <- "Sand Lance"
    x[x == "Lanternfishes"] <- "Lanternfish"
    x[x == "Common Angler"] <- "Monkfish"
    x[x == "Scyphozoan (Jellyfish)"] <- "Jellyfish"
    x[x == "Cephalopod Rossi.(B.a.sq"] <- "Rossi spp."
    x[x == "Octopus   Octo."] <- "Octopus spp."
    names(x) <- common_spp$spec
    common_spp <- x

    ## Aggregate Redfish, Wolffish, and Skate

    redfish_spp <- common_spp[grep("Redfish", common_spp)]
    redfish <- "Redfish spp."
    names(redfish) <- paste0(names(redfish_spp), collapse = "")

    wolf_spp <- common_spp[grep("Wolffish", common_spp)]
    wolf <- "Wolffish spp."
    names(wolf) <- paste0(names(wolf_spp), collapse = "")

    skate_spp <- common_spp[grep("Skate", common_spp)]
    skate <- "Skate spp."
    names(skate) <- paste0(names(skate_spp), collapse = "")

    common_spp <- common_spp[!common_spp %in% c(redfish_spp, wolf_spp, skate_spp)]
    common_spp <- c(common_spp, redfish, wolf, skate)

    totals$spec[totals$spec %in% as.numeric(names(redfish_spp))] <- as.numeric(names(redfish))
    totals$spec[totals$spec %in% as.numeric(names(wolf_spp))] <- as.numeric(names(wolf))
    totals$spec[totals$spec %in% as.numeric(names(skate_spp))] <- as.numeric(names(skate))
    totals$common.name <- common_spp[as.character(totals$spec)]

    totals %>%
        filter(!is.na(common.name)) %>%
        plot_ly(x = ~survey.year, y = ~total_weight, color = ~common.name) %>%
        add_bars() %>%
        layout(barmode = "stack")

    ## Replace species specific codes in setdet with redfish, wolffish, and skate spp. "codes"
    setdet$spec[setdet$spec %in% as.numeric(names(redfish_spp))] <- as.numeric(names(redfish))
    setdet$spec[setdet$spec %in% as.numeric(names(wolf_spp))] <- as.numeric(names(wolf))
    setdet$spec[setdet$spec %in% as.numeric(names(skate_spp))] <- as.numeric(names(skate))

    setdet <- Rstrap::setdet %>%
        filter(strat %in% core_strat,
               (!is.na(season) & season == "spring" & survey.year %in% spring_years) |
               (!is.na(season) & season == "fall" & survey.year %in% fall_years),
               (rec == 5 | spec %in% as.numeric(names(common_spp))))

    setdet$common.name <- common_spp[as.character(setdet$spec)]

    setdet

}



## Note: index converted to kt

one_strat <- function(setdet, years, season, series, NAFOdiv, species, species_name, region) {

    out <- strat.fun(setdet = setdet, program = "strat2",
                     data.series = series, species = species,
                     survey.year = years,
                     season = season, NAFOdiv = NAFOdiv,
                     export = NULL, plot.results = FALSE)
    tab <- out$strat2$biomass$summary[, c("survey.year", "total")]

    tab <- data.frame(species = species_name, region = region, gear = series,
                      season = Hmisc::capitalize(season),
                      year = tab$survey.year, index = tab$total / 1000000)
    tab
}

stack_strat <- function(setdet, NAFOdiv, species, species_name, region) {

    spring_years <- sort(unique(setdet$survey.year[setdet$season == "spring"]))
    fall_years <- sort(unique(setdet$survey.year[setdet$season == "fall"]))

    spring_yankee <- try(one_strat(setdet, spring_years[spring_years <= 1982],
                                   "spring", "Yankee", NAFOdiv, species, species_name, region))
    spring_engel <- try(one_strat(setdet, spring_years[spring_years >= 1982 & spring_years <= 1995],
                                  "spring", "Engel", NAFOdiv, species, species_name, region))
    spring_campelen <- try(one_strat(setdet, spring_years[spring_years >= 1996],
                                     "spring", "Campelen", NAFOdiv, species, species_name, region))

    fall_yankee <- try(one_strat(setdet, fall_years[fall_years <= 1982],
                                 "fall", "Yankee", NAFOdiv, species, species_name, region))
    fall_engel <- try(one_strat(setdet, fall_years[fall_years >= 1982 & fall_years <= 1995],
                                "fall", "Engel", NAFOdiv, species, species_name, region))
    fall_campelen <- try(one_strat(setdet, fall_years[fall_years >= 1996],
                                   "fall", "Campelen", NAFOdiv, species, species_name, region))

    dat <- list()
    for (nm in c("spring_yankee", "spring_engel", "spring_campelen",
           "fall_yankee", "fall_engel", "fall_campelen")) {
        if (class(get(nm)) != "try-error") dat[[nm]] <- get(nm)
    }

    do.call(rbind, dat)

}



library(Rstrap)
library(dplyr)
library(plotly)

## Prep data

region_data <- function(divisions) {

    sub_setdet <- Rstrap::setdet %>%
        filter(NAFOdiv %in% divisions)

    ## Treat winter as spring survey for 3Ps as is done in the assessment
    ## Gear change roughly corresponds to change in survey timing
    sub_setdet$season[!is.na(sub_setdet$season) &
                          sub_setdet$NAFOdiv == "3P" &
                          sub_setdet$season == "winter"] <- "spring"

    ## Tidy up common name
    sp_names <- unique(na.omit(sub_setdet[, c("spec", "common.name")]))
    x <- sp_names$common.name
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
    names(x) <- sp_names$spec
    common_name <- x

    ## Aggregate Redfish, Wolffish, and Skate species
    ## because these species are not differentiated in the landings

    redfish_spp <- common_name[grep("Redfish", common_name)]
    redfish <- "Redfish spp."
    names(redfish) <- 794.999

    wolf_spp <- common_name[grep("Wolffish", common_name)]
    wolf <- "Wolffish spp."
    names(wolf) <- 700.999

    skate_spp <- common_name[grep("Skate", common_name)]
    skate <- "Skate spp."
    names(skate) <- 90.999

    common_name <- common_name[!common_name %in% c(redfish_spp, wolf_spp, skate_spp)]
    common_name <- c(common_name, redfish, wolf, skate)

    ## Replace species specific codes in setdet with redfish, wolffish, and skate spp. "codes"
    sub_setdet$spec[sub_setdet$spec %in% as.numeric(names(redfish_spp))] <- as.numeric(names(redfish))
    sub_setdet$spec[sub_setdet$spec %in% as.numeric(names(wolf_spp))] <- as.numeric(names(wolf))
    sub_setdet$spec[sub_setdet$spec %in% as.numeric(names(skate_spp))] <- as.numeric(names(skate))

    sub_setdet$common.name <- common_name[as.character(sub_setdet$spec)]

    sub_setdet

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


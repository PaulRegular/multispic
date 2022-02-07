
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



## Under one_strat
## If more than 5 years of data are present for a survey, the core strata are identified (those
## covered in > 90% of unique survey years) and years where > 10% of the biomass were potentially
## missing because of poor survey coverage are excluded.
## Note: index converted to kt

one_strat <- function(setdet, season, series, NAFOdiv, species, species_name, region) {

    ## Identify years with data across all divisions
    ind <- setdet$which.survey == "multispecies" & setdet$rec == 5 &
        setdet$NAFOdiv %in% NAFOdiv & setdet$data.series %in% series & setdet$season == season
    year_div_table <- table(setdet$survey.year[ind], setdet$NAFOdiv[ind])
    year_ind <- rowSums(year_div_table > 0) == length(NAFOdiv)
    years <- as.numeric(names(year_ind[year_ind]))

    if (length(years) >= 5) {

        ## First run analysis using all available data
        all_strat <- strat.fun(setdet = setdet, program = "strat2",
                               data.series = series, species = species,
                               survey.year = years,
                               season = season, NAFOdiv = NAFOdiv,
                               export = NULL, plot.results = FALSE)

        all_setdet <- all_strat$raw.data$set.details
        all_means <- all_strat$strat2$biomass$details

        if (all(all_means$totals == 0)) {

            tab <- data.frame(species = species_name, region = region, gear = series,
                              season = Hmisc::capitalize(season),
                              year = years, index = 0, cv = 0)

        } else {

            ## Calculate percent biomass in each strata using grand totals
            ## Limit to strata covered across more than 90% of survey years (i.e., keep core strata)
            strat_percents <- all_means %>%
                group_by(strat) %>%
                summarise(strat_total = sum(totals), n = n(), percent_years = n() / length(years)) %>%
                ungroup() %>%
                mutate(grand_total = sum(strat_total), mean_percent = strat_total / sum(strat_total)) %>%
                filter(percent_years > 0.9)
            keep_strat <- strat_percents$strat

            ## Identify years where more than 10% of the biomass was likely to be missing because of poor coverage
            coverage <- table(all_setdet$strat, all_setdet$survey.year) == 0
            coverage <- coverage[rownames(coverage) %in% as.character(keep_strat), ]
            percents <- replicate(ncol(coverage), matrix(strat_percents$mean_percent, ncol = 1), simplify = TRUE)
            coverage <- coverage * percents
            percent_missing <- colSums(coverage)
            drop_years <- percent_missing > 0.1

            keep_years <- as.numeric(names(percent_missing)[!drop_years])

            sub_strat <- strat.fun(setdet = setdet, program = "strat2",
                                   data.series = series, species = species,
                                   survey.year = keep_years, strat = keep_strat,
                                   season = season, NAFOdiv = NAFOdiv,
                                   export = NULL, plot.results = FALSE)

            tab <- sub_strat$strat2$biomass$summary[, c("survey.year", "var", "mean", "total")]

            tab <- data.frame(species = species_name, region = region, gear = series,
                              season = Hmisc::capitalize(season),
                              year = tab$survey.year, index = tab$total / 1000000,
                              cv = sqrt(tab$var) / tab$mean)
            tab

        }

    } else {

        warning(paste0("There were insufficient survey data to produce the following index: ", species_name, " - ", region, " - ", season, " - ", series))
        tab <- NULL

    }

    tab

}

stack_strat <- function(setdet, NAFOdiv, species, species_name, region) {

    spring_yankee <- one_strat(setdet, "spring", "Yankee", NAFOdiv, species, species_name, region)
    spring_engel <- one_strat(setdet, "spring", "Engel", NAFOdiv, species, species_name, region)
    spring_campelen <- one_strat(setdet, "spring", "Campelen", NAFOdiv, species, species_name, region)

    fall_yankee <- one_strat(setdet, "fall", "Yankee", NAFOdiv, species, species_name, region)
    fall_engel <- one_strat(setdet, "fall", "Engel", NAFOdiv, species, species_name, region)
    fall_campelen <- one_strat(setdet, "fall", "Campelen", NAFOdiv, species, species_name, region)

    dat <- rbind(spring_yankee, spring_engel, spring_campelen, fall_yankee, fall_engel, fall_campelen)
    dat

}


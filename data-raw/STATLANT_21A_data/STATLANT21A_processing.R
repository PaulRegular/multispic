
library(dplyr)

## Import STATLANT 21A landings data from https://www.nafo.int/Data/STATLANT-21A
## Note the landings data are in tonnes
landings <- read.csv("data-raw/STATLANT_21A_data/STATLANT21A_Extraction.csv")
names(landings) <- c("year", "country", "division", "species", "landings")

## Filter to focal area
landings <- landings %>%
    filter(division %in% c("2J", "3K", "3L", "3N", "3O", "3P")) %>%
    mutate(region = case_when(
        division %in% c("2J", "3K") ~ "2J3K",
        division %in% c("3L", "3N", "3O") ~ "3LNO",
        division %in% c("3P") ~ "3P",
        TRUE ~ "ERROR"
    ))

## Demersal fish species with all time reported landings of > 1000 tonnes
totals <- landings %>%
    group_by(species) %>%
    summarise(n_years = length(unique(year)), total = sum(landings)) %>%
    arrange(-total)
data.frame(totals)

keep_sp <- c("ATLANTIC COD - COD" = "Atlantic Cod",
             "ATLANTIC REDFISHES (NS) - RED" = "Redfish spp.",
             "AMERICAN PLAICE - PLA" = "American Plaice",
             "GREENLAND HALIBUT - GHL" = "Greenland Halibut",
             "YELLOWTAIL FLOUNDER - YEL" = "Yellowtail Flounder",
             "WITCH FLOUNDER - WIT" = "Witch Flounder",
             "SKATES (NS) - SKA" = "Skate spp.",
             "ROUNDNOSE GRENADIER - RNG" = "Roundnose Grenadier",
             "HADDOCK - HAD" = "Haddock",
             "WHITE HAKE - HKW" = "White Hake",
             "WOLFFISHES (CATFISH) (NS) - CAT" = "Wolffish spp.",
             "ROUGHHEAD GRENADIER - RHG" = "Roughhead Grenadier",
             "ATLANTIC HALIBUT - HAL" = "Atlantic Halibut",
             "AMERICAN ANGLER - ANG" = "Monkfish",
             "RED HAKE - HKR" = "Red Hake",
             "WINTER FLOUNDER - FLW" = "Winter Flounder", # may not be present in the survey data b/c inshore
             "SILVER HAKE - HKS" = "Silver Hake",
             "BEAKED REDFISH(DEEP-WATER) - REB" = "Redfish spp.",
             "ATLANTIC WOLFFISH - CAA" = "Wolffish spp.",
             "THORNY SKATE (STARRY RAY) - RJR" = "Skate spp.",
             "GOLDEN REDFISH - REG" = "Redfish spp.",
             "NORTHERN WOLFFISH - CAB" = "Wolffish spp.",
             "SPOTTED WOLFFISH - CAS" = "Wolffish spp.")

## TODO: Consider options for dealing with general FINFISHES (NS), GROUNDFISHES (NS), FLATFISHES (NS) categories

## Subset to demersal fish species with cumulative catch of > 1000 tonnes
landings <- landings[landings$species %in% names(keep_sp), ]
landings$species <- keep_sp[as.character(landings$species)] # simplify names

## Sum by region and species
landings <- landings %>%
    group_by(year, region, species) %>%
    summarise(landings = sum(landings) / 1000) %>%
    as.data.frame()

## Fill missing values with zero
grd <- expand.grid(year = unique(landings$year),
                   region = unique(landings$region),
                   species = unique(landings$species))
landings <- merge(landings, grd, by = c("year", "region", "species"), all = TRUE)
landings$landings[is.na(landings$landings)] <- 0

## Export
landings <- landings[, c("year", "species", "region", "landings")]
landings <- landings[order(landings$year, landings$region, landings$species), ]
names(landings) <- c("Year", "Species", "Region", "Landings (kt)")
write.csv(landings, file = "data-raw/landings.csv", row.names = FALSE)




library(dplyr)

## Import STATLANT 21A landings data from 3LNO
## Note the landings data are in tonnes
landings <- read.csv("data-raw/STATLANT_21A_data/STATLANT21A_Extraction.csv")
names(landings) <- c("year", "country", "division", "species", "landings")

## Groundfish species with all time reported landings of > 1000 tonnes and > 10 years
totals <- landings %>%
    group_by(species) %>%
    summarise(n_years = length(unique(year)), total = sum(landings)) %>%
    arrange(-total)
data.frame(totals)

keep_sp <- c("ATLANTIC COD - COD" = "Atlantic Cod",
             "AMERICAN PLAICE - PLA" = "American Plaice",
             "ATLANTIC REDFISHES (NS) - RED" = "Redfish spp.",
             "CAPELIN - CAP" = "Capelin",
             "YELLOWTAIL FLOUNDER - YEL" = "Yellowtail Flounder",
             "GREENLAND HALIBUT - GHL" = "Greenland Halibut",
             "SKATES (NS) - SKA" = "Skate spp.",
             "HADDOCK - HAD" = "Haddock",
             "WITCH FLOUNDER - WIT" = "Witch Flounder",
             "ATLANTIC HERRING - HER" = "Atlantic Herring",
             "WHITE HAKE - HKW" = "White Hake",
             "WOLFFISHES (CATFISH) (NS) - CAT" = "Wolffish spp.",
             "ROUGHHEAD GRENADIER - RHG" = "Roughhead Grenadier",
             "ATLANTIC HALIBUT - HAL" = "Atlantic Halibut",
             "ROUNDNOSE GRENADIER - RNG" = "Roundnose Grenadier",
             "AMERICAN ANGLER - ANG" = "Monkfish",
             "SILVER HAKE - HKS" = "Silver Hake")

landings <- landings[landings$species %in% names(keep_sp), ]

## Sum to 3LNO
landings <- landings %>%
    group_by(species, year) %>%
    summarise(landings = sum(landings) / 1000) %>%
    as.data.frame()

## Simplify names
landings$species <- keep_sp[as.character(landings$species)]
landings$stock <- "3LNO"

## Export
landings <- landings[, c("species", "stock", "year", "landings")]
names(landings) <- c("Species", "Stock", "Year", "Landings (kt)")
write.csv(landings, file = "data-raw/landings.csv", row.names = FALSE)



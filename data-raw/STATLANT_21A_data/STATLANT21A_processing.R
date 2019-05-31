
library(dplyr)

## Import STATLANT 21A landings data from 3LNO
## Note the landinigs data are in tonnes
landings <- read.csv("data-raw/STATLANT_21A_data/STATLANT21A_Extraction.csv", skip = 10)
names(landings) <- c("year", "division", "species", "abbrev", "landings")

sub_landings <- landings[landings$year == 1984, ]
head(sub_landings[order(sub_landings$landings, decreasing = TRUE), ], 100)

keep_sp <- c("Atlantic Cod" = "Cod",
             "American Plaice" = "Plaice",
             "Atlantic Redfishes (ns)" = "Redfish",
             "Yellowtail Flounder" = "Yellowtail",
             "White Hake" = "Hake",
             "Witch Flounder" = "Witch",
             "Haddock" = "Haddock",
             "Skates (ns)" = "Skate")

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



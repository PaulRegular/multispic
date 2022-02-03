
index <- read.csv("data-raw/index.csv")

landings <- read.csv("data-raw/landings.csv")
names(landings) <- c("year", "region", "species", "landings")

covariates <- read.csv("data-raw/covariates.csv")
landings <- merge(landings, covariates, by = "year", all.x = TRUE)

## Impose common end year
terminal_year <- min(max(index$year), max(landings$year))
index <- index[index$year <= terminal_year, ]
landings <- landings[landings$year <= terminal_year, ]

save(index, file = "data/index.RData")
save(landings, file = "data/landings.RData")


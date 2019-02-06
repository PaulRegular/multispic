
library(units)

index <- read.csv("data-raw/index.csv")
index$index <- set_units(index$index, "kt")
save(index, file = "data/index.RData")

landings <- read.csv("data-raw/landings.csv")
names(landings) <- c("species", "stock", "year", "landings")
landings$landings <- set_units(landings$landings, "kt")
save(landings, file = "data/landings.RData")

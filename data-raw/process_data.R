
library(units)

index <- read.csv("data-raw/index.csv")
names(index) <- c("year", "survey", "index")
index$index <- set_units(index$index, "kt")
save(index, file = "data/index.RData")

landings <- read.csv("data-raw/landings.csv")
names(landings) <- c("year", "tac", "landings")
landings$landings <- set_units(landings$landings, "kt")
save(landings, file = "data/landings.RData")

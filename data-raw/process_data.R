
index <- read.csv("data-raw/index.csv")
save(index, file = "data/index.RData")

landings <- read.csv("data-raw/landings.csv")
names(landings) <- c("year", "region", "species", "landings")
save(landings, file = "data/landings.RData")

covariates <- read.csv("data-raw/covariates.csv")
save(covariates, file = "data/covariates.RData")

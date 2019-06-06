
index <- read.csv("data-raw/index.csv")
save(index, file = "data/index.RData")

landings <- read.csv("data-raw/landings.csv")
names(landings) <- c("species", "stock", "year", "landings")
save(landings, file = "data/landings.RData")

covariates <- read.csv("data-raw/covariates.csv")
covariates <- covariates[order(covariates$year), ]
save(covariates, file = "data/covariates.RData")

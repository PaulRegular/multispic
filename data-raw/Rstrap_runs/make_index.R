
source("data-raw/Rstrap_runs/index_helpers.R")

regions <- list("2J3K" = c("2J", "3K"),
                "3LNO" = c("3L", "3N", "3O"),
                "3Ps" = c("3P"))

index <- lapply(seq_along(regions), function(i) {

    setdet <- region_data(regions[[i]])

    common_spp <- unique(na.omit(setdet[, c("spec", "common.name")]))

    ## Generate index across common species
    region_index <- lapply(seq(nrow(common_spp)), function(j) {
        stack_strat(setdet,
                    regions[[i]],
                    common_spp$spec[j],
                    common_spp$common.name[j],
                    names(regions[i]))
    })
    region_index <- do.call(rbind, region_index)
    rownames(region_index) <- NULL

    region_index

})
index <- do.call(rbind, index)

index[index$index == 0, ]
table(index$species[index$index == 0])

## Zeros are an occasional problem.
## They are more commonly observed for rare and non-commercial species.
## Hard to know if they are true zeros or simply missing data.

write.csv(index, file = "data-raw/index.csv", row.names = FALSE)


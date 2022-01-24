
source("data-raw/index/index_helpers.R")

## Load landings file to limit to focal species (see data-raw/landings/make_landings.R)
landings <- read.csv("data-raw/landings.csv")
keep_sp <- sort(unique(landings$Species))


regions <- list("2J3K" = c("2J", "3K"),
                "3LNO" = c("3L", "3N", "3O"),
                "3Ps" = c("3P"))

index <- lapply(seq_along(regions), function(i) {

    sub_setdet <- region_data(regions[[i]])

    focal_sp <- unique(na.omit(sub_setdet[sub_setdet$common.name %in% keep_sp, c("spec", "common.name")]))

    ## Generate index across common species
    region_index <- lapply(seq(nrow(focal_sp)), function(j) {
        stack_strat(sub_setdet,
                    regions[[i]],
                    focal_sp$spec[j],
                    focal_sp$common.name[j],
                    names(regions[i]))
    })
    region_index <- do.call(rbind, region_index)
    rownames(region_index) <- NULL

    ## Drop spring Yankee-Engel series for Grenadier because deep water coverage was poor
    ## (Kept by the algorithm but index strata are not representative because)
    if (names(regions[i]) == "3LNO") {
        region_index <- region_index %>%
            filter(!species %in% c("Roughhead Grenadier", "Roundnose Grenadier"))
    }

    ## Drop species with a partial series (less than 30 years of data)
    ## Zeros are also an occasional problem for species like Roundnose Grenadier, especially
    ## early in the time series. It is hard to know if these are true zeros or simply because
    ## they were not sampled.
    # region_index %>% group_by(species) %>% summarise(n_years = length(unique(year)))
    region_index <- region_index %>%
        filter(index > 0) %>%
        group_by(species) %>%
        mutate(n_years = length(unique(year))) %>%
        filter(n_years > 30)

    region_index

})
index <- do.call(rbind, index)

## Early values are not applicable to many species. Manually define start of each series.
plot_ly(index, x = ~year, y = ~index, color = ~species, frame = ~region) %>% add_markers()
index <- index %>%
    filter(region == "2J3K" & year >= 1978 |
               region %in% c("3LNO", "3Ps") & year >= 1983)
plot_ly(index, x = ~year, y = ~index, color = ~species, frame = ~region) %>% add_markers()

index$n_years <- NULL

write.csv(index, file = "data-raw/index.csv", row.names = FALSE)


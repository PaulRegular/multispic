
source("data-raw/index/index_helpers.R")

## Load landings file to limit to focal species (see data-raw/landings/make_landings.R)
landings <- read.csv("data-raw/landings.csv")
keep_sp <- sort(unique(landings$Species))

## Load strat areas file to calculate total area in each region
library(vroom)
sa <- vroom::vroom_fwf("data-raw/index/stratum_areas.pc4",
                       col_positions = fwf_widths(widths = c(3, 4, 4, 3),
                                                  col_names = c("strat", "area", "max_depth", "div")),
                       col_types = cols(strat = "i", area = "i", max_depth = "c", div = "c"))
sa$div[sa$div == "3PS"] <- "3P"


regions <- list("2J3K" = c("2J", "3K"),
                "3LNO" = c("3L", "3N", "3O"),
                "3Ps" = c("3P"))

index <- lapply(seq_along(regions), function(i) {

    sub_setdet <- region_data(regions[[i]])

    total_area <- sum(sa$area[sa$div %in% regions[[i]]])

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

    region_index$total_area <- total_area
    region_index$coverage <- region_index$survey_area / region_index$total_area

    region_index

})
index <- do.call(rbind, index)

index %>%
    group_by(gear, species, season) %>%
    plot_ly(x = ~year, y = ~index, color = ~species, frame = ~region,
            legendgroup = ~species, hoverinfo = "x+y+text", hovertext = ~paste(season, gear)) %>%
    add_lines() %>%
    layout(title = "Pre-processed indices")

## Manually unify the start year for 2J3KL
index <- index %>%
    filter(region == "2J3K" & year >= 1978 | region != "2J3K")

## Drop Grenadier as there has been poor coverage of this deep water species
index <- index %>%
    filter(!species %in% c("Roughhead Grenadier", "Roundnose Grenadier"))

## Drop species with a partial series
## (less than 30 years of data or if species isn't present across each series)
## Zeros are also an occasional problem, especially early in the time series.
## It is hard to know if these are true zeros or simply because they were not sampled.
# index %>% group_by(region, species) %>% summarise(n_years = length(unique(year))) %>% as.data.frame()
index <- index %>%
    filter(index > 0) %>%
    group_by(region, species) %>%
    mutate(n_years = length(unique(year)),
           n_series = length(unique(gear))) %>%
    filter(n_years > 30,
           (region == "2J3K" & n_series == 2 |
                region != "2J3K" & n_series == 3)) %>%
    as.data.frame()
index$n_years <- NULL
index$n_series <- NULL

p <- index %>%
    group_by(gear, species, season) %>%
    plot_ly(x = ~year, y = ~index, color = ~species, frame = ~region,
            legendgroup = ~species, hoverinfo = "x+y+text", hovertext = ~paste(season, gear)) %>%
    add_lines()  %>%
    layout(title = "Processed indices")
p
p %>% layout(yaxis = list(type = "log"))

write.csv(index, file = "data-raw/index.csv", row.names = FALSE)


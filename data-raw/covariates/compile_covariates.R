
library(data.table)
library(dplyr)

## Calculate annual NAO index
monthly_nao <- fread("data-raw/covariates/monthly_NAO.txt")
names(monthly_nao) <- c("year", "month", "index")

nao <- monthly_nao %>%
    filter(month %in% c(12, 1, 2)) %>%
    group_by(year) %>%
    summarise(nao = mean(index))

## Import core cil data
core_cil <- fread("data-raw/covariates/core_cil.csv")

## Import composite environmental index
cei <- fread("data-raw/covariates/cei.csv")
names(cei) <- c("year", "cei")

## Import timing of ice retreat data
tice <- fread("data-raw/covariates/tice.csv") %>%
    select(year, tice)

## Merge and export
covariates <- nao %>%
    full_join(core_cil, by = "year") %>%
    full_join(cei, by = "year") %>%
    full_join(tice, by = "year") %>%
    arrange(year)

write.csv(covariates, file = "data-raw/covariates.csv", row.names = FALSE)

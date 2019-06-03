
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

## Merge and export
covariates <- merge(nao, core_cil, by = c("year"), all = TRUE)
write.csv(covariates, file = "data-raw/covariates.csv", row.names = FALSE)

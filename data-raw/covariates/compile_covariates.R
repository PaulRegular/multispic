
library(data.table)
library(dplyr)

## Import monthly NAO index
monthly_nao <- fread("https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii.table")
names(monthly_nao) <- c("year", month.abb)

## Lag Dec to pair with Jan for winter NAO index
monthly_nao$Dec <- c(NA, monthly_nao$Dec[-1])

## Calculate seasonal NAO indices
nao <- data.frame(year = monthly_nao$year,
                  winter_nao = rowMeans(monthly_nao[, c("Dec", "Jan", "Feb")]),
                  spring_nao = rowMeans(monthly_nao[, c("Mar", "Apr", "May")]),
                  summer_nao = rowMeans(monthly_nao[, c("Jun", "Jul", "Aug")]),
                  fall_nao = rowMeans(monthly_nao[, c("Sep", "Oct", "Nov")]))
covariates <- nao

write.csv(covariates, file = "data-raw/covariates.csv", row.names = FALSE)

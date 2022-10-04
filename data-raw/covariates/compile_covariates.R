
library(data.table)
library(dplyr)

## NAO ---------------------------------------------------------------------------------------------

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


## NL climate index --------------------------------------------------------------------------------

nlci <- fread("data-raw/covariates/climate/NL_climate_index.csv") %>% as.data.frame()
names(nlci) <- c("year", "nlci")

## The cumulative sum of the NLCI represents the effect of consecutive warm or cold years
## Hypothesis is that higher trophic level species may be responding to the cumulative effects
## of sequential warm or cold years
nlci$cumsum_nlci <- replace(nlci$nlci, is.na(nlci$nlci), 0) |> cumsum()
plot(cumsum_nlci ~ year, data = nlci, type = "l")

## Likewise, a moving average may capture a similar effect; select 3 years as recruitment to the
## survey is belated by 3 years or more for several of the focal species
nlci$ma_nlci <- zoo::rollmean(nlci$nlci, 3, na.pad = TRUE, align = "right")
plot(ma_nlci ~ year, data = nlci, type = "l", xlim = c(1978, 2021))


## Merge and export --------------------------------------------------------------------------------

covariates <- merge(nao, nlci, by = "year", all = TRUE)

## Add 1991 breakpoint covariate - this is when the system change.
## Hypothesis is that energy flow to higher trophic levels have been hampered since
## the collapse of mid-trophic level forage fish such as capelin.
covariates$shift <- ifelse(covariates$year < 1991, "pre-1991", "post-1991")
covariates$shift <- factor(covariates$shift, levels = c("pre-1991", "post-1991"))

write.csv(covariates, file = "data-raw/covariates.csv", row.names = FALSE)



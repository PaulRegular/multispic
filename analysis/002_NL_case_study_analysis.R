
source("analysis/001_NL_case_study_helpers.R")

for (r in c("2J3K", "3LNO", "3Ps")) {

    full <- nl_multispic(region = r, species = NULL)

    no_nao <- update(full, pe_covariates = NULL)
    just_nao <- update(full, species_cor = "none", temporal_cor = "none")
    one_species_cor <- update(no_nao, species_cor = "one")
    no_species_cor <- update(one_species_cor, species_cor = "none")
    no_temporal_cor <- update(no_species_cor, temporal_cor = "none")

    full$loo <- run_loo(full)
    no_nao$loo<- run_loo(no_nao)
    just_nao$loo <- run_loo(just_nao)
    one_species_cor$loo  <- run_loo(one_species_cor)
    no_species_cor$loo  <- run_loo(no_species_cor)
    no_temporal_cor$loo  <- run_loo(no_temporal_cor)

    fits <- mget(c("full", "no_nao", "just_nao", "one_species_cor",
                   "no_species_cor", "no_temporal_cor"))

    saveRDS(fits, file = paste0("analysis/exports/fits_", r, ".rds"))

}




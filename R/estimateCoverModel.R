## Cover-presence probability of each ecoregionGroup x speciesCode, for makeSpeciesEcoregion().
##
## Rows where the species is present in every pixel (coverPres == coverNum) make the binomial
## model fail, so they are left out of the fit and given probability 1. When every row is like
## that -- a study area with a single tree species, for instance -- nothing is left to fit:
## fitting the empty table stopped with "contrasts can be applied only to factors with 2 or
## more levels" (ELF 3.1.2, 2026-09-10).
##
## `fitCover` is called with the rows to fit and returns the `statsModel()` list.
## Returns a list: `modelCover` for makeSpeciesEcoregion() (the fitted list, or one probability
## per row of `cohortDataShort` when any row was left out), `model` (the fitted list, `NULL` when
## nothing was fitted) and `cohortDataShort` (with a `pred` column when probabilities were set).
estimateCoverModel <- function(cohortDataShort, fitCover) {
  cdsWh <- cohortDataShort$coverPres == cohortDataShort$coverNum
  cds <- Copy(cohortDataShort)
  cds <- cds[!cdsWh]

  if (!NROW(cds)) {
    return(list(modelCover = rep(1, NROW(cohortDataShort)), model = NULL,
                cohortDataShort = cohortDataShort))
  }

  model <- fitCover(cds)
  modelCover <- model
  if (isTRUE(any(cdsWh))) {
    cds[, pred := fitted(model$mod, response = "response")]
    cohortDataShort <- cds[, -c("coverPres", "coverNum")][cohortDataShort,
                                                          on = c("ecoregionGroup", "speciesCode"), nomatch = NA]
    cohortDataShort[is.na(pred), pred := 1]
    modelCover <- cohortDataShort$pred
  }
  list(modelCover = modelCover, model = model, cohortDataShort = cohortDataShort)
}

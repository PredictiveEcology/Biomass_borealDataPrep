## maxB can be no lower than the biomass a species was actually observed at in that group.
##
## makeSpeciesEcoregion() predicts maxB from the biomass model at `cover = 100` and old age, and
## clamps a negative prediction to 0. A model that under-predicts -- typically a group where a
## species is rare -- then says the species cannot grow where it demonstrably grows. The largest
## observed cohort biomass is a lower bound on what the species reaches there, so use it.
##
## `cohortData` is the estimation data (all rows, not the model's subsample), with `B`,
## `speciesCode` and `ecoregionGroup`. maxANPP is recomputed for raised rows with the rule
## makeSpeciesEcoregion() uses, `asInteger(maxB / 30)`. The number of raised rows is returned
## as attribute "nRaised".
floorMaxBAtObserved <- function(speciesEcoregion, cohortData) {
  obs <- cohortData[!is.na(B), list(obsMaxB = max(B)),
                    by = list(eg = as.character(ecoregionGroup), sc = as.character(speciesCode))]
  se <- data.table::copy(data.table::as.data.table(speciesEcoregion))
  se[, `:=`(eg = as.character(ecoregionGroup), sc = as.character(speciesCode))]
  se[obs, obsMaxB := i.obsMaxB, on = c("eg", "sc")]
  raise <- which(!is.na(se$obsMaxB) & se$maxB < se$obsMaxB)
  if (length(raise)) {
    newMaxB <- se$obsMaxB[raise]
    if (is.integer(se$maxB)) newMaxB <- as.integer(newMaxB)
    data.table::set(se, raise, "maxB", newMaxB)
    if ("maxANPP" %in% names(se)) {
      newANPP <- LandR::asInteger(se$maxB[raise] / 30)
      if (!is.integer(se$maxANPP)) newANPP <- as.numeric(newANPP)
      data.table::set(se, raise, "maxANPP", newANPP)
    }
  }
  se[, c("eg", "sc", "obsMaxB") := NULL]
  data.table::setattr(se, "nRaised", length(raise))
  se[]
}

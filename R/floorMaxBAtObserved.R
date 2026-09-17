## A clamped maxB is replaced by the largest biomass the species was observed at in its group.
##
## makeSpeciesEcoregion() predicts maxB from the biomass model at `cover = 100` and old age, and
## clamps a negative prediction to 0. Where a species is rare the fit can go negative, and the
## clamp then says the species cannot grow where it demonstrably grows. Only those rows -- maxB
## of 0 for a species that was observed in the group -- are changed.
##
## Deliberately NOT a floor on every row: the largest of thousands of observed cohorts is an
## outlier statistic, and flooring every row at it overrode genuine fits (black spruce upland
## 5,541 -> 10,800 on a 60 km test window). A fitted value is kept whatever it is.
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
  raise <- which(!is.na(se$obsMaxB) & se$obsMaxB > 0 & se$maxB <= 0)
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

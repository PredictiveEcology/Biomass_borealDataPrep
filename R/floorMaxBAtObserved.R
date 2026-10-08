## A clamped maxB is replaced by a high quantile of the biomass the species was observed at in
## its group -- the 95th percentile by default, not the maximum, so one freak cohort cannot set
## the ceiling for a whole ecoregion x species.
##
## makeSpeciesEcoregion() predicts maxB from the biomass model at `cover = 100` and old age, and
## clamps a negative prediction to 0. Where a species is rare the fit can go negative, and the
## clamp then says the species cannot grow where it demonstrably grows. Only those rows -- maxB
## of 0 for a species that was observed in the group -- are changed.
##
## The quantile follows what the module and LandR already do with observed extremes rather than
## inventing a rule: longevity is `quantile(age, 0.99) * 1.3` (Biomass_borealDataPrep.R:724), and
## LandR's maxB estimation summarises its predictions at the 0.05/0.25/0.5/0.75/0.95 quantiles
## (LandR/R/maxBestimation.R:959). `type = 7`, R's default, so this is the ordinary sample
## quantile; with few observations it sits close to the maximum, which is the sensible limit.
##
## Deliberately NOT a floor on every row: an upper quantile of thousands of observed cohorts is
## still an extreme statistic, and flooring every row at it overrode genuine fits (black spruce
## upland 5,541 -> 10,800 on a 60 km test window). A fitted value is kept whatever it is.
##
## `cohortData` is the estimation data (all rows, not the model's subsample), with `B`,
## `speciesCode` and `ecoregionGroup`. maxANPP is recomputed for raised rows with the rule
## makeSpeciesEcoregion() uses, `asInteger(maxB / 30)`. The number of raised rows is returned
## as attribute "nRaised".
floorMaxBAtObserved <- function(speciesEcoregion, cohortData, probs = 0.95) {
  stopifnot(length(probs) == 1, probs > 0, probs <= 1)
  obs <- cohortData[!is.na(B),
                    list(obsMaxB = as.numeric(stats::quantile(B, probs = probs, names = FALSE))),
                    by = list(eg = as.character(ecoregionGroup), sc = as.character(speciesCode))]
  se <- data.table::copy(data.table::as.data.table(speciesEcoregion))
  se[, `:=`(eg = as.character(ecoregionGroup), sc = as.character(speciesCode))]
  se[obs, obsMaxB := i.obsMaxB, on = c("eg", "sc")]
  raise <- which(!is.na(se$obsMaxB) & se$obsMaxB > 0 & se$maxB <= 0)
  if (length(raise)) {
    newMaxB <- se$obsMaxB[raise]
    if (is.integer(se$maxB)) newMaxB <- LandR::asInteger(newMaxB)
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

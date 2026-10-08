## Land-cover class of a pixel from its species cover, for the class-240 fill.
##
## 210, 220 and 230 are statements about composition -- a pixel cannot be pine-dominated and
## carry a broadleaf class -- so species cover can decide them. The reverse is not true: a map
## that says "broadleaf" cannot say which species. That is why the fill takes composition from
## here rather than from another year's land cover or from a neighbouring pixel's class.
##
## Site is NOT decided here. Whether a pixel is wet (80/81) is a property of the ground, and
## no amount of species cover expresses it; the caller keeps site from the land-cover record
## and uses this only to choose the composition class within an upland pixel.
##
## The rule has the same shape as `LandR::vegTypeMapGenerator(mixedType = 2)` -- conifer share
## of tree cover at or above the threshold is conifer-leading, at or below its complement is
## deciduous-leading, anything between is mixed -- and the tests verify that agreement directly
## at 0.8. The DEFAULT here is 0.75 rather than 0.8, because this function answers a different
## question: which NTEMS legend code does this pixel carry, not which LandR vegetation type is
## it. See the provenance below.
##
## NTEMS itself drew the line at 75%, not 80%, and on a different quantity. Hermosilla et al.
## (2018) build the VLCE classes on the NFI land-cover scheme via the EOSD legend (Wulder &
## Nelson 2003), which defines them as:
##   Coniferous  coniferous trees are 75% or more of total basal area
##   Broadleaf   broadleaf trees are 75% or more of total basal area
##   Mixed Wood  neither coniferous nor broadleaf accounts for 75% or more of total basal area
## and the NFI Photo Plot Data Dictionary (v5.2 and v6.1, identical wording) says the same for
## TC/TB/TM on total tree VOLUME, with >= 10% crown cover to be treed at all. Crown closure in
## EOSD is a separate axis (dense >60%, open 26-60%, sparse 10-25%) -- it sets density, never
## composition.
##
## So the default is 0.75, matching the product whose codes we write. It deliberately does NOT
## track the module's `vegLeadingProportion` parameter: that one is shared across modules and
## `paramCheckOtherMods(sim, "vegLeadingProportion", ifSetButDifferent = "error")` makes it an
## error for two modules to disagree on it, so it cannot be bent to 0.75 without changing what
## every other module means by leading vegetation. Composition on crown cover is still not the
## basal area / volume the definition is written on; `deciduousCoverWeight` closes part of
## that gap by making deciduous cover conifer-equivalent before the ratio is taken.
##
## `coverDT` has one row per pixel and one numeric column per species, in percent. Column names
## are species codes in the `sppEquivCol` convention, optionally prefixed "cover." as
## `makePixelTable()` leaves them. Species absent from `sppEquiv` are ignored.
##
## Returns an integer vector, one per row: 210, 220, 230, or NA where the pixel has no tree
## cover at all and therefore nothing to infer from -- those pixels are left for the caller's
## fallback rather than being forced into a class.
## `deciduousCoverWeight` makes deciduous cover comparable with coniferous before the ratio
## is taken: hardwoods carry a much wider canopy than softwoods, so the same crown cover is
## less stem for a deciduous species. `partitionBiomass()` applies it the same way -- deciduous
## cover is MULTIPLIED by the weight (`cover * c(1, x)[decid + 1]`); equal 50/50 cover then
## yields conifer/deciduous biomass of 1.1878 = 1/0.8419 at the old default. Dividing would
## invert it. The weight is not bounded above by 1: the module estimates it per study area
## (`R/deciduousCoverWeight.R`) and some landscapes come back above 1.
speciesLeadingClass <- function(coverDT, sppEquiv, sppEquivCol = "LandR",
                                vegLeadingProportion = getOption(
                                  "NTEMS.mixedwoodProp",
                                  getOption("LandR.lccLeadingProportion", 0.75)
                                ),
                                deciduousCoverWeight = 1) {
  stopifnot(vegLeadingProportion > 0.5, vegLeadingProportion <= 1,
            deciduousCoverWeight > 0)

  sppCols <- grep("^cover\\.", colnames(coverDT), value = TRUE)
  bare <- if (length(sppCols)) sub("^cover\\.", "", sppCols) else colnames(coverDT)
  if (!length(sppCols)) {
    sppCols <- bare[bare %in% sppEquiv[[sppEquivCol]]]
    bare <- sppCols
  }
  keep <- bare %in% sppEquiv[[sppEquivCol]]
  sppCols <- sppCols[keep]
  bare <- bare[keep]
  if (!length(sppCols)) {
    return(rep(NA_integer_, NROW(coverDT)))
  }

  ## Conifer or deciduous, from the species table rather than a hard-coded list.
  type <- sppEquiv[["Type"]][match(bare, sppEquiv[[sppEquivCol]])]
  isConifer <- !is.na(type) & type == "Conifer"
  isDecid <- !is.na(type) & type == "Deciduous"

  m <- as.matrix(coverDT[, sppCols, with = FALSE])
  m[is.na(m)] <- 0
  conCover <- if (any(isConifer)) rowSums(m[, isConifer, drop = FALSE]) else rep(0, nrow(m))
  decCover <- if (any(isDecid)) rowSums(m[, isDecid, drop = FALSE]) else rep(0, nrow(m))
  decCover <- decCover * deciduousCoverWeight   ## conifer-equivalent; see the note above
  total <- conCover + decCover

  out <- rep(NA_integer_, nrow(m))
  has <- total > 0
  frac <- conCover[has] / total[has]

  ## Both thresholds are inclusive, matching `vegTypeMapGenerator()`: measured against it, a
  ## conifer share of exactly 0.80 is conifer-leading and exactly 0.20 is deciduous-leading,
  ## with mixed strictly between. Compare with a tolerance rather than exactly, because
  ## `1 - 0.8` is 0.19999999999999996 while `20/100` is 0.20000000000000001, so an exact
  ## `<=` sent a true 20/80 split to mixed instead of deciduous.
  eps <- sqrt(.Machine$double.eps)
  out[has] <- data.table::fifelse(
    frac >= vegLeadingProportion - eps, 210L,
    data.table::fifelse(frac <= (1 - vegLeadingProportion) + eps, 220L, 230L)
  )
  out
}

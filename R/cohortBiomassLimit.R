## The largest cohort B, in g/m2, that createBiomass_coreInputs() can make from a raw biomass map
## whose maximum is `maxRawB` (t/ha * 100), for its check that every cohort is within the map.
##
## A cohort's B is rounded half up twice: makePixelTable() sets the pixel's `totalBiomass` to
## `asInteger(t/ha * 100)`, partitionBiomass() splits it among the pixel's cohorts, and each share
## is rounded to `pixelGroupBiomassClass`. Both roundings are monotone and a share is at most
## `totalBiomass`, so no cohort exceeds `maxRawB` rounded the same way. The check used
## `round(maxRawB, -2)`, one half-even rounding to a fixed 100, which can be lower:
## NorthernAlbertaUplands' maximum, 329.4970703 t/ha, gives `totalBiomass` 32950 and, in a
## single-species pixel, a cohort of 33000 against a limit of 32900.
## Young cohorts re-estimated inside fire perimeters are not rounded to the class, and
## youngBiomassLimit() caps them at `maxRawB`, which can be above the rounded value; hence `max()`.
cohortBiomassLimit <- function(maxRawB, pixelGroupBiomassClass) {
  rounded <- LandR::asInteger(LandR::asInteger(maxRawB) / pixelGroupBiomassClass) *
    pixelGroupBiomassClass
  max(rounded, maxRawB)
}

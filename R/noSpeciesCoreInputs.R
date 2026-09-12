#' Core inputs for a study area with no tree species
#'
#' Everything `createBiomass_coreInputs()` does after the input rasters are checked estimates
#' tree traits (cover, biomass, `speciesEcoregion`) from species cover layers. With zero species
#' layers there is nothing to estimate, so the outputs that define trees are empty. Fire models
#' fitted on non-forest fuel classes only still need them to exist, with the columns and geometry
#' their consumers address by name.
#'
#' @param rasterToMatch the `SpatRaster` defining the geometry of the simulation outputs
#'
#' @return a list with `cohortData` (a 0-row `data.table`) and `pixelGroupMap` (a `SpatRaster`
#'   with `rasterToMatch`'s geometry, 0 -- no cohorts -- wherever `rasterToMatch` has data)
noSpeciesCoreInputs <- function(rasterToMatch) {
  cohortData <- data.table(
    pixelGroup = integer(0),
    ecoregionGroup = factor(),
    speciesCode = factor(),
    age = integer(0),
    B = integer(0),
    totalBiomass = integer(0)
  )

  ## 0 is how the with-species path codes pixels that have data but no cohorts
  pixelGroupMap <- terra::rast(rasterToMatch)
  pixelGroupMap[!is.na(terra::values(mat = FALSE, rasterToMatch))] <- 0L

  list(cohortData = cohortData, pixelGroupMap = pixelGroupMap)
}

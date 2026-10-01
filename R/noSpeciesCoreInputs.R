#' Core inputs for a study area with no tree species
#'
#' Everything `createBiomass_coreInputs()` does after the input rasters are checked estimates
#' tree traits (cover, biomass, `speciesEcoregion`) from species cover layers. With zero species
#' layers there is nothing to estimate, so the outputs that define trees are empty. Fire models
#' fitted on non-forest fuel classes only still need them to exist, with the columns and geometry
#' their consumers address by name.
#'
#' The tables are empty but well-formed rather than absent, because `suppliedElsewhere()` reports
#' TRUE for anything this module DECLARES in `createsOutput()`, whether or not it was assigned.
#' Biomass_regeneration's `.inputObjects` fallbacks for `species` and `sufficientLight` are
#' therefore suppressed, and it reads them unguarded (Biomass_regeneration.R:249-252). A
#' promised-but-unset output is worse than no output at all. `species` is not built here: the
#' caller gets it from the same `prepSpeciesTable()` call as the with-species path, which
#' returns the 0-row table with the full column set when `sppEquiv` has no rows.
#'
#' @param rasterToMatch the `SpatRaster` defining the geometry of the simulation outputs
#'
#' @return a named list of the outputs `createBiomass_coreInputs()` would otherwise build:
#'   `cohortData`, `pixelGroupMap`, `sufficientLight` and `speciesEcoregion`
#'
#' @importFrom data.table data.table
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

  ## not species-dependent: the same LANDIS-test table the with-species path uses
  sufficientLight <- data.frame(
    speciesshadetolerance = 1:5,
    X0 = c(rep(1, 4), 0),
    X1 = c(0, rep(1, 3), 0),
    X2 = c(0, 0, rep(1, 3)),
    X3 = c(rep(0, 3), rep(1, 2)),
    X4 = c(rep(0, 4), 1),
    X5 = c(rep(0, 4), 1)
  )

  ## the columns makeSpeciesEcoregion() returns, with no rows
  speciesEcoregion <- data.table(
    ecoregionGroup = factor(),
    speciesCode = factor(),
    establishprob = numeric(0),
    maxB = integer(0),
    maxANPP = integer(0),
    year = integer(0)
  )

  list(cohortData = cohortData, pixelGroupMap = pixelGroupMap,
       sufficientLight = sufficientLight, speciesEcoregion = speciesEcoregion)
}

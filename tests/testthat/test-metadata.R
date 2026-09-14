## The module's metadata is its public contract: a project using this module binds
## to these object names and classes. Renaming or retyping one breaks every caller,
## which is exactly the class of change the raster -> terra migration makes, so it is
## worth asserting here rather than discovering downstream.
##
## When a change is deliberate, update this file in the same commit and bump the
## module version to match: removed, renamed or retyped is a MAJOR bump.

test_that("module metadata parses", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_type(md, "list")
  expect_identical(md$name, moduleName)
})

test_that("inputs are the expected names and classes", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  inputs <- stats::setNames(md$inputObjects$objectClass, md$inputObjects$objectName)
  expect_identical(
    inputs[order(names(inputs))],
    c(cloudFolderID              = "character",
      columnsForPixelGroups      = "character",
      ecoregionLayer             = "sf",
      ecoregionRst               = "SpatRaster",
      firePerimeters             = "SpatRaster",
      imputedPixID               = "integer",
      rasterToMatch              = "SpatRaster",
      rasterToMatch_biomassParam = "SpatRaster",
      rawBiomassMap              = "SpatRaster",
      rstLCC                     = "SpatRaster",
      speciesLayers              = "SpatRaster",
      speciesTable               = "data.table",
      sppColorVect               = "character",
      sppEquiv                   = "data.table",
      sppNameVector              = "character",
      standAgeMap                = "SpatRaster",
      studyArea                  = "sf",
      studyArea_biomassParam     = "sf")
  )
})

test_that("outputs are the expected names and classes", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  outputs <- stats::setNames(md$outputObjects$objectClass, md$outputObjects$objectName)
  expect_identical(
    outputs[order(names(outputs))],
    c(biomassMap       = "SpatRaster",
      cohortData       = "data.table",
      ecoregion        = "data.table",
      ecoregionMap     = "SpatRaster",
      firePerimeters   = "SpatRaster",
      imputedPixID     = "integer",
      minRelativeB     = "data.frame",
      modelBiomass     = "data.frame",
      modelCover       = "data.frame",
      pixelFateDT      = "data.table",
      pixelGroupMap    = "SpatRaster",
      rstLCC           = "SpatRaster",
      species          = "data.table",
      speciesEcoregion = "data.table",
      speciesLayers    = "SpatRaster",
      standAgeMap      = "SpatRaster",
      sufficientLight  = "data.frame")
  )
})

test_that("parameters are the expected names", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_identical(
    sort(md$parameters$paramName),
    sort(c(".plotInitialTime", ".plotInterval", ".plots", ".saveInitialTime",
           ".saveInterval", ".seed", ".sslVerify", ".studyAreaName", ".useCache",
           ".useCacheArgs", ".useCloud", "adjustAgeAndLongevity", "biomassModel",
           "coverModel", "coverPctToBiomassPctModel", "dataSource", "dataYear",
           "deciduousCoverDiscount", "earliestFireYear", "ecoregionLayerField",
           "exportModels", "fitDeciduousCoverDiscount", "fixModelBiomass",
           "forestedLCCClasses", "imputeBadAgeModel", "landis", "LCCClassesToReplaceNN",
           "LCCClassesToReplaceNNMethod", "minCoverThreshold", "minRelativeBFunction",
           "omitNonTreedPixels", "overrideAgeInFires", "overrideBiomassInFires",
           "pixelGroupAgeClass", "pixelGroupBiomassClass", "rmImputedPix",
           "speciesTableAreas", "speciesUpdateFunction", "sppEquivCol",
           "subsetDataAgeModel", "subsetDataAttempts", "subsetDataBiomassModel",
           "successionTimestep", "useCloudCacheForStats", "vegLeadingProportion"))
  )
})

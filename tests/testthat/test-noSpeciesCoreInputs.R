## Four ELFs (3.2.1, 3.2.4, 3.2.5, 3.3.2) have no tree species at all, so `sim$speciesLayers`
## arrives as a zero-layer SpatRaster (the shape fixed by the no-tree-species contract). Every
## step of createBiomass_coreInputs() after the raster checks estimates tree traits from species
## cover, and stopped on the empty input -- most memorably at "No trait values were found for .",
## where `0 == 0` made the check fire and there was nothing to name. These pin the no-species
## branch and the checks it must not weaken.

emptySpeciesLayers <- function(template) {
  ## LandR::.emptySpatRaster() is the one place that owns this (terra exports no way to build a
  ## zero-layer SpatRaster); repeated here only as a test fixture, for LandR versions without it.
  if (exists(".emptySpatRaster", asNamespace("LandR"))) {
    return(get(".emptySpatRaster", asNamespace("LandR"))(template))
  }
  out <- methods::new("SpatRaster")
  out@pntr <- terra:::SpatRaster$new(c(nrow(template), ncol(template), 0),
                                     as.vector(terra::ext(template)), terra::crs(template))
  out
}

testRasterToMatch <- function() {
  rtm <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2400, ymin = 0, ymax = 2400,
                     crs = "EPSG:3978", vals = 1)
  rtm[1:5] <- NA ## some pixels outside the study area
  rtm
}

moduleExprs <- function() {
  parse(testthat::test_path("..", "..", "Biomass_borealDataPrep.R"), keep.source = FALSE)
}

moduleFunctionBody <- function(name) {
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name(name)), moduleExprs())
  stopifnot(length(def) == 1L)
  body(eval(def[[1]][[3]]))
}

test_that("noSpeciesCoreInputs returns the contracted cohortData and pixelGroupMap", {
  rtm <- testRasterToMatch()
  out <- noSpeciesCoreInputs(rtm)

  expect_s3_class(out$cohortData, "data.table")
  expect_identical(nrow(out$cohortData), 0L)
  expect_identical(names(out$cohortData),
                   c("pixelGroup", "ecoregionGroup", "speciesCode", "age", "B", "totalBiomass"))
  expect_true(is.factor(out$cohortData$ecoregionGroup))
  expect_true(is.factor(out$cohortData$speciesCode))
  expect_true(is.integer(out$cohortData$pixelGroup))

  expect_s4_class(out$pixelGroupMap, "SpatRaster")
  expect_true(terra::compareGeom(out$pixelGroupMap, rtm))
  ## no tree pixel groups, and NA exactly where rtm is NA
  expect_identical(is.na(terra::values(out$pixelGroupMap, mat = FALSE)),
                   is.na(terra::values(rtm, mat = FALSE)))
  expect_true(all(terra::values(out$pixelGroupMap, mat = FALSE) == 0, na.rm = TRUE))
})

test_that("createBiomass_coreInputs takes the no-species branch and skips the trait estimation", {
  rtm <- testRasterToMatch()
  sim <- new.env()
  sim$speciesLayers <- emptySpeciesLayers(rtm)
  sim$rasterToMatch <- rtm

  ## the branch, extracted from the module source so the test tracks the real code
  bod <- as.list(moduleFunctionBody("createBiomass_coreInputs"))
  branch <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
                     identical(deparse(x[[2]]), "nlyr(sim$speciesLayers) == 0L"), bod)
  expect_length(branch, 1L)

  e <- new.env()
  e$sim <- sim
  e$nlyr <- terra::nlyr
  e$noSpeciesCoreInputs <- noSpeciesCoreInputs
  expect_no_error(withCallingHandlers(eval(branch[[1]], e), message = function(m) invokeRestart("muffleMessage")))

  expect_identical(nrow(sim$cohortData), 0L)
  expect_s4_class(sim$pixelGroupMap, "SpatRaster")
  expect_true(terra::compareGeom(sim$pixelGroupMap, rtm))
})

test_that("a NULL speciesLayers still stops, and layers that exist are not diverted", {
  bod <- as.list(moduleFunctionBody("createBiomass_coreInputs"))

  ## mis-ordered modules (NULL) is a different, still-fatal condition
  nullCheck <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
                        identical(deparse(x[[2]]), "is.null(sim$speciesLayers)"), bod)
  expect_length(nullCheck, 1L)
  e <- new.env()
  e$sim <- list(speciesLayers = NULL)
  expect_error(eval(nullCheck[[1]], e), "speciesLayers' are missing")

  ## with real layers the branch condition is FALSE, so the normal path is unchanged
  rtm <- testRasterToMatch()
  e2 <- new.env()
  e2$sim <- list(speciesLayers = c(rtm, rtm))
  e2$nlyr <- terra::nlyr
  expect_false(eval(quote(nlyr(sim$speciesLayers) == 0L), e2))
})

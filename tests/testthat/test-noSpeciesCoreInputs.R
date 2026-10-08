## Four ELFs (3.2.1, 3.2.4, 3.2.5, 3.3.2) have no tree species at all: `sim$sppEquiv` has zero
## rows and `sim$speciesLayers` arrives as NULL (a zero-layer SpatRaster was tried first and
## rejected -- terra cannot wrap(), unwrap() or write one, so it did not survive Cache). Every
## step of createBiomass_coreInputs() after the raster checks estimates tree traits from species
## cover, and stopped on the empty input -- most memorably at "No trait values were found for .",
## where `0 == 0` made the check fire and there was nothing to name. These pin the no-species
## branch and the checks it must not weaken.

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
  sim$speciesLayers <- NULL
  sim$sppEquiv <- data.table::data.table(LandR = character(0), FuelClass = character(0))
  sim$rasterToMatch <- rtm

  ## the branch, extracted from the module source so the test tracks the real code
  bod <- as.list(moduleFunctionBody("createBiomass_coreInputs"))
  branch <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
                     identical(deparse(x[[2]]), "noSpecies"), bod)
  expect_length(branch, 1L)

  e <- new.env()
  e$sim <- sim
  e$noSpecies <- TRUE
  e$noSpeciesCoreInputs <- noSpeciesCoreInputs
  ## prepSpeciesTable() with a 0-row sppEquiv is LandR's; here it is stubbed to what it returns
  e$prepSpeciesTable <- function(...) data.table::data.table(species = character(0))
  e$Cache <- identity
  e$P <- function(sim) list(speciesTableAreas = "BSW", sppEquivCol = "LandR")
  expect_no_error(withCallingHandlers(eval(branch[[1]], e), message = function(m) invokeRestart("muffleMessage")))

  expect_identical(nrow(sim$cohortData), 0L)
  expect_s4_class(sim$pixelGroupMap, "SpatRaster")
  expect_true(terra::compareGeom(sim$pixelGroupMap, rtm))
  ## declared outputs must be assigned, not left NULL (see the last test)
  expect_identical(nrow(sim$species), 0L)
  expect_identical(nrow(sim$sufficientLight), 5L)
  expect_identical(nrow(sim$speciesEcoregion), 0L)
})

test_that("the no-species branch builds species with the same call as the with-species path", {
  bod <- as.list(moduleFunctionBody("createBiomass_coreInputs"))
  branch <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
                     identical(deparse(x[[2]]), "noSpecies"), bod)[[1]]
  isSpeciesAssign <- function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
    identical(deparse(x[[2]]), "sim$species")
  inBranch <- Filter(isSpeciesAssign, as.list(branch[[3]]))
  atTopLevel <- Filter(isSpeciesAssign, bod)
  expect_length(inBranch, 1L)
  expect_gte(length(atTopLevel), 1L)
  expect_identical(inBranch[[1]][[3]], atTopLevel[[1]][[3]])
})

test_that("a NULL speciesLayers still stops when there ARE species, and is not an error when there are none", {
  bod <- as.list(moduleFunctionBody("createBiomass_coreInputs"))

  ## mis-ordered modules (NULL layers, species expected) is a different, still-fatal condition
  nullCheck <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
                        identical(deparse(x[[2]]), "is.null(sim$speciesLayers) && !noSpecies"), bod)
  expect_length(nullCheck, 1L)
  e <- new.env()
  e$sim <- list(speciesLayers = NULL)
  e$noSpecies <- FALSE
  expect_error(eval(nullCheck[[1]], e), "speciesLayers' are missing")
  e$noSpecies <- TRUE
  expect_no_error(eval(nullCheck[[1]], e))

  ## `noSpecies` is decided by sppEquiv, the single place "no tree species" is established
  noSpp <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                    identical(x[[2]], as.name("noSpecies")), bod)
  expect_length(noSpp, 1L)
  e3 <- new.env()
  e3$sim <- list(sppEquiv = data.table::data.table(LandR = character(0)))
  expect_true(eval(noSpp[[1]][[3]], e3))
  e3$sim <- list(sppEquiv = data.table::data.table(LandR = "Pice_mar"))
  expect_false(eval(noSpp[[1]][[3]], e3))
  e3$sim <- list(sppEquiv = NULL) ## not yet supplied is not "no species"
  expect_false(eval(noSpp[[1]][[3]], e3))
})

## sufficientLight and speciesEcoregion are empty (or constant) but MUST be present:
## suppliedElsewhere() reports TRUE for anything this module declares in createsOutput(), so
## Biomass_regeneration's .inputObjects fallbacks are suppressed and it reads
## sim$species / sufficientLight / speciesEcoregion unguarded (Biomass_regeneration.R:249-252).

test_that("the no-species branch also supplies sufficientLight and speciesEcoregion", {
  out <- noSpeciesCoreInputs(testRasterToMatch())

  expect_named(out, c("cohortData", "pixelGroupMap", "sufficientLight", "speciesEcoregion"))

  ## the columns makeSpeciesEcoregion() returns, with no rows
  expect_identical(nrow(out$speciesEcoregion), 0L)
  expect_identical(names(out$speciesEcoregion),
                   c("ecoregionGroup", "speciesCode", "establishprob", "maxB", "maxANPP", "year"))

  ## sufficientLight is not species-dependent: it is the full LANDIS-test table, not empty
  expect_identical(nrow(out$sufficientLight), 5L)
  expect_identical(out$sufficientLight$speciesshadetolerance, 1:5)

  expect_false(any(vapply(out, is.null, logical(1))))
})

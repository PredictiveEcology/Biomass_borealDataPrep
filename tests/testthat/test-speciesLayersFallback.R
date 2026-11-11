## Without Biomass_speciesData, .inputObjects makes `speciesLayers` itself from SCANFI. In an ELF
## with no tree species, prepSpeciesLayers_SCANFI() returns NULL, which must pass through
## (init's no-species branch expects it) rather than reach NAcover2zero(), which errors on NULL.

speciesLayersFallback <- function() {
  exprs <- parse(testthat::test_path("..", "..", "Biomass_borealDataPrep.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name(".inputObjects")), exprs)
  stopifnot(length(def) == 1L)
  bod <- as.list(body(eval(def[[1]][[3]])))
  branch <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
                     identical(deparse(x[[2]]), "!suppliedElsewhere(\"speciesLayers\", sim)"), bod)
  stopifnot(length(branch) == 1L)
  branch[[1]]
}

test_that("a NULL SCANFI result (no tree species) passes through the speciesLayers fallback", {
  rtm <- terra::rast(nrows = 4, ncols = 4, crs = "EPSG:3978", vals = 1)
  sim <- new.env()
  sim$rasterToMatch_biomassParam <- rtm
  sim$sppEquiv <- data.table::data.table(LandR = character(0))

  e <- new.env()
  e$sim <- sim
  e$dPath <- tempdir()
  e$cacheTags <- "test"
  e$suppliedElsewhere <- function(...) FALSE
  e$P <- function(sim) list(.sslVerify = TRUE, .studyAreaName = "test", sppEquivCol = "LandR",
                            dataYear = 2020)
  e$prepSpeciesLayers_SCANFI <- function(...) NULL
  e$Cache <- function(x, ...) x

  expect_no_error(eval(speciesLayersFallback(), envir = e))
  expect_null(sim$speciesLayers)
})

## When `studyArea` is smaller than `studyArea_biomassParam`, createBiomass_coreInputs() renumbers
## `pixelIndex` from the `rasterToMatch_biomassParam` grid to the `rasterToMatch` grid, but
## `sim$firePerimeters` stays on `rasterToMatch_biomassParam`. The young-cohort spin-up looked up
## fire status with the renumbered index, so cohorts were picked by the fire status of other cells.

## rasterToMatch_biomassParam is 6 x 10 cells; rasterToMatch is the 4 x 5 window in its rows 2-5
## and columns 3-7, so new index n (window row r, column c) is old index 10 * r + c + 2.
## Fire is in old cells 24 and 35 (new 7 and 13) and in old cell 5, outside the window.
testGrids <- function() {
  bp <- terra::rast(nrows = 6, ncols = 10, xmin = 0, xmax = 2400, ymin = 0, ymax = 1440,
                    crs = "EPSG:3978", vals = 1)
  bp[c(1, 60)] <- NA ## outside the window, so the two grids' NA counts differ
  rtm <- terra::rast(nrows = 4, ncols = 5, xmin = 480, xmax = 1680, ymin = 240, ymax = 1200,
                     crs = "EPSG:3978", vals = 1)
  fire <- terra::rast(bp)
  fire[] <- NA_real_
  fire[c(5, 24, 35)] <- c(2001, 2010, 2015)
  list(bp = bp, rtm = rtm, fire = fire)
}

moduleFunctionBody <- function(name) {
  exprs <- parse(testthat::test_path("..", "..", "Biomass_borealDataPrep.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name(name)), exprs)
  stopifnot(length(def) == 1L)
  body(eval(def[[1]][[3]]))
}

## every `name <- ...` anywhere in `expr`
findAssignments <- function(expr, name) {
  out <- list()
  if (identical(expr[[1]], as.name("<-")) && identical(expr[[2]], as.name(name))) {
    out <- list(expr)
  }
  for (i in seq_along(expr)[-1]) {
    if (is.call(expr[[i]])) out <- c(out, findAssignments(expr[[i]], name))
  }
  out
}

test_that("inFirePerimeters() looks up renumbered pixels in their own cells", {
  g <- testGrids()
  idx <- terra::crop(terra::setValues(g$bp, seq_len(terra::ncell(g$bp))), g$rtm)
  map <- data.table::data.table(pixelIndex = terra::values(idx, mat = FALSE),
                                newPixelIndex = seq_len(terra::ncell(idx)))
  expect_equal(map$pixelIndex[c(1, 5, 7, 13, 18)], c(13, 17, 24, 35, 45))

  young <- c(5L, 7L, 13L, 18L) ## on rasterToMatch
  expected <- c(FALSE, TRUE, TRUE, FALSE)
  expect_identical(inFirePerimeters(young, g$fire, map), expected)
  ## the fixture can tell the grids apart: read straight off firePerimeters, the new indices
  ## land in other cells
  expect_false(identical(!is.na(terra::values(g$fire, mat = FALSE)[young]), expected))

  ## nothing renumbered: pixelIndex is already on the firePerimeters grid
  expect_identical(inFirePerimeters(c(5L, 17L, 24L, 45L), g$fire), c(TRUE, FALSE, TRUE, FALSE))

  expect_error(inFirePerimeters(99L, g$fire, map), "not in pixelIndexMap")
})

test_that("createBiomass_coreInputs() picks young cohorts in fire cells after subsetting to studyArea", {
  g <- testGrids()
  bod <- moduleFunctionBody("createBiomass_coreInputs")

  ## the module's own subsetting block and fire lookup, so the test tracks the real code
  subsetBlock <- Filter(function(x) {
    is.call(x) && identical(x[[1]], as.name("if")) &&
      length(findAssignments(x[[3]], "newPixelIndexDT")) > 0
  }, as.list(bod))
  expect_length(subsetBlock, 1L)
  fireLookup <- findAssignments(bod, "youngRows2")
  expect_length(fireLookup, 1L)

  sim <- new.env()
  sim$rasterToMatch <- g$rtm
  sim$rasterToMatch_biomassParam <- g$bp
  sim$firePerimeters <- g$fire
  sim$imputedPixID <- integer(0)

  ## one young cohort in old cells 13, 17, 24, 35, 45 (in the window) and 5 (outside, dropped)
  oldIdx <- c(5L, 13L, 17L, 24L, 35L, 45L)
  e <- new.env()
  e$sim <- sim
  e$cacheTags <- "test"
  e$Cache <- function(x, ...) x
  ## attached by SpaDES when the module runs; the package rendition does not import LandR
  e$.compareRas <- LandR::.compareRas
  e$pixelCohortData <- data.table::data.table(pixelIndex = oldIdx, oldIdx = oldIdx,
                                              speciesCode = "Pice_mar", age = 20L, B = 100L)
  suppressMessages(eval(subsetBlock[[1]], e))

  ## renumbered to rasterToMatch
  pcd <- e$pixelCohortData[order(oldIdx)]
  expect_identical(pcd$oldIdx, c(13L, 17L, 24L, 35L, 45L))
  expect_identical(pcd$pixelIndex, c(1L, 5L, 7L, 13L, 18L))

  e$young <- e$pixelCohortData
  eval(fireLookup[[1]], e)
  expect_identical(sort(e$young$oldIdx[e$youngRows2]), c(24L, 35L))
})

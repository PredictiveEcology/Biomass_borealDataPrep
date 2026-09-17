## .inputObjects stopped unless dataYear was 2000, 2010 or 2020, although SCANFI V2 provides every
## 5 years from 1985 to 2025, so 1985 and 1990 could not be used. The check is now left to
## LandR::prepRawBiomassMap(), which knows which years each source provides; this pins that it
## still stops, before downloading anything, on a year the source does not have.
test_that("prepRawBiomassMap stops on a year SCANFI V2 does not provide", {
  skip_if_not_installed("LandR")
  expect_error(
    LandR::prepRawBiomassMap(dataSource = "SCANFI", dataYear = "1987", dataVersion = "V2",
                             destinationPath = withr::local_tempdir()),
    "SCANFI V2 data is currently available for 1985, 1990"
  )
})

## `prepSpeciesLayers_SCANFI()`'s formal is `dataYear`. Partial matching does not save a
## `year =` argument, because "year" is not a prefix of "dataYear": it falls into `...` and is
## silently discarded, so species cover was always the 2020 set while rawBiomassMap and
## standAgeMap honoured dataYear -- a run with dataYear != 2020 mixed years without saying so.
## The identical bug was fixed for Biomass_speciesData in its #47; this module was missed in
## that sweep, so it is pinned here. Static, so it costs nothing and needs no Drive access.
test_that("the SCANFI species-layer call passes dataYear, not year", {
  code <- parse(file.path(moduleRoot, paste0(moduleName, ".R")))
  found <- list()
  walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("prepSpeciesLayers_SCANFI"))) {
        found[[length(found) + 1L]] <<- x
      }
      lapply(as.list(x), walk)
    }
    invisible(NULL)
  }
  lapply(as.list(code), walk)
  expect_length(found, 1L)

  argNames <- names(as.list(found[[1L]]))[-1L]
  expect_true("dataYear" %in% argNames)
  expect_false("year" %in% argNames)
})

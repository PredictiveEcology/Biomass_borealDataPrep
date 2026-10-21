#' Which pixels lie inside a fire perimeter
#'
#' `firePerimeters` is on the `rasterToMatch_biomassParam` grid. When `studyArea` is smaller than
#' `studyArea_biomassParam`, `createBiomass_coreInputs()` renumbers `pixelIndex` to the
#' `rasterToMatch` grid. The two grids have different widths, so the same index names a different
#' cell on each; `pixelIndexMap` maps the renumbered indices back before the lookup.
#'
#' @param pixelIndex integer vector of pixel indices, e.g. `pixelCohortData$pixelIndex`.
#' @param firePerimeters `SpatRaster` of fire years on the `rasterToMatch_biomassParam` grid,
#'   `NA` outside fire perimeters.
#' @param pixelIndexMap `NULL` if `pixelIndex` is on the `firePerimeters` grid. Otherwise a
#'   `data.table` with columns `pixelIndex` (on `rasterToMatch_biomassParam`) and
#'   `newPixelIndex` (on `rasterToMatch`), as built when subsetting to `studyArea`.
#'
#' @return logical vector the length of `pixelIndex`, `TRUE` inside a fire perimeter.
#' @noRd
inFirePerimeters <- function(pixelIndex, firePerimeters, pixelIndexMap = NULL) {
  if (!is.null(pixelIndexMap)) {
    pixelIndex <- pixelIndexMap$pixelIndex[match(pixelIndex, pixelIndexMap$newPixelIndex)]
    if (anyNA(pixelIndex)) {
      stop("inFirePerimeters(): some pixelIndex values are not in pixelIndexMap$newPixelIndex")
    }
  }
  !is.na(terra::values(firePerimeters, mat = FALSE)[pixelIndex])
}

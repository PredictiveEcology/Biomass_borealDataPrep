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

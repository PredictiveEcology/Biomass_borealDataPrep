library(data.table)

test_that("young-cohort limit passes a pioneer 17% over the age-based limit but not one above maxRawB", {
  maxRawB <- 30000
  longevity <- c(Lari_lar = 400, Pice_mar = 250)
  maxAgeHighQualityData <- 40
  oldLimit <- 2.8 * maxRawB / min(longevity / maxAgeHighQualityData)
  lim <- youngBiomassLimit(maxRawB, longevity, maxAgeHighQualityData)
  expect_lt(lim, maxRawB)
  expect_true(all(1.17 * oldLimit <= lim))
  expect_false(all(1.17 * oldLimit <= oldLimit))
  ## capped at the raw map maximum
  expect_equal(youngBiomassLimit(maxRawB, longevity, maxAgeHighQualityData = 400), maxRawB)
  expect_false(all(maxRawB + 1 <= youngBiomassLimit(maxRawB, longevity, 400)))
})

test_that("spin-up starts every cohort at the new-cohort biomass, not zero", {
  cd <- data.table::data.table(pixelIndex = 1:2, speciesCode = c("Pice_mar", "Betu_pap"),
                               age = c(12L, 7L), B = c(500L, 800L))
  out <- spinUpStartCohorts(cd, initialB = 10, speciesEcoregion = NULL)
  expect_identical(out$B, c(10L, 10L))
  expect_identical(out$age, c(1L, 1L))
  expect_identical(cd$B, c(500L, 800L)) # input untouched
})

test_that("spin-up with NA initialB starts cohorts at maxANPP from speciesEcoregion", {
  cd <- data.table::data.table(pixelIndex = 1:3, speciesCode = c("Pice_mar", "Betu_pap", "Pice_mar"),
                               ecoregionGroup = c("a", "a", "b"), age = c(12L, 7L, 9L), B = c(500L, 800L, 5L))
  se <- data.table::data.table(speciesCode = c("Pice_mar", "Betu_pap", "Pice_mar"),
                               ecoregionGroup = c("a", "a", "b"), maxANPP = c(107L, 150L, 0L), maxB = 3000L)
  out <- spinUpStartCohorts(cd, initialB = NA, speciesEcoregion = se)
  expect_identical(out$B, c(107L, 150L, 0L))
  expect_identical(out$age, c(1L, 1L, 1L))
})

## Biomass_core stops on a cohort whose speciesCode x ecoregionGroup has no speciesEcoregion row
## (it used to drop it silently). Such cohorts (zero-biomass or imputed pixels) reach the spin-up,
## so they are left out of it and come back with B = NA, to be filled from the biomass model.
test_that("spin-up leaves out cohorts with no speciesEcoregion row and does not change the others", {
  skip_if_not_installed("LandR")
  skip_if_not_installed("SpaDES.core")
  ## the module's code runs with its reqdPkgs attached; the package rendition does not import them all
  for (pkg in c("SpaDES.core", "LandR", "terra")) withr::local_package(pkg)
  bcPath <- file.path(modulePath, "Biomass_core")
  skip_if_not(dir.exists(bcPath), "Biomass_core not next to the module")
  skip_if_not(isTRUE(SpaDES.core::moduleVersion("Biomass_core", modulePath) >= "2.1.0.9002"),
              "Biomass_core without the speciesEcoregion check")
  mp <- withr::local_tempdir()
  file.symlink(moduleRoot, file.path(mp, moduleName))
  file.symlink(bcPath, file.path(mp, "Biomass_core"))
  paths <- list(cachePath = tempfile(), inputPath = tempfile(), modulePath = mp, outputPath = tempfile())
  for (p in paths) dir.create(p)

  sppEquiv <- data.table::as.data.table(LandR::sppEquivalencies_CA)[LandR %in% c("Pice_mar", "Abie_bal")]
  sppEquiv[, EN_generic_short := LandR]
  sppCodes <- factor(c("Pice_mar", "Abie_bal"))
  species <- data.table::data.table(
    species = as.character(sppCodes), Area = "BSW", longevity = c(250L, 200L), sexualmature = c(30L, 25L),
    shadetolerance = c(4, 5), firetolerance = 1L, seeddistance_eff = c(80L, 25L),
    seeddistance_max = c(200L, 160L), resproutprob = 0, resproutage_min = 0L, resproutage_max = 0L,
    postfireregen = "none", leaflongevity = 3L, wooddecayrate = 0.1, mortalityshape = 15L,
    growthcurve = 0, leafLignin = 0.2, hardsoft = "soft", speciesCode = sppCodes)
  ecoregions <- factor(c("a", "b"))
  cd <- data.table::data.table(
    pixelIndex = 1:4, speciesCode = sppCodes[c(1, 1, 2, 2)], ecoregionGroup = ecoregions[c(1, 1, 1, 2)],
    age = c(5L, 8L, 6L, 7L), B = c(10L, 20L, 30L, 0L))
  ## ecoregion "b" has no speciesEcoregion row
  se <- data.table::data.table(
    speciesCode = sppCodes, ecoregionGroup = ecoregions[1], establishprob = 0.5,
    maxB = c(5000L, 4000L), maxANPP = c(150L, 120L), year = 0)
  spinUp <- function(x) {
    spinUpPartial(
      x, se, maxAge = 10,
      minRelativeB = data.frame(ecoregionGroup = ecoregions, X1 = 0.2, X2 = 0.4, X3 = 0.5, X4 = 0.7, X5 = 0.9),
      species = species, sppEquiv = sppEquiv, sppEquivCol = "LandR",
      sppColorVect = c(Pice_mar = "red", Abie_bal = "blue"), paths = paths,
      currentModule = moduleName, modules = "Biomass_core", initialB = 10)
  }
  msgs <- character()
  capture.output(withCallingHandlers(
    out <- spinUp(cd),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
  expect_true(any(grepl("Leaving out of the spin-up 1 cohort\\(s\\) in 1 .* e\\.g\\. Abie_bal x b", msgs)))
  expect_true(is.na(out$B[out$pixelIndex == 4]))
  expect_false(anyNA(out$B[out$pixelIndex != 4]))
  expect_identical(out$pixelIndex, cd$pixelIndex)
  ## the cohorts that had traits get what they get without the left-out cohort
  capture.output(ref <- spinUp(cd[pixelIndex != 4]))
  expect_equal(out[pixelIndex != 4], ref)
})

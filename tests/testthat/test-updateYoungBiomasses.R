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
  out <- spinUpStartCohorts(cd, initialB = 10)
  expect_identical(out$B, c(10L, 10L))
  expect_identical(out$age, c(1L, 1L))
  expect_identical(cd$B, c(500L, 800L)) # input untouched
})

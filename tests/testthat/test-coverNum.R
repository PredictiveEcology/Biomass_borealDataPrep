## coverNum is the denominator of the establishment (cover presence) model: the number of
## PIXELS in each ecoregionGroup. It used to count cohort rows, so a pixel with three species
## counted three times and every presence probability was deflated.
library(data.table)

test_that("coverNum counts pixels, not cohorts", {
  cd <- data.table(pixelIndex = c(1, 1, 1, 2, 2, 3),
                   speciesCode = c("a", "b", "c", "a", "b", "a"),
                   ecoregionGroup = factor("1_210"))
  pt <- data.table(pixelIndex = 1:3)
  out <- coverNumByGroup(cd, pt)
  expect_identical(out$ecoregionGroup, "1_210")
  expect_identical(out$coverNum, 3L)        # the old count was 6
})

test_that("coverNum is per group, only for pixels in pixelTable, and skips NA groups", {
  cd <- data.table(pixelIndex = c(1, 1, 2, 3, 3, 4, 5),
                   speciesCode = c("a", "b", "a", "a", "b", "a", "a"),
                   ecoregionGroup = c("1_210", "1_210", "1_210", "2_220", "2_220", "2_220", NA))
  pt <- data.table(pixelIndex = c(1, 2, 3, 5))   # pixel 4 is not in pixelTable
  out <- coverNumByGroup(cd, pt)
  expect_identical(out$ecoregionGroup, c("1_210", "2_220"))
  expect_identical(out$coverNum, c(2L, 1L))
})

test_that("a species in every pixel has presence probability 1", {
  cd <- data.table(pixelIndex = rep(1:4, each = 2), speciesCode = rep(c("a", "b"), 4),
                   cover = c(10, 5, 20, 0, 30, 0, 40, 5), ecoregionGroup = "1_210")
  pres <- cd[, list(coverPres = sum(cover > 0)), by = c("ecoregionGroup", "speciesCode")]
  num <- coverNumByGroup(cd, data.table(pixelIndex = 1:4))
  p <- pres[num, on = "ecoregionGroup"][, coverPres / coverNum]
  expect_identical(p, c(1, 0.5))
})

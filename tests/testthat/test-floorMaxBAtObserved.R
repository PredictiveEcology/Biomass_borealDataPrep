## maxB is floored at the largest biomass the species was observed at in its group: a fit clamped
## to 0 would otherwise stop a species growing where it demonstrably grows.
library(data.table)

se <- data.table(ecoregionGroup = factor(c("1_081", "1_081", "1_210")),
                 speciesCode = factor(c("Betu_pap", "Pice_mar", "Pice_mar")),
                 maxB = c(0L, 3200L, 5400L), maxANPP = c(0L, 107L, 180L), establishprob = 0.5)
cd <- data.table(ecoregionGroup = c("1_081", "1_081", "1_081", "1_210"),
                 speciesCode = c("Betu_pap", "Betu_pap", "Pice_mar", "Pice_mar"),
                 B = c(300L, 700L, 2900L, 6000L))

test_that("a maxB below the observed maximum is raised to it, with maxANPP to match", {
  out <- floorMaxBAtObserved(se, cd)
  expect_identical(out$maxB, c(700L, 3200L, 6000L))
  expect_identical(out$maxANPP, c(23L, 107L, 200L))   # asInteger(700/30), untouched, 6000/30
  expect_identical(attr(out, "nRaised"), 2L)
})

test_that("rows already above the observations, and the other columns, are untouched", {
  out <- floorMaxBAtObserved(se, cd)
  expect_identical(out[speciesCode == "Pice_mar" & ecoregionGroup == "1_081", maxB], 3200L)
  expect_identical(out$establishprob, se$establishprob)
  expect_identical(names(out), names(se))
  expect_identical(levels(out$ecoregionGroup), levels(se$ecoregionGroup))
})

test_that("a group-species with no observations keeps its fitted value", {
  out <- floorMaxBAtObserved(se, cd[speciesCode != "Betu_pap"])
  expect_identical(out[speciesCode == "Betu_pap", maxB], 0L)
})

test_that("the input table is not modified", {
  before <- copy(se)
  invisible(floorMaxBAtObserved(se, cd))
  expect_identical(se, before)
})

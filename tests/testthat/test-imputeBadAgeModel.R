## Biomass_borealDataPrep#131: `imputeBadAgeModel`'s old default fit age directly, which
## could predict a negative age for a young, high-cover, low-biomass stand; that got clamped
## to age 0 while biomass/cover stayed positive, which CBMutils::cumPoolsCreateAGB() rejects.
## The default now comes from LandR, whose response is log(age), so this checks the module
## did not drift back to a locally-written formula.

test_that("imputeBadAgeModel default is LandR::imputeBadAgeModelDefault()", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  default <- md$parameters[md$parameters$paramName == "imputeBadAgeModel", "default"][[1]]
  expect_identical(default, LandR::imputeBadAgeModelDefault())
})

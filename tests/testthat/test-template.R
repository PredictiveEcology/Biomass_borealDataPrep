
# Please do three things to ensure this template is correctly modified:
# 1. Rename this file based on the content you are testing using
#    `test-functionName.R` format so that your can directly call `moduleCoverage`
#    to calculate module coverage information.
#    `functionName` is a function's name in your module (e.g., `Biomass_coreDataPrepEvent1`).
# 2. Copy this file to the tests folder (i.e., `~/GitHub/LandWeb/Biomass_borealDataPrep/tests/testthat`).

# 3. Modify the test description based on the content you are testing:
test_that("test Event1 and Event2", {
  origLibPaths <- .libPaths()
  on.exit(.libPaths(origLibPaths))
  origDir <- getwd()
  td <- reproducible::tempdir2()
  moduleFile <- dir(file.path("..", ".."), pattern = "\\.R$", full.names = TRUE)
  moduleName <- basename(gsub("\\.R", "", moduleFile))
  modulePathIn <- file.path(td, moduleName)
  reproducible::checkPath(modulePathIn, create = TRUE)
  modulePath <- dirname(modulePathIn)
  file.copy(moduleFile, file.path(modulePathIn, basename(moduleFile)))
  withr::local_options("SpaDES.project.updateRprofile" = FALSE)

  warns <- capture_warning(
    out <- SpaDES.project::setupProject(paths = list(projectPath = td,
                                                     modulePath = modulePath,
                                                     packagePath = .libPaths()[1]),
                                        modules = moduleName,
                                        times = list(start = 0, end = 1))
  )
  browser()
  sa <- setupStudyArea(list(NAME_1 = "Alberta", "NAME_2" = "Division No. 17", level = 2))
  terra
  expect_error(outFinal <- SpaDES.core::simInitAndSpades2(out), regexp = "Please provide a.*polygon")
  out$studyArea <- sa
  expect_error(outFinal <- SpaDES.core::simInitAndSpades2(out), regexp = "Please provide a.*Large.*polygon")
  out$studyAreaLarge <- terra::buffer(sa, 1000)
  outFinal <- SpaDES.core::simInitAndSpades2(out)

  skip("rest is not yet tested")

  # module <- list(moduleName)
  # path <- list(modulePath = modulePath,
  #              outputPath = file.path(td, "outputs"))
  # parameters <- list(
  #   #.progress = list(type = "graphical", interval = 1),
  #   .globals = list(verbose = FALSE),
  #   Biomass_borealDataPrep = list(.saveInitialTime = NA)
  # )
  # times <- list(start = 0, end = 1)

  # If your test function contains `time(sim)`, you can test the function at a
  # particular simulation time by defining the start time above.
  # object1 <- "object1" # please specify
  # object2 <- "object2" # please specify
  # objects <- list("object1" = object1, "object2" = object2)
  #
  # mySim <- simInit(times = times,
  #                  params = parameters,
  #                  modules = module,
  #                  objects = objects,
  #                  paths = path)

  # if (exists("Biomass_coreDataPrepEvent1", envir = .GlobalEnv)) {
  #   simOutput <- Biomass_coreDataPrepEvent1(mySim)
  # } else {
  #   simOutput <- mySim$Biomass_coreDataPrepEvent1(mySim)
  # }
  #
  # expectedOutputEvent1Test1 <- " this is test for event 1. " # please define your expection of your output
  # expect_is(class(simOutput$event1Test1), "character")
  # expect_equal(simOutput$event1Test1, expectedOutputEvent1Test1) # or other expect function in testthat package.
  # expect_equal(simOutput$event1Test2, as.numeric(999)) # or other expect function in testthat package.
  #
  # if (exists("Biomass_coreDataPrepEvent2", envir = .GlobalEnv)) {
  #   simOutput <- Biomass_coreDataPrepEvent2(mySim)
  # } else {
  #   simOutput <- mySim$Biomass_coreDataPrepEvent2(mySim)
  # }
  #
  # expectedOutputEvent2Test1 <- " this is test for event 2. " # please define your expection of your output
  # expect_is(class(simOutput$event2Test1), "character")
  # expect_equal(simOutput$event2Test1, expectedOutputEvent2Test1) # or other expect function in testthat package.
  # expect_equal(simOutput$event2Test2, as.numeric(777)) # or other expect function in testthat package.
})

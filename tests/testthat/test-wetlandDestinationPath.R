## The CWIM wetland layer is a 5.4 GB national download. Since LandR 1.2.0.9060 prepInputs_CWIM()
## takes `destinationPath` and, when it is NULL, falls back to the working directory, so the
## .inputObjects call must pass the module's input path (`dPath`) like its neighbouring prepInputs
## calls do. Otherwise a run with no shared copy downloads the file into the project folder.

wetlandBranch <- function() {
  exprs <- parse(testthat::test_path("..", "..", "Biomass_borealDataPrep.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name(".inputObjects")), exprs)
  stopifnot(length(def) == 1L)
  bod <- as.list(body(eval(def[[1]][[3]])))
  branch <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
                     grepl("suppliedElsewhere(\"rstWetland\", sim)", deparse(x[[2]])[1], fixed = TRUE),
                   bod)
  stopifnot(length(branch) == 1L)
  branch[[1]]
}

test_that("the CWIM wetland layer is downloaded into the module's input path", {
  seen <- new.env()
  testthat::local_mocked_bindings(
    prepInputs_CWIM = function(to, destinationPath = NULL, ...) {
      seen$destinationPath <- destinationPath
      "rstWetland"
    },
    .package = "LandR"
  )

  e <- new.env()
  e$sim <- new.env()
  e$sim$rasterToMatch_biomassParam <- "rtm"
  e$dPath <- file.path(testPaths$inputPath, "wetland")
  e$suppliedElsewhere <- function(...) FALSE
  e$P <- function(sim) list(wetlandSource = "CWIM", .studyAreaName = "test")
  e$currentModule <- function(sim) "Biomass_borealDataPrep"
  e$Cache <- function(x, ...) x

  eval(wetlandBranch(), envir = e)
  expect_identical(e$sim$rstWetland, "rstWetland")
  expect_identical(seen$destinationPath, e$dPath)
})

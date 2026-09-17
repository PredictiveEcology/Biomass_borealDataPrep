## Site x composition strata. Codes stay three digits wide so the "<eco>_<code>" string keeps
## its padding and the `_240` assertions keep matching.
library(data.table)

test_that("landcover mode is the NTEMS one-axis map: wet treed 81, wet other 80", {
  skip_if_not_installed("LandR")
  skip_if_not("wetlandToLCC" %in% getNamespaceExports("LandR"))
  comp <- c(210, 220, 230, 240, 210, 50, 20, NA)
  wet  <- c(1,   1,   1,   1,   0,   1,  1,  1)
  expect_identical(siteCompositionCodes(comp, wet, "landcover"),
                   c(81, 81, 81, 81, 210, 80, 20, NA))
})

test_that("siteComposition mode keeps composition under wet ground", {
  comp <- c(210, 220, 230, 240, 210, 220, 50, 20, NA)
  wet  <- c(1,   1,   1,   1,   0,   NA,  1,  1,  1)
  expect_identical(siteCompositionCodes(comp, wet, "siteComposition"),
                   c(810, 820, 830, 840, 210, 220, 80, 20, NA))
})

test_that("every stratum code is three digits wide", {
  codes <- c(siteCompositionCodes(c(210, 220, 230, 240), c(1, 1, 1, 1), "siteComposition"),
             210, 220, 230, 240, 290, 890, 990)
  expect_true(all(nchar(as.character(codes)) == 3))
})

test_that("the classes left to convertUnwantedLCC include the wet unresolved code", {
  expect_identical(unresolvedClasses(240, "landcover"), 240)
  expect_identical(unresolvedClasses(240, "siteComposition"), c(240, 840))
})

test_that("thin strata lose composition first, then site", {
  counts <- data.table(
    eco =     c("1",  "1",  "1",  "1",  "1",  "2",  "2"),
    stratum = c(210L, 220L, 230L, 810L, 820L, 210L, 810L),
    N =       c(500L, 40L,  30L,  300L, 10L,  20L,  15L)
  )
  ## minN = 60: 220 (40) and 230 (30) are thin; pooled, 290 holds 70 >= 60 and stops there
  m <- collapseThinStrata(counts, minN = 60)
  get <- function(e, s) m[eco == e & stratum == s, newStratum]
  expect_identical(get("1", 210L), 210L)    # plenty
  expect_identical(get("1", 220L), 290L)
  expect_identical(get("1", 230L), 290L)
  expect_identical(get("1", 810L), 810L)
  expect_identical(get("1", 820L), 990L)    # alone in the wet pool, which holds 10 < 60
  ## minN = 100: the upland pool (70) is itself thin, so it loses its site as well
  m2 <- collapseThinStrata(counts, minN = 100)
  expect_identical(m2[eco == "1" & stratum == 220L, newStratum], 990L)
  expect_identical(m2[eco == "1" & stratum == 230L, newStratum], 990L)
  expect_identical(m2[eco == "1" & stratum == 210L, newStratum], 210L)
})

test_that("a thin wet pool also loses its site", {
  counts <- data.table(eco = "1", stratum = c(810L, 820L), N = c(300L, 10L))
  m <- collapseThinStrata(counts, minN = 60)
  expect_identical(m[stratum == 820L, newStratum], 990L)   # 890 holds only 10
  expect_identical(m[stratum == 810L, newStratum], 810L)
})

test_that("a pooled class that is still thin loses its site", {
  counts <- data.table(eco = c("2", "2"), stratum = c(210L, 810L), N = c(20L, 15L))
  m <- collapseThinStrata(counts, minN = 100)
  expect_identical(m$newStratum, c(990L, 990L))
})

test_that("minN = 0 pools nothing", {
  counts <- data.table(eco = "1", stratum = c(210L, 820L), N = c(1L, 1L))
  expect_identical(collapseThinStrata(counts, minN = 0)$newStratum, c(210L, 820L))
})

test_that("inferred pixels take their own stratum, then the pooled ones, then go back unresolved", {
  mapping <- data.table(eco = c("1", "1", "2"), stratum = c(210L, 220L, 810L),
                        N = c(500L, 5L, 500L), newStratum = c(210L, 290L, 810L))
  got <- stratumForInferred(
    eco     = c("1",  "1",  "1",  "2",  "2",  "3"),
    stratum = c(210L, 220L, 230L, 820L, 230L, 210L),
    mapping = mapping
  )
  expect_identical(got, c(
    210L,  # its own stratum exists
    290L,  # its own stratum was pooled
    290L,  # no 230 in eco 1, but eco 1 has an upland pool
    840L,  # wet, no 820 and no 890/990 in eco 2 -> back to the wet unresolved code
    240L,  # upland, eco 2 has neither 230 nor 290 nor 990
    240L   # an ecoregion with no estimation data at all
  ))
})

test_that("without pooling an inferred pixel either keeps its stratum or goes back to 240", {
  mapping <- data.table(eco = "1", stratum = 81L, N = 10L, newStratum = 81L)
  expect_identical(
    stratumForInferred(c("1", "1", "2"), c(81L, 210L, 81L), mapping, allowPooling = FALSE),
    c(81L, 240L, 240L)
  )
})

test_that("mismatched lengths are refused", {
  expect_error(siteCompositionCodes(1:2, 1), "same length")
  expect_error(stratumForInferred("1", 1:2, data.table()), "same length")
})

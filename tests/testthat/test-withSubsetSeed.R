## Parameter estimation subsamples are seeded, so identical runs give identical maxB -- without
## touching the random number stream anything else uses.

test_that("a seeded expression is reproducible", {
  a <- withSubsetSeed(7L, sample(1000, 5))
  b <- withSubsetSeed(7L, sample(1000, 5))
  expect_identical(a, b)
  expect_false(identical(a, withSubsetSeed(8L, sample(1000, 5))))
})

test_that("the RNG state is restored, so later draws are unaffected", {
  set.seed(123); ref <- { runif(1); runif(3) }
  set.seed(123); runif(1)
  invisible(withSubsetSeed(99L, runif(50)))
  expect_identical(runif(3), ref)   # same as if the seeded call had never run
})

test_that("NA or NULL leaves the expression unseeded and the stream moving", {
  set.seed(1); x1 <- withSubsetSeed(NA, runif(1)); x2 <- runif(1)
  set.seed(1); expect_identical(c(runif(1), runif(1)), c(x1, x2))
  expect_identical(withSubsetSeed(NULL, 42), 42)
})

test_that("with no RNG state yet, none is left behind", {
  had <- exists(".Random.seed", envir = globalenv())
  if (had) {
    old <- get(".Random.seed", envir = globalenv())
    rm(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()))
  }
  invisible(withSubsetSeed(3L, runif(1)))
  expect_false(exists(".Random.seed", envir = globalenv()))
})

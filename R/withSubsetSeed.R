## Evaluate `expr` with the random number generator seeded, then put the generator back exactly
## as it was, so seeding a subsample changes no random draw that follows it. `seed` NULL or NA
## evaluates `expr` unseeded.
withSubsetSeed <- function(seed, expr) {
  if (is.null(seed) || length(seed) != 1L || is.na(seed)) {
    return(expr)
  }
  hadSeed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (hadSeed) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  on.exit({
    if (hadSeed) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  expr
}

## `dPath` is `asPath(inputPath(sim), 1)`, and inputPath differs for every run/rep, so a Cache()
## key that holds `destinationPath` (a download folder) changes with every rep and recomputes
## identical results. Every Cache() call that wraps a call passing `destinationPath` (or
## `outputPath = dPath`) must omit it. Static, so it needs no Drive access.
test_that("Cache() calls omit the download folders of the call they wrap", {
  code <- parse(file.path(moduleRoot, paste0(moduleName, ".R")))
  found <- list()
  walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("Cache")) && length(x) > 1 && is.call(x[[2]]))
        found[[length(found) + 1L]] <<- x
      lapply(as.list(x), walk)
    }
    invisible(NULL)
  }
  lapply(as.list(code), walk)

  folderArgs <- c("destinationPath", "outputPath")
  checked <- 0L
  for (cl in found) {
    wrapped <- cl[[2]]
    passed <- intersect(folderArgs, names(as.list(wrapped)))
    ## only folders given as `dPath` (inputPath-derived); a user's `outputPath(sim)` is not touched here
    passed <- passed[vapply(passed, function(a) identical(wrapped[[a]], as.name("dPath")), logical(1))]
    if (!length(passed)) next
    checked <- checked + 1L
    omitted <- tryCatch(eval(as.list(cl)$omitArgs), error = function(e) character())
    expect_true(all(passed %in% omitted),
                info = paste("Cache() of", deparse(wrapped[[1]]), "must omit", paste(passed, collapse = ", ")))
  }
  expect_gte(checked, 5L)
})

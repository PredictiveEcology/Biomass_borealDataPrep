#!/usr/bin/env Rscript
## Compare the SCANFI arms of one area (see run-all-scanfi.sh).
##
## Usage: Rscript compare-scanfi.R <area> [arm ...]
## Writes <ROOT>/runs/compare-<area>.md and prints it.

suppressPackageStartupMessages({ library(terra); library(data.table) })

.args <- commandArgs(trailingOnly = TRUE)
if (!length(.args)) stop("usage: compare-scanfi.R <area> [arm ...]")
area <- .args[[1]]
arms <- if (length(.args) > 1) .args[-1] else c("baseline", "forestland", "landcover", "sitecomp")
ROOT <- Sys.getenv("CLASS240_ROOT", "/mnt/fast/class240-regression")

lab <- c("0" = "unclassified", "20" = "water", "30" = "rock/barren", "31" = "snow_ice",
         "32" = "rock_rubble", "33" = "exposed_barren", "40" = "bryoids", "50" = "shrubs",
         "80" = "wetland", "81" = "wetland_treed", "100" = "herbs", "210" = "coniferous",
         "220" = "broadleaf", "230" = "mixedwood", "240" = "disturbed",
         "290" = "upland pooled", "810" = "wet conifer", "820" = "wet broadleaf",
         "830" = "wet mixed", "840" = "wet unresolved", "890" = "wet pooled", "990" = "any pooled")

runDir <- function(arm) file.path(ROOT, "runs", paste0(area, "__", arm))
have <- arms[file.exists(file.path(vapply(arms, runDir, ""), "meta.csv"))]
if (!length(have)) stop("no finished runs for ", area)
missing <- setdiff(arms, have)

rd <- function(arm, nm) {
  d <- runDir(arm)
  tif <- file.path(d, paste0(nm, ".tif"))
  if (file.exists(tif)) return(rast(tif))
  f <- file.path(d, paste0(nm, ".rds"))
  if (file.exists(f)) readRDS(f) else NULL
}
freqTab <- function(r) {
  if (is.null(r)) return(data.table(value = numeric(), count = numeric()))
  f <- as.data.table(freq(r))[, .(value, count)]
  f
}
wide <- function(dt, valueCol) {
  dcast(dt, code + class ~ arm, value.var = valueCol, fill = 0)
}
md <- character()
out <- function(...) md <<- c(md, paste0(...))
tbl <- function(dt) {
  dt <- as.data.frame(dt)
  out("| ", paste(names(dt), collapse = " | "), " |")
  out("|", paste(rep("---", ncol(dt)), collapse = "|"), "|")
  for (i in seq_len(nrow(dt))) {
    v <- vapply(dt[i, ], function(x) {
      if (is.numeric(x)) format(round(x, 3), big.mark = ",", scientific = FALSE, trim = TRUE)
      else as.character(x)
    }, "")
    out("| ", paste(v, collapse = " | "), " |")
  }
  out("")
}

out("# SCANFI 2020 arms -- ", area, "")
out("")
metas <- rbindlist(lapply(have, function(a) fread(file.path(runDir(a), "meta.csv"))), fill = TRUE)
tbl(metas[, .(label, landrSha, moduleSha, stratum, wetland, elapsedMin)])
if (length(missing)) out("Not finished: ", paste(missing, collapse = ", "), "\n")

## 1. land cover in and out -----------------------------------------------------------
lcIn <- rbindlist(lapply(have, function(a) freqTab(rd(a, "rstLCC-input"))[, arm := a]))
lcOut <- rbindlist(lapply(have, function(a) freqTab(rd(a, "rstLCC"))[, arm := a]))
for (x in list(list("Land cover built by the arm's LandR (input, pixels)", lcIn),
               list("Land cover the module returns (sim$rstLCC, pixels)", lcOut))) {
  d <- x[[2]][, .(code = value, class = lab[as.character(value)], count, arm)]
  out("## ", x[[1]], "\n")
  w <- wide(d, "count")
  setcolorder(w, c("code", "class", intersect(arms, names(w))))
  tbl(w)
}

## 2. headline numbers ------------------------------------------------------------------
head1 <- rbindlist(lapply(have, function(a) {
  pgm <- rd(a, "pixelGroupMap"); cd <- rd(a, "cohortData"); se <- rd(a, "speciesEcoregion")
  imp <- rd(a, "imputedPixID"); wet <- rd(a, "rstWetland")
  data.table(
    arm = a,
    simulatedPixels = if (is.null(pgm)) NA_real_ else sum(!is.na(values(pgm, mat = FALSE))),
    wetPixels = if (is.null(wet)) NA_real_ else sum(values(wet, mat = FALSE) %in% 1),
    imputedPixels = length(imp),
    cohortRows = NROW(cd),
    pixelGroups = if (is.null(cd)) NA_real_ else uniqueN(cd$pixelGroup),
    ecoregionGroups = if (is.null(se)) NA_real_ else uniqueN(se$ecoregionGroup),
    speciesEcoregionRows = NROW(se)
  )
}))
out("## Headline\n")
tbl(head1)

## 3. strata: pixels per stratum code --------------------------------------------------
strata <- rbindlist(lapply(have, function(a) {
  em <- rd(a, "ecoregionMap")
  lvF <- file.path(runDir(a), "ecoregionMap-levels.rds")
  if (is.null(em) || !file.exists(lvF)) return(NULL)
  lv <- readRDS(lvF)
  levels(em) <- NULL   ## count raw IDs; freq() on a categorical raster returns the labels
  f <- as.data.table(freq(em))[, .(value = as.integer(value), count)]
  idCol <- names(lv)[1]
  lv[[idCol]] <- as.integer(lv[[idCol]])
  f <- lv[f, on = setNames("value", idCol)]
  f[, grp := if ("ecoregionGroup" %in% names(f)) as.character(ecoregionGroup) else as.character(f[[2]])]
  f[, code := as.integer(sub(".*_", "", grp))]
  f[, .(pixels = sum(count), groups = uniqueN(grp)), by = code][, arm := a]
}))
if (NROW(strata)) {
  strata[, class := lab[as.character(code)]]
  out("## Pixels per stratum code (from ecoregionMap)\n")
  w <- dcast(strata, code + class ~ arm, value.var = "pixels", fill = 0)
  setcolorder(w, c("code", "class", intersect(arms, names(w))))
  tbl(w)
  out("Groups (ecoregion x code) per arm:\n")
  tbl(strata[, .(groups = sum(groups)), by = arm])
}

## 4. parameters -------------------------------------------------------------------------
par <- rbindlist(lapply(have, function(a) {
  se <- rd(a, "speciesEcoregion")
  if (is.null(se)) return(NULL)
  se <- as.data.table(se)
  se[, code := as.integer(sub(".*_", "", as.character(ecoregionGroup)))]
  ## 810-890 and 81 are wet; 990 pools across sites, so it is neither
  se[, site := fifelse(code %in% c(81, 810, 820, 830, 840, 890), "wet",
                       fifelse(code == 990, "pooled (any site)", "upland"))]
  se[, arm := a]
  se[, .(arm, speciesCode, site, maxB, maxANPP, establishprob)]
}))
if (NROW(par)) {
  out("## maxB by species (mean over ecoregion groups)\n")
  w <- dcast(par[, .(maxB = mean(maxB)), by = .(arm, speciesCode)], speciesCode ~ arm,
             value.var = "maxB")
  setcolorder(w, c("speciesCode", intersect(arms, names(w))))
  tbl(w)
  out("## establishprob by species (mean over ecoregion groups)\n")
  w <- dcast(par[, .(ep = mean(establishprob)), by = .(arm, speciesCode)], speciesCode ~ arm,
             value.var = "ep")
  setcolorder(w, c("speciesCode", intersect(arms, names(w))))
  tbl(w)
  out("## maxB by species and site (arms that distinguish site)\n")
  s2 <- par[, .(maxB = mean(maxB), groups = .N), by = .(arm, speciesCode, site)]
  s2 <- s2[arm %in% s2[site != "upland", unique(arm)]]
  if (nrow(s2)) {
    tbl(dcast(s2, speciesCode + site ~ arm, value.var = "maxB")[order(speciesCode, site)])
  } else {
    out("(no arm has wet strata)\n")
  }
}

f <- file.path(ROOT, "runs", paste0("compare-", area, ".md"))
writeLines(md, f)
cat(md, sep = "\n")
message("\nwritten: ", f)

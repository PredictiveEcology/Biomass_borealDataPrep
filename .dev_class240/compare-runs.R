#!/usr/bin/env Rscript
## Stage 3 of the class-240 regression harness: diff two runs produced by run-one.R.
##
## The question a diff has to answer is not "did anything change" -- something always does --
## but "what changed, where, and is it the thing we meant to change". So the report is by
## object, and for the parameter tables it is per ecoregion group and species.
##
## Usage:
##   Rscript .dev_class240/compare-runs.R <areaA__labelA> <areaB__labelB>
## e.g.
##   Rscript .dev_class240/compare-runs.R big__baseline big__forestland

suppressPackageStartupMessages({
  library(terra); library(data.table); library(digest)
})

.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) < 2) stop("usage: compare-runs.R <runA> <runB>")
ROOT <- Sys.getenv("CLASS240_ROOT", "/mnt/fast/class240-regression")
dirA <- file.path(ROOT, "runs", .args[[1]])
dirB <- file.path(ROOT, "runs", .args[[2]])
for (d in c(dirA, dirB)) if (!dir.exists(d)) stop("no such run: ", d)

readObj <- function(dir, nm) {
  f <- file.path(dir, paste0(nm, ".rds"))
  if (!file.exists(f)) return(NULL)
  obj <- readRDS(f)
  if (inherits(obj, "PackedSpatRaster")) terra::unwrap(obj) else obj
}

hr <- function(title) message("\n", strrep("-", 72), "\n", title, "\n", strrep("-", 72))

metaA <- fread(file.path(dirA, "meta.csv")); metaB <- fread(file.path(dirB, "meta.csv"))
hr("runs")
print(rbind(metaA, metaB)[, .(label, area, landrVersion, moduleSha, elapsedMin)])

## 1. Which saved objects are identical, byte for byte -------------------------------------
hr("objects: identical or not (sha256 of the saved object)")
mA <- fread(file.path(dirA, "manifest.csv")); mB <- fread(file.path(dirB, "manifest.csv"))
cmp <- merge(mA[, .(object, shaA = sha256, bytesA = bytes)],
             mB[, .(object, shaB = sha256, bytesB = bytes)], by = "object", all = TRUE)
cmp[, identical := shaA == shaB]
print(cmp[, .(object, identical, bytesA, bytesB)])

## 2. The land-cover map: how many pixels of each class, 240 above all ----------------------
hr("rstLCC class counts")
lccA <- readObj(dirA, "rstLCC"); lccB <- readObj(dirB, "rstLCC")
if (!is.null(lccA) && !is.null(lccB)) {
  fA <- as.data.table(terra::freq(lccA))[, .(class = value, A = count)]
  fB <- as.data.table(terra::freq(lccB))[, .(class = value, B = count)]
  f <- merge(fA, fB, by = "class", all = TRUE)
  f[is.na(A), A := 0L][is.na(B), B := 0L][, delta := B - A]
  print(f)
}

## 3. The parameters themselves ------------------------------------------------------------
hr("speciesEcoregion: maxB / maxANPP / establishprob")
seA <- readObj(dirA, "speciesEcoregion"); seB <- readObj(dirB, "speciesEcoregion")
if (!is.null(seA) && !is.null(seB)) {
  keyCols <- intersect(c("ecoregionGroup", "speciesCode", "year"), names(seA))
  valCols <- intersect(c("maxB", "maxANPP", "establishprob"), names(seA))
  j <- merge(as.data.table(seA)[, c(keyCols, valCols), with = FALSE],
             as.data.table(seB)[, c(keyCols, valCols), with = FALSE],
             by = keyCols, all = TRUE, suffixes = c(".A", ".B"))
  message("groups only in A: ", j[is.na(get(paste0(valCols[1], ".B"))), .N],
          " | only in B: ", j[is.na(get(paste0(valCols[1], ".A"))), .N],
          " | in both: ", j[!is.na(get(paste0(valCols[1], ".A"))) &
                              !is.na(get(paste0(valCols[1], ".B"))), .N])
  for (v in valCols) {
    a <- j[[paste0(v, ".A")]]; b <- j[[paste0(v, ".B")]]
    ok <- !is.na(a) & !is.na(b)
    if (!any(ok)) next
    d <- b[ok] - a[ok]
    message(sprintf("  %-14s changed in %d of %d rows | mean delta %+.2f | max |delta| %.2f | median A %.1f B %.1f",
                    v, sum(abs(d) > 1e-8), sum(ok), mean(d), max(abs(d)),
                    median(a[ok]), median(b[ok])))
  }
  ## the groups that moved most, which is where to look first
  v1 <- valCols[[1]]
  j[, delta := get(paste0(v1, ".B")) - get(paste0(v1, ".A"))]
  worst <- j[!is.na(delta)][order(-abs(delta))][1:min(10, .N)]
  message("\n  largest ", v1, " changes:")
  print(worst[, c(keyCols, paste0(v1, c(".A", ".B")), "delta"), with = FALSE])
}

## 4. Which pixels were simulated at all ----------------------------------------------------
hr("cohortData and pixel coverage")
cdA <- readObj(dirA, "cohortData"); cdB <- readObj(dirB, "cohortData")
if (!is.null(cdA) && !is.null(cdB)) {
  message("cohort rows      A: ", nrow(cdA), "  B: ", nrow(cdB), "  (", nrow(cdB) - nrow(cdA), ")")
  message("pixelGroups      A: ", uniqueN(cdA$pixelGroup), "  B: ", uniqueN(cdB$pixelGroup))
  message("ecoregionGroups  A: ", uniqueN(cdA$ecoregionGroup), "  B: ", uniqueN(cdB$ecoregionGroup))
}
pgA <- readObj(dirA, "pixelGroupMap"); pgB <- readObj(dirB, "pixelGroupMap")
if (!is.null(pgA) && !is.null(pgB)) {
  nA <- sum(!is.na(terra::values(pgA))); nB <- sum(!is.na(terra::values(pgB)))
  message("simulated pixels A: ", nA, "  B: ", nB, "  (", nB - nA, ", ",
          sprintf("%+.1f%%", 100 * (nB - nA) / max(nA, 1)), ")")
}

## 5. Pixel fates, the module's own accounting ----------------------------------------------
hr("pixelFateDT")
pfA <- readObj(dirA, "pixelFateDT"); pfB <- readObj(dirB, "pixelFateDT")
if (!is.null(pfA) && !is.null(pfB)) {
  fa <- as.data.table(pfA)[, .(fate, A = pixelsRemoved)]
  fb <- as.data.table(pfB)[, .(fate, B = pixelsRemoved)]
  print(merge(fa, fb, by = "fate", all = TRUE))
}

## 6. Imputed / inferred pixels -------------------------------------------------------------
hr("imputedPixID")
ipA <- readObj(dirA, "imputedPixID"); ipB <- readObj(dirB, "imputedPixID")
if (!is.null(ipA) && !is.null(ipB)) {
  message("imputed pixels   A: ", length(ipA), "  B: ", length(ipB))
  message("only in A: ", length(setdiff(ipA, ipB)), "  only in B: ", length(setdiff(ipB, ipA)))
}

message("\n", strrep("=", 72), "\ndone\n")

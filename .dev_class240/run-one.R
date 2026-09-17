#!/usr/bin/env Rscript
## Stage 2 of the class-240 regression harness: run the module once, against a named LandR
## build and a named module worktree, reading the snapshotted inputs from stage 1.
##
## Everything that could differ between two runs is pinned here: the inputs are files, the
## cache is off, the parameters are explicit, and the LandR build is chosen by path. What is
## left to vary is the code under test.
##
## Usage:
##   Rscript .dev_class240/run-one.R <label> <area> <landrPath> [modulePath]
## e.g.
##   Rscript .dev_class240/run-one.R baseline big ~/GitHub/LandR
##   Rscript .dev_class240/run-one.R forestland big ~/GitHub/worktrees/landr-221-forestland

suppressPackageStartupMessages({
  library(terra); library(sf); library(data.table); library(digest)
})

.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) < 3) stop("usage: run-one.R <label> <area> <landrPath> [modulePath]")
label     <- .args[[1]]
areaLabel <- .args[[2]]
landrPath <- normalizePath(.args[[3]], mustWork = TRUE)
modulePath <- if (length(.args) >= 4) normalizePath(.args[[4]], mustWork = TRUE) else normalizePath(file.path(getwd(), ".."), mustWork = TRUE)

ROOT <- Sys.getenv("CLASS240_ROOT", "/mnt/fast/class240-regression")
fixDir <- file.path(ROOT, "fixtures", areaLabel)
if (!dir.exists(fixDir)) stop("no fixtures for area '", areaLabel, "'; run prepare-inputs.R first")
outDir <- file.path(ROOT, "runs", paste0(areaLabel, "__", label))
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

options(
  spades.useRequire = FALSE,
  spades.moduleCodeChecks = FALSE,
  reproducible.useCache = FALSE,   ## the comparison must not read another run's cache
  reproducible.useMemoise = FALSE,
  reproducible.verbose = 0,
  LandR.assertions = TRUE
)

suppressPackageStartupMessages(pkgload::load_all(landrPath, quiet = TRUE, export_all = FALSE))
landrVersion <- as.character(packageVersion("LandR"))
message("LandR   : ", landrPath, "  (", landrVersion, ")")
message("module  : ", modulePath)
message("area    : ", areaLabel, "   label: ", label)

## Rasters are GeoTIFFs beside the manifest; everything else is an .rds. See the note in
## prepare-inputs.R: `wrap()` + `saveRDS` silently stores a reference to a temp file once a
## raster is disk-backed, and the fixture is dead as soon as that temp directory is.
readObj <- function(nm) {
  tif <- file.path(fixDir, paste0(nm, ".tif"))
  if (file.exists(tif)) return(terra::rast(tif))
  f <- file.path(fixDir, paste0(nm, ".rds"))
  if (!file.exists(f)) return(NULL)
  obj <- readRDS(f)
  if (inherits(obj, "PackedSpatRaster")) terra::unwrap(obj) else obj
}

manifest <- fread(file.path(fixDir, "manifest.csv"))
objects <- setNames(lapply(manifest$object, readObj), manifest$object)
objects <- objects[!vapply(objects, is.null, logical(1))]
message("inputs  : ", paste(names(objects), collapse = ", "))

## `rstLCC` is built HERE, by the LandR build under test, not read from the fixtures: for the
## #221 arms it is the treatment, since the rule lives in `prepInputs_*_LCC_FAO()`. NTEMS is
## used rather than SCANFI because the SCANFI land-cover ids are not publicly retrievable, and
## because NTEMS keeps classes 80/81. Arguments that exist only on the #221 branch are passed
## only when that build has them, so one script drives every arm.
ntemsYear <- as.integer(Sys.getenv("CLASS240_NTEMS_YEAR", "2000"))
lccArgs <- list(year = ntemsYear, to = objects$rasterToMatch, disturbedCode = 240,
                destinationPath = file.path(ROOT, "inputs"))
hasForestLand <- "forestLandFrom" %in% names(formals(LandR::prepInputs_NTEMS_LCC_FAO))
if (hasForestLand) {
  lccArgs$forestLandFrom <- Sys.getenv("CLASS240_FORESTLAND_FROM", "both")
  lccArgs$faoYear <- as.integer(Sys.getenv("CLASS240_FAO_YEAR", "2022"))
}
message("rstLCC  : NTEMS ", ntemsYear,
        if (hasForestLand) paste0(" | forestLandFrom=", lccArgs$forestLandFrom,
                                  " faoYear=", lccArgs$faoYear) else " | baseline rule (FAO 2019 code 2)")
objects$rstLCC <- do.call(LandR::prepInputs_NTEMS_LCC_FAO, lccArgs)

## Strip the category table. `prepInputs_NTEMS_LCC_FAO()` attaches one (prepInputs_NTEMS.R:126,
## `levels(out) <- cls`) while `prepInputs_SCANFI_LCC_FAO()` returns plain numeric codes, and
## `ecoregionProducer()` branches on `is.factor()` (LandR ecoregions.R:66) into
## `raster::factorValues()`, which errors on a terra-native categorical raster. The land-cover
## half of `ecoregion_lcc` then comes back NA and every group collapses to `<ecoregion>_NA`,
## i.e. the module silently estimates maxB / maxANPP per ecoregion with no LCC stratum at all.
## Stripping the levels hands over exactly what the SCANFI path would have.
## `levels(x) <- NULL` restores the land-cover CODES (20, 81, 230, 240 ...).
## Not `terra::catalyze()`: that returns a layer named "label" holding the category INDICES
## (1, 2, 3, 4), which is numeric and so passes the `is.factor()` branch, but then feeds wrong
## values into `ecoregionProducer()` -- `paddedFloatToChar()` stopped with
## "x%%1: non-numeric argument to binary operator". The replacement form is not exported under
## a qualified name, hence `library(terra)` at the top of this script.
levels(objects$rstLCC) <- NULL
saveRDS(terra::wrap(objects$rstLCC), file.path(outDir, "rstLCC-input.rds"))

t0 <- Sys.time()
sim <- SpaDES.core::simInit(
  times = list(start = 0, end = 1),
  modules = "Biomass_borealDataPrep",
  paths = list(
    modulePath = modulePath,
    inputPath = file.path(ROOT, "inputs"),
    cachePath = file.path(ROOT, "cache-runs", paste0(areaLabel, "__", label)),
    outputPath = outDir
  ),
  objects = objects,
  params = list(Biomass_borealDataPrep = list(
    ## Same source as the fixtures were built from. Every input is supplied, so nothing should
    ## be fetched here, but leaving these at the SCANFI default would let any unsupplied input
    ## reach for Drive ids that 404.
    dataSource = Sys.getenv("CLASS240_DATA_SOURCE", "KNN"),
    dataYear = as.integer(Sys.getenv("CLASS240_DATA_YEAR", "2001")),
    sppEquivCol = Sys.getenv("CLASS240_SPP_EQUIV_COL", "LandR"),
    .studyAreaName = paste0("class240_", areaLabel),
    .plots = NA,
    .useCache = FALSE,
    exportModels = "none"
  ))
)
sim <- SpaDES.core::spades(sim, .plotInitialTime = NA)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
message("ran in ", round(elapsed, 1), " min")

## Outputs worth diffing. The models themselves are deliberately not saved: they are large,
## and what matters is the parameters they produce.
OUT_OBJECTS <- c("speciesEcoregion", "cohortData", "ecoregion", "ecoregionMap", "pixelGroupMap",
                 "imputedPixID", "pixelFateDT", "minRelativeB", "species", "rstLCC",
                 "standAgeMap", "biomassMap")

saveObj <- function(obj, path) {
  if (inherits(obj, "SpatRaster")) saveRDS(terra::wrap(obj), path) else saveRDS(obj, path)
}

outManifest <- rbindlist(lapply(OUT_OBJECTS, function(nm) {
  obj <- tryCatch(sim[[nm]], error = function(e) NULL)
  if (is.null(obj)) return(NULL)
  f <- file.path(outDir, paste0(nm, ".rds"))
  saveObj(obj, f)
  data.table(object = nm, class = class(obj)[1], bytes = file.size(f),
             sha256 = digest::digest(file = f, algo = "sha256"))
}))

meta <- data.table(label = label, area = areaLabel, landrPath = landrPath,
                   landrVersion = landrVersion, modulePath = modulePath,
                   moduleSha = system(paste("git -C", shQuote(modulePath), "rev-parse --short HEAD"),
                                      intern = TRUE),
                   elapsedMin = round(elapsed, 2), when = format(Sys.time()))
fwrite(outManifest, file.path(outDir, "manifest.csv"))
fwrite(meta, file.path(outDir, "meta.csv"))
print(meta)
print(outManifest[, .(object, class, bytes)])
message("\nSaved to ", outDir)

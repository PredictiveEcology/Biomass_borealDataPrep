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

readObj <- function(nm) {
  f <- file.path(fixDir, paste0(nm, ".rds"))
  if (!file.exists(f)) return(NULL)
  obj <- readRDS(f)
  if (inherits(obj, "PackedSpatRaster")) terra::unwrap(obj) else obj
}

manifest <- fread(file.path(fixDir, "manifest.csv"))
objects <- setNames(lapply(manifest$object, readObj), manifest$object)
objects <- objects[!vapply(objects, is.null, logical(1))]
message("inputs  : ", paste(names(objects), collapse = ", "))

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

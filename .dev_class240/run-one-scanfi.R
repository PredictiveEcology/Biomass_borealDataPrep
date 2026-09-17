#!/usr/bin/env Rscript
## Stage 2, SCANFI edition: run the module once against a named LandR build and module
## worktree, on the SCANFI 2020 fixtures from prepare-inputs-scanfi.R.
##
## As in run-one.R, `rstLCC` is built here by the LandR build under test -- the forest-land
## rule lives in prepInputs_SCANFI_LCC_FAO() -- and every other input is a fixture. The SCANFI
## land cover itself comes through reproducible.urlRemap (arbutus), as in stage 1.
##
## Usage:
##   Rscript run-one-scanfi.R <label> <area> <landrPath> <modulePath>
## Environment:
##   CLASS240_STRATUM   "" (module default / not supported), "landcover" or "siteComposition"
##   CLASS240_WETLAND   "1" to supply rstWetland from CWIM3A (needs LandR::prepInputs_CWIM)
##   CLASS240_FORESTLAND_FROM, CLASS240_FAO_YEAR   passed when the LandR build has them

suppressPackageStartupMessages({
  library(terra); library(sf); library(data.table); library(digest)
})

.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) < 4) stop("usage: run-one-scanfi.R <label> <area> <landrPath> <modulePath>")
label      <- .args[[1]]
areaLabel  <- .args[[2]]
landrPath  <- normalizePath(.args[[3]], mustWork = TRUE)
modulePath <- normalizePath(.args[[4]], mustWork = TRUE)

ROOT <- Sys.getenv("CLASS240_ROOT", "/mnt/fast/class240-regression")
MANIFEST <- Sys.getenv("CLASS240_REMAP_MANIFEST", "~/arbutus_manifest_SCANFI_v2_clean.csv")
stratum <- Sys.getenv("CLASS240_STRATUM", "")
useWetland <- identical(Sys.getenv("CLASS240_WETLAND", "0"), "1")

fixDir <- file.path(ROOT, "fixtures", areaLabel)
if (!dir.exists(fixDir)) stop("no fixtures for area '", areaLabel, "'")
outDir <- file.path(ROOT, "runs", paste0(areaLabel, "__", label))
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

options(
  spades.useRequire = FALSE,
  spades.moduleCodeChecks = FALSE,
  reproducible.useCache = FALSE,
  reproducible.useMemoise = FALSE,
  reproducible.verbose = 0,
  reproducible.destinationPath = file.path(ROOT, "inputs"),
  ## empty manifest: no remap here, so the LandR build under test must supply its own
  reproducible.urlRemap = if (nzchar(MANIFEST)) reproducible::makeUrlRemap(utils::read.csv(path.expand(MANIFEST))),
  reproducible.useCOG = FALSE,
  reproducible.gdriveNoAuth = TRUE,
  LandR.assertions = TRUE
)
terra::terraOptions(memmax = as.numeric(Sys.getenv("CLASS240_MEMMAX_GB", "4")), todisk = TRUE)

suppressPackageStartupMessages(pkgload::load_all(landrPath, quiet = TRUE, export_all = FALSE))
landrVersion <- as.character(packageVersion("LandR"))
message("LandR   : ", landrPath, "  (", landrVersion, ")")
message("module  : ", modulePath)
message("area    : ", areaLabel, "   label: ", label, "   stratum: ",
        if (nzchar(stratum)) stratum else "(module default)", "   wetland: ", useWetland)

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

## rstLCC: the LandR build under test, SCANFI 2020, the module's own crop/mask/project form.
lccArgs <- list(year = 2020, disturbedCode = 240,
                maskTo = objects$studyArea_biomassParam,
                cropTo = objects$rasterToMatch_biomassParam,
                projectTo = objects$rasterToMatch_biomassParam,
                destinationPath = file.path(ROOT, "inputs"))
hasForestLand <- "forestLandFrom" %in% names(formals(LandR::prepInputs_SCANFI_LCC_FAO))
if (hasForestLand) {
  lccArgs$forestLandFrom <- Sys.getenv("CLASS240_FORESTLAND_FROM", "both")
  lccArgs$faoYear <- as.integer(Sys.getenv("CLASS240_FAO_YEAR", "2022"))
}
message("rstLCC  : SCANFI 2020 | ",
        if (hasForestLand) paste0("forestLandFrom=", lccArgs$forestLandFrom, " faoYear=", lccArgs$faoYear)
        else "baseline rule (FAO 2019 code 2)")
objects$rstLCC <- do.call(LandR::prepInputs_SCANFI_LCC_FAO, lccArgs)
levels(objects$rstLCC) <- NULL
terra::writeRaster(objects$rstLCC, file.path(outDir, "rstLCC-input.tif"), overwrite = TRUE)

if (useWetland) {
  if (!"prepInputs_CWIM" %in% getNamespaceExports("LandR")) {
    stop("CLASS240_WETLAND=1 but this LandR build has no prepInputs_CWIM()")
  }
  objects$rstWetland <- LandR::prepInputs_CWIM(to = objects$rasterToMatch_biomassParam)
  terra::writeRaster(objects$rstWetland, file.path(outDir, "rstWetland-input.tif"), overwrite = TRUE)
}

modParams <- list(
  dataSource = "SCANFI",
  dataYear = 2020L,
  sppEquivCol = "LandR",
  .studyAreaName = paste0("class240_", areaLabel),
  .plots = NA,
  .useCache = FALSE,
  exportModels = "none"
)
## Deterministic mode. Three things draw random numbers in this module: the 50-per-group
## subsample for the biomass model, the same for the age-imputation model, and
## convertUnwantedLCC()'s "nearestRandom". Unseeded, maxB varies run to run with a median CV of
## 10% (up to 62%) on the 60 km window, which swamps most between-arm differences. Seeding alone
## does not fix a comparison, because different arms subsample different data; so use all rows,
## and the deterministic neighbour rule.
deterministic <- identical(Sys.getenv("CLASS240_DETERMINISTIC", "0"), "1")
if (deterministic) {
  modParams$subsetDataBiomassModel <- 100000L
  modParams$subsetDataAgeModel <- 100000L
  modParams$LCCClassesToReplaceNNMethod <- "nearestWeighted"
  set.seed(1)
}
if (nzchar(stratum)) modParams$stratumType <- stratum
if (!useWetland && nzchar(stratum)) modParams$wetlandSource <- "none"

t0 <- Sys.time()
sim <- SpaDES.core::simInit(
  times = list(start = 0, end = 1),
  modules = "Biomass_borealDataPrep",
  paths = list(
    modulePath = dirname(modulePath),
    inputPath = file.path(ROOT, "inputs"),
    cachePath = file.path(ROOT, "cache-runs", paste0(areaLabel, "__", label)),
    outputPath = outDir
  ),
  objects = objects,
  params = list(Biomass_borealDataPrep = modParams)
)
sim <- SpaDES.core::spades(sim, .plotInitialTime = NA)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
message("ran in ", round(elapsed, 1), " min")

OUT_OBJECTS <- c("speciesEcoregion", "cohortData", "ecoregion", "ecoregionMap", "pixelGroupMap",
                 "imputedPixID", "pixelFateDT", "minRelativeB", "species", "rstLCC",
                 "rstWetland", "standAgeMap", "biomassMap")
saveObj <- function(obj, path) {
  if (inherits(obj, "SpatRaster")) {
    terra::writeRaster(obj, sub("\\.rds$", ".tif", path), overwrite = TRUE)
    ## keep the category table of ecoregionMap: it is what maps codes to ecoregionGroup
    if (!is.null(terra::levels(obj)[[1]]) && NROW(terra::levels(obj)[[1]])) {
      saveRDS(as.data.table(terra::levels(obj)[[1]]), sub("\\.rds$", "-levels.rds", path))
    }
    sub("\\.rds$", ".tif", path)
  } else {
    saveRDS(obj, path)
    path
  }
}
outManifest <- rbindlist(lapply(OUT_OBJECTS, function(nm) {
  obj <- tryCatch(sim[[nm]], error = function(e) NULL)
  if (is.null(obj)) return(NULL)
  f <- saveObj(obj, file.path(outDir, paste0(nm, ".rds")))
  data.table(object = nm, class = class(obj)[1], file = basename(f), bytes = file.size(f))
}))
moduleDir <- modulePath
meta <- data.table(label = label, area = areaLabel, landrPath = landrPath,
                   landrVersion = landrVersion,
                   landrSha = system(paste("git -C", shQuote(landrPath), "rev-parse --short HEAD"), intern = TRUE),
                   modulePath = moduleDir,
                   moduleSha = system(paste("git -C", shQuote(moduleDir), "rev-parse --short HEAD"), intern = TRUE),
                   stratum = stratum, wetland = useWetland, deterministic = deterministic,
                   elapsedMin = round(elapsed, 2), when = format(Sys.time()))
fwrite(outManifest, file.path(outDir, "manifest.csv"))
fwrite(meta, file.path(outDir, "meta.csv"))
print(meta)
message("\nSaved to ", outDir)

#!/usr/bin/env Rscript
## Stage 1, SCANFI edition: materialise the module's inputs from SCANFI v2 (2020) ONCE and
## snapshot them, exactly as prepare-inputs.R does for kNN 2001.
##
## Why a second script rather than a flag: the kNN script routes AROUND the module's default
## inputs (it supplies land cover and species cover itself) because the SCANFI Drive ids 404.
## Those same ids are mirrored on arbutus, and `reproducible.urlRemap` redirects them there, so
## here the module's OWN `.inputObjects` builds every input on its default SCANFI path. Nothing
## is routed around except the species list, which is named so the parameter fits have a real
## set of species to work with.
##
## The remap manifest is ~/arbutus_manifest_SCANFI_v2_clean.csv (read-only). NOT the copy
## FireSenseTesting/global.R reads (PredictiveEcology.org/scripts/, untracked): that one has no
## `type` column and so none of the nine `dir` rows, and without them listGoogleDriveFolder()
## cannot remap the SCANFI species FOLDER and falls back to an authenticated drive_ls(). Every
## id the two files share maps to the same URL.
##
## Usage:
##   Rscript .dev_class240/prepare-inputs-scanfi.R [big_s2020|small_s2020|all]

suppressPackageStartupMessages({
  library(terra); library(sf); library(data.table); library(digest)
})

.args <- commandArgs(trailingOnly = TRUE)
which_area <- if (length(.args)) .args[[1]] else "all"

ROOT <- Sys.getenv("CLASS240_ROOT", "/mnt/fast/class240-regression")
MODULE_PATH <- Sys.getenv("CLASS240_MODULE_PATH", unset = normalizePath(file.path(getwd(), ".."), mustWork = TRUE))
if (!dir.exists(file.path(MODULE_PATH, "Biomass_borealDataPrep"))) {
  stop("modulePath must be the directory containing Biomass_borealDataPrep/; got: ", MODULE_PATH)
}
LANDR_BASE <- Sys.getenv("CLASS240_LANDR_BASE", "~/GitHub/worktrees/landr-baseline-dev")
MANIFEST <- Sys.getenv("CLASS240_REMAP_MANIFEST", "~/arbutus_manifest_SCANFI_v2_clean.csv")

REF_RAS <- "/mnt/fast/inputs/CA_FAO_forest_2019/CA_FAO_forest_2019.tif"
CENTRE_LONLAT <- c(-89.20, 53.88)
DATA_YEAR <- 2020L
SPP_EQUIV_COL <- "LandR"
SPP_NAMES <- strsplit(Sys.getenv(
  "CLASS240_SPP_NAMES", "Pice_mar,Pice_gla,Pinu_ban,Popu_tre,Betu_pap,Abie_bal,Lari_lar"
), ",")[[1]]

## Same squares as the kNN fixtures, so the two editions are directly comparable.
## CLASS240_FIXTURE_SUFFIX keeps a new fixture set beside an old one (e.g. "_int").
.sfx <- Sys.getenv("CLASS240_FIXTURE_SUFFIX", "")
AREAS <- setNames(list(60000, 25000), paste0(c("big_s2020", "small_s2020"), .sfx))

## The ids this edition depends on. Checked up front so a stale manifest fails here, loudly,
## rather than as an anonymous Drive 404 an hour into the run.
## An empty CLASS240_REMAP_MANIFEST sets no remap here, so a LandR that installs its own SCANFI
## mirror on load (LandR#230) is what serves the files.
useManifest <- nzchar(MANIFEST)
remapTable <- if (useManifest) utils::read.csv(path.expand(MANIFEST)) else data.frame(id = character())
needIds <- c(landcover2020 = "1EGp7LUA7cXMR6KpXDmu617xsjwGM6aIx",
             age2020 = "1nXPS3bpFUESYieNfXO25OKlZJEgqtRnD",
             speciesDir2020 = "15T4HIFeqzwp0TuOuxmYoexuXdLFnCZBi")
missingIds <- needIds[!needIds %in% remapTable$id]
if (useManifest && length(missingIds)) {
  stop("remap manifest lacks: ", paste(names(missingIds), collapse = ", "), " (", MANIFEST, ")")
}

options(
  spades.useRequire = FALSE,
  spades.moduleCodeChecks = FALSE,
  reproducible.verbose = 1,
  reproducible.useMemoise = FALSE,
  reproducible.destinationPath = file.path(ROOT, "inputs"),
  reproducible.cachePath = file.path(ROOT, "cache-s2020"),
  reproducible.urlRemap = if (useManifest) reproducible::makeUrlRemap(remapTable),
  reproducible.useCOG = FALSE,
  ## With the remap every id the module asks for is served by arbutus; never prompt for a login.
  reproducible.gdriveNoAuth = TRUE
)
dir.create(file.path(ROOT, "inputs"), recursive = TRUE, showWarnings = FALSE)
terra::terraOptions(memmax = as.numeric(Sys.getenv("CLASS240_MEMMAX_GB", "4")), todisk = TRUE)

INPUT_OBJECTS <- c(
  "studyArea", "studyArea_biomassParam", "rasterToMatch", "rasterToMatch_biomassParam",
  "rawBiomassMap", "ecoregionLayer", "firePerimeters", "standAgeMap",
  "speciesLayers", "speciesTable", "columnsForPixelGroups", "sppEquiv", "sppColorVect",
  "sppNameVector"
)

RTM_RES <- 240
makeRTM <- function(sa) {
  saV <- terra::vect(sa)
  terra::mask(terra::rast(saV, res = c(RTM_RES, RTM_RES), vals = 1), mask = saV)
}
makeStudyArea <- function(side) {
  ref <- terra::rast(REF_RAS)
  ctr <- terra::project(terra::vect(matrix(CENTRE_LONLAT, ncol = 2), crs = "EPSG:4326"),
                        terra::crs(ref))
  xy <- terra::crds(ctr)
  e <- terra::ext(xy[1] - side / 2, xy[1] + side / 2, xy[2] - side / 2, xy[2] + side / 2)
  sf::st_as_sf(terra::as.polygons(e, crs = terra::crs(ref)))
}
## GeoTIFF, never wrap()+saveRDS -- see prepare-inputs.R for why.
saveObj <- function(obj, path) {
  if (inherits(obj, "SpatRaster")) {
    tif <- sub("\\.rds$", ".tif", path)
    terra::writeRaster(obj, tif, overwrite = TRUE, gdal = c("COMPRESS=ZSTD", "TILED=YES"))
    tif
  } else {
    saveRDS(obj, path)
    path
  }
}

prepareOne <- function(areaLabel) {
  side <- AREAS[[areaLabel]]
  message("\n=== preparing SCANFI ", DATA_YEAR, " inputs for '", areaLabel, "' (", side / 1000, " km)")
  outDir <- file.path(ROOT, "fixtures", areaLabel)
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

  suppressPackageStartupMessages(pkgload::load_all(LANDR_BASE, quiet = TRUE, export_all = FALSE))
  message("    LandR used for input derivation: ", LANDR_BASE, " (",
          as.character(packageVersion("LandR")), ")")

  sa <- makeStudyArea(side)
  rtm <- makeRTM(sa)

  sim <- SpaDES.core::simInit(
    times = list(start = 0, end = 1),
    modules = "Biomass_borealDataPrep",
    paths = list(
      modulePath = MODULE_PATH,
      inputPath = file.path(ROOT, "inputs"),
      cachePath = file.path(ROOT, "cache-s2020"),
      outputPath = file.path(ROOT, "scratch", areaLabel)
    ),
    objects = list(studyArea = sa, studyArea_biomassParam = sa,
                   rasterToMatch = rtm, rasterToMatch_biomassParam = rtm,
                   sppNameVector = SPP_NAMES),
    params = list(Biomass_borealDataPrep = list(
      dataSource = "SCANFI",
      dataYear = DATA_YEAR,
      sppEquivCol = SPP_EQUIV_COL,
      .studyAreaName = paste0("class240_", areaLabel),
      .plots = NA,
      .useCache = ".inputObjects"
    ))
  )

  manifest <- rbindlist(lapply(INPUT_OBJECTS, function(nm) {
    obj <- tryCatch(sim[[nm]], error = function(e) NULL)
    if (is.null(obj)) {
      message("    - ", nm, ": absent, skipped")
      return(NULL)
    }
    f <- saveObj(obj, file.path(outDir, paste0(nm, ".rds")))
    data.table(object = nm, class = class(obj)[1], file = basename(f),
               bytes = file.size(f), sha256 = digest::digest(file = f, algo = "sha256"))
  }))
  fwrite(manifest, file.path(outDir, "manifest.csv"))
  print(manifest[, .(object, class, bytes)])
  if (!is.null(sim$speciesLayers)) {
    message("    species layers: ", paste(names(sim$speciesLayers), collapse = ", "))
  }
  invisible(manifest)
}

todo <- if (identical(which_area, "all")) names(AREAS) else which_area
for (a in todo) prepareOne(a)
message("\nDone. Fixtures under ", file.path(ROOT, "fixtures"))

#!/usr/bin/env Rscript
## Stage 1 of the class-240 regression harness: materialise the module's inputs ONCE and
## snapshot them, so that every later run reads byte-identical inputs and touches no network.
##
## Why snapshot rather than let `.inputObjects` run each time: the module builds twelve inputs
## from remote sources (SCANFI, NTEMS, FAO, NFDB, ecodistricts). Re-deriving them per run makes
## the comparison depend on downloads, caches and the FAO/land-cover rule that is itself under
## test -- `rstLCC` is built by `prepInputs_SCANFI_LCC_FAO()`, the function LandR #221 changes.
## Snapshotting pins them as data.
##
## Study area: northwestern Ontario boreal (lon -89.20, lat 53.88), chosen by probing the
## national FAO 2019 and VLCE2 2000 rasters already on disk. In a 30 x 30 km window there:
## 76% treed, 269,992 pixels of class 81 (treed wetland), 36,914 pixels that the current rule
## makes 240 (3.7%) and 80,615 under the proposed rule (8.1%).
##
## Two sizes, because the module only runs its SCANFI year-fill loop when more than 1000 pixels
## are 240 (Biomass_borealDataPrep.R:809). At 240 m:
##   "big"   60 km -> ~62,500 pixels, ~2,300 of them 240 -> the loop fires
##   "small" 25 km -> ~10,800 pixels,   ~400 of them 240 -> it does not
## The difference between those two paths is the defect under test, so both must be covered.
##
## Usage:
##   Rscript .dev_class240/prepare-inputs.R [big|small|all]

suppressPackageStartupMessages({
  library(terra); library(sf); library(data.table); library(digest)
})

.args <- commandArgs(trailingOnly = TRUE)
which_area <- if (length(.args)) .args[[1]] else "all"

ROOT <- Sys.getenv("CLASS240_ROOT", "/mnt/fast/class240-regression")
## SpaDES resolves a module by directory name, so `modulePath` is the directory *containing*
## `Biomass_borealDataPrep/`. Run this from the module root, or set CLASS240_MODULE_PATH.
MODULE_PATH <- Sys.getenv("CLASS240_MODULE_PATH", unset = normalizePath(file.path(getwd(), ".."), mustWork = TRUE))
if (!dir.exists(file.path(MODULE_PATH, "Biomass_borealDataPrep"))) {
  stop("modulePath must be the directory containing Biomass_borealDataPrep/; got: ", MODULE_PATH)
}
LANDR_BASE <- Sys.getenv("CLASS240_LANDR_BASE", "~/GitHub/LandR")

## The national grid the probe used; the study areas are defined in its CRS.
REF_RAS <- "/mnt/fast/inputs/CA_FAO_forest_2019/CA_FAO_forest_2019.tif"
CENTRE_LONLAT <- c(-89.20, 53.88)

AREAS <- list(big = 60000, small = 25000)

options(
  spades.useRequire = FALSE,
  spades.moduleCodeChecks = FALSE,
  reproducible.verbose = 1,
  reproducible.useMemoise = FALSE,
  reproducible.destinationPath = file.path(ROOT, "inputs")
)

dir.create(file.path(ROOT, "inputs"), recursive = TRUE, showWarnings = FALSE)

## Objects the module would otherwise derive itself; each is snapshotted.
INPUT_OBJECTS <- c(
  "studyArea", "studyArea_biomassParam", "rasterToMatch", "rasterToMatch_biomassParam",
  "rawBiomassMap", "rstLCC", "ecoregionLayer", "firePerimeters", "standAgeMap",
  "speciesLayers", "speciesTable", "columnsForPixelGroups", "sppEquiv", "sppColorVect",
  "sppNameVector"
)

makeStudyArea <- function(side) {
  ref <- terra::rast(REF_RAS)
  ctr <- terra::project(terra::vect(matrix(CENTRE_LONLAT, ncol = 2), crs = "EPSG:4326"),
                        terra::crs(ref))
  xy <- terra::crds(ctr)
  e <- terra::ext(xy[1] - side / 2, xy[1] + side / 2, xy[2] - side / 2, xy[2] + side / 2)
  sf::st_as_sf(terra::as.polygons(e, crs = terra::crs(ref)))
}

saveObj <- function(obj, path) {
  if (inherits(obj, "SpatRaster")) {
    saveRDS(terra::wrap(obj), path)   ## wrap(): a SpatRaster does not survive a plain saveRDS
  } else {
    saveRDS(obj, path)
  }
}

prepareOne <- function(areaLabel) {
  side <- AREAS[[areaLabel]]
  message("\n=== preparing inputs for '", areaLabel, "' (", side / 1000, " km square)")
  outDir <- file.path(ROOT, "fixtures", areaLabel)
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

  sa <- makeStudyArea(side)

  ## Inputs are derived by the *baseline* LandR, so that the snapshot is the status quo.
  suppressPackageStartupMessages(pkgload::load_all(LANDR_BASE, quiet = TRUE, export_all = FALSE))
  message("    LandR used for input derivation: ", as.character(packageVersion("LandR")))

  sim <- SpaDES.core::simInit(
    times = list(start = 0, end = 1),
    modules = "Biomass_borealDataPrep",
    paths = list(
      modulePath = MODULE_PATH,
      inputPath = file.path(ROOT, "inputs"),
      cachePath = file.path(ROOT, "cache"),
      outputPath = file.path(ROOT, "scratch", areaLabel)
    ),
    objects = list(studyArea = sa, studyArea_biomassParam = sa),
    params = list(Biomass_borealDataPrep = list(
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
    f <- file.path(outDir, paste0(nm, ".rds"))
    saveObj(obj, f)
    data.table(object = nm, class = class(obj)[1], file = basename(f),
               bytes = file.size(f), sha256 = digest::digest(file = f, algo = "sha256"))
  }))
  fwrite(manifest, file.path(outDir, "manifest.csv"))
  print(manifest[, .(object, class, bytes)])

  ## What the land-cover map looks like going in: the baseline against which the
  ## #221 rule change is read.
  lccCounts <- as.data.table(terra::freq(sim$rstLCC))
  fwrite(lccCounts, file.path(outDir, "rstLCC-classCounts.csv"))
  message("    rstLCC class counts:")
  print(lccCounts)

  invisible(manifest)
}

todo <- if (identical(which_area, "all")) names(AREAS) else which_area
for (a in todo) prepareOne(a)
message("\nDone. Fixtures under ", file.path(ROOT, "fixtures"))

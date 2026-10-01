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

## NTEMS land-cover year used for the placeholder rstLCC. 2000 and 2001 are the VLCE2 years
## already on local disk, so nothing downloads.
NTEMS_YEAR <- as.integer(Sys.getenv("CLASS240_NTEMS_YEAR", "2000"))

## Source for biomass, stand age and species cover. kNN, not the module's SCANFI default:
## two SCANFI 2020 Drive ids return 404 (stand age, prepInputObjects.R:601; land cover,
## maps.R's id table), while the kNN branches fetch over public NFIS FTP and the kNN
## stand-age raster is already on local disk. kNN offers 2001 and 2011 only.
DATA_SOURCE <- Sys.getenv("CLASS240_DATA_SOURCE", "KNN")
DATA_YEAR <- as.integer(Sys.getenv("CLASS240_DATA_YEAR", "2001"))

## Naming convention for the species table; the module's own default (Biomass_borealDataPrep.R:193).
SPP_EQUIV_COL <- Sys.getenv("CLASS240_SPP_EQUIV_COL", "LandR")

## Species are named explicitly rather than derived from the species-presence raster. Derived,
## this window came back with a single species (Pice_mar), which leaves a parameter comparison
## almost nothing to move: the module fits maxB / maxANPP per species x ecoregion group. This
## is the usual boreal set, and every one has a kNN raster on NFIS FTP.
SPP_NAMES <- strsplit(Sys.getenv(
  "CLASS240_SPP_NAMES",
  "Pice_mar,Pice_gla,Pinu_ban,Popu_tre,Betu_pap,Abie_bal,Lari_lar"
), ",")[[1]]

AREAS <- list(big = 60000, small = 25000)

options(
  spades.useRequire = FALSE,
  spades.moduleCodeChecks = FALSE,
  reproducible.verbose = 1,
  reproducible.useMemoise = FALSE,
  reproducible.destinationPath = file.path(ROOT, "inputs"),
  ## `loadkNNSpeciesLayers()` (LandR maps.R:1096-1099) falls back to this option when no
  ## `cachePath` is in its dots, then calls `basename()` on it unguarded (:1207, :1214). In a
  ## fresh Rscript session the option is NULL, so the species step dies with
  ## "basename(cachePath): a character vector argument expected".
  reproducible.cachePath = file.path(ROOT, "cache")
)

dir.create(file.path(ROOT, "inputs"), recursive = TRUE, showWarnings = FALSE)

## Cap terra's memory. This machine is shared with long-running simulations, and an earlier
## attempt was killed by the OOM reclaim while cropping national 250 m species rasters. The
## work is a sequence of windowed reads, so it needs very little resident memory; LandR's own
## `prepInputs_NTEMS_LCC_FAO()` sets the same knobs for the same reason.
terra::terraOptions(memmax = as.numeric(Sys.getenv("CLASS240_MEMMAX_GB", "4")), todisk = TRUE)

## Objects the module would otherwise derive itself; each is snapshotted.
##
## `rstLCC` is deliberately NOT here. It is built by `prepInputs_*_LCC_FAO()`, which is the
## function LandR #221 changes, so it is the treatment, not a fixture: pinning it would make
## that comparison vacuous. Each run builds its own `rstLCC` with its own LandR build (see
## run-one.R) and both the land-cover map and the downstream outputs are compared.
INPUT_OBJECTS <- c(
  "studyArea", "studyArea_biomassParam", "rasterToMatch", "rasterToMatch_biomassParam",
  "rawBiomassMap", "ecoregionLayer", "firePerimeters", "standAgeMap",
  "speciesLayers", "speciesTable", "columnsForPixelGroups", "sppEquiv", "sppColorVect",
  "sppNameVector"
)

## The module's own rasterToMatch construction (Biomass_borealDataPrep.R:1575-1588), copied
## because supplying a land-cover raster means supplying the grid it is aligned to.
RTM_RES <- 240
makeRTM <- function(sa) {
  saV <- terra::vect(sa)
  if (terra::is.lonlat(saV)) {
    saV <- terra::project(saV, paste("+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77",
                                     "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))
  }
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

## Rasters are written as real GeoTIFFs, not `wrap()`ed into .rds.
##
## `terra::wrap()` embeds cell values only while the raster is small enough to sit in memory.
## With `todisk = TRUE` (set above, to keep this job inside its memory budget) every computed
## raster spills to a file under `tempdir()`, and `wrap()` then stores a ~3 KB *reference* to
## that file. The .rds looks fine, weighs nothing, and is unusable the moment the session's
## temp directory is gone -- which is exactly how an earlier snapshot produced fixtures whose
## `readRDS` failed with "[rast] file does not exist: /tmp/Rtmp.../spat_....tif".
##
## A GeoTIFF beside the manifest has none of that fragility, and can be inspected with gdalinfo.
saveObj <- function(obj, path) {
  if (inherits(obj, "SpatRaster")) {
    tif <- sub("\\.rds$", ".tif", path)
    terra::writeRaster(obj, tif, overwrite = TRUE,
                       gdal = c("COMPRESS=ZSTD", "TILED=YES"))
    tif
  } else {
    saveRDS(obj, path)
    path
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

  rtm <- makeRTM(sa)

  ## A land-cover raster has to be supplied here, but only to stop `.inputObjects` deriving one:
  ## the module's default route is `prepInputs_SCANFI_LCC_FAO(year = 2020)`, whose Drive id is
  ## not publicly retrievable (404), and which would also drop classes 80/81. This NTEMS build
  ## comes from rasters already on local disk. It is NOT snapshotted -- each run builds its own.
  rstLCCTmp <- LandR::prepInputs_NTEMS_LCC_FAO(
    year = NTEMS_YEAR, to = rtm, disturbedCode = 240,
    destinationPath = file.path(ROOT, "inputs")
  )

  ## Species cover must be supplied too. The module's species block calls
  ## `prepSpeciesLayers_SCANFI()` unconditionally (Biomass_borealDataPrep.R:1731-1744) --
  ## `dataSource` is consulted only for biomass (:1602) and stand age (:1687) -- and the SCANFI
  ## per-species ids 404 like the rest of the SCANFI 2020 set. kNN 2001 has 75 per-species
  ## rasters on public NFIS FTP, so build from those instead.
  ##
  ## `sppHarmonize()` is what the module itself uses (:1720) to derive the species table and
  ## the name and colour vectors; calling it here with NULLs reproduces that derivation, and
  ## all four objects are supplied together so the module uses exactly what the layers were
  ## built from.
  sppOuts <- LandR::sppHarmonize(
    sppEquiv = NULL, sppNameVector = SPP_NAMES, sppEquivCol = SPP_EQUIV_COL,
    sppColorVect = NULL, vegLeadingProportion = 0.8, studyArea = sa,
    dPath = file.path(ROOT, "inputs")
  )
  message("    species: ", paste(sppOuts$sppNameVector, collapse = ", "))
  speciesLayersTmp <- LandR::prepSpeciesLayers_KNN(
    destinationPath = file.path(ROOT, "inputs"),
    outputPath = file.path(ROOT, "inputs"),
    studyArea = sa,
    rasterToMatch = rtm,
    sppEquiv = sppOuts$sppEquiv,
    sppEquivCol = sppOuts$sppEquivCol,
    thresh = 10,
    year = DATA_YEAR
  )
  speciesLayersTmp <- LandR::NAcover2zero(speciesLayersTmp, rtm)

  sim <- SpaDES.core::simInit(
    times = list(start = 0, end = 1),
    modules = "Biomass_borealDataPrep",
    paths = list(
      modulePath = MODULE_PATH,
      inputPath = file.path(ROOT, "inputs"),
      cachePath = file.path(ROOT, "cache"),
      outputPath = file.path(ROOT, "scratch", areaLabel)
    ),
    objects = list(studyArea = sa, studyArea_biomassParam = sa,
                   rasterToMatch = rtm, rasterToMatch_biomassParam = rtm,
                   rstLCC = rstLCCTmp,
                   speciesLayers = speciesLayersTmp,
                   sppEquiv = sppOuts$sppEquiv,
                   sppNameVector = sppOuts$sppNameVector,
                   sppColorVect = sppOuts$sppColorVect),
    params = list(Biomass_borealDataPrep = list(
      ## kNN rather than the SCANFI default: the SCANFI 2020 stand-age and land-cover Drive
      ## ids both return 404 (prepInputObjects.R:601 and maps.R's LCC table), whereas the kNN
      ## branches fetch biomass and stand age over public NFIS FTP, and the kNN stand-age
      ## raster is already on local disk.
      dataSource = DATA_SOURCE,
      dataYear = DATA_YEAR,
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

  ## No rstLCC counts here: the land-cover map is built per run, by the LandR build under
  ## test, and compared there (see compare-runs.R).

  invisible(manifest)
}

todo <- if (identical(which_area, "all")) names(AREAS) else which_area
for (a in todo) prepareOne(a)
message("\nDone. Fixtures under ", file.path(ROOT, "fixtures"))

Known issues: <https://github.com/PredictiveEcology/Biomass_borealDataPrep/issues>

version 1.5.13
=============

## bug fixes
* A study area in which every species is present in every pixel of its `ecoregionGroup` (e.g., one
  tree species, cover 100 everywhere) no longer stops with "contrasts can be applied only to factors
  with 2 or more levels". Rows with 100% presence are left out of `coverModel` and given probability
  1; when that was every row, the model was fitted to an empty table. The cover model is now skipped
  in that case and every row gets probability 1, as those rows already did. The step is in
  `estimateCoverModel()` (`R/estimateCoverModel.R`), with tests.

version 1.5.12
=============

## dependency changes
* requires `LandR (>= 1.2.0.9005)`, which added the `method` argument to
  `convertUnwantedLCC()` and removed that function's deterministic lowest-class tie-break.
  1.2.0.9004 had also dropped the `newPossLCC` column this module uses to write
  replacement classes back into `rstLCCAdj`, where the `is.null()` guard silently stopped
  firing and left `rstLCCAdj` (hence `ecoregionMap`) showing the un-replaced classes;
  1.2.0.9005 restores it.

## new features
* new parameter `LCCClassesToReplaceNNMethod`, passed to `LandR::convertUnwantedLCC()` as
  its `method`. Default `"nearestRandom"` preserves the stochastic allocation this module
  has always had. `"nearestWeighted"` gives the same abundance weighting but keys the draw
  on the pixel's ground position, so it is deterministic without a seed and a grid-aligned
  crop of the study area agrees with the full extent.

# Biomass_borealDataPrep 1.5.11 (2026-06-02)

* New `landis` mode parameter (default `FALSE`): when enabled, the forested LCC classes are collapsed to a single class in `createBiomass_coreInputs()` (before `prepEcoregions`/`makePixelTable`) so `ecoregionGroup` is defined by ecoregion only, not ecoregion x LCC; `maxB`, `maxANPP`, and species establishment probability are then estimated per ecoregion, as expected by LANDIS-II Biomass Succession. Default `FALSE` preserves the standard ecoregion x LCC behaviour.
* Support study areas with more than one polygon: allow `sim$studyArea` with `NROW > 1` (only the local `studyArea` object needs a single feature).
* Switch serialization/caching from `qs` to `qs2`, with an accompanying refactor of `updateYoungBiomasses()`; add `googledrive` to `reqdPkgs` and use cloud storage for the module event.
* Address multiple:1 mapping in `sppEquiv`, and bump the minimum `LandR` version several times to pull in upstream bugfixes (`prepEcoregions`, fire effects).
* Fixes: use `terra::aggregate`; drop the invalid `startYear` argument from `prepInputsStandAgeMap`; move `sim` assignment inside its block in `.inputObjects`; complete the remaining `crayon` to `cli` conversions.

# Biomass_borealDataPrep 1.5.10 (2025-10-09)

* Add stand-age adjustment: cap pixel ages to species longevity via a new function, plus a diagnostic plot comparing age distributions before and after adjustment.
* Merge SCANFI support (updated data source with corresponding manual/Rmd rebuild and bibliography updates).
* Package metadata: switch from `crayon` to `cli`; general cleanup.

# Biomass_borealDataPrep 1.5.9 (2025-09-11)

* Overhaul the land-cover (LCC) adjustment function: fix the adjustment logic, remove hard-coded 34:36 class references, and improve handling of empty pixels.
* Fixes to the `dataYear` logic and to the `lastYrOnNTEMS` handling (NTEMS annual data).
* Caching improvements: add `.functionName` to several `Cache` calls for a speedup; address a `data.table` warning; update `sf` metadata.

# Biomass_borealDataPrep 1.5.8 (2025-05-27)

* Refactor biomass parameterization into a dedicated parameterization object (`biomassParam`), including guarding against negative biomass at low cover values (via `pmax`), fixing backwards and lat/lon tests, and ensuring the rasterToMatch is projected in metres.
* Fire behavior corrections and cleanup.
* Land cover: add treed wetland (class 81) to the forested classes; default `sppEquivCol` to `"LandR"`.
* Reprojection safety: disable `gdalwarp` for `postProcess` of raster indices to avoid unintended reprojection, and fix use of the temporary `reproducible.gdalwarp` option.
* Misc: `studyArea` is no longer a module output; skip redundant subsetting when `subsetDataBiomassModel` is `FALSE`; bump `LandR` for a `prepRTM` bugfix; replace deprecated `filename2`.

# Biomass_borealDataPrep 1.5.7 (2024-06-06)

* Use CRAN versions of package dependencies.
* Spinup: install `Biomass_core` via `SpaDES.project` (guarded to versions < 1.3.5) and qualify `getModule` with `SpaDES.project::`.
* Fixes: correct `initialEcoregionGroup` to `initialEcoregionCode`; stop forcing `reproducible.useTerra` before `prepInputsStandAgeMap`; fix a `.plots` bug; add missing `lme4::` namespace; drop use of the removed `SpaDES.core` `dotSeed`.
* Extensive GitHub Actions workflow and manual (Rmd) rebuild maintenance.

# Biomass_borealDataPrep 1.5.4 (2022-05-27)

* Baseline release (pre-2023): `terra` added as a dependency; comment fixes. Earlier history not detailed here.

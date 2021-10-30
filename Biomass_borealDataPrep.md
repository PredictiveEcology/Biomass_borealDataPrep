---
title: "Biomass_borealDataPrep"
author: "Yong Luo, Eliot McIntire, Ceres Barros, Alex Chubaty"
date: "24 September 2019, updated 26 Jul 2021"
output:
  html_document:
    keep_md: yes
editor_options:
  chunk_output_type: console
---



[![Gitter](https://badges.gitter.im/PredictiveEcology/LandR_Biomass.svg)](https://gitter.im/PredictiveEcology/LandR_Biomass?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

# Overview

This module converts open datasets that are available for all of Canada's forests, into the input requirements for LandR Biomass, a forest landscape succession model derived from the Landis-II Biomass Succession Extension model (v3.2; hereafter LBSE). 
This has been designed and tested for some parts of the Western Boreal Forest. 

Specifically, it takes the stand biomass, stand age (defaulting to the kNN biomass/age maps), land-cover (Land Cover of Canada map by default) and ecological zonation maps of Canada (ecodistricts by default), as well as species specific % cover maps of Canada (defaulting to the kNN maps) and derives LandR Biomass parameter estimates from these data.
Species traits are taken from those used by Dominic Cyr for LANDIS-II simulations, with some exceptions (see below).

Keeping data preparation outside of the LandR `Biomass_core` module maintains the modularity of the LandR modules.

# Functioning

After downloading all the necessary data, the module proceeds to prepare the necessary objects and parameters for the simulation (see 'Data dependencies' below).
Depending on the objects, some are parameterized using empirical models, others based on literature (e.g. longevity values for western boreal species taken from Burton & Cumming (1995) - see `?LandR::speciesTableUpdate`), or expert knowledge (e.g., `sufficientLight` values adjusted to reflect western boreal forest succession dynamics)

## Filling data gaps

* mismatches between stand age, stand biomass and species cover are dealt with by trusting species % cover first.
  If `cover > 0` and `age == 0`, `age` is empirically estimated using a linear mixed effects model:

    `age ~ B * species + (1 | initialEcoregionCode) + cover`

  If `cover == 0`, `age == 0` and `biomass > 0`, `biomass` is set to `0`.

* species `longevity` parameters are taken from published numbers from closest regions (Burton & Cumming 1995)

* the user can choose to replace certain land cover classes with others (parameter `LCCClassesToReplaceNN`). For example, in a "Natural Range of Variation" project, urban pixels can be changed to forested. Since the currently available simulation modules do not drive land cover change, species-ecoregion traits will be spatially non-varying. This means that initial land cover classes that are "transient", such as "recent burn" are likely not appropriate to have "permanent" traits associated with them. This module offers the ability to estimate the parameters for these pixels *as if they were a neighbouring land cover class*. In these two cases, the module can replace these classes with a neighbouring forest type (defined by the parameter `forestedLCCClasses` and probabilistically chosen in proportion to its neighbourhood abundance up to a radius of 5 times the pixel resolution).


## Parameterization

LandR-Biomass, like LBSE, has species-level traits (invariant species tratis) and species-ecolocation level traits (spatially varying traits). Ecolocations are defined, by default, as *land cover class* X *ecodistrict* (Canadian, federal definition) combinations, and the traits are estimated by *species* X *ecolocation* combinations. This can mean that there are some species-ecolocation traits that are unique at the pixel level as there may be rare combinations of, say, a land cover type, for a given species, in a given ecodistrict. 
**Note that `ecolocations` are called `ecoregionGroup`'s across LandR Biomass objects, and species `speciesCode`**. For simplicity, below we refer to these as "ecolocation" and "species", but bear in mind that when passing formulas to the modules LandR object nomenclature will need to be used. 

### Spatially varying species traits

* establishment probabilities: species establishment probabilities (or sometimes `SEP` in LANDIS) by ecolocation are empirically estimated using species cover data (converted to presence/absence) using a GLMEM defined as (the equation is a parameter in the module that can be changed by the user):
  + `prob. presence ~ species + (1|ecolocation)`

* `maxB` is estimated empirically, using stand age, stand biomass per species, per ecolocation data, using an LMEM. Because the following equation represents a curve whose positive slope is decreasing with age (via `log(age)`), we estimate `maxB` as the expected `B` when `cover = 100` and `logAge = log(longevity)` for a given species in a given `ecolocation`. This equation is a parameter in the module that can be changed by the user.
  + `B ~ logAge * species + cover * species + (logAge + cover + species | ecolocation)`

* `maxANPP`: is defined as maxB/30 following LANDIS-II. But see below for updates.

All empirically estimated parameters *can* be estimated using data from a larger study area (`studyAreaLarge`) than the one used to run the simulation (`studyArea`), if the user provides such a polygon.

### Invariant species traits:

*  `growthcurve` and `mortalityshape` are by default taken from LANDIS-II parameters. However, they can be estimated from data by using the `Biomass_speciesParameters` module (see below for a brief explanation).

#### Updating species growth/mortality traits `Biomass_speciesParameters`
Briefly:
1.	We run ~41,000,000 hypothetical species with full factorial combinations of longevity, ratio of `maxANPP` to `maxBiomass`, `growthcurve`, `mortalityshape`

2.	We take the closest permanent and temporary sample plots (PSP) in or near the study area and find the hypothetical species in previous step that most closely matches the growth dynamics in the PSPs. This gives us the `growthcurve`, `mortalityshape`, and ratio of maximum biomass (`maxB`) to maximum ANPP (`maxANPP`) for each species in our study area

3.	We introduce a new parameter, `actualMaxBiomass`. We recognize that `maxB`, as obtained empirically above, cannot be easily reached in simulations because all reasonable values of `growthcurve`, `mortalityshape` and `longevity` prevent the equation from reaching `maxB` (it acts as an asymptote that is never approached). The `actualMaxBiomass` is then obtained by multiplying the empirically estimated `maxB` by the ratio between the `maxBiomass` parameter used for the simulations in step 1 and the maximum simulated biomass *actually* achieved in the simulations (of step 1).  We use this `actualMaxBiomass` so that the resulting non-linear growth curves will hit the the empirically estimated `maxB`. This adjustment effectively allows the growth equation's "maximum biomass" parameter to actually be the empirically estimated maximum biomass.  

4.	Species-specific `maxANPP` is estimated by multiplying the empirically estimated `maxB` (spatial) above and the ratio of the simulated `maxANPP` parameter (point 1) to the maximum simulated biomass (step 1) at the species level. 

## Get the module

See [SpaDES-modules repository](https://github.com/PredictiveEcology/SpaDES-modules) to see how to download this and other SpaDES modules.
Alternatively, it can be forked or cloned from github.com directly.

# Load libraries


```r
library(magrittr) # for %>% pipe
library(SpaDES)
```

# Set up paths

```r
moduleName <- "Biomass_borealDataPrep"
spadesModulesDirectory <- ".." # where the module will be located -- this is correct, if this module
                               # is an Rstudio project, i.e., one up from the project

inputPath <- file.path(dirname(spadesModulesDirectory), "inputs") %>% checkPath(create = TRUE)
outputPath <- file.path(dirname(spadesModulesDirectory), "outputs") 
cachePath = file.path(outputPath, "cache")
         
setPaths(cachePath = cachePath,
         modulePath = spadesModulesDirectory,
         inputPath = inputPath,
         outputPath = outputPath)
paths <- getPaths()
```

# Choose a study area


```r
library(raster)
# modulePath <- Cache(readline, paste0("Where is the module path? (e.g., ~/module, with no quotes).\n",
#                                      "Press Enter to accept the path in getPaths()$modulePath: "),
#                     cacheRepo = cachePath)
# setPaths(cachePath = cachePath, modulePath = modulePath)

## do you want to hand-draw a map or use defaults?
# - note that large areas will take longer to compute
handDrawMap <- FALSE

if (handDrawMap) {
  dev()
  clearPlot()
  canadaMap <- Cache(getData, 'GADM', country = 'CAN', level = 1, path = Paths$inputPath,
                     cacheRepo = getPaths()$cachePath, quick = FALSE)
  Plot(canadaMap, speedup = 5, visualSqueeze = 0.9) # 5 seemed optimal
  
  ## hand-drawn study area
  if (!exists("studyAreaLarge")) {
    message("Since there is no object called 'studyAreaLarge', please draw a study area with 10 points")
    severalrandompoints <- Cache(clickCoordinates, 10)
    # if(startsWith(attr(severalrandompoints, "tags"), "cache")) message("Taking studyAreaLarge from Cache")
    studyAreaLarge <- SpatialPolygons(list(Polygons(list(Polygon(severalrandompoints$coords)), ID = 1)),
                                          proj4string = crs(canadaMap))
  }
  Plot(studyAreaLarge, addTo = "canadaMap", col = "red")
  studyArea <- raster::buffer(studyAreaLarge, buffer = -200) ## make a smaller area within the first
} else {
  studyAreaLarge <- Cache(randomStudyArea, size = 1e8)
  studyArea <- raster::buffer(studyAreaLarge, buffer = -200)   ## make a smaller area within the first
}

times <- list(start = 0, end = 10)
modules <- list("Biomass_borealDataPrep")
objects <- list("studyAreaLarge" = studyAreaLarge,
                "studyArea" = studyArea) 

mySim <- simInit(times = times, #params = parameters, 
                 modules = modules, #, "Biomass_core"),
                 objects = objects, paths = getPaths())
```

# Run `spades`

This module is about data preparation, so there is no stochastic elements.
The `spades` call will only cause one event to occur (the `init` event)


```r
simOut <- spades(mySim, debug = TRUE)
```

# Visualize

The `Plot` function will visualize all known .quickPlot type objects, which includes `Raster*` and `SpatialPolygons*` objects.
After running this module, these are the outputs, which would likely be used as inputs to `Biomass_core`.


```r
dev()
clearPlot()

Plot(simOut)
```

# Downloads

During the `simInit` call, if the user does not provide alternatives for the expected inputs, the module will download 3 large `.tar` files (~2 GB each) and 1 `.zip` file (~45 MB) from the internet.

# Data dependencies

**NOTE:** all raster _inputs_ should be at the scale of `rasterToMatchLarge`/`studyAreaLarge` and all raster _outputs_ will be at the scale of `rasterToMatch`/`studyArea.`

## Module parameters


```
## defineParameter: 'speciesUpdateFunction' is not of specified type 'list'.
```

```
## defineParameter: 'pixelGroupAgeClass' is not of specified type 'numeric'.
```

```
## defineParameter: 'speciesUpdateFunction' is not of specified type 'list'.
```

```
## defineParameter: '.plotInitialTime' is not of specified type 'numeric'.
```



|paramName                 |paramClass |default      |min   |max  |paramDesc                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
|:-------------------------|:----------|:------------|:-----|:----|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|biomassModel              |call       |lme4::lm.... |NA    |NA   |Model and formula for estimating biomass (B) from ecoregionGroup (currently ecoregionLayer LandCoverClass), speciesCode, logAge (gives a downward curving relationship), and cover. Defaults to a LMEM, which can be slow if dealing with very large datasets (e.g. 36 000 points take 20min). For faster fitting try P(sim)$subsetDataBiomassModel == TRUE, or quote(RcppArmadillo::fastLm(formula = B ~ logAge speciesCode ecoregionGroup + cover speciesCode ecoregionGroup)). A custom model call can also be provided, as long as the 'data' argument is NOT included.                                                                                  |
|coverModel                |call       |glm, cbi.... |NA    |NA   |Model and formula used for estimating cover from ecoregion and speciesCode and potentially others. Defaults to a GLMEM if there are > 1 grouping levels. A custom model call can also be provided, as long as the 'data' argument is NOT included                                                                                                                                                                                                                                                                                                                                                                                                            |
|coverPctToBiomassPctModel |call       |glm, I(l.... |NA    |NA   |Model to estimate the relationship between % cover and % biomass, referred to as deciduousCoverDiscount. It is a number between 0 and 1 that translates %cover, as provided in several databases, to %biomass. It is assumed that all hardwoods are equivalent and all softwoods are equivalent and that %cover of hardwoods will be an overesimate of the %biomass of hardwoods. E.g., 30% cover of hardwoods might translate to 20% biomass of hardwoods. The reason this discount exists is because hardwoods in Canada have a much wider canopy than softwoods.                                                                                          |
|deciduousCoverDiscount    |numeric    |0.8418911    |NA    |NA   |This was estimated with data from NWT on March 18, 2020 and may or may not be universal. Will not be used if P(sim)$fitDeciduousCoverDiscount is TRUE                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
|fitDeciduousCoverDiscount |logical    |FALSE        |NA    |NA   |If TRUE, this will re-estimate deciduousCoverDiscount. This may be unstable and is not recommended currently. If FALSE, will use the current default                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
|dataYear                  |numeric    |2001         |NA    |NA   |Used to override the default 'sourceURL' of KNN datasets (species cover, stand biomass and stand age), which point to 2001 data, to fetch KNN data for another year. Currently, the only other possible year is 2011.                                                                                                                                                                                                                                                                                                                                                                                                                                        |
|ecoregionLayerField       |character  |             |NA    |NA   |the name of the field used to distinguish ecoregions, if supplying a polygon. Defaults to NULL and tries to use 'ECODISTRIC' where available (for legacy reasons), or the row numbers of sim$ecoregionLayer. If this field is not numeric, it will be coerced to numeric                                                                                                                                                                                                                                                                                                                                                                                     |
|exportModels              |character  |none         |NA    |NA   |Controls whether models used to estimate maximum B/ANPP ('biomassModel') and species establishment ('coverModel') probabilities are exported for posterior analyses or not. This may be important when models fail to converge or hit singularity (but can still be used to make predictions) and the user wants to investigate them further. Can be set to 'none' (no models are exported), 'all' (both are exported), 'biomassModel' or 'coverModel'. BEWARE: because this is intended for posterior model inspection, the models will be exported with data, which may mean very large simList(s)!                                                        |
|fireURL                   |character  |https://.... |NA    |NA   |A url to a fire database, such as the Canadian National Fire Database, that is a zipped shapefile with fire polygons, an attribute (i.e., a column) named 'Year'. If supplied (omitted with NULL or NA), this will be used to 'update' age pixels on standAgeMap with 'time since fire' as derived from this fire polygons map. Biomass is also updated in these pixels, when the last fire is more recent than 1986. If NULL or NAm, no age and biomass imputation will be done in these pixels.                                                                                                                                                            |
|fixModelBiomass           |logical    |FALSE        |NA    |NA   |should modelBiomass be fixed in the case of non-convergence? Only scaling of variables is implemented at this time                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
|forestedLCCClasses        |numeric    |1, 2, 3,.... |0     |NA   |The classes in the rstLCC layer that are 'treed' and will therefore be run in Biomass_core. Defaults to forested classes in LCC2010 map.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
|imputeBadAgeModel         |call       |lme4::lm.... |NA    |NA   |Model and formula used for imputing ages that are either missing or do not match well with Biomass or Cover. Specifically, if Biomass or Cover is 0, but age is not, then age will be imputed. Similarly, if Age is 0 and either Biomass or Cover is not, then age will be imputed                                                                                                                                                                                                                                                                                                                                                                           |
|LCCClassesToReplaceNN     |numeric    |             |NA    |NA   |This will replace these classes on the landscape with the closest forest class P(sim)$forestedLCCClasses. If the user is using the default 2005 data product for rstLCC, then users may wish to include 36 (cities -- if running a historic range of variation project), and 34:35 (burns) Since this is about estimating parameters for growth, it doesn't make any sense to have unique estimates for transient classes in most cases. If no classes are to be replaced, pass 'LCCClassesToReplaceNN' = numeric(0) when supplying parameters.                                                                                                              |
|minCoverThreshold         |numeric    |5            |0     |100  |Cover that is equal to or below this number will be omitted from the dataset                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|minRelativeBFunction      |call       |LandR::m.... |NA    |NA   |A quoted function that makes the table of min. relative B determining a stand shade level for each ecoregionGroup. Using the internal object `pixelCohortData` is advisable to access/use the list of ecoregionGroups per pixel. The function must output a data.frame with 6 columns, named 'ecoregionGroup' and 'X1' to 'X5', with one line per ecoregionGroup code ('ecoregionGroup'), and the min. relative biomass for each stand shade level X1-5. The default function uses values from LANDIS-II available at: https://github.com/dcyr/LANDIS-II_IA_generalUseFiles/blob/master/LandisInputs/BSW/biomass-succession-main-inputs_BSW_Baseline.txt%7E. |
|omitNonTreedPixels        |logical    |TRUE         |FALSE |TRUE |Should this module use only treed pixels, as identified by P(sim)$forestedLCCClasses?                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
|pixelGroupAgeClass        |numeric    |params(s.... |NA    |NA   |When assigning pixelGroup membership, this defines the resolution of ages that will be considered 'the same pixelGroup', e.g., if it is 10, then 6 and 14 will be the same                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|pixelGroupBiomassClass    |numeric    |100          |NA    |NA   |When assigning pixelGroup membership, this defines the resolution of biomass that will be considered 'the same pixelGroup', e.g., if it is 100, then 5160 and 5240 will be the same                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
|rmImputedPix              |logical    |FALSE        |NA    |NA   |Should sim$imputedPixID be removed from the simulation?                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
|speciesUpdateFunction     |list       |LandR::s.... |NA    |NA   |Unnamed list of (one or more) quoted functions that updates species table to customize values. By default, 'LandR::speciesTableUpdate' is used to change longevity and shade tolerance values, using values appropriate to Boreal Shield West (BSW), Boreal Plains (BP) and Montane Cordillera (MC) ecoprovinces (see ?LandR::speciesTableUpdate for details). Set to NULL if default trait values from 'speciesTable' are to be kept instead                                                                                                                                                                                                                |
|sppEquivCol               |character  |Boreal       |NA    |NA   |The column in sim$specieEquivalency data.table to use as a naming convention                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|speciesTableAreas         |character  |BSW, BP, MC  |NA    |NA   |One or more of the Ecoprovince short forms that are in the `speciesTable` file, e.g., BSW, MC etc. Default is good for Alberta and maybe other places.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
|subsetDataAgeModel        |numeric    |50           |NA    |NA   |the number of samples to use when subsampling the age data model and when fitting DeciduousCoverDiscount; Can be TRUE/FALSE/NULL or numeric; if TRUE, uses 50. If FALSE/NULL no subsetting is done.                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
|subsetDataBiomassModel    |numeric    |             |NA    |NA   |the number of samples to use when subsampling the biomass data model; Can be TRUE/FALSE/NULL or numeric; if TRUE, uses 50. If FALSE/NULL no subsetting is done.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
|successionTimestep        |numeric    |10           |NA    |NA   |defines the simulation time step, default is 10 years                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
|useCloudCacheForStats     |logical    |TRUE         |NA    |NA   |Some of the statistical models take long (at least 30 minutes, likely longer). If this is TRUE, then it will try to get previous cached runs from googledrive.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|.plotInitialTime          |numeric    |start(sim)   |NA    |NA   |This is here for backwards compatibility. Please use .plots                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
|.plots                    |character  |NA           |NA    |NA   |This describes the type of 'plotting' to do. See ?Plots for possible types. To omit, set to NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|.plotInterval             |numeric    |NA           |NA    |NA   |This describes the simulation time interval between plot events                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
|.saveInitialTime          |numeric    |NA           |NA    |NA   |This describes the simulation time at which the first save event should occur                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
|.saveInterval             |numeric    |NA           |NA    |NA   |This describes the simulation time interval between save events                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
|.seed                     |list       |             |NA    |NA   |Named list of seeds to use for each event (names). E.g., list('init' = 123) will set.seed(123) at the start of the init event and unset it at the end. Defaults to NULL, meaning that no seeds will be set                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|.studyAreaName            |character  |NA           |NA    |NA   |Human-readable name for the study area used. If NA, a hash of studyArea will be used.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
|.useCache                 |character  |.inputOb.... |NA    |NA   |Internal. Can be names of events or the whole module name; these will be cached by SpaDES                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |

## Inputs

This module has several input requirements. 
One is a study area, which should be provided as a `SpatialPolygonsDataFrame`, and named `studyAreaLarge`.
This should be inside the boundaries of the boreal forest of Canada. 
When first running the code in this `.Rmd` file, you will be prompted to draw a polygon if none is provided as an input.


```
## defineParameter: 'speciesUpdateFunction' is not of specified type 'list'.
```

```
## defineParameter: 'pixelGroupAgeClass' is not of specified type 'numeric'.
```

```
## defineParameter: 'speciesUpdateFunction' is not of specified type 'list'.
```

```
## defineParameter: '.plotInitialTime' is not of specified type 'numeric'.
```



|objectName            |objectClass              |desc                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |sourceURL                                                                                                                                                                                                      |
|:---------------------|:------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|cloudFolderID         |character                |The google drive location where cloudCache will store large statistical objects                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |NA                                                                                                                                                                                                             |
|columnsForPixelGroups |character                |The names of the columns in cohortData that define unique pixelGroups. Default is c('ecoregionGroup', 'speciesCode', 'age', 'B')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |NA                                                                                                                                                                                                             |
|ecoregionLayer        |SpatialPolygonsDataFrame |A SpatialPolygonsDataFrame that characterizes the unique ecological regions used to parameterize the biomass, cover, and species establishment probability models. It will be overlaid with landcover to generate classes for every ecoregion/LCC combination. It must have same extent and crs as studyAreaLarge if suppplied by user. It is superseded by sim$ecoregionRst if that object is supplied by the user                                                                                                                                                                                                                                                                                                                                                     |http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip                                                                                                                                         |
|ecoregionRst          |RasterLayer              |A raster that characterizes the unique ecological regions used to parameterize the biomass, cover, and species establishment probability models. If this object is provided, it will supersede sim$ecoregionLayer. It will be overlaid with landcover to generate classes for every ecoregion/LCC combination. It must have same extent and crs as rasterToMatchLarge if suppplied by user - use reproducible::postProcess. If it uses an attribute table, it must contain the field 'ecoregion' to represent raster values                                                                                                                                                                                                                                             |NA                                                                                                                                                                                                             |
|rstLCC                |RasterLayer              |A land classification map in study area. It must be 'corrected', in the sense that: 1) Every class must not conflict with any other map in this module (e.g., speciesLayers should not have data in LCC classes that are non-treed); 2) It can have treed and non-treed classes. The non-treed will be removed within this module if P(sim)$omitNonTreedPixels is TRUE; 3) It can have transient pixels, such as 'young fire'. These will be converted to a the nearest non-transient class, probabilistically if there is more than 1 nearest neighbour class, based on P(sim)$LCCClassesToReplaceNN. The default layer used, if not supplied, is Canada national land classification in 2010. The metadata (res, proj, ext, origin) need to match rasterToMatchLarge. |NA                                                                                                                                                                                                             |
|rasterToMatch         |RasterLayer              |A raster of the studyArea in the same resolution and projection as rawBiomassMap. This is the scale used for all outputs for use in the simulation. If not supplied will be forced to match the default rawBiomassMap.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |NA                                                                                                                                                                                                             |
|rasterToMatchLarge    |RasterLayer              |A raster of the studyAreaLarge in the same resolution and projection as rawBiomassMap. This is the scale used for all inputs for use in the simulation. If not supplied will be forced to match the default rawBiomassMap.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |NA                                                                                                                                                                                                             |
|rawBiomassMap         |RasterLayer              |total biomass raster layer in study area. Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived total aboveground biomass map from 2001 (in tonnes/ha), unless 'dataYear' != 2001. If necessary, biomass values are rescaled to match changes in resolution. See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata.                                                                                                                                                                                                                                                                                                                                                                                  |http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/canada-forests-attributes_attributs-forests-canada/2001-attributes_attributs-2001/NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif |
|speciesLayers         |RasterStack              |cover percentage raster layers by species in Canada species map. Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived species cover maps from 2001 using a cover threshold of 10 - see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata                                                                                                                                                                                                                                                                                                                                                                                                                                                            |http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/canada-forests-attributes_attributs-forests-canada/2001-attributes_attributs-2001/                                                                     |
|speciesTable          |data.table               |species attributes table, default is from Dominic Cyr and Yan Boulanger's project                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv                                                                                                                   |
|sppColorVect          |character                |named character vector of hex colour codes corresponding to each species                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |NA                                                                                                                                                                                                             |
|sppEquiv              |data.table               |table of species equivalencies. See LandR::sppEquivalencies_CA.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |NA                                                                                                                                                                                                             |
|sppNameVector         |character                |an optional vector of species names to be pulled from sppEquiv. If not provided, then species will be taken from the entire P(sim)$sppEquivCol in sppEquiv. See LandR::sppEquivalencies_CA.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |NA                                                                                                                                                                                                             |
|standAgeMap           |RasterLayer              |stand age map in study area. Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived biomass map from 2001, unless 'dataYear' != 2001. See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/canada-forests-attributes_attributs-forests-canada/2001-attributes_attributs-2001/NFI_MODIS250m_2001_kNN_Structure_Stand_Age_v1.tif                    |
|studyArea             |SpatialPolygonsDataFrame |Polygon to use as the study area. Defaults to an area in Southwestern Alberta, Canada.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |NA                                                                                                                                                                                                             |
|studyAreaLarge        |SpatialPolygonsDataFrame |multipolygon (larger area than studyArea) used for parameter estimation, with attribute LTHFC describing the fire return interval. Defaults to a square shapefile in Southwestern Alberta, Canada.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |NA                                                                                                                                                                                                             |

### Creates Inputs

Most of the inputs will be created automatically, if they are not provided by the user. 
The automatic creation will work for western boreal forests of Canada.
These are zip files and tar files that are available from various Natural Resources Canada web pages. 

## Outputs

This will show the outputs of this module, which can be used directly as the inputs for Biomass_core:


```
## defineParameter: 'speciesUpdateFunction' is not of specified type 'list'.
```

```
## defineParameter: 'pixelGroupAgeClass' is not of specified type 'numeric'.
```

```
## defineParameter: 'speciesUpdateFunction' is not of specified type 'list'.
```

```
## defineParameter: '.plotInitialTime' is not of specified type 'numeric'.
```



|objectName       |objectClass |desc                                                                                                                                                                                                                                                                                                                                        |
|:----------------|:-----------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|biomassMap       |RasterLayer |total biomass raster layer in study area, filtered for pixels covered by cohortData. Units in g/m2                                                                                                                                                                                                                                          |
|cohortData       |data.table  |initial community table, created from available biomass (g/m2), age and species cover data, as well as eco zonation information                                                                                                                                                                                                             |
|ecoregion        |data.table  |ecoregion look up table                                                                                                                                                                                                                                                                                                                     |
|ecoregionMap     |RasterLayer |ecoregion map that has mapcodes match ecoregion table and speciesEcoregion table                                                                                                                                                                                                                                                            |
|imputedPixID     |integer     |A vector of pixel IDs - matching rasterMatch IDs - that suffered data imputation. Data imputation may be in age (to match last fire event post 1950s, or 0 cover), biomass (to match fire-related imputed ages, correct for missing values or for 0 age/cover), land cover (to convert non-forested classes into to nearest forested class) |
|pixelGroupMap    |RasterLayer |initial community map that has mapcodes match initial community table                                                                                                                                                                                                                                                                       |
|pixelFateDT      |data.table  |A small table that keeps track of the pixel removals and cause. This may help diagnose issues related to understanding the creation of cohortData                                                                                                                                                                                           |
|minRelativeB     |data.frame  |define the cut points to classify stand shadeness                                                                                                                                                                                                                                                                                           |
|modelCover       |data.frame  |If P(sim)$exportModels is 'all', or 'cover', fitted biomass model, as defined by P(sim)$coverModel                                                                                                                                                                                                                                          |
|modelBiomass     |data.frame  |If P(sim)$exportModels is 'all', or 'biomass', fitted biomass model, as defined by P(sim)$biomassModel                                                                                                                                                                                                                                      |
|rawBiomassMap    |RasterLayer |total biomass raster layer in study area. Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived total aboveground biomass map (in tonnes/ha) from 2001. See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata                                                            |
|species          |data.table  |a table that has species traits such as longevity...                                                                                                                                                                                                                                                                                        |
|speciesEcoregion |data.table  |define the maxANPP, maxB and establishprob change with both ecoregion and simulation time                                                                                                                                                                                                                                                   |
|studyArea        |            |Polygon to use as the study area. Defaults to an area in Southwestern Alberta, Canada.                                                                                                                                                                                                                                                      |
|sufficientLight  |data.frame  |define how the species with different shade tolerance respond to stand shadeness. Table values follow LANDIS-II test traits available at: https://raw.githubusercontent.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/biomass-succession_test.txt                                   |


```r
## species table
simOut$speciesTable
```


```r
Plot(simOut$biomassMap)
simOut$studyAreaLarge <- spTransform(simOut$studyAreaLarge, crs(simOut$biomassMap))
Plot(simOut$studyAreaLarge, addTo = "simOut$biomassMap")
```

## Getting help

- <https://gitter.im/PredictiveEcology/LandR_Biomass>

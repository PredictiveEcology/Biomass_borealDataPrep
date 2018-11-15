## ------------------------------------------------------------------
## Function to load kNN species layers from online data repository
## ------------------------------------------------------------------

## dPath: directory to data folder
## rasterToMatch: passed to prepInputs
## studyArea: passed to prepInputs
## species is either a character vector of species names to download,
##    or a two-column matrix with the species names to download and final names, with column names = c("speciesNamesRaw", "speciesNamesEnd")
##    should two raw species names share the same final name, their biomass data will be considered as the "same species"
## thresh: is the minimum number of pixels where the species must have biomass > 0 to be considered present in the study area.
##    Defaults to 1
## url: is the source url for the data, passed to prepInputs.

loadkNNSpeciesLayers <- function(dPath, rasterToMatch, studyArea,
                                 speciesList = NULL, thresh = 1, url, cachePath, ...) {


  if(class(speciesList) == "matrix") {
    ## check column names
    if(!setequal(colnames(speciesList), c("speciesNamesRaw", "speciesNamesEnd")))
      stop("names(species) must be c('speciesNamesRaw', 'speciesNamesEnd'), for raw species names and final species names respectively")
  }

  # Changed by Eliot Oct 20 2018 -- can't start with untar because tar file may not be present
  suffix <- if (basename(cachePath) == "cache") paste0(as.character(ncell(rasterToMatch)),"px") else
    basename(cachePath)
  suffix <- paste0("_", suffix)

  ## Make sure raw names are compatible with kNN names
  kNNnames <- lapply(strsplit(speciesList[,1], "_"), function(x) {
    x[1] <- substring(x[1], 1, 4)
    x[2] <- paste0(toupper(substring(x[2], 1, 1)), substring(x[2], 2, 3))
    x
  })
  kNNnames <- sapply(kNNnames, function(x) paste(x, collapse = "_"))
  speciesList[, 1] <- kNNnames

  species1 <- Cache(loadFun, url = url, spp = speciesList, #[, "speciesNamesRaw"],
                    #loadFun,
                    dPath = dPath,
                    suffix = suffix,
                    studyArea = studyArea, rasterToMatch = rasterToMatch,
                    userTags = "kNN_SppLoad")
  # species1 <- Cache(lapply, seq_len(NROW(speciesList)),
  #                   spp = speciesList, #[, "speciesNamesRaw"],
  #                   loadFun, url = url, dPath = dPath,
  #                   suffix = suffix,
  #                   studyArea = studyArea, rasterToMatch = rasterToMatch,
  #                   userTags = "kNN_SppLoad")

  ## get all kNN species
  if (FALSE) { #TODO This no longer does all species }
    allSpp <- Cache(untar, tarfile = file.path(dPath, "kNN-Species.tar"), list = TRUE)
    allSpp <- allSpp %>%
      grep(".zip", ., value = TRUE) %>%
      sub("_v0.zip", "", .) %>%
      sub(".*Species_", "", .)


    ## check for missing species
    if(any(!speciesList[,1] %in% allSpp)) {
      warning("Some species not present in kNN database.
            /n  Check if this is correct")
      speciesList <- speciesList[speciesList[, 1] %in% allSpp,]
    }
  }

  names(species1) <- speciesList[, "speciesNamesRaw"]

  ## Sum species that share same final name
  if(any(duplicated(speciesList[, 2]))) {
    dubs <- unique(speciesList[duplicated(speciesList[, 2]), 2])   ## get the duplicated final names

    ## make a list of species that will be summed (those with duplicated final names)
    spp2sum <- lapply(dubs, FUN = function(x) {
      speciesList[speciesList[, 2] %in% x, 1]
    })

    names(spp2sum) = dubs

    for(i in 1:length(spp2sum)) {
      sumSpecies <- spp2sum[[i]]
      newLayerName <- names(spp2sum)[i]

      fname <- .suffix(file.path(dPath, paste0("KNN", newLayerName, ".tif")), suffix)
      a <- Cache(sumRastersBySpecies,
                 speciesLayers = species1[sumSpecies],
                 newLayerName = newLayerName,
                 filenameToSave = asPath(fname),
                 ...)
      a <- raster(fname) ## ensure a gets a filename

      ## replace spp rasters by the summed one
      species1[sumSpecies] <- NULL
      species1[[newLayerName]] <- a
    }
  }

  ## Rename species layers - note: merged species were renamed already
  nameReplace <- as.matrix(speciesList[,2])
  rownames(nameReplace) = speciesList[, 1]

  toReplace <- names(species1)[names(species1) %in% rownames(nameReplace)]
  names(species1)[names(species1) %in% toReplace] <- nameReplace[toReplace, 1]

  ## remove layers that have less data than thresh (i.e. spp absent in study area)
  ## count no. of pixels that have biomass
  layerData <- Cache(sapply, X = species1, function(x) sum(x[] > 0, na.rm = TRUE))

  ## remove layers that had < thresh pixels with biomass
  belowThresh <- layerData < thresh
  if (any(belowThresh))
    species1[belowThresh] <- NULL

  ## return stack and final species matrix
  list(specieslayers = stack(species1), speciesList = speciesList)
}


## ------------------------------------------------------------------
## Function to sum rasters of species layers
## ------------------------------------------------------------------

## speciesLayers: stack of species layers rasters
## layersToSum: names/indices of layers to be summed - optional
## filenameToSave: file path to save output raster
## newLayerName: name of the output raster layer

sumRastersBySpecies <- function(speciesLayers, layersToSum,
                                filenameToSave, newLayerName) {
  ras_out <- raster::calc(raster::stack(speciesLayers[layersToSum]), sum)
  names(ras_out) <- newLayerName
  writeRaster(ras_out, filename = filenameToSave, datatype = "INT2U", overwrite = TRUE)
  ras_out # Work around for Cache
}

loadFun <- function(speciesListIndex, spp, suffix, url, dPath,
                    studyArea, rasterToMatch) {

  if (is.null(spp)) {
    knownSp <- c("Abie_Ama", "Abie_Bal", "Abie_Gra", "Abie_Las", "Abie_Spp",
                 "Acer_Cir", "Acer_Mac", "Acer_Neg", "Acer_Pen", "Acer_Rub", "Acer_Sac",
                 "Acer_Sah", "Acer_Spi", "Acer_Spp", "Alnu_Inc_Rug", "Alnu_Inc_Ten",
                 "Alnu_Inc", "Alnu_Rub", "Alnu_Spp", "Arbu_Men", "Asim_Tri", "Betu_All",
                 "Betu_Pap", "Betu_Pop", "Betu_Spp", "Carp_Car", "Cary_Cor", "Cast_Den",
                 "Cham_Noo", "Crat_Spp", "Fagu_Gra", "Frax_Ame", "Frax_Nig", "Frax_Pen_Sub",
                 "Frax_Pen", "Frax_Spp", "Generic_BroadLeaf_Spp", "Generic_NeedleLeaf_Spp",
                 "Gled_Tri", "Jugl_Cin", "Jugl_Nig", "Juni_Vir", "Lari_Kae", "Lari_Lar",
                 "Lari_Lya", "Lari_Occ", "Lari_Spp", "Malu_Fus", "Malu_Spp", "Ostr_Vir",
                 "Pice_Abi", "Pice_Eng_Gla", "Pice_Eng", "Pice_Gla", "Pice_Mar",
                 "Pice_Rub", "Pice_Sit", "Pice_Spp", "Pinu_Alb", "Pinu_Ban", "Pinu_Con_Lat",
                 "Pinu_Con", "Pinu_Fle", "Pinu_Mon", "Pinu_Pon", "Pinu_Res", "Pinu_Rig",
                 "Pinu_Spp", "Pinu_Str", "Pinu_Syl", "Plat_Occ", "Popu_Bal", "Popu_Del",
                 "Popu_Gra", "Popu_Spp", "Popu_Tre", "Popu_Tri", "Prun_Pen", "Prun_Ser",
                 "Prun_Vir", "Pseu_Men_Gla", "Pseu_Men_Men", "Pseu_Men", "Quer_Alb",
                 "Quer_Bic", "Quer_Gar", "Quer_Mac", "Quer_Rub", "Robi_Pse", "Sali_Beb",
                 "Sali_Nig", "Sali_Spp", "Sass_Alb", "Sorb_Ame", "Sorb_Dec", "Sorb_Spp",
                 "Thuj_Occ", "Thuj_Pli", "Thuj_Spp", "Tili_Ame", "Tsug_Can", "Tsug_Het",
                 "Tsug_Mer_Het", "Tsug_Mer", "Tsug_Spp", "Ulmu_Ame", "Ulmu_Rub",
                 "Ulmu_Spp", "Ulmu_Tho")
    stop("This loadFun has not been tested for all species. Please specify the actual species desired by name",
         " Known species are:\n", paste(knownSp, collapse = "\n"))
  }
  archive <- asPath("kNN-Species.tar")
  ## check if species is a vector/matrix
  if (is.null(spp)) {
      ## set to NULL so prepInputs extracts all of them
    targetFile <- NULL

    # just get tar file, no crop/reproject etc. Too many
    tarFile <- prepInputs(
      targetFile = targetFile,
      url = url,
      archive = archive,
      destinationPath = asPath(dPath),
      fun = "raster::raster")#,
      #studyArea = studyArea,
      #rasterToMatch = rasterToMatch,
      #method = "bilinear",
      #datatype = "INT2U",
      #filename2 = postProcessedFilename


    ## make a matrix of raw and final species names
    spp <-  matrix(data = rep(spp, 2),
                   nrow = length(spp), ncol = 2, byrow = FALSE)
    colnames(spp) = c("speciesNamesRaw", "speciesNamesEnd")

  } else if (class(spp) == "matrix") {
    ## check column names
    if(!setequal(colnames(spp), c("speciesNamesRaw", "speciesNamesEnd")))
      stop("names(species) must be c('speciesNamesRaw', 'speciesNamesEnd'), for raw species names and final species names respectively")
    targetFiles <- paste0("NFI_MODIS250m_kNN_Species_", spp[, "speciesNamesRaw"], "_v0.tif")
    names(targetFiles) <- targetFiles
    archives <- cbind(archive1 = archive, archive2 = paste0("NFI_MODIS250m_kNN_Species_", spp[, "speciesNamesRaw"], "_v0.zip"))
    archives <- split(archives, archives[, "archive2"])
  } else stop("species must be a character vector or a two-column matrix")

  postProcessedFilenames <- .suffix(targetFiles, suffix = suffix)


  species1 <- Map(targetFile = targetFiles, archive = archives,
                  filename2 = postProcessedFilenames,
                  MoreArgs = list(url = url,
                                  destinationPath = asPath(dPath),
                                  fun = "raster::raster",
                                  studyArea = studyArea,
                                  rasterToMatch = rasterToMatch,
                                  method = "bilinear",
                                  datatype = "INT2U"
                  ),
                  prepInputs)

  # species1 <- prepInputs(
  #   targetFile = targetFile,
  #   url = url,
  #   archive = archive,
  #   destinationPath = asPath(dPath),
  #   fun = "raster::raster",
  #   studyArea = studyArea,
  #   rasterToMatch = rasterToMatch,
  #   method = "bilinear",
  #   datatype = "INT2U",
  #   filename2 = postProcessedFilename
  #   )

  names(species1) <- spp[, "speciesNamesRaw"]
  return(species1)
}

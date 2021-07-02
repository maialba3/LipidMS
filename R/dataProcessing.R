# dataProcessing
#' Process mzXML files individually: peakpicking and isotope annotation
#'
#' Process mzXML files individually: peakpicking and isotope anotation
#'
#' @param file file path.
#' @param acquisitionmode character value: MS, DIA or DDA.
#' @param polarity character value: negative or positive.
#' @param dmzagglom mz tolerance (in ppm) used for partitioning and clustering.
#' @param drtagglom RT window used for partitioning (in seconds). 
#' @param drtclust RT window used for clustering (in seconds). 
#' @param minpeak minimum number of measurements required for a peak. 
#' @param drtgap maximum RT gap length to be filled (in seconds). 
#' @param drtminpeak minimum RT width of a peak (in seconds). At least minpeak 
#' within the drtminpeak window are required to define a peak.
#' @param drtmaxpeak maximum RT width of a single peak (in seconds).
#' @param recurs maximum number of peaks within one EIC.
#' @param sb signal-to-base ratio.
#' @param sn signal-to-noise ratio.
#' @param minint minimum intensity of a peak. 
#' @param weight weight for assigning measurements to a peak.
#' @param dmzIso mass tolerance for isotope matching.
#' @param drtIso time window for isotope matching.
#'
#' @return an msobject that contains metadata of the mzXML file, raw data and
#' extracted peaks.
#'
#' @details It is important that mzXML files are centroided. 
#' 
#' This function executes 2 steps: 1) peak-picking based on enviPick
#' package and 2) isotope annotation.
#' 
#' Numeric arguments accept one or two values for MS1 and MS2, respectively.
#' 
#' 
#' @seealso \link{batchdataProcessing} and \link{setmsbatch}
#'
#' @examples
#' \dontrun{
#' msobject <- dataProcessing("input_file.mzXML", acquisitionmode="DIA", polarity,
#' dmzagglom = 25, drtagglom = 500, drtclust = 60, minpeak = c(5, 3),
#' drtgap = 5, drtminpeak = 20, drtmaxpeak = 100, recurs = 5, sb = c(3, 2),
#' sn = 2, minint = c(1000, 100), weight = 2, dmzIso = 10, drtIso = 5)
#' }
#'
#' @references Peak-picking algorithm has been imported from enviPick R-package:
#' https://cran.r-project.org/web/packages/enviPick/index.html
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@iislafe.es>
dataProcessing <- function(file, 
                           acquisitionmode, 
                           polarity,
                           dmzagglom = 15,
                           drtagglom = 500,
                           drtclust = 100,
                           minpeak = c(5, 3),
                           drtgap = 10,
                           drtminpeak = c(15, 15),
                           drtmaxpeak = c(100, 200),
                           recurs = 5,
                           sb = c(3, 2),
                           sn = 2,
                           minint = c(1000, 100),
                           weight = c(2, 3),
                           dmzIso = 5,
                           drtIso = 5){
  #============================================================================#
  # check arguments
  #============================================================================#
  if (length(file) > 1){
    stop("dataProcessing just admits one file to process. 
         \nIn case you want to process multiple files, use batchdataProcessing instead.")
  }
  if (!file.exists(file)){
    stop("File doesn't exist")
  }
  if (missing(acquisitionmode)){
    stop("acquisitionmode argument is required: it must be \"MS\" \"DIA\" or \"DDA\"")
  }
  if (acquisitionmode %in% c("ms", "dia", "dda", "MS", "DIA", "DDA")){
    acquisitionmode <- toupper(acquisitionmode)
  }
  if (!tolower(polarity) %in% c("positive", "negative")){
    stop("Polarity must be set to positive or negative")
  } else {
    polarity <- tolower(polarity)
  }
  if(length(dmzagglom)>1){dmzagglom1<-dmzagglom[1];dmzagglom2<-dmzagglom[2]}else{dmzagglom1<-dmzagglom2<-dmzagglom}
  if(length(drtagglom)>1){drtagglom1<-drtagglom[1];drtagglom2<-drtagglom[2]}else{drtagglom1<-drtagglom2<-drtagglom}
  if(length(drtclust)>1){drtclust1<-drtclust[1];drtclust2<-drtclust[2]}else{drtclust1<-drtclust2<-drtclust}
  if(length(minpeak)>1){minpeak1<-minpeak[1];minpeak2<-minpeak[2]}else{minpeak1<-minpeak2<-minpeak}
  if(length(drtgap)>1){drtgap1<-drtgap[1];drtgap2<-drtgap[2]}else{drtgap1<-drtgap2<-drtgap}
  if(length(drtminpeak)>1){drtminpeak1<-drtminpeak[1];drtminpeak2<-drtminpeak[2]}else{drtminpeak1<-drtminpeak2<-drtminpeak}
  if(length(drtmaxpeak)>1){drtmaxpeak1<-drtmaxpeak[1];drtmaxpeak2<-drtmaxpeak[2]}else{drtmaxpeak1<-drtmaxpeak2<-drtmaxpeak}
  if(length(recurs)>1){recurs1<-recurs[1];recurs2<-recurs[2]}else{recurs1<-recurs2<-recurs}
  if(length(sb)>1){sb1<-sb[1];sb2<-sb[2]}else{sb1<-sb2<-sb}
  if(length(sn)>1){sn1<-sn[1];sn2<-sn[2]}else{sn1<-sn2<-sn}
  if(length(minint)>1){minint1<-minint[1];minint2<-minint[2]}else{minint1<-minint2<-minint}
  if(length(weight)>1){weight1<-weight[1];weight2<-weight[2]}else{weight1<-weight2<-weight}
  if(length(dmzIso)>1){dmzIso1<-dmzIso[1];dmzIso2<-dmzIso[2]}else{dmzIso1<-dmzIso2<-dmzIso}
  if(length(drtIso)>1){drtIso1<-drtIso[1];drtIso2<-drtIso[2]}else{drtIso1<-drtIso2<-drtIso}
  
  
  #============================================================================#
  # read file and filter scans by polarity if required
  #============================================================================#
  cat(paste(c("\n", file), collapse=""))
  cat("\n Reading MS file...")
  msobject <- readMSfile(file, polarity)
  msobject$metaData$generalMetadata$acquisitionmode <- acquisitionmode
  cat("OK")
  
  #============================================================================#
  # Peak-picking: based on enviPick algorithm
  #============================================================================#
  cat("\n Searching for features...")
  ##############################################################################
  # msLevel 1
  if("MS1" %in% names(msobject$rawData)){
    cat("\n   Processing MS1...")
    for (cE in names(msobject$rawData$MS1)){
      cat("\n     partitioning...")
      msobject <- partitioning(msobject, dmzagglom = dmzagglom1, drtagglom = drtagglom1,
                               minpeak = minpeak1, mslevel = "MS1", cE = cE)
      cat("OK")
      cat("\n     clustering...")
      msobject <- clustering(msobject, dmzagglom = dmzagglom1, drtclust = drtclust1,
                             minpeak = minpeak1, mslevel = "MS1", cE = cE)
      cat("OK")
      cat("\n     detecting peaks...")
      msobject <- peakdetection(msobject, minpeak = minpeak1, drtminpeak = drtminpeak1,
                                drtmaxpeak = drtmaxpeak1, drtgap = drtgap1,
                                recurs = recurs1, weight = weight1,
                                sb = sb1, sn = sn1, minint = minint1,
                                ended = 2, mslevel = "MS1", cE = cE)
      cat("OK")
      msobject$peaklist$MS1 <- do.call(rbind, msobject$peaklist[["MS1"]])
      rownames(msobject$peaklist$MS1) <- msobject$peaklist$MS1$peakID
      msobject$rawData$MS1 <- do.call(rbind, msobject$rawData$MS1)
      msobject$rawData$MS1 <- msobject$rawData$MS1[,c("mz", "RT", "int", "peak", "Scan")]
      colnames(msobject$rawData$MS1) <- c("mz", "RT", "int", "peakID", "Scan")
      # msobject$rawData$MS1 <- msobject$rawData$MS1[!grepl("_0$", msobject$rawData$MS1$peakID),] # keep all raw data
    }
  }
  ##############################################################################
  # msLevel 2
  # if acquired in DIA, MS2 is processed as msLevel 1
  if (acquisitionmode == "DIA"){
    if("MS2" %in% names(msobject$rawData)){
      cat("\n   Processing MS2...")
      for (cE in names(msobject$rawData$MS2)){
        cat("\n     partitioning...")
        msobject <- partitioning(msobject, dmzagglom = dmzagglom2, drtagglom = drtagglom2,
                                 minpeak = minpeak2, mslevel = "MS2", cE = cE)
        cat("OK")
        cat("\n     clustering...")
        msobject <- clustering(msobject, dmzagglom = dmzagglom2, drtclust = drtclust2,
                               minpeak = minpeak2, mslevel = "MS2", cE = cE)
        cat("OK")
        cat("\n     detecting peaks...")
        msobject <- peakdetection(msobject, minpeak = minpeak2, drtminpeak = drtminpeak2,
                                  drtmaxpeak = drtmaxpeak2, drtgap = drtgap2,
                                  recurs = recurs2, weight = weight2,
                                  sb = sb2, sn = sn2, minint = minint2,
                                  ended = 2, mslevel = "MS2", cE = cE)
        cat("OK")
      }
      msobject$peaklist$MS2 <- do.call(rbind, msobject$peaklist[["MS2"]])
      rownames(msobject$peaklist$MS2) <- msobject$peaklist$MS2$peakID
      msobject$rawData$MS2 <- do.call(rbind, msobject$rawData$MS2)
      msobject$rawData$MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peak", "Scan")]
      colnames(msobject$rawData$MS2) <- c("mz", "RT", "int", "peakID", "Scan")
      # msobject$rawData$MS2 <- msobject$rawData$MS2[!grepl("_0$", msobject$rawData$MS2$peakID),] # keep all raw data
    }
  } else if (acquisitionmode == "DDA"){
    ############################################################################
    # if acquired in DDA, scans from MS2 are extracted directly
    if ("MS2" %in% names(msobject$rawData)){
      cat("\n   Processing MS2...")
      for (cE in names(msobject$rawData$MS2)){
        msobject$rawData$MS2[[cE]]$peakID <- paste(paste("MS2_", cE, sep=""), 
                                                   msobject$rawData$MS2[[cE]]$Scan, sep="_")
        
        cat("OK")
      }
      msobject$rawData$MS2 <- do.call(rbind, msobject$rawData$MS2)
      msobject$rawData$MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID", "Scan")]
      colnames(msobject$rawData$MS2) <- c("mz", "RT", "int", "peakID", "Scan")
      msobject$rawData$MS2 <- msobject$rawData$MS2[msobject$rawData$MS2$int >= minint2,]
    }
  }
  #============================================================================#
  # Search for isotopes
  #============================================================================#
  cat("\n Searching for isotopes...")
  ##############################################################################
  # msLevel 1
  if("MS1" %in% names(msobject$rawData)){
    cat("\n   MS1...")
    peaklistIso <- annotateIsotopes(peaklist = msobject$peaklist$MS1,
                                    rawScans = msobject$rawData$MS1,
                                    dmz = dmzIso1, 
                                    drt = drtIso1,
                                    massdiff = 1.003355,
                                    charge = 1,
                                    isotopeAb = 0.01109,
                                    m0mass = 12,
                                    corThr = 0.8,
                                    checkInt = TRUE,
                                    checkCor = TRUE)
    msobject$peaklist$MS1 <- peaklistIso
    cat("OK")
  }
  ##############################################################################
  # msLevel 2
  if("MS2" %in% names(msobject$rawData)){
    if (acquisitionmode == "DIA"){
      cat("\n   MS2...")
      peaklistIso <- annotateIsotopes(peaklist = msobject$peaklist$MS2,
                                      rawScans = msobject$rawData$MS2,
                                      dmz = dmzIso2, 
                                      drt = drtIso2,
                                      massdiff = 1.003355,
                                      charge = 1,
                                      isotopeAb = 0.01109,
                                      m0mass = 12,
                                      corThr = 0.8,
                                      checkInt = TRUE,
                                      checkCor = TRUE)
      msobject$peaklist$MS2 <- peaklistIso
      cat("OK")
    }
  }
  cat("\n")
  
  if("MS1" %in% names(msobject$rawData)){
    msobject$processing$MS1$partIndex <- NULL
    msobject$processing$MS1$clustIndex <- NULL
    msobject$processing$MS1$peakIndex <- NULL
  }
  if("MS2" %in% names(msobject$rawData)){
    msobject$processing$MS2$partIndex <- NULL
    msobject$processing$MS2$clustIndex <- NULL
    msobject$processing$MS2$peakIndex <- NULL
  }
  msobject$annotation <- list()
    
  return(msobject)
}

# setmsbatch
#' Create msbatch for batch processing.
#'
#' Create msbatch from a list of msobjects to build an msbatch.
#'
#' @param msobjectlist list of msobjects.
#' @param metadata sample metadata. Optional. It can be a csv file or a data.frame 
#' with 3 columns (sample, acquistionmode and sampletype).
#' 
#' @return msbatch
#' 
#' @details samples are sorted following the metadata data.frame.
#' 
#' @seealso \link{dataProcessing} and \link{batchdataProcessing}
#'
#' @examples
#' \dontrun{
#' msbatch <- setmsbatch(msobjectlist)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@iislafe.es>
setmsbatch <- function(msobjectlist, 
                       metadata){
  
  #============================================================================#
  # Extract metadata from csv or build it fro msobjects list
  #============================================================================#
  if (missing(metadata)){
    sample <- unlist(lapply(msobjectlist, function(x) x$metaData$generalMetadata$file))
    sampletype <- rep("sample", length(msobjectlist))
    acquisitionmode <- unlist(lapply(msobjectlist, function(x) x$metaData$generalMetadata$acquisitionmode))
    metadata <- data.frame(sample, acquisitionmode, sampletype)
  } else {
    if (is.character(metadata)){
      metadata <- read.csv(metadata, header = TRUE)
    }
    if(is.data.frame(metadata)){
      if(!all(c("sample", "acquisitionmode", "sampletype") %in% colnames(metadata))){
        stop("metadata must have at least 3 columns sample, acquisitionmode and sampletype")
      }
      if (!(all(file.exists(metadata$sample)) | all(file.exists(paste(metadata$sample, ".mzXML", sep = ""))))){
        stop("samples in metadata must have the same name than their mzXML files")
      }
      if (!all(grepl("mzXML", metadata$sample))){
        metadata$sample <- paste(metadata$sample, ".mzXML", sep = "")
      }
      if (!all(metadata$acquisitionmode %in% c("ms", "dia", "dda", "MS", "DIA", "DDA"))){
        stop("acquisitionmode in metadata can be MS, DIA or DDA")
      } else {
        metadata$acquisitionmode <- toupper(metadata$acquisitionmode)
      }
    }
  }
  
  #============================================================================#
  # Check that all samples are in the msobject list and sort them 
  #============================================================================#
  msfiles <- unlist(lapply(msobjectlist, function(x) x$metaData$generalMetadata$file))
  inmsobjectlist <- unlist(sapply(metadata$sample, function(x) x %in% msfiles))
  metadata <- metadata[inmsobjectlist,]
  
  msobjectlistorder <- unlist(sapply(metadata$sample, function(x) match(x, msfiles)))
  if (length(msobjectlist) != length(msobjectlistorder)){
    stop("file names in metadata don't match msobjects file names")
  } else {
    msobjectlist <- msobjectlist[msobjectlistorder]
  }
  
  #============================================================================#
  # Create msbatch
  #============================================================================#
  msbatch <- list()
  msbatch$metaData <- metadata
  msbatch$msobjects <- msobjectlist
  msbatch$alignment <- list()
  msbatch$alignment$aligned <- FALSE
  msbatch$grouping <- list()
  msbatch$grouping$grouped <- FALSE
  msbatch$features <- data.frame()
  
  return(msbatch)
}

# batchdataProcessing
#' Process several mzXML files (peakpicking and isotope annotation) and create a 
#' msbatch for batch processing.
#'
#' Process several mzXML files (peakpicking and isotope annotation) and create a 
#' msbatch for batch processing.
#'
#' @param files  file paths of the mzXML files. Optional.
#' @param metadata csv file or data.frame with 3 columns: sample (samples named as the 
#' mzXML files), acquisitionmode (MS, DIA or DDA) and groups (i.e. blank, QC, sample).
#' DIA, DDA and MS files are allowed, but only DIA and DDA files will be used 
#' for lipid annotation.
#' @param polarity character value: negative or positive.
#' @param dmzagglom mz tolerance (in ppm) used for partitioning and clustering.
#' @param drtagglom rt window used for partitioning (in seconds).
#' @param drtclust rt window used for clustering (in seconds).
#' @param minpeak minimum number of measurements required for a peak.
#' @param drtgap maximum RT gap length to be filled (in seconds).
#' @param drtminpeak minimum RT width of a peak (in seconds). At least minpeak 
#' within the drtminpeak window are required to define a peak.
#' @param drtmaxpeak maximum RT width of a single peak (in seconds).
#' @param recurs maximum number of peaks within one EIC.
#' @param sb signal-to-base ratio.
#' @param sn signal-to-noise ratio.
#' @param minint minimum intensity of a peak.
#' @param weight weight for assigning measurements to a peak.
#' @param dmzIso mass tolerance for isotope matching.
#' @param drtIso time windows for isotope matching. 
#' @param parallel logical.
#' @param ncores number of cores to be used in case parallel is TRUE.
#'
#' @return msbatch
#'
#' @details This function executes 2 steps: 1) creates an msobject for each 
#' sample (using the \link{dataProcessing} function) and 2) sets an msbatch 
#' (\link{setmsbatch} function).
#' 
#' Numeric arguments accept one or two values for MS1 and MS2, respectively.
#' 
#' @seealso \link{dataProcessing} and \link{setmsbatch}
#'
#' @examples
#' \dontrun{
#' # if metadata is a data frame:
#' msbatch <- batchdataProcessing(metadata$sample, metadata, polarity = "positive",
#' dmzagglom = 25, drtagglom = 500, drtclust = 60, minpeak = c(5, 3),
#' drtgap = 5, drtminpeak = 20, drtmaxpeak = 100, recurs = 5, sb = c(3, 2),
#' sn = 2, minint = c(1000, 100), weight = 2, dmzIso = 10, drtIso = 5)
#' 
#' # if metadata is a csv file:
#' msbatch <- batchdataProcessing(metadata = "metadata.csv", polarity = "positive",
#' dmzagglom = 25, drtagglom = 500, drtclust = 60, minpeak = c(5, 3),
#' drtgap = 5, drtminpeak = 20, drtmaxpeak = 100, recurs = 5, sb = c(3, 2),
#' sn = 2, minint = c(1000, 100), weight = 2, dmzIso = 10, drtIso = 5)
#' }
#' 
#' @references Peak-picking algorithm has been imported from enviPick R-package:
#' https://cran.r-project.org/web/packages/enviPick/index.html
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@iislafe.es>
batchdataProcessing <- function(files, 
                                metadata, 
                                polarity,
                                dmzagglom = 15,
                                drtagglom = 500,
                                drtclust = 100,
                                minpeak = c(5, 3),
                                drtgap = 10,
                                drtminpeak = 15,
                                drtmaxpeak = c(100, 200),
                                recurs = 5,
                                sb = c(3, 2),
                                sn = 2,
                                minint = c(1000, 100),
                                weight = c(2, 3),
                                dmzIso = 10,
                                drtIso = 5, 
                                parallel = FALSE,
                                ncores){
  #============================================================================#
  # check arguments
  #============================================================================#
  if (missing(metadata)){
    stop("metadata must be a data.frame with at least 3 columns 
         (sample, acquisitionmode, sampletype) or a csv file")
  } else {
    if (is.character(metadata)){
      metadata <- read.csv(metadata, header = TRUE)
    }
    if(is.data.frame(metadata)){
      if(!all(c("sample", "acquisitionmode", "sampletype") %in% colnames(metadata))){
        stop("metadata must have at least 3 columns sample, acquisitionmode and sampletype")
      }
      if (!all(grepl("mzXML", metadata$sample))){
        metadata$sample[!grepl("mzXML", metadata$sample)] <- 
          paste(metadata$sample[!grepl("mzXML", metadata$sample)], ".mzXML", sep="")
      }
      if (!all(metadata$acquisitionmode %in% c("ms", "dia", "dda", "MS", "DIA", "DDA"))){
        stop("acquisitionmode in metadata can be MS, DIA or DDA")
      } else {
        metadata$acquisitionmode <- toupper(metadata$acquisitionmode)
      }
    }
  }
  if (missing(files)){
    files <- metadata$sample
  }
  
  if (!(all(files %in% metadata$sample))){
    message(paste("File",  files[!files %in% metadata$sample], 
                  "can't be found in metadata and will be omitted.\n"))
    files <- files[files %in% metadata$sample]
    if (length(files) < 2){
      stop("samples in metadata must have the same name than their mzXML files")
    }
  }
  if (!all(file.exists(metadata$sample))){
    message(paste("mzXML file for",  metadata$sample[!file.exists(metadata$sample)], 
                  "sample can't be found and will be omitted.\n"))
    metadata <- metadata[file.exists(metadata$sample),]
    files <- files[files %in% metadata$sample]
  }
  if (!all(file.exists(files))){
    stop(paste("File",  files[!file.exists(files)], "doesn't exist.\n"))
  }
  if (!all(metadata$sample %in% files)){
    metadata <- metadata[metadata$sample %in% files,]
  }
  if (!tolower(polarity) %in% c("positive", "negative")){
    stop("Polarity must be set to positive or negative")
  } else {
    polarity <- tolower(polarity)
  }
  if (parallel){
    if (missing(ncores)){
      stop("ncores argument is required if parallel is TRUE")
    }
    if (ncores > parallel::detectCores()){
      ncores <- parallel::detectCores() - 1
      message("ncores is greater than available cores. ", ncores, " will be used.")
    }
  }
  if (nrow(metadata) == 0){
    stop("No samples to be processed.")
  } else if (nrow(metadata) == 1){
    stop("For single mzXML files, use dataProcessing function instead.")
  }
  
  #============================================================================#
  # Process samples
  #============================================================================#
  if (parallel) {
    cl <- makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    `%d%` <- `%dopar%`
  } else {
    `%d%` <- `%do%`
  }
  f <- c()
  msobjects <- foreach::foreach(f = 1:nrow(metadata)) %d% {
    dataProcessing(metadata$sample[f],
                   metadata$acquisitionmode[f],
                   polarity = polarity,
                   dmzagglom = dmzagglom,
                   drtagglom = drtagglom,
                   drtclust = drtclust,
                   minpeak = minpeak,
                   drtgap = drtgap,
                   drtminpeak = drtminpeak,
                   drtmaxpeak = drtmaxpeak,
                   recurs = recurs,
                   sb = sb,
                   sn = sn,
                   minint = minint,
                   weight = weight,
                   dmzIso = dmzIso,
                   drtIso = drtIso)
  }
  if (parallel){
    parallel::stopCluster(cl)
  }
  msbatch <- setmsbatch(msobjects, metadata)
  
  return(msbatch)
}

# alignmsbatch
#' Align samples from an msbatch
#'
#' Align samples from an msbatch
#'
#' @param msbatch msbatch obtained from the \link{setmsbatch} function.
#' @param dmz mass tolerance between peak groups in ppm.
#' @param drt maximum rt distance between peaks for alignment in seconds.
#' @param minsamples minimum number of samples represented in each cluster 
#' used for the alignment.
#' @param minsamplesfrac minimum samples fraction represented in each cluster 
#' used for the alignment. Used to calculate minsamples in case it is missing.
#' @param span span parameter for loess rt deviation smoothing.
#' @param parallel logical. If TRUE, parallel processing will be performed.
#' @param ncores number of cores to be used in case parallel is TRUE.
#' 
#' @return aligned msbatch
#'
#' @examples
#' \dontrun{
#' msbatch <- setmsbatch(msbatch)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@iislafe.es>
alignmsbatch <- function(msbatch, 
                         dmz = 5, 
                         drt = 30, 
                         minsamples, 
                         minsamplesfrac = 0.75, 
                         span = 0.4, 
                         parallel = FALSE, 
                         ncores){
  #============================================================================#
  # Check arguments
  #============================================================================#
  ##############################################################################
  # check msbatch structure
  if (!is.list(msbatch) | !all(names(msbatch) %in% c("metaData", "msobjects", "alignment", "grouping", "features")) | 
      !is.data.frame(msbatch$metaData) | !is.list(msbatch$msobjects) | !is.list(msbatch$alignment) | 
      !is.list(msbatch$grouping) | !is.data.frame(msbatch$features)){
    stop("Wrong msbatch format")
  }
  ##############################################################################
  # check that all msobjects have an mslevel 1
  whichmslevel1 <- which(unlist(lapply(msbatch$msobjects, function(x) 
    1 %in% unique(x$metaData$scansMetadata$msLevel))))
  if (length(whichmslevel1) != nrow(msbatch$metaData)){
    warning("Removing samples with no MS1 level for alignment")
    msbatch$metaData <- msbatch$metaData[whichmslevel1,]
    msbatch$msobjects <- msbatch$msobjects[whichmslevel1]
  }
  ##############################################################################
  # set minsamples for alignment
  if (missing(minsamples)){
    minsamples <- floor(minsamplesfrac * length(msbatch$msobjects))
  }
  if (minsamples < 1){minsamples <- 1}
  ##############################################################################
  # Check parallel
  if (parallel){
    if (missing(ncores)){
      stop("ncores argument is required if parallel is TRUE")
    }
    if (ncores > parallel::detectCores()){
      ncores <- parallel::detectCores() - 1
      message("ncores is greater than availables cores. ", ncores, " will be used.")
    }
  }
  
  #============================================================================#
  # Extract peaks from all samples
  #============================================================================#
  peaks <- getallpeaks(msbatch)
  
  #============================================================================#
  # Create mz partitions based on dmz and drt
  #============================================================================#
  cat("\nCreating m.z partitions...")
  part <- .Call("agglom", as.numeric(peaks$mz),
                as.numeric(peaks$RT), as.integer(1),
                as.numeric(dmz), as.numeric(drt),
                PACKAGE = "LipidMS")
  peaks <- peaks[order(part, decreasing = FALSE),] 
  part <- part[order(part, decreasing = FALSE)]
  
  ##############################################################################
  # index partitions
  partIndex <- indexrtpart(peaks, part, minsamples)
  peaks$partID <- partIndex$idvector
  msbatch$alignment$partIndex <- partIndex$index
  cat("OK")
  
  #============================================================================#
  # Create RT clusters for each mz partition. Duplicate samples are not allowed 
  # in the same cluster
  #============================================================================#
  cat("\nClustering peaks by RT...")
  ##############################################################################
  # Clusterize (in parallel if required)
  if (parallel) {
    cl <- makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    `%d%` <- `%dopar%`
  } else {
    `%d%` <- `%do%`
  }
  clus <- foreach::foreach(p = 1:nrow(msbatch$alignment$partIndex)) %d% {
    start <- msbatch$alignment$partIndex[p, 1]
    end <- msbatch$alignment$partIndex[p, 2]
    measures <- peaks[start:end,]
    clusters <- clust(values = measures$RT, 
                      mins = measures$minRT, 
                      maxs = measures$maxRT, 
                      samples = measures$sample,
                      unique.samples = TRUE,
                      maxdist = drt,
                      ppm = FALSE)
    return(list(start = start, end = end, clusters = clusters))
  }
  if (parallel){
    parallel::stopCluster(cl)
  }
  
  ##############################################################################
  # Merge clust results
  startat <- 0
  clusts <- rep(0, nrow(peaks))
  roworder <- 1:nrow(peaks)
  for (p in 1:nrow(msbatch$alignment$partIndex)){
    start <- clus[[p]]$start
    end <- clus[[p]]$end
    clusters <- clus[[p]]$clusters + startat
    clusts[start:end] <- clusters
    roworder[start:end] <- roworder[start:end][order(clusters, decreasing = FALSE)]
    startat <- max(clusters)
  }
  peaks <- peaks[roworder,]
  clusts <- clusts[roworder]
  
  ##############################################################################
  # index clusters
  clustIndex <- indexrtpart(peaks, clusts, minsamples)
  peaks$clustID <- clustIndex$idvector
  msbatch$alignment$clustIndex <- clustIndex$index
  cat("OK")
  
  #============================================================================#
  # Create a RT matrix for each group (rows) and sample (columns)
  #============================================================================#
  cat("\nEstimating RT deviation...")
  ##############################################################################
  # rtgroupsMatrix
  rtgroupsMatrix <- matrix(nrow = nrow(msbatch$alignment$clustIndex), 
                           ncol = length(msbatch$msobjects))
  for (c in 1:nrow(msbatch$alignment$clustIndex)){
    start <- msbatch$alignment$clustIndex[c,1]
    end <- msbatch$alignment$clustIndex[c,2]
    cluster <- peaks[start:end,]
    for (s in as.numeric(unique(cluster$sample))){
      rtgroupsMatrix[c, s] <- as.numeric(cluster$RT[cluster$sample == s])
    }
  }
  
  ##############################################################################
  # sort by median rt
  rtmedian <- apply(rtgroupsMatrix, 1, median, na.rm = TRUE)
  rtgroupsMatrix <- rtgroupsMatrix[order(rtmedian, decreasing = FALSE),]
  rtmedian <- sort(rtmedian, decreasing = FALSE)
  
  ##############################################################################
  # rtdevMatrix: differences between median RT and individual samples RT
  rtdevMatrix <- rtgroupsMatrix - rtmedian
  cat("OK")
  
  #============================================================================#
  # RT correction
  #============================================================================#
  cat("\nAligning samples...")
  rtdevcorrected <- list()
  rtmodels <- list()
  for (i in 1:length(msbatch$msobjects)){
    ############################################################################
    # Adjust rt dev with loess
    rtlo <- loess(rtdevMatrix[,i] ~ rtmedian, span = span, degree = 1, 
                  family = "gaussian")
    rtmodels[[i]] <- rtlo
    
    ############################################################################
    # correct raw scans RT (in metaData)
    rt <- msbatch$msobjects[[i]]$metaData$scansMetadata$RT
    msbatch$msobjects[[i]]$metaData$scansMetadata$RT <- rtcorrection(rt, rtlo)
    rtdevcorrected[[i]] <- list(RT = rt, 
                                RTdev = rt - msbatch$msobjects[[i]]$metaData$scansMetadata$RT) 
    
    ############################################################################
    # correct raw scans RT (in MS1 and MS2)
    mslevels <- c("MS1", "MS2")[which(c("MS1", "MS2") %in% names(msbatch$msobjects[[i]]$rawData))]
    for (mslevel in mslevels){
      rt <- msbatch$msobjects[[i]]$rawData[[mslevel]]$RT
      msbatch$msobjects[[i]]$rawData[[mslevel]]$RT <- rtcorrection(rt, rtlo)
    }
    
    ############################################################################
    # correct peaklist RT (in MS1 and MS2)
    mslevels <- c("MS1", "MS2")[which(c("MS1", "MS2") %in% names(msbatch$msobjects[[i]]$peaklist))]
    for (mslevel in mslevels){
      rt <- msbatch$msobjects[[i]]$peaklist[[mslevel]]$RT
      msbatch$msobjects[[i]]$peaklist[[mslevel]]$RT <- rtcorrection(rt, rtlo)
    }
  }
  cat("OK\n")
  
  #============================================================================#
  # Save results in msbatch
  #============================================================================
  msbatch$alignment$aligned <- TRUE
  msbatch$alignment$peaks <- peaks
  msbatch$alignment$parameters$dmz <- dmz
  msbatch$alignment$parameters$drt <- drt
  msbatch$alignment$parameters$minsamplesfrac <- minsamplesfrac
  msbatch$alignment$parameters$span <- span
  msbatch$alignment$partIndex <- NULL
  msbatch$alignment$clustIndex <- NULL
  msbatch$alignment$rtdevMatrix <- NULL
  # msbatch$alignment$rtmodels <- rtmodels
  msbatch$alignment$rtdevcorrected <- rtdevcorrected
  msbatch$grouping <- list()
  msbatch$features <- data.frame()
  
  return(msbatch)
}

# groupmsbatch
#' Group features from an msbatch
#'
#' Group features from an msbatch
#'
#' @param msbatch msbatch obtained from \link{setmsbatch} or \link{alignmsbatch} 
#' functions.
#' @param dmz mass tolerance between peak groups for grouping in ppm.
#' @param drtagglom rt window for mz partitioning.
#' @param drt rt window for peaks clustering.
#' @param minsamples minimum number of samples represented in clusters 
#' used for grouping.
#' @param minsamplesfrac minimum samples fraction represented in each cluster 
#' used for grouping. Used to calculate minsamples in case it is missing.
#' @param parallel logical. If TRUE, parallel processing is performed.
#' @param ncores number of cores to be used in case parallel is TRUE.
#' 
#' @return grouped msbatch
#'
#' @examples
#' \dontrun{
#' msbatch <- groupmsbatch(msbatch)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@iislafe.es>
groupmsbatch <- function(msbatch, 
                         dmz = 5, 
                         drtagglom = 30, 
                         drt = 15, 
                         minsamples, 
                         minsamplesfrac = 0.25,
                         parallel = FALSE,
                         ncores){
  
  #============================================================================#
  # Check arguments
  #============================================================================#
  ##############################################################################
  # check msbatch structure
  if (!is.list(msbatch) | !all(names(msbatch) %in% c("metaData", "msobjects", "alignment", "grouping", "features")) | 
      !is.data.frame(msbatch$metaData) | !is.list(msbatch$msobjects) | !is.list(msbatch$alignment) | 
      !is.list(msbatch$grouping) | !is.data.frame(msbatch$features)){
    stop("Wrong msbatch format")
  }
  ##############################################################################
  # check that all msobjects have an mslevel 1
  whichmslevel1 <- which(unlist(lapply(msbatch$msobjects, function(x) 
    1 %in% unique(x$metaData$scansMetadata$msLevel))))
  if (length(whichmslevel1) != nrow(msbatch$metaData)){
    warning("Removing samples with no MS1 level for alignment")
    msbatch$metaData <- msbatch$metaData[whichmslevel1,]
    msbatch$msobjects <- msbatch$msobjects[whichmslevel1]
  }
  ##############################################################################
  # set minsamples for alignment
  if (missing(minsamples)){
    minsamples <- floor(minsamplesfrac * length(msbatch$msobjects))
  }
  if (minsamples < 1){minsamples <- 1}
  ##############################################################################
  # Check parallel
  if (parallel){
    if (missing(ncores)){
      stop("ncores argument is required if parallel is TRUE")
    }
    if (ncores > parallel::detectCores()){
      ncores <- parallel::detectCores() - 1
      message("ncores is greater than available cores. ", ncores, " will be used.")
    }
  }
  
  #============================================================================#
  # Extract peaks from all samples
  #============================================================================#
  peaks <- getallpeaks(msbatch)
  
  #============================================================================#
  # Create mz partitions based on dmz and drt
  #============================================================================#
  cat("\nCreating m.z partitions...")
  part <- .Call("agglom", as.numeric(peaks$mz),
                as.numeric(peaks$RT), as.integer(1),
                as.numeric(dmz), as.numeric(drtagglom),
                PACKAGE = "LipidMS")
  peaks <- peaks[order(part, decreasing = FALSE),] 
  part <- part[order(part, decreasing = FALSE)]
  
  ##############################################################################
  # index partitions
  partIndex <- indexrtpart(peaks, part, minsamples)
  peaks$partID <- partIndex$idvector
  msbatch$grouping$partIndex <- partIndex$index
  cat("OK")
  
  #============================================================================#
  # Create mz clusters for each partition.
  #============================================================================#
  cat("\nClustering peaks by m.z...")
  ##############################################################################
  # Clusterize (in parallel if required)
  if (parallel) {
    cl <- makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    `%d%` <- `%dopar%`
  } else {
    `%d%` <- `%do%`
  }
  clus <- foreach::foreach(p = 1:nrow(msbatch$grouping$partIndex)) %d% {
    start <- msbatch$grouping$partIndex[p, 1]
    end <- msbatch$grouping$partIndex[p, 2]
    measures <- peaks[start:end,]
    clusters <- clust(values = measures$mz,
                      mins = measures$mz,
                      maxs = measures$mz,
                      samples = measures$sample,
                      unique.samples = FALSE,
                      maxdist = dmz,
                      ppm = TRUE)
    return(list(start = start, end = end, clusters = clusters))
  }
  if (parallel){
    parallel::stopCluster(cl)
  }
  
  ##############################################################################
  # merge clust results
  startat <- 0
  clusts <- rep(0, nrow(peaks))
  roworder <- 1:nrow(peaks)
  for (p in 1:nrow(msbatch$grouping$partIndex)){
    start <- clus[[p]]$start
    end <- clus[[p]]$end
    clusters <- clus[[p]]$clusters + startat
    clusts[start:end] <- clusters
    roworder[start:end] <- roworder[start:end][order(clusters, decreasing = FALSE)]
    startat <- max(clusters)
  }
  peaks <- peaks[roworder,]
  clusts <- clusts[roworder]
  
  ##############################################################################
  # index clusters
  clustIndex <- indexrtpart(peaks, clusts, minsamples)
  peaks$clustID <- clustIndex$idvector
  msbatch$grouping$clustIndex <- clustIndex$index
  cat("OK")
  
  #============================================================================#
  # Create RT clusters for each mz cluster
  #============================================================================#
  cat("\nGrouping peaks by RT...")
  ##############################################################################
  # Clusterize (in parallel if required)
  if (parallel) {
    cl <- makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    `%d%` <- `%dopar%`
  } else {
    `%d%` <- `%do%`
  }
  gr <- foreach::foreach(c = 1:nrow(msbatch$grouping$clustIndex)) %d% {
    start <- msbatch$grouping$clustIndex[c, 1]
    end <- msbatch$grouping$clustIndex[c, 2]
    measures <- peaks[start:end,]
    groups <- clust(values = measures$RT,
                    mins = measures$RT,
                    maxs = measures$RT,
                    samples = measures$sample,
                    unique.samples = TRUE,
                    maxdist = drt,
                    ppm = FALSE)
    return(list(start = start, end = end, groups = groups))
  }
  if (parallel){
    parallel::stopCluster(cl)
  }
  
  ##############################################################################
  # Merge clust results
  startat <- 0
  group <- rep(0, nrow(peaks))
  roworder <- 1:nrow(peaks)
  for (c in 1:nrow(msbatch$grouping$clustIndex)){
    start <- gr[[c]]$start
    end <- gr[[c]]$end
    groups <- gr[[c]]$groups + startat
    group[start:end] <- groups
    roworder[start:end] <- roworder[start:end][order(groups, decreasing = FALSE)]
    startat <- max(groups)
  }
  peaks <- peaks[roworder,]
  group <- group[roworder]
  
  ##############################################################################
  # index clusters
  groupIndex <- indexrtpart(peaks, group, minsamples)
  peaks$groupID <- groupIndex$idvector
  msbatch$grouping$groupIndex <- groupIndex$index
  cat("OK")
  
  #============================================================================#
  # Save grouping results in msbatch
  #============================================================================#
  msbatch$grouping$grouped <- TRUE
  msbatch$grouping$peaks <- peaks
  msbatch$grouping$parameters$dmz <- dmz
  msbatch$grouping$parameters$drtagglom <- drtagglom
  msbatch$grouping$parameters$drt <- drt
  msbatch$grouping$parameters$minsamplesfrac <- minsamplesfrac
  msbatch$grouping$parameters$minsamples <- minsamples
  
  #============================================================================#
  # Create feature table
  #============================================================================#
  cat("\nBuilding data matrix...")
  msbatch <- getfeaturestable(msbatch)
  cat("OK\n")
  
  #============================================================================#
  # Remove unnecessary data
  #============================================================================#
  msbatch$grouping$partIndex <- NULL
  msbatch$grouping$clustIndex <- NULL
  msbatch$grouping$groupIndex <- NULL
  
  return(msbatch)
}

# fillpeaksmsbatch
#' Fill peaks from a grouped msbatch
#'
#' Use grouping results to target all peaks from the msbatch in each sample and 
#' refill intensities at the features table.
#'
#' @param msbatch msbatch obtained from the \link{groupmsbatch} function.
#' 
#' @return msbatch
#'
#' @examples
#' \dontrun{
#' msbatch <- fillpeaksmsbatch(msbatch)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@iislafe.es>
fillpeaksmsbatch <- function(msbatch){
  #============================================================================#
  # Check arguments
  #============================================================================#
  ##############################################################################
  # check msbatch structure
  if (!msbatch$grouping$grouped){
    stop("msbatch needs to be grouped before filling peaks. Use groupmsbatch function.")
  }
  if (!is.list(msbatch) | !all(names(msbatch) %in% c("metaData", "msobjects", "alignment", "grouping", "features")) | 
      !is.data.frame(msbatch$metaData) | !is.list(msbatch$msobjects) | !is.list(msbatch$alignment) | 
      !is.list(msbatch$grouping) | !is.data.frame(msbatch$features)){
    stop("Wrong msbatch format")
  }
  ##############################################################################
  # check that all msobjects have an mslevel 1
  whichmslevel1 <- which(unlist(lapply(msbatch$msobjects, function(x) 
    1 %in% unique(x$metaData$scansMetadata$msLevel))))
  if (length(whichmslevel1) != nrow(msbatch$metaData)){
    warning("Removing samples with no MS1 level for alignment")
    msbatch$metaData <- msbatch$metaData[whichmslevel1,]
    msbatch$msobjects <- msbatch$msobjects[whichmslevel1]
  }
  
  #============================================================================#
  # Extract features table and additional info required for filling peaks
  #============================================================================#
  features <- msbatch$features[,!colnames(msbatch$features) %in% make.names(msbatch$metaData$sample)]
  fmatrix <- msbatch$features[,colnames(msbatch$features) %in% make.names(msbatch$metaData$sample)]
  # peaks <- msbatch$grouping$peaks
  # 
  dmz <- msbatch$grouping$parameters$dmz
  
  #============================================================================#
  # Fill peaks
  #============================================================================#
  ##############################################################################
  # for each sample
  for (s in 1:ncol(fmatrix)){
    ############################################################################
    # Make an index for raw MS data to speed up the function
    MS1 <- msbatch$msobjects[[s]]$rawData$MS1
    drt <- msbatch$msobjects[[s]]$metaData$generalMetadata$endTime - 
      msbatch$msobjects[[s]]$metaData$generalMetadata$startTime
    mz <- MS1$mz
    ord <- order(mz)
    mz <- mz[ord]
    rt <- MS1$RT[ord]
    int <- MS1$int[ord]
    part <- .Call("agglom", mz, rt, as.integer(1), dmz*2, drt, PACKAGE = "LipidMS")
    part <- part[order(part, decreasing = FALSE)]
    maxit <- max(part)
    index <- .Call("indexed", as.integer(part), int, 0, max(int),
                   as.integer(maxit), PACKAGE = "LipidMS")
    index <- cbind(index, mz[index[,1]])
    inimz <- index[,4]
    
    ############################################################################
    # for each NA 
    for(m in 1:nrow(fmatrix)){
      ##########################################################################
      # get mz and rt limits
      minmz <- features$mz[m] - features$mz[m]*dmz/1e6
      maxmz <- features$mz[m] + features$mz[m]*dmz/1e6
      minRT <- features$iniRT[m]
      maxRT <- features$endRT[m]
      ##########################################################################
      # Search into the rawData and fill up fmatrix with recalculated intensities
      p <- which(inimz < maxmz)
      if (length(p) > 0){
        p <- p[length(p)]
        start <- index[p, 1]
        end <- index[p, 2]
        subset <- (start:end)[mz[start:end] >= minmz & mz[start:end] <= maxmz & 
                                rt[start:end] >= minRT & rt[start:end] <= maxRT]
        if (length(subset) > 0){
          fmatrix[m,s] <- sum(int[subset], na.rm = TRUE)
        } else {
          fmatrix[m,s] <- 0
        }
      } else {
        fmatrix[m,s] <- 0
      }
    }
  }
  
  #============================================================================#
  # Save results in msbatch
  #============================================================================#
  msbatch$features <- cbind(features, fmatrix)
  
  return(msbatch)
}



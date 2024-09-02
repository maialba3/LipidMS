# readMSfile
#' Read mzXML file and initiate msobject
#'
#' Read mzXML file and initiate msobject
#'
#' @param file file path for a mzXML file
#' @param polarity character value: negative or positive.
#'
#' @return msobject
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
readMSfile <- function(file, polarity){
  # 1. read data with readMzXmlFile
  ms <- readMzXmlData::readMzXmlFile(file.path(file))

  # 2. Extract metaData (general and scan by scan)
  if (is.null(ms[[1]]$metaData$startTime)){
    startTime <- min(unlist(lapply(ms, function(x) x$metaData$retentionTime)))
    endTime <- max(unlist(lapply(ms, function(x) x$metaData$retentionTime)))
  } else {
    startTime <- ms[[1]]$metaData$startTime
    endTime <- ms[[1]]$metaData$endTime
  }
  
  collisionEnergies <- sort(unique(unlist(lapply(ms, function(x) if(!is.null(x$metaData$collisionEnergy)){x$metaData$collisionEnergy} else {0}))))
  
  generalMetadata <- list(file = file, scans = length(ms),
                          startTime = startTime,
                          endTime = endTime,
                          collisionEnergies = collisionEnergies)
  
  scansMetadata <-
    data.frame(msLevel = unlist(lapply(ms, function(x) if(!is.null(x$metaData$msLevel)){x$metaData$msLevel} else {NA})),
               polarity = unlist(lapply(ms, function(x) if(!is.null(x$metaData$polarity)){x$metaData$polarity} else {NA})),
               scanType = unlist(lapply(ms, function(x) if(!is.null(x$metaData$scanType)){x$metaData$scanType} else {NA})),
               centroided = unlist(lapply(ms, function(x) if(!is.null(x$metaData$centroided)){x$metaData$centroided} else {NA})),
               RT = unlist(lapply(ms, function(x) if(!is.null(x$metaData$retentionTime)){x$metaData$retentionTime} else {NA})),
               peaksCount = unlist(lapply(ms, function(x) if(!is.null(x$metaData$peaksCount)){x$metaData$peaksCount} else {NA})),
               lowMz = unlist(lapply(ms, function(x) if(!is.null(x$metaData$lowMz)){x$metaData$lowMz} else {NA})),
               highMz = unlist(lapply(ms, function(x) if(!is.null(x$metaData$highMz)){x$metaData$highMz} else {NA})),
               basePeakMz = unlist(lapply(ms, function(x) if(!is.null(x$metaData$basePeakMz)){x$metaData$basePeakMz} else {NA})),
               basePeakInt = unlist(lapply(ms, function(x) if(!is.null(x$metaData$basePeakInt)){x$metaData$basePeakInt} else {NA})),
               totIonCurrent = unlist(lapply(ms, function(x) if(!is.null(x$metaData$totIonCurrent)){x$metaData$totIonCurrent} else {NA})),
               precursor = unlist(lapply(ms, function(x) if(!is.null(x$metaData$precursorMz)){x$metaData$precursorMz} else {NA})),
               collisionEnergy = unlist(lapply(ms, function(x) if(!is.null(x$metaData$collisionEnergy)){x$metaData$collisionEnergy} else {0})),
               stringsAsFactors = FALSE)
  scanOrder <- rep(0,nrow(scansMetadata))
  for (l in unique(scansMetadata$msLevel)){
    for (c in unique(scansMetadata$collisionEnergy))
      scanOrder[scansMetadata$msLevel == l & scansMetadata$collisionEnergy == c] <-
        as.numeric(factor(scansMetadata$RT[scansMetadata$msLevel == l & scansMetadata$collisionEnergy == c]))
  }
  scansMetadata$Scan <- scanOrder

  # 3. Extract scans
  mz <- unlist(lapply(ms, function(x) x$spectrum$mass))
  int <- unlist(lapply(ms, function(x) x$spectrum$intensity))
  RT <- unlist(mapply(rep, scansMetadata$RT, scansMetadata$peaksCount))
  mslevel <- unlist(mapply(rep, scansMetadata$msLevel, scansMetadata$peaksCount))
  collisionEnergy <- unlist(mapply(rep, scansMetadata$collisionEnergy, scansMetadata$peaksCount))
  scannum <- unlist(mapply(rep, scansMetadata$Scan, scansMetadata$peaksCount))
  scans <- data.frame(mz = mz, int = int, RT = RT, mslevel = mslevel,
                      collisionEnergy = collisionEnergy,
                      part = 0, clust = 0, peak = 0, Scan = scannum)
  
  if (polarity == "positive") {pol <- "+"} else {pol <- "-"}
  keepPolarity <- scansMetadata$Scan[which(scansMetadata$polarity == pol)]
  scansMetadata <- scansMetadata[scansMetadata$Scan %in% keepPolarity,]
  scans <- scans[scans$Scan %in% keepPolarity,]
  generalMetadata$polarity <- polarity
  
  if (nrow(scans) == 0){
    stop(paste(c("No scans were found for ESI", pol), collapse = ""))
  }
  
  # 4. Generate msobject
  msobject <- list()
  msobject$metaData <- list(generalMetadata = generalMetadata,
                            scansMetadata = scansMetadata)
  msobject$processing <- list()
  
  if (all(c(1, 2) %in% unique(scansMetadata$msLevel))){
    MS1 <-  scans[scans$mslevel == 1,]
    MS1 <- split(MS1, MS1$collisionEnergy)
    # Add MS1 to msobject
    msobject$rawData$MS1 <- MS1
    
    MS2 <- scans[scans$mslevel == 2,]
    MS2 <- split(MS2, MS2$collisionEnergy)
    # Add MS2 to msobject
    msobject$rawData$MS2 <- MS2
  } else { # exception for mzXML files from equipments without explicit DIA mode (such as Agilent)
    if (any(scansMetadata$collisionEnergy == 0)){
      MS1 <-  scans[scans$collisionEnergy == 0,]
      MS1 <- split(MS1, MS1$collisionEnergy)
      # Add MS1 to msobject
      msobject$rawData$MS1 <- MS1
    }
    if (any(scansMetadata$collisionEnergy > 0)){
      MS2 <- scans[scans$collisionEnergy > 0,]
      MS2 <- split(MS2, MS2$collisionEnergy)
      scansMetadata$msLevel[scansMetadata$collisionEnergy > 0] <- 2
      # Add MS2 to msobject
      msobject$rawData$MS2 <- MS2
    }
  }
  return(msobject)
}

# partitioning
#' agglomarative partitioning for LC-HRMS data based on enviPick algorithm
#'
#' agglomarative partitioning for LC-HRMS data based on enviPick algorithm
#'
#' @param msobject msobject generated by \link{readMSfile}
#' @param dmzagglom mz tolerance for partitions
#' @param drtagglom RT window for partitions
#' @param minpeak minimum number of measures to define a peak
#' @param mslevel MS level information to access msobject
#' @param cE collision energy information to access msobject
#'
#' @return msobject
#'
#' @keywords internal
#' 
#' @references Peak-picking algorithm has been imported from enviPick R-package:
#' https://cran.r-project.org/web/packages/enviPick/index.html
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
partitioning <- function(msobject,
                         dmzagglom,
                         drtagglom,
                         minpeak,
                         mslevel,
                         cE){
  # save parameters
  msobject$processing[[mslevel]]$parameters$maxint <- max(msobject$rawData[[mslevel]][[cE]]$int)
  msobject$processing[[mslevel]]$parameters$dmzagglom <- dmzagglom
  msobject$processing[[mslevel]]$parameters$drtagglom <- drtagglom
  msobject$processing[[mslevel]]$parameters$minpeak <- minpeak

  # order ms measures by increasing mz
  msobject$rawData[[mslevel]][[cE]] <- msobject$rawData[[mslevel]][[cE]][order(msobject$rawData[[mslevel]][[cE]]$mz,
                                                   decreasing = FALSE),]

  # Agglomerative partitioning: agglom function from enviPick package
  part <- .Call("agglom", 
                as.numeric(msobject$rawData[[mslevel]][[cE]]$mz),
                as.numeric(msobject$rawData[[mslevel]][[cE]]$RT), 
                as.integer(1),
                as.numeric(dmzagglom), 
                as.numeric(drtagglom),
                PACKAGE = "LipidMS")

  # Index of partitions: indexed function from enviPick package
    # order ms measures by partition order
  msobject$rawData[[mslevel]][[cE]] <- msobject$rawData[[mslevel]][[cE]][order(part,
                                                               decreasing = FALSE),]
  part <- part[order(part, decreasing = FALSE)]
  maxit <- max(part)
  index <- .Call("indexed", 
                 as.integer(part), 
                 as.numeric(msobject$rawData[[mslevel]][[cE]]$int),
                 as.integer(minpeak),
                 as.numeric(msobject$processing[[mslevel]]$parameters$maxint),
                 as.integer(maxit),
                 PACKAGE = "LipidMS")
  index <- index[index[,2] != 0,,drop = FALSE]
  colnames(index) <- c("start", "end", "length")

  # Assign partition ID: partID function from enviPick
  partID <- .Call("partID", 
                  as.integer(index),
                  as.integer(nrow(msobject$rawData[[mslevel]][[cE]])),
                  PACKAGE = "LipidMS")

  # save partitions
  msobject$rawData[[mslevel]][[cE]]$part <- partID
  msobject$processing[[mslevel]]$partIndex[[cE]] <- index

  return(msobject)
}

# clustering
#' EIC extraction based on previous partitions generated by \link{partitioning}
#'
#' EIC extraction based on previous partitions generated by \link{partitioning}
#'
#' @param msobject msobject generated by \link{partitioning}
#' @param dmzagglom mz tolerance for clusters
#' @param drtclust RT window for clusters
#' @param minpeak minimum number of measures to define a peak
#' @param mslevel info to access msobject
#' @param cE info to access msobject
#'
#' @return msobject
#'
#' @keywords internal
#' 
#' @references Peak-picking algorithm has been imported from enviPick R-package:
#' https://cran.r-project.org/web/packages/enviPick/index.html
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
clustering <- function(msobject,
                       dmzagglom,
                       drtclust,
                       minpeak,
                       mslevel,
                       cE){

  # save parameters
  msobject$processing[[mslevel]]$parameters$drtclust <- drtclust

  startat <- 0
  roworder <- 1:nrow(msobject$rawData[[mslevel]][[cE]])
  for (k in 1:nrow(msobject$processing[[mslevel]]$partIndex[[cE]])) {
    # get EICs using getEIC function from enviPick
    start <- msobject$processing[[mslevel]]$partIndex[[cE]][k,1]
    end <- msobject$processing[[mslevel]]$partIndex[[cE]][k,2]
    if ((end - (start+1)) > 1){
      clusters <- .Call("getEIC",
                        as.numeric(msobject$rawData[[mslevel]][[cE]]$mz[start:end]),
                        as.numeric(msobject$rawData[[mslevel]][[cE]]$RT[start:end]),
                        as.numeric(msobject$rawData[[mslevel]][[cE]]$int[start:end]),
                        as.integer(order(msobject$rawData[[mslevel]][[cE]]$int[start:end], decreasing = TRUE)),
                        as.integer(order(msobject$rawData[[mslevel]][[cE]]$RT[start:end], decreasing = FALSE)),
                        as.numeric(dmzagglom), 
                        as.integer(1), 
                        as.numeric(drtclust),
                        as.integer(1), 
                        PACKAGE = "LipidMS")
      clust <- clusters[, 10] + startat
      msobject$rawData[[mslevel]][[cE]]$clust[start:end] <- clust
      roworder[start:end] <- roworder[start:end][order(clust, decreasing = FALSE)]
      startat <- max(clust)
    } else {
      msobject$rawData[[mslevel]][[cE]]$clust[start:end] <- startat
      startat <- startat + 1
    }
  }
  msobject$rawData[[mslevel]][[cE]] <- msobject$rawData[[mslevel]][[cE]][roworder,]

  # Index of clusters: indexed function from enviPick
  maxit <- max(msobject$rawData[[mslevel]][[cE]]$clust)
  index <- .Call("indexed",
                 as.integer(msobject$rawData[[mslevel]][[cE]]$clust),
                 as.numeric(msobject$rawData[[mslevel]][[cE]]$int),
                 as.integer(minpeak),
                 as.numeric(msobject$processing[[mslevel]]$parameters$maxint),
                 as.integer(maxit),
                 PACKAGE = "LipidMS")
  index <- index[index[,2] != 0,,drop = FALSE]
  colnames(index) <- c("start", "end", "length")

  # Assign cluster ID: partID function from enviPick
  clustID <- .Call("partID", 
                   as.integer(index),
                   as.integer(nrow(msobject$rawData[[mslevel]][[cE]])),
                   PACKAGE = "LipidMS")

  # save clusters
  msobject$rawData[[mslevel]][[cE]]$clust <- clustID
  msobject$processing[[mslevel]]$clustIndex[[cE]] <- index

  return(msobject)
}

# peakdetection
#' peak-pick based on previous EIC clusters generated by \link{clustering}
#'
#' peak-pick based on previous EIC clusters generated by \link{clustering}
#'
#' @param msobject msobject generated by \link{clustering}
#' @param minpeak minimum number of measures to define a peak
#' @param drtminpeak minimum RT length of a peak
#' @param drtmaxpeak maximum RT length of a peak
#' @param drtgap maximum RT gap to be filled
#' @param recurs maximum number of peaks for a EIC
#' @param weight weight for assigning measurements to a peak
#' @param ended number of failures allowed when detecting peaks
#' @param sb signal-to-base ration
#' @param sn signal-to-noise ratio
#' @param minint minimum intensity
#' @param mslevel info to access msobject
#' @param cE info to access msobject
#'
#' @return msobject
#'
#' @keywords internal
#' 
#' @references Peak-picking algorithm has been imported from enviPick R-package:
#' https://cran.r-project.org/web/packages/enviPick/index.html
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
peakdetection <- function(msobject,
                          minpeak,
                          drtminpeak,
                          drtmaxpeak,
                          drtgap,
                          recurs,
                          weight,
                          ended,
                          sb,
                          sn,
                          minint,
                          mslevel,
                          cE){

  # save parameters
  msobject$processing[[mslevel]]$parameters$drtminpeak <- drtminpeak
  msobject$processing[[mslevel]]$parameters$drtmaxpeak <- drtmaxpeak
  msobject$processing[[mslevel]]$parameters$drtgap <- drtgap
  msobject$processing[[mslevel]]$parameters$recurs <- recurs
  msobject$processing[[mslevel]]$parameters$weight <- weight
  msobject$processing[[mslevel]]$parameters$ended <- ended
  msobject$processing[[mslevel]]$parameters$sb <- sb
  msobject$processing[[mslevel]]$parameters$sn <- sn
  msobject$processing[[mslevel]]$parameters$minint <- minint

  msobject$rawData[[mslevel]][[cE]]$id <- 1:nrow(msobject$rawData[[mslevel]][[cE]])
  level <- as.numeric(gsub("MS", "", mslevel))

  startat <- 0
  npeaks <- 0
  areas <- c()
  roworder <- 1:nrow(msobject$rawData[[mslevel]][[cE]])
  for (k in 1:nrow(msobject$processing[[mslevel]]$clustIndex[[cE]])) {
    if (msobject$processing[[mslevel]]$clustIndex[[cE]][k, 3] >= minpeak){
      start <- msobject$processing[[mslevel]]$clustIndex[[cE]][k,1]
      end <- msobject$processing[[mslevel]]$clustIndex[[cE]][k,2]
      # Fill RT gaps < drtgap: gapfill function from enviPick
      out1 <- .Call("gapfill",
                    as.numeric(msobject$rawData[[mslevel]][[cE]]$RT[start:end]),
                    as.numeric(msobject$rawData[[mslevel]][[cE]]$int[start:end]),
                    as.integer(order(msobject$rawData[[mslevel]][[cE]]$RT[start:end], decreasing = FALSE)),
                    as.numeric(msobject$rawData[[mslevel]][[cE]]$mz[start:end]),
                    as.numeric(msobject$rawData[[mslevel]][[cE]]$id[start:end]),
                    as.numeric(msobject$metaData$scansMetadata$RT[msobject$metaData$scansMetadata$msLevel == level & msobject$metaData$scansMetadata$collisionEnergy == cE]),
                    as.numeric(drtgap),
                    PACKAGE = "LipidMS")
      out1 <- matrix(out1,ncol=10)
      colnames(out1)<-c("mz","intens","RT","index","intens_filt","1pick","pickcrit","baseline","intens_corr","2pick")
      # Filter step: yet to be implemented
      out1[, 5] <- out1[, 2]
      # Peak detection, baseline calculation and 2nd peak detection
      out2 <- .Call("pickpeak",
                    as.numeric(out1),
                    as.numeric(drtminpeak),
                    as.numeric(drtmaxpeak),
                    as.integer(minpeak),
                    as.integer(recurs),
                    as.numeric(weight),
                    as.numeric(sb),
                    as.numeric(sn), 
                    as.numeric(minint),
                    as.numeric(msobject$processing[[mslevel]]$parameters$maxint),
                    as.integer(ended),
                    as.integer(2),
                    PACKAGE = "LipidMS")
      out2 <- matrix(out2, ncol = 10)
      colnames(out2) <- c("m/z", "intens", "RT", "index", "intens_filt", "1pick", 
                          "pickcrit", "baseline", "intens_corr", "2pick");
      if(!all(out2[,10] == 0)){
        npeaks <- npeaks + length(unique(out2[,10]))
        out2[,10] <- out2[,10] + startat
        out2 <- out2[out2[,10] != startat,]
        area <- tapply(out2[,8], out2[,10], sum) # sum int of filled peak
        areas <- c(areas, area)
        peak <- as.numeric(sapply(msobject$rawData[[mslevel]][[cE]]$id[start:end], 
                                  function(x) if(x %in% out2[,4]){out2[out2[,4] == x,10]} else {0}))
        msobject$rawData[[mslevel]][[cE]]$peak[start:end] <- peak
        roworder[start:end] <- roworder[start:end][order(peak, decreasing = FALSE)]
        startat <- c(max(out2[,10]))
      }
    }
  }
  msobject$rawData[[mslevel]][[cE]] <- msobject$rawData[[mslevel]][[cE]][roworder,]

  # assign peakID
  # Index of peaks: indexed function from enviPick
  maxit <- max(msobject$rawData[[mslevel]][[cE]]$peak)
  if(maxit > 0){
    index <- .Call("indexed",
                   as.integer(msobject$rawData[[mslevel]][[cE]]$peak),
                   as.numeric(msobject$rawData[[mslevel]][[cE]]$int),
                   as.integer(minpeak),
                   as.numeric(msobject$processing[[mslevel]]$parameters$maxint),
                   as.integer(maxit),
                   PACKAGE="LipidMS")
    if(any(index[,2]!=0)){
      index <- index[index[,2] != 0,,drop=FALSE];
      # Assign peakID: partID function from enviPick
      peakID <- .Call("partID",
                      as.integer(index),
                      as.integer(length(msobject$rawData[[mslevel]][[cE]]$peak)),
                      PACKAGE = "LipidMS")
      colnames(index) <- c("start","end","length")
      # save peaks
      msobject$rawData[[mslevel]][[cE]]$peak <- peakID
      msobject$processing[[mslevel]]$peakIndex[[cE]] <- index
    }
  }

  # create peaklist
  maxit <- max(msobject$rawData[[mslevel]][[cE]]$peak)
  peaklist <- data.frame()
  if (maxit > 0){
    for (p in 1:nrow( msobject$processing[[mslevel]]$peakIndex[[cE]])){
      start <- msobject$processing[[mslevel]]$peakIndex[[cE]][p,1]
      end <- msobject$processing[[mslevel]]$peakIndex[[cE]][p,2]

      mz <- mean(msobject$rawData[[mslevel]][[cE]]$mz[start:end])
      # mz <- weighted.mean(msobject$rawData[[mslevel]][[cE]]$mz[start:end],
      #                     msobject$rawData[[mslevel]][[cE]]$int[start:end])
      maxint <- max(msobject$rawData[[mslevel]][[cE]]$int[start:end])
      ints <- msobject$rawData[[mslevel]][[cE]]$int[start:end]
      ints <- ints - min(ints)
      sumint <- sum(ints)
      area <- areas[p]
      RT <- mean(msobject$rawData[[mslevel]][[cE]]$RT[start:end][msobject$rawData[[mslevel]][[cE]]$int[start:end] == maxint])
      minRT <- min(msobject$rawData[[mslevel]][[cE]]$RT[start:end])
      maxRT <- max(msobject$rawData[[mslevel]][[cE]]$RT[start:end])
      peakid <- p

      peaklist <- rbind(peaklist,
                        data.frame(mz = mz, RT = RT, int = sumint,
                                   minRT = minRT, maxRT = maxRT,
                                   peakID = peakid))
    }
    peaklist <- peaklist[order(peaklist$int, decreasing = TRUE),]
    peaklist$peakID <- paste(paste(mslevel, cE, sep="_"), peaklist$peakID, sep="_")
    msobject$rawData[[mslevel]][[cE]]$peak <- paste(paste(mslevel, cE, sep="_"), 
                                                    msobject$rawData[[mslevel]][[cE]]$peak, sep="_")
  } else {
    warning("No peaks found")
    peaklist <- data.frame()
  }
  msobject$peaklist[[mslevel]][[cE]] <- peaklist
  return(msobject)
}

# annotateIsotopes
#' Annotate isotopes
#'
#' Annotate isotopes based on mass differences, retention time and peak
#' correlation if required.
#'
#' @param peaklist extracted peaks. Data.frame with 4 columns (mz, RT, int
#' and peakID).
#' @param rawScans raw scan data. Data.frame with 5 columns (mz, RT, int,
#' peakID and Scan).
#' @param dmz mass tolerance in ppm.
#' @param drt RT windows with the same units used in peaklist.
#' @param massdiff mass difference.
#' @param charge charge.
#' @param isotopeAb isotope natural abundance.
#' @param m0mass mass of the most abundant naturally occurring stable isotope.
#' @param corThr peak correlation threshold.
#' @param checkInt logical. If TRUE, relative intensity is checked.
#' @param checkCor logical. If TRUE, peaks correlation is checked.
#'
#' @return peaklist with 6 columns (mz, RT, int, peakID, isotope and isoGroup).
#'
#' @keywords internal
#' 
#' @references Isotope annotation has been adapted from CAMERA algorithm: 
#' Kuhl C, Tautenhahn R, Boettcher C, Larson TR, Neumann S (2012). “CAMERA: an 
#' integrated strategy for compound spectra extraction and annotation of liquid 
#' chromatography/mass spectrometry data sets.” Analytical Chemistry, 84, 283–289. 
#' http://pubs.acs.org/doi/abs/10.1021/ac202450g.
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@iislafe.es>
annotateIsotopes <- function(peaklist, rawScans, dmz, drt,
                             massdiff, charge, isotopeAb, m0mass, 
                             corThr, checkInt, checkCor){
  
  peaklist <- peaklist[order(peaklist$mz, decreasing = F),]
  newpeaklist <- c()
  cluster <- 1 # counter
  
  # move isotope groups from peaklist to newpeaklist until peaklist is empty
  while (nrow(peaklist) > 0){
    # start from the first feature in peaklist (sorted by increasing mz)
    f <- peaklist[1,]
    
    # Search for coeluting features with greater mz than f
    ss <- peaklist[abs(peaklist$RT - f$RT) < drt & peaklist$mz >= f$mz,]
    
    if (nrow(ss) > 1){
      dm <- (ss$mz-f$mz)
      ss$isotope <- round(dm/massdiff/charge, 0)
      ##########################################################################
      # Check mass differences
      errors <- abs(dm-ss$isotope*(massdiff/charge))*1e6/(f$mz+round(dm, 0)*(massdiff/charge))
      group <- cbind(ss[errors < dmz, c(colnames(peaklist), "isotope"),])
      
      ##########################################################################
      # Check intensity ratios
      if (nrow(group) > 1){
        if(checkInt){
          keep <- rep(FALSE, nrow(group))
          keep[1] <- TRUE
          for (i in 2:nrow(group)){
            prev <- which(group$isotope - group$isotope[i] == -round(massdiff, 0)/charge)
            if(length(prev) > 0 & group$isotope[i] == round(massdiff, 0)/charge){
              intcheck <- any((group$int[i]/group$int[prev]) > isotopeAb &
                                (group$int[i]/group$int[prev]) < group$mz[prev[1]]*isotopeAb/m0mass)
              keep[i] <- intcheck
            } else if(length(prev) > 0 & group$isotope[i] > 1){
              keep[i] <- any(group$isotope - group$isotope[i] == -round(massdiff, 0)/charge)
            } else {
              break
            }
          }
          group <- group[keep,]
        }
        
        ##########################################################################
        # Check correlation
        if(checkCor){
          if(nrow(group) > 1){
            keep <- rep(FALSE, nrow(group))
            keep[1] <- TRUE
            group$cor <- sapply(group$peakID[1:nrow(group)], 
                                coelutionScore, group$peakID[1], rawScans)
            group <- group[group$cor >= corThr,]
          }
        } else {
          group$cor <- 0
        }
        
        ##########################################################################
        # Check overlaps. Keep the most correlated peak
        if (nrow(group) > 1){
          keep <- rep(TRUE, nrow(group))
          check_overlaps <- which(table(group$isotope) > 1)
          if (length(check_overlaps) > 0){
            repeated <- as.numeric(names(check_overlaps))
            for (rep in repeated){
              whichrep <- which(group$isotope == rep)
              keep[whichrep] <- FALSE
              closest <- whichrep[which.min(abs(f$RT - group$RT[group$isotope == rep]))]
              keep[closest] <- TRUE
            }
          }
          group <- group[keep,]
        } else {
          group$cor <- 0
        }
      } else {
        group$cor <- 0
      }
      
      if (nrow(group) > 1){
        group$isoGroup <- cluster
        cluster <- cluster + 1
      } else {
        group$isoGroup <- 0
      }
    } else {
      ss$isotope <- 0
      ss$cor <- 0
      group <- ss
      group$isoGroup <- 0
    }
    newpeaklist <- rbind(newpeaklist, group)
    peaklist <- peaklist[which(!peaklist$peakID %in% group$peakID),]
  }
  
  newpeaklist$isotope <- paste("[M+", newpeaklist$isotope, "]", sep="")
  newpeaklist$isotope[newpeaklist$isoGroup == 0] <- ""
  newpeaklist <- newpeaklist[order(newpeaklist$mz),]
  newpeaklist$cor <- NULL
  return(newpeaklist)
}

# getallpeaks
#' Extract peaks from all msobjects in a msbatch.
#' 
#' Extract peaks from all MS1 peaklists of the msobjects in a msbatch.
#'
#' @param msbatch msbatch object
#'
#' @return data.frame with 6 columns (mz, RT, int, peakID, isotope and isoGroup).
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@iislafe.es>
getallpeaks <- function(msbatch){
  # extract peaks
  peakspersample <- unlist(lapply(msbatch$msobjects, function(x) nrow(x$peaklist$MS1)))
  peaks <- data.frame(do.call(rbind, lapply(msbatch$msobjects, function(x) x$peaklist$MS1)))
  peaks$sample <- as.numeric(unlist(mapply(rep, 1:length(msbatch$msobjects), peakspersample)))
  peaks <- peaks[order(peaks$mz, peaks$RT, decreasing = FALSE),]
  return(peaks)
}

# indexrtpart
#' Index partitions or clusters assigned during alignment.
#' 
#' Index partitions or clusters assigned during alignment.
#'
#' @param peaks peaks obtained using \link{getallpeaks}.
#' @param part numeric vector with the partition/cluster assigned to each peak.
#' @param minsamples minimum number of samples represented in each partition.
#'
#' @return list with two elements: index and vector with the assigned partitions.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@iislafe.es>
indexrtpart <- function(peaks, part, minsamples){
  maxit <- max(part)
  
  i <- rle(part)
  end <- cumsum(i$lengths)
  start <- c(1, end[-length(end)] + 1)
  ind <- data.frame(start, end, length = i$lengths, value = i$values)
  ind <- ind[ind$value != 0,]
  keep <- c()
  for (i in 1:nrow(ind)){
    # keep only partitions which meet the alignment requirements
    s <- peaks$sample[ind$start[i]:ind$end[i]]
    if(length(unique(s)) >= minsamples){
      keep <- append(keep, TRUE)
    } else {
      keep <- append(keep, FALSE)
    }
  }
  ind$value[!keep] <- 0
  ind <- ind[keep,]
  ind$value <- as.numeric(as.factor(ind$value))
  
  pID <- rep(0, nrow(peaks))
  for(x in 1:nrow(ind)){
    pID[ind$start[x]:ind$end[x]] <- ind$value[x]
  }
  
  return(list(index = ind, idvector = pID))
}

# clustdist
#' Calculate max distance between clusters.
#' 
#' Calculate max distance between clusters.
#' 
#' @param mins lower bound for each cluster.
#' @param maxs higher bound for each cluster.
#' 
#' @return numeric vector with the assigned clusters
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@iislafe.es>
clustdist <- function(mins, maxs){
  # calculate max distance between 2 clusters
  cdiff <- matrix(nrow = length(mins), ncol = length(mins))
  for (x in 1:(length(mins)-1)){
    for(y in (x+1):length(maxs)){
      cdiff[y, x] <- max(abs(mins[x] - maxs[y]), abs(mins[y] - maxs[x]))
    }
  }
  return(cdiff)
}

# clust
#' Clustering for MS peaks based on mz or RT.
#' 
#' Clustering for MS peaks based on mz or RT.
#'
#' @param values values to clusterize (mz or RT).
#' @param mins lower bound for each value.
#' @param maxs higher bound for each value.
#' @param samples numeric vector that indicates to which sample each value belongs.
#' @param unique.samples logical. FALSE if several measures from the same sample 
#' can be clusterized together.
#' @param maxdist maximum distance allowed within a cluster.
#' @param ppm logical. TRUE if maxdist is given in ppm.
#'
#' @return numeric vector with the assigned clusters
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@iislafe.es>
clust <- function(values, mins, maxs, samples, unique.samples, maxdist, ppm){
  # values: vector of values to cluterize
  # mins: vector of minimum values
  # maxs: vector of maximum values
  # samples: vector indicating to which sample/cluster belongs each value from mins and maxs
  # values, mins, maxs and samples have the same length
  # unique.samples can be TRUE or FALSE (whether or not a cluster can contain different values from the same sample)
  # maxdist: maximum distance allowed
  # ppm: TRUE or FALSE if maxdist is in ppm
  
  if (missing(mins) | missing(maxs)){
    mins <- maxs <- values
  }
  if (length(mins) > 1){
    clust <- rep(0, length(mins)) # cluster id assigned
    n <- rep(1, length(mins)) # n peaks assigned to each cluster. Initialize the algortihm with as many clusters as points
    at <- list()
    at <- lapply(1:length(samples), function(x) at[[x]] <- samples[x]) # samples represented within each cluster
    atclust <- list()
    atclust <- lapply(1:length(samples), function(x) atclust[[x]] <- x) # values assigned to each cluster
    
    distmatrix <- clustdist(mins, maxs) # vector of distances calculated with clustdist
    distmatrix[distmatrix == -1] <- NA
    
    mindist <- which.min(distmatrix)
    n1 <- ceiling(mindist/length(mins))
    n2 <- mindist-(length(mins)*(n1-1))
    if (unique.samples){
      do <- FALSE
      while(!do){
        if (any(at[[n1]] %in% at[[n2]])){
          distmatrix[mindist] <- NA
          if (any(!is.na(distmatrix))){
            mindist <- which.min(distmatrix)
            n1 <- ceiling(mindist/length(mins))
            n2 <- mindist-(length(mins)*(n1-1))
          } else {
            do <- TRUE
          }
        } else {
          do <- TRUE
        }
      }
    }
    
    while(any(!is.na(distmatrix))){
      # condition 1 to join two clusters: dist n2-n1 is the minimum distance between n2 and any other cluster
      cond1 <- order(distmatrix[(length(mins)*(n1-1)+1):(length(mins)*n1)])[1] == n2
      
      # condition 2: dist n2-n1 is below maxdist 
      if (ppm == TRUE){ # if maxdist is in ppm
        dist <- abs(values[n2] - values[n1]) * 1e6 / values[n1]
        cond2 <- dist  <= maxdist 
      } else {
        cond2 <- abs(values[n2] - values[n1])  <= maxdist 
      }
      
      if(cond1 & cond2){ # if both conditions are true, join clusters y remove n2
        mins[n1] <- min(mins[n1], mins[n2])
        maxs[n1] <- max(maxs[n1], maxs[n2])
        values[n1] <- (values[n1] * n[n1] + values[n2] * n[n2])/(n[n1] + n[n2]) # mean value
        n[n1] <- n[n1] + n[n2]
        mins <- mins[-n2]
        maxs <- maxs[-n2]
        values <- values[-n2]
        samples <- samples[-n2]
        at[[n1]] <- append(at[[n1]], at[[n2]])
        at[[n2]] <- NULL
        atclust[[n1]] <- append(atclust[[n1]], atclust[[n2]])
        atclust[[n2]] <- NULL
        
        if (length(mins) > 1){ # update distances between clusters
          distmatrix <- clustdist(mins, maxs)
          
          mindist <- which.min(distmatrix)
          n1 <- ceiling(mindist/length(mins))
          n2 <- mindist-(length(mins)*(n1-1))
          any(at[[n1]] %in% at[[n2]])
          if (unique.samples){
            do <- FALSE
            while(!do){
              if (any(at[[n1]] %in% at[[n2]])){
                distmatrix[mindist] <- NA
                if (any(!is.na(distmatrix))){
                  mindist <- which.min(distmatrix)
                  n1 <- ceiling(mindist/length(mins))
                  n2 <- mindist-(length(mins)*(n1-1))
                } else {
                  do <- TRUE
                }
              } else {
                do <- TRUE
              }
            }
          }
        } else {
          distmatrix[mindist] <- NA
        }
      } else {
        distmatrix[mindist] <- NA
        if(any(!is.na(distmatrix))){
          mindist <- which.min(distmatrix)
          n1 <- ceiling(mindist/length(mins))
          n2 <- mindist-(length(mins)*(n1-1))
          if (unique.samples){
            do <- FALSE
            while(!do){
              if (any(at[[n1]] %in% at[[n2]])){
                distmatrix[mindist] <- NA
                if (any(!is.na(distmatrix))){
                  mindist <- which.min(distmatrix)
                  n1 <- ceiling(mindist/length(mins))
                  n2 <- mindist-(length(mins)*(n1-1))
                } else {
                  do <- TRUE
                }
              } else {
                do <- TRUE
              }
            }
          }
        }
      }
    }
    for (c in 1:length(atclust)){
      pos <- atclust[[c]]
      clust[pos] <- c
    }
  } else {
    clust <- 1
  }
  return(clust)
}

# rtcorrection
#' Correct RT based on a rtmodel.
#' 
#' Correct RT based on a rtmodel.
#'
#' @param rt rt vector to correct.
#' @param rtmodel rt model.
#'
#' @return corrected rt vector.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@iislafe.es>
rtcorrection <- function(rt, rtmodel){
  rtdevsmoothed <- predict(rtmodel, rt)
  rtdevsmoothed[is.na(rtdevsmoothed)] <- 0
  rt <- rt - rtdevsmoothed
  return(rt)
}

# getfeaturestable
#' Write features table based on groups
#' 
#' Write features table based on groups
#'
#' @param msbatch msbatch object
#'
#' @return data.frame
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@iislafe.es>
getfeaturestable <- function(msbatch){
  # check msbatch structure
  if (!is.list(msbatch) | !all(names(msbatch) %in% c("metaData", "msobjects", "alignment", "grouping", "features")) | 
      !is.data.frame(msbatch$metaData) | !is.list(msbatch$msobjects) | !is.list(msbatch$alignment) | 
      !is.list(msbatch$grouping) | !is.data.frame(msbatch$features)){
    stop("Wrong msbatch format")
  }
  # check that all msobjects have an mslevel 1
  whichmslevel1 <- which(unlist(lapply(msbatch$msobjects, function(x) 
    1 %in% unique(x$metaData$scansMetadata$msLevel))))
  if (length(whichmslevel1) != nrow(msbatch$metaData)){
    warning("Removing samples with no MS1 level for alignment")
    msbatch$metaData <- msbatch$metaData[whichmslevel1,]
    msbatch$msobjects <- msbatch$msobjects[whichmslevel1]
  }
  
  if(msbatch$grouping$grouped != TRUE | nrow(msbatch$grouping$groupIndex) <= 1){
    stop("Grouping must be performed before running getFeatureTable function
         \n(use groupmsbatch(msbatch)) or not enough groups have been identified")
  }
  
  featureMatrix <- matrix(nrow = nrow(msbatch$grouping$groupIndex), 
                          ncol = length(msbatch$msobjects))
  mz <- rep(0, nrow(featureMatrix))
  minmz <- rep(0, nrow(featureMatrix))
  maxmz <- rep(0, nrow(featureMatrix))
  RT <- rep(0, nrow(featureMatrix))
  minRT <- rep(0, nrow(featureMatrix))
  maxRT <- rep(0, nrow(featureMatrix))
  iniRT <- rep(0, nrow(featureMatrix))
  endRT <- rep(0, nrow(featureMatrix))
  n <- rep(0, nrow(featureMatrix))
  group <- 1:nrow(featureMatrix)
  # isotope <- rep("", nrow(featureMatrix))
  
  for (g in 1:nrow(msbatch$grouping$groupIndex)){
    gr <- msbatch$grouping$peaks[msbatch$grouping$peaks$groupID == g,]
    for (s in as.numeric(unique(gr$sample))){
      featureMatrix[g, s] <- sum(as.numeric(gr$int[gr$sample == s])) # sum() added
      n[g] <- n[g] + 1
    }
    mz[g] <- mean(gr$mz)
    minmz[g] <- min(gr$mz)
    maxmz[g] <- max(gr$mz)
    RT[g] <- mean(gr$RT)
    minRT[g] <- min(gr$RT, na.rm = TRUE)
    maxRT[g] <- max(gr$RT, na.rm = TRUE)
    iniRT[g] <- median(gr$minRT, na.rm = TRUE)
    endRT[g] <- median(gr$maxRT, na.rm = TRUE)
    # iso <- table(gr$isotope)
    # isotope[g] <- names(which.max(iso[!is.na(iso)]))
  }
  
  colnames(featureMatrix) <- msbatch$metaData$sample
  
  # Isotopes and isotopegroups
  index <- matrix(nrow = length(mz), 
                  ncol = length(msbatch$msobjects))
  peaks <- msbatch$grouping$peaks
  
  for (g in 1:length(mz)){
    gr <- which(peaks$groupID == g)
    for (s in as.numeric(unique(peaks$sample[gr]))){
      index[g, s] <- max(gr[which(peaks$sample[gr] == s)]) # max() added 
    }
  }
  
  # isotope
  isotopes <- t(apply(index, 1, function(x) peaks$isotope[as.numeric(x)]))
  isotopes[isotopes == ""] <- NA
  isotopes <- apply(isotopes, 1, function(x){
    i <- as.numeric(gsub("\\[|\\]|M|\\+", "", x))
    if (any(!is.na(i))){
      return(paste("[M+", min(i, na.rm = TRUE), "]", sep=""))
    } else {
      return("")
    }
  })
  
  # groups
  gind <- t(apply(index, 1, function(x) peaks$isoGroup[as.numeric(x)]))
  gind[gind == 0] <- NA
  
  isogroup <- rep(0, length(mz))
  counter <- 1
  
  for (s in 1:ncol(gind)){
    groups <- unique(gind[!is.na(gind[,s]),s])
    if (length(groups) > 0){
      for (g in groups){
        ig <- which(gind[,s] == g)
        if (any(isogroup[ig] != 0)){
          pg <- unique(isogroup[ig][isogroup[ig] != 0])
          if (length(pg) > 1){
            ig <- unique(c(ig, which(isogroup %in% pg)))
            pg <- min(pg)
          }
          isogroup[ig] <- pg
        } else {
          isogroup[ig] <- counter
          counter <- counter + 1
        }
      }
    }
  }
  
  # save results in msbatch
  rtminutes <- round(RT/60, 2)
  featureTable <- data.frame(mz, minmz, maxmz, RT, minRT, maxRT, iniRT, endRT, 
                             rtminutes, npeaks = n, group, isotopes, isogroup, 
                             featureMatrix)
  rownames(featureTable) <- make.names(paste(round(mz, 3), round(RT, 0), sep="_"), 
                                       unique = TRUE)
  
  colnames(featureTable)[1:13] <- c("mz", "minmz", "maxmz", "RT", "minRT", "maxRT",
                                    "iniRT", "endRT", "RTminutes", "npeaks", 
                                    "group", "isotope", "isoGroup")
  
  msbatch$features <- featureTable
  
  return(msbatch)
}

# removeduplicatedpeaks
#' Remove duplicated features after grouping step
#' 
#' Remove duplicated features after grouping step
#'
#' @param msbatch msbatch object 
#' @param ppm mz tolerance in ppm
#'
#' @return msbatch
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@iislafe.es>
removeduplicatedpeaks <- function(msbatch, 
                                  ppm, 
                                  dmz = 5, 
                                  drt = 30,
                                  thr_overlap = 0.7){
  
  peaks <- msbatch$features
  peaks[is.na(peaks)] <- 0
  idpeaks <- 1:nrow(peaks)
  
  #============================================================================#
  # Create mz partitions based on dmz and drt
  #============================================================================#
  part <- .Call("agglom", as.numeric(peaks$mz),
                as.numeric(peaks$RT), as.integer(1),
                as.numeric(dmz), as.numeric(drt),
                PACKAGE = "LipidMS")
  
  #============================================================================#
  # For partitions with more than one feature, check if there is RT overlapping
  #============================================================================#
  dup <- as.numeric(names(table(part))[which(table(part) > 1)])
  
  count <- 0
  toremove <- c()
  tokeep <- c()
  for (d in dup){
    ss <- peaks[part == d,]
    ss <- ss[order(ss$iniRT, decreasing = FALSE),]
    # ss[,1:11]
    # check overlapping
    end <- FALSE
    last <- list()
    while (nrow(ss) > 1 & !end){
      tocompare <- sapply(2:nrow(ss), function(x){
        ss$group[which(ss$iniRT[x] < ss$endRT[1:(x-1)])]
      }, simplify = FALSE)
      ncomp <- sum(sapply(tocompare, length) > 0)
      icomp <- 0
      if (length(last) == length(tocompare)){
        if(all(unlist(last) == unlist(tocompare))){
          keepgoing <- FALSE
        } else {
          keepgoing <- TRUE
        }
      } else {
        keepgoing <- TRUE
      }
      if (any(sapply(tocompare, length) > 0) & !end & keepgoing){
        f <- which(sapply(tocompare, length) > 0)
        
        # overlapping %
        for (i in f){
          ss2 <- ss[!is.na(ss$mz),]
          g1 <- tocompare[[i]][1]
          g2 <- ss$group[i+1]
          
          if(g1 %in% ss2$group & g2 %in% ss2$group & 
             !(g1 %in% toremove) & !(g2 %in% toremove)){
            p1 <- which(ss2$group == g1)
            p2 <- which(ss2$group == g2)
            width_p1 <- ss2$endRT[p1] - ss2$iniRT[p1]
            width_p2 <- ss2$endRT[p2] - ss2$iniRT[p2]
            time_overlap <- (min(ss2$endRT[p1], ss2$endRT[p2]) -
                               max(ss2$iniRT[p1], ss2$iniRT[p2]))
            perc_overlap <- time_overlap / c(width_p1, width_p2)
            
            if (any(perc_overlap > thr_overlap)){
              rem <- c(g1, g2)[which.max(perc_overlap)]
              keep <- c(g1, g2)[which.min(perc_overlap)]
              
              # reorder raw data
              msbatch$grouping$peaks$groupID[msbatch$grouping$groupIndex$start[rem]:
                                               msbatch$grouping$groupIndex$end[rem]] <- keep
              
              toremove <- c(toremove, rem)
              tokeep <- c(tokeep, keep)
              count <- count + 1
              
              ss[which(ss$group == rem), 1:11] <- NA
              icomp <- icomp + 1
              
              if (rem == g1 & nrow(ss[(!is.na(ss$mz)) & (!ss$group %in% toremove),, 
                                      drop=FALSE]) > 1){
                ss <- ss[!is.na(ss$mz),,drop=FALSE]
                ss <- ss[!ss$group %in% toremove,,drop=FALSE]
                break
              }
              
            } else {
              icomp <- icomp + 1
            }
          } else {
            icomp <- icomp + 1
          }
        }
        last <- tocompare
        if (icomp == ncomp){
          end <- TRUE
        }
      } else {
        end <- TRUE
      }
    }
  }
  
  # reorder peaks
  msbatch$grouping$peaks <- msbatch$grouping$peaks[order(msbatch$grouping$peaks$groupID),]
  
  # create new index and feature table
  i <- rle(msbatch$grouping$peaks$groupID)
  end <- cumsum(i$lengths)
  start <- c(1, end[-length(end)] + 1)
  ind <- data.frame(start, end, length = i$lengths, value = i$values)
  ind <- ind[ind$value != 0,]
  ind$value <- as.numeric(as.factor(ind$value))
  pID <- rep(0, nrow(msbatch$grouping$peaks))
  for(x in 1:nrow(ind)){
    pID[ind$start[x]:ind$end[x]] <- ind$value[x]
  }
  msbatch$grouping$peaks$groupID <- pID
  msbatch$grouping$groupIndex <- ind
  
  msbatch <- getfeaturestable(msbatch)
  
  return(msbatch)
}

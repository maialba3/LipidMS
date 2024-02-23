# idPOS
#' Lipids annotation for ESI+
#'
#' Lipids annotation based on fragmentation patterns for LC-MS/MS DIA or DDA data
#' acquired in positive mode. This function compiles all functions written for
#' ESI+ annotations.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 5 seconds.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param lipidClasses classes of interest to run the identification functions.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification); and the annotatedPeaklist element shows the original 
#' MS1 peaklist with the annotations on it.
#'
#' @examples
#' \dontrun{
#' msobject <- idPOS(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idPOS <- function(msobject,
                  ppm_precursor = 5,
                  ppm_products = 10,
                  rttol = 5,
                  coelCutoff = 0.8,
                  lipidClasses = c("MG", "LPC", "LPE", "PC", "PCo", "PCp", "PE", 
                                   "PEo", "PEp", "PG", "PI", "Sph", "SphP", "Cer", 
                                   "AcylCer", "CerP", "SM", "Carnitines", "CE", 
                                   "DG", "TG"),
                  dbs,
                  verbose = TRUE){

  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if ("results" %in% names(msobject$annotation)){
    if(verbose){cat("\n Removing previous results...")}
    msobject$annotation$results <- NULL
    msobject$annotation$detailsAnnotation <- NULL
    msobject$annotation$annotatedPeaklist <- NULL
    if(verbose){cat("OK")}
  }
  if (!all(lipidClasses %in% c("MG", "LPC", "LPE", "PC", "PCo", "PCp", "PE", 
                               "PEo", "PEp", "PG", "PI", "Sph", "SphP", "Cer", 
                               "AcylCer", "CerP", "SM", "Carnitines", "CE", "DG", 
                               "TG"))){
    stop("Lipid classes allowed for positive annotation are: MG, LPC, LPE, PC,
          PCo, PCp, PE, PEo, PEp, PG, PI, Sph, SphP, Cer, CerP, AcylCer, SM, 
          Carnitines, CE, DG and TG")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }

  if(verbose){cat("\n Starting annotation...")}
  if ("MG" %in% lipidClasses){
    if(verbose){cat("\n  Searching for MG...")}
    msobject <-  idMGpos(msobject = msobject, ppm_precursor = ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("LPC" %in% lipidClasses){
    if(verbose){cat("\n  Searching for LPC...")}
    msobject <-  idLPCpos(msobject = msobject, ppm_precursor = ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("LPE" %in% lipidClasses){
    if(verbose){cat("\n  Searching for LPE...")}
    msobject <-  idLPEpos(msobject = msobject, ppm_precursor = ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
    
  }
  if ("PC" %in% lipidClasses){
    if(verbose){cat("\n  Searching for PC...")}
    msobject <-  idPCpos(msobject = msobject, ppm_precursor = ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("PCo" %in% lipidClasses){
    if(verbose){cat("\n  Searching for PCo...")}
    msobject <-  idPCopos(msobject = msobject, ppm_precursor = ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("PCp" %in% lipidClasses){
    if(verbose){cat("\n  Searching for PCp...")}
    msobject <-  idPCppos(msobject = msobject, ppm_precursor = ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("PE" %in% lipidClasses){
    if(verbose){cat("\n  Searching for PE...")}
    msobject <-  idPEpos(msobject = msobject, ppm_precursor = ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
    
  }
  if ("PEo" %in% lipidClasses){
    if(verbose){cat("\n  Searching for PEo...")}
    msobject <-  idPEopos(msobject = msobject, ppm_precursor = ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("PEp" %in% lipidClasses){
    if(verbose){cat("\n  Searching for PEp...")}
    msobject <-  idPEppos(msobject = msobject, ppm_precursor = ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("PG" %in% lipidClasses){
    if(verbose){cat("\n  Searching for PG...")}
    msobject <-  idPGpos(msobject = msobject, ppm_precursor = ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("PI" %in% lipidClasses){
    if(verbose){cat("\n  Searching for PI...")}
    msobject <-  idPIpos(msobject = msobject, ppm_precursor = ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("Sph" %in% lipidClasses){
    if(verbose){cat("\n  Searching for Sph...")}
    msobject <-  idSphpos(msobject = msobject, ppm_precursor = ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("SphP" %in% lipidClasses){
    if(verbose){cat("\n  Searching for SphP...")}
    msobject <-  idSphPpos(msobject = msobject, ppm_precursor = ppm_precursor,
                           ppm_products = ppm_products, rttol = rttol,
                           coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("Cer" %in% lipidClasses){
    if(verbose){cat("\n  Searching for Cer...")}
    msobject <-  idCerpos(msobject = msobject, ppm_precursor = ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("CerP" %in% lipidClasses){
    if(verbose){cat("\n  Searching for CerP...")}
    msobject <-  idCerPpos(msobject = msobject, ppm_precursor = ppm_precursor,
                           ppm_products = ppm_products, rttol = rttol,
                           coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("AcylCer" %in% lipidClasses){
    if(verbose){cat("\n  Searching for AcylCer...")}
    msobject <-  idAcylCerpos(msobject = msobject, ppm_precursor = ppm_precursor,
                              ppm_products = ppm_products, rttol = rttol,
                              coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("SM" %in% lipidClasses){
    if(verbose){cat("\n  Searching for SM...")}
    msobject <-  idSMpos(msobject = msobject, ppm_precursor = ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("Carnitines" %in% lipidClasses){
    if(verbose){cat("\n  Searching for Carnitines...")}
    msobject <-  idCarpos(msobject = msobject, ppm_precursor = ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("CE" %in% lipidClasses){
    if(verbose){cat("\n  Searching for CE...")}
    msobject <-  idCEpos(msobject = msobject, ppm_precursor = ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("DG" %in% lipidClasses){
    if(verbose){cat("\n  Searching for DG...")}
    msobject <-  idDGpos(msobject = msobject, ppm_precursor = ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if ("TG" %in% lipidClasses){
    if(verbose){cat("\n  Searching for TG...")}
    msobject <-  idTGpos(msobject = msobject, ppm_precursor = ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs, verbose = verbose)
    if(verbose){cat("OK")}
  }
  if(verbose){cat("\n Preparing output...")}
  msobject <- crossTables(msobject,
                          ppm = ppm_precursor, 
                          rttol = rttol,
                          dbs = dbs)
  if(verbose){cat("OK\n")}
  return(msobject)
}

# idMGpos
#' Monoacylglycerol (MG) annotation for ESI+
#'
#' MG identification based on fragmentation patterns for LC-MS/MS DIA and DDA data
#' acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for MG in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idMGpos} function involves 2 steps. 1) FullMS-based
#' identification of candidate MG as M+H-H2O, M+NH4 and M+Na. 2) Search of
#' MG class fragments if any is assigned.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (in this
#' case, just MS-only or Subclass level (if any class fragment is defined) are
#' possible) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idMGpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idMGpos <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("M+H-H2O", "M+NH4", "M+Na"),
                    clfrags = c(),
                    clrequired = c(),
                    ftype = c(),
                    coelCutoff = 0.8,
                    dbs,
                    verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous MG annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "MG",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("MG" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious MG annotations removed")}
      msobject$annotation$detailsAnnotation$MG <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$mgdb, ppm = ppm_precursor,
                               rt = rt, adducts = adducts, rttol = rttol,
                               dbs = dbs, rawData = rawData,
                               coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    # isolation of coeluting fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb = list(),
                           intrules  = c(), intConf = list(), nchains = 0,
                           class="MG",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$MG <- list()
    msobject$annotation$detailsAnnotation$MG$candidates <- candidates
    msobject$annotation$detailsAnnotation$MG$classfragments <- classConf$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$MG$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$MG <- list()
  }
  return(msobject)
}

# idLPCpos
#' Lysophosphocholines (LPC) annotation for ESI+
#'
#' LPC identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for LPC in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments. See \link{chainFrags} for details.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idLPCpos} function involves 3 steps. 1) FullMS-based
#' identification of candidate LPC as M+H and M+Na. 2) Search of
#' LPC class fragments: 104.1075 and 184.0739 coeluting with the
#' precursor ion. 3) Search of specific fragments that confirm chain
#' composition (MG as M+H-H2O).
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (in this
#' case, as LPC only have one chain, only Subclass and FA level are possible)
#' and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idLPCpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idLPCpos <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M+H", "M+Na"),
                     clfrags = c(104.1075, 184.0739),
                     clrequired = c(F, F),
                     ftype = c("F", "F"),
                     chainfrags_sn1 = c("mg_M+H-H2O"),
                     coelCutoff = 0.8,
                     dbs,
                     verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous LPC annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "LPC",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("LPC" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious LPC annotations removed")}
      msobject$annotation$detailsAnnotation$LPC <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$lysopcdb, ppm = ppm_precursor,
                               rt = rt, adducts = adducts, rttol = rttol,
                               dbs = dbs, rawData = rawData,
                               coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=1, sn1)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules = c(), rates = c(),
                                   intrequired = c(), nchains=1,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb,
                           intrules  = c(), intConf, nchains = 1, class="LPC",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$LPC <- list()
    msobject$annotation$detailsAnnotation$LPC$candidates <- candidates
    msobject$annotation$detailsAnnotation$LPC$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$LPC$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$LPC$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$LPC <- list()
  }
  return(msobject)
}

# idLPEpos
#' Lysophosphoethanolamines (LPE) annotation for ESI+
#'
#' LPE identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for LPE in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments. See \link{chainFrags} for details.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idLPEpos} function involves 3 steps. 1) FullMS-based
#' identification of candidate LPE as M+H and M+Na. 2) Search of
#' LPE class fragments: neutral loss of 141.01909 coeluting with the
#' precursor ion. 3) Search of specific fragments that confirm chain
#' composition in sn1 (MG as M+H-H2O).
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (in this
#' case, as LPE only have one chain, only Subclass and FA level are possible)
#' and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idLPEpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idLPEpos <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M+H", "M+Na"),
                     clfrags = c(141.01909),
                     clrequired = c(F),
                     ftype = c("NL"),
                     chainfrags_sn1 = c("mg_M+H-H2O"),
                     coelCutoff = 0.8,
                     dbs,
                     verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "LPE",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("LPE" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious LPE annotations removed")}
      msobject$annotation$detailsAnnotation$LPE <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$lysopedb, ppm = ppm_precursor,
                               rt = rt, adducts = adducts, rttol = rttol,
                               dbs = dbs, rawData = rawData,
                               coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=1, sn1)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules = c(), rates = c(),
                                   intrequired = c(), nchains=1,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb,
                           intrules  = c(), intConf, nchains = 1, class="LPE",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$LPE <- list()
    msobject$annotation$detailsAnnotation$LPE$candidates <- candidates
    msobject$annotation$detailsAnnotation$LPE$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$LPE$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$LPE$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$LPE <- list()
  }
  return(msobject)
}

# idPCpos
#' Phosphocholines (PC) annotation for ESI+
#'
#' PC identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for PC in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param rates character vector with the expected rates between fragments given
#' as a string (e.g. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idPCpos} function involves 5 steps. 1) FullMS-based
#' identification of candidate PC as M+H and M+Na. 2) Search of PC class
#' fragments: 104.1075, 184.0739 and neutral loss of 183.06604 coeluting with
#' the precursor ion. 3) Search of specific fragments that inform about chain
#' composition in sn1 (lysoPC as M+H or M+H-H2O resulting from the loss of the
#' FA chain at sn2) and sn2 (lysoPC as M+H or M+H-H2O resulting from the loss of
#' the FA chain at sn1 or the difference between precursor and sn1 chain
#' fragments). 4) Look for possible chains structure based on the combination of
#' chain fragments. 5) Check intensity rules to confirm chains position. In this
#' case, lysoPC from sn1 is at least twice more intense than lysoPC from sn2.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idPCpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idPCpos <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("M+H", "M+Na"),
                    clfrags = c(104.1075, 184.0739, 183.06604),
                    clrequired = c(F, F, F),
                    ftype = c("F", "F", "NL"),
                    chainfrags_sn1 = c("lysopc_M+H", "lysopc_M+H-H2O"),
                    chainfrags_sn2 = c("lysopc_M+H", "lysopc_M+H-H2O", ""),
                    intrules = c("lysopc_sn1/lysopc_sn2"),
                    rates = c("2/1"),
                    intrequired = c(T),
                    coelCutoff = 0.8,
                    dbs,
                    verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous PC annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "PC",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("PC" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious PC annotations removed")}
      msobject$annotation$detailsAnnotation$PC <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$pcdb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=2,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="PC",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PC <- list()
    msobject$annotation$detailsAnnotation$PC$candidates <- candidates
    msobject$annotation$detailsAnnotation$PC$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$PC$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$PC$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PC <- list()
  }
  return(msobject)
}

# idPCopos
#' Plasmanyl Phosphocholines (PCo) annotation for ESI+
#'
#' PCo identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for PC in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param rates character vector with the expected rates between fragments given
#' as a string (e.g. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idPCopos} function involves 5 steps. 1) FullMS-based
#' identification of candidate PCo as M+H and M+Na. 2) Search of PC class
#' fragments: 104.1075, 184.0739 and neutral loss of 183.06604 coeluting with
#' the precursor ion. 3) Search of specific fragments that inform about chain
#' composition in sn1 (LPCo as M+H or M+H-H2O resulting from the loss of the
#' FA chain at sn2) and sn2 (LPC as M+H-H2O resulting from the loss of
#' the FA chain at sn1 or the difference between precursor and sn1 chain
#' fragments). 4) Look for possible chains structure based on the combination of
#' chain fragments. 5) Check intensity rules to confirm chains position. In this
#' case, LPCo from sn1 is at least twice more intense than LPC from sn2.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idPCopos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idPCopos <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M+H", "M+Na"),
                     clfrags = c(104.1075, 184.0739, 183.06604),
                     clrequired = c(F, F, F),
                     ftype = c("F", "F", "NL"),
                     chainfrags_sn1 = c("lysopco_M+H", "lysopco_M+H-H2O"),
                     chainfrags_sn2 = c("lysopc_M+H", "lysopc_M+H-H2O", ""),
                     intrules = c("lysopco_sn1/lysopc_sn2"),
                     rates = c("2/1"),
                     intrequired = c(T),
                     coelCutoff = 0.8,
                     dbs,
                     verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous PC annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "PCo",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("PCo" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious PCo annotations removed")}
      msobject$annotation$detailsAnnotation$PCo <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$pcodb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=2,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="PCo",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PCo <- list()
    msobject$annotation$detailsAnnotation$PCo$candidates <- candidates
    msobject$annotation$detailsAnnotation$PCo$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$PCo$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$PCo$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PCo <- list()
  }
  return(msobject)
}

# idPCppos
#' Plasmenyl Phosphocholines (PCp) annotation for ESI+
#'
#' PCp identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for PC in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param rates character vector with the expected rates between fragments given
#' as a string (e.g. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idPCppos} function involves 5 steps. 1) FullMS-based
#' identification of candidate PC as M+H and M+Na. 2) Search of PC class
#' fragments: 104.1075, 184.0739 and neutral loss of 183.06604 coeluting with
#' the precursor ion. 3) Search of specific fragments that inform about chain
#' composition in sn1 (LPCp as M+H or M+H-H2O resulting from the loss of the
#' FA chain at sn2) and sn2 (LPC as M+H-H2O resulting from the loss of
#' the FA chain at sn1 or the difference between precursor and sn1 chain
#' fragments). 4) Look for possible chains structure based on the combination of
#' chain fragments. 5) Check intensity rules to confirm chains position. In this
#' case, LPC from sn2 is at least twice more intense than LPCo from sn1.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idPCppos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idPCppos <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M+H", "M+Na"),
                     clfrags = c(104.1075, 184.0739, 183.06604),
                     clrequired = c(F, F, F),
                     ftype = c("F", "F", "NL"),
                     chainfrags_sn1 = c("lysopcp_M+H", "lysopcp_M+H-H2O"),
                     chainfrags_sn2 = c("lysopc_M+H-H2O", ""),
                     intrules = c("lysopcp_sn1/lysopc_sn2"),
                     rates = c("1/2"),
                     intrequired = c(T),
                     coelCutoff = 0.8,
                     dbs,
                     verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous PC annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "PCp",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("PCp" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious PCp annotations removed")}
      msobject$annotation$detailsAnnotation$PCp <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$pcpdb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=2,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="PCp",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PCp <- list()
    msobject$annotation$detailsAnnotation$PCp$candidates <- candidates
    msobject$annotation$detailsAnnotation$PCp$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$PCp$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$PCp$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PCp <- list()
  }
  return(msobject)
}

# idPEpos
#' Phosphoethanolamines (PE) annotation for ESI+
#'
#' PE identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for PE in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param rates character vector with the expected rates between fragments given
#' as a string (e.g. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See
#' \link{createLipidDB} and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idPEpos} function involves 5 steps. 1) FullMS-based
#' identification of candidate PE as M+H and M+Na. 2) Search of PE class
#' fragments: loss of head group (DG as M+H-H2O) coeluting with the precursor
#' ion. 3) Search of specific fragments that inform about chain composition at
#' sn1 (MG as M+H-H2O resulting from the loss of the FA chain at sn2 and
#' the head group or LPE as M+H-H2O resulting just from the loss of the FA chain)
#' and sn2 (MG as M+H-H2O resulting from the loss of the head group and FA chain 
#' from sn2). 4) Look for possible chains structure based on the combination of 
#' chain fragments. 5) Check intensity rules to confirm chains position. 
#' LPE or MG from sn1 is at least 3 times more intense than the ones from sn2.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idPEpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idPEpos <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("M+H", "M+Na"),
                    clfrags = c("dg_M+H-H2O"),
                    clrequired = c(F),
                    ftype = c("BB"),
                    chainfrags_sn1 = c("lysope_M+H-H2O", "mg_M+H-H2O"),
                    chainfrags_sn2 = c("mg_M+H-H2O"),
                    intrules = c("lysope_sn1/lysope_sn1", "mg_sn1/mg_sn2"),
                    rates = c("3/1", "1/2"),
                    intrequired = c(F, F),
                    coelCutoff = 0.8,
                    dbs,
                    verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "PE",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("PE" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious PE annotations removed")}
      msobject$annotation$detailsAnnotation$PE <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$pedb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=2,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="PE",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PE <- list()
    msobject$annotation$detailsAnnotation$PE$candidates <- candidates
    msobject$annotation$detailsAnnotation$PE$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$PE$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$PE$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PE <- list()
  }
  return(msobject)
}

# idPEopos
#' Plasmanyl Phosphoethanolamines (PEo) annotation for ESI+
#'
#' PEo identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for PE in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param rates character vector with the expected rates between fragments given
#' as a string (e.g. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See
#' \link{createLipidDB} and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idPEopos} function involves 5 steps. 1) FullMS-based
#' identification of candidate PE as M+H and M+Na. 2) Search of PE class
#' fragments: loss of head group (NL of 141.0193) coeluting with the precursor
#' ion. 3) Search of specific fragments that inform about chain composition at
#' sn1 (LPEo as M+H or M+H-H2O resulting from the loss of the FA chain at sn2)
#' and sn2 (MG as M+H-H2O resulting just from the loss of the head group and the 
#' FA chain at sn1). 4) Look for possible chains structure based on the
#' combination of chain fragments. 5) Check intensity rules to confirm chains
#' position. LPEo from sn1 is at least 2 times more intense than MG from sn2.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idPEopos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idPEopos <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M+H", "M+Na"),
                     clfrags = c(141.0193),
                     clrequired = c(F),
                     ftype = c("NL"),
                     chainfrags_sn1 = c("lysopeo_M+H", "lysopeo_M+H-H2O"),
                     chainfrags_sn2 = c("mg_M+H-H2O"),
                     intrules = c("lysopeo_sn1/mg_sn2"),
                     rates = c("2/1"),
                     intrequired = c(T),
                     coelCutoff = 0.8,
                     dbs,
                     verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "PEo",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("PEo" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious PEo annotations removed")}
      msobject$annotation$detailsAnnotation$PEo <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$peodb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=2,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="PEo",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PEo <- list()
    msobject$annotation$detailsAnnotation$PEo$candidates <- candidates
    msobject$annotation$detailsAnnotation$PEo$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$PEo$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$PEo$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PEo <- list()
  }
  return(msobject)
}

# idPEppos
#' Plasmenyl Phosphoethanolamines (PEp) annotation for ESI+
#'
#' PEp identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for PE in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param rates character vector with the expected rates between fragments given
#' as a string (e.g. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See
#' \link{createLipidDB} and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idPEppos} function involves 5 steps. 1) FullMS-based
#' identification of candidate PE as M+H and M+Na. 2) Search of PE class
#' fragments: loss of head group (NL of 140.012) coeluting with the precursor
#' ion. 3) Search of specific fragments that inform about chain composition at
#' sn1 (LPEp as M+H or M+H-H2O resulting from the loss of the FA chain at sn2)
#' and sn2 (MG as M+H-H2O from sn2 resulting from the loss of the FA chain at sn1). 
#' 4) Look for possible chains structure based on the combination of chain 
#' fragments. 5) Check intensity rules to confirm chains position. MG from sn2 
#' is at least 3 times more intense than LPEp from sn1.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idPEppos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idPEppos <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M+H", "M+Na"),
                     clfrags = c(140.012),
                     clrequired = c(F),
                     ftype = c("NL"),
                     chainfrags_sn1 = c("lysopep_M+H", "lysopep_M+H-H2O"),
                     chainfrags_sn2 = c("mg_M+H-H2O"),
                     intrules = c("lysopep_sn1/mg_sn2"),
                     rates = c("1/3"),
                     intrequired = c(T),
                     coelCutoff = 0.8,
                     dbs, 
                     verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "PEp",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("PEp" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious PEp annotations removed")}
      msobject$annotation$detailsAnnotation$PEp <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$pepdb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=2,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="PEp",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PEp <- list()
    msobject$annotation$detailsAnnotation$PEp$candidates <- candidates
    msobject$annotation$detailsAnnotation$PEp$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$PEp$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$PEp$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PEp <- list()
  }
  return(msobject)
}

# idPGpos
#' Phosphoglycerols (PG) annotation for ESI+
#'
#' PG identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for PE in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param rates character vector with the expected rates between fragments given
#' as a string (e.g. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See
#' \link{createLipidDB} and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idPGpos} function involves 5 steps. 1) FullMS-based
#' identification of candidate PG as M+H, M+NH4 and M+Na. 2) Search of PG class
#' fragments: loss of head group (DG as M+H-H2O) coeluting with the precursor
#' ion. 3) Search of specific fragments that inform about chain composition at
#' sn1 (MG as M+H-H2O resulting from the loss of the FA chain at sn2)
#' and sn2 (MG as M+H-H2O resulting from the loss of the FA chain at sn1).
#' 4) Look for possible chains structure based on the combination of chain
#' fragments. 5) Check intensity rules to confirm chains position. MG from sn2
#' is at least twice more intense than the one from sn1.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idPGpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idPGpos <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("M+H", "M+NH4", "M+Na"),
                    clfrags = c("dg_M+H-H2O"),
                    clrequired = c(F),
                    ftype = c("BB"),
                    chainfrags_sn1 = c("mg_M+H-H2O"),
                    chainfrags_sn2 = c("mg_M+H-H2O"),
                    intrules = c("mg_sn1/mg_sn2"),
                    rates = c("1/2"),
                    intrequired = c(F),
                    coelCutoff = 0.8,
                    dbs, 
                    verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "PG",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("PG" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious PG annotations removed")}
      msobject$annotation$detailsAnnotation$PG <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$pgdb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=2,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="PG",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PG <- list()
    msobject$annotation$detailsAnnotation$PG$candidates <- candidates
    msobject$annotation$detailsAnnotation$PG$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$PG$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$PG$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PG <- list()
  }
  return(msobject)
}

# idPIpos
#' Phosphoinositols (PI) annotation for ESI+
#'
#' PI identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for PE in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param rates character vector with the expected rates between fragments given
#' as a string (e.g. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See
#' \link{createLipidDB} and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idPIpos} function involves 5 steps. 1) FullMS-based
#' identification of candidate PI as M+H, M+NH4 and M+Na. 2) Search of PI class
#' fragments: loss of head group (DG as M+H-H2O) coeluting with the precursor
#' ion. 3) Search of specific fragments that inform about chain composition at
#' sn1 (MG as M+H-H2O or LPI as M+H-H2O resulting from the loss of the FA chain 
#' at sn2) and sn2 (MG as M+H-H2O or LPI as M+H-H2O resulting from the loss of 
#' the FA chain at sn1). 4) Look for possible chains structure based on the 
#' combination of chain fragments. 5) Check intensity rules to confirm chains 
#' position. MG or LPI from sn1 are at least twice more intense than the ones 
#' from sn2.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idPIpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idPIpos <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("M+H", "M+NH4", "M+Na"),
                    clfrags = c("dg_M+H-H2O"),
                    clrequired = c(F),
                    ftype = c("BB"),
                    chainfrags_sn1 = c("mg_M+H-H2O", "lysopi_M+H-H2O"),
                    chainfrags_sn2 = c("mg_M+H-H2O", "lysopi_M+H-H2O"),
                    intrules = c("mg_sn1/mg_sn2", "lysopi_sn1/lysopi_sn2"),
                    rates = c("2/1", "2/1"),
                    intrequired = c(F, F),
                    coelCutoff = 0.8,
                    dbs, 
                    verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "PI",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("PI" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious PI annotations removed")}
      msobject$annotation$detailsAnnotation$PI <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$pidb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=2,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="PI",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PI <- list()
    msobject$annotation$detailsAnnotation$PI$candidates <- candidates
    msobject$annotation$detailsAnnotation$PI$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$PI$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$PI$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$PI <- list()
  }
  return(msobject)
}

# idSphpos
#' Sphingoid bases (Sph) annotation for ESI-
#'
#' Sph identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursors and product
#' ions. By default, 3 seconds.
#' @param rt rt window where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for Sph in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idSphpos} function involves 2 steps. 1) FullMS-based
#' identification of candidate Sph as M+H. 2) Search of Sph class fragments:
#' neutral loss of 1 or 2 H2O molecules.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (in this
#' case, as Sph only have one chain, only Subclass and FA level are possible)
#' and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idSphpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idSphpos <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M+H"),
                     clfrags = c("sph_M+H-H2O", "sph_M+H-2H2O"),
                     clrequired = c(F, F),
                     ftype = c("BB", "BB"),
                     coelCutoff = 0.8,
                     dbs,
                     verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "Sph",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("Sph" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious Sph annotations removed")}
      msobject$annotation$detailsAnnotation$Sph <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$sphdb, ppm = ppm_precursor,
                               rt = rt, adducts = adducts, rttol = rttol,
                               dbs = dbs, rawData = rawData,
                               coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb = list(),
                           intrules  = c(), intConf = list(), nchains = 0,
                           class="Sph",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$Sph <- list()
    msobject$annotation$detailsAnnotation$Sph$candidates <- candidates
    msobject$annotation$detailsAnnotation$Sph$classfragments <- classConf$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$Sph$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$Sph <- list()
  }
  return(msobject)
}

# idSphPpos
#' Sphingoid bases phosphate (SphP) annotation for ESI+
#'
#' SphP identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursors and product
#' ions. By default, 3 seconds.
#' @param rt rt window where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for Sph in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idSphPpos} function involves 2 steps. 1) FullMS-based
#' identification of candidate SphP as M+H. 2) Search of SphP class fragments:
#' neutral loss of 1 or 2 H2O molecules, or H2O and NH4.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (in this
#' case, as SphP only have one chain, only Subclass and FA level are possible).
#' and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idSphPpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idSphPpos <- function(msobject,
                      ppm_precursor = 5,
                      ppm_products = 10,
                      rttol = 3,
                      rt,
                      adducts = c("M+H"),
                      clfrags = c("sphP_M+H-H2O", "sphP_M+H-2H2O", "sphP_M+H-H2O-NH4"),
                      clrequired = c(F, F, F),
                      ftype = c("BB", "BB", "BB"),
                      coelCutoff = 0.7,
                      dbs,
                      verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "SphP",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("SphP" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious SphP annotations removed")}
      msobject$annotation$detailsAnnotation$SphP <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$sphPdb, ppm = ppm_precursor,
                               rt = rt, adducts = adducts, rttol = rttol,
                               dbs = dbs, rawData = rawData,
                               coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # prepare output
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb = list(),
                           intrules  = c(), intConf = list(), nchains = 0,
                           class="SphP",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$SphP <- list()
    msobject$annotation$detailsAnnotation$SphP$candidates <- candidates
    msobject$annotation$detailsAnnotation$SphP$classfragments <- classConf$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$SphP$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$SphP <- list()
  }
  return(msobject)
}

# idCerpos
#' Ceramides (Cer) annotation for ESI+
#'
#' Ceramides identification based on fragmentation patterns for LC-MS/MS DIA or
#' DDA data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for Cer in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param rates character vector with the expected rates between fragments given
#' as a string (e.g. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between peaks (adducts, parent and
#' fragment ions...). Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idCerpos} function involves 5 steps. 1) FullMS-based
#' identification of candidate Cer as M+H, M+H-H2O and M+Na. 2) Search of Cer
#' class fragments: there isn't any class fragment by default. 3) Search of
#' specific fragments that inform about the sphingoid base (Sph as M+H-2H2O
#' resulting from the loss of the FA chain) and the FA chain (by default it is
#' calculated using the difference between precursor and sph fragments).
#' 4) Look for possible chains structure based on the combination of chain
#' fragments. 5) Check intensity rules to confirm chains position. In this case,
#' there are no intensity rules by default.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idCerpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idCerpos <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M+H-H2O", "M+Na", "M+H"),
                     clfrags = c(),
                     clrequired = c(),
                     ftype = c(),
                     chainfrags_sn1 = c("sph_M+H-2H2O"),
                     chainfrags_sn2 = c(""),
                     intrules = c(),
                     rates = c(),
                     intrequired = c(),
                     coelCutoff = 0.8,
                     dbs,
                     verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "Cer",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("Cer" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious Ceramide annotations removed")}
      msobject$annotation$detailsAnnotation$Cer <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, db = dbs$cerdb,
                               ppm = ppm_precursor, rt = rt, adducts = adducts,
                               rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=2,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="Cer",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$Cer <- list()
    msobject$annotation$detailsAnnotation$Cer$candidates <- candidates
    msobject$annotation$detailsAnnotation$Cer$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$Cer$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$Cer$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$Cer <- list()
  }
  return(msobject)
}

# idCerPpos
#' Ceramides phosphate (CerP) annotation for ESI+
#'
#' CerP identification based on fragmentation patterns for LC-MS/MS DIA or
#' DDA data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for Cer in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param rates character vector with the expected rates between fragments given
#' as a string (e.g. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between peaks (adducts, parent and
#' fragment ions...). Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idCerPpos} function involves 5 steps. 1) FullMS-based
#' identification of candidate CerP as M+H. 2) Search of Cer class fragments: 
#' Cer as M+H-H2O and M+H-2H2O resulting from the loss of the phosphate group and
#' 1 or 2 H2O molecules. 3) Search of specific fragments that inform about the 
#' sphingoid base (Sph as M+H-2H2O resulting from the loss of the FA chain and 
#' the phosphate group) and the FA chain (by default it is calculated using the 
#' difference between precursor and sph fragments). 4) Look for possible chains 
#' structure based on the combination of chain fragments. 5) Check intensity 
#' rules to confirm chains position. In this case, there are no intensity rules 
#' by default.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idCerPpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idCerPpos <- function(msobject,
                      ppm_precursor = 5,
                      ppm_products = 10,
                      rttol = 3,
                      rt,
                      adducts = c("M+H"),
                      clfrags = c("cer_M+H-H2O", "cer_M+H-2H2O"),
                      clrequired = c(F, F),
                      ftype = c("BB", "BB"),
                      chainfrags_sn1 = c("sph_M+H-2H2O"),
                      chainfrags_sn2 = c(""),
                      intrules = c(),
                      rates = c(),
                      intrequired = c(),
                      coelCutoff = 0.8,
                      dbs,
                      verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "CerP",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("CerP" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious CerP annotations removed")}
      msobject$annotation$detailsAnnotation$CerP <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, db = dbs$cerPdb,
                               ppm = ppm_precursor, rt = rt, adducts = adducts,
                               rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=2,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="CerP",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$CerP <- list()
    msobject$annotation$detailsAnnotation$CerP$candidates <- candidates
    msobject$annotation$detailsAnnotation$CerP$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$CerP$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$CerP$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$CerP <- list()
  }
  return(msobject)
}

# idSMpos
#' Sphyngomyelines (SM) annotation for ESI+
#'
#' SM identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for SM in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param rates character vector with the expected rates between fragments given
#' as a string (e.g. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idSMpos} function involves 5 steps. 1) FullMS-based
#' identification of candidate SM as M+H and M+Na. 2) Search of SM class
#' fragments: 104.1075, 184.0739 and neutral loss of 183.06604 coeluting with
#' the precursor ion. 3) Search of specific fragments that inform about the
#' composition of the sphingoid base (Sph as M+H-2H2O resulting from the loss of
#' the FA chain) and the FA chain (by default it is calculated using the
#' difference between precursor and sph chain fragments). 4) Look for possible
#' chains structure based on the combination of chain fragments. 5) Check
#' intensity rules to confirm chains position. In this case, there are no
#' intensity rules by default as FA chain is unlikely to be detected.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idSMpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idSMpos <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("M+H", "M+Na"),
                    clfrags = c(104.1075, 184.0739, 183.06604),
                    clrequired = c(F, F, F),
                    ftype = c("F", "F", "NL"),
                    chainfrags_sn1 = c("sph_M+H-2H2O"),
                    chainfrags_sn2 = c(""),
                    intrules = c(),
                    rates = c(),
                    intrequired = c(),
                    coelCutoff = 0.8,
                    dbs,
                    verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "SM",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("SM" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious SM annotations removed")}
      msobject$annotation$detailsAnnotation$SM <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$smdb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=2,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="SM",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$SM <- list()
    msobject$annotation$detailsAnnotation$SM$candidates <- candidates
    msobject$annotation$detailsAnnotation$SM$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$SM$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$SM$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$SM <- list()
  }
  return(msobject)
}

# idCarpos
#' Acylcarnitine annotation for ESI+
#'
#' Acylcarnitines identification based on fragmentation patterns for LC-MS/MS
#' DIA or DDA data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for Carnitines in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments. See \link{chainFrags} for details.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idCarpos} function involves 3 steps. 1) FullMS-based
#' identification of candidate carnitines as M+H and M+Na. 2) Search of
#' carnitine class fragments: 60.0807 and 85.0295 or its loss (FA as M+H-H20)
#' coeluting with the precursor ion. 3) Search of specific fragments coming from
#' the FA chain (FA as M+H-H2O).
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (in this
#' case, as Carnitines only have one chain, only Subclass and FA level are
#' possible) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idCarpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idCarpos <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M+H", "M+Na"),
                     clfrags = c(60.0807, 85.0295, "fa_M+H-H2O"),
                     clrequired = c(F, F, F),
                     ftype = c("F", "F", "BB"),
                     chainfrags_sn1 = c("fa_M+H-H2O"),
                     coelCutoff = 0.8,
                     dbs,
                     verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "Carnitine",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("Carnitine" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious Carnitine annotations removed")}
      msobject$annotation$detailsAnnotation$Carnitine <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$carnitinedb, ppm = ppm_precursor,
                               rt = rt, adducts = adducts, rttol = rttol,
                               dbs = dbs, rawData = rawData,
                               coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=1, sn1)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules = c(), rates = c(), intrequired = c(),
                                   nchains=1, chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb,
                           intrules  = c(), intConf, nchains = 1,
                           class="Carnitine",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$Carnitine <- list()
    msobject$annotation$detailsAnnotation$Carnitine$candidates <- candidates
    msobject$annotation$detailsAnnotation$Carnitine$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$Carnitine$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$Carnitine$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$Carnitine <- list()
  }
  return(msobject)
}

# idCEpos
#' Cholesteryl Esters (CE) annotation for ESI+
#'
#' CE identification based on fragmentation patterns for LC-MS/MS DIA or DDA data
#' acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for CE in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments. See \link{chainFrags} for details.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idCEpos} function involves 3 steps. 1) FullMS-based
#' identification of candidate CE as 2M+NH4, 2M+Na, M+NH4 and M+Na. 2) Search of
#' CE class fragments: 369.3516 or its loss (FA as M+H-H20) coeluting with the
#' precursor ion. 3) Search of specific fragments that confirm chain composition
#' (FA as M+H-H2O).
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (in this
#' case, as CE only have one chain, only Subclass and FA level are possible)
#' and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idCEpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idCEpos <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("2M+NH4", "2M+Na", "M+NH4", "M+Na"),
                    clfrags = c(369.3516, "fa_M+H-H2O"),
                    clrequired = c(F, F),
                    ftype = c("F", "BB"),
                    chainfrags_sn1 = c("fa_M+H-H2O"),
                    coelCutoff = 0.8,
                    dbs,
                    verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "CE",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("CE" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious CE annotations removed")}
      msobject$annotation$detailsAnnotation$CE <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$CEdb, ppm = ppm_precursor,
                               rt = rt, adducts = adducts, rttol = rttol,
                               dbs = dbs, rawData = rawData,
                               coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=1, sn1)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules = c(), rates = c(),
                                   intrequired = c(), nchains=1,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb,
                           intrules = c(),
                           intConf, nchains = 1, class="CE",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$CE <- list()
    msobject$annotation$detailsAnnotation$CE$candidates <- candidates
    msobject$annotation$detailsAnnotation$CE$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$CE$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$CE$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$CE <- list()
  }
  return(msobject)
}

# idDGpos
#' Diacylglycerols (DG) annotation for ESI+
#'
#' DG identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for DG in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param rates character vector with the expected rates between fragments given
#' as a string (e.g. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idDGpos} function involves 5 steps. 1) FullMS-based
#' identification of candidate DG as M+H-H2O, M+NH4 and M+Na. 2) Search of DG
#' class fragments: there are no class fragment by default. 3) Search of
#' specific fragments that inform about the FA chains (MGs as M+H-H2O
#' resulting from the loss of the FA chains). 4) Look for possible chains
#' structure based on the combination of chain fragments. 5) Check intensity
#' rules to confirm chains position: MG coming from the loss of the sn2 chain is
#' more intense than the one coming from the loss of sn1.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idDGpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idDGpos <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("M+H-H2O", "M+NH4", "M+Na"),
                    clfrags = c(),
                    clrequired = c(),
                    ftype = c(),
                    chainfrags_sn1 = c("mg_M+H-H2O"),
                    chainfrags_sn2 = c("mg_M+H-H2O"),
                    intrules = c("mg_sn1/mg_sn2"),
                    rates = c("1"),
                    intrequired = c(T),
                    coelCutoff = 0.8,
                    dbs,
                    verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "DG",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("DG" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious DG annotations removed")}
      msobject$annotation$detailsAnnotation$DG <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$dgdb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=2,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="DG",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$DG <- list()
    msobject$annotation$detailsAnnotation$DG$candidates <- candidates
    msobject$annotation$detailsAnnotation$DG$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$DG$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$DG$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$DG <- list()
  }
  return(msobject)
}

# idTGpos
#' Triacylglycerols (TG) annotation for ESI+
#'
#' TG identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for TG in ESI+. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param chainfrags_sn3 character vector containing the fragmentation rules for
#' the chain fragments in sn3 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn2 chains.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}. If some intensity rules should be employed to
#' identify the chains position but they are't known yet, use "Unknown". If it
#' isn't required, leave an empty vector.
#' @param rates character vector with the expected rates between fragments given
#' as a string (e.g. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idTGpos} function involves 5 steps. 1) FullMS-based
#' identification of candidate TG as M+NH4 and M+Na. 2) Search of TG class
#' fragments: there are no class fragment by default. 3) Search of specific
#' fragments that inform about the FA chains: DGs resulting from the loss of FA
#' chains as M+H-H2O.  4) Look for possible chains structure based on the
#' combination of chain fragments. 5) Check intensity rules to confirm chains
#' position. In the case of TG, DG resulting from the loss of sn2 if the most
#' intense, followed by the loss of sn1 and sn3, but this FA position level
#' still needs to be improved due to the high level of coelution for TG.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been written based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idTGpos(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idTGpos <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("M+NH4", "M+Na"),
                    clfrags = c(),
                    clrequired = c(),
                    ftype = c(),
                    chainfrags_sn1 = c("cbdiff-dg_M+H-H2O"),
                    chainfrags_sn2 = c("cbdiff-dg_M+H-H2O"),
                    chainfrags_sn3 = c("cbdiff-dg_M+H-H2O"),
                    intrules = c("cbdiff-dg_sn2/cbdiff-dg_sn1",
                                 "cbdiff-dg_sn2/cbdiff-dg_sn3",
                                 "cbdiff-dg_sn1/cbdiff-dg_sn3"),
                    rates = c("1", "1", "1"),
                    intrequired = c(T, T, T),
                    coelCutoff = 0.8,
                    dbs,
                    verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "TG",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("TG" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious TG annotations removed")}
      msobject$annotation$detailsAnnotation$TG <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$tgdb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products,
                      candidates = candidates, sn1, dbs)
    sn3 <- chainFrags(coelfrags, chainfrags_sn3, ppm_products, candidates, sn1,
                      dbs)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=3, sn1, sn2, sn3)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=3,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb,
                           intrules, intConf, nchains = 3, class="TG",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$TG <- list()
    msobject$annotation$detailsAnnotation$TG$candidates <- candidates
    msobject$annotation$detailsAnnotation$TG$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$TG$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$TG$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$TG <- list()
  }
  return(msobject)
}

# idAcylCerpos
#' Acylceramides (AcylCer) annotation for ESI+
#'
#' AcylCer identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in positive mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for Cer in ESI-. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the sphingoid base. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param chainfrags_sn3 character vector containing the fragmentation rules for
#' the acyl chain. See \link{chainFrags} for details.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param rates character vector with the expected ratesbetween fragments given
#' as a string (e.g. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param verbose print information messages.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, lipid class, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), mz, RT (in seconds), I (intensity), Adducts, ppm (mz error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and Score (parent-fragment coelution 
#' score mean in DIA data or relative sum intensity in DDA of all fragments used 
#' for the identification).
#'
#' @details \code{idAcylCerpos} function involves 5 steps. 1) FullMS-based
#' identification of candidate AcylCer as M+H, M+H-H2O and M+Na. 2) Search of 
#' AcylCer class fragments: there are no class fragments by default. 3) Search of 
#' specific fragments that inform about the acyl chain (Cer as M+H, M+H-H2O or
#' M+H-2H2H), the sphingoid base (Sph as M+H-H2O or M+H-2H2O) and the FA chain 
#' (FA as M+H but with a N intead of an O, what results in a mass difference of 
#' 0.02329 with the Mn of the FA chain). 4) Look for possible chains structure 
#' based on the combination of chain fragments. 5) Check intensity rules to 
#' confirm chains position. In this case, Sph fragment must be twice more
#' intense than the loss of the acyl chain and at least 5 times more intense than
#' the FA chain from sn3.
#'
#' Results data frame shows: ID, lipid class, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (mz error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and Score (parent-fragment coelution score mean in DIA data or relative 
#' sum intensity in DDA of all fragments used for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Synapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' msobject <- idCerPneg(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idAcylCerpos <- function(msobject,
                         ppm_precursor = 5,
                         ppm_products = 10,
                         rttol = 3,
                         rt,
                         adducts = c("M+H", "M+H-H2O", "M+Na"),
                         clfrags = c(),
                         clrequired = c(),
                         ftype = c(),
                         chainfrags_sn1 = c("cbdiff-cer_M+H", "cbdiff-cer_M+H-H2O", "cbdiff-cer_M+H-2H2O"),
                         chainfrags_sn2 = c("sph_M+H-H2O", "sph_M+H-2H2O"),
                         chainfrags_sn3 = c("fa_Mn+0.02329"),
                         intrules = c("sph_sn2/cbdiff-cer_sn1", "sph_sn2/fa_sn3"),
                         rates = c("2/1", "5/1"),
                         intrequired = c(T, T),
                         coelCutoff = 0.8,
                         dbs,
                         verbose = TRUE){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "positive"){
    stop("Data wasn't acquired in positive mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "rawData", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!all(c("MS1", "MS2") %in% names(msobject$rawData))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the adductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "isoGroup")]
  # Peaklist MS2: 
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    MS2 <- msobject$rawData$MS2[,c("mz", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("mz", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$rawData$MS1, msobject$rawData$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("RT", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject$annotation)){
    if (nrow(msobject$annotation$results) > 0){
      msobject$annotation$results <- msobject$annotation$results[msobject$annotation$results$Class != "AcylCer",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject$annotation)){
    if("AcylCer" %in% names(msobject$annotation$detailsAnnotation)){
      if(verbose){cat("\nPrevious AcylCer annotations removed")}
      msobject$annotation$detailsAnnotation$AcylCer <- list()
    }
  }
  ##############################################################################
  # set rt limits
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }
  ##############################################################################
  # Start identification steps
  
  # candidates search
  candidates <- findCandidates(MS1, dbs$acylcerdb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)
  
  if (nrow(candidates) > 0){
    if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }
    
    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)
    
    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    sn3 <- chainFrags(coelfrags, chainfrags_sn3, ppm_products, candidates, sn1,
                      dbs)
    
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=3, sn1, sn2, sn3)
    
    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=2,
                                   chainsComb)
    
    # prepare output
    res <- organizeResults(candidates, coelfrags, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 3, class="AcylCer",
                           acquisitionmode = msobject$metaData$generalMetadata$acquisitionmode)
    
    # update msobject
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$AcylCer <- list()
    msobject$annotation$detailsAnnotation$AcylCer$candidates <- candidates
    msobject$annotation$detailsAnnotation$AcylCer$classfragments <- classConf$fragments
    msobject$annotation$detailsAnnotation$AcylCer$chainfragments <- chainsComb$fragments
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      msobject$annotation$detailsAnnotation$AcylCer$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject$annotation)){
      msobject$annotation$results <- rbind(msobject$annotation$results, res)
    } else {
      msobject$annotation$results <- res
    }
    msobject$annotation$detailsAnnotation$AcylCer <- list()
  }
  return(msobject)
}


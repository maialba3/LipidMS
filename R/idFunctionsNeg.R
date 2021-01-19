# idNEG
#' Lipids annotation for ESI-
#'
#' Lipids annotation based on fragmentation patterns for LC-MS/MS DIA or DDA data
#' acquired in negative mode. This function compiles all functions writen for
#' ESI- annotations.
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idNEG(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idNEG <- function(msobject,
                  ppm_precursor = 5,
                  ppm_products = 10,
                  rttol = 5,
                  coelCutoff = 0.8,
                  lipidClasses = c("FA", "FAHFA", "LPC", "LPE", "LPG", "LPI",
                                   "LPS", "PC", "PE", "PG", "PI", "PS", "Sph",
                                   "SphP", "Cer", "CL", "BA"),
                  dbs){

  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if ("results" %in% names(msobject)){
    cat("\n Removing previous results...")
    msobject$results <- NULL
    msobject$detailsAnnotation <- NULL
    msobject$annotatedPeaklist <- NULL
    cat("OK")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(lipidClasses %in% c("FA", "FAHFA", "LPC", "LPE", "LPG", "LPI",
                               "LPS", "PC", "PE", "PG", "PI", "PS", "Sph",
                               "SphP", "Cer", "CL", "BA"))){
    stop("Lipid classes allowed for negative annotation are: FA, FAHFA, LPC, LPE,
         LPG, LPI, LPS, PC, PE, PG, PI, PS, Sph, SphP, Cer, CL and BA")
  }

  cat("\n Starting annotation...")
  if ("FA" %in% lipidClasses){
    cat("\n  Searching for FA...")
    msobject <-  idFAneg(msobject = msobject, ppm_precursor= ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("FAHFA" %in% lipidClasses){
    cat("\n  Searching for FAHFA...")
    msobject <-  idFAHFAneg(msobject = msobject, ppm_precursor= ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("LPC" %in% lipidClasses){
    cat("\n  Searching for LPC...")
    msobject <-  idLPCneg(msobject = msobject, ppm_precursor= ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("LPE" %in% lipidClasses){
    cat("\n  Searching for LPE...")
    msobject <-  idLPEneg(msobject = msobject, ppm_precursor= ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("LPG" %in% lipidClasses){
    cat("\n  Searching for LPG...")
    msobject <-  idLPGneg(msobject = msobject, ppm_precursor= ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("LPI" %in% lipidClasses){
    cat("\n  Searching for LPI...")
    msobject <-  idLPIneg(msobject = msobject, ppm_precursor= ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("LPS" %in% lipidClasses){
    cat("\n  Searching for LPS...")
    msobject <-  idLPSneg(msobject = msobject, ppm_precursor= ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("PC" %in% lipidClasses){
    cat("\n  Searching for PC...")
    msobject <-  idPCneg(msobject = msobject, ppm_precursor= ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("PE" %in% lipidClasses){
    cat("\n  Searching for PE...")
    msobject <-  idPEneg(msobject = msobject, ppm_precursor= ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("PG" %in% lipidClasses){
    cat("\n  Searching for PG...")
    msobject <-  idPGneg(msobject = msobject, ppm_precursor= ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("PI" %in% lipidClasses){
    cat("\n  Searching for PI...")
    msobject <-  idPIneg(msobject = msobject, ppm_precursor= ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("PS" %in% lipidClasses){
    cat("\n  Searching for PS...")
    msobject <-  idPSneg(msobject = msobject, ppm_precursor= ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("Sph" %in% lipidClasses){
    cat("\n  Searching for Sph...")
    msobject <-  idSphneg(msobject = msobject, ppm_precursor= ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("SphP" %in% lipidClasses){
    cat("\n  Searching for SphP...")
    msobject <-  idSphPneg(msobject = msobject, ppm_precursor= ppm_precursor,
                           ppm_products = ppm_products, rttol = rttol,
                           coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("Cer" %in% lipidClasses){
    cat("\n  Searching for Cer...")
    msobject <-  idCerneg(msobject = msobject, ppm_precursor= ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("CL" %in% lipidClasses){
    cat("\n  Searching for CL...")
    msobject <-  idCLneg(msobject = msobject, ppm_precursor= ppm_precursor,
                         ppm_products = ppm_products, rttol = rttol,
                         coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }
  if ("BA" %in% lipidClasses){
    cat("\n  Searching for Bile acids...")
    msobject <-  idBAneg(msobject = msobject, ppm_precursor= ppm_precursor,
                          ppm_products = ppm_products, rttol = rttol,
                          coelCutoff = coelCutoff, dbs = dbs)
    cat("OK")
  }

  cat("\n Preparing output...")
  if (nrow(msobject$results) > 0){
    annotatedPeaklist <- crossTables(msobject$peaklist$MS1, msobject$results,
                                     ppm = ppm_precursor, rttol = rttol,
                                     dbs = dbs)
  } else {
    annotatedPeaklist <- "No results were found"
  }
  msobject$annotatedPeaklist <- annotatedPeaklist
  cat("OK\n")
  return(msobject)
}

# idFAneg
#' Fatty Acids (FA) annotation for ESI-
#'
#' FA identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for FA in ESI-. Adducts allowed can
#' be modified in addutcsTable (dbs argument).
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idFAneg} function involves 2 steps. 1) FullMS-based
#' identification of candidate FA as M-H or 2M-H. 2) Search of FA class
#' fragments: neutral loss of H2O coeluting with the precursor ion or the
#' molecular ion.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (in this
#' case, just MS-only or Subclass level (if any class fragment is defined) are
#' possible) and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idFAneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
idFAneg <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("M-H", "2M-H"),
                    clfrags = c("fa_M-H", "fa_M-H-H2O"),
                    clrequired = c(FALSE, FALSE),
                    ftype = c("BB", "BB"),
                    coelCutoff = 0.8,
                    dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "FA",]
    }
  }
  if ("lipidClasses" %in% names(msobject)){
    if("FA" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious FA annotations removed")
      msobject$detailsAnnotation$FA <- list()
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
  candidates <- findCandidates(MS1, dbs$fadb, ppm = ppm_precursor,
                               rt = rt, adducts = adducts, rttol = rttol,
                               dbs = dbs, rawData = rawData,
                               coelCutoff = coelCutoff)

  if (nrow(candidates) > 0){
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }

    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)

    # prepare output
    res <- organizeResults(candidates, clfrags, classConf, chainsComb = list(),
                           intrules  = c(), intConf = list(), nchains = 0,
                           class="FA")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$FA <- list()
    msobject$detailsAnnotation$FA$candidates <- candidates
    msobject$detailsAnnotation$FA$classfragments <- classConf$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$FA$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$FA <- list()
  }
  return(msobject)
}

# idFAHFAneg
#' FAHFA annotation for ESI-
#'
#' FAHFA identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for FAHFA in ESI-. Adducts allowed can
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idFAHFAneg} function involves 5 steps. 1) FullMS-based
#' identification of candidate FAHFA as M-H. 2) Search of FAHFA class fragments:
#' there is't any class fragment by default. 3) Search of specific fragments
#' that inform about chain composition in sn1 (HFA as M-H resulting from the
#' loss of the FA chain) and sn2 (FA chain as M-H). 4) Look for possible
#' chains structure based on the combination of chain fragments. 5) Check
#' intensity rules to confirm chains position. In this case, HFA intensity has
#' to be higher than FA.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idFAHFAneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idFAHFAneg <- function(msobject,
                       ppm_precursor = 5,
                       ppm_products = 10,
                       rttol = 3,
                       rt,
                       adducts = c("M-H"),
                       clfrags = c(),
                       clrequired = c(),
                       ftype = c(),
                       chainfrags_sn1 = c("hfa_M-H"),
                       chainfrags_sn2 = c("fa_M-H"),
                       intrules = c("hfa_sn1/fa_sn2"),
                       rates = c("3/1"),
                       intrequired = c(T),
                       coelCutoff = 0.8,
                       dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous FAHFA annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "FAHFA",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("FAHFA" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious FAHFA annotations removed")
      msobject$detailsAnnotation$FAHFA <- list()
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
  candidates <- findCandidates(MS1, dbs$fahfadb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)

  if (nrow(candidates) > 0){
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
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
    res <- organizeResults(candidates, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="FAHFA")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$FAHFA <- list()
    msobject$detailsAnnotation$FAHFA$candidates <- candidates
    msobject$detailsAnnotation$FAHFA$classfragments <- classConf$fragments
    msobject$detailsAnnotation$FAHFA$chainfragments <- chainsComb$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$FAHFA$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$FAHFA <- list()
  }
  return(msobject)
}

# idLPCneg
#' Lysophosphocholines (LPC) annotation for ESI-
#'
#' LPC identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for LPC in ESI-. Adducts allowed can
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idLPCneg} function involves 3 steps. 1) FullMS-based
#' identification of candidate LPC as M+CH3COO, M-CH3 and M+CH3COO-CH3. To avoid
#' incorrect annotations of PE as PC, candidates which are present just as M-CH3
#' will be ignored. 2) Search of LPC class fragments: 168.0426, 224.0688, lysoPA
#' as M-H or lysoPC as M-CH3 coeluting with the precursor ion. 3) Search of
#' specific fragments that confirm chain composition (FA as M-H).
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (in this
#' case, as LPC only have one chain, only Subclass and FA level are possible)
#' and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idLPCneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idLPCneg <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M+CH3COO", "M-CH3", "M+CH3COO-CH3"),
                     clfrags = c(168.0426, 224.0688, "lysopa_M-H", "lysopc_M-CH3"),
                     clrequired = c(F, F, F, F),
                     ftype = c("F", "F", "BB", "BB"),
                     chainfrags_sn1 = c("fa_M-H"),
                     coelCutoff = 0.8,
                     dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous LPC annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "LPC",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("LPC" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious LPC annotations removed")
      msobject$detailsAnnotation$LPC <- list()
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
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
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
    res <- organizeResults(candidates, clfrags, classConf, chainsComb,
                           intrules  = c(), intConf, nchains = 1, class="LPC")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$LPC <- list()
    msobject$detailsAnnotation$LPC$candidates <- candidates
    msobject$detailsAnnotation$LPC$classfragments <- classConf$fragments
    msobject$detailsAnnotation$LPC$chainfragments <- chainsComb$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$LPC$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$LPC <- list()
  }
  return(msobject)
}

# idLPEneg
#' Lysophosphoethanolamines (LPE) annotation for ESI-
#'
#' LPE identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for LPE in ESI-. Adducts allowed can
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idLPEneg} function involves 3 steps. 1) FullMS-based
#' identification of candidate LPE as M-H. 2) Search of
#' LPE class fragments: 140.0115, 196.038 and 214.048 coeluting with the
#' precursor ion. If a loss of CH3 group is found coeluting with any candidate,
#' this will be excluded as it is a characteristic fragment of LPC.3) Search of
#' specific fragments that confirm chain composition (FA as M-H).
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (in this
#' case, as LPE only have one chain, only Subclass and FA level are possible)
#' and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idLPEneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idLPEneg <- function(msobject, ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M-H"),
                     clfrags = c(140.0115, 196.038, 214.048, "lysope_M-CH3"),
                     clrequired = c(F, F, F, "excluding"),
                     ftype = c("F", "F", "F", "BB"),
                     chainfrags_sn1 = c("fa_M-H"),
                     coelCutoff = 0.8,
                     dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "LPE",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("LPE" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious LPE annotations removed")
      msobject$detailsAnnotation$LPE <- list()
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
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
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
    res <- organizeResults(candidates, clfrags, classConf, chainsComb,
                           intrules  = c(), intConf, nchains = 1, class="LPE")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$LPE <- list()
    msobject$detailsAnnotation$LPE$candidates <- candidates
    msobject$detailsAnnotation$LPE$classfragments <- classConf$fragments
    msobject$detailsAnnotation$LPE$chainfragments <- chainsComb$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$LPE$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$LPE <- list()
  }
  return(msobject)
}

# idLPGneg
#' Lysophosphoglycerols (LPG) annotation for ESI-
#'
#' LPG identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for LPG in ESI-. Adducts allowed can
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idLPGneg} function involves 3 steps. 1) FullMS-based
#' identification of candidate LPG as M-H. 2) Search of LPG class fragments:
#' 152.9958, 227.0326, 209.022 and neutral loss of 74.0359 coeluting with the
#' precursor ion. 3) Search of specific fragments that confirm chain composition
#' (FA as M-H).
#'
#' Results data frame shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (in this
#' case, as LPG only have one chain, only Subclass and FA level are possible)
#' and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idLPGneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idLPGneg <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M-H"),
                     clfrags = c(152.9958, 227.0326, 209.022, 74.0359),
                     clrequired = c(F, F, F, F),
                     ftype = c("F", "F", "F", "NL"),
                     chainfrags_sn1 = c("fa_M-H"),
                     coelCutoff = 0.8,
                     dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "LPG",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("LPG" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious LPG annotations removed")
      msobject$detailsAnnotation$LPG <- list()
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
  candidates <- findCandidates(MS1, dbs$lysopgdb, ppm = ppm_precursor,
                               rt = rt, adducts = adducts, rttol = rttol,
                               dbs = dbs, rawData = rawData,
                               coelCutoff = coelCutoff)

  if (nrow(candidates) > 0){
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
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
    res <- organizeResults(candidates, clfrags, classConf, chainsComb,
                           intrules  = c(), intConf, nchains = 1, class="LPG")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$LPG <- list()
    msobject$detailsAnnotation$LPG$candidates <- candidates
    msobject$detailsAnnotation$LPG$classfragments <- classConf$fragments
    msobject$detailsAnnotation$LPG$chainfragments <- chainsComb$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$LPG$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$LPG <- list()
  }
  return(msobject)
}

# idLPIneg
#' Lysophosphoinositols (LPI) annotation for ESI-
#'
#' LPI identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for LPI in ESI-. Adducts allowed can
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idLPIneg} function involves 3 steps. 1) FullMS-based
#' identification of candidate LPI as M-H. 2) Search of
#' LPI class fragments: 241.0115, 223.0008, 259.0219 and 297.0375 coeluting
#' with the precursor ion. 3) Search of specific fragments that confirm chain
#' composition (FA as M-H).
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (in this
#' case, as LPI only have one chain, only Subclass and FA level are possible)
#' and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idLPIneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idLPIneg <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M-H"),
                     clfrags = c(241.0115, 223.0008, 259.0219, 297.0375),
                     clrequired = c(F, F, F, F),
                     ftype = c("F", "F", "F", "F"),
                     chainfrags_sn1 = c("fa_M-H"),
                     coelCutoff = 0.8,
                     dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "LPI",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("LPI" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious LPI annotations removed")
      msobject$detailsAnnotation$LPI <- list()
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
  candidates <- findCandidates(MS1, dbs$lysopidb, ppm = ppm_precursor,
                               rt = rt, adducts = adducts, rttol = rttol,
                               dbs = dbs, rawData = rawData,
                               coelCutoff = coelCutoff)

  if (nrow(candidates) > 0){
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
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
    res <- organizeResults(candidates, clfrags, classConf, chainsComb,
                           intrules  = c(), intConf, nchains = 1, class="LPI")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$LPI <- list()
    msobject$detailsAnnotation$LPI$candidates <- candidates
    msobject$detailsAnnotation$LPI$classfragments <- classConf$fragments
    msobject$detailsAnnotation$LPI$chainfragments <- chainsComb$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$LPI$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$LPI <- list()
  }
  return(msobject)
}

# idLPSneg
#' Lysophosphoserines (LPS) annotation for ESI-
#'
#' LPS identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for LPS in ESI-. Adducts allowed can
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idLPSneg} function involves 3 steps. 1) FullMS-based
#' identification of candidate LPS as M-H and M+Na-2H. 2) Search of
#' LPS class fragments: neutral loss of 87.032 coeluting with the precursor ion.
#' 3) Search of specific fragments that confirm chain composition (FA as M-H).
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (in this
#' case, as LPS only have one chain, only Subclass and FA level are possible)
#' and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idLPSneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idLPSneg <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M-H", "M+Na-2H"),
                     clfrags = c(87.032),
                     clrequired = c(F),
                     ftype = c("NL"),
                     chainfrags_sn1 = c("fa_M-H"),
                     coelCutoff = 0.8,
                     dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "LPS",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("LPS" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious LPS annotations removed")
      msobject$detailsAnnotation$LPS <- list()
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
  candidates <- findCandidates(MS1, dbs$lysopsdb, ppm = ppm_precursor,
                               rt = rt, adducts = adducts, rttol = rttol,
                               dbs = dbs, rawData = rawData,
                               coelCutoff = coelCutoff)

  if (nrow(candidates) > 0){
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
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
    res <- organizeResults(candidates, clfrags, classConf, chainsComb,
                           intrules  = c(), intConf, nchains = 1, class="LPS")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$LPS <- list()
    msobject$detailsAnnotation$LPS$candidates <- candidates
    msobject$detailsAnnotation$LPS$classfragments <- classConf$fragments
    msobject$detailsAnnotation$LPS$chainfragments <- chainsComb$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$LPS$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$LPS <- list()
  }
  return(msobject)
}

# idPCneg
#' Phosphocholines (PC) annotation for ESI-
#'
#' PC identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for PC in ESI-. Adducts allowed can
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idPCneg} function involves 5 steps. 1) FullMS-based
#' identification of candidate PC as M+CH3COO, M-CH3 or M+CH3COO-CH3. To avoid
#' incorrect annotations of PE as PC, candidates which are present just as M-CH3
#' will be ignored. 2) Search of PC class fragments: 168.0426, 224.0688 or loss
#' of CH3 coeluting with the precursor ion. 3) Search of specific fragments that
#' inform about chain composition in sn1 (lysoPC as M-CH3 resulting from the
#' loss of the FA chain at sn2) and sn2 (lysoPC as M-CH3 resulting from the loss
#' of sn1 or FA as M-H). 4) Look for possible chains structure based on the
#' combination of chain fragments. 5) Check intensity rules to confirm chains
#' position. In this case, lysoPC from sn1 is at least 3 times more intense than
#' lysoPC from sn2.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idPCneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idPCneg <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("M+CH3COO", "M-CH3", "M+CH3COO-CH3"),
                    clfrags = c(168.0426, 224.0688, "pc_M-CH3"),
                    clrequired = c(F, F, F),
                    ftype = c("F", "F", "BB"),
                    chainfrags_sn1 = c("lysopc_M-CH3"),
                    chainfrags_sn2 = c("fa_M-H", "lysopc_M-CH3"),
                    intrules = c("lysopc_sn1/lysopc_sn2"),
                    rates = c("3/1"),
                    intrequired = c(T),
                    coelCutoff = 0.8,
                    dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "PC",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("PC" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious PC annotations removed")
      msobject$detailsAnnotation$PC <- list()
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
  # remove PC which ony appear as M-CH3
  if(length(adducts) > 1 & "M-CH3" %in% adducts){
    candidates <- candidates[candidates$adducts != "M-CH3",]
  }

  if (nrow(candidates) > 0){
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
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
    res <- organizeResults(candidates, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="PC")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$PC <- list()
    msobject$detailsAnnotation$PC$candidates <- candidates
    msobject$detailsAnnotation$PC$classfragments <- classConf$fragments
    msobject$detailsAnnotation$PC$chainfragments <- chainsComb$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$PC$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$PC <- list()
  }
  return(msobject)
}

# idPEneg
#' Phosphoethanolamines (PE) annotation for ESI-
#'
#' PE identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for PE in ESI-. Adducts allowed can
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idPEneg} function involves 5 steps. 1) FullMS-based
#' identification of candidate PE as M-H. 2) Search of PE class fragments:
#' 140.0115, 196.038, 214.048 ion coeluting with the precursor ion. If a loss of
#' CH3 group is found coeluting with any candidate, this will be excluded as it
#' is a characteristic fragment of PC. 3) Search of specific fragments that
#' inform about chain composition in sn1 (lysoPE as M-H resulting from the loss
#' of the FA chain at sn2) and sn2 (lysoPE as M-H resulting from the loss
#' of the FA chain at sn1 or FA chain as M-H). 4) Look for possible
#' chains structure based on the combination of chain fragments. 5) Check
#' intensity rules to confirm chains position. In this case, lysoPE from sn1 is
#' at least 3 times more intense than lysoPE from sn2.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idPEneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idPEneg <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 5,
                    rt,
                    adducts = c("M-H"),
                    clfrags = c(140.0118, 196.038, 214.048, "pe_M-CH3"),
                    clrequired = c(F, F, F, "excluding"),
                    ftype = c("F", "F", "F", "BB"),
                    chainfrags_sn1 = c("lysope_M-H"),
                    chainfrags_sn2 = c("lysope_M-H", "fa_M-H"),
                    intrules = c("lysope_sn1/lysope_sn2"),
                    rates = c("3/1"),
                    intrequired = c(T),
                    coelCutoff = 0.8,
                    dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "PE",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("PE" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious PE annotations removed")
      msobject$detailsAnnotation$PE <- list()
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
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
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
    res <- organizeResults(candidates, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="PE")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$PE <- list()
    msobject$detailsAnnotation$PE$candidates <- candidates
    msobject$detailsAnnotation$PE$classfragments <- classConf$fragments
    msobject$detailsAnnotation$PE$chainfragments <- chainsComb$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$PE$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$PE <- list()
  }
  return(msobject)
}

# idPGneg
#' Phosphoglycerols (PG) annotation for ESI-
#'
#' PG identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for PG in ESI-. Adducts allowed can
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idPGneg} function involves 5 steps. 1) FullMS-based
#' identification of candidate PG as M-H. 2) Search of PG class fragments:
#' 152.9958, 227.0326, 209.022 and neutral loss of 74.0359 coeluting with the
#' precursor ion. 3) Search of specific fragments that inform about chain
#' composition at sn1 (lysoPG as M-H resulting from the loss of the FA chain
#' at sn2) and sn2 (lysoPG as M-H resulting from the loss of the FA chain
#' at sn1 or FA chain as M-H). 4) Look for possible chains structure
#' based on the combination of chain fragments. 5) Check intensity rules to
#' confirm chains position. In this case, lysoPG from sn1 is at least 3 times
#' more intense than lysoPG from sn2.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idPGneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idPGneg <- function(msobject, ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("M-H"),
                    clfrags = c(152.9958, 227.0326, 209.022, 74.0359),
                    clrequired = c(F, F, F, F),
                    ftype = c("F", "F", "F", "NL"),
                    chainfrags_sn1 = c("lysopg_M-H"),
                    chainfrags_sn2 = c("lysopg_M-H", "fa_M-H"),
                    intrules = c("lysopg_sn1/lysopg_sn2"),
                    rates = c("2/1"),
                    intrequired = c(T),
                    coelCutoff = 0.8,
                    dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "PG",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("PG" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious PG annotations removed")
      msobject$detailsAnnotation$PG <- list()
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
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
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
    res <- organizeResults(candidates, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="PG")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$PG <- list()
    msobject$detailsAnnotation$PG$candidates <- candidates
    msobject$detailsAnnotation$PG$classfragments <- classConf$fragments
    msobject$detailsAnnotation$PG$chainfragments <- chainsComb$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$PG$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$PG <- list()
  }
  return(msobject)
}

# idPIneg
#' Phosphoinositols (PI) annotation for ESI-
#'
#' PI identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for PI in ESI-. Adducts allowed can
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idPIneg} function involves 5 steps. 1) FullMS-based
#' identification of candidate PI as M-H. 2) Search of PI class fragments:
#' 241.0115, 223.0008, 259.0219 and 297.0375 coeluting with the precursor
#' ion. 3) Search of specific fragments that inform about chain composition at
#' sn1 (lysoPI as M-H resulting from the loss of the FA chain at sn2 or lysoPA
#' as M-H if it also losses the head group) and sn2 (lysoPI or lysoPA as M-H
#' resulting from the loss of the FA chain at sn1 or FA chain as M-H). 4) Look
#' for possible chains structure based on the combination of chain fragments.
#' 5) Check intensity rules to confirm chains position. In this case, lysoPI or
#' lysoPA from sn1 is at least 3 times more intense than lysoPI or lysoPA from
#'  sn2.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idPIneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idPIneg <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("M-H"),
                    clfrags = c(241.0115, 223.0008, 259.0219, 297.0375),
                    clrequired = c(F, F, F, F),
                    ftype = c("F", "F", "F", "F"),
                    chainfrags_sn1 = c("lysopi_M-H", "lysopa_M-H"),
                    chainfrags_sn2 = c("lysopi_M-H", "lysopa_M-H", "fa_M-H"),
                    intrules = c("lysopi_sn1/lysopi_sn2", "lysopa_sn1/lysopa_sn2"),
                    rates = c("3/1", "3/1"),
                    intrequired = c(F, F),
                    coelCutoff = 0.8,
                    dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "PI",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("PI" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious PI annotations removed")
      msobject$detailsAnnotation$Cer <- list()
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
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
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
    res <- organizeResults(candidates, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="PI")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$PI <- list()
    msobject$detailsAnnotation$PI$candidates <- candidates
    msobject$detailsAnnotation$PI$classfragments <- classConf$fragments
    msobject$detailsAnnotation$PI$chainfragments <- chainsComb$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$PI$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$PI <- list()
  }
  return(msobject)
}

# idPSneg
#' Phosphoserines (PS) annotation for ESI-
#'
#' PS identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for PS in ESI-. Adducts allowed can
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idPSneg} function involves 5 steps. 1) FullMS-based
#' identification of candidate PS as M-H or M+Na-2H. 2) Search of PS class
#' fragments: neutral loss of 87.032 (serine) coeluting with the precursor ion.
#' 3) Search of specific fragments that inform about chain composition at sn1
#' (lysoPA as M-H or M-H-H2O resulting from the loss of the FA chain at sn2 and
#' the head group) and sn2 (lysoPA as M-H or M-H-H2O resulting from the loss of
#' the FA chain at sn1 and the head group or FA chain as M-H). 4) Look for
#' possible chains structure based on the combination of chain fragments.
#' 5) Check intensity rules to confirm chains position. In this case, lysoPA
#' from sn1 is at least 3 times more intense than lysoPA from sn2.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idPSneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idPSneg <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("M-H", "M+Na-2H"),
                    clfrags = c(87.032, 152.9958),
                    clrequired = c(F, F),
                    ftype = c("NL", "F"),
                    chainfrags_sn1 = c("lysopa_M-H", "lysopa_M-H-H2O"),
                    chainfrags_sn2 = c("lysopa_M-H", "lysopa_M-H-H2O", "fa_M-H"),
                    intrules = c("lysopa_sn1/lysopa_sn2"),
                    rates = c("3/1"),
                    intrequired = c(T),
                    coelCutoff = 0.8,
                    dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "PS",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("PS" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious PS annotations removed")
      msobject$detailsAnnotation$PS <- list()
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
  candidates <- findCandidates(MS1, dbs$psdb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)

  if (nrow(candidates) > 0){
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
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
    res <- organizeResults(candidates, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="PS")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$PS <- list()
    msobject$detailsAnnotation$PS$candidates <- candidates
    msobject$detailsAnnotation$PS$classfragments <- classConf$fragments
    msobject$detailsAnnotation$PS$chainfragments <- chainsComb$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$PS$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$PS <- list()
  }
  return(msobject)
}

# idSphneg
#' Sphingoid bases (Sph) annotation for ESI-
#'
#' Sph identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for Sph in ESI-. Adducts allowed can
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idSphneg} function involves 2 steps. 1) FullMS-based
#' identification of candidate Sph as M-H. 2) Search of Sph class fragments:
#' neutral loss of 1 or 2 H2O molecules.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (in this
#' case, as Sph only have one chain, only Subclass and FA level are possible)
#' and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idSphneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idSphneg <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3, rt,
                     adducts = c("M-H"),
                     clfrags = c("sph_M-H-H2O", "sph_M-H-2H2O"),
                     clrequired = c(F, F),
                     ftype = c("BB", "BB"),
                     coelCutoff = 0.8,
                     dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "Sph",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("Sph" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious Sph annotations removed")
      msobject$detailsAnnotation$Sph <- list()
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
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }

    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)

    # prepare output
    res <- organizeResults(candidates, clfrags, classConf, chainsComb = list(),
                           intrules  = c(), intConf = list(), nchains = 0,
                           class="Sph")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$Sph <- list()
    msobject$detailsAnnotation$Sph$candidates <- candidates
    msobject$detailsAnnotation$Sph$classfragments <- classConf$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$Sph$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$Sph <- list()
  }
  return(msobject)
}

# idSphPpneg
#' Sphingoid bases phosphate (SphP) annotation for ESI-
#'
#' SphP identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for SphP in ESI-. Adducts allowed can
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idSphpos} function involves 2 steps. 1) FullMS-based
#' identification of candidate SphP as M-H. 2) Search of SphP class fragments:
#' 78.9585, 96.969 or neutral loss of 1 H2O molecule.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (in this
#' case, as SphP only have one chain, only Subclass and FA level are possible)
#' and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idSphPneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idSphPneg <- function(msobject,
                      ppm_precursor = 5,
                      ppm_products = 10,
                      rttol = 3,
                      rt,
                      adducts = c("M-H"),
                      clfrags = c(78.9585, 96.9691, "sphP_M-H-H2O"),
                      clrequired = c(F, F, F),
                      ftype = c("F", "F", "BB"),
                      coelCutoff = 0.8,
                      dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "SphP",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("SphP" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious SphP annotations removed")
      msobject$detailsAnnotation$SphP <- list()
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
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }

    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)

    # prepare output
    res <- organizeResults(candidates, clfrags, classConf, chainsComb = list(),
                           intrules  = c(), intConf = list(), nchains = 0,
                           class="SphP")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$SphP <- list()
    msobject$detailsAnnotation$SphP$candidates <- candidates
    msobject$detailsAnnotation$SphP$classfragments <- classConf$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$SphP$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$SphP <- list()
  }
  return(msobject)
}

# idCerneg
#' Ceramides (Cer) annotation for ESI-
#'
#' Cer identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
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
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idCerneg} function involves 5 steps. 1) FullMS-based
#' identification of candidate Cer as M-H and M+CH3COO. 2) Search of Cer class
#' fragments: there are no class fragment by default. 3) Search of specific
#' fragments that inform about the sphingoid base (Sph as M-H-2H2O resulting
#' from the loss of the FA chain or loss of part of the sphingoid base) and the
#' FA chain (FA as M-H but with a N intead of an O, what means a mass difference
#' of 1.9918 from the exact mass of the FA). 4) Look for possible chains
#' structure based on the combination of chain fragments. 5) Check intensity
#' rules to confirm chains position. In this case, there are no intensity rules
#' by default.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idCerneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idCerneg <- function(msobject,
                     ppm_precursor = 5,
                     ppm_products = 10,
                     rttol = 3,
                     rt,
                     adducts = c("M-H","M+CH3COO"),
                     clfrags = c(),
                     clrequired = c(),
                     ftype = c(),
                     chainfrags_sn1 = c("NL-nlsph_M-H", "sph_M-H-2H2O", "sph_M-H-H2O"),
                     chainfrags_sn2 = c("fa_Mn-1.9918"),
                     intrules = c(),
                     rates = c(),
                     intrequired = c(),
                     coelCutoff = 0.8,
                     dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "Cer",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("Cer" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious Ceramide annotations removed")
      msobject$detailsAnnotation$Cer <- list()
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
  candidates <- findCandidates(MS1, dbs$cerdb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)

  if (nrow(candidates) > 0){
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
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
    res <- organizeResults(candidates, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="Cer")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$Cer <- list()
    msobject$detailsAnnotation$Cer$candidates <- candidates
    msobject$detailsAnnotation$Cer$classfragments <- classConf$fragments
    msobject$detailsAnnotation$Cer$chainfragments <- chainsComb$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$Cer$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$Cer <- list()
  }
  return(msobject)
}

# idCLneg
#' Cardiolipines (CL) annotation for ESI-
#'
#' CL identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for CL in ESI-. Adducts allowed can
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
#' the chain fragments in sn2 position. See \link{chainFrags} for details.
#' @param chainfrags_sn3 character vector containing the fragmentation rules for
#' the chain fragments in sn3 position. See \link{chainFrags} for details.
#' @param chainfrags_sn4 character vector containing the fragmentation rules for
#' the chain fragments in sn4 position. See \link{chainFrags} for details.
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
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idCLneg} function involves 5 steps. 1) FullMS-based
#' identification of candidate CL as M-H or M-2H. 2) Search of CL class fragments:
#' no class fragments are searched by defaults as they use to have bad coelution
#' scores. 3) Search of specific fragments that inform about chain composition
#' at sn1 (lysoPA as M-H-H2O), sn2 (lysoPA as M-H-H2O), sn3 (lysoPA as M-H-H2O)
#' and sn4 (lysoPA as M-H-H2O). 4) Look for possible chains structure based on
#' the combination of chain fragments. 5) Check intensity rules to confirm
#' chains position. For CL there are no intensity rules by default.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idCLneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idCLneg <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 5,
                    rt,
                    adducts = c("M-H", "M+Na-2H"),
                    clfrags = c(),
                    clrequired = c(),
                    ftype = c(),
                    chainfrags_sn1 = c("lysopa_M-H-H2O"),
                    chainfrags_sn2 = c("lysopa_M-H-H2O"),
                    chainfrags_sn3 = c("lysopa_M-H-H2O"),
                    chainfrags_sn4 = c("lysopa_M-H-H2O"),
                    intrules = c("Unknown"),
                    rates = c(),
                    intrequired = c(),
                    coelCutoff = 0.8,
                    dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
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
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "CL",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("CL" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious CL annotations removed")
      msobject$detailsAnnotation$CL <- list()
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
  candidates <- findCandidates(MS1, dbs$cldb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)

  if (nrow(candidates) > 0){
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
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
    sn4 <- chainFrags(coelfrags, chainfrags_sn4, ppm_products, candidates, sn1,
                      dbs)

    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=4, sn1, sn2, sn3, sn4)

    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=4,
                                   chainsComb)

    # prepare output
    res <- organizeResults(candidates, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 4, class="CL")


    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$CL <- list()
    msobject$detailsAnnotation$CL$candidates <- candidates
    msobject$detailsAnnotation$CL$classfragments <- classConf$fragments
    msobject$detailsAnnotation$CL$chainfragments <- chainsComb$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$CL$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$CL <- list()
  }
  return(msobject)
}

# idBAneg
#' Bile Acids (BA) annotation for ESI-
#'
#' BA identification based on fragmentation patterns for LC-MS/MS DIA or DDA
#' data acquired in negative mode.
#'
#' @param msobject an msobject returned by \link{dataProcessing}.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for BA in ESI-. Adducts allowed can
#' be modified in the adducsTable (dbs argument).
#' @param conjfrag character vector containing the fragmentation rules for
#' the BA-conjugates. By default just taurine and glycine are considered,
#' but baconjdb can be modified to add more possible conjugates.
#' See \link{chainFrags} for details. It can also be an empty vector.
#' @param bafrag character vector containing the fragmentation rules for
#' other BA fragments. See \link{chainFrags} for details. It can be an empty
#' vector.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#'
#' @return annotated msobject (list with several elements). The results element
#' is a data frame that shows: ID, class of lipid, CDB (total number of carbons
#' and double bounds), FA composition (specific chains composition if it has
#' been confirmed), m.z, RT (in seconds), I (intensity), Adducts, ppm (m.z error),
#' confidenceLevel (Subclass, FA level, where chains are known but not their
#' positions, or FA position level), peakID, and PFCS (parent-fragment coelution
#' score mean of all fragments used for the identification); and the
#' annotatedPeaklist element shows the original MS1 peaklist with the annotations
#' on it.
#'
#' @details \code{idBAneg} function involves 3 steps. 1) FullMS-based
#' identification of candidate BA as M-H. 2) Search of BA-conjugate fragments if
#' required. 3) Search of fragments coming from the loss of H2O.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (MS-only
#' if no rules are defined, or Subclass level if they are supported by fragments)
#' and PFCS (parent-fragment coelution score mean of all fragments used for the
#' identification).
#'
#' @note This function has been writen based on fragmentation patterns
#' observed for three different platforms (QTOF 6550 from Agilent, Sinapt G2-Si
#' from Waters and Q-exactive from Thermo), but it may need to be customized for
#' other platforms or acquisition settings.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idBAneg(LipidMSdata2::msobjectDIAneg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idBAneg <- function(msobject,
                    ppm_precursor = 5,
                    ppm_products = 10,
                    rttol = 3,
                    rt,
                    adducts = c("M-H"),
                    conjfrag = c("baconj_M-H"),
                    bafrag = c("ba_M-H-H2O", "ba_M-H-2H2O"),
                    coelCutoff = 0.8,
                    dbs){
  ##############################################################################
  # check arguments
  if (msobject$metaData$generalMetadata$polarity != "negative"){
    stop("Data wasn't acquired in negative mode")
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (!all(c("metaData", "processing", "MS1", "MS2", "peaklist") %in% names(msobject))){
    stop("Wrong msobject format")
  }
  if (!msobject$metaData$acquisitionmode %in% c("DIA", "DDA")){
    stop("Acquisition mode must be DIA or DDA")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
  }
  ##############################################################################
  # extract data from msobject
  # Peaklist MS1: remove isotopes
  MS1 <- msobject$peaklist$MS1
  MS1 <- MS1[MS1$isotope %in% c("[M+0]"),
             !colnames(MS1) %in% c("isotope", "group")]
  # Peaklist MS2: remove isotopes
  if (msobject$metaData$acquisitionmode == "DDA"){
    MS2 <- msobject$MS2[,c("m.z", "RT", "int", "peakID")]
  } else {
    MS2 <- msobject$peaklist$MS2[,c("m.z", "RT", "int", "peakID")]
  }
  rawData <- rbind(msobject$MS1, msobject$MS2)
  # if acquisition mode is DDA, extract precursors
  if (msobject$metaData$acquisitionmode == "DDA"){
    precursors <- msobject$metaData$scansMetadata[msobject$metaData$scansMetadata$collisionEnergy > 0 &
                                                    msobject$metaData$scansMetadata$msLevel == 2,
                                                  c("retentionTime", "precursor", "Scan")]
  }
  ##############################################################################
  # Remove previous ceramide annotations
  if ("results" %in% names(msobject)){
    if (nrow(msobject$results) > 0){
      msobject$results <- msobject$results[msobject$results$Class != "BA",]
    }
  }
  if ("detailsAnnotation" %in% names(msobject)){
    if("BA" %in% names(msobject$detailsAnnotation)){
      cat("\nPrevious BA annotations removed")
      msobject$detailsAnnotation$BA <- list()
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
  candidates <- findCandidates(MS1, dbs$badb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)


  if (nrow(candidates) > 0){
    if (msobject$metaData$acquisitionmode == "DIA"){
      if (nrow(rawData) == 0){
        coelCutoff <- 0 # if no rawData is supplied, coelution score between precursors and fragments will be ignored
      }
      # isolation of coeluting fragments
      coelfrags <- coelutingFrags(candidates, MS2, rttol, rawData,
                                  coelCutoff = coelCutoff)
    } else if (msobject$metaData$acquisitionmode == "DDA"){
      coelCutoff <- 0
      coelfrags <- ddaFrags(candidates, precursors, rawData, ppm = ppm_products)
    }

    # search fragments
    check <- rep(list(vector()), nrow(candidates))
    conjugates <- rep(list(vector()), nrow(candidates))
    bas <- rep(list(vector()), nrow(candidates))
    if (length(conjfrag) > 0){
      conjugates <- chainFrags(coelfrags, conjfrag, ppm_products, dbs = dbs,
                               candidates = candidates)
      for (c in 1:nrow(candidates)){
        if (nrow(conjugates[[c]]) > 0){
          if (dbs$badb[dbs$badb$total ==
                       candidates$cb[c], "conjugate"] %in%
              conjugates[[c]]$cb | dbs$badb[dbs$badb$total ==
                                            candidates$cb[c],
                                            "conjugate"] == ""){
            check[[c]] <- append(check[[c]], TRUE)
          } else {
            check[[c]] <- append(check[[c]], FALSE)
          }
        } else if (dbs$badb[dbs$badb$total ==
                            candidates$cb[c],
                            "conjugate"] == ""){
          check[[c]] <- append(check[[c]], TRUE)
          conjugates[[c]] <- data.frame(0,0,0,0,0,0,0,0)
          colnames(conjugates[[c]]) <- c("cb", "m.z", "RT", "int",
                                         "peakID", "coelScore", "db", "adduct")
        } else {
          check[[c]] <- append(check[[c]], FALSE)
          conjugates[[c]] <- data.frame(0,0,0,0,0,0,0,0)
          colnames(conjugates[[c]]) <- c("cb", "m.z", "RT", "int",
                                         "peakID", "coelScore", "db", "adduct")
        }
      }
      classConf <- list()
    }
    if (length(bafrag) > 0){
      bas <- chainFrags(coelfrags, bafrag, ppm_products, dbs = dbs,
                        candidates = candidates)
      for (c in 1:nrow(candidates)){
        if(nrow(bas[[c]]) > 0){
          if (candidates$cb[c] %in% bas[[c]]$cb ||
              dbs$badb[dbs$badb$total == candidates$cb[c], "base"]
              %in% bas[[c]]$cb){
            check[[c]] <- append(check[[c]], TRUE)
          } else {
            check[[c]] <- append(check[[c]], FALSE)
          }
        } else {
          check[[c]] <- append(check[[c]], FALSE)
          bas[[c]] <- data.frame(0,0,0,0,0,NA,0,0)
          colnames(bas[[c]]) <- c("cb", "m.z", "RT", "int",
                                  "peakID", "coelScore", "db", "adduct")
        }
      }
    }
    check <- matrix(unlist(check), nrow = nrow(candidates), byrow = T)
    classConf$presence <- check
    classConf$passed <- unlist(apply(check, 1, sum)) > 0
    classConf$fragments <- Map(rbind, conjugates, bas)

    # prepare output
    if (length(bafrag) > 0 | length(conjfrag) > 0){
      clfrags <- c("fragment")
    } else {
      clfrags <- c()
    }

    res <- organizeResults(candidates, clfrags = clfrags, classConf = classConf,
                           chainsComb = c(), intrules = c(),
                           intConf = c(), nchains = 0, class="BA")

    # update msobject
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$BA <- list()
    msobject$detailsAnnotation$BA$candidates <- candidates
    msobject$detailsAnnotation$BA$classfragments <- classConf$fragments
    msobject$detailsAnnotation$BA$chainfragments <- bas$fragments
    if (msobject$metaData$acquisitionmode == "DDA"){
      msobject$detailsAnnotation$BA$coelfrags <- coelfrags
    }
  } else {
    res <- data.frame()
    if ("results" %in% names(msobject)){
      msobject$results <- rbind(msobject$results, res)
    } else {
      msobject$results <- res
    }
    msobject$detailsAnnotation$BA <- list()
  }
  return(msobject)
}

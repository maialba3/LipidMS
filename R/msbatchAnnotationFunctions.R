# annotatemsbatch
#' Lipid annotation for an msbatch
#'
#' Summarize annotation results of an msbatch into the feature table
#'
#' @param msbatch msbatch
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 5 seconds.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param lipidClassesPos classes of interest in ESI+.
#' @param lipidClassesNeg classes of interest in ESI-.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#' @param simplifyAnnotations logical. If TRUE, only the most frequent id will be 
#' kept (recommended when only pool samples have been acquired in DIA or DDA). If 
#' FALSE, all annotations will be shown.
#' @param parallel logical.
#' @param ncores number of cores to be used in case parallel is TRUE.
#' @param global_gb numeric. Gigabytes to set as future.globals.maxSize **inside** the function.
#' Defaults to `getOption("LipidMS.future.globals.maxSizeGB", Inf)`.
#'
#' @return msbatch
#' 
#' @examples
#' \dontrun{
#' msbatch <- annotatemsbatch(msbatch)
#' 
#' msbatch$features
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
annotatemsbatch <- function(msbatch,
                            ppm_precursor = 5,
                            ppm_products = 10,
                            rttol = 5,
                            coelCutoff = 0.8,
                            lipidClassesPos = c("MG", "LPC", "LPE", "PC", "PCo",
                                                "PCp", "PE", "PEo", "PEp", "PG",
                                                "PI", "Sph", "SphP", "Cer", 
                                                "CerP", "AcylCer", "SM", 
                                                "Carnitines", "CE", "DG", "TG"),
                            lipidClassesNeg = c("FA", "FAHFA", "LPC", "LPE", 
                                                "LPG", "LPI", "LPS", "PC", "PCo", 
                                                "PCp", "PE", "PEo", "PEp", "PG", 
                                                "PI", "PS", "Sph", "SphP", 
                                                "Cer", "CerP", "AcylCer", "SM", 
                                                "CL", "BA"),
                            dbs,
                            simplifyAnnotations = FALSE,
                            parallel = FALSE,
                            ncores,
                            global_gb = getOption("LipidMS.future.globals.maxSizeGB", Inf)){
  
  
  ##############################################################################
  # Check arguments
  if (parallel){
    if (missing(ncores)){
      stop("ncores argument is required if parallel is TRUE")
    }
    if (ncores > parallel::detectCores()){
      ncores <- parallel::detectCores() - 1
      message("ncores is greater than availables cores. ", ncores, " will be used.")
    }
  }
  if (missing(dbs)){
    dbs <- assignDB()
  }
  
  ##############################################################################
  # lipid annotation in samples acquired in DIA or DDA mode
  toannotate <- which(msbatch$metaData$acquisitionmode %in% c("DIA", "DDA"))
  
  if (length(toannotate) > 0){
    if (!parallel) ncores <- 1L
    
    res <- .with_future_plan(
      workers = ncores,
      globals_gb = global_gb,
      expr = {
        future.apply::future_lapply(
          toannotate,
          function(ix){
            mso <- msbatch$msobjects[[ix]]
            if (mso$metaData$generalMetadata$polarity == "positive"){
              idPOS(msobject = mso, 
                    ppm_precursor = ppm_precursor,
                    ppm_products = ppm_products,
                    rttol = rttol,
                    coelCutoff = coelCutoff,
                    lipidClasses = lipidClassesPos,
                    dbs = dbs)
            } else {
              idNEG(msobject = mso,
                    ppm_precursor = ppm_precursor,
                    ppm_products = ppm_products,
                    rttol = rttol,
                    coelCutoff = coelCutoff,
                    lipidClasses = lipidClassesNeg,
                    dbs = dbs)
            }
          },
          future.seed = TRUE
        )
      }
    )
    
    for (k in seq_along(toannotate)) {
      msbatch$msobjects[[toannotate[k]]] <- res[[k]]
    }
    
    
    # remove previous annotations and write results on the features table
    msbatch$features <- msbatch$features[, !colnames(msbatch$features) %in% 
                                           c("LipidMSid", "Adduct", 
                                             "confidenceLevel", "Score", 
                                             "ScoreInt", "nsamples")]
    msbatch <- joinAnnotationResults(msbatch, 
                                     simplifyAnnotations = simplifyAnnotations)
  } else {
    message("No available files to search for lipid annotations (DDA or DIA acquired samples are required)")
  }
  
  return(msbatch)
}


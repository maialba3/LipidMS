#' Annotate lipids in LC-MS/MS data
#'
#' @description
#' Performs lipid annotation on an [xcms::XcmsExperiment] containing
#' chromatographic peaks and MS/MS data. LipidMS combines accurate-mass
#' matching at the MS1 level with adduct grouping and rule-based evaluation
#' of MS2 spectra to generate lipid annotations.
#'
#' The function supports data acquired using data-dependent acquisition (DDA)
#' and data-independent acquisition (DIA/SWATH). Full-scan samples can also be
#' present in the input object and contribute to MS1 candidate detection and
#' isotope annotation, while MS2 annotation is performed using DDA and
#' DIA/SWATH samples.
#'
#' @details
#' `annotateLipids()` implements the main LipidMS annotation workflow.
#'
#' The workflow consists of the following steps:
#'
#' 1. Validation of the input [xcms::XcmsExperiment], including sample
#'    metadata, chromatographic peaks, acquisition modes and MS2-capable
#'    samples.
#' 2. Annotation of MS1 isotopes using [annotateIsotopes()] if isotope
#'    information is not already present in `chromPeakData(x)`.
#' 3. Generation of MS1 lipid candidates by matching chromatographic peak
#'    masses against the LipidMS lipid databases and supported adducts.
#' 4. Grouping and scoring of co-eluting adduct candidates corresponding to
#'    the same neutral lipid species.
#' 5. Extraction of experimental MS2 spectra. Acquired spectra are associated
#'    with precursor chromatographic peaks for DDA samples, whereas MS2 spectra
#'    are reconstructed from DIA/SWATH data.
#' 6. Evaluation of the MS2 spectra against the LipidMS fragmentation rules.
#' 7. Storage of peak-level and, when available, feature-level lipid
#'    annotations in the returned [xcms::XcmsExperiment].
#' 8. Standarize lipid names and identifiers with RefMet package.
#'
#' The input object must contain MS1 chromatographic peaks and a non-empty
#' `sampleData()` table. Acquisition mode must be specified in the
#' `acquisition_mode` column of `sampleData(x)`. Accepted acquisition modes
#' include full-scan MS, DDA, DIA and SWATH. At least one DDA or DIA/SWATH
#' sample must be present.
#'
#' @param x An [xcms::XcmsExperiment] containing detected MS1 chromatographic
#'   peaks and MS data. The object must contain non-empty sample metadata with
#'   an `acquisition_mode` column. At least one sample must contain DDA or
#'   DIA/SWATH MS2 data.
#' @param paramannot A named `list` containing the parameters used for MS1 and
#'   MS2 lipid annotation. The default parameter structure can be generated
#'   with [defaultAnnotationParameters()]. Important parameters include
#'   precursor and product ion mass tolerances, retention-time tolerance,
#'   isotope annotation settings and minimum MS2 intensity.
#' @param paramdia A named `list` containing parameters used for reconstruction
#'   of DIA/SWATH spectra. These parameters are required only when DIA or SWATH
#'   samples are present. They can be generated with
#'   [defaultDIAExtractionParameters()]. For datasets containing only DDA
#'   samples, an empty `list()` can be supplied.
#' @param dbs Lipid database object used for MS1 candidate generation.
#'   By default, databases are generated with [createLipidDbs()].
#' @param rules Fragmentation-rule object used to evaluate MS2 spectra.
#'   By default, the rules returned by [defaultFragmentRules()] are used.
#' @param adduct_table A `data.frame` defining the supported ion adducts and
#'   their mass transformations. By default, [defaultAdductTable()] is used.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object controlling
#'   parallel processing. Defaults to [BiocParallel::SerialParam()].
#' @param verbose `logical(1)` indicating whether progress information and
#'   timing for the different annotation steps should be printed.
#'   Defaults to `TRUE`.
#'
#' @return
#' An [xcms::XcmsExperiment] containing the lipid annotation results.
#' Peak-level annotation information is stored in `chromPeakData(x)`.
#' When grouped chromatographic features are available, feature-level
#' annotation summaries are also added to `featureDefinitions(x)`.
#' Additional LipidMS results, parameters and annotation metadata are stored
#' in the object's metadata.
#'
#' @seealso
#' [defaultAnnotationParameters()],
#' [defaultDIAExtractionParameters()],
#' [annotateIsotopes()],
#' [findMs1LipidCandidates()],
#' [groupCoelutingAdducts()],
#' [extractMs2Spectra()],
#' [matchFragmentRules()]
#'
#' @examples
#' data(xcmsexperiment_negative)
#'
#' param <- defaultAnnotationParameters("negative")
#'
#' result <- annotateLipids(
#'     xcmsexperiment_negative,
#'     paramannot = param,
#'     paramdia = list(),
#'     verbose = FALSE
#' )
#'
#' result
#'
#' @export
annotateLipids <- function(x, paramannot, paramdia, dbs = createLipidDbs(),
                           rules = defaultFragmentRules(),
                           adduct_table = defaultAdductTable(),
                           BPPARAM = BiocParallel::SerialParam(),
                           verbose = TRUE){


  validation <- .validate_annotate_lipids_input(x)

  .validate_paramannot(paramannot)

  acquisition <- validation$acquisition
  sample_names <- validation$sample_names
  idx_dda <- validation$idx_dda
  idx_dia <- validation$idx_dia
  idx_ms2 <- validation$idx_ms2
  has_grouped_features <- validation$has_grouped_features

  if (length(idx_dia)) {
    .validate_paramdia(paramdia)
  }

  if (isTRUE(verbose)) {
    cp_full <- xcms::chromPeaks(x)
    sample_index_full <- .get_chrom_peak_sample_idx(cp_full)
    ms2_peak_ids_full <- rownames(cp_full)[sample_index_full %in% idx_ms2]

    if (isTRUE(verbose)){
      cat("\nAnnotation subset strategy:\n")
      cat("  Full samples: ", length(sample_names), "\n", sep = "")
      cat("  DDA samples: ", length(idx_dda), "\n", sep = "")
      cat("  DIA/SWATH samples: ", length(idx_dia), "\n", sep = "")
      cat("  MS2-capable samples used for annotation subset: ", length(idx_ms2), "\n", sep = "")
      cat("  Full chromPeaks: ", nrow(cp_full), "\n", sep = "")
      cat("  Full chromPeaks in DDA/DIA/SWATH samples: ", length(ms2_peak_ids_full), "\n", sep = "")
      cat("  MS2 sample names: ", paste(sample_names[idx_ms2], collapse = ", "), "\n", sep = "")
    }

  }

  if (!"is_isotope" %in% colnames(xcms::chromPeakData(x))) {
    x <- .step_time("annotateIsotopes", {
      if (isTRUE(verbose)) {
        cat("\nDiscarding isotopes using the full dataset before MS2 subsetting...\n")
      }
      annotateIsotopes(
        x,
        mslevel = 1L,
        param = list(
          dmzprecursor = paramannot$dmzprecursor,
          mzabs = paramannot$isotope_mzabs,
          maxcharge = paramannot$isotope_maxcharge,
          maxiso = paramannot$isotope_maxiso,
          minfrac = paramannot$isotope_minfrac,
          intval = paramannot$intval,
          rt_tolerance = paramannot$isotope_rttol
        ),
        BPPARAM = .get_param_value(
          paramannot,
          "BPPARAM_isotopes",
          BiocParallel::SerialParam()
        )
      )
    }, verbose = verbose)
  } else {
    if (isTRUE(verbose)) {
      cat("\nIsotopes already identified in full object.\n")
    }
  }

  ms1_candidates <- .step_time("findMs1LipidCandidates", {
    if (isTRUE(verbose)) {
      cat("\nSearching for lipid candidates on the full dataset...\n")
    }
    findMs1LipidCandidates(
      x = x,
      dbs = dbs,
      paramannot = paramannot,
      adduct_table = adduct_table
    )
  }, verbose = verbose)

  ms1_candidates <- .step_time("groupCoelutingAdducts", {
    groupCoelutingAdducts(
      candidates = ms1_candidates,
      paramannot = paramannot
    )
  }, verbose = verbose)

  if (!"chrom_peak_id" %in% colnames(ms1_candidates)) {
    stop("'ms1_candidates' must contain column 'chrom_peak_id'.", call. = FALSE)
  }

  x_annot <- .step_time("subset_ms2_data", {
    .subset_xcms_experiment_samples(x, idx_ms2)
  }, verbose = verbose)

  # Keep only the globally detected MS1 candidates that belong to MS2-capable
  # samples for the MS2 extraction/matching steps.
  ms2_subset_peak_ids <- rownames(xcms::chromPeaks(x_annot))
  ms1_candidates_ms2 <- ms1_candidates[
    ms1_candidates$chrom_peak_id %in% ms2_subset_peak_ids,
    ,
    drop = FALSE
  ]

  requested_peak_ids <- unique(ms1_candidates_ms2$chrom_peak_id)

  if (isTRUE(verbose)) {
    cat("\nMS1 candidates in full dataset: ", nrow(ms1_candidates), "\n", sep = "")
    cat("MS1 candidates retained for MS2 extraction: ",
        nrow(ms1_candidates_ms2), "\n", sep = "")
    cat("Requested MS2 peak IDs: ", length(requested_peak_ids), "\n", sep = "")
  }

  if (!length(requested_peak_ids)) {
    warning("No full-dataset MS1 candidates were present in DDA/DIA/SWATH samples; returning unannotated object.")
    ms2_data <- list(samples = list(), x_dia = NULL)
    ms2_annotation <- list(
      annotations = data.frame(),
      evidence = data.frame(),
      all_annotations = data.frame()
    )
  } else {
    ms2_data <- .step_time("extractMs2Spectra", {
      if (isTRUE(verbose)) {
        cat("\nExtracting MS2 data for lipid candidates from DDA/DIA/SWATH samples...\n")
      }
      extractMs2Spectra(
        x = x_annot,
        paramannot = paramannot,
        paramdia = paramdia,
        requested_peak_ids = requested_peak_ids,
        BPPARAM = BPPARAM,
        verbose = verbose
      )
    }, verbose = verbose)

    ms2_annotation <- .step_time("matchFragmentRules", {
      if (isTRUE(verbose)) {
        cat("\nSearching expected fragments in MS2 spectra...\n")
      }
      matchFragmentRules(
        candidates = ms1_candidates_ms2,
        ms2_data = ms2_data,
        rules = rules,
        dbs = dbs,
        paramannot = paramannot
      )
    }, verbose = verbose)
  }

  x_annot <- .step_time("addMs2AnnotationsToSubsetChromPeakData", {
    if (isTRUE(verbose)) {
      cat("\nAdding lipid annotations to subset chromPeakData...\n")
    }

    addMs2AnnotationsToChromPeakData(
      x = x_annot,
      ms2_annotation = ms2_annotation,
      candidates = ms1_candidates_ms2
    )
  }, verbose = verbose)

  x <- .step_time("merge_subset_annotations_into_full_dataset", {
    .merge_chrom_peak_data_columns_by_id(
      x_full = x,
      x_subset = x_annot,
      verbose = verbose
    )
  }, verbose = verbose)

  if (has_grouped_features) {
    x <- .step_time("addFeatureLipidAnnotations", {
      if (isTRUE(verbose)) {
        cat(
          "\nAdding lipid annotations to featureDefinitions...\n"
        )
      }

      addFeatureLipidAnnotations(x)
    }, verbose = verbose)
  } else {
    warning(
      paste(
        "Lipid annotations were stored in chromPeakData(x),",
        "but featureDefinitions(x) was not updated because",
        "valid grouped features were not found.",
        "Run chromatographic peak grouping and then call",
        "addFeatureLipidAnnotations(x) if feature-level summaries",
        "are required."
      ),
      call. = FALSE
    )
  }

  x <- .step_time("addRefMetAnnotations", {

    if (isTRUE(verbose)) {
      cat(
        "\nAdding standarized ids from RefMet to annotations...\n"
      )
    }

    .add_refmet_annotations(x)

  }, verbose = verbose)

  ms2_data <- annotateAcceptedFragmentsInMs2Data(
    ms2_data = ms2_data,
    ms2_annotation = ms2_annotation,
    x = x,
    ppm = paramannot$dmzproducts
  )

  x <- addLipidMSMetadata(
    x = x,
    ms1_candidates = ms1_candidates,
    ms2_data = ms2_data,
    ms2_annotation = ms2_annotation,
    params = paramannot,
    rules = rules
  )

  if (isTRUE(verbose)) {
    cat("\nLipidMS annotation workflow finished.\n")
  }

  return(x)
}



#' Validate an XcmsExperiment for lipid annotation
#'
#' Validate the structure of an `XcmsExperiment` before lipid annotation.
#'
#' The function checks:
#'
#' - that the object is an `XcmsExperiment`
#' - that sample metadata are available
#' - that MS1 chromatographic peaks have been detected
#' - that `chromPeaks(x)` and `chromPeakData(x)` are aligned
#' - that chromatographic peak sample indices are valid
#' - that acquisition modes are available and valid
#' - that at least one DDA, DIA or SWATH sample is present
#' - whether valid grouped features are available
#'
#' @param x `XcmsExperiment`.
#'
#' @return Named list containing:
#'
#' - `acquisition`
#' - `sample_names`
#' - `idx_dda`
#' - `idx_dia`
#' - `idx_ms2`
#' - `has_grouped_features`
#'
#' @noRd
#'
#' @keywords internal
.validate_annotate_lipids_input <- function(x) {
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("Package 'xcms' is required.", call. = FALSE)
  }

  if (!requireNamespace("MsExperiment", quietly = TRUE)) {
    stop("Package 'MsExperiment' is required.", call. = FALSE)
  }

  if (!inherits(x, "XcmsExperiment")) {
    stop(
      "'x' must be an XcmsExperiment.",
      call. = FALSE
    )
  }

  sample_data <- tryCatch(
    as.data.frame(MsExperiment::sampleData(x)),
    error = function(e) NULL
  )

  if (is.null(sample_data) || !nrow(sample_data)) {
    stop(
      "'x' must contain non-empty sampleData.",
      call. = FALSE
    )
  }

  n_samples <- nrow(sample_data)

  chrom_peaks <- tryCatch(
    xcms::chromPeaks(x),
    error = function(e) NULL
  )

  if (is.null(chrom_peaks) || !nrow(chrom_peaks)) {
    stop(
      paste(
        "'x' does not contain chromatographic peaks.",
        "Run MS1 chromatographic peak detection before annotateLipids()."
      ),
      call. = FALSE
    )
  }

  if (is.null(rownames(chrom_peaks))) {
    stop(
      "'chromPeaks(x)' must have row names.",
      call. = FALSE
    )
  }

  required_peak_columns <- c(
    "mz",
    "mzmin",
    "mzmax",
    "rt",
    "rtmin",
    "rtmax"
  )

  missing_peak_columns <- setdiff(
    required_peak_columns,
    colnames(chrom_peaks)
  )

  if (length(missing_peak_columns)) {
    stop(
      "'chromPeaks(x)' is missing required columns: ",
      paste(missing_peak_columns, collapse = ", "),
      call. = FALSE
    )
  }

  chrom_peak_data <- tryCatch(
    as.data.frame(xcms::chromPeakData(x)),
    error = function(e) NULL
  )

  if (is.null(chrom_peak_data)) {
    stop(
      "Could not extract chromPeakData(x).",
      call. = FALSE
    )
  }

  if (nrow(chrom_peak_data) != nrow(chrom_peaks)) {
    stop(
      "chromPeakData(x) and chromPeaks(x) are not aligned.",
      call. = FALSE
    )
  }

  sample_index <- .get_chrom_peak_sample_idx(chrom_peaks)

  if (length(sample_index) != nrow(chrom_peaks)) {
    stop(
      "Could not assign one sample index to each chromatographic peak.",
      call. = FALSE
    )
  }

  if (anyNA(sample_index)) {
    stop(
      "Some chromatographic peaks have missing sample indices.",
      call. = FALSE
    )
  }

  if (any(sample_index < 1L | sample_index > n_samples)) {
    stop(
      paste(
        "Sample indices in chromPeaks(x) are not compatible",
        "with the number of rows in sampleData(x)."
      ),
      call. = FALSE
    )
  }

  ms_level_column <- intersect(
    c("ms_level", "msLevel", "mslevel"),
    colnames(chrom_peak_data)
  )

  if (length(ms_level_column)) {
    ms_level <- as.integer(
      chrom_peak_data[[ms_level_column[[1L]]]]
    )

    if (!any(ms_level == 1L, na.rm = TRUE)) {
      stop(
        paste(
          "'x' does not contain MS1 chromatographic peaks.",
          "Run MS1 peak detection before annotateLipids()."
        ),
        call. = FALSE
      )
    }
  }

  acquisition <- .validate_and_normalize_acquisition_mode(
    sample_data = sample_data
  )

  idx_dda <- which(acquisition == "DDA")
  idx_dia <- which(acquisition == "DIA")
  idx_ms2 <- sort(unique(c(idx_dda, idx_dia)))

  if (!length(idx_ms2)) {
    stop(
      paste(
        "No DDA, DIA or SWATH samples were found.",
        "At least one MS2-capable sample is required for lipid annotation."
      ),
      call. = FALSE
    )
  }

  sample_names <- .get_sample_names(x)

  if (length(sample_names) != n_samples) {
    stop(
      paste(
        "The number of sample names does not match",
        "the number of rows in sampleData(x)."
      ),
      call. = FALSE
    )
  }

  has_grouped_features <- .has_valid_grouped_features(x)

  invisible(list(
    acquisition = acquisition,
    sample_names = sample_names,
    idx_dda = idx_dda,
    idx_dia = idx_dia,
    idx_ms2 = idx_ms2,
    has_grouped_features = has_grouped_features
  ))
}

#' Validate and normalize acquisition modes
#'
#' Validate the acquisition mode stored for each sample and normalize accepted
#' values to `"FULLSCAN"`, `"DDA"` or `"DIA"`.
#'
#' @param sample_data `data.frame` obtained from `sampleData(x)`.
#'
#' @return Character vector with one normalized acquisition mode per sample.
#'
#' @noRd
#'
#' @keywords internal
.validate_and_normalize_acquisition_mode <- function(sample_data) {
  if (!is.data.frame(sample_data) || !nrow(sample_data)) {
    stop(
      "'sample_data' must be a non-empty data.frame.",
      call. = FALSE
    )
  }

  if (!"acquisition_mode" %in% colnames(sample_data)) {
    stop(
      paste(
        "sampleData(x) must contain column 'acquisition_mode'.",
        "Accepted values include FULLSCAN/MS, DDA, DIA and SWATH."
      ),
      call. = FALSE
    )
  }

  acquisition_raw <- as.character(sample_data$acquisition_mode)
  acquisition_raw <- trimws(acquisition_raw)

  invalid_empty <- is.na(acquisition_raw) | !nzchar(acquisition_raw)

  if (any(invalid_empty)) {
    stop(
      "'sampleData(x)$acquisition_mode' contains missing or empty values.",
      call. = FALSE
    )
  }

  acquisition <- toupper(acquisition_raw)

  acquisition_compact <- gsub(
    pattern = "[[:space:]_-]+",
    replacement = "",
    x = acquisition
  )

  acquisition[
    acquisition_compact %in% c(
      "FULLSCAN",
      "FULLSCAN/MS",
      "FULLSCANMS",
      "FULLMS",
      "MS",
      "MS1"
    )
  ] <- "FULLSCAN"

  acquisition[
    acquisition_compact %in% c(
      "DDA",
      "DATADEPENDENT",
      "DATADEPENDENTACQUISITION"
    )
  ] <- "DDA"

  acquisition[
    acquisition_compact %in% c(
      "DIA",
      "SWATH",
      "SWATHMS",
      "DATAINDEPENDENT",
      "DATAINDEPENDENTACQUISITION"
    )
  ] <- "DIA"

  allowed <- c("FULLSCAN", "DDA", "DIA")
  invalid <- !acquisition %in% allowed

  if (any(invalid)) {
    invalid_values <- unique(acquisition_raw[invalid])

    stop(
      "Unsupported acquisition_mode value",
      if (length(invalid_values) > 1L) "s" else "",
      ": ",
      paste(invalid_values, collapse = ", "),
      ". Accepted modes include FULLSCAN/MS, DDA, DIA and SWATH.",
      call. = FALSE
    )
  }

  acquisition
}

#' Validate lipid annotation parameters
#'
#' @param paramannot Named annotation parameter list.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @noRd
#'
#' @keywords internal
.validate_paramannot <- function(paramannot) {
  if (!is.list(paramannot)) {
    stop("'paramannot' must be a named list.", call. = FALSE)
  }

  required <- c(
    "polarity",
    "dmzprecursor",
    "dmzproducts",
    "rttol",
    "isotope_mzabs",
    "isotope_maxcharge",
    "isotope_maxiso",
    "isotope_minfrac",
    "isotope_rttol",
    "intval",
    "min_ms2_intensity",
    "remove_isotopes"
  )

  missing <- setdiff(required, names(paramannot))

  if (length(missing)) {
    stop(
      "'paramannot' is missing required elements: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  polarity <- tolower(as.character(paramannot$polarity))

  if (
    length(polarity) != 1L ||
    is.na(polarity) ||
    !polarity %in% c("positive", "negative")
  ) {
    stop(
      "'paramannot$polarity' must be 'positive' or 'negative'.",
      call. = FALSE
    )
  }

  positive_numeric <- c(
    "dmzprecursor",
    "dmzproducts",
    "rttol",
    "isotope_mzabs",
    "isotope_rttol"
  )

  for (nm in positive_numeric) {
    value <- paramannot[[nm]]

    if (
      !is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value) ||
      value <= 0
    ) {
      stop(
        "'paramannot$", nm,
        "' must be a finite numeric value greater than zero.",
        call. = FALSE
      )
    }
  }

  non_negative_numeric <- "min_ms2_intensity"

  for (nm in non_negative_numeric) {
    value <- paramannot[[nm]]

    if (
      !is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value) ||
      value < 0
    ) {
      stop(
        "'paramannot$", nm,
        "' must be a finite non-negative numeric value.",
        call. = FALSE
      )
    }
  }

  positive_integer <- c(
    "isotope_maxcharge",
    "isotope_maxiso"
  )

  for (nm in positive_integer) {
    value <- paramannot[[nm]]

    if (
      !is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value) ||
      value < 1 ||
      value != as.integer(value)
    ) {
      stop(
        "'paramannot$", nm,
        "' must be a positive integer.",
        call. = FALSE
      )
    }
  }

  minfrac <- paramannot$isotope_minfrac

  if (
    !is.numeric(minfrac) ||
    length(minfrac) != 1L ||
    is.na(minfrac) ||
    !is.finite(minfrac) ||
    minfrac < 0 ||
    minfrac > 1
  ) {
    stop(
      "'paramannot$isotope_minfrac' must be between 0 and 1.",
      call. = FALSE
    )
  }

  intval <- paramannot$intval

  if (
    !is.character(intval) ||
    length(intval) != 1L ||
    is.na(intval) ||
    !intval %in% c("into", "maxo", "intb")
  ) {
    stop(
      paste(
        "'paramannot$intval' must be one of",
        "'into', 'maxo' or 'intb'."
      ),
      call. = FALSE
    )
  }

  remove_isotopes <- paramannot$remove_isotopes

  if (
    !is.logical(remove_isotopes) ||
    length(remove_isotopes) != 1L ||
    is.na(remove_isotopes)
  ) {
    stop(
      "'paramannot$remove_isotopes' must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Validate DIA extraction parameters
#'
#' @param paramdia Named DIA parameter list.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @noRd
#'
#' @keywords internal
.validate_paramdia <- function(paramdia) {
  if (!is.list(paramdia)) {
    stop(
      paste(
        "'paramdia' must be a named list when DIA or SWATH",
        "samples are present."
      ),
      call. = FALSE
    )
  }

  required <- c(
    "findChromPeaksIsolationWindowParam",
    "expandRt",
    "minCor",
    "diffRt",
    "intensity",
    "minint",
    "BPPARAM"
  )

  missing <- setdiff(required, names(paramdia))

  if (length(missing)) {
    stop(
      "'paramdia' is missing required elements: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  cwp <- paramdia$findChromPeaksIsolationWindowParam

  if (is.null(cwp)) {
    stop(
      paste(
        "'paramdia$findChromPeaksIsolationWindowParam'",
        "must not be NULL for DIA/SWATH samples."
      ),
      call. = FALSE
    )
  }

  non_negative <- c("expandRt", "diffRt", "minint")

  for (nm in non_negative) {
    value <- paramdia[[nm]]

    if (
      !is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value) ||
      value < 0
    ) {
      stop(
        "'paramdia$", nm,
        "' must be a finite non-negative numeric value.",
        call. = FALSE
      )
    }
  }

  min_cor <- paramdia$minCor

  if (
    !is.numeric(min_cor) ||
    length(min_cor) != 1L ||
    is.na(min_cor) ||
    !is.finite(min_cor) ||
    min_cor < 0 ||
    min_cor > 1
  ) {
    stop(
      "'paramdia$minCor' must be between 0 and 1.",
      call. = FALSE
    )
  }

  intensity <- paramdia$intensity

  if (
    !is.character(intensity) ||
    length(intensity) != 1L ||
    is.na(intensity) ||
    !intensity %in% c("into", "maxo", "intb")
  ) {
    stop(
      paste(
        "'paramdia$intensity' must be one of",
        "'into', 'maxo' or 'intb'."
      ),
      call. = FALSE
    )
  }

  if (!inherits(paramdia$BPPARAM, "BiocParallelParam")) {
    stop(
      "'paramdia$BPPARAM' must inherit from 'BiocParallelParam'.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Check whether an XcmsExperiment contains valid grouped features
#'
#' Check whether `featureDefinitions(x)` contains grouped chromatographic
#' features with valid peak indices.
#'
#' @param x `XcmsExperiment`.
#'
#' @return `logical(1)`.
#'
#' @noRd
#'
#' @keywords internal
.has_valid_grouped_features <- function(x) {
  if (!inherits(x, "XcmsExperiment")) {
    return(FALSE)
  }

  chrom_peaks <- tryCatch(
    xcms::chromPeaks(x),
    error = function(e) NULL
  )

  if (is.null(chrom_peaks) || !nrow(chrom_peaks)) {
    return(FALSE)
  }

  features <- tryCatch(
    as.data.frame(xcms::featureDefinitions(x)),
    error = function(e) NULL
  )

  if (
    is.null(features) ||
    !nrow(features) ||
    !"peakidx" %in% colnames(features)
  ) {
    return(FALSE)
  }

  valid_peakidx <- vapply(
    features$peakidx,
    function(idx) {
      if (is.null(idx) || !length(idx)) {
        return(FALSE)
      }

      if (!is.numeric(idx) && !is.integer(idx)) {
        return(FALSE)
      }

      idx <- as.integer(idx)

      all(
        !is.na(idx) &
          idx >= 1L &
          idx <= nrow(chrom_peaks)
      )
    },
    logical(1)
  )

  all(valid_peakidx)
}

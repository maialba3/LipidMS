# ==============================================================================
# Export lipid annotation results
# ==============================================================================

#' Write an annotation table to Excel
#'
#' @param x `data.frame`.
#' @param file `character(1)` output path.
#' @param sheet `character(1)` worksheet name.
#'
#' @return Invisibly returns `file`.
#'
#' @noRd
#'
#' @keywords internal
.write_annotation_xlsx <- function(x,
                                   file,
                                   sheet = "annotations") {
  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop(
      "Package 'writexl' is required to export annotation tables to Excel.",
      call. = FALSE
    )
  }

  if (
    !is.character(file) ||
    length(file) != 1L ||
    is.na(file) ||
    !nzchar(file)
  ) {
    stop(
      "'file' must be a non-empty character(1).",
      call. = FALSE
    )
  }

  if (!grepl("\\.xlsx$", file, ignore.case = TRUE)) {
    file <- paste0(file, ".xlsx")
  }

  writexl::write_xlsx(
    stats::setNames(list(x), sheet),
    path = file
  )

  invisible(file)
}

#' Export grouped features with lipid annotations
#'
#' Export a flat feature-level annotation table from an `XcmsExperiment`.
#'
#' Feature definitions are combined with the feature intensity matrix obtained
#' with `featureValues(x, value = "into", filled = TRUE)`.
#'
#' If `paramannot$only_confirmed_monoisotopic` is `TRUE`, only features
#' containing a confirmed monoisotopic chromatographic peak are exported.
#'
#' If grouped features are unavailable, or confirmed monoisotopic features
#' cannot be identified when required, an empty table is returned and,
#' when `file` is provided, an empty Excel file is written.
#'
#' @param x `XcmsExperiment`.
#' @param paramannot Annotation parameter list.
#' @param file Optional `character(1)` output path.
#'
#' @return `data.frame`.
#'
#' @examples
#' data(xcmsexperiment_negative)
#'
#' paramannot <- defaultAnnotationParameters("negative")
#'
#' out <- exportAnnotatedFeatures(
#'     xcmsexperiment_negative,
#'     paramannot = paramannot
#' )
#'
#' head(out)
#'
#' @export
exportAnnotatedFeatures <- function(x,
                                    paramannot,
                                    file = NULL) {
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("Package 'xcms' is required.", call. = FALSE)
  }

  if (!inherits(x, "XcmsExperiment")) {
    stop("'x' must be an XcmsExperiment.", call. = FALSE)
  }

  empty_out <- data.frame(
    feature_id = character(),
    stringsAsFactors = FALSE
  )

  fd <- tryCatch(
    as.data.frame(xcms::featureDefinitions(x)),
    error = function(e) NULL
  )

  if (is.null(fd) || !nrow(fd)) {
    message(
      paste(
        "No grouped features are available.",
        "An empty feature annotation table will be exported."
      )
    )

    if (!is.null(file)) {
      .write_annotation_xlsx(
        empty_out,
        file = file,
        sheet = "features"
      )
    }

    return(empty_out)
  }

  fv <- tryCatch(
    as.data.frame(
      xcms::featureValues(
        x,
        value = "into",
        filled = TRUE
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fv) || nrow(fv) != nrow(fd)) {
    message(
      paste(
        "Feature intensity values could not be retrieved or are not",
        "aligned with featureDefinitions(x).",
        "An empty feature annotation table will be exported."
      )
    )

    if (!is.null(file)) {
      .write_annotation_xlsx(
        empty_out,
        file = file,
        sheet = "features"
      )
    }

    return(empty_out)
  }

  feature_id <- rownames(fd)

  duplicate_cols <- intersect(
    colnames(fd),
    colnames(fv)
  )

  if (length(duplicate_cols)) {
    fv <- fv[
      ,
      setdiff(colnames(fv), duplicate_cols),
      drop = FALSE
    ]
  }

  y <- cbind(
    feature_id = feature_id,
    fd,
    fv
  )

  rownames(y) <- NULL

  only_mono <- isTRUE(
    paramannot$only_confirmed_monoisotopic
  )

  if (only_mono) {
    if (!"has_monoisotopic_peak" %in% colnames(y)) {
      message(
        paste(
          "'only_confirmed_monoisotopic = TRUE', but",
          "'has_monoisotopic_peak' is not available in",
          "featureDefinitions(x). An empty table will be exported."
        )
      )

      y <- empty_out

    } else {
      keep <- !is.na(y$has_monoisotopic_peak) &
        y$has_monoisotopic_peak

      y <- y[keep, , drop = FALSE]
    }
  }

  remove_cols <- intersect(
    c(
      "peakidx",
      "monoisotopic_peak_ids"
    ),
    colnames(y)
  )

  if (length(remove_cols)) {
    y <- y[
      ,
      setdiff(colnames(y), remove_cols),
      drop = FALSE
    ]
  }

  if (!nrow(y)) {
    message(
      "No features fulfilled the export criteria."
    )
  }

  if (!is.null(file)) {
    .write_annotation_xlsx(
      y,
      file = file,
      sheet = "features"
    )
  }

  y
}

#' Export chromatographic peaks with lipid annotations and metadata
#'
#' Export chromatographic peak information together with peak-level metadata
#' and lipid annotations.
#'
#' If `paramannot$only_confirmed_monoisotopic` is `TRUE`, only chromatographic
#' peaks with `isotope_rank == 0` are exported.
#'
#' If chromatographic peaks are unavailable, or confirmed monoisotopic peaks
#' cannot be identified when required, an empty table is returned and,
#' when `file` is provided, an empty Excel file is written.
#'
#' @param x `XcmsExperiment`.
#' @param paramannot Annotation parameter list.
#' @param file Optional `character(1)` output path.
#'
#' @return `data.frame`.
#'
#' @examples
#' data(xcmsexperiment_negative)
#'
#' paramannot <- defaultAnnotationParameters("negative")
#'
#' out <- exportAnnotatedChromPeaks(
#'     xcmsexperiment_negative,
#'     paramannot = paramannot
#' )
#'
#' head(out)
#'
#' @export
exportAnnotatedChromPeaks <- function(x,
                                      paramannot,
                                      file = NULL) {
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("Package 'xcms' is required.", call. = FALSE)
  }

  if (!inherits(x, "XcmsExperiment")) {
    stop("'x' must be an XcmsExperiment.", call. = FALSE)
  }

  empty_out <- data.frame(
    chrom_peak_id = character(),
    stringsAsFactors = FALSE
  )

  peak_data <- tryCatch(
    as.data.frame(xcms::chromPeakData(x)),
    error = function(e) NULL
  )

  peaks <- tryCatch(
    as.data.frame(xcms::chromPeaks(x)),
    error = function(e) NULL
  )

  if (is.null(peaks) || !nrow(peaks)) {
    message(
      paste(
        "No chromatographic peaks are available.",
        "An empty chromatographic peak table will be exported."
      )
    )

    if (!is.null(file)) {
      .write_annotation_xlsx(
        empty_out,
        file = file,
        sheet = "chromPeaks"
      )
    }

    return(empty_out)
  }

  if (
    is.null(peak_data) ||
    nrow(peak_data) != nrow(peaks)
  ) {
    message(
      paste(
        "chromPeakData(x) is unavailable or is not aligned with",
        "chromPeaks(x). An empty chromatographic peak table",
        "will be exported."
      )
    )

    if (!is.null(file)) {
      .write_annotation_xlsx(
        empty_out,
        file = file,
        sheet = "chromPeaks"
      )
    }

    return(empty_out)
  }

  peak_ids <- rownames(peaks)

  if (
    !is.null(rownames(peak_data)) &&
    !identical(
      rownames(peak_data),
      peak_ids
    )
  ) {
    idx <- match(
      peak_ids,
      rownames(peak_data)
    )

    if (anyNA(idx)) {
      message(
        paste(
          "chromPeakData(x) could not be matched to",
          "chromPeaks(x) by chromatographic peak ID.",
          "An empty table will be exported."
        )
      )

      if (!is.null(file)) {
        .write_annotation_xlsx(
          empty_out,
          file = file,
          sheet = "chromPeaks"
        )
      }

      return(empty_out)
    }

    peak_data <- peak_data[
      idx,
      ,
      drop = FALSE
    ]
  }

  duplicate_cols <- intersect(
    colnames(peak_data),
    colnames(peaks)
  )

  if (length(duplicate_cols)) {
    peaks <- peaks[
      ,
      setdiff(
        colnames(peaks),
        duplicate_cols
      ),
      drop = FALSE
    ]
  }

  y <- cbind(
    chrom_peak_id = peak_ids,
    peak_data,
    peaks
  )

  rownames(y) <- NULL

  only_mono <- isTRUE(
    paramannot$only_confirmed_monoisotopic
  )

  if (only_mono) {
    if (!"isotope_rank" %in% colnames(y)) {
      message(
        paste(
          "'only_confirmed_monoisotopic = TRUE', but",
          "'isotope_rank' is not available in chromPeakData(x).",
          "An empty table will be exported."
        )
      )

      y <- empty_out

    } else {
      keep <- !is.na(y$isotope_rank) &
        y$isotope_rank == 0L

      y <- y[keep, , drop = FALSE]
    }
  }

  remove_cols <- intersect(
    "isotope_rank",
    colnames(y)
  )

  if (length(remove_cols)) {
    y <- y[
      ,
      setdiff(
        colnames(y),
        remove_cols
      ),
      drop = FALSE
    ]
  }

  if (nrow(y)) {
    y <- tryCatch(
      .add_sample_metadata_to_peak_table(
        out = y,
        x = x
      ),
      error = function(e) {
        message(
          "Sample metadata could not be added to the exported table."
        )

        y
      }
    )
  }

  if (!nrow(y)) {
    message(
      "No chromatographic peaks fulfilled the export criteria."
    )
  }

  if (!is.null(file)) {
    .write_annotation_xlsx(
      y,
      file = file,
      sheet = "chromPeaks"
    )
  }

  y
}

# ==============================================================================
# Isotope annotation utilities
# ==============================================================================

#' Annotate isotope peaks
#'
#' Annotate isotopologues and store the results in `chromPeakData()`.
#'
#' This function is designed to be called internally by `annotateLipids()`.
#' For `mslevel = 1L`, isotope annotation can be performed using CAMERA.
#' For `mslevel = 2L`, the interface is already prepared but the method is not
#' implemented yet.
#'
#' @param x `XcmsExperiment` or compatible xcms object.
#' @param mslevel `integer(1)`. MS level to annotate. Currently only `1L` is
#' implemented.
#' @param param Isotope annotation parameter list. Usually generated or passed
#' internally from `annotateLipids()`.
#' @param method `character(1)`. Isotope annotation backend. Currently
#' `"camera"` is implemented for MS1.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object controlling
#'   parallel processing.
#'
#' @examples
#' data(xcmsexperiment_negative)
#'
#' param <- defaultAnnotationParameters("negative")
#'
#' x <- annotateIsotopes(
#'     xcmsexperiment_negative,
#'     param = param
#' )
#'
#' head(xcms::chromPeakData(x))
#'
#' @return The input object with isotope annotation columns added to
#' `chromPeakData()`.
#'
#' @export
annotateIsotopes <- function(x,
                             mslevel = 1L,
                             param = list(),
                             method = c("camera"),
                             BPPARAM = BiocParallel::SerialParam()) {
  method <- match.arg(method)

  if (!mslevel %in% c(1L, 2L)) {
    stop("'mslevel' must be either 1L or 2L.", call. = FALSE)
  }

  param <- .prepare_isotope_param(
    param = param,
    mslevel = mslevel
  )

  if (mslevel == 1L) {
    return(.annotate_ms1_isotopes(
      x = x,
      param = param,
      method = method,
      BPPARAM = BPPARAM
    ))
  }

  .annotate_ms2_isotopes(
    x = x,
    param = param,
    method = method
  )
}

#' Prepare isotope annotation parameters
#'
#' Normalize isotope annotation parameters. The mass tolerance in ppm is expected
#' to be inherited from `dmzprecursor` for MS1 and from `dmzproducts` for MS2
#' when called from `annotateLipids()`.
#'
#' @param param Named list of isotope annotation parameters.
#' @param mslevel `integer(1)`. MS level.
#'
#' @return Named list of completed parameters.
#'
#' @noRd
#'
#' @keywords internal
.prepare_isotope_param <- function(param, mslevel = 1L) {
  if (is.null(param)) {
    param <- list()
  }

  if (!is.list(param)) {
    stop("'param' must be a list.", call. = FALSE)
  }

  ppm <- param$ppm
  if (is.null(ppm)) {
    ppm <- if (mslevel == 1L) {
      param$dmzprecursor
    } else {
      param$dmzproducts
    }
  }

  if (is.null(ppm)) {
    ppm <- 5
  }

  list(
    ppm = ppm,
    mzabs = .get_param_value(param, "mzabs", 0.005),
    maxcharge = .get_param_value(param, "maxcharge", 3L),
    maxiso = .get_param_value(param, "maxiso", 4L),
    intval = .get_param_value(param, "intval", "into"),
    minfrac = .get_param_value(param, "minfrac", 0.5),
    rt_tolerance = .get_param_value(param, "rt_tolerance", 1),
    group_fwhm = .get_param_value(param, "group_fwhm", TRUE)
  )
}

#' Get parameter value with default
#'
#' @param param Named list.
#' @param name Parameter name.
#' @param default Default value.
#'
#' @return Parameter value.
#'
#' @noRd
#'
#' @keywords internal
.get_param_value <- function(param, name, default) {
  value <- param[[name]]
  if (is.null(value)) {
    default
  } else {
    value
  }
}

#' Annotate MS1 isotopes
#'
#' Annotate MS1 isotopologues and store CAMERA-derived information in
#' `chromPeakData()`.
#'
#' @param x `XcmsExperiment` or compatible xcms object.
#' @param param Prepared isotope parameter list.
#' @param method `character(1)`. Backend method.
#'
#' @return Object with updated `chromPeakData()`.
#'
#' @noRd
#'
#' @keywords internal
.annotate_ms1_isotopes <- function(x,
                                   param,
                                   method = "camera",
                                   BPPARAM = BiocParallel::SerialParam()) {
  if (!identical(method, "camera")) {
    stop("Only method = 'camera' is implemented for MS1.", call. = FALSE)
  }

  .annotate_ms1_isotopes_camera(
    x = x,
    param = param,
    BPPARAM = BPPARAM
  )
}

#' Annotate MS2 isotopes
#'
#' Placeholder for future MS2 isotope annotation.
#'
#' The function is intentionally present so that `annotateLipids()` can call the
#' same isotope annotation interface for MS1 and MS2. In the future, the ppm
#' tolerance should be inherited from `dmzproducts`.
#'
#' @param x `XcmsExperiment` or compatible xcms object.
#' @param param Prepared isotope parameter list.
#' @param method `character(1)`. Backend method.
#'
#' @return Currently stops because MS2 isotope annotation is not implemented.
#'
#' @noRd
#'
#' @keywords internal
.annotate_ms2_isotopes <- function(x,
                                   param,
                                   method = "camera") {
  stop("MS2 isotope annotation is not implemented yet.", call. = FALSE)
}

#' Annotate MS1 isotopes using CAMERA by sample
#'
#' @param x `XcmsExperiment` or compatible xcms object.
#' @param param Prepared isotope parameter list.
#'
#' @return Object with isotope annotation columns added to `chromPeakData()`.
#'
#' @noRd
#'
#' @keywords internal
.annotate_ms1_isotopes_camera <- function(x,
                                          param,
                                          BPPARAM = BiocParallel::SerialParam()) {
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("Package 'xcms' is required.", call. = FALSE)
  }

  if (!requireNamespace("CAMERA", quietly = TRUE)) {
    stop("Package 'CAMERA' is required for isotope annotation.", call. = FALSE)
  }

  chrom_peaks <- xcms::chromPeaks(x)
  chrom_peak_data <- xcms::chromPeakData(x)

  if (!nrow(chrom_peaks)) {
    stop("'x' does not contain chromatographic peaks.", call. = FALSE)
  }

  peak_ids <- rownames(chrom_peaks)
  if (is.null(peak_ids)) {
    stop("'chromPeaks(x)' must have row names.", call. = FALSE)
  }

  chrom_peak_data <- .initialize_isotope_chrom_peak_data(chrom_peak_data)

  sample_col <- .get_chrom_peak_sample_column(chrom_peaks, chrom_peak_data)
  sample_ids <- sample_col$sample_ids

  sample_index <- split(seq_along(sample_ids), sample_ids)

  sample_names <- names(sample_index)

  isotope_list <- BiocParallel::bplapply(
    sample_names,
    function(sample_name) {
      peak_idx <- sample_index[[sample_name]]

      if (length(peak_idx) < 2L) {
        return(NULL)
      }

      isotope_info <- .annotate_ms1_isotopes_camera_one_sample(
        x = x,
        peak_idx = peak_idx,
        peak_ids = peak_ids,
        chrom_peaks = chrom_peaks,
        param = param
      )

      if (!nrow(isotope_info)) {
        return(NULL)
      }

      isotope_info$sample_name <- sample_name
      isotope_info
    },
    BPPARAM = BPPARAM
  )

  isotope_list <- isotope_list[!vapply(isotope_list, is.null, logical(1))]

  for (isotope_info in isotope_list) {
    sample_name <- isotope_info$sample_name[[1L]]
    isotope_info$sample_name <- NULL

    idx <- match(isotope_info$peak_id, peak_ids)
    keep <- !is.na(idx)

    idx <- idx[keep]
    isotope_info <- isotope_info[keep, , drop = FALSE]

    if (!length(idx)) {
      next
    }

    chrom_peak_data$is_isotope[idx] <- isotope_info$is_isotope
    chrom_peak_data$isotope_group[idx] <- paste(
      sample_name,
      isotope_info$isotope_group,
      sep = ":"
    )
    chrom_peak_data$isotope_rank[idx] <- isotope_info$isotope_rank
    chrom_peak_data$isotope_charge[idx] <- isotope_info$isotope_charge
    chrom_peak_data$monoisotopic_peak_id[idx] <-
      isotope_info$monoisotopic_peak_id
    chrom_peak_data$monoisotopic_mz[idx] <- isotope_info$monoisotopic_mz
    chrom_peak_data$isotope_annotation[idx] <-
      isotope_info$isotope_annotation
    chrom_peak_data$isotope_annotation_source[idx] <- "CAMERA"
  }

  xcms::chromPeakData(x) <- as.data.frame(chrom_peak_data)
  x
}

#' Initialize isotope columns in chromPeakData
#'
#' @param chrom_peak_data `DataFrame` returned by `chromPeakData()`.
#'
#' @return Updated `DataFrame`.
#'
#' @noRd
#'
#' @keywords internal
.initialize_isotope_chrom_peak_data <- function(chrom_peak_data) {
  chrom_peak_data$is_isotope <- FALSE
  chrom_peak_data$isotope_group <- NA_character_
  chrom_peak_data$isotope_rank <- NA_integer_
  chrom_peak_data$isotope_charge <- NA_integer_
  chrom_peak_data$monoisotopic_peak_id <- NA_character_
  chrom_peak_data$monoisotopic_mz <- NA_real_
  chrom_peak_data$isotope_annotation <- NA_character_
  chrom_peak_data$isotope_annotation_source <- NA_character_

  chrom_peak_data
}

#' Convert xcms object to CAMERA xsAnnotate
#'
#' Convert an `XcmsExperiment` or `xcmsSet` to a CAMERA `xsAnnotate`
#' object suitable for per-sample isotope annotation.
#'
#' Group information inherited from xcms is removed when the object contains
#' a single sample. CAMERA isotope annotation is performed on chromatographic
#' peaks and pseudospectra generated internally with `groupFWHM()`, so existing
#' xcms feature grouping must not be used as CAMERA group information.
#'
#' @param x `XcmsExperiment`, `xcmsSet`, or `xsAnnotate`.
#'
#' @return CAMERA `xsAnnotate` object.
#'
#' @noRd
#'
#' @keywords internal
.as_camera_xsannotate <- function(x) {
  if (inherits(x, "xsAnnotate")) {
    return(x)
  }

  if (inherits(x, "xcmsSet")) {
    xcms_set <- x
  } else {
    xcms_set <- tryCatch(
      methods::as(x, "xcmsSet"),
      error = function(e) NULL
    )

    if (is.null(xcms_set)) {
      stop(
        paste(
          "Could not coerce 'x' to 'xcmsSet'.",
          "CAMERA requires an object compatible with xsAnnotate."
        ),
        call. = FALSE
      )
    }
  }

  n_samples <- length(xcms::sampnames(xcms_set))

  if (n_samples == 1L) {
    xcms_set@groups <- matrix(
      numeric(0),
      nrow = 0L,
      ncol = ncol(xcms_set@groups),
      dimnames = list(
        NULL,
        colnames(xcms_set@groups)
      )
    )

    xcms_set@groupidx <- list()
  }

  CAMERA::xsAnnotate(
    xs = xcms_set,
    sample = if (n_samples == 1L) 1L else NA
  )
}

#' Parse CAMERA isotope annotations
#'
#' @param peaklist CAMERA peak list.
#' @param peak_ids Chromatographic peak identifiers.
#' @param chrom_peaks Matrix returned by `chromPeaks()`.
#' @param ppm Mass tolerance in ppm.
#' @param rt_tolerance Retention time tolerance.
#'
#' @return `data.frame` with parsed isotope annotations.
#'
#' @noRd
#'
#' @keywords internal
.parse_camera_isotopes <- function(peaklist,
                                   peak_ids,
                                   chrom_peaks,
                                   ppm = 5,
                                   rt_tolerance = 1) {
  iso <- peaklist$isotopes
  has_iso <- !is.na(iso) & nzchar(iso)

  if (!any(has_iso)) {
    return(.empty_isotope_table())
  }

  peaklist_idx <- which(has_iso)

  out <- lapply(peaklist_idx, function(i) {
    annotation <- iso[[i]]
    parsed <- .parse_camera_isotope_string(annotation)

    peak_id <- .map_camera_peak_to_chrom_peak(
      peaklist_row = peaklist[i, , drop = FALSE],
      peak_ids = peak_ids,
      chrom_peaks = chrom_peaks,
      ppm = ppm,
      rt_tolerance = rt_tolerance
    )

    if (is.na(peak_id)) {
      return(NULL)
    }

    is_mono <- !is.na(parsed$isotope_rank) && parsed$isotope_rank == 0L

    data.frame(
      peak_id = peak_id,
      is_isotope = !is_mono,
      isotope_group = parsed$isotope_group,
      isotope_rank = parsed$isotope_rank,
      isotope_charge = parsed$isotope_charge,
      monoisotopic_peak_id = if (is_mono) peak_id else NA_character_,
      monoisotopic_mz = if (is_mono) chrom_peaks[peak_id, "mz"] else NA_real_,
      isotope_annotation = annotation,
      stringsAsFactors = FALSE
    )
  })

  out <- out[!vapply(out, is.null, logical(1))]

  if (!length(out)) {
    return(.empty_isotope_table())
  }

  out <- do.call(rbind, out)
  .fill_camera_monoisotopic_links(out, chrom_peaks)
}

#' Return empty isotope annotation table
#'
#' @return Empty `data.frame` with isotope annotation columns.
#'
#' @noRd
#'
#' @keywords internal
.empty_isotope_table <- function() {
  data.frame(
    peak_id = character(),
    is_isotope = logical(),
    isotope_group = character(),
    isotope_rank = integer(),
    isotope_charge = integer(),
    monoisotopic_peak_id = character(),
    monoisotopic_mz = numeric(),
    isotope_annotation = character(),
    stringsAsFactors = FALSE
  )
}

#' Parse a CAMERA isotope annotation string
#'
#' @param x CAMERA isotope annotation.
#'
#' @return Named list with isotope group, rank and charge.
#'
#' @noRd
#'
#' @keywords internal
.parse_camera_isotope_string <- function(x) {
  isotope_group <- NA_character_
  isotope_rank <- NA_integer_
  isotope_charge <- NA_integer_

  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    return(list(
      isotope_group = isotope_group,
      isotope_rank = isotope_rank,
      isotope_charge = isotope_charge
    ))
  }

  group_match <- regexpr("\\[[^]]+\\]", x, perl = TRUE)
  if (group_match[1L] != -1L) {
    isotope_group <- regmatches(x, group_match)
    isotope_group <- gsub("^\\[|\\]$", "", isotope_group)
  }

  rank_match <- regexpr("M\\+[0-9]+", x, perl = TRUE)
  if (rank_match[1L] != -1L) {
    isotope_rank <- as.integer(sub("M\\+", "", regmatches(x, rank_match)))
  } else if (grepl("\\[M\\]", x, perl = TRUE) || grepl("^M($|[^+])", x)) {
    isotope_rank <- 0L
  }

  charge_match <- regexpr("[0-9]+\\+$", x, perl = TRUE)
  if (charge_match[1L] != -1L) {
    isotope_charge <- as.integer(sub("\\+$", "", regmatches(x, charge_match)))
  }

  list(
    isotope_group = isotope_group,
    isotope_rank = isotope_rank,
    isotope_charge = isotope_charge
  )
}

#' Map CAMERA peak list row to chromPeak ID
#'
#' @param peaklist_row One-row CAMERA peak list.
#' @param peak_ids Chromatographic peak identifiers.
#' @param chrom_peaks Matrix returned by `chromPeaks()`.
#' @param ppm Mass tolerance in ppm.
#' @param rt_tolerance Retention time tolerance.
#'
#' @return `character(1)` chromPeak ID or `NA_character_`.
#'
#' @noRd
#'
#' @keywords internal
.map_camera_peak_to_chrom_peak <- function(peaklist_row,
                                           peak_ids,
                                           chrom_peaks,
                                           ppm = 5,
                                           rt_tolerance = 1) {
  if ("peakidx" %in% colnames(peaklist_row)) {
    idx <- as.integer(peaklist_row$peakidx[[1]])
    if (!is.na(idx) && idx >= 1L && idx <= length(peak_ids)) {
      return(peak_ids[[idx]])
    }
  }

  if (!all(c("mz", "rt") %in% colnames(peaklist_row))) {
    return(NA_character_)
  }

  mz <- as.numeric(peaklist_row$mz[[1]])
  rt <- as.numeric(peaklist_row$rt[[1]])

  mz_tol <- mz * ppm / 1e6

  candidates <- which(
    abs(chrom_peaks[, "mz"] - mz) <= mz_tol &
      abs(chrom_peaks[, "rt"] - rt) <= rt_tolerance
  )

  if (!length(candidates)) {
    return(NA_character_)
  }

  candidates <- candidates[order(
    abs(chrom_peaks[candidates, "mz"] - mz) +
      abs(chrom_peaks[candidates, "rt"] - rt)
  )]

  peak_ids[[candidates[[1L]]]]
}

#' Fill monoisotopic links within CAMERA isotope groups
#'
#' @param isotope_info Parsed isotope annotation table.
#' @param chrom_peaks Matrix returned by `chromPeaks()`.
#'
#' @return Updated isotope annotation table.
#'
#' @noRd
#'
#' @keywords internal
.fill_camera_monoisotopic_links <- function(isotope_info, chrom_peaks) {
  groups <- unique(stats::na.omit(isotope_info$isotope_group))

  for (grp in groups) {
    idx <- which(isotope_info$isotope_group == grp)

    if (!length(idx)) {
      next
    }

    mono_idx <- idx[which.min(isotope_info$isotope_rank[idx])]
    mono_peak_id <- isotope_info$peak_id[[mono_idx]]
    mono_mz <- chrom_peaks[mono_peak_id, "mz"]

    isotope_info$monoisotopic_peak_id[idx] <- mono_peak_id
    isotope_info$monoisotopic_mz[idx] <- mono_mz
  }

  isotope_info
}

#' Get sample identifier for each chromPeak
#'
#' @param chrom_peaks Matrix returned by `chromPeaks()`.
#' @param chrom_peak_data `DataFrame` returned by `chromPeakData()`.
#'
#' @return List with element `sample_ids`.
#'
#' @noRd
#'
#' @keywords internal
.get_chrom_peak_sample_column <- function(chrom_peaks, chrom_peak_data) {
  if ("sample" %in% colnames(chrom_peak_data)) {
    return(list(sample_ids = as.character(chrom_peak_data$sample)))
  }

  if ("sample" %in% colnames(chrom_peaks)) {
    return(list(sample_ids = as.character(chrom_peaks[, "sample"])))
  }

  if ("sampleidx" %in% colnames(chrom_peaks)) {
    return(list(sample_ids = as.character(chrom_peaks[, "sampleidx"])))
  }

  stop(
    paste(
      "Could not determine sample for each chromPeak.",
      "Expected 'sample' in chromPeakData() or 'sample'/'sampleidx'",
      "in chromPeaks()."
    ),
    call. = FALSE
  )
}

#' Annotate MS1 isotopes in one sample using CAMERA
#'
#' @param x Original xcms object.
#' @param peak_idx Integer indices of chromPeaks belonging to one sample.
#' @param peak_ids Chromatographic peak identifiers.
#' @param chrom_peaks Full chromPeaks matrix.
#' @param param Prepared isotope parameter list.
#'
#' @return Parsed isotope annotation table.
#'
#' @noRd
#'
#' @keywords internal
.annotate_ms1_isotopes_camera_one_sample <- function(
    x,
    peak_idx,
    peak_ids,
    chrom_peaks,
    param
) {
  n_samples <- nrow(
    as.data.frame(MsExperiment::sampleData(x))
  )

  if (n_samples == 1L) {
    x_sample <- x
  } else {
    x_sample <- .subset_xcms_object_to_chrom_peaks(
      x = x,
      peak_idx = peak_idx
    )
  }

  camera_object <- .as_camera_xsannotate(
    x_sample
  )

  if (isTRUE(param$group_fwhm)) {
    camera_object <- CAMERA::groupFWHM(
      camera_object,
      intval = param$intval
    )
  }

  camera_object <- CAMERA::findIsotopes(
    camera_object,
    maxcharge = param$maxcharge,
    maxiso = param$maxiso,
    ppm = param$ppm,
    mzabs = param$mzabs,
    intval = param$intval,
    minfrac = param$minfrac
  )

  peaklist <- CAMERA::getPeaklist(
    camera_object,
    intval = param$intval
  )

  if (!"isotopes" %in% colnames(peaklist)) {
    return(.empty_isotope_table())
  }

  sample_peak_ids <- peak_ids[peak_idx]

  sample_chrom_peaks <- chrom_peaks[
    peak_idx,
    ,
    drop = FALSE
  ]

  .parse_camera_isotopes(
    peaklist = peaklist,
    peak_ids = sample_peak_ids,
    chrom_peaks = sample_chrom_peaks,
    ppm = param$ppm,
    rt_tolerance = param$rt_tolerance
  )
}

#' Subset xcms object to selected chromPeaks
#'
#' @param x xcms object.
#' @param peak_idx Integer indices of chromPeaks.
#'
#' @return xcms object subset to selected chromPeaks when possible.
#'
#' @noRd
#'
#' @keywords internal
.subset_xcms_object_to_chrom_peaks <- function(x, peak_idx) {
  peak_ids <- rownames(xcms::chromPeaks(x))[peak_idx]

  out <- tryCatch(
    xcms::filterChromPeaks(x, keep = peak_ids),
    error = function(e) NULL
  )

  if (!is.null(out)) {
    return(out)
  }

  stop(
    paste(
      "Could not subset object by chromPeaks.",
      "A per-sample CAMERA annotation requires filterChromPeaks()."
    ),
    call. = FALSE
  )
}

# .regroup_annotation_subset_for_camera <- function(x_subset, paramannot, verbose = TRUE) {
#   n_samples <- length(.get_sample_names(x_subset))
#
#   if (n_samples <= 1L) {
#     if (isTRUE(verbose)) {
#       cat("Annotation subset contains one sample; regrouping is not required for CAMERA.\n")
#     }
#     return(x_subset)
#   }
#
#   bw <- .get_param_value(paramannot, "annotation_subset_group_bw", 5)
#   min_fraction <- .get_param_value(paramannot, "annotation_subset_group_min_fraction", 0.5)
#   bin_size <- .get_param_value(paramannot, "annotation_subset_group_bin_size", 0.01)
#
#   if (isTRUE(verbose)) {
#     cat(
#       "Regrouping annotation subset for CAMERA/xsAnnotate compatibility: ",
#       "sampleGroups=1 for ", n_samples, " samples; ",
#       "bw=", bw, "; minFraction=", min_fraction, "; binSize=", bin_size, "\n",
#       sep = ""
#     )
#   }
#
#   xcms::groupChromPeaks(
#     x_subset,
#     param = xcms::PeakDensityParam(
#       sampleGroups = rep(1L, n_samples),
#       bw = bw,
#       minFraction = min_fraction,
#       binSize = bin_size
#     )
#   )
# }

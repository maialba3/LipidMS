# ==============================================================================
# MS2 spectra extraction utilities for lipid annotation
# ==============================================================================

#' Extract MS2 spectra for lipid annotation
#'
#' Extract MS2 information per sample and per MS1 chromatographic peak.
#'
#' For DDA samples, MS2 spectra are matched to MS1 chromPeaks using precursor
#' m/z tolerance and retention time limits. If multiple MS2 scans match a peak,
#' scans within `rttol` of the MS1 apex are retained.
#'
#' For DIA/SWATH samples, MS2 chromatographic peaks are detected with
#' `findChromPeaksIsolationWindow()` and spectra are reconstructed with
#' `reconstructChromPeakSpectra()`. Reconstructed fragment signals are linked
#' back to MS2 chromPeak IDs when possible.
#'
#' @param x `XcmsExperiment`.
#' @param requested_peak_ids Optional character vector with chromatographic
#'   peak identifiers for which MS2 spectra should be extracted. If `NULL`,
#'   all chromatographic peaks are considered.
#' @param paramannot Named list with annotation parameters. Required elements:
#' `dmzprecursor`, `dmzproducts`, and `rttol`.
#' @param paramdia Named list with DIA extraction parameters.
#' @param BPPARAM BiocParallel parameter.
#' @param verbose `logical(1)`. Whether to print progress information.
#'
#' @examples
#' data(xcmsexperiment_negative)
#'
#' param <- defaultAnnotationParameters("negative")
#' dbs <- createLipidDbs()
#'
#' candidates <- findMs1LipidCandidates(
#'     xcmsexperiment_negative,
#'     dbs,
#'     param
#' )
#'
#' ms2_data <- extractMs2Spectra(
#'     x = xcmsexperiment_negative,
#'     requested_peak_ids = unique(candidates$chrom_peak_id),
#'     paramannot = param,
#'     paramdia = list(),
#'     verbose = FALSE
#' )
#'
#' names(ms2_data)
#'
#' @return Named list with elements:
#' - `samples`: per-sample MS2 spectra information.
#'
#' @export
extractMs2Spectra <- function(x,
                              requested_peak_ids = NULL,
                              paramannot,
                              paramdia = list(),
                              BPPARAM = BiocParallel::SerialParam(),
                              verbose = TRUE) {
  .check_ms2_extraction_inputs(x, paramannot)

  acquisition <- .get_acquisition_mode(x)
  sample_names <- .get_sample_names(x)

  out <- vector("list", length(sample_names))
  names(out) <- sample_names

  idx_dda <- which(acquisition == "DDA")
  idx_dia <- which(acquisition %in% c("DIA", "SWATH"))

  if (isTRUE(verbose)) {
    message(
      "extractMs2Spectra: samples = ", length(sample_names),
      "; DDA = ", length(idx_dda),
      "; DIA/SWATH = ", length(idx_dia)
    )
    message(
      "extractMs2Spectra: requested_peak_ids = ",
      if (is.null(requested_peak_ids)) "ALL" else length(requested_peak_ids)
    )
  }

  # DDA: subset the XcmsExperiment to DDA samples before matching already
  # acquired MS2 spectra. This reduces the object copied to workers.
  if (length(idx_dda)) {
    BPPARAM_DDA <- .cap_bpparam_workers(BPPARAM, length(idx_dda))

    BPPARAM_DDA <- .snow_to_multicore_for_sourced_helpers(
      BPPARAM_DDA,
      verbose = verbose
    )

    if (isTRUE(verbose)) {
      message(
        "extractMs2Spectra: extracting DDA spectra from ", length(idx_dda),
        " sample(s); workers = ", .bpparam_nworkers_safe(BPPARAM_DDA)
      )
    }

    x_dda <- .subset_xcms_experiment_samples(x, idx_dda)
    dda_sample_names <- sample_names[idx_dda]

    dda_res <- BiocParallel::bplapply(
      seq_along(idx_dda),
      function(j) {
        .extract_dda_sample_spectra(
          x = x_dda,
          sample_idx = j,
          sample_name = dda_sample_names[[j]],
          requested_peak_ids = requested_peak_ids,
          dmzprecursor = paramannot$dmzprecursor,
          rttol = paramannot$rttol,
          min_intensity = .get_param_value(paramannot, "min_ms2_intensity", 0)
        )
      },
      BPPARAM = BPPARAM_DDA
    )

    for (j in seq_along(idx_dda)) {
      out[[dda_sample_names[[j]]]] <- dda_res[[j]]
    }

    rm(x_dda)
    invisible(gc())
  }

  # DIA/SWATH: subset to DIA/SWATH samples before detecting MS2 chromPeaks.
  # The original per-sample reconstruction logic is kept unchanged.
  if (length(idx_dia)) {
    # BPPARAM_DIA <- .cap_bpparam_workers(BPPARAM, length(idx_dia))
    BPPARAM_DIA <- .cap_bpparam_workers(
      .get_param_value(paramdia, "BPPARAM", BPPARAM),
      length(idx_dia)
    )
    BPPARAM_DIA <- .snow_to_multicore_for_sourced_helpers(
      BPPARAM_DIA,
      verbose = verbose
    )
    paramdia_dia <- paramdia
    paramdia_dia$BPPARAM <- BiocParallel::SerialParam()

    if (isTRUE(verbose)) {
      message(
        "extractMs2Spectra: extracting DIA/SWATH spectra from ", length(idx_dia),
        " sample(s); reconstructChromPeakSpectra() workers = serial; outer DIA loop = ",
        .bpparam_nworkers_safe(BPPARAM_DIA)
      )
    }

    x_dia <- .subset_xcms_experiment_samples(x, idx_dia)

    x_dia <- .prepare_dia_chrom_peaks(
      x = x_dia,
      paramdia = paramdia,
      BPPARAM = BPPARAM_DIA
    )

    dia_sample_names <- sample_names[idx_dia]

    dia_res <- BiocParallel::bplapply(
      seq_along(idx_dia),
      function(j) {
        .extract_dia_sample_spectra(
          x_dia = x_dia,
          sample_idx = j,
          sample_name = dia_sample_names[[j]],
          requested_peak_ids = requested_peak_ids,
          dmzproducts = paramannot$dmzproducts,
          rttol = paramannot$rttol,
          paramdia = paramdia,
          BPPARAM = paramdia_dia$BPPARAM
        )
      },
      BPPARAM = BPPARAM_DIA
    )

    for (j in seq_along(idx_dia)) {
      out[[dia_sample_names[[j]]]] <- dia_res[[j]]
    }

    rm(x_dia)
    invisible(gc())
  }

  if (isTRUE(verbose)) {
    n_spectra <- sum(vapply(out, function(z) {
      if (is.null(z) || is.null(z$spectra)) {
        return(0L)
      }
      sum(vapply(z$spectra, function(sp) {
        is.data.frame(sp) && nrow(sp) > 0L
      }, logical(1)))
    }, integer(1)))

    message(
      "extractMs2Spectra: finished; MS1 peaks with extracted/reconstructed MS2 = ",
      n_spectra
    )
  }

  list(samples = out)
}

# ------------------------------------------------------------------------------
# Validation and metadata helpers
# ------------------------------------------------------------------------------

.check_ms2_extraction_inputs <- function(x, paramannot) {
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("Package 'xcms' is required.", call. = FALSE)
  }

  if (!requireNamespace("Spectra", quietly = TRUE)) {
    stop("Package 'Spectra' is required.", call. = FALSE)
  }

  if (!inherits(x, "XcmsExperiment")) {
    stop("'x' must be an XcmsExperiment.", call. = FALSE)
  }

  cp <- xcms::chromPeaks(x)
  if (!nrow(cp)) {
    stop("'x' must contain MS1 chromPeaks.", call. = FALSE)
  }

  required <- c("dmzprecursor", "dmzproducts", "rttol")
  missing <- setdiff(required, names(paramannot))
  if (length(missing)) {
    stop(
      "Missing required elements in 'paramannot': ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.get_param_value <- function(x, name, default) {
  value <- x[[name]]
  if (is.null(value)) default else value
}

.get_acquisition_mode <- function(x) {
  sample_data <- .get_sample_data(x)

  .validate_and_normalize_acquisition_mode(sample_data)
}

.get_sample_data <- function(x) {
  if (!requireNamespace("MsExperiment", quietly = TRUE)) {
    stop("Package 'MsExperiment' is required.", call. = FALSE)
  }

  as.data.frame(MsExperiment::sampleData(x))
}

.get_sample_names <- function(x) {
  sd <- .get_sample_data(x)

  if ("sample" %in% colnames(sd)) {
    return(as.character(sd$sample))
  }

  basename(rownames(sd))
}

.get_chrom_peak_sample_idx <- function(chrom_peaks) {
  if ("sample" %in% colnames(chrom_peaks)) {
    return(as.integer(chrom_peaks[, "sample"]))
  }

  if ("sampleidx" %in% colnames(chrom_peaks)) {
    return(as.integer(chrom_peaks[, "sampleidx"]))
  }

  stop(
    "Could not find sample index in chromPeaks().",
    call. = FALSE
  )
}

.empty_ms2_spectrum_table <- function() {
  data.frame(
    ms1_chrom_peak_id = character(),
    ms2_spectrum_id = character(),
    ms2_chrom_peak_id = character(),
    mz = numeric(),
    intensity = numeric(),
    rtime = numeric(),
    stringsAsFactors = FALSE
  )
}

.get_spectra <- function(x) {
  if (!requireNamespace("MsExperiment", quietly = TRUE)) {
    stop("Package 'MsExperiment' is required.", call. = FALSE)
  }

  MsExperiment::spectra(x)
}

#' Limit BiocParallel workers to the number of tasks
#'
#' The user-selected backend is preserved whenever possible. On Windows, or when
#' using SnowParam from a sourced script, non-exported helper functions may not be
#' visible to remote workers; in that case we fall back to SerialParam for the
#' outer LipidMS loops to keep behaviour robust across systems.
#'
#' @param BPPARAM BiocParallel parameter supplied by the user.
#' @param n Number of independent tasks.
#'
#' @return BiocParallel parameter with an appropriate number of workers.
#'
#' @noRd
#'
#' @keywords internal
.cap_bpparam_workers <- function(BPPARAM, n) {
  n <- max(1L, as.integer(n))

  current_workers <- tryCatch(
    BiocParallel::bpnworkers(BPPARAM),
    error = function(e) n
  )

  workers <- min(current_workers, n)

  if (workers <= 1L || inherits(BPPARAM, "SerialParam")) {
    return(BiocParallel::SerialParam())
  }

  if (inherits(BPPARAM, "MulticoreParam")) {
    return(BiocParallel::MulticoreParam(workers = workers))
  }

  if (inherits(BPPARAM, "SnowParam")) {
    return(BiocParallel::SnowParam(workers = workers))
  }

  BPPARAM
}

#' Safely report the number of BiocParallel workers
#'
#' @param BPPARAM BiocParallel parameter.
#'
#' @return Integer worker count or NA.
#'
#' @noRd
#'
#' @keywords internal
.bpparam_nworkers_safe <- function(BPPARAM) {
  tryCatch(
    BiocParallel::bpnworkers(BPPARAM),
    error = function(e) NA_integer_
  )
}

.snow_to_multicore_for_sourced_helpers <- function(BPPARAM, verbose = TRUE) {
  if (!inherits(BPPARAM, "SnowParam")) {
    return(BPPARAM)
  }

  workers <- .bpparam_nworkers_safe(BPPARAM)
  if (is.na(workers) || workers < 1L) {
    workers <- 1L
  }

  if (.Platform$OS.type == "windows") {
    if (isTRUE(verbose)) {
      message(
        "extractMs2Spectra: SnowParam with sourced helpers detected on Windows; ",
        "using SerialParam to avoid missing helper functions."
      )
    }
    return(BiocParallel::SerialParam())
  }

  if (isTRUE(verbose)) {
    message(
      "extractMs2Spectra: converting SnowParam to MulticoreParam(workers = ",
      workers,
      ") for sourced helper visibility."
    )
  }

  BiocParallel::MulticoreParam(workers = workers)
}



# ------------------------------------------------------------------------------
# DDA extraction
# ------------------------------------------------------------------------------

.extract_dda_sample_spectra <- function(x,
                                        sample_idx,
                                        sample_name,
                                        requested_peak_ids,
                                        dmzprecursor,
                                        rttol,
                                        min_intensity = 0) {
  cp <- xcms::chromPeaks(x)
  sample_index <- .get_chrom_peak_sample_idx(cp)

  keep_ms1 <- sample_index == sample_idx
  cp_ms1 <- cp[keep_ms1, , drop = FALSE]
  if (!is.null(requested_peak_ids)) {
    cp_ms1 <- cp_ms1[rownames(cp_ms1) %in% requested_peak_ids, , drop = FALSE]
  }

  if (!nrow(cp_ms1)) {
    return(list(
      sample = sample_name,
      acquisition = "DDA",
      spectra = list()
    ))
  }

  sps <- .get_spectra(x)
  sps <- .filter_spectra_to_sample(sps, sample_idx = sample_idx,
                                   sample_name = sample_name)
  sps <- Spectra::filterMsLevel(sps, msLevel = 2L)

  spectra_by_peak <- stats::setNames(
    replicate(nrow(cp_ms1), .empty_ms2_spectrum_table(), simplify = FALSE),
    rownames(cp_ms1)
  )

  if (!length(sps)) {
    return(list(
      sample = sample_name,
      acquisition = "DDA",
      spectra = spectra_by_peak
    ))
  }

  precursor_mz <- .get_spectra_precursor_mz(sps)
  spec_rt <- Spectra::rtime(sps)
  spec_ids <- .get_spectrum_ids(sps)

  for (i in seq_len(nrow(cp_ms1))) {
    peak_id <- rownames(cp_ms1)[[i]]

    mz <- cp_ms1[i, "mz"]
    rt <- cp_ms1[i, "rt"]
    rtmin <- cp_ms1[i, "rtmin"]
    rtmax <- cp_ms1[i, "rtmax"]

    mz_tol <- mz * dmzprecursor / 1e6

    keep <- !is.na(precursor_mz) &
      abs(precursor_mz - mz) <= mz_tol &
      spec_rt >= rtmin &
      spec_rt <= rtmax &
      abs(spec_rt - rt) <= rttol

    if (!any(keep)) {
      next
    }

    spectra_by_peak[[peak_id]] <- .spectra_to_fragment_table(
      sps = sps[keep],
      ms1_chrom_peak_id = peak_id,
      ms2_chrom_peak_id = NA_character_,
      min_intensity = min_intensity
    )
  }

  list(
    sample = sample_name,
    acquisition = "DDA",
    spectra = spectra_by_peak
  )
}

.filter_spectra_to_sample <- function(sps, sample_idx, sample_name) {
  if ("sample_id" %in% Spectra::spectraVariables(sps)) {
    sample_id <- Spectra::spectraData(sps)$sample_id
    keep <- sample_id == sample_idx | sample_id == sample_name
    if (any(keep, na.rm = TRUE)) {
      return(sps[keep])
    }
  }

  data_origin <- tryCatch(
    basename(Spectra::dataOrigin(sps)),
    error = function(e) NULL
  )

  if (!is.null(data_origin)) {
    keep <- data_origin == basename(sample_name)
    if (any(keep, na.rm = TRUE)) {
      return(sps[keep])
    }
  }

  sps
}

.get_spectra_precursor_mz <- function(sps) {
  vars <- Spectra::spectraVariables(sps)

  if ("precursorMz" %in% vars) {
    return(Spectra::precursorMz(sps))
  }

  if ("isolationWindowTargetMz" %in% vars) {
    return(Spectra::spectraData(sps)$isolationWindowTargetMz)
  }

  rep(NA_real_, length(sps))
}

.get_spectrum_ids <- function(sps) {
  vars <- Spectra::spectraVariables(sps)

  if ("spectrumId" %in% vars) {
    return(as.character(Spectra::spectraData(sps)$spectrumId))
  }

  if ("acquisitionNum" %in% vars) {
    return(as.character(Spectra::spectraData(sps)$acquisitionNum))
  }

  if ("scanIndex" %in% vars) {
    return(as.character(Spectra::spectraData(sps)$scanIndex))
  }

  as.character(seq_along(sps))
}

#' Convert MS2 spectra to a fragment table
#'
#' @param sps `Spectra` object containing MS2 spectra.
#' @param ms1_chrom_peak_id `character(1)`. Linked MS1 chromatographic peak ID.
#' @param ms2_chrom_peak_id Optional MS2 chromatographic peak ID.
#' @param min_intensity Minimum fragment intensity.
#'
#' @return `data.frame` with fragment m/z, intensity and peak identifiers.
#'
#' @noRd
#'
#' @keywords internal
.spectra_to_fragment_table <- function(sps,
                                       ms1_chrom_peak_id = NA_character_,
                                       ms2_chrom_peak_id = NA_character_,
                                       min_intensity = 0) {
  if (!inherits(sps, "Spectra")) {
    stop("'sps' must be a Spectra object.", call. = FALSE)
  }

  if (!length(sps)) {
    return(data.frame(
      mz = numeric(),
      intensity = numeric(),
      chrom_peak_id = character(),
      ms1_chrom_peak_id = character(),
      ms2_chrom_peak_id = character(),
      spectrum_id = character(),
      stringsAsFactors = FALSE
    ))
  }

  pd <- Spectra::peaksData(sps)
  out <- vector("list", length(pd))

  spectrum_id <- seq_along(pd)

  if ("acquisitionNum" %in% Spectra::spectraVariables(sps)) {
    spectrum_id <- Spectra::acquisitionNum(sps)
  }

  for (i in seq_along(pd)) {
    p <- pd[[i]]

    if (is.null(p) || !nrow(p)) {
      next
    }

    cn <- colnames(p)
    cn_low <- tolower(cn)

    mz_col <- cn[match("mz", cn_low)]
    int_col <- if ("intensity" %in% cn_low) {
      cn[match("intensity", cn_low)]
    } else if ("int" %in% cn_low) {
      cn[match("int", cn_low)]
    } else {
      cn[min(2L, ncol(p))]
    }

    frag <- data.frame(
      mz = as.numeric(p[, mz_col]),
      intensity = as.numeric(p[, int_col]),
      chrom_peak_id = ms1_chrom_peak_id,
      ms1_chrom_peak_id = ms1_chrom_peak_id,
      ms2_chrom_peak_id = ms2_chrom_peak_id,
      spectrum_id = as.character(spectrum_id[[i]]),
      stringsAsFactors = FALSE
    )

    frag <- frag[is.finite(frag$mz) &
                   is.finite(frag$intensity) &
                   frag$intensity >= min_intensity, ,
                 drop = FALSE]

    out[[i]] <- frag
  }

  out <- out[!vapply(out, is.null, logical(1))]

  if (!length(out)) {
    return(data.frame(
      mz = numeric(),
      intensity = numeric(),
      chrom_peak_id = character(),
      ms1_chrom_peak_id = character(),
      ms2_chrom_peak_id = character(),
      spectrum_id = character(),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------------------------
# DIA extraction
# ------------------------------------------------------------------------------

#' Subset an XcmsExperiment by sample index
#'
#' @param x `XcmsExperiment`.
#' @param sample_idx Integer sample indices.
#'
#' @return Subset `XcmsExperiment`.
#'
#' @noRd
#'
#' @keywords internal
.subset_xcms_experiment_samples <- function(x, sample_idx) {
  if (!length(sample_idx)) {
    stop("'sample_idx' must contain at least one sample index.", call. = FALSE)
  }

  out <- tryCatch(
    x[, sample_idx],
    error = function(e) NULL
  )

  if (!is.null(out)) {
    return(out)
  }

  out <- tryCatch(
    x[sample_idx],
    error = function(e) NULL
  )

  if (!is.null(out)) {
    return(out)
  }

  stop("Could not subset 'x' by sample index.", call. = FALSE)
}

.prepare_dia_chrom_peaks <- function(x, paramdia, BPPARAM) {
  if (!"findChromPeaksIsolationWindow" %in% getNamespaceExports("xcms")) {
    stop(
      "xcms::findChromPeaksIsolationWindow() is required for DIA extraction.",
      call. = FALSE
    )
  }

  param <- paramdia$findChromPeaksIsolationWindowParam
  if (is.null(param)) {
    param <- paramdia$param
  }

  if (is.null(param)) {
    stop(
      paste(
        "DIA extraction requires",
        "'paramdia$findChromPeaksIsolationWindowParam' or 'paramdia$param'."
      ),
      call. = FALSE
    )
  }

  xcms::findChromPeaksIsolationWindow(
    x,
    param = param,
    BPPARAM = BPPARAM
  )
}

.extract_dia_sample_spectra <- function(x_dia,
                                        sample_idx,
                                        sample_name,
                                        requested_peak_ids,
                                        dmzproducts,
                                        rttol,
                                        paramdia,
                                        BPPARAM =  BiocParallel::SerialParam()) {
  cp <- xcms::chromPeaks(x_dia)
  cpd <- as.data.frame(xcms::chromPeakData(x_dia))
  sample_index <- .get_chrom_peak_sample_idx(cp)

  ms_level <- .get_chrom_peak_ms_level(cpd)

  keep_ms1 <- sample_index == sample_idx & ms_level == 1L
  keep_ms2 <- sample_index == sample_idx & ms_level == 2L

  cp_ms1 <- cp[keep_ms1, , drop = FALSE]
  if (!is.null(requested_peak_ids)) {
    cp_ms1 <- cp_ms1[rownames(cp_ms1) %in% requested_peak_ids, , drop = FALSE]
  }
  cp_ms2 <- cp[keep_ms2, , drop = FALSE]
  cpd_ms2 <- cpd[keep_ms2, , drop = FALSE]

  spectra_by_peak <- stats::setNames(
    replicate(nrow(cp_ms1), .empty_ms2_spectrum_table(), simplify = FALSE),
    rownames(cp_ms1)
  )

  if (!nrow(cp_ms1)) {
    return(list(
      sample = sample_name,
      acquisition = "DIA",
      spectra = spectra_by_peak
    ))
  }

  requested_peak_ids <- rownames(cp_ms1)

  spectra_dia <- tryCatch(
    xcms::reconstructChromPeakSpectra(
      x_dia,
      expandRt = .get_param_value(paramdia, "expandRt", 1),
      minCor = .get_param_value(paramdia, "minCor", 0.8),
      diffRt = .get_param_value(paramdia, "diffRt", 3),
      intensity = .get_param_value(paramdia, "intensity", "into"),
      peakId = requested_peak_ids,
      BPPARAM = BPPARAM
    ),
    error = function(e) NULL
  )

  if (is.null(spectra_dia)) {
    return(list(
      sample = sample_name,
      acquisition = "DIA",
      spectra = spectra_by_peak
    ))
  }

  returned_peak_ids <- .extract_reconstructed_peak_ids(
    spectra_obj = spectra_dia,
    requested_peak_ids = requested_peak_ids
  )

  pd <- Spectra::peaksData(spectra_dia)
  n_use <- min(length(pd), length(returned_peak_ids))

  ms2_peak_table <- .build_ms2_peak_table(
    cp_ms2 = cp_ms2,
    cpd_ms2 = cpd_ms2
  )

  for (i in seq_len(n_use)) {
    pid <- returned_peak_ids[[i]]

    if (is.na(pid) || !pid %in% requested_peak_ids) {
      next
    }

    sp <- pd[[i]]
    if (is.null(sp) || !nrow(sp)) {
      next
    }

    frag_table <- .reconstructed_peaks_to_fragment_table(
      peaks = sp,
      ms1_chrom_peak_id = pid,
      ms1_mz = cp_ms1[pid, "mz"],
      ms1_rt = cp_ms1[pid, "rt"],
      ms2_peak_table = ms2_peak_table,
      dmzproducts = dmzproducts,
      rttol = rttol,
      min_intensity = .get_param_value(paramdia, "minint", 0)
    )

    spectra_by_peak[[pid]] <- frag_table
  }

  list(
    sample = sample_name,
    acquisition = "DIA",
    spectra = spectra_by_peak
  )
}

.get_chrom_peak_ms_level <- function(chrom_peak_data) {
  candidate <- c("ms_level", "msLevel", "mslevel")
  hit <- intersect(candidate, colnames(chrom_peak_data))

  if (!length(hit)) {
    rep(1L, nrow(chrom_peak_data))
  } else {
    as.integer(chrom_peak_data[[hit[[1L]]]])
  }
}

.build_ms2_peak_table <- function(cp_ms2, cpd_ms2) {
  if (!nrow(cp_ms2)) {
    return(data.frame(
      ms2_chrom_peak_id = character(),
      mz = numeric(),
      mzmin = numeric(),
      mzmax = numeric(),
      rtime = numeric(),
      rtmin = numeric(),
      rtmax = numeric(),
      isolationWindowTargetMz = numeric(),
      isolationWindowLowerMz = numeric(),
      isolationWindowUpperMz = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  out <- data.frame(
    ms2_chrom_peak_id = rownames(cp_ms2),
    mz = cp_ms2[, "mz"],
    mzmin = cp_ms2[, "mzmin"],
    mzmax = cp_ms2[, "mzmax"],
    rtime = cp_ms2[, "rt"],
    rtmin = cp_ms2[, "rtmin"],
    rtmax = cp_ms2[, "rtmax"],
    stringsAsFactors = FALSE
  )

  out$isolationWindowTargetMz <- .get_optional_col(
    cpd_ms2, "isolationWindowTargetMz", NA_real_
  )
  out$isolationWindowLowerMz <- .get_optional_col(
    cpd_ms2, "isolationWindowLowerMz", NA_real_
  )
  out$isolationWindowUpperMz <- .get_optional_col(
    cpd_ms2, "isolationWindowUpperMz", NA_real_
  )

  out
}

.get_optional_col <- function(x, name, default) {
  if (name %in% colnames(x)) {
    x[[name]]
  } else {
    rep(default, nrow(x))
  }
}

.extract_reconstructed_peak_ids <- function(spectra_obj, requested_peak_ids) {
  sd <- tryCatch(
    as.data.frame(Spectra::spectraData(spectra_obj)),
    error = function(e) NULL
  )

  candidate_cols <- c(
    "peak_id", "peakId", "chrom_peak_id", "chromPeakId",
    "chromPeakID", "ms1_peak_id", "ms1PeakId",
    "chrom_peak_ids", "chromPeakIds"
  )

  if (!is.null(sd)) {
    hit <- intersect(candidate_cols, colnames(sd))

    if (length(hit)) {
      ids <- sd[[hit[[1L]]]]

      if (is.list(ids)) {
        ids <- vapply(ids, function(z) {
          if (!length(z) || all(is.na(z))) {
            NA_character_
          } else {
            as.character(z[[1L]])
          }
        }, character(1))
      } else {
        ids <- as.character(ids)
      }

      return(ids)
    }
  }

  pd <- Spectra::peaksData(spectra_obj)
  if (length(pd) == length(requested_peak_ids)) {
    return(requested_peak_ids)
  }

  stop("Unable to map reconstructed spectra back to MS1 chromPeak IDs.",
       call. = FALSE)
}

.reconstructed_peaks_to_fragment_table <- function(peaks,
                                                   ms1_chrom_peak_id,
                                                   ms1_mz,
                                                   ms1_rt,
                                                   ms2_peak_table,
                                                   dmzproducts,
                                                   rttol,
                                                   min_intensity = 0) {
  peaks <- as.data.frame(peaks)

  cn <- colnames(peaks)
  cn_low <- tolower(cn)

  mz_col <- cn[match("mz", cn_low)]
  int_col <- if ("intensity" %in% cn_low) {
    cn[match("intensity", cn_low)]
  } else if ("int" %in% cn_low) {
    cn[match("int", cn_low)]
  } else {
    cn[min(2L, ncol(peaks))]
  }

  rt_col <- if ("rt" %in% cn_low) {
    cn[match("rt", cn_low)]
  } else {
    NA_character_
  }

  frag_mz <- as.numeric(peaks[[mz_col]])
  frag_int <- as.numeric(peaks[[int_col]])
  frag_rt <- if (!is.na(rt_col)) {
    as.numeric(peaks[[rt_col]])
  } else {
    rep(ms1_rt, length(frag_mz))
  }

  ms2_chrom_peak_id <- vapply(
    seq_along(frag_mz),
    function(i) {
      .match_dia_fragment_to_ms2_peak(
        frag_mz = frag_mz[[i]],
        frag_rt = frag_rt[[i]],
        ms1_mz = ms1_mz,
        ms2_peak_table = ms2_peak_table,
        ppm = dmzproducts,
        rt_tolerance = rttol
      )
    },
    character(1)
  )

  out <- data.frame(
    ms1_chrom_peak_id = ms1_chrom_peak_id,
    ms2_spectrum_id = ms1_chrom_peak_id,
    ms2_chrom_peak_id = ms2_chrom_peak_id,
    mz = frag_mz,
    intensity = frag_int,
    rtime = frag_rt,
    stringsAsFactors = FALSE
  )

  out <- out[out$intensity > min_intensity, , drop = FALSE]

  if (!nrow(out)) {
    return(.empty_ms2_spectrum_table())
  }

  rownames(out) <- NULL
  out
}

.match_dia_fragment_to_ms2_peak <- function(frag_mz,
                                            frag_rt,
                                            ms1_mz,
                                            ms2_peak_table,
                                            ppm,
                                            rt_tolerance) {
  if (!nrow(ms2_peak_table)) {
    return(NA_character_)
  }

  keep <- abs(ms2_peak_table$rtime - frag_rt) <= rt_tolerance

  has_iw <- all(c(
    "isolationWindowLowerMz",
    "isolationWindowUpperMz"
  ) %in% colnames(ms2_peak_table))

  if (has_iw) {
    keep <- keep &
      (
        is.na(ms2_peak_table$isolationWindowLowerMz) |
          is.na(ms2_peak_table$isolationWindowUpperMz) |
          (
            ms2_peak_table$isolationWindowLowerMz <= ms1_mz &
              ms2_peak_table$isolationWindowUpperMz >= ms1_mz
          )
      )
  }

  cand <- ms2_peak_table[keep, , drop = FALSE]

  if (!nrow(cand)) {
    return(NA_character_)
  }

  mz_tol <- frag_mz * ppm / 1e6

  keep_mz <- cand$mzmin <= frag_mz + mz_tol &
    cand$mzmax >= frag_mz - mz_tol

  if (any(keep_mz)) {
    cand <- cand[keep_mz, , drop = FALSE]
  }

  score <- abs(cand$mz - frag_mz) / max(mz_tol, 1e-6) +
    abs(cand$rtime - frag_rt) / max(rt_tolerance, 1e-6)

  cand$ms2_chrom_peak_id[[which.min(score)]]
}

.add_dia_ms2_links_to_chrom_peak_data <- function(x_dia, samples) {
  cpd <- as.data.frame(xcms::chromPeakData(x_dia))

  if (!"linked_ms1_chrom_peak_id" %in% colnames(cpd)) {
    cpd$linked_ms1_chrom_peak_id <- NA_character_
  }

  all_links <- unlist(lapply(samples, function(sample_res) {
    if (is.null(sample_res) || !identical(sample_res$acquisition, "DIA")) {
      return(character())
    }

    spec <- sample_res$spectra
    unlist(lapply(spec, function(z) {
      if (is.null(z) || !nrow(z)) {
        return(character())
      }

      ids <- z$ms2_chrom_peak_id
      ids <- ids[!is.na(ids) & nzchar(ids)]
      stats::setNames(z$ms1_chrom_peak_id[match(ids, z$ms2_chrom_peak_id)], ids)
    }))
  }))

  if (length(all_links)) {
    idx <- match(names(all_links), rownames(xcms::chromPeaks(x_dia)))
    keep <- !is.na(idx)
    cpd$linked_ms1_chrom_peak_id[idx[keep]] <- unname(all_links[keep])
  }

  xcms::chromPeakData(x_dia) <- as.data.frame(cpd)
  x_dia
}

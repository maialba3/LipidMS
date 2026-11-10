# ==============================================================================
# Organize annotation results
# ==============================================================================

#' Add accepted MS2 lipid annotations to chromPeakData
#'
#' Add compact lipid annotation results obtained after MS2 fragment-rule
#' matching to the `chromPeakData()` slot of an `XcmsExperiment`.
#'
#' Only annotations with `status == "accepted"` are exported.
#'
#' Each accepted annotation is propagated to all chromatographic peaks belonging
#' to the corresponding `adduct_group_id`, lipid class and sum composition.
#'
#' The original adduct assigned to each chromatographic peak is preserved. This
#' means that peaks belonging to the same `adduct_group_id` keep their own
#' `adduct` value instead of receiving a collapsed adduct string.
#'
#' If several lipid annotations are associated with the same chromatographic
#' peak, values are collapsed into pipe-separated (`"|"`) strings after
#' filtering redundant annotations and sorting them by confidence.
#'
#' The following columns are added to `chromPeakData()`:
#' - `has_lipid_ms2_annotation`
#' - `n_lipid_ms2_annotations`
#' - `lipid_class`
#' - `sum_composition`
#' - `carbon`
#' - `double_bond`
#' - `oxygen`
#' - `adduct`
#' - `adduct_group_id`
#' - `annotation_level`
#' - `annotation`
#'
#' @param x `XcmsExperiment`.
#' @param ms2_annotation Named list returned by [matchFragmentRules()].
#' Must contain element `all_annotations`.
#' @param candidates Candidate table returned by [groupCoelutingAdducts()].
#'
#' @examples
#' data(xcmsexperiment_negative)
#'
#' param <- defaultAnnotationParameters("negative")
#' dbs <- createLipidDbs()
#' rules <- defaultFragmentRules()
#'
#' candidates <- findMs1LipidCandidates(
#'     xcmsexperiment_negative,
#'     dbs,
#'     param
#' )
#' candidates <- groupCoelutingAdducts(candidates, param)
#'
#' ms2_data <- extractMs2Spectra(
#'     xcmsexperiment_negative,
#'     requested_peak_ids = unique(candidates$chrom_peak_id),
#'     paramannot = param,
#'     paramdia = list(),
#'     verbose = FALSE
#' )
#'
#' ms2_annotation <- matchFragmentRules(
#'     candidates,
#'     ms2_data,
#'     rules,
#'     dbs,
#'     param
#' )
#'
#' x <- addMs2AnnotationsToChromPeakData(
#'     xcmsexperiment_negative,
#'     ms2_annotation,
#'     candidates
#' )
#'
#' head(xcms::chromPeakData(x))
#'
#' @return Updated `XcmsExperiment` with lipid annotation columns added to
#' `chromPeakData()`.
#'
#' @export
addMs2AnnotationsToChromPeakData <- function(x, ms2_annotation, candidates) {
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("Package 'xcms' is required.", call. = FALSE)
  }

  cpd <- as.data.frame(xcms::chromPeakData(x))
  peak_ids <- rownames(xcms::chromPeaks(x))

  out_cols <- c(
    "lipid_class", "sum_composition", "carbon", "double_bond", "oxygen",
    "adduct", "adduct_group_id", "adduct_group_score",
    "adduct_group_n_adducts", "annotation_level", "annotation"
  )

  for (col in out_cols) {
    cpd[[col]] <- NA_character_
  }

  cpd$has_lipid_ms2_annotation <- FALSE
  cpd$n_lipid_ms2_annotations <- 0L

  if (
    is.null(ms2_annotation) ||
    is.null(ms2_annotation$all_annotations) ||
    !nrow(ms2_annotation$all_annotations) ||
    !nrow(candidates)
  ) {
    xcms::chromPeakData(x) <- as.data.frame(cpd)
    return(x)
  }

  ann <- ms2_annotation$all_annotations
  ann <- ann[ann$status == "accepted", , drop = FALSE]

  if (!nrow(ann)) {
    xcms::chromPeakData(x) <- as.data.frame(cpd)
    return(x)
  }

  required_ann <- c(
    "adduct_group_id", "lipid_class", "sum_composition",
    "annotation_level", "annotation", "status"
  )

  required_cand <- c(
    "chrom_peak_id", "adduct_group_id", "class", "sum_composition",
    "carbon", "double_bond", "oxygen", "adduct"
  )

  missing_ann <- setdiff(required_ann, colnames(ann))
  missing_cand <- setdiff(required_cand, colnames(candidates))

  if (length(missing_ann)) {
    stop(
      "'ms2_annotation$all_annotations' is missing columns: ",
      paste(missing_ann, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(missing_cand)) {
    stop(
      "'candidates' is missing columns: ",
      paste(missing_cand, collapse = ", "),
      call. = FALSE
    )
  }

  peak_candidates <- .prepare_peak_level_annotation_candidates(candidates)

  ann_expanded <- merge(
    ann[, c(
      "adduct_group_id", "lipid_class", "sum_composition",
      "annotation_level", "annotation"
    ), drop = FALSE],
    peak_candidates,
    by.x = c("adduct_group_id", "lipid_class", "sum_composition"),
    by.y = c("adduct_group_id", "class", "sum_composition"),
    all = FALSE,
    sort = FALSE
  )

  if (!nrow(ann_expanded)) {
    xcms::chromPeakData(x) <- as.data.frame(cpd)
    return(x)
  }

  ann_expanded <- unique(ann_expanded)

  by_peak <- split(ann_expanded, ann_expanded$chrom_peak_id)

  for (peak_id in names(by_peak)) {
    peak_idx <- match(peak_id, peak_ids)

    if (is.na(peak_idx)) {
      next
    }

    tab <- by_peak[[peak_id]]
    tab <- .keep_best_annotations_per_sum_composition(tab)
    tab <- .sort_ms2_annotations_for_chrom_peak(tab)

    cpd$has_lipid_ms2_annotation[peak_idx] <- TRUE
    cpd$n_lipid_ms2_annotations[peak_idx] <- nrow(tab)

    for (col in out_cols) {
      cpd[[col]][peak_idx] <- paste(tab[[col]], collapse = "|")
    }
  }

  xcms::chromPeakData(x) <- as.data.frame(cpd)
  x
}

#' Prepare peak-level candidates for annotation export
#'
#' Build one compact candidate row per chromatographic peak, adduct group, lipid
#' class and sum composition, while preserving the original peak-level adduct.
#'
#' The number of unique adducts observed in each `adduct_group_id`, lipid class
#' and sum composition is calculated and propagated to every chromatographic
#' peak in the group as `adduct_group_n_adducts`.
#'
#' Unlike group-level summaries, this helper does not collapse the `adduct`
#' column. Each chromatographic peak keeps its own adduct assignment so that
#' the adduct-to-peak relationship remains explicit.
#'
#' @param candidates Candidate table returned by [groupCoelutingAdducts()].
#'
#' @return `data.frame` with one row per unique combination of:
#' `chrom_peak_id`, `adduct_group_id`, `class`, `sum_composition`, and `adduct`.
#'
#' The output contains:
#' - `chrom_peak_id`
#' - `adduct_group_id`
#' - `class`
#' - `sum_composition`
#' - `carbon`
#' - `double_bond`
#' - `oxygen`
#' - `adduct`
#' - `adduct_group_score`
#' - `adduct_group_n_adducts`
#'
#' @noRd
#'
#' @keywords internal
.prepare_peak_level_annotation_candidates <- function(candidates) {
  if (!nrow(candidates)) {
    return(data.frame(
      chrom_peak_id = character(),
      adduct_group_id = character(),
      class = character(),
      sum_composition = character(),
      carbon = numeric(),
      double_bond = numeric(),
      oxygen = numeric(),
      adduct = character(),
      adduct_group_score = numeric(),
      adduct_group_n_adducts = integer(),
      stringsAsFactors = FALSE
    ))
  }

  if (!"adduct_group_n_adducts" %in% colnames(candidates)) {
    group_key <- paste(
      candidates$adduct_group_id,
      candidates$class,
      candidates$sum_composition,
      sep = "|"
    )

    candidates$adduct_group_n_adducts <- stats::ave(
      candidates$adduct,
      group_key,
      FUN = function(z) length(unique(z))
    )
  }

  if (!"adduct_group_score" %in% colnames(candidates)) {
    if ("adduct_score" %in% colnames(candidates)) {
      group_key <- paste(
        candidates$adduct_group_id,
        candidates$class,
        candidates$sum_composition,
        sep = "|"
      )

      candidates$adduct_group_score <- round(stats::ave(
        candidates$adduct_score,
        group_key,
        FUN = function(z) sum(z, na.rm = TRUE),3)
      )
    } else {
      candidates$adduct_group_score <- NA_real_
    }
  }

  out <- unique(candidates[, c(
    "chrom_peak_id", "adduct_group_id", "class", "sum_composition",
    "carbon", "double_bond", "oxygen", "adduct",
    "adduct_group_score", "adduct_group_n_adducts"
  ), drop = FALSE])

  rownames(out) <- NULL
  out
}

#' Sort MS2 annotations for compact peak-level export
#'
#' Sort annotation rows before collapsing them into pipe-separated strings
#' for storage in `chromPeakData()`.
#'
#' Annotations are ordered by:
#'
#' 1. Highest annotation confidence level:
#'    `FA position > FA > Subclass > MS-only > rejected`
#' 2. Highest number of supporting adducts within the adduct group.
#' 3. Lipid class name.
#' 4. Sum composition.
#'
#' This ensures that the most reliable annotations appear first when several
#' annotations are associated with the same chromatographic peak.
#'
#' @param x `data.frame` with annotation rows. Must contain columns
#' `annotation_level` and `adduct_group_n_adducts`.
#'
#' @return Sorted `data.frame`.
#'
#' @noRd
#'
#' @keywords internal
.sort_ms2_annotations_for_chrom_peak <- function(x) {
  if (!nrow(x)) {
    return(x)
  }

  if (!"adduct_group_n_adducts" %in% colnames(x)) {
    x$adduct_group_n_adducts <- 1L
  }

  level_rank <- c(
    rejected = 0L,
    `MS-only` = 1L,
    Subclass = 2L,
    FA = 3L,
    `FA position` = 4L
  )

  rank <- unname(level_rank[x$annotation_level])
  rank[is.na(rank)] <- 0L

  x <- x[
    order(
      -rank,
      -x$adduct_group_n_adducts,
      x$lipid_class,
      x$sum_composition
    ),
    ,
    drop = FALSE
  ]

  rownames(x) <- NULL
  x
}


#' Keep best annotation per class and sum composition
#'
#' Collapse redundant annotations for the same lipid class and total
#' composition.
#'
#' For each `lipid_class` and `sum_composition`, only the annotation with the
#' highest `annotation_level` is kept. If the best level is `"FA position"`,
#' all unique positional annotations at that level are kept.
#'
#' This means that duplicate lower-confidence annotations for the same total
#' composition are removed, while alternative positional assignments such as
#' `PC(16:0/18:1)` and `PC(18:1/16:0)` are retained when both are supported at
#' `"FA position"` level.
#'
#' @param x `data.frame` with annotation rows. Must contain columns
#' `lipid_class`, `sum_composition`, `annotation_level`, and `annotation`.
#'
#' @return Filtered `data.frame`.
#'
#' @noRd
#'
#' @keywords internal
.keep_best_annotations_per_sum_composition <- function(x) {
  if (!nrow(x)) {
    return(x)
  }

  required <- c(
    "lipid_class", "sum_composition", "adduct_group_score", "adduct_group_n_adducts", "annotation_level", "annotation"
  )

  missing <- setdiff(required, colnames(x))
  if (length(missing)) {
    stop(
      "'x' is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  level_rank <- c(
    rejected = 0L,
    `MS-only` = 1L,
    Subclass = 2L,
    FA = 3L,
    `FA position` = 4L
  )

  x$annotation_level_rank <- unname(level_rank[x$annotation_level])
  x$annotation_level_rank[is.na(x$annotation_level_rank)] <- 0L

  split_key <- paste(
    x$lipid_class,
    x$sum_composition,
    sep = "|"
  )

  groups <- split(seq_len(nrow(x)), split_key)

  out <- lapply(groups, function(idx) {
    tab <- x[idx, , drop = FALSE]

    best_rank <- max(tab$annotation_level_rank, na.rm = TRUE)
    tab <- tab[tab$annotation_level_rank == best_rank, , drop = FALSE]

    if (best_rank < level_rank[["FA position"]]) {
      tab <- tab[1L, , drop = FALSE]
    } else {
      tab <- tab[!duplicated(tab$annotation), , drop = FALSE]
    }

    tab
  })

  out <- do.call(rbind, out)
  out$annotation_level_rank <- NULL
  rownames(out) <- NULL
  out
}

#' Export chromatographic peaks with lipid annotations and metadata
#'
#' Build a flat peak-level annotation table from an `XcmsExperiment`.
#'
#' The output contains one row per chromatographic peak, combining columns from
#' `chromPeaks(x)` and `chromPeakData(x)`. If sample metadata are available in
#' `phenoData(x)`, they are appended using the peak sample index.
#'
#' @param x `XcmsExperiment`.
#' @param file Optional `character(1)` output path. If `NULL`, no file is
#' written.
#' @param sep Field separator used when `file` is provided.
#'
#' @return `data.frame` with chromatographic peak information, peak metadata and
#' lipid annotation columns.
#'
#' @noRd
#'
#' @keywords internal
.exportAnnotatedChromPeaks <- function(x, file = NULL, sep = ",") {
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("Package 'xcms' is required.", call. = FALSE)
  }

  peaks <- as.data.frame(xcms::chromPeaks(x))
  peak_data <- as.data.frame(xcms::chromPeakData(x))

  if (nrow(peaks) != nrow(peak_data)) {
    stop("chromPeaks and chromPeakData are not aligned.", call. = FALSE)
  }

  peaks$chrom_peak_id <- rownames(peaks)
  peak_data$chrom_peak_id <- rownames(peak_data)

  out <- merge(
    peaks,
    peak_data,
    by = "chrom_peak_id",
    all.x = TRUE,
    sort = FALSE
  )

  out <- .add_sample_metadata_to_peak_table(out = out, x = x)

  if (!is.null(file)) {
    .write_annotation_table(out, file = file, sep = sep)
  }

  out
}

#' Store LipidMS auxiliary objects inside an XcmsExperiment
#'
#' Store LipidMS annotation metadata and auxiliary objects inside the
#' `metadata(x)` slot of an `XcmsExperiment`.
#'
#' This function provides a centralized location for storing intermediate and
#' final LipidMS annotation results together with processing parameters and
#' package metadata, enabling reproducible and self-contained annotated
#' `XcmsExperiment` objects.
#'
#' Stored objects are placed inside:
#'
#' `metadata(x)$LipidMS`
#'
#' Typical stored components include:
#'
#' - MS1 lipid candidates
#' - reconstructed or extracted MS2 data
#' - annotation parameters
#' - fragmentation rules
#' - package version information
#' - processing timestamps
#'
#' @param x `XcmsExperiment`.
#' @param ms1_candidates Optional `data.frame` containing MS1 lipid candidates.
#' @param ms2_data Optional object containing MS2 spectra, reconstructed DIA
#' spectra or fragment information.
#' @param ms2_annotation MS2 annotation result returned by
#'   [matchFragmentRules()].
#' @param params Optional `list` with annotation parameters.
#' @param rules Optional `list` with fragmentation or annotation rules.
#' @param version Optional `character(1)` specifying the LipidMS version used.
#' If `NULL`, the installed package version is used when available.
#'
#' @return
#' Updated `XcmsExperiment` with LipidMS information stored in
#' `metadata(x)$LipidMS`.
#'
#' @details
#' This function is intended to preserve annotation provenance and intermediate
#' results within the `XcmsExperiment` object itself.
#'
#' Peak-level annotations should generally be stored in:
#'
#' `chromPeakData(x)`
#'
#' while feature-level summaries should be stored in:
#'
#' `featureDefinitions(x)`
#'
#' Large or hierarchical objects such as reconstructed MS2 spectra are better
#' stored in `metadata(x)` rather than duplicated across chromatographic peak
#' tables.
#'
#' The resulting object can be serialized with:
#'
#' `saveRDS(x, file)`
#'
#' preserving the complete annotation workflow state.
#'
#' @seealso
#' [addMs2AnnotationsToChromPeakData()],
#' [addFeatureLipidAnnotations()]
#'
#' @examples
#' data(xcmsexperiment_negative)
#'
#' param <- defaultAnnotationParameters("negative")
#'
#' x <- addLipidMSMetadata(
#'     xcmsexperiment_negative,
#'     params = list(annotation = param)
#' )
#'
#' S4Vectors::metadata(x)$LipidMS
#'
#' @export
addLipidMSMetadata <- function(x,
                                ms1_candidates = NULL,
                                ms2_data = NULL,
                                ms2_annotation = NULL,
                                params = NULL,
                                rules = NULL,
                                version = NULL) {
  md <- S4Vectors::metadata(x)

  if (is.null(md$LipidMS)) {
    md$LipidMS <- list()
  }

  if (!is.null(ms1_candidates)) {
    md$LipidMS$ms1_candidates <- ms1_candidates
  }

  if (!is.null(ms2_data)) {
    md$LipidMS$ms2_data <- ms2_data
  }

  if (!is.null(ms2_annotation)) {
    md$LipidMS$ms2_annotation <- ms2_annotation
  }

  if (!is.null(params)) {
    md$LipidMS$params <- params
  }

  if (!is.null(rules)) {
    md$LipidMS$rules <- rules
  }

  if (is.null(version)) {
    version <- tryCatch(
      as.character(utils::packageVersion("LipidMS")),
      error = function(e) NA_character_
    )
  }

  md$LipidMS$version <- version
  md$LipidMS$timestamp <- Sys.time()

  S4Vectors::metadata(x) <- md
  x
}

#' Add lipid and isotope annotation summaries to featureDefinitions
#'
#' Propagate peak-level lipid annotations and isotope summaries from
#' `chromPeakData(x)` to grouped feature-level metadata in
#' `featureDefinitions(x)`.
#'
#' For each grouped feature, chromatographic peaks listed in
#' `featureDefinitions(x)$peakidx` are collected. Lipid annotations stored at
#' peak level are expanded, filtered, ranked and collapsed into compact
#' pipe-separated feature-level summaries. Isotope information is also
#' summarized across all chromatographic peaks assigned to each feature.
#'
#' The resulting columns are stored directly in `featureDefinitions(x)`, making
#' the annotated `XcmsExperiment` self-contained and ready for export or
#' conversion to downstream objects such as `SummarizedExperiment`.
#'
#' @param x `XcmsExperiment`.
#'
#' @return
#' Updated `XcmsExperiment` with lipid annotation and isotope summary columns
#' added to `featureDefinitions(x)`.
#'
#' @details
#' Peak-level lipid annotations are expected to be already present in
#' `chromPeakData(x)`, usually after running
#' [addMs2AnnotationsToChromPeakData()].
#'
#' Lipid annotation summaries are generated using:
#'
#' - `.collapse_feature_annotations()`
#' - `.expand_pipe_separated_annotations()`
#' - `.keep_best_annotations_per_sum_composition()`
#' - `.keep_best_feature_annotations_per_class()`
#' - `.sort_ms2_annotations_for_chrom_peak()`
#'
#' Isotope summaries are generated using:
#'
#' - `.summarize_feature_isotopes()`
#'
#' Existing feature-level columns with the same names as the newly generated
#' annotation or isotope summary columns are replaced.
#'
#' Typical added lipid annotation columns include:
#'
#' - `has_monoisotopic_peak`
#' - `n_monoisotopic_peaks`
#' - `monoisotopic_peak_ids`
#' - `monoisotopic_isotope_groups`
#' - `isotope_pattern_summary`
#' - `has_lipid_ms2_annotation`
#' - `n_lipid_ms2_annotations`
#' - `lipid_class`
#' - `sum_composition`
#' - `carbon`
#' - `double_bond`
#' - `oxygen`
#' - `adduct`
#' - `adduct_group_id`
#' - `adduct_group_score`
#' - `adduct_group_n_adducts`
#' - `annotation_level`
#' - `annotation`
#'
#' @seealso
#' [addMs2AnnotationsToChromPeakData()],
#' [exportAnnotatedFeatures()]
#'
#' @examples
#' data(xcmsexperiment_negative)
#'
#' param <- defaultAnnotationParameters("negative")
#'
#' x <- annotateLipids(
#'     xcmsexperiment_negative,
#'     paramannot = param,
#'     paramdia = list(),
#'     verbose = FALSE
#' )
#'
#' x <- addFeatureLipidAnnotations(x)
#'
#' head(xcms::featureDefinitions(x))
#'
#' @export
addFeatureLipidAnnotations <- function(x) {
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("Package 'xcms' is required.", call. = FALSE)
  }

  features <- as.data.frame(xcms::featureDefinitions(x))
  peak_data <- as.data.frame(xcms::chromPeakData(x))
  peaks <- as.data.frame(xcms::chromPeaks(x))

  if (!nrow(features)) {
    stop("'x' does not contain grouped features.", call. = FALSE)
  }

  if (nrow(peaks) != nrow(peak_data)) {
    stop("chromPeaks and chromPeakData are not aligned.", call. = FALSE)
  }

  if (!"peakidx" %in% colnames(features)) {
    stop("featureDefinitions(x) does not contain 'peakidx'.", call. = FALSE)
  }

  .validate_lipid_annotation_columns(peak_data)

  annotation_summary <- lapply(seq_len(nrow(features)), function(i) {
    peak_idx <- features$peakidx[[i]]

    if (is.null(peak_idx) || !length(peak_idx)) {
      return(.empty_feature_annotation_summary())
    }

    peak_idx <- peak_idx[
      peak_idx >= 1L &
        peak_idx <= nrow(peak_data)
    ]

    if (!length(peak_idx)) {
      return(.empty_feature_annotation_summary())
    }

    tab <- peak_data[
      peak_idx,
      .lipid_annotation_export_columns(),
      drop = FALSE
    ]

    .collapse_feature_annotations(tab)
  })

  isotope_summary <- lapply(seq_len(nrow(features)), function(i) {
    peak_idx <- features$peakidx[[i]]

    if (is.null(peak_idx) || !length(peak_idx)) {
      return(.empty_feature_isotope_summary())
    }

    peak_idx <- peak_idx[
      peak_idx >= 1L &
        peak_idx <= nrow(peak_data)
    ]

    if (!length(peak_idx)) {
      return(.empty_feature_isotope_summary())
    }

    tab <- peak_data[peak_idx, , drop = FALSE]

    .summarize_feature_isotopes(tab)
  })

  annotation_summary <- do.call(rbind, annotation_summary)
  isotope_summary <- do.call(rbind, isotope_summary)

  rownames(annotation_summary) <- rownames(features)
  rownames(isotope_summary) <- rownames(features)

  new_data <- cbind(
    isotope_summary,
    annotation_summary
  )

  duplicated_cols <- intersect(
    colnames(features),
    colnames(new_data)
  )

  if (length(duplicated_cols)) {
    features[, duplicated_cols] <- NULL
  }

  features <- cbind(
    features,
    new_data
  )

  xcms::featureDefinitions(x) <- as.data.frame(features)

  x
}

#' Add sample metadata to a peak-level annotation table
#'
#' Append sample metadata from `phenoData(x)` when available.
#'
#' @param out Peak-level table containing either `sample` or `sampleidx`.
#' @param x `XcmsExperiment`.
#'
#' @return `data.frame` with sample metadata appended when possible.
#'
#' @noRd
#'
#' @keywords internal
.add_sample_metadata_to_peak_table <- function(out, x) {
  sample_col <- intersect(c("sample", "sampleidx"), colnames(out))

  if (!length(sample_col)) {
    return(out)
  }

  sample_meta <- .get_xcms_sample_metadata(x)

  if (is.null(sample_meta) || !nrow(sample_meta)) {
    return(out)
  }

  sample_meta$sample <- seq_len(nrow(sample_meta))

  merge(
    out,
    sample_meta,
    by.x = sample_col[[1L]],
    by.y = "sample",
    all.x = TRUE,
    sort = FALSE
  )
}

#' Get sample metadata from an XcmsExperiment
#'
#' Extract sample metadata using `phenoData(x)` when available.
#'
#' @param x `XcmsExperiment`.
#'
#' @return `data.frame` with sample metadata, or `NULL` if unavailable.
#'
#' @noRd
#'
#' @keywords internal
.get_xcms_sample_metadata <- function(x) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    return(NULL)
  }

  if (!"phenoData" %in% methods::slotNames(x)) {
    return(NULL)
  }

  tryCatch(
    as.data.frame(Biobase::pData(Biobase::phenoData(x))),
    error = function(e) NULL
  )
}


#' Return lipid annotation columns exported from chromPeakData
#'
#' @return Character vector with lipid annotation column names.
#'
#' @noRd
#'
#' @keywords internal
.lipid_annotation_export_columns <- function() {
  c(
    "has_lipid_ms2_annotation",
    "n_lipid_ms2_annotations",
    "lipid_class",
    "sum_composition",
    "carbon",
    "double_bond",
    "oxygen",
    "adduct",
    "adduct_group_id",
    "adduct_group_score",
    "adduct_group_n_adducts",
    "annotation_level",
    "annotation"
  )
}


#' Validate lipid annotation columns in chromPeakData
#'
#' @param peak_data `data.frame` returned by `chromPeakData(x)`.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @noRd
#'
#' @keywords internal
.validate_lipid_annotation_columns <- function(peak_data) {
  missing <- setdiff(
    .lipid_annotation_export_columns(),
    colnames(peak_data)
  )

  if (length(missing)) {
    stop(
      "chromPeakData(x) is missing lipid annotation columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Collapse feature-level lipid annotations
#'
#' Collapse lipid annotations from all chromatographic peaks assigned to one
#' grouped feature.
#'
#' Peak-level annotation columns may contain pipe-separated values. This helper
#' expands them back into rows, removes duplicates, sorts annotations by
#' confidence and collapses them again into feature-level pipe-separated values.
#'
#' @param tab `data.frame` with lipid annotation columns from `chromPeakData(x)`.
#'
#' @return One-row `data.frame` with feature-level annotation summary.
#'
#' @noRd
#'
#' @keywords internal
.collapse_feature_annotations <- function(tab) {
  if (!nrow(tab)) {
    return(.empty_feature_annotation_summary())
  }

  tab <- tab[tab$has_lipid_ms2_annotation %in% TRUE, , drop = FALSE]

  if (!nrow(tab)) {
    return(.empty_feature_annotation_summary())
  }

  expanded <- .expand_pipe_separated_annotations(tab)

  if (!nrow(expanded)) {
    return(.empty_feature_annotation_summary())
  }

  expanded <- unique(expanded)

  if (!"adduct_group_n_adducts" %in% colnames(expanded)) {
    expanded$adduct_group_n_adducts <- 1L
  }

  expanded <- .keep_best_annotations_per_sum_composition(expanded)
  expanded <- .keep_best_feature_annotations_per_class(expanded)
  expanded <- .sort_ms2_annotations_for_chrom_peak(expanded)

  data.frame(
    has_lipid_ms2_annotation = TRUE,
    n_lipid_ms2_annotations = nrow(expanded),
    lipid_class = paste(expanded$lipid_class, collapse = "|"),
    sum_composition = paste(expanded$sum_composition, collapse = "|"),
    carbon = paste(expanded$carbon, collapse = "|"),
    double_bond = paste(expanded$double_bond, collapse = "|"),
    oxygen = paste(expanded$oxygen, collapse = "|"),
    adduct = paste(expanded$adduct, collapse = "|"),
    adduct_group_id = paste(expanded$adduct_group_id, collapse = "|"),
    adduct_group_score = paste(round(expanded$adduct_group_score, 3), collapse = "|"),
    adduct_group_n_adducts = paste(expanded$adduct_group_n_adducts, collapse = "|"),
    annotation_level = paste(expanded$annotation_level, collapse = "|"),
    annotation = paste(expanded$annotation, collapse = "|"),
    stringsAsFactors = FALSE
  )
}


#' Expand pipe-separated peak annotations
#'
#' Convert compact peak-level annotation strings into one row per annotation.
#'
#' @param tab `data.frame` with lipid annotation columns from `chromPeakData(x)`.
#'
#' @return `data.frame` with one row per expanded lipid annotation.
#'
#' @noRd
#'
#' @keywords internal
.expand_pipe_separated_annotations <- function(tab) {
  cols <- c(
    "lipid_class",
    "sum_composition",
    "carbon",
    "double_bond",
    "oxygen",
    "adduct",
    "adduct_group_id",
    "adduct_group_score",
    "adduct_group_n_adducts",
    "annotation_level",
    "annotation"
  )

  out <- vector("list", nrow(tab))
  k <- 1L

  for (i in seq_len(nrow(tab))) {
    values <- lapply(cols, function(col) {
      value <- as.character(tab[[col]][[i]])

      if (is.na(value) || !nzchar(value)) {
        return(NA_character_)
      }

      strsplit(value, "\\|", perl = TRUE)[[1L]]
    })

    names(values) <- cols
    n <- max(vapply(values, length, integer(1)))

    values <- lapply(values, function(z) rep_len(z, n))

    out[[k]] <- as.data.frame(values, stringsAsFactors = FALSE)
    k <- k + 1L
  }

  out <- do.call(rbind, out[seq_len(k - 1L)])
  out <- out[!is.na(out$annotation) & nzchar(out$annotation), , drop = FALSE]

  out$adduct_group_score <- as.numeric(out$adduct_group_score)

  out$adduct_group_n_adducts <- as.integer(out$adduct_group_n_adducts)

  rownames(out) <- NULL
  out
}


#' Return an empty feature-level annotation summary
#'
#' @return One-row `data.frame` with empty feature annotation fields.
#'
#' @noRd
#'
#' @keywords internal
.empty_feature_annotation_summary <- function() {
  data.frame(
    has_lipid_ms2_annotation = FALSE,
    n_lipid_ms2_annotations = 0L,
    lipid_class = NA_character_,
    sum_composition = NA_character_,
    carbon = NA_character_,
    double_bond = NA_character_,
    oxygen = NA_character_,
    adduct = NA_character_,
    adduct_group_id = NA_character_,
    adduct_group_score = NA_character_,
    adduct_group_n_adducts = NA_character_,
    annotation_level = NA_character_,
    annotation = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Write an annotation table
#'
#' @param x `data.frame`.
#' @param file `character(1)` output path.
#' @param sep Field separator.
#'
#' @return Invisibly returns `file`.
#'
#' @noRd
#'
#' @keywords internal
.write_annotation_table <- function(x, file, sep = ",") {

  if (!is.character(file) || length(file) != 1L ||
      is.na(file) || !nzchar(file)) {
    stop("'file' must be a non-empty character(1).", call. = FALSE)
  }

  ext <- tolower(tools::file_ext(file))

  if (ext == "xlsx") {

    if (!requireNamespace("writexl", quietly = TRUE)) {
      stop(
        "Package 'writexl' is required to write .xlsx files.",
        call. = FALSE
      )
    }

    writexl::write_xlsx(x, path = file)

  } else {

    utils::write.table(
      x,
      file = file,
      sep = sep,
      quote = TRUE,
      row.names = FALSE,
      col.names = TRUE
    )

  }

  invisible(file)
}

#' Summarize isotope confirmation for one grouped feature
#'
#' Summarize isotope-related information across all chromatographic peaks
#' belonging to one grouped feature, focusing on confirmation of monoisotopic
#' (`M+0`) peaks.
#'
#' The function uses `is_isotope`, `isotope_group`, and `isotope_rank` from
#' `chromPeakData(x)`. Peaks with `isotope_rank == 0` are considered
#' monoisotopic peaks (`M+0`), regardless of whether they are flagged as isotope
#' peaks or not.
#'
#' The main exported fields indicate whether the feature contains at least one
#' confirmed `M+0` peak and which chromatographic peaks and isotope groups
#' support that confirmation.
#'
#' A compact `isotope_pattern_summary` field is also returned to keep the result
#' flexible for later inspection of the full isotope pattern. This field stores
#' one entry per isotope group in the format:
#'
#' `group=>rank1,rank2,...`
#'
#' For example:
#'
#' `sample1:grp12=>0,1,2|sample2:grp8=>0,1`
#'
#' @param tab `data.frame` containing `chromPeakData(x)` rows for one feature.
#'
#' @return One-row `data.frame` with feature-level isotope summary.
#'
#' @noRd
#'
#' @keywords internal
.summarize_feature_isotopes <- function(tab) {
  if (!nrow(tab) || !"isotope_rank" %in% colnames(tab)) {
    return(.empty_feature_isotope_summary())
  }

  rank <- tab$isotope_rank

  if (is.factor(rank)) {
    rank <- as.character(rank)
  }

  if (is.character(rank)) {
    valid <- is.na(rank) | grepl("^[0-9]+$", rank)

    if (!all(valid)) {
      stop(
        "'isotope_rank' contains non-integer values.",
        call. = FALSE
      )
    }
  }

  rank <- as.integer(rank)

  is_m0 <- !is.na(rank) & rank == 0L
  m0_rows <- tab[is_m0, , drop = FALSE]

  out <- data.frame(
    has_monoisotopic_peak = any(is_m0),
    n_monoisotopic_peaks = sum(is_m0),
    monoisotopic_peak_ids = NA_character_,
    monoisotopic_isotope_groups = NA_character_,
    isotope_pattern_summary = NA_character_,
    stringsAsFactors = FALSE
  )

  if (nrow(m0_rows)) {
    if ("chrom_peak_id" %in% colnames(m0_rows)) {
      out$monoisotopic_peak_ids <- paste(
        unique(m0_rows$chrom_peak_id),
        collapse = "|"
      )
    }

    if ("isotope_group" %in% colnames(m0_rows)) {
      vals <- unique(as.character(m0_rows$isotope_group))
      vals <- vals[!is.na(vals) & nzchar(vals)]

      if (length(vals)) {
        out$monoisotopic_isotope_groups <- paste(vals, collapse = "|")
      }
    }
  }

  out$isotope_pattern_summary <- .summarize_feature_isotope_pattern(tab)

  out
}

#' Return an empty feature-level isotope summary
#'
#' @return One-row `data.frame` with empty isotope summary fields.
#'
#' @noRd
#'
#' @keywords internal
.empty_feature_isotope_summary <- function() {
  data.frame(
    has_monoisotopic_peak = FALSE,
    n_monoisotopic_peaks = 0L,
    monoisotopic_peak_ids = NA_character_,
    monoisotopic_isotope_groups = NA_character_,
    isotope_pattern_summary = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Summarize isotope pattern by isotope group
#'
#' Build a compact feature-level summary of isotope ranks observed in each
#' isotope group.
#'
#' The output contains one entry per isotope group in the format:
#'
#' `group=>rank1,rank2,...`
#'
#' Entries are collapsed with `"|"`.
#'
#' @param tab `data.frame` containing `chromPeakData(x)` rows for one feature.
#'
#' @return `character(1)` with the summarized isotope pattern, or
#' `NA_character_` if unavailable.
#'
#' @noRd
#'
#' @keywords internal
.summarize_feature_isotope_pattern <- function(tab) {
  required <- c("isotope_group", "isotope_rank")

  if (!all(required %in% colnames(tab)) || !nrow(tab)) {
    return(NA_character_)
  }

  group <- as.character(tab$isotope_group)
  rank <- tab$isotope_rank

  if (is.factor(rank)) {
    rank <- as.character(rank)
  }

  if (is.character(rank)) {
    valid <- is.na(rank) | grepl("^[0-9]+$", rank)

    if (!all(valid)) {
      stop(
        "'isotope_rank' contains non-integer values.",
        call. = FALSE
      )
    }
  }

  rank <- as.integer(rank)

  keep <- !is.na(group) & nzchar(group) & !is.na(rank)

  if (!any(keep)) {
    return(NA_character_)
  }

  tab2 <- data.frame(
    isotope_group = group[keep],
    isotope_rank = rank[keep],
    stringsAsFactors = FALSE
  )

  groups <- split(tab2$isotope_rank, tab2$isotope_group)

  out <- vapply(names(groups), function(grp) {
    vals <- sort(unique(groups[[grp]]))
    paste0(grp, "=>", paste(vals, collapse = ","))
  }, character(1), USE.NAMES = FALSE)

  paste(out, collapse = "|")
}

#' Keep best feature-level annotations per lipid class
#'
#' Collapse feature-level annotations across samples when several annotations
#' share the same lipid class but differ in total composition.
#'
#' For each `lipid_class`, only annotations with the highest
#' `annotation_level` are kept. If several annotations share the same best
#' level, all unique annotations are retained.
#'
#' This helper is intended for feature-level summaries, where annotations from
#' multiple chromatographic peaks and samples are merged.
#'
#' @param x `data.frame` with annotation rows. Must contain columns
#' `lipid_class`, `annotation_level`, and `annotation`.
#'
#' @return Filtered `data.frame`.
#'
#' @noRd
#'
#' @keywords internal
.keep_best_feature_annotations_per_class <- function(x) {
  if (!nrow(x)) {
    return(x)
  }

  required <- c("lipid_class", "adduct_group_score", "adduct_group_n_adducts", "annotation_level", "annotation")

  missing <- setdiff(required, colnames(x))
  if (length(missing)) {
    stop(
      "'x' is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  level_rank <- c(
    rejected = 0L,
    `MS-only` = 1L,
    Subclass = 2L,
    FA = 3L,
    `FA position` = 4L
  )

  x$annotation_level_rank <- unname(level_rank[x$annotation_level])
  x$annotation_level_rank[is.na(x$annotation_level_rank)] <- 0L

  if (!"adduct_group_score" %in% colnames(x)) {
    x$adduct_group_score <- 0
  }

  if (!"adduct_group_n_adducts" %in% colnames(x)) {
    x$adduct_group_n_adducts <- 1L
  }

  x$adduct_group_score <- as.numeric(x$adduct_group_score)
  x$adduct_group_score[is.na(x$adduct_group_score)] <- 0

  x$adduct_group_n_adducts <- as.integer(x$adduct_group_n_adducts)
  x$adduct_group_n_adducts[is.na(x$adduct_group_n_adducts)] <- 1L

  groups <- split(seq_len(nrow(x)), x$lipid_class)

  out <- lapply(groups, function(idx) {
    tab <- x[idx, , drop = FALSE]

    best_score <- max(tab$adduct_group_score, na.rm = TRUE)
    tab <- tab[tab$adduct_group_score == best_score, , drop = FALSE]

    best_rank <- max(tab$annotation_level_rank, na.rm = TRUE)
    tab <- tab[tab$annotation_level_rank == best_rank, , drop = FALSE]

    best_n <- max(tab$adduct_group_n_adducts, na.rm = TRUE)
    tab <- tab[tab$adduct_group_n_adducts == best_n, , drop = FALSE]

    tab <- tab[!duplicated(tab$annotation), , drop = FALSE]
    tab
  })

  out <- do.call(rbind, out)

  out$annotation_level_rank <- NULL
  rownames(out) <- NULL
  out
}

#' Add accepted MS2 fragment annotations to extracted MS2 data
#'
#' Adds lipid annotation evidence to the fragment tables stored in an
#' extracted MS2 data object. Only evidence supporting lipid annotations that
#' were accepted and written to `chromPeakData(x)` is propagated.
#'
#' Fragment-level annotations include the evidence scope, fragmentation rule
#' identifier and, when available, chain-specific information.
#'
#' @param ms2_data `list`. Extracted MS2 data object containing a `samples`
#'   element with spectra grouped by chromatographic peak identifier.
#' @param ms2_annotation `list`. MS2 annotation results containing an
#'   `evidence` element.
#' @param x `XCMSnExp` or `XcmsExperiment` object. Object containing lipid
#'   annotations in `chromPeakData(x)`.
#' @param ppm `numeric(1)`. Mass tolerance, in ppm, used to match evidence
#'   fragments to extracted fragment peaks.
#'
#' @return A modified `ms2_data` object with fragment-level annotation columns
#'   added to each spectrum table:
#'   `has_lipid_fragment_annotation`, `evidence_scope`, `rule_id`,
#'   `chain_role`, `chain_composition`, `chain_carbon`,
#'   `chain_double_bond` and `chain_oxygen`.
#'
#' @noRd
#'
#' @keywords internal
annotateAcceptedFragmentsInMs2Data <- function(ms2_data,
                                               ms2_annotation,
                                               x,
                                               ppm = 10) {
  if (is.null(ms2_data) || !is.list(ms2_data) ||
      !"samples" %in% names(ms2_data)) {
    stop("'ms2_data' must contain element 'samples'.", call. = FALSE)
  }

  if (is.null(ms2_annotation$evidence) ||
      !nrow(ms2_annotation$evidence)) {
    return(ms2_data)
  }

  accepted_keys <- .get_accepted_annotation_keys_from_chrom_peak_data(x)

  if (!nrow(accepted_keys)) {
    return(ms2_data)
  }

  evidence <- ms2_annotation$evidence
  evidence <- evidence[evidence$supports_annotation, , drop = FALSE]

  evidence <- merge(
    evidence,
    accepted_keys,
    by = c("adduct_group_id", "lipid_class", "sum_composition"),
    all = FALSE,
    sort = FALSE
  )

  if (!nrow(evidence)) {
    return(ms2_data)
  }

  for (sample_name in names(ms2_data$samples)) {
    sample_i <- ms2_data$samples[[sample_name]]

    if (is.null(sample_i$spectra) || !length(sample_i$spectra)) {
      next
    }

    for (peak_id in names(sample_i$spectra)) {
      fragments <- sample_i$spectra[[peak_id]]

      if (is.null(fragments) || !nrow(fragments)) {
        next
      }

      sample_i$spectra[[peak_id]] <- .annotate_fragment_table_with_evidence(
        fragments = fragments,
        evidence = evidence,
        peak_id = peak_id,
        ppm = ppm
      )
    }

    ms2_data$samples[[sample_name]] <- sample_i
  }

  ms2_data
}

#' Get accepted lipid annotation keys from chromPeakData
#'
#' Extracts the accepted lipid annotation keys written to `chromPeakData(x)`.
#' Multi-annotation fields separated by `"|"` are expanded so that each
#' accepted lipid annotation is represented by one row.
#'
#' @param x `XCMSnExp` or `XcmsExperiment` object containing lipid annotation
#'   columns in `chromPeakData(x)`.
#'
#' @return A `data.frame` with unique accepted annotation keys:
#'   `adduct_group_id`, `lipid_class` and `sum_composition`.
#'
#' @noRd
#'
#' @keywords internal
.get_accepted_annotation_keys_from_chrom_peak_data <- function(x) {
  cpd <- as.data.frame(xcms::chromPeakData(x))

  required <- c(
    "has_lipid_ms2_annotation",
    "adduct_group_id",
    "lipid_class",
    "sum_composition"
  )

  missing <- setdiff(required, colnames(cpd))

  if (length(missing)) {
    stop(
      "'chromPeakData(x)' is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  cpd <- cpd[cpd$has_lipid_ms2_annotation, , drop = FALSE]

  if (!nrow(cpd)) {
    return(data.frame(
      adduct_group_id = character(),
      lipid_class = character(),
      sum_composition = character(),
      stringsAsFactors = FALSE
    ))
  }

  out <- lapply(seq_len(nrow(cpd)), function(i) {
    data.frame(
      adduct_group_id = strsplit(cpd$adduct_group_id[[i]], "|", fixed = TRUE)[[1L]],
      lipid_class = strsplit(cpd$lipid_class[[i]], "|", fixed = TRUE)[[1L]],
      sum_composition = strsplit(cpd$sum_composition[[i]], "|", fixed = TRUE)[[1L]],
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  out <- unique(out)
  rownames(out) <- NULL
  out
}

#' Annotate one fragment table with lipid evidence
#'
#' Matches observed fragment peaks in a single MS2 fragment table against
#' lipid annotation evidence and adds fragment-level annotation columns.
#'
#' When several evidence rows match the same fragment peak, values are collapsed
#' using `"|"`.
#'
#' @param fragments `data.frame`. Fragment table for one chromatographic peak.
#'   Must contain an `mz` column.
#' @param evidence `data.frame`. Lipid annotation evidence table.
#' @param peak_id `character(1)`. Chromatographic peak identifier used to
#'   subset evidence by `ms2_source_peak_id`.
#' @param ppm `numeric(1)`. Mass tolerance, in ppm.
#'
#' @return The input `fragments` table with additional fragment annotation
#'   columns.
#'
#' @noRd
#'
#' @keywords internal
.annotate_fragment_table_with_evidence <- function(fragments,
                                                   evidence,
                                                   peak_id,
                                                   ppm) {
  out_cols <- c(
    "evidence_scope",
    "rule_id",
    "chain_role",
    "chain_composition",
    "chain_carbon",
    "chain_double_bond",
    "chain_oxygen"
  )

  for (col in out_cols) {
    fragments[[col]] <- NA_character_
  }

  fragments$has_lipid_fragment_annotation <- FALSE

  evidence <- evidence[
    evidence$ms2_source_peak_id == peak_id,
    ,
    drop = FALSE
  ]

  if (!nrow(evidence)) {
    return(fragments)
  }

  for (i in seq_len(nrow(fragments))) {

    mz_i <- fragments$mz[[i]]

    if (is.na(mz_i)) {
      next
    }

    tol <- mz_i * ppm * 1e-6

    hit <- abs(evidence$observed_mz - mz_i) <= tol
    hit[is.na(hit)] <- FALSE

    if ("spectrum_id" %in% colnames(fragments)) {

      sid <- as.character(fragments$spectrum_id[[i]])

      sid_hit <- (
        is.na(evidence$ms2_source_spectrum_id) |
          evidence$ms2_source_spectrum_id == sid
      )

      sid_hit[is.na(sid_hit)] <- FALSE

      hit <- hit & sid_hit
    }

    if (!any(hit, na.rm = TRUE)) {
      next
    }

    ev_i <- evidence[hit, , drop = FALSE]

    fragments$has_lipid_fragment_annotation[[i]] <- TRUE

    for (col in out_cols) {

      value <- unique(as.character(ev_i[[col]]))

      value <- value[!is.na(value) & nzchar(value)]

      if (length(value)) {
        fragments[[col]][[i]] <- paste(value, collapse = "|")
      }
    }
  }

  fragments
}

.merge_chrom_peak_data_columns_by_id <- function(x_full, x_subset, verbose = TRUE) {
  cpd_full <- as.data.frame(xcms::chromPeakData(x_full))
  cpd_subset <- as.data.frame(xcms::chromPeakData(x_subset))

  ids <- intersect(rownames(cpd_subset), rownames(cpd_full))

  if (!length(ids)) {
    stop("Could not merge chromPeakData: no shared chromPeak IDs between full and subset objects.",
         call. = FALSE)
  }

  for (cn in colnames(cpd_subset)) {
    if (!cn %in% colnames(cpd_full)) {
      proto <- cpd_subset[[cn]]

      if (is.logical(proto)) {
        # Annotation flags should be FALSE for non-annotated/fullscan peaks.
        cpd_full[[cn]] <- rep(FALSE, nrow(cpd_full))
      } else if (is.integer(proto)) {
        cpd_full[[cn]] <- rep(NA_integer_, nrow(cpd_full))
      } else if (is.numeric(proto)) {
        cpd_full[[cn]] <- rep(NA_real_, nrow(cpd_full))
      } else if (is.factor(proto)) {
        cpd_full[[cn]] <- factor(rep(NA_character_, nrow(cpd_full)),
                                 levels = levels(proto))
      } else if (is.character(proto)) {
        cpd_full[[cn]] <- rep(NA_character_, nrow(cpd_full))
      } else if (is.list(proto)) {
        cpd_full[[cn]] <- I(vector("list", nrow(cpd_full)))
      } else {
        cpd_full[[cn]] <- rep(NA, nrow(cpd_full))
      }
    }

    cpd_full[ids, cn] <- cpd_subset[ids, cn]
  }

  xcms::chromPeakData(x_full) <- as.data.frame(cpd_full)

  if (isTRUE(verbose)) {
    cat("Merged chromPeakData columns from annotation subset into full object. Peaks merged: ",
        length(ids), "\n", sep = "")
  }

  x_full
}

#' Add RefMet information to LipidMS annotations
#'
#' Add RefMet standardized names, molecular formulas and external database
#' identifiers to lipid annotations stored in chromPeakData() and
#' featureDefinitions().
#'
#' @param x An annotated XcmsExperiment.
#'
#' @return Updated XcmsExperiment.
#'
#' @keywords internal
.add_refmet_annotations <- function(x) {

  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("Package 'xcms' is required.", call. = FALSE)
  }

  cpd <- as.data.frame(xcms::chromPeakData(x))

  if (!"annotation" %in% colnames(cpd)) {
    return(x)
  }

  fd <- tryCatch(
    as.data.frame(xcms::featureDefinitions(x)),
    error = function(e) NULL
  )

  annotations <- .extract_refmet_annotations(cpd$annotation)

  if (
    !is.null(fd) &&
    nrow(fd) &&
    "annotation" %in% colnames(fd)
  ) {
    annotations <- unique(c(
      annotations,
      .extract_refmet_annotations(fd$annotation)
    ))
  }

  ## No annotations present
  if (!length(annotations)) {
    cpd <- .initialize_refmet_columns(cpd)
    xcms::chromPeakData(x) <- cpd

    if (!is.null(fd) && nrow(fd)) {
      fd <- .initialize_refmet_columns(fd)
      xcms::featureDefinitions(x) <- fd
    }

    return(x)
  }

  refmet <- RefMet::refmet_map_df(annotations)

  if (nrow(refmet) != length(annotations)) {
    stop(
      "refmet_map_df() did not return one row per input annotation.",
      call. = FALSE
    )
  }

  refmet_map <- data.frame(
    annotation = annotations,
    formula = as.character(refmet$Formula),
    refmet_name = as.character(refmet$Standardized.name),
    refmet_id = as.character(refmet$RefMet_ID),
    lipidmaps_id = as.character(refmet$LM_ID),
    hmdb_id = as.character(refmet$HMDB_ID),
    kegg_id = as.character(refmet$KEGG_ID),
    stringsAsFactors = FALSE
  )

  cpd <- .add_refmet_to_table(
    tab = cpd,
    refmet_map = refmet_map
  )

  xcms::chromPeakData(x) <- cpd

  if (
    !is.null(fd) &&
    nrow(fd) &&
    "annotation" %in% colnames(fd)
  ) {
    fd <- .add_refmet_to_table(
      tab = fd,
      refmet_map = refmet_map
    )

    xcms::featureDefinitions(x) <- fd
  }

  x
}

#' Extract unique annotations for RefMet mapping
#'
#' @noRd
#'
#' @keywords internal
.extract_refmet_annotations <- function(x) {

  x <- as.character(x)

  x <- x[
    !is.na(x) &
      nzchar(trimws(x))
  ]

  if (!length(x)) {
    return(character())
  }

  annotations <- unlist(
    strsplit(x, "\\|", perl = TRUE),
    use.names = FALSE
  )

  annotations <- trimws(annotations)

  annotations <- annotations[
    !is.na(annotations) &
      nzchar(annotations)
  ]

  unique(annotations)
}

#' Initialize RefMet annotation columns
#'
#' @noRd
#'
#' @keywords internal
.initialize_refmet_columns <- function(tab) {

  cols <- c(
    "formula",
    "refmet_name",
    "refmet_id",
    "lipidmaps_id",
    "hmdb_id",
    "kegg_id"
  )

  for (col in cols) {
    tab[[col]] <- NA_character_
  }

  tab
}

#' Add RefMet fields to an annotation table
#'
#' @noRd
#'
#' @keywords internal
.add_refmet_to_table <- function(tab, refmet_map) {

  tab <- .initialize_refmet_columns(tab)

  if (
    !nrow(tab) ||
    !"annotation" %in% colnames(tab) ||
    !nrow(refmet_map)
  ) {
    return(tab)
  }

  valid <- !is.na(tab$annotation) &
    nzchar(trimws(as.character(tab$annotation)))

  if (!any(valid)) {
    return(tab)
  }

  refmet_cols <- c(
    "formula",
    "refmet_name",
    "refmet_id",
    "lipidmaps_id",
    "hmdb_id",
    "kegg_id"
  )

  for (i in which(valid)) {

    annotations <- strsplit(
      as.character(tab$annotation[[i]]),
      "\\|",
      perl = TRUE
    )[[1L]]

    annotations <- trimws(annotations)

    idx <- match(
      annotations,
      refmet_map$annotation
    )

    for (col in refmet_cols) {

      values <- refmet_map[[col]][idx]

      tab[[col]][[i]] <- .collapse_refmet_values(values)
    }
  }

  tab
}

#' Collapse RefMet values preserving annotation positions
#'
#' @noRd
#'
#' @keywords internal
.collapse_refmet_values <- function(x) {

  x <- as.character(x)

  if (!length(x) || all(is.na(x) | !nzchar(x))) {
    return(NA_character_)
  }

  x[is.na(x)] <- ""

  paste(x, collapse = "|")
}

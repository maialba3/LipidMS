# ==============================================================================
# MS2 fragment-rule matching
# ==============================================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Match MS2 fragment rules to lipid candidates
#'
#' Match class-level and chain-level MS2 fragmentation rules against all MS2
#' fragments available for each `adduct_group_id`.
#'
#' @param candidates `data.frame` returned by `findMs1LipidCandidates()`.
#' @param ms2_data Output from `extractMs2Spectra()`.
#' @param rules Named list with elements `class_fragments`,
#' `chain_fragments`, optionally `chain_inference`, and optionally
#' `relative_intensity`.
#' @param dbs Object returned by `createLipidDbs()`.
#' @param paramannot Annotation parameter list. Requires `polarity` and
#' `dmzproducts`.
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
#' annotations <- matchFragmentRules(
#'     candidates,
#'     ms2_data,
#'     rules,
#'     dbs,
#'     param
#' )
#'
#' annotations$annotations
#'
#' @return Named list with:
#' - `annotations`: best annotation per `adduct_group_id`
#' - `evidence`: matched MS2 fragment evidence
#' - `all_annotations`: candidate-level annotation attempts
#'
#' @export
matchFragmentRules <- function(candidates,
                               ms2_data,
                               rules,
                               dbs,
                               paramannot) {
  .validate_fragment_matching_inputs(candidates, paramannot)

  if (!nrow(candidates)) {
    return(list(
      annotations = .empty_ms2_annotation_table(),
      evidence = .empty_ms2_evidence_table(),
      all_annotations = .empty_ms2_annotation_table()
    ))
  }

  candidates <- candidates[
    !is.na(candidates$adduct_group_id) &
      nzchar(candidates$adduct_group_id),
    ,
    drop = FALSE
  ]

  if (!nrow(candidates)) {
    return(list(
      annotations = .empty_ms2_annotation_table(),
      evidence = .empty_ms2_evidence_table(),
      all_annotations = .empty_ms2_annotation_table()
    ))
  }

  if (is.null(rules$class_fragments)) {
    rules$class_fragments <- .empty_class_rule_table()
  }

  if (is.null(rules$chain_fragments)) {
    rules$chain_fragments <- .empty_chain_rule_table()
  }

  if (is.null(rules$chain_inference_rules)) {
    rules$chain_inference_rules <- defaultChainInferenceRules()
  }

  if (is.null(rules$combination_rules)) {
    rules$combination_rules <- defaultChainCombinationRules()
  }

  if (is.null(rules$relative_intensity)) {
    rules$relative_intensity <- .empty_intensity_rule_table()
  }

  groups <- split(seq_len(nrow(candidates)), candidates$adduct_group_id)

  ann_out <- vector("list", length(groups))
  ev_out <- vector("list", length(groups))
  k <- 1L

  for (group_id in names(groups)) {
    group_candidates <- candidates[groups[[group_id]], , drop = FALSE]

    group_spectra <- collectAdductGroupSpectra(
      group_candidates = group_candidates,
      ms2_data = ms2_data
    )

    res <- .match_fragment_rules_one_group(
      group_candidates = group_candidates,
      group_spectra = group_spectra,
      rules = rules,
      dbs = dbs,
      paramannot = paramannot
    )

    ann_out[[k]] <- res$annotations
    ev_out[[k]] <- res$evidence
    k <- k + 1L
  }

  all_annotations <- .rbind_or_empty(
    ann_out,
    .empty_ms2_annotation_table()
  )

  evidence <- .rbind_or_empty(
    ev_out,
    .empty_ms2_evidence_table()
  )

  annotations <- .keep_best_group_annotation(all_annotations)

  list(
    annotations = annotations,
    evidence = evidence,
    all_annotations = all_annotations
  )
}

#' Match MS2 rules for one adduct group
#'
#' @param group_candidates Candidate rows from one adduct group.
#' @param group_spectra MS2 fragment table for the adduct group.
#' @param rules Fragment rule list.
#' @param dbs Lipid database object.
#' @param paramannot Annotation parameter list.
#'
#' @return Named list with `annotations` and `evidence`.
#'
#' @noRd
#'
#' @keywords internal
.match_fragment_rules_one_group <- function(group_candidates,
                                            group_spectra,
                                            rules,
                                            dbs,
                                            paramannot) {
  ann_out <- vector("list", nrow(group_candidates))
  ev_out <- vector("list", nrow(group_candidates))

  for (i in seq_len(nrow(group_candidates))) {
    candidate <- group_candidates[i, , drop = FALSE]

    class_evidence <- matchClassRules(
      candidate = candidate,
      group_spectra = group_spectra,
      class_rules = rules$class_fragments,
      dbs = dbs,
      paramannot = paramannot
    )

    if (.has_exclude_evidence(class_evidence)) {
      ann_out[[i]] <- .build_ms2_annotation_row(
        candidate = candidate,
        annotation_level = "rejected",
        annotation = NA_character_,
        score_ms2 = 0,
        status = "excluded_by_class_rule"
      )
      ev_out[[i]] <- class_evidence
      next
    }

    chain_evidence <- matchChainRules(
      candidate = candidate,
      group_spectra = group_spectra,
      chain_rules = rules$chain_fragments,
      dbs = dbs,
      paramannot = paramannot
    )

    chain_sets <- buildChainCombinations(
      candidate = candidate,
      chain_evidence = chain_evidence,
      inference_rules = rules$chain_inference_rules,
      combination_rules = rules$combination_rules,
      paramannot = paramannot,
      dbs = dbs
    )

    if (.needs_relative_intensity_matching(chain_sets)) {
      intensity_evidence <- matchRelativeIntensityRules(
        candidate = candidate,
        chain_sets = chain_sets,
        chain_evidence = chain_evidence,
        intensity_rules = rules$relative_intensity,
        paramannot = paramannot
      )
    } else {
      intensity_evidence <- .empty_ms2_evidence_table()
    }

    ann_out[[i]] <- .rank_candidate_annotation(
      candidate = candidate,
      class_evidence = class_evidence,
      chain_evidence = chain_evidence,
      chain_sets = chain_sets,
      intensity_evidence = intensity_evidence
    )

    ev_out[[i]] <- .rbind_or_empty(
      list(class_evidence, chain_evidence, intensity_evidence),
      .empty_ms2_evidence_table()
    )
  }

  list(
    annotations = .rbind_or_empty(ann_out, .empty_ms2_annotation_table()),
    evidence = .rbind_or_empty(ev_out, .empty_ms2_evidence_table())
  )
}

#' Collect MS2 spectra for one adduct group
#'
#' @param group_candidates Candidate rows from one adduct group.
#' @param ms2_data Output from `extractMs2Spectra()`.
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
#' group_ids <- unique(stats::na.omit(candidates$adduct_group_id))
#'
#' if (length(group_ids)) {
#'     group_candidates <- candidates[
#'         candidates$adduct_group_id == group_ids[1],
#'         ,
#'         drop = FALSE
#'     ]
#'
#'     spectra <- collectAdductGroupSpectra(
#'         group_candidates,
#'         ms2_data
#'     )
#'
#'     head(spectra)
#' }
#'
#' @return `data.frame` with MS2 fragments from all candidate peaks in the
#' adduct group.
#'
#' @export
collectAdductGroupSpectra <- function(group_candidates, ms2_data) {
  peak_ids <- unique(group_candidates$chrom_peak_id)

  out <- lapply(peak_ids, function(pid) {
    .extract_ms2_fragments_for_peak(ms2_data = ms2_data, peak_id = pid)
  })

  .rbind_or_empty(out, .empty_group_spectra_table())
}

#' Extract MS2 fragments for one chromatographic peak
#'
#' @param ms2_data Output from `extractMs2Spectra()`.
#' @param peak_id MS1 chromatographic peak ID.
#'
#' @return `data.frame` with columns `mz`, `intensity`, `chrom_peak_id`,
#' `ms1_chrom_peak_id`, `ms2_chrom_peak_id`, and `spectrum_id`.
#'
#' @noRd
#'
#' @keywords internal
.extract_ms2_fragments_for_peak <- function(ms2_data, peak_id) {
  if (is.null(ms2_data) || !is.list(ms2_data) ||
      !"samples" %in% names(ms2_data)) {
    stop(
      "'ms2_data' must be the output from extractMs2Spectra() with element ",
      "'samples'.",
      call. = FALSE
    )
  }

  out <- vector("list", length(ms2_data$samples))
  k <- 1L

  for (sample_i in ms2_data$samples) {
    if (is.null(sample_i) || !"spectra" %in% names(sample_i)) {
      next
    }

    spectra_i <- sample_i$spectra

    if (!is.list(spectra_i) || !peak_id %in% names(spectra_i)) {
      next
    }

    frag_i <- spectra_i[[peak_id]]

    if (is.null(frag_i) || !is.data.frame(frag_i) || !nrow(frag_i)) {
      next
    }

    mz_col <- intersect(c("mz", "fragment_mz"), colnames(frag_i))
    int_col <- intersect(c("intensity", "int", "fragment_intensity"),
                         colnames(frag_i))

    if (!length(mz_col) || !length(int_col)) {
      next
    }

    tmp <- data.frame(
      mz = as.numeric(frag_i[[mz_col[[1L]]]]),
      intensity = as.numeric(frag_i[[int_col[[1L]]]]),
      chrom_peak_id = if ("chrom_peak_id" %in% colnames(frag_i)) {
        as.character(frag_i$chrom_peak_id)
      } else {
        rep(peak_id, nrow(frag_i))
      },
      ms1_chrom_peak_id = if ("ms1_chrom_peak_id" %in% colnames(frag_i)) {
        as.character(frag_i$ms1_chrom_peak_id)
      } else {
        rep(peak_id, nrow(frag_i))
      },
      ms2_chrom_peak_id = if ("ms2_chrom_peak_id" %in% colnames(frag_i)) {
        as.character(frag_i$ms2_chrom_peak_id)
      } else {
        rep(NA_character_, nrow(frag_i))
      },
      spectrum_id = if ("spectrum_id" %in% colnames(frag_i)) {
        as.character(frag_i$spectrum_id)
      } else if ("ms2_spectrum_id" %in% colnames(frag_i)) {
        as.character(frag_i$ms2_spectrum_id)
      } else {
        rep(NA_character_, nrow(frag_i))
      },
      stringsAsFactors = FALSE
    )

    tmp <- tmp[is.finite(tmp$mz) & is.finite(tmp$intensity), , drop = FALSE]

    if (!nrow(tmp)) {
      next
    }

    out[[k]] <- tmp
    k <- k + 1L
  }

  if (k == 1L) {
    return(.empty_group_spectra_table())
  }

  out <- do.call(rbind, out[seq_len(k - 1L)])
  rownames(out) <- NULL
  out
}

#' Match class-level rules
#'
#' @param candidate One-row candidate table.
#' @param group_spectra MS2 fragment table for the adduct group.
#' @param class_rules Class-level fragmentation rules.
#' @param dbs Lipid database object.
#' @param paramannot Annotation parameter list.
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
#' group_ids <- unique(stats::na.omit(candidates$adduct_group_id))
#'
#' if (length(group_ids)) {
#'     group_candidates <- candidates[
#'         candidates$adduct_group_id == group_ids[1],
#'         ,
#'         drop = FALSE
#'     ]
#'
#'     group_spectra <- collectAdductGroupSpectra(
#'         group_candidates,
#'         ms2_data
#'     )
#'
#'     evidence <- matchClassRules(
#'         candidate = group_candidates[1, , drop = FALSE],
#'         group_spectra = group_spectra,
#'         class_rules = rules$class_fragments,
#'         dbs = dbs,
#'         paramannot = param
#'     )
#'
#'     evidence
#' }
#'
#' @return MS2 evidence table.
#'
#' @export
matchClassRules <- function(candidate,
                            group_spectra,
                            class_rules,
                            dbs,
                            paramannot) {
  if (!nrow(group_spectra) || is.null(class_rules) || !nrow(class_rules)) {
    return(.empty_ms2_evidence_table())
  }

  rules_i <- class_rules[
    class_rules$lipid_class == candidate$class[[1L]] &
      class_rules$polarity == paramannot$polarity,
    ,
    drop = FALSE
  ]

  if (!nrow(rules_i)) {
    return(.empty_ms2_evidence_table())
  }

  expected <- .build_expected_class_fragments(
    candidate = candidate,
    rules = rules_i,
    dbs = dbs
  )

  .match_expected_fragments(
    candidate = candidate,
    expected = expected,
    observed = group_spectra,
    ppm = paramannot$dmzproducts,
    evidence_scope = "class"
  )
}

#' Match chain-level rules
#'
#' @param candidate One-row candidate table.
#' @param group_spectra MS2 fragment table for the adduct group.
#' @param chain_rules Chain-level fragmentation rules.
#' @param dbs Lipid database object.
#' @param paramannot Annotation parameter list.
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
#' group_ids <- unique(stats::na.omit(candidates$adduct_group_id))
#'
#' if (length(group_ids)) {
#'     group_candidates <- candidates[
#'         candidates$adduct_group_id == group_ids[1],
#'         ,
#'         drop = FALSE
#'     ]
#'
#'     group_spectra <- collectAdductGroupSpectra(
#'         group_candidates,
#'         ms2_data
#'     )
#'
#'     evidence <- matchChainRules(
#'         candidate = group_candidates[1, , drop = FALSE],
#'         group_spectra = group_spectra,
#'         chain_rules = rules$chain_fragments,
#'         dbs = dbs,
#'         paramannot = param
#'     )
#'
#'     evidence
#' }
#'
#' @return MS2 evidence table.
#'
#' @export
matchChainRules <- function(candidate,
                            group_spectra,
                            chain_rules,
                            dbs,
                            paramannot) {
  if (!nrow(group_spectra) || is.null(chain_rules) || !nrow(chain_rules)) {
    return(.empty_ms2_evidence_table())
  }

  rules_i <- chain_rules[
    chain_rules$lipid_class == candidate$class[[1L]] &
      chain_rules$polarity == paramannot$polarity,
    ,
    drop = FALSE
  ]

  if (!nrow(rules_i)) {
    return(.empty_ms2_evidence_table())
  }

  expected <- .build_expected_chain_fragments(
    candidate = candidate,
    rules = rules_i,
    dbs = dbs,
    paramannot = paramannot
  )

  .match_expected_fragments(
    candidate = candidate,
    expected = expected,
    observed = group_spectra,
    ppm = paramannot$dmzproducts,
    evidence_scope = "chain"
  )
}

#' Build expected class-level fragments
#'
#' @param candidate One-row candidate table.
#' @param rules Class-level rules.
#' @param dbs Lipid database object.
#'
#' @return Expected fragment table.
#'
#' @noRd
#'
#' @keywords internal
.build_expected_class_fragments <- function(candidate, rules, dbs) {
  out <- vector("list", nrow(rules))
  k <- 1L

  for (i in seq_len(nrow(rules))) {
    rule <- rules[i, , drop = FALSE]
    rule_type <- rule$rule_type[[1L]]

    expected_mz <- .calculate_expected_rule_mz(
      candidate = candidate,
      chain = .candidate_total_composition(candidate),
      rule = rule,
      rule_type_col = "rule_type",
      dbs = dbs
    )

    if (length(expected_mz) != 1L || !is.finite(expected_mz)) {
      next
    }

    out[[k]] <- .expected_fragment_row(
      rule_id = rule$rule_id[[1L]],
      evidence_scope = "class",
      rule_type = rule_type,
      effect = rule$effect[[1L]],
      expected_mz = expected_mz,
      chain_role = NA_character_,
      chain_composition = NA_character_,
      chain_carbon = NA_real_,
      chain_double_bond = NA_real_,
      chain_oxygen = NA_real_,
      required = rule$required[[1L]],
      score_weight = rule$score_weight[[1L]]
    )

    k <- k + 1L
  }

  .rbind_or_empty(out[seq_len(k - 1L)], .empty_expected_fragment_table())
}

#' Build total composition pseudo-chain
#'
#' Create a pseudo-chain representation from the candidate total
#' composition (`sum_composition`).
#'
#' @param candidate One-row candidate table.
#'
#' @return One-row pseudo-chain table.
#'
#' @noRd
#'
#' @keywords internal
.candidate_total_composition <- function(candidate) {
  data.frame(
    chain_role = "sum",
    chain_composition = candidate$sum_composition[[1L]],
    carbon = candidate$carbon[[1L]],
    double_bond = candidate$double_bond[[1L]],
    oxygen = candidate$oxygen[[1L]],
    stringsAsFactors = FALSE
  )
}

#' Build expected chain-level fragments
#'
#' Build expected MS2 fragments from chain-level rules. For `building_block`
#' rules, expected m/z values are calculated from the corresponding lipid
#' database in `dbs$dbs[[building_block_class]]`. For
#' `composition_difference` rules, expected m/z values are calculated from the
#' reference-class fragment. The chain composition assigned to the evidence is
#' the difference between the candidate chain composition and the reference
#' fragment composition.
#'
#' @param candidate One-row MS1 candidate table.
#' @param rules Chain-level fragmentation rules.
#' @param dbs Object returned by [createLipidDbs()].
#' @param paramannot Annotation parameter list.
#'
#' @return Expected fragment table.
#'
#' @noRd
#'
#' @keywords internal
.build_expected_chain_fragments <- function(candidate,
                                            rules,
                                            dbs,
                                            paramannot = list()) {
  if (!nrow(rules)) {
    return(.empty_expected_fragment_table())
  }

  out <- vector("list", nrow(rules))
  k <- 1L

  for (i in seq_len(nrow(rules))) {
    rule <- rules[i, , drop = FALSE]

    expected_i <- .build_expected_chain_fragments_one_rule(
      candidate = candidate,
      rule = rule,
      dbs = dbs
    )

    if (!nrow(expected_i)) {
      next
    }

    out[[k]] <- expected_i
    k <- k + 1L
  }

  .rbind_or_empty(out[seq_len(k - 1L)], .empty_expected_fragment_table())
}

#' Build expected chain fragments for one rule
#'
#' @param candidate One-row MS1 candidate table.
#' @param rule One-row chain-level rule.
#' @param dbs Object returned by [createLipidDbs()].
#'
#' @return Expected fragment table.
#'
#' @noRd
#'
#' @keywords internal
.build_expected_chain_fragments_one_rule <- function(candidate, rule, dbs) {
  evidence_type <- rule$evidence_type[[1L]]
  mass_shift <- .get_rule_numeric(rule, "mass_shift", 0)

  if (identical(evidence_type, "fragment")) {
    mz <- .get_rule_numeric(rule, "mz", NA_real_) + mass_shift

    if (!is.finite(mz)) {
      return(.empty_expected_fragment_table())
    }

    return(.expected_fragment_row(
      rule_id = rule$rule_id[[1L]],
      evidence_scope = "chain",
      rule_type = "fragment",
      effect = "support",
      expected_mz = mz,
      chain_role = rule$chain_role[[1L]],
      chain_composition = NA_character_,
      chain_carbon = NA_real_,
      chain_double_bond = NA_real_,
      chain_oxygen = NA_real_,
      required = rule$required[[1L]],
      score_weight = rule$score_weight[[1L]]
    ))
  }

  if (identical(evidence_type, "neutral_loss")) {
    return(.build_expected_chain_neutral_loss(
      candidate = candidate,
      rule = rule,
      dbs = dbs,
      mass_shift = mass_shift
    ))
  }

  if (identical(evidence_type, "building_block")) {
    return(.build_expected_chain_building_block(
      candidate = candidate,
      rule = rule,
      dbs = dbs,
      mass_shift = mass_shift
    ))
  }

  if (identical(evidence_type, "composition_difference")) {
    return(.build_expected_chain_composition_difference(
      candidate = candidate,
      rule = rule,
      dbs = dbs,
      mass_shift = mass_shift
    ))
  }

  .empty_expected_fragment_table()
}


#' Build expected fragments from building-block rules
#'
#' @param candidate One-row MS1 candidate table.
#' @param rule One-row chain-level rule.
#' @param dbs Object returned by [createLipidDbs()].
#' @param mass_shift Numeric mass shift.
#'
#' @return Expected fragment table.
#'
#' @noRd
#'
#' @keywords internal
.build_expected_chain_building_block <- function(candidate,
                                                 rule,
                                                 dbs,
                                                 mass_shift = 0) {
  bb <- .get_lipid_db_table(dbs, rule$building_block_class[[1L]])

  if (!nrow(bb)) {
    return(.empty_expected_fragment_table())
  }

  bb <- .filter_composition_subset(
    x = bb,
    total = candidate[, c("carbon", "double_bond", "oxygen"), drop = FALSE]
  )

  if (!nrow(bb)) {
    return(.empty_expected_fragment_table())
  }

  expected_mz <- calculateAdductMz(
    neutral_exact_mass = bb$neutral_exact_mass,
    adduct = rule$building_block_adduct[[1L]]
  ) + mass_shift

  keep <- is.finite(expected_mz)

  if (!any(keep)) {
    return(.empty_expected_fragment_table())
  }

  bb <- bb[keep, , drop = FALSE]
  expected_mz <- expected_mz[keep]

  data.frame(
    rule_id = rule$rule_id[[1L]],
    evidence_scope = "chain",
    rule_type = "building_block",
    effect = "support",
    expected_mz = expected_mz,
    chain_role = rule$chain_role[[1L]],
    chain_composition = bb$sum_composition,
    chain_carbon = bb$carbon,
    chain_double_bond = bb$double_bond,
    chain_oxygen = bb$oxygen,
    required = rule$required[[1L]],
    score_weight = rule$score_weight[[1L]],
    stringsAsFactors = FALSE
  )
}


#' Build expected fragments from composition-difference rules
#'
#' @param candidate One-row MS1 candidate table.
#' @param rule One-row chain-level rule.
#' @param dbs Object returned by [createLipidDbs()].
#' @param mass_shift Numeric mass shift.
#'
#' @return Expected fragment table.
#'
#' @noRd
#'
#' @keywords internal
.build_expected_chain_composition_difference <- function(candidate,
                                                         rule,
                                                         dbs,
                                                         mass_shift = 0) {
  ref <- .get_lipid_db_table(dbs, rule$reference_class[[1L]])

  if (!nrow(ref)) {
    return(.empty_expected_fragment_table())
  }

  fa <- .get_lipid_db_table(dbs, "FA")

  if (!nrow(fa)) {
    return(.empty_expected_fragment_table())
  }

  ref <- .filter_composition_subset(
    x = ref,
    total = candidate[, c("carbon", "double_bond", "oxygen"), drop = FALSE]
  )

  if (!nrow(ref)) {
    return(.empty_expected_fragment_table())
  }

  diff <- data.frame(
    carbon = candidate$carbon[[1L]] - ref$carbon,
    double_bond = candidate$double_bond[[1L]] - ref$double_bond,
    oxygen = candidate$oxygen[[1L]] - ref$oxygen,
    stringsAsFactors = FALSE
  )

  diff$chain_composition <- .chain_composition_to_string(diff)

  keep <- diff$carbon >= 0 &
    diff$double_bond >= 0 &
    diff$oxygen >= 0 &
    diff$chain_composition %in% fa$sum_composition

  if (!any(keep)) {
    return(.empty_expected_fragment_table())
  }

  ref <- ref[keep, , drop = FALSE]
  diff <- diff[keep, , drop = FALSE]

  expected_mz <- calculateAdductMz(
    neutral_exact_mass = ref$neutral_exact_mass,
    adduct = rule$building_block_adduct[[1L]]
  ) + mass_shift

  keep <- is.finite(expected_mz)

  if (!any(keep)) {
    return(.empty_expected_fragment_table())
  }

  diff <- diff[keep, , drop = FALSE]
  expected_mz <- expected_mz[keep]

  data.frame(
    rule_id = rule$rule_id[[1L]],
    evidence_scope = "chain",
    rule_type = "composition_difference",
    effect = "support",
    expected_mz = expected_mz,
    chain_role = rule$chain_role[[1L]],
    chain_composition = diff$chain_composition,
    chain_carbon = diff$carbon,
    chain_double_bond = diff$double_bond,
    chain_oxygen = diff$oxygen,
    required = rule$required[[1L]],
    score_weight = rule$score_weight[[1L]],
    stringsAsFactors = FALSE
  )
}

#' Build expected fragments from neutral-loss rules
#'
#' @param candidate One-row MS1 candidate table.
#' @param rule One-row chain-level rule.
#' @param dbs Object returned by [createLipidDbs()].
#' @param mass_shift Numeric mass shift.
#'
#' @return Expected fragment table.
#'
#' @noRd
#'
#' @keywords internal
.build_expected_chain_neutral_loss <- function(candidate,
                                               rule,
                                               dbs,
                                               mass_shift = 0) {
  neutral_loss <- .get_rule_numeric(rule, "neutral_loss", NA_real_)

  if (is.finite(neutral_loss)) {
    expected_mz <- candidate$theoretical_mz[[1L]] - neutral_loss + mass_shift

    if (!is.finite(expected_mz)) {
      return(.empty_expected_fragment_table())
    }

    return(.expected_fragment_row(
      rule_id = rule$rule_id[[1L]],
      evidence_scope = "chain",
      rule_type = "neutral_loss",
      effect = "support",
      expected_mz = expected_mz,
      chain_role = rule$chain_role[[1L]],
      chain_composition = NA_character_,
      chain_carbon = NA_real_,
      chain_double_bond = NA_real_,
      chain_oxygen = NA_real_,
      required = rule$required[[1L]],
      score_weight = rule$score_weight[[1L]]
    ))
  }

  bb_class <- rule$building_block_class[[1L]]

  if (is.na(bb_class) || !nzchar(bb_class)) {
    return(.empty_expected_fragment_table())
  }

  bb <- .get_lipid_db_table(dbs, bb_class)

  if (!nrow(bb)) {
    return(.empty_expected_fragment_table())
  }

  bb <- .filter_composition_subset(
    x = bb,
    total = candidate[, c("carbon", "double_bond", "oxygen"), drop = FALSE]
  )

  if (!nrow(bb)) {
    return(.empty_expected_fragment_table())
  }

  expected_mz <- candidate$theoretical_mz[[1L]] -
    bb$neutral_exact_mass + mass_shift

  keep <- is.finite(expected_mz)

  if (!any(keep)) {
    return(.empty_expected_fragment_table())
  }

  bb <- bb[keep, , drop = FALSE]
  expected_mz <- expected_mz[keep]

  data.frame(
    rule_id = rule$rule_id[[1L]],
    evidence_scope = "chain",
    rule_type = "neutral_loss",
    effect = "support",
    expected_mz = expected_mz,
    chain_role = rule$chain_role[[1L]],
    chain_composition = bb$sum_composition,
    chain_carbon = bb$carbon,
    chain_double_bond = bb$double_bond,
    chain_oxygen = bb$oxygen,
    required = rule$required[[1L]],
    score_weight = rule$score_weight[[1L]],
    stringsAsFactors = FALSE
  )
}

#' Combine observed chains
#'
#' Generate combinations across chain roles and keep combinations whose total
#' carbon, double bond and oxygen counts match the candidate sum composition.
#'
#' @param observed_chains Observed chain table.
#' @param total One-row total composition table.
#'
#' @return Chain-combination table.
#'
#' @noRd
#'
#' @keywords internal
.combine_observed_chains <- function(observed_chains, total) {
  if (!nrow(observed_chains)) {
    return(.empty_chain_annotation_table())
  }

  roles <- unique(observed_chains$chain_role)
  by_role <- split(observed_chains, observed_chains$chain_role)

  if (!length(by_role)) {
    return(.empty_chain_annotation_table())
  }

  idx_grid <- expand.grid(
    lapply(by_role, function(z) seq_len(nrow(z))),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  if (!nrow(idx_grid)) {
    return(.empty_chain_annotation_table())
  }

  out <- vector("list", nrow(idx_grid))
  k <- 1L

  for (i in seq_len(nrow(idx_grid))) {
    combo <- do.call(rbind, Map(
      function(tab, idx) tab[idx, , drop = FALSE],
      by_role,
      as.list(idx_grid[i, , drop = TRUE])
    ))

    if (sum(combo$carbon) != total$carbon[[1L]]) {
      next
    }

    if (sum(combo$double_bond) != total$double_bond[[1L]]) {
      next
    }

    if (sum(combo$oxygen) != total$oxygen[[1L]]) {
      next
    }

    out[[k]] <- data.frame(
      chain_role = combo$chain_role,
      chain_composition = combo$chain_composition,
      carbon = combo$carbon,
      double_bond = combo$double_bond,
      oxygen = combo$oxygen,
      evidence_status = "observed",
      stringsAsFactors = FALSE
    )

    k <- k + 1L
  }

  .rbind_or_empty(out[seq_len(k - 1L)], .empty_chain_annotation_table())
}

#' Combine chains by required roles
#'
#' Generate combinations with exactly one chain per required role and keep only
#' combinations whose summed composition equals the candidate total.
#'
#' @param chains Chain table with `chain_role`, `chain_composition`, `carbon`,
#' `double_bond`, and `oxygen`.
#' @param roles Required roles in output order.
#' @param total One-row total composition table.
#'
#' @return Chain annotation table. Each valid combination is represented by
#' one row per chain and linked by `combination_id`.
#'
#' @noRd
#'
#' @keywords internal
.combine_chains_by_roles <- function(chains, roles, total) {
  chains <- chains[chains$chain_role %in% roles, , drop = FALSE]

  if (!nrow(chains)) {
    return(.empty_chain_annotation_table())
  }

  by_role <- lapply(roles, function(role) {
    chains[chains$chain_role == role, , drop = FALSE]
  })
  names(by_role) <- roles

  if (any(vapply(by_role, nrow, integer(1)) == 0L)) {
    return(.empty_chain_annotation_table())
  }

  idx_grid <- expand.grid(
    lapply(by_role, function(z) seq_len(nrow(z))),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  out <- vector("list", nrow(idx_grid))
  k <- 1L

  for (i in seq_len(nrow(idx_grid))) {
    combo <- do.call(rbind, Map(
      function(tab, idx) tab[idx, , drop = FALSE],
      by_role,
      as.list(idx_grid[i, , drop = TRUE])
    ))

    combo <- combo[match(roles, combo$chain_role), , drop = FALSE]

    if (sum(combo$carbon) != total$carbon[[1L]]) next
    if (sum(combo$double_bond) != total$double_bond[[1L]]) next
    if (sum(combo$oxygen) != total$oxygen[[1L]]) next

    positional_id <- paste(combo$chain_composition, collapse = "/")
    unordered_id <- paste(sort(combo$chain_composition), collapse = "_")

    out[[k]] <- data.frame(
      chain_role = combo$chain_role,
      chain_composition = combo$chain_composition,
      carbon = combo$carbon,
      double_bond = combo$double_bond,
      oxygen = combo$oxygen,
      evidence_status = "observed",
      combination_id = paste0("C", k),
      positional_annotation = positional_id,
      unordered_annotation = unordered_id,
      stringsAsFactors = FALSE
    )

    k <- k + 1L
  }

  out <- .rbind_or_empty(out[seq_len(k - 1L)], .empty_chain_annotation_table())

  if (!nrow(out)) {
    return(out)
  }

  out
}


#' Get a lipid database table
#'
#' Extract and standardize one lipid database from the `createLipidDbs()` output.
#'
#' @param dbs Object returned by [createLipidDbs()].
#' @param lipid_class Lipid class name.
#'
#' @return Standardized lipid database table.
#'
#' @noRd
#'
#' @keywords internal
.get_lipid_db_table <- function(dbs, lipid_class) {
  if (is.null(lipid_class) || is.na(lipid_class) || !nzchar(lipid_class)) {
    return(.empty_lipid_db_table())
  }

  if (!is.null(dbs$dbs) && lipid_class %in% names(dbs$dbs)) {
    x <- dbs$dbs[[lipid_class]]
  } else if (lipid_class %in% names(dbs)) {
    x <- dbs[[lipid_class]]
  } else {
    return(.empty_lipid_db_table())
  }

  if (is.null(x) || !is.data.frame(x) || !nrow(x)) {
    return(.empty_lipid_db_table())
  }

  .standardize_lipid_db_columns(x)
}


#' Standardize lipid database column names
#'
#' @param x Lipid database table.
#'
#' @return Standardized lipid database table.
#'
#' @noRd
#'
#' @keywords internal
.standardize_lipid_db_columns <- function(x) {
  if ("total" %in% colnames(x) && !"sum_composition" %in% colnames(x)) {
    x$sum_composition <- x$total
  }

  if ("mass" %in% colnames(x) && !"neutral_exact_mass" %in% colnames(x)) {
    x$neutral_exact_mass <- x$mass
  }

  if ("formula" %in% colnames(x) && !"molecular_formula" %in% colnames(x)) {
    x$molecular_formula <- x$formula
  }

  if ("n_carbon" %in% colnames(x) && !"carbon" %in% colnames(x)) {
    x$carbon <- x$n_carbon
  }

  if ("n_double_bond" %in% colnames(x) && !"double_bond" %in% colnames(x)) {
    x$double_bond <- x$n_double_bond
  }

  if ("n_double_bonds" %in% colnames(x) && !"double_bond" %in% colnames(x)) {
    x$double_bond <- x$n_double_bonds
  }

  if ("n_oxygen" %in% colnames(x) && !"oxygen" %in% colnames(x)) {
    x$oxygen <- x$n_oxygen
  }

  required <- c(
    "class", "sum_composition", "carbon", "double_bond", "oxygen",
    "molecular_formula", "neutral_exact_mass"
  )

  missing <- setdiff(required, colnames(x))

  if (length(missing)) {
    stop(
      "Lipid database table is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  x[, required, drop = FALSE]
}


#' Filter a composition table to rows not exceeding a total composition
#'
#' @param x Lipid database table.
#' @param total One-row total composition table.
#'
#' @return Filtered table.
#'
#' @noRd
#'
#' @keywords internal
.filter_composition_subset <- function(x, total) {
  x[
    x$carbon <= total$carbon[[1L]] &
      x$double_bond <= total$double_bond[[1L]] &
      x$oxygen <= total$oxygen[[1L]],
    ,
    drop = FALSE
  ]
}


#' Empty lipid database table
#'
#' @return Empty standardized lipid database table.
#'
#' @noRd
#'
#' @keywords internal
.empty_lipid_db_table <- function() {
  data.frame(
    class = character(),
    sum_composition = character(),
    carbon = numeric(),
    double_bond = numeric(),
    oxygen = numeric(),
    molecular_formula = character(),
    neutral_exact_mass = numeric(),
    stringsAsFactors = FALSE
  )
}

#' Calculate expected m/z for a rule
#'
#' @param candidate One-row candidate table.
#' @param chain Optional one-row chain composition table.
#' @param rule One-row rule table.
#' @param rule_type_col Name of the rule type column.
#' @param dbs Lipid database object.
#'
#' @return `numeric(1)`.
#'
#' @noRd
#'
#' @keywords internal
.calculate_expected_rule_mz <- function(candidate,
                                        chain,
                                        rule,
                                        rule_type_col,
                                        dbs) {
  rule_type <- rule[[rule_type_col]][[1L]]
  mass_shift <- .get_rule_numeric(rule, "mass_shift", 0)

  out <- NA_real_

  if (identical(rule_type, "fragment")) {
    out <- .get_rule_numeric(rule, "mz", NA_real_) + mass_shift
  } else if (identical(rule_type, "neutral_loss")) {
    out <- candidate$theoretical_mz[[1L]] -
      .get_rule_numeric(rule, "neutral_loss", NA_real_) +
      mass_shift
  } else if (identical(rule_type, "building_block")) {
    neutral_mass <- .get_building_block_neutral_mass(
      lipid_class = rule$building_block_class[[1L]],
      chain = chain,
      dbs = dbs
    )

    if (length(neutral_mass) == 1L && is.finite(neutral_mass)) {
      out <- calculateAdductMz(
        neutral_exact_mass = neutral_mass,
        adduct = rule$building_block_adduct[[1L]]
      ) + mass_shift
    }
  } else if (identical(rule_type, "composition_difference")) {
    reference_mass <- .get_building_block_neutral_mass(
      lipid_class = rule$reference_class[[1L]],
      chain = chain,
      dbs = dbs
    )

    if (length(reference_mass) == 1L && is.finite(reference_mass)) {
      neutral_mass <- candidate$neutral_exact_mass[[1L]] - reference_mass

      out <- calculateAdductMz(
        neutral_exact_mass = neutral_mass,
        adduct = rule$building_block_adduct[[1L]]
      ) + mass_shift
    }
  }

  if (length(out) != 1L) {
    return(NA_real_)
  }

  out
}

#' Get neutral mass for a building block
#'
#' @param lipid_class Building-block lipid class.
#' @param chain Optional one-row chain composition table.
#' @param dbs Lipid database object.
#'
#' @return `numeric(1)`.
#'
#' @noRd
#'
#' @keywords internal
.get_building_block_neutral_mass <- function(lipid_class, chain, dbs) {
  if (is.null(chain) || !nrow(chain)) {
    return(NA_real_)
  }

  total <- .chain_composition_to_string(chain)

  registry <- lipidMassModels()
  model <- .get_mass_model(lipid_class, registry = registry)
  res <- model(total)

  if (!"neutral_exact_mass" %in% names(res)) {
    return(NA_real_)
  }

  res$neutral_exact_mass
}

#' Get candidate-compatible chain pool
#'
#' @param candidate One-row candidate table.
#' @param rules Chain-level rules.
#' @param paramannot Annotation parameter list.
#'
#' @return `data.frame` with chain compositions.
#'
#' @noRd
#'
#' @keywords internal
.get_candidate_chain_pool <- function(candidate, rules, paramannot) {
  chain_pools <- paramannot$chain_pools %||% defaultChainPools()

  pool_names <- unique(c(
    "fa",
    "sph",
    "phyto",
    "alkyl",
    "vinyl"
  ))

  pools <- chain_pools[intersect(pool_names, names(chain_pools))]
  pools <- pools[vapply(pools, length, integer(1)) > 0L]

  if (!length(pools)) {
    return(.empty_chain_pool_table())
  }

  parsed <- lapply(pools, function(x) {
    m <- .parse_chain_vector_local(x)
    data.frame(
      chain_composition = .chain_composition_to_string(m),
      carbon = m$carbon,
      double_bond = m$double_bond,
      oxygen = m$oxygen,
      stringsAsFactors = FALSE
    )
  })

  out <- .rbind_or_empty(parsed, .empty_chain_pool_table())

  out <- out[
    out$carbon <= candidate$carbon[[1L]] &
      out$double_bond <= candidate$double_bond[[1L]] &
      out$oxygen <= candidate$oxygen[[1L]],
    ,
    drop = FALSE
  ]

  out <- unique(out)
  rownames(out) <- NULL
  out
}

#' Match expected fragments against observed fragments
#'
#' @param candidate One-row candidate table.
#' @param expected Expected fragment table.
#' @param observed Observed fragment table.
#' @param ppm Product ion tolerance in ppm.
#' @param evidence_scope Evidence scope.
#'
#' @return MS2 evidence table.
#'
#' @noRd
#'
#' @keywords internal
.match_expected_fragments <- function(candidate,
                                      expected,
                                      observed,
                                      ppm,
                                      evidence_scope) {
  if (!nrow(expected) || !nrow(observed)) {
    return(.empty_ms2_evidence_table())
  }

  out <- vector("list", nrow(expected))
  k <- 1L

  for (i in seq_len(nrow(expected))) {
    mz_i <- expected$expected_mz[[i]]

    if (!is.finite(mz_i)) {
      next
    }

    tol <- mz_i * ppm * 1e-6
    idx <- which(abs(observed$mz - mz_i) <= tol)

    if (!length(idx)) {
      next
    }

    best <- idx[which.max(observed$intensity[idx])]

    out[[k]] <- data.frame(
      evidence_id = paste0(
        candidate$adduct_group_id[[1L]], "_",
        expected$rule_id[[i]], "_",
        k
      ),
      adduct_group_id = candidate$adduct_group_id[[1L]],
      candidate_peak_id = candidate$chrom_peak_id[[1L]],
      ms2_source_peak_id = observed$chrom_peak_id[[best]],
      ms2_source_spectrum_id = observed$spectrum_id[[best]],
      sample = candidate$sample[[1L]],
      lipid_class = candidate$class[[1L]],
      sum_composition = candidate$sum_composition[[1L]],
      evidence_scope = evidence_scope,
      rule_id = expected$rule_id[[i]],
      rule_type = expected$rule_type[[i]],
      effect = expected$effect[[i]],
      expected_mz = mz_i,
      observed_mz = observed$mz[[best]],
      ppm_error = (observed$mz[[best]] - mz_i) / mz_i * 1e6,
      intensity = observed$intensity[[best]],
      chain_role = expected$chain_role[[i]],
      chain_composition = expected$chain_composition[[i]],
      chain_carbon = expected$chain_carbon[[i]],
      chain_double_bond = expected$chain_double_bond[[i]],
      chain_oxygen = expected$chain_oxygen[[i]],
      supports_annotation = expected$effect[[i]] == "support",
      required = expected$required[[i]],
      score_weight = expected$score_weight[[i]],
      stringsAsFactors = FALSE
    )

    k <- k + 1L
  }

  .rbind_or_empty(out[seq_len(k - 1L)], .empty_ms2_evidence_table())
}

#' Build valid chain combinations
#'
#' Build chain combinations using the expected chain-role architecture for the
#' candidate class and polarity. Missing roles are inferred only if explicitly
#' allowed by `inference_rules`.
#'
#' @param candidate One-row MS1 candidate table.
#' @param chain_evidence MS2 chain evidence table.
#' @param inference_rules Chain-inference rules.
#' @param combination_rules Chain-combination rules.
#' @param paramannot Annotation parameter list.
#' @param dbs Optional lipid database object used to constrain inferred
#'   chains. Usually returned by [createLipidDbs()].
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
#' group_ids <- unique(stats::na.omit(candidates$adduct_group_id))
#'
#' if (length(group_ids)) {
#'     group_candidates <- candidates[
#'         candidates$adduct_group_id == group_ids[1],
#'         ,
#'         drop = FALSE
#'     ]
#'
#'     group_spectra <- collectAdductGroupSpectra(
#'         group_candidates,
#'         ms2_data
#'     )
#'
#'     candidate <- group_candidates[1, , drop = FALSE]
#'
#'     chain_evidence <- matchChainRules(
#'         candidate,
#'         group_spectra,
#'         rules$chain_fragments,
#'         dbs,
#'         param
#'     )
#'
#'     combinations <- buildChainCombinations(
#'         candidate,
#'         chain_evidence,
#'         paramannot = param,
#'         dbs = dbs
#'     )
#'
#'     combinations
#' }
#'
#' @return `data.frame` with valid chain combinations.
#'
#' @export
buildChainCombinations <- function(candidate,
                                   chain_evidence,
                                   inference_rules = defaultChainInferenceRules(),
                                   combination_rules = defaultChainCombinationRules(),
                                   paramannot = list(),
                                   dbs = NULL) {
  roles <- .get_combination_roles(
    lipid_class = candidate$class[[1L]],
    polarity = paramannot$polarity,
    combination_rules = combination_rules
  )

  if (!length(roles) || !nrow(chain_evidence)) {
    return(.empty_chain_annotation_table())
  }

  observed <- chain_evidence[
    chain_evidence$evidence_scope == "chain" &
      !is.na(chain_evidence$chain_composition) &
      nzchar(chain_evidence$chain_composition),
    ,
    drop = FALSE
  ]

  if (!nrow(observed)) {
    return(.empty_chain_annotation_table())
  }

  observed <- unique(observed[, c(
    "chain_role", "chain_composition",
    "chain_carbon", "chain_double_bond", "chain_oxygen"
  ), drop = FALSE])

  colnames(observed) <- c(
    "chain_role", "chain_composition",
    "carbon", "double_bond", "oxygen"
  )

  observed <- observed[observed$chain_role %in% roles, , drop = FALSE]

  total <- candidate[, c("carbon", "double_bond", "oxygen"), drop = FALSE]

  direct <- .combine_chains_by_roles(
    chains = observed,
    roles = roles,
    total = total
  )

  if (nrow(direct)) {
    direct$evidence_status <- "observed"
    return(direct)
  }

  inferred <- inferMissingChains(
    lipid_class = candidate$class[[1L]],
    polarity = paramannot$polarity,
    total = total,
    observed_chains = observed,
    inference_rules = inference_rules,
    dbs = dbs
  )

  if (!nrow(inferred)) {
    return(.empty_chain_annotation_table())
  }

  inferred$chain_composition <- .chain_composition_to_string(inferred)

  chains <- unique(rbind(
    observed[, c(
      "chain_role", "chain_composition", "carbon",
      "double_bond", "oxygen"
    ), drop = FALSE],
    inferred[, c(
      "chain_role", "chain_composition", "carbon",
      "double_bond", "oxygen"
    ), drop = FALSE]
  ))

  out <- .combine_chains_by_roles(
    chains = chains,
    roles = roles,
    total = total
  )

  if (nrow(out)) {
    out$evidence_status <- "observed_or_inferred"
  }

  out
}

#' Match relative intensity rules
#'
#' Match positional relative-intensity rules against observed chain fragments.
#'
#' @param candidate One-row MS1 lipid candidate table.
#' @param chain_sets `data.frame` containing valid chain combinations generated
#'   by [buildChainCombinations()].
#' @param chain_evidence MS2 evidence table containing observed chain fragments.
#' @param intensity_rules Relative-intensity rules, typically returned by
#'   [defaultRelativeIntensityRules()].
#' @param paramannot Annotation parameter list. Must contain the ionization
#'   `polarity`.
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
#' group_ids <- unique(stats::na.omit(candidates$adduct_group_id))
#'
#' if (length(group_ids)) {
#'     group_candidates <- candidates[
#'         candidates$adduct_group_id == group_ids[1],
#'         ,
#'         drop = FALSE
#'     ]
#'
#'     candidate <- group_candidates[1, , drop = FALSE]
#'
#'     group_spectra <- collectAdductGroupSpectra(
#'         group_candidates,
#'         ms2_data
#'     )
#'
#'     chain_evidence <- matchChainRules(
#'         candidate,
#'         group_spectra,
#'         rules$chain_fragments,
#'         dbs,
#'         param
#'     )
#'
#'     chain_sets <- buildChainCombinations(
#'         candidate,
#'         chain_evidence,
#'         paramannot = param,
#'         dbs = dbs
#'     )
#'
#'     intensity_evidence <- matchRelativeIntensityRules(
#'         candidate,
#'         chain_sets,
#'         chain_evidence,
#'         rules$relative_intensity,
#'         param
#'     )
#'
#'     intensity_evidence
#' }
#'
#' @return MS2 evidence table with one row per supported positional
#' combination.
#'
#' @export
matchRelativeIntensityRules <- function(candidate,
                                        chain_sets,
                                        chain_evidence,
                                        intensity_rules,
                                        paramannot = list()) {
  if (is.null(intensity_rules) || !nrow(intensity_rules) ||
      !nrow(chain_sets) || !nrow(chain_evidence)) {
    return(.empty_ms2_evidence_table())
  }

  rules_i <- intensity_rules[
    intensity_rules$lipid_class == candidate$class[[1L]] &
      intensity_rules$polarity == paramannot$polarity,
    ,
    drop = FALSE
  ]

  if (!nrow(rules_i)) {
    return(.empty_ms2_evidence_table())
  }

  observed <- .prepare_intensity_chain_evidence(chain_evidence)

  if (!nrow(observed)) {
    return(.empty_ms2_evidence_table())
  }

  combo_ids <- unique(chain_sets$combination_id)
  combo_ids <- combo_ids[nzchar(combo_ids)]

  out <- vector("list", nrow(rules_i) * length(combo_ids))
  k <- 1L

  for (i in seq_len(nrow(rules_i))) {
    rule <- rules_i[i, , drop = FALSE]
    numerator <- .parse_intensity_operand(rule$numerator[[1L]])
    denominator <- .parse_intensity_operand(rule$denominator[[1L]])

    if (!nrow(numerator) || !nrow(denominator)) {
      next
    }

    for (combo_id in combo_ids) {
      combo <- chain_sets[chain_sets$combination_id == combo_id, , drop = FALSE]

      num_chain <- .get_combo_chain(combo, numerator$chain_role[[1L]])
      den_chain <- .get_combo_chain(combo, denominator$chain_role[[1L]])

      if (!nrow(num_chain) || !nrow(den_chain)) {
        next
      }

      num_ev <- .get_operand_evidence(observed, numerator, num_chain)
      den_ev <- .get_operand_evidence(observed, denominator, den_chain)

      if (!nrow(num_ev) || !nrow(den_ev)) {
        next
      }

      ratio <- num_ev$intensity[[1L]] / den_ev$intensity[[1L]]

      if (!is.finite(ratio) || ratio < rule$min_ratio[[1L]]) {
        next
      }

      out[[k]] <- .intensity_evidence_row(
        candidate = candidate,
        rule = rule,
        combo_id = combo_id,
        ratio = ratio,
        numerator_evidence = num_ev
      )
      k <- k + 1L
    }
  }

  .rbind_or_empty(out[seq_len(k - 1L)], .empty_ms2_evidence_table())
}

#' Prepare chain evidence for relative-intensity matching
#'
#' @param chain_evidence chain fragments found
#'
#' @noRd
#'
#' @keywords internal
.prepare_intensity_chain_evidence <- function(chain_evidence) {
  observed <- chain_evidence[
    chain_evidence$evidence_scope == "chain" &
      chain_evidence$effect == "support" &
      is.finite(chain_evidence$intensity) &
      chain_evidence$intensity > 0 &
      !is.na(chain_evidence$chain_role) &
      nzchar(chain_evidence$chain_role) &
      !is.na(chain_evidence$chain_composition) &
      nzchar(chain_evidence$chain_composition),
    ,
    drop = FALSE
  ]

  if (!nrow(observed)) {
    return(observed)
  }

  observed$fragment_label <- .fragment_label_from_chain_rule_id(
    observed$rule_id
  )

  observed[
    !is.na(observed$fragment_label) & nzchar(observed$fragment_label),
    ,
    drop = FALSE
  ]
}

#' Extract the fragment label encoded in a chain rule id
#'
#' @noRd
#'
#' @keywords internal
.fragment_label_from_chain_rule_id <- function(rule_id) {
  vapply(rule_id, function(id) {
    parts <- strsplit(id, "_", fixed = TRUE)[[1L]]

    if (length(parts) >= 5L &&
        identical(parts[[2L]], "building") &&
        identical(parts[[3L]], "block")) {
      return(parts[[5L]])
    }

    NA_character_
  }, character(1), USE.NAMES = FALSE)
}

#' Parse a relative-intensity rule operand
#'
#' @noRd
#'
#' @keywords internal
.parse_intensity_operand <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    return(.empty_intensity_operand_table())
  }

  parts <- strsplit(x, "_", fixed = TRUE)[[1L]]

  if (length(parts) < 2L) {
    return(.empty_intensity_operand_table())
  }

  data.frame(
    fragment_label = paste(parts[-length(parts)], collapse = "_"),
    chain_role = parts[[length(parts)]],
    stringsAsFactors = FALSE
  )
}

#' Get the chain assigned to a role in one combination
#'
#' @noRd
#'
#' @keywords internal
.get_combo_chain <- function(combo, chain_role) {
  out <- combo[combo$chain_role == chain_role, , drop = FALSE]

  if (!nrow(out)) {
    return(out)
  }

  out[1L, , drop = FALSE]
}

#' Get aggregated evidence for an intensity-rule operand
#'
#' Retrieve all evidence rows matching the fragment type, chain role and chain
#' composition specified by an intensity-rule operand.
#'
#' Intensities are summed across matching evidence rows. This allows different
#' variants of the same fragment type, spectra or precursor adducts to
#' contribute jointly to the intensity assigned to one chain, while preventing
#' unrelated fragment types such as neutral losses, fatty-acid fragments or
#' lysophospholipid fragments from being combined.
#'
#' The most intense matching evidence row is retained as the representative
#' source row, and its `intensity` value is replaced by the total summed
#' intensity.
#'
#' @param observed Prepared chain-level evidence returned by
#' `.prepare_intensity_chain_evidence()`.
#' @param operand One-row intensity operand table returned by
#' `.parse_intensity_operand()`.
#' @param chain One-row chain-combination table.
#'
#' @return One-row evidence table with the summed intensity, or an empty table
#' if no matching evidence is available.
#'
#' @noRd
#'
#' @keywords internal
.get_operand_evidence <- function(observed, operand, chain) {
  out <- observed[
    observed$fragment_label == operand$fragment_label[[1L]] &
      observed$chain_role == operand$chain_role[[1L]] &
      observed$chain_composition == chain$chain_composition[[1L]],
    ,
    drop = FALSE
  ]

  if (!nrow(out)) {
    return(out)
  }

  keep <- is.finite(out$intensity) & out$intensity > 0
  out <- out[keep, , drop = FALSE]

  if (!nrow(out)) {
    return(out)
  }

  total_intensity <- sum(out$intensity, na.rm = TRUE)

  representative <- out[
    which.max(out$intensity),
    ,
    drop = FALSE
  ]

  representative$intensity <- total_intensity
  representative
}

#' Build relative-intensity evidence row
#'
#' @noRd
#'
#' @keywords internal
.intensity_evidence_row <- function(candidate,
                                    rule,
                                    combo_id,
                                    ratio,
                                    numerator_evidence) {
  data.frame(
    evidence_id = paste0(
      candidate$adduct_group_id[[1L]], "_",
      rule$rule_id[[1L]], "_",
      combo_id
    ),
    adduct_group_id = candidate$adduct_group_id[[1L]],
    candidate_peak_id = candidate$chrom_peak_id[[1L]],
    ms2_source_peak_id = numerator_evidence$ms2_source_peak_id[[1L]],
    ms2_source_spectrum_id = numerator_evidence$ms2_source_spectrum_id[[1L]],
    sample = candidate$sample[[1L]],
    lipid_class = candidate$class[[1L]],
    sum_composition = candidate$sum_composition[[1L]],
    evidence_scope = "relative_intensity",
    rule_id = rule$rule_id[[1L]],
    rule_type = "relative_intensity",
    effect = "support",
    expected_mz = NA_real_,
    observed_mz = NA_real_,
    ppm_error = NA_real_,
    intensity = ratio,
    chain_role = NA_character_,
    chain_composition = combo_id,
    chain_carbon = NA_real_,
    chain_double_bond = NA_real_,
    chain_oxygen = NA_real_,
    supports_annotation = TRUE,
    required = rule$required[[1L]],
    score_weight = rule$score_weight[[1L]],
    stringsAsFactors = FALSE
  )
}

#' Empty relative-intensity operand table
#'
#' @noRd
#'
#' @keywords internal
.empty_intensity_operand_table <- function() {
  data.frame(
    fragment_label = character(),
    chain_role = character(),
    stringsAsFactors = FALSE
  )
}

#' Select a uniquely determined positional chain combination
#'
#' Determine whether the available chain combinations correspond to a single
#' unique positional annotation.
#'
#' This helper is intended for cases where one or more chains have been inferred
#' from the lipid sum composition. If only one positional combination is
#' possible after combining observed and inferred chains, the annotation can be
#' reported using positional notation (e.g. `PC(16:0/18:1)`) even in the absence
#' of relative-intensity evidence, because no alternative positional isomer
#' remains.
#'
#' If multiple positional combinations are possible (e.g.
#' `PC(16:0/18:1)` and `PC(18:1/16:0)`), this function returns an empty table and
#' the ambiguity must instead be resolved using relative-intensity rules or
#' reported as a non-positional annotation.
#'
#' @param chain_sets `data.frame` returned by [buildChainCombinations()].
#'
#' @return A subset of `chain_sets` corresponding to the unique positional
#'   combination if exactly one exists. Otherwise, an empty `data.frame` with
#'   the same columns.
#'
#' @noRd
#'
#' @keywords internal
.select_unique_positional_chain_set <- function(chain_sets) {
  if (is.null(chain_sets) || !nrow(chain_sets)) {
    return(chain_sets[FALSE, , drop = FALSE])
  }

  required <- c(
    "combination_id",
    "positional_annotation",
    "unordered_annotation"
  )

  if (!all(required %in% colnames(chain_sets))) {
    return(chain_sets[FALSE, , drop = FALSE])
  }

  combo_summary <- unique(chain_sets[, required, drop = FALSE])

  combo_summary <- combo_summary[
    !is.na(combo_summary$combination_id) &
      nzchar(combo_summary$combination_id) &
      !is.na(combo_summary$positional_annotation) &
      nzchar(combo_summary$positional_annotation),
    ,
    drop = FALSE
  ]

  if (!nrow(combo_summary)) {
    return(chain_sets[FALSE, , drop = FALSE])
  }

  positional_annotations <- unique(combo_summary$positional_annotation)

  if (length(positional_annotations) != 1L) {
    return(chain_sets[FALSE, , drop = FALSE])
  }

  selected_id <- combo_summary$combination_id[[1L]]

  chain_sets[
    chain_sets$combination_id == selected_id,
    ,
    drop = FALSE
  ]
}

#' Rank candidate annotation from MS2 evidence
#'
#' Assign the annotation level according to the available MS2 evidence.
#'
#' Annotation levels are:
#' - `"MS-only"`: candidate supported only by MS1/adduct evidence.
#' - `"Subclass"`: class-level MS2 fragments support the lipid class.
#' - `"FA"`: chain composition is supported by observed or allowed inferred
#'   chain evidence.
#' - `"FA position"`: chain composition is supported by relative-intensity
#'   evidence, by chain-role-specific fragments that directly define the
#'   position, or because all equivalent acyl positions contain identical
#'   chains and positional exchange therefore produces the same structure.
#'
#' @param candidate One-row candidate table.
#' @param class_evidence Class-level MS2 evidence table.
#' @param chain_evidence Chain-level MS2 evidence table.
#' @param chain_sets Chain-combination table.
#' @param intensity_evidence Relative-intensity evidence table.
#'
#' @return One-row MS2 annotation table.
#'
#' @noRd
#'
#' @keywords internal
.rank_candidate_annotation <- function(candidate,
                                       class_evidence,
                                       chain_evidence,
                                       chain_sets,
                                       intensity_evidence) {
  class_support <- nrow(class_evidence) > 0L &&
    any(class_evidence$effect == "support", na.rm = TRUE)

  chain_support <- nrow(chain_sets) > 0L

  single_chain_support <- .has_single_chain_support(
    candidate = candidate,
    class_evidence = class_evidence,
    chain_evidence = chain_evidence
  )

  intensity_support <- nrow(intensity_evidence) > 0L &&
    any(intensity_evidence$effect == "support", na.rm = TRUE)

  direct_position_sets <- .select_direct_positional_chain_sets(
    candidate = candidate,
    chain_sets = chain_sets,
    chain_evidence = chain_evidence
  )

  unique_position_sets <- .select_unique_positional_chain_set(chain_sets)

  unique_position_support <- nrow(unique_position_sets) > 0L

  direct_position_support <- nrow(direct_position_sets) > 0L

  score <- sum(
    class_evidence$score_weight[class_evidence$effect == "support"],
    na.rm = TRUE
  ) +
    sum(
      chain_evidence$score_weight[chain_evidence$effect == "support"],
      na.rm = TRUE
    ) +
    sum(
      intensity_evidence$score_weight[intensity_evidence$effect == "support"],
      na.rm = TRUE
    )

  if (chain_support) {
    if (intensity_support) {
      supported_combo_ids <- unique(
        intensity_evidence$chain_composition[
          intensity_evidence$effect == "support" &
            intensity_evidence$evidence_scope == "relative_intensity"
        ]
      )

      supported_combo_ids <- supported_combo_ids[
        !is.na(supported_combo_ids) &
          nzchar(supported_combo_ids)
      ]

      selected_sets <- chain_sets[
        chain_sets$combination_id %in% supported_combo_ids,
        ,
        drop = FALSE
      ]

      level <- "FA position"
    } else if (direct_position_support) {
      selected_sets <- direct_position_sets
      level <- "FA position"
    } else if (unique_position_support) {
      selected_sets <- unique_position_sets
      level <- "FA position"
    } else {
      selected_sets <- .collapse_unordered_chain_sets(chain_sets)
      level <- "FA"
    }

    annotation <- .format_chain_set_annotations(
      lipid_class = candidate$class[[1L]],
      chain_sets = selected_sets,
      positional = identical(level, "FA position")
    )

    return(.build_ms2_annotation_row(
      candidate = candidate,
      annotation_level = level,
      annotation = annotation,
      score_ms2 = score,
      status = "accepted"
    ))
  }

  if (single_chain_support) {
    return(.build_ms2_annotation_row(
      candidate = candidate,
      annotation_level = "FA",
      annotation = .format_single_chain_annotation(
        candidate = candidate,
        chain_evidence = chain_evidence
      ),
      score_ms2 = score,
      status = "accepted"
    ))
  }

  if (class_support) {
    return(.build_ms2_annotation_row(
      candidate = candidate,
      annotation_level = "Subclass",
      annotation = .format_sum_annotation(
        candidate$class[[1L]],
        candidate$sum_composition[[1L]]
      ),
      score_ms2 = score,
      status = "accepted"
    ))
  }

  .build_ms2_annotation_row(
    candidate = candidate,
    annotation_level = "MS-only",
    annotation = .format_sum_annotation(
      candidate$class[[1L]],
      candidate$sum_composition[[1L]]
    ),
    score_ms2 = 0,
    status = "ms1_only"
  )
}

#' Test whether relative-intensity rules are needed
#'
#' Relative-intensity rules are only informative when several positional
#' permutations share the same unordered chain composition and all positional
#' roles are equivalent acyl positions (`sn1`, `sn2`, ...). Classes with
#' structurally different roles, such as plasmalogens and sphingolipids, are
#' positioned by their role-specific fragments instead.
#'
#' @noRd
#'
#' @keywords internal
.needs_relative_intensity_matching <- function(chain_sets) {
  if (!nrow(chain_sets)) {
    return(FALSE)
  }

  roles <- unique(chain_sets$chain_role)

  if (length(roles) < 2L || any(!grepl("^sn[0-9]+$", roles))) {
    return(FALSE)
  }

  by_combo <- split(chain_sets, chain_sets$combination_id)

  unordered <- vapply(
    by_combo,
    function(x) x$unordered_annotation[[1L]],
    character(1)
  )

  any(duplicated(unordered))
}

#' Select chain sets directly defining positional annotation
#'
#' Sphingolipids and plasmalogens can be positioned without relative-intensity
#' rules because their chain roles are chemically different. For example,
#' `sph|fa`, `phyto|fa`, `alkyl|sn2` and `vinyl|sn2` define the chain order by
#' fragment identity.
#'
#' @noRd
#'
#' @keywords internal
.select_direct_positional_chain_sets <- function(candidate,
                                                 chain_sets,
                                                 chain_evidence) {
  if (!nrow(chain_sets)) {
    return(chain_sets)
  }

  lipid_class <- candidate$class[[1L]]

  if (!.is_role_defined_position_class(lipid_class)) {
    return(chain_sets[FALSE, , drop = FALSE])
  }

  observed <- chain_evidence[
    chain_evidence$evidence_scope == "chain" &
      chain_evidence$effect == "support" &
      !is.na(chain_evidence$chain_role) &
      nzchar(chain_evidence$chain_role) &
      !is.na(chain_evidence$chain_composition) &
      nzchar(chain_evidence$chain_composition),
    ,
    drop = FALSE
  ]

  if (!nrow(observed)) {
    return(chain_sets[FALSE, , drop = FALSE])
  }

  observed_key <- paste(observed$chain_role, observed$chain_composition)
  combo_ids <- unique(chain_sets$combination_id)

  keep_ids <- vapply(combo_ids, function(combo_id) {
    combo <- chain_sets[chain_sets$combination_id == combo_id, , drop = FALSE]
    combo_key <- paste(combo$chain_role, combo$chain_composition)

    if (.is_sphingolipid_class(lipid_class)) {
      position_roles <- combo$chain_role %in% c("sph", "phyto", "fa")
      position_key <- combo_key[position_roles]
      return(length(position_key) > 0L && all(position_key %in% observed_key))
    }

    if (.is_plasmalogen_class(lipid_class)) {
      ether_roles <- combo$chain_role %in% c("alkyl", "vinyl")
      ether_key <- combo_key[ether_roles]
      return(length(ether_key) > 0L && all(ether_key %in% observed_key))
    }

    FALSE
  }, logical(1), USE.NAMES = FALSE)

  chain_sets[chain_sets$combination_id %in% combo_ids[keep_ids], , drop = FALSE]
}

#' Test whether a class has role-defined chain positions
#'
#' @noRd
#'
#' @keywords internal
.is_role_defined_position_class <- function(lipid_class) {
  .is_sphingolipid_class(lipid_class) || .is_plasmalogen_class(lipid_class)
}

#' Test whether chain fragments directly define positional annotation
#'
#' @noRd
#'
#' @keywords internal
.has_direct_positional_chain_support <- function(candidate,
                                                 chain_sets,
                                                 chain_evidence) {
  nrow(.select_direct_positional_chain_sets(
    candidate = candidate,
    chain_sets = chain_sets,
    chain_evidence = chain_evidence
  )) > 0L
}


#' Test whether a single-chain lipid has chain-level evidence
#'
#' @noRd
#'
#' @keywords internal
.has_single_chain_support <- function(candidate, class_evidence, chain_evidence) {
  class_support <- nrow(class_evidence) > 0L &&
    any(class_evidence$effect == "support", na.rm = TRUE)

  if (!class_support || !nrow(chain_evidence)) {
    return(FALSE)
  }

  observed <- .get_single_chain_evidence(candidate, chain_evidence)
  nrow(observed) > 0L
}

#' Extract chain-level evidence matching the candidate sum composition
#'
#' @noRd
#'
#' @keywords internal
.get_single_chain_evidence <- function(candidate, chain_evidence) {
  if (!nrow(chain_evidence)) {
    return(chain_evidence[FALSE, , drop = FALSE])
  }

  out <- chain_evidence[
    chain_evidence$evidence_scope == "chain" &
      chain_evidence$effect == "support" &
      !is.na(chain_evidence$chain_composition) &
      nzchar(chain_evidence$chain_composition) &
      chain_evidence$chain_carbon == candidate$carbon[[1L]] &
      chain_evidence$chain_double_bond == candidate$double_bond[[1L]] &
      chain_evidence$chain_oxygen == candidate$oxygen[[1L]],
    ,
    drop = FALSE
  ]

  out
}

#' Format single-chain annotations from chain evidence
#'
#' @noRd
#'
#' @keywords internal
.format_single_chain_annotation <- function(candidate, chain_evidence) {
  observed <- .get_single_chain_evidence(candidate, chain_evidence)

  if (!nrow(observed)) {
    return(.format_sum_annotation(
      candidate$class[[1L]],
      candidate$sum_composition[[1L]]
    ))
  }

  chains <- unique(observed$chain_composition)
  chains <- sort(chains[nzchar(chains)])

  paste0(
    .display_lipid_class(candidate$class[[1L]]),
    "(",
    paste(chains, collapse = "|"),
    ")"
  )
}

#' Collapse positional chain sets to one row per unordered annotation
#'
#' @noRd
#'
#' @keywords internal
.collapse_unordered_chain_sets <- function(chain_sets) {
  if (!nrow(chain_sets)) {
    return(chain_sets)
  }

  keep_ids <- vapply(
    split(chain_sets, chain_sets$combination_id),
    function(x) x$combination_id[[1L]],
    character(1)
  )

  unordered <- vapply(
    split(chain_sets, chain_sets$combination_id),
    function(x) x$unordered_annotation[[1L]],
    character(1)
  )

  keep_ids <- keep_ids[!duplicated(unordered)]
  chain_sets[chain_sets$combination_id %in% keep_ids, , drop = FALSE]
}

#' Format annotations from chain-combination rows
#'
#' @noRd
#'
#' @keywords internal
.format_chain_set_annotations <- function(lipid_class,
                                          chain_sets,
                                          positional) {
  if (!nrow(chain_sets)) {
    return(NA_character_)
  }

  combo_ids <- unique(chain_sets$combination_id)

  annotations <- vapply(combo_ids, function(combo_id) {
    combo <- chain_sets[chain_sets$combination_id == combo_id, , drop = FALSE]
    combo <- .format_combo_chains(lipid_class, combo, positional)
    paste0(combo$class, "(", combo$chains, ")")
  }, character(1), USE.NAMES = FALSE)

  annotations <- unique(annotations[nzchar(annotations)])
  annotations <- sort(annotations)

  paste(annotations, collapse = "|")
}

#' Format one chain combination
#'
#' Format one valid chain combination using positional or non-positional
#' separators.
#'
#' When positional information is available, chains are separated with `"/"`.
#' Without positional evidence, chains are normally separated with `"_"`.
#'
#' As a special case, combinations consisting exclusively of equivalent acyl
#' positions (`sn1`, `sn2`, ...) with identical chain compositions are formatted
#' with `"/"` even when the overall annotation level is `"FA"`. In this case,
#' exchanging positions does not produce a different structural alternative.
#'
#' @param lipid_class `character(1)`. Internal lipid class name.
#' @param combo One chain combination.
#' @param positional `logical(1)`. Whether positional information is supported.
#'
#' @return Named list with formatted class and chain annotation.
#'
#' @noRd
#'
#' @keywords internal
.format_combo_chains <- function(lipid_class, combo, positional) {
  roles <- as.character(combo$chain_role)
  chains <- as.character(combo$chain_composition)
  names(chains) <- roles

  class_out <- .display_lipid_class(lipid_class)

  if (.is_acylcer_class(lipid_class)) {
    chains <- .format_sphingolipid_chains(
      lipid_class = lipid_class,
      roles = roles,
      chains = chains
    )

    sn1 <- chains[roles == "sn1"][[1L]]
    backbone <- paste(chains[roles != "sn1"], collapse = "/")

    return(list(
      class = class_out,
      chains = paste(sn1, backbone, sep = "-")
    ))
  }

  if (.is_sphingolipid_class(lipid_class)) {
    chains <- .format_sphingolipid_chains(
      lipid_class = lipid_class,
      roles = roles,
      chains = chains
    )

    return(list(
      class = class_out,
      chains = paste(chains, collapse = "/")
    ))
  }

  if (.is_plasmalogen_class(lipid_class)) {
    chains <- .format_plasmalogen_chains(
      lipid_class = lipid_class,
      roles = roles,
      chains = chains
    )

    return(list(
      class = class_out,
      chains = paste(chains, collapse = "/")
    ))
  }

  identical_equivalent_acyl_positions <-
    length(chains) > 1L &&
    all(grepl("^sn[0-9]+$", roles)) &&
    length(unique(chains)) == 1L

  separator <- if (
    isTRUE(positional) ||
    identical_equivalent_acyl_positions
  ) {
    "/"
  } else {
    "_"
  }

  list(
    class = class_out,
    chains = paste(chains, collapse = separator)
  )
}

#' Format sphingolipid chains with sphingoid-base prefixes
#'
#' @noRd
#'
#' @keywords internal
.format_sphingolipid_chains <- function(lipid_class, roles, chains) {
  out <- chains
  base_prefix <- if (grepl("_t$", lipid_class)) "t" else "d"
  base_roles <- roles %in% c("sph", "phyto")
  out[base_roles] <- paste0(base_prefix, out[base_roles])
  unname(out)
}

#' Format plasmalogen chains with O-/P- inside parentheses
#'
#' @noRd
#'
#' @keywords internal
.format_plasmalogen_chains <- function(lipid_class, roles, chains) {
  out <- chains
  prefix <- if (grepl("-P$", lipid_class)) "P-" else "O-"
  ether_roles <- roles %in% c("alkyl", "vinyl")
  out[ether_roles] <- paste0(prefix, out[ether_roles])
  unname(out)
}

#' Format a sum-composition annotation
#'
#' @noRd
#'
#' @keywords internal
.format_sum_annotation <- function(lipid_class, sum_composition) {
  class_out <- .display_lipid_class(lipid_class)

  if (.is_plasmalogen_class(lipid_class)) {
    prefix <- if (grepl("-P$", lipid_class)) "P-" else "O-"
    return(paste0(class_out, "(", prefix, sum_composition, ")"))
  }

  if (.is_sphingolipid_class(lipid_class)) {
    prefix <- if (grepl("_t$", lipid_class)) "t" else "d"
    return(paste0(class_out, "(", prefix, sum_composition, ")"))
  }

  paste0(class_out, "(", sum_composition, ")")
}

#' Return display class without internal chain-type suffixes
#'
#' @noRd
#'
#' @keywords internal
.display_lipid_class <- function(lipid_class) {
  lipid_class <- sub("_[dt]$", "", lipid_class)
  sub("-[OP]$", "", lipid_class)
}

#' Test whether a class is a sphingolipid class
#'
#' @noRd
#'
#' @keywords internal
.is_sphingolipid_class <- function(lipid_class) {
  grepl("^(Cer|CerP|SM|AcylCer|HexCer|AHexCer|SHexCer)_[dt]$", lipid_class)
}

#' Test whether a class is an acylceramide class
#'
#' @noRd
#'
#' @keywords internal
.is_acylcer_class <- function(lipid_class) {
  grepl("^(AcylCer|AHexCer)_[dt]$", lipid_class)
}

#' Test whether a class is an ether or vinyl ether phospholipid class
#'
#' @noRd
#'
#' @keywords internal
.is_plasmalogen_class <- function(lipid_class) {
  grepl("-[OP]$", lipid_class)
}
#' Keep best annotation per adduct group
#'
#' @param annotations Annotation table.
#'
#' @return `data.frame`.
#'
#' @noRd
#'
#' @keywords internal
.keep_best_group_annotation <- function(annotations) {
  if (!nrow(annotations)) {
    return(annotations)
  }

  rank <- c(
    rejected = 0L,
    `MS-only` = 1L,
    Subclass = 2L,
    FA = 3L,
    `FA position` = 4L
  )

  annotations$level_rank <- unname(rank[annotations$annotation_level])
  annotations$level_rank[is.na(annotations$level_rank)] <- 0L

  annotations <- annotations[
    order(
      annotations$adduct_group_id,
      -annotations$level_rank,
      -annotations$score_ms2
    ),
    ,
    drop = FALSE
  ]

  out <- annotations[!duplicated(annotations$adduct_group_id), , drop = FALSE]
  out$level_rank <- NULL
  rownames(out) <- NULL
  out
}

#' Build one MS2 annotation row
#'
#' @return One-row annotation table.
#'
#' @noRd
#'
#' @keywords internal
.build_ms2_annotation_row <- function(candidate,
                                      annotation_level,
                                      annotation,
                                      score_ms2,
                                      status) {
  data.frame(
    adduct_group_id = candidate$adduct_group_id[[1L]],
    sample = candidate$sample[[1L]],
    lipid_class = candidate$class[[1L]],
    sum_composition = candidate$sum_composition[[1L]],
    annotation_level = annotation_level,
    annotation = annotation,
    score_ms2 = score_ms2,
    best_candidate_peak_id = candidate$chrom_peak_id[[1L]],
    status = status,
    stringsAsFactors = FALSE
  )
}

#' Check whether exclusion evidence is present
#'
#' @noRd
#'
#' @keywords internal
.has_exclude_evidence <- function(evidence) {
  nrow(evidence) > 0L && any(evidence$effect == "exclude")
}

#' Check whether required class evidence is present
#'
#' @noRd
#'
#' @keywords internal
.has_required_class_evidence <- function(evidence,
                                         class_rules,
                                         candidate,
                                         paramannot) {
  rules_i <- class_rules[
    class_rules$lipid_class == candidate$class[[1L]] &
      class_rules$polarity == paramannot$polarity &
      class_rules$effect == "support",
    ,
    drop = FALSE
  ]

  if (!nrow(rules_i)) {
    return(TRUE)
  }

  required_rules <- rules_i[rules_i$required, , drop = FALSE]

  if (nrow(required_rules)) {
    return(all(required_rules$rule_id %in% evidence$rule_id))
  }

  any(evidence$effect == "support")
}

#' Convert chain composition to string
#'
#' @noRd
#'
#' @keywords internal
.chain_composition_to_string <- function(x) {
  if (!"oxygen" %in% colnames(x)) {
    x$oxygen <- 0
  }

  out <- paste(x$carbon, x$double_bond, sep = ":")
  has_oxygen <- x$oxygen > 0

  out[has_oxygen] <- paste0(out[has_oxygen], ";O", x$oxygen[has_oxygen])
  out
}

#' Parse chain vector locally
#'
#' @noRd
#'
#' @keywords internal
.parse_chain_vector_local <- function(x) {
  parsed <- lapply(x, function(z) {
    z <- trimws(z)
    m <- regexec("^([0-9]+):([0-9]+)(?:;O([0-9]*))?$", z, perl = TRUE)
    p <- regmatches(z, m)[[1L]]

    if (!length(p)) {
      stop("Invalid chain composition: '", z, "'.", call. = FALSE)
    }

    oxygen <- 0
    if (length(p) >= 4L && !is.na(p[4L]) && nzchar(p[4L])) {
      oxygen <- as.numeric(p[4L])
    } else if (grepl(";O$", z, perl = TRUE)) {
      oxygen <- 1
    }

    c(
      carbon = as.numeric(p[2L]),
      double_bond = as.numeric(p[3L]),
      oxygen = oxygen
    )
  })

  out <- as.data.frame(do.call(rbind, parsed), stringsAsFactors = FALSE)
  out$carbon <- as.numeric(out$carbon)
  out$double_bond <- as.numeric(out$double_bond)
  out$oxygen <- as.numeric(out$oxygen)
  out
}

#' Get numeric rule value
#'
#' @noRd
#'
#' @keywords internal
.get_rule_numeric <- function(rule, column, default) {
  if (!column %in% colnames(rule)) {
    return(default)
  }

  value <- rule[[column]][[1L]]

  if (is.na(value)) {
    return(default)
  }

  as.numeric(value)
}

#' Build one expected fragment row
#'
#' @noRd
#'
#' @keywords internal
.expected_fragment_row <- function(rule_id,
                                   evidence_scope,
                                   rule_type,
                                   effect,
                                   expected_mz,
                                   chain_role,
                                   chain_composition,
                                   chain_carbon,
                                   chain_double_bond,
                                   chain_oxygen,
                                   required,
                                   score_weight) {
  data.frame(
    rule_id = rule_id,
    evidence_scope = evidence_scope,
    rule_type = rule_type,
    effect = effect,
    expected_mz = expected_mz,
    chain_role = chain_role,
    chain_composition = chain_composition,
    chain_carbon = chain_carbon,
    chain_double_bond = chain_double_bond,
    chain_oxygen = chain_oxygen,
    required = required,
    score_weight = score_weight,
    stringsAsFactors = FALSE
  )
}

#' Validate inputs for fragment matching
#'
#' @noRd
#'
#' @keywords internal
.validate_fragment_matching_inputs <- function(candidates, paramannot) {
  required <- c(
    "adduct_group_id", "chrom_peak_id", "sample", "class",
    "sum_composition", "carbon", "double_bond", "oxygen",
    "theoretical_mz", "neutral_exact_mass"
  )

  missing <- setdiff(required, colnames(candidates))
  if (length(missing)) {
    stop(
      "'candidates' is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(paramannot$polarity)) {
    stop("'paramannot$polarity' is required.", call. = FALSE)
  }

  if (is.null(paramannot$dmzproducts)) {
    stop("'paramannot$dmzproducts' is required.", call. = FALSE)
  }

  invisible(TRUE)
}

#' Bind rows or return an empty table
#'
#' @noRd
#'
#' @keywords internal
.rbind_or_empty <- function(x, empty) {
  x <- x[!vapply(x, is.null, logical(1))]
  x <- x[vapply(x, nrow, integer(1)) > 0L]

  if (!length(x)) {
    return(empty)
  }

  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

#' Empty group spectra table
#'
#' @noRd
#'
#' @keywords internal
.empty_group_spectra_table <- function() {
  data.frame(
    mz = numeric(),
    intensity = numeric(),
    chrom_peak_id = character(),
    ms1_chrom_peak_id = character(),
    ms2_chrom_peak_id = character(),
    spectrum_id = character(),
    stringsAsFactors = FALSE
  )
}

#' Empty expected fragment table
#'
#' @noRd
#'
#' @keywords internal
.empty_expected_fragment_table <- function() {
  data.frame(
    rule_id = character(),
    evidence_scope = character(),
    rule_type = character(),
    effect = character(),
    expected_mz = numeric(),
    chain_role = character(),
    chain_composition = character(),
    chain_carbon = numeric(),
    chain_double_bond = numeric(),
    chain_oxygen = numeric(),
    required = logical(),
    score_weight = numeric(),
    stringsAsFactors = FALSE
  )
}

#' Empty MS2 evidence table
#'
#' @noRd
#'
#' @keywords internal
.empty_ms2_evidence_table <- function() {
  data.frame(
    evidence_id = character(),
    adduct_group_id = character(),
    candidate_peak_id = character(),
    ms2_source_peak_id = character(),
    ms2_source_spectrum_id = character(),
    sample = integer(),
    lipid_class = character(),
    sum_composition = character(),
    evidence_scope = character(),
    rule_id = character(),
    rule_type = character(),
    effect = character(),
    expected_mz = numeric(),
    observed_mz = numeric(),
    ppm_error = numeric(),
    intensity = numeric(),
    chain_role = character(),
    chain_composition = character(),
    chain_carbon = numeric(),
    chain_double_bond = numeric(),
    chain_oxygen = numeric(),
    supports_annotation = logical(),
    required = logical(),
    score_weight = numeric(),
    stringsAsFactors = FALSE
  )
}

#' Empty chain annotation table
#'
#' @return Empty chain annotation table.
#'
#' @noRd
#'
#' @keywords internal
.empty_chain_annotation_table <- function() {
  data.frame(
    chain_role = character(),
    chain_composition = character(),
    carbon = numeric(),
    double_bond = numeric(),
    oxygen = numeric(),
    evidence_status = character(),
    combination_id = character(),
    positional_annotation = character(),
    unordered_annotation = character(),
    stringsAsFactors = FALSE
  )
}

#' Empty MS2 annotation table
#'
#' @noRd
#'
#' @keywords internal
.empty_ms2_annotation_table <- function() {
  data.frame(
    adduct_group_id = character(),
    sample = integer(),
    lipid_class = character(),
    sum_composition = character(),
    annotation_level = character(),
    annotation = character(),
    score_ms2 = numeric(),
    best_candidate_peak_id = character(),
    status = character(),
    stringsAsFactors = FALSE
  )
}

#' Empty chain pool table
#'
#' @noRd
#'
#' @keywords internal
.empty_chain_pool_table <- function() {
  data.frame(
    chain_composition = character(),
    carbon = numeric(),
    double_bond = numeric(),
    oxygen = numeric(),
    stringsAsFactors = FALSE
  )
}

#' Empty class rule table
#'
#' @noRd
#'
#' @keywords internal
.empty_class_rule_table <- function() {
  data.frame(
    lipid_class = character(),
    polarity = character(),
    rule_type = character(),
    rule_id = character(),
    mz = numeric(),
    neutral_loss = numeric(),
    building_block_class = character(),
    building_block_adduct = character(),
    effect = character(),
    required = logical(),
    score_weight = numeric(),
    stringsAsFactors = FALSE
  )
}

#' Empty chain rule table
#'
#' @noRd
#'
#' @keywords internal
.empty_chain_rule_table <- function() {
  data.frame(
    lipid_class = character(),
    polarity = character(),
    rule_id = character(),
    chain_role = character(),
    evidence_type = character(),
    mz = numeric(),
    neutral_loss = numeric(),
    building_block_class = character(),
    building_block_adduct = character(),
    reference_class = character(),
    mass_shift = numeric(),
    required = logical(),
    score_weight = numeric(),
    stringsAsFactors = FALSE
  )
}

#' Empty intensity rule table
#'
#' @noRd
#'
#' @keywords internal
.empty_intensity_rule_table <- function() {
  data.frame(stringsAsFactors = FALSE)
}

#' Test whether identical chains make positional ambiguity irrelevant
#'
#' Return `TRUE` when all valid chain combinations contain only equivalent
#' acyl positions (`sn1`, `sn2`, ...) and all chains within each combination
#' have exactly the same composition.
#'
#' For example, `PC(16:0/16:0)` has no alternative positional isomer because
#' exchanging sn1 and sn2 produces the same annotation.
#'
#' This rule is restricted to equivalent acyl positions. Chemically distinct
#' roles such as `sph`, `phyto`, `alkyl`, `vinyl` or `fa` are not evaluated
#' using this shortcut.
#'
#' @param chain_sets Chain-combination table returned by
#' [buildChainCombinations()].
#'
#' @return `logical(1)`.
#'
#' @noRd
#'
#' @keywords internal
.has_unambiguous_identical_acyl_positions <- function(chain_sets) {
  if (is.null(chain_sets) || !nrow(chain_sets)) {
    return(FALSE)
  }

  required <- c(
    "combination_id",
    "chain_role",
    "chain_composition"
  )

  if (!all(required %in% colnames(chain_sets))) {
    return(FALSE)
  }

  combo_ids <- unique(as.character(chain_sets$combination_id))
  combo_ids <- combo_ids[
    !is.na(combo_ids) &
      nzchar(combo_ids)
  ]

  if (!length(combo_ids)) {
    return(FALSE)
  }

  all(vapply(
    combo_ids,
    function(combo_id) {
      combo <- chain_sets[
        chain_sets$combination_id == combo_id,
        ,
        drop = FALSE
      ]

      roles <- as.character(combo$chain_role)
      chains <- as.character(combo$chain_composition)

      valid <- !is.na(roles) &
        nzchar(roles) &
        !is.na(chains) &
        nzchar(chains)

      roles <- roles[valid]
      chains <- chains[valid]

      if (length(chains) < 2L) {
        return(FALSE)
      }

      if (length(chains) != nrow(combo)) {
        return(FALSE)
      }

      all(grepl("^sn[0-9]+$", roles)) &&
        length(unique(chains)) == 1L
    },
    logical(1),
    USE.NAMES = FALSE
  ))
}

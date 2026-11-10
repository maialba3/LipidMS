# ==============================================================================
# Fragmentation rule utilities
# ==============================================================================

# ==============================================================================
# Class-level fragmentation and chain inference rules
# ==============================================================================

#' Return default fragmentation rules
#'
#' @examples
#' rules <- defaultFragmentRules()
#' names(rules)
#'
#' @return Named list with class, chain and relative-intensity rules.
#'
#' @export
defaultFragmentRules <- function() {
  list(
    class_fragments = defaultClassFragmentRules(),
    chain_fragments = defaultChainFragmentRules(),
    combination_rules = defaultChainCombinationRules(),
    chain_inference_rules = defaultChainInferenceRules(),
    relative_intensity = defaultRelativeIntensityRules()
  )
}

#' Make a class fragment rule identifier
#'
#' @param rule_type Character. One of `"fragment"`, `"neutral_loss"` or
#' `"building_block"`.
#' @param mz Numeric. Fragment m/z.
#' @param neutral_loss Numeric. Neutral loss mass.
#' @param building_block_class Character. Building block lipid class.
#' @param building_block_adduct Character. Building block adduct.
#' @param effect Character vector. Either `"support"` or `"exclude"`.
#'
#' @return `character(1)`.
#'
#' @noRd
#'
#' @keywords internal
.make_fragment_rule_id <- function(rule_type,
                                   mz = NA_real_,
                                   neutral_loss = NA_real_,
                                   building_block_class = NA_character_,
                                   building_block_adduct = NA_character_,
                                   effect = "support") {

  prefix <- if (effect == "exclude") "excl" else "frag"

  if (identical(rule_type, "fragment")) {
    return(paste0(
      prefix, "_",
      format(mz, digits = 8, trim = TRUE)
    ))
  }

  if (identical(rule_type, "neutral_loss")) {
    return(paste0(
      prefix, "_nl_",
      format(neutral_loss, digits = 8, trim = TRUE)
    ))
  }

  if (identical(rule_type, "building_block")) {
    return(paste(
      prefix,
      "bb",
      building_block_class,
      building_block_adduct,
      sep = "_"
    ))
  }

  stop("Unsupported rule_type: ", rule_type, call. = FALSE)
}

#' Build class-level MS2 fragmentation rules
#'
#' Build a standardized data frame for class-level fragmentation rules. Rules
#' are defined at lipid class and polarity level, independently of precursor
#' adduct. The `effect` column defines whether a matched rule supports or
#' excludes the candidate class.
#'
#' @param lipid_class Character vector with candidate lipid class names.
#' @param polarity Character vector with ionization polarity.
#' @param rule_type Character vector. One of `"fragment"`, `"neutral_loss"` or
#' `"building_block"`.
#' @param mz Numeric vector with diagnostic fragment m/z values.
#' @param neutral_loss Numeric vector with neutral loss masses.
#' @param building_block_class Character vector with building block lipid class.
#' @param building_block_adduct Character vector with building block adduct.
#' @param effect Character vector. Either `"support"` or `"exclude"`.
#' @param required Logical vector indicating whether the rule is mandatory.
#' @param score_weight Numeric vector with rule score weight.
#'
#' @return `data.frame` with standardized class-level MS2 rules.
#'
#' @noRd
#'
#' @keywords internal
.class_rules_df <- function(lipid_class,
                            polarity,
                            rule_type,
                            mz = NA_real_,
                            neutral_loss = NA_real_,
                            building_block_class = NA_character_,
                            building_block_adduct = NA_character_,
                            effect = "support",
                            required = FALSE,
                            score_weight = 1) {
  n <- max(
    length(lipid_class),
    length(polarity),
    length(rule_type),
    length(mz),
    length(neutral_loss),
    length(building_block_class),
    length(building_block_adduct),
    length(effect),
    length(required),
    length(score_weight)
  )

  out <- data.frame(
    lipid_class = rep_len(lipid_class, n),
    polarity = rep_len(polarity, n),
    rule_type = rep_len(rule_type, n),
    mz = rep_len(mz, n),
    neutral_loss = rep_len(neutral_loss, n),
    building_block_class = rep_len(building_block_class, n),
    building_block_adduct = rep_len(building_block_adduct, n),
    effect = rep_len(effect, n),
    required = rep_len(required, n),
    score_weight = rep_len(score_weight, n),
    stringsAsFactors = FALSE
  )

  allowed_rule_types <- c("fragment", "neutral_loss", "building_block")
  if (!all(out$rule_type %in% allowed_rule_types)) {
    stop(
      "'rule_type' must be one of: ",
      paste(allowed_rule_types, collapse = ", "),
      call. = FALSE
    )
  }

  allowed_effects <- c("support", "exclude")
  if (!all(out$effect %in% allowed_effects)) {
    stop(
      "'effect' must be one of: ",
      paste(allowed_effects, collapse = ", "),
      call. = FALSE
    )
  }

  invalid_fragment <- out$rule_type == "fragment" & !is.finite(out$mz)
  if (any(invalid_fragment)) {
    stop("Fragment rules require a finite 'mz'.", call. = FALSE)
  }

  invalid_neutral_loss <- out$rule_type == "neutral_loss" &
    !is.finite(out$neutral_loss)
  if (any(invalid_neutral_loss)) {
    stop("Neutral loss rules require a finite 'neutral_loss'.", call. = FALSE)
  }

  invalid_bb <- out$rule_type == "building_block" &
    (
      is.na(out$building_block_class) |
        !nzchar(out$building_block_class) |
        is.na(out$building_block_adduct) |
        !nzchar(out$building_block_adduct)
    )
  if (any(invalid_bb)) {
    stop(
      "Building block rules require 'building_block_class' and ",
      "'building_block_adduct'.",
      call. = FALSE
    )
  }

  out$rule_id <- mapply(
    .make_fragment_rule_id,
    rule_type = out$rule_type,
    mz = out$mz,
    neutral_loss = out$neutral_loss,
    building_block_class = out$building_block_class,
    building_block_adduct = out$building_block_adduct,
    effect = out$effect,
    USE.NAMES = FALSE
  )

  out[, c(
    "lipid_class", "polarity", "rule_type", "rule_id",
    "mz", "neutral_loss", "building_block_class",
    "building_block_adduct", "effect", "required", "score_weight"
  )]
}

#' Return default class fragment rules
#'
#' Class rules are defined at lipid-class and polarity level, independently of
#' the precursor adduct. Evidence from all adducts in the same
#' `adduct_group_id` should be considered during scoring.
#'
#' @examples
#' rules <- defaultClassFragmentRules()
#'
#' head(rules)
#' unique(rules$lipid_class)
#'
#' @return `data.frame`.
#'
#' @export
defaultClassFragmentRules <- function() {
  rules <- list(

    # FA and MG
    .class_rules_df("FA", "negative", "building_block",
                    building_block_class = c("FA"),
                    building_block_adduct = c("M-H", "M-H-H2O"), required = TRUE),
    .class_rules_df("MG", "positive", "building_block",
                    building_block_class = c("MG"),
                    building_block_adduct = c("M+H-H2O"), required = TRUE),

    # phospholipids: LPL/PL
    .class_rules_df("LPC", "positive", "fragment",
                    mz = c(104.1075, 184.0739), required = FALSE),
    .class_rules_df("LPC", "negative", "fragment",
                    mz = c(168.0426, 224.0688), required = FALSE),
    .class_rules_df("LPC", "negative", "building_block",
                    building_block_class = c("LPA"),
                    building_block_adduct = c("M-H"), required = FALSE),
    .class_rules_df("LPC", "negative", "building_block",
                    building_block_class = c("LPC"),
                    building_block_adduct = c("M-CH3"), required = FALSE),
    .class_rules_df("LPE", "positive", "neutral_loss",
                    neutral_loss = 141.01909, required = FALSE),
    .class_rules_df("LPE", "positive", "fragment",
                    mz = 184.0739, effect = "exclude", required = FALSE),
    .class_rules_df("LPE", "negative", "fragment",
                    mz = c(140.0115, 196.038, 214.048), required = FALSE),
    .class_rules_df("LPE", "negative", "building_block",
                    building_block_class = "LPE",
                    building_block_adduct = "M-CH3", effect = "exclude", required = FALSE),
    .class_rules_df("LPG", "negative", "fragment",
                    mz = c(152.9958, 227.0326, 209.022), required = FALSE),
    .class_rules_df("LPG", "negative", "neutral_loss",
                    neutral_loss = 74.0359, required = FALSE),
    .class_rules_df("LPI", "negative", "fragment",
                    mz = c(241.0115, 223.0008, 259.0219, 297.0375), required = FALSE),
    .class_rules_df("LPS", "negative", "neutral_loss",
                    neutral_loss = 87.032, required = FALSE),

    .class_rules_df("PC", "positive", "fragment",
                    mz = c(104.1075, 184.0739), required = FALSE),
    .class_rules_df("PC", "positive", "neutral_loss",
                    neutral_loss = 183.06604, required = FALSE),
    .class_rules_df("PC", "positive", "neutral_loss",
                    neutral_loss = 141.0193, effect = "exclude"),
    .class_rules_df("PC", "positive", "neutral_loss",
                    neutral_loss = 260.0297, effect = "exclude"),
    # .class_rules_df("PC", "positive", "building_block",
    #                 building_block_class = "DG",
    #                 building_block_adduct = "M+H-H2O", required = FALSE),
    .class_rules_df("PC", "negative", "fragment",
                    mz = c(168.0426, 224.0688), required = FALSE),
    .class_rules_df("PC", "negative", "building_block",
                    building_block_class = c("PC"),
                    building_block_adduct = c("M-CH3"), required = FALSE),
    .class_rules_df("PC", "negative", "neutral_loss",
                    neutral_loss = c(87.032), effect = "exclude"),

    .class_rules_df("PC-O", "positive", "fragment",
                    mz = c(104.1075, 184.0739), required = FALSE),
    .class_rules_df("PC-O", "positive", "neutral_loss",
                    neutral_loss = 141.0193, effect = "exclude"),
    # .class_rules_df("PC-O", "positive", "neutral_loss",
    #                 neutral_loss = c(183.06604), required = FALSE),
    .class_rules_df("PC-O", "negative", "fragment",
                    mz = c(168.0426, 224.0688), required = FALSE),
    .class_rules_df("PC-O", "negative", "building_block",
                    building_block_class = c("PC"),
                    building_block_adduct = c("M-CH3"), required = FALSE),

    .class_rules_df("PC-P", "positive", "fragment",
                    mz = c(104.1075, 184.0739), required = FALSE),
    .class_rules_df("PC-P", "positive", "neutral_loss",
                    neutral_loss = 141.0193, effect = "exclude"),
    # .class_rules_df("PC-P", "positive", "neutral_loss",
    #                 neutral_loss = c(183.06604), required = FALSE),
    .class_rules_df("PC-P", "negative", "fragment",
                    mz = c(168.0426, 224.0688), required = FALSE),
    .class_rules_df("PC-P", "negative", "building_block",
                    building_block_class = c("PC"),
                    building_block_adduct = c("M-CH3"), required = FALSE),

    .class_rules_df("PE", "positive", "building_block",
                    building_block_class = "DG",
                    building_block_adduct = "M+H-H2O", required = FALSE),
    .class_rules_df("PE", "positive", "neutral_loss",
                    neutral_loss = 141.0193, required = FALSE),
    .class_rules_df("PE", "positive", "fragment",
                    mz = c(104.1075, 184.0739), effect = "exclude"),
    .class_rules_df("PE", "positive", "neutral_loss",
                    neutral_loss = 183.06604, effect = "exclude"),
    .class_rules_df("PE", "negative", "fragment",
                    mz = c(140.0115, 196.038, 214.048), required = FALSE),
    .class_rules_df("PE", "negative", "building_block",
                    building_block_class = "PE",
                    building_block_adduct = "M-CH3", effect = "exclude", required = FALSE),

    .class_rules_df("PE-O", "positive", "neutral_loss",
                    neutral_loss = 141.0193, required = FALSE),
    .class_rules_df("PE-O", "positive", "fragment",
                    mz = c(104.1075, 184.0739), effect = "exclude", required = FALSE),
    .class_rules_df("PE-O", "negative", "fragment",
                    mz = c(140.0115, 196.038, 214.048), required = FALSE),
    .class_rules_df("PE-O", "negative", "building_block",
                    building_block_class = "PE-O",
                    building_block_adduct = "M-CH3", effect = "exclude", required = FALSE),

    .class_rules_df("PE-P", "positive", "neutral_loss",
                    neutral_loss = 141.0193, required = FALSE),
    .class_rules_df("PE-P", "positive", "fragment",
                    mz = c(104.1075, 184.0739), effect = "exclude", required = FALSE),
    .class_rules_df("PE-P", "negative", "fragment",
                    mz = c(140.0115, 196.038, 214.048), required = FALSE),
    .class_rules_df("PE-P", "negative", "building_block",
                    building_block_class = "PE-P",
                    building_block_adduct = "M-CH3", effect = "exclude", required = FALSE),

    .class_rules_df("PG", "positive", "building_block",
                    building_block_class = "DG",
                    building_block_adduct = "M+H-H2O", required = TRUE),
    .class_rules_df("PG", "negative", "fragment",
                    mz = c(227.0326, 209.022), required = FALSE), #152.9958,
    .class_rules_df("PG", "negative", "neutral_loss",
                    neutral_loss = 74.0359, required = TRUE),

    .class_rules_df("PI", "positive", "building_block",
                    building_block_class = "DG",
                    building_block_adduct = "M+H-H2O", required = FALSE),
    .class_rules_df("PI", "positive", "neutral_loss",
                    neutral_loss = 260.0297, required = FALSE),
    .class_rules_df("PI", "positive", "neutral_loss",
                    neutral_loss = 183.06604, effect = "exclude"),
    .class_rules_df("PI", "negative", "fragment",
                    mz = c(241.0115, 223.0008, 259.0219, 297.0375), required = FALSE),

    .class_rules_df("PS", "negative", "neutral_loss",
                    neutral_loss = c(87.032), required = FALSE),
    .class_rules_df("PS", "negative", "fragment",
                    mz = c(168.0426, 224.0688), effect = "exclude", required = FALSE),
    .class_rules_df("PS", "negative", "building_block", building_block_class = "PC",
                    building_block_adduct = "M-CH3", effect = "exclude"),

    # sphingolipids
    .class_rules_df("Sph_d", "positive", "building_block",
                    building_block_class = "Sph_d",
                    building_block_adduct = c("M+H-H2O", "M+H-2H2O"),
                    required = FALSE),
    .class_rules_df("Sph_d", "negative", "building_block",
                    building_block_class = "Sph_d",
                    building_block_adduct = c("M-H-H2O", "M-H-2H2O"),
                    required = FALSE),

    .class_rules_df("SphP_d", "positive", "building_block",
                    building_block_class = "SphP_d",
                    building_block_adduct = c("M+H-H2O", "M+H-2H2O", "M+H-H2O-NH4"),
                    required = FALSE),
    .class_rules_df("SphP_d", "negative", "fragment",
                    mz = c(78.9585, 96.9691), required = FALSE),
    .class_rules_df("SphP_d", "negative", "building_block",
                    building_block_class = "SphP_d",
                    building_block_adduct = "M-H-H2O",
                    required = FALSE),

    .class_rules_df("SM_d", "positive", "fragment",
                    mz = c(104.1075, 184.0739), required = FALSE),
    .class_rules_df("SM_d", "positive", "neutral_loss",
                    neutral_loss = c(183.0604), required = FALSE),
    .class_rules_df("SM_d", "negative", "fragment",
                    mz = c(168.0426, 224.0688), required = FALSE),
    .class_rules_df("SM_d", "negative", "building_block",
                    building_block_class = "SM_d",
                    building_block_adduct = "M-CH3",
                    required = FALSE),

    .class_rules_df("CerP_d", "positive", "building_block",
                    building_block_class = "Cer_d",
                    building_block_adduct = c("M+H-H2O", "M+H-2H2O"),
                    required = FALSE),
    .class_rules_df("CerP_d", "negative", "fragment",
                    mz = c(78.9585, 96.9691), required = FALSE),

    .class_rules_df("HexCer_d", "positive", "building_block",
                    building_block_class = "Cer_d",
                    building_block_adduct = c("M+H-H2O", "M+H-2H2O"),
                    required = FALSE),

    .class_rules_df("AHexCer_d", "positive", "neutral_loss",
                    neutral_loss = c(162.0528), required = FALSE),
    .class_rules_df("HexCer_d", "positive", "neutral_loss",
                    neutral_loss = c(162.0528), required = FALSE),

    # fatty acyls
    .class_rules_df("Carnitine", "positive", "fragment",
                    mz = c(60.0807, 85.0295), required = FALSE),
    .class_rules_df("Carnitine", "positive", "building_block",
                    building_block_class = "FA",
                    building_block_adduct = "M+H-H2O"),
    .class_rules_df("CE", "positive", "fragment",
                    mz = c(369.3516), required = FALSE),
    .class_rules_df("CE", "positive", "building_block",
                    building_block_class = "FA",
                    building_block_adduct = "M+H-H2O", required = FALSE)
  )

  do.call(rbind, rules)
}

# ==============================================================================
# Chain-level fragmentation rule constructors
# ==============================================================================

#' Build chain-level MS2 fragmentation rules
#'
#' Build a standardized data frame describing which MS2 signals provide direct
#' evidence for lipid chains.
#'
#' Supported evidence types are:
#' `"building_block"`, `"fragment"`, `"neutral_loss"` and
#' `"composition_difference"`.
#'
#' @param lipid_class Character vector with candidate lipid class names.
#' @param polarity Character vector with ionization polarity.
#' @param chain_role Character vector describing the role assigned to the
#' observed chain evidence, e.g. `"fa"`, `"sn1"`, `"sn2"`, `"sph"`,
#' `"phyto"`, `"alkyl"` or `"vinyl"`.
#' @param evidence_type Character vector. One of `"building_block"`,
#' `"fragment"`, `"neutral_loss"` or `"composition_difference"`.
#' @param mz Numeric vector with diagnostic fragment m/z values. Required for
#' `"fragment"` rules.
#' @param neutral_loss Numeric vector with neutral loss masses. Required for
#' `"neutral_loss"` rules.
#' @param building_block_class Character vector with the lipid class used as a
#' building block, e.g. `"FA"`, `"LPC"`, `"LPE"`, `"MG"` or `"Sph"`.
#' Required for `"building_block"` rules.
#' @param building_block_adduct Character vector with the adduct used to
#' calculate expected m/z, e.g. `"M-H"`, `"M+H"` or `"M+H-H2O"`.
#' Required for `"building_block"` and `"composition_difference"` rules.
#' @param reference_class Character vector with reference lipid class or core
#' used for `"composition_difference"` rules, e.g. `"Cer_d"` or `"Cer_t"`.
#' @param mass_shift Numeric vector with additional mass shift applied to the
#' expected m/z. Defaults to `0`.
#' @param required Logical vector indicating whether the rule is mandatory.
#' @param score_weight Numeric vector with rule score weight.
#'
#' @return `data.frame` with standardized chain-level MS2 rules.
#'
#' @noRd
#'
#' @keywords internal
.chain_rules_df <- function(lipid_class,
                            polarity,
                            chain_role,
                            evidence_type = "building_block",
                            mz = NA_real_,
                            neutral_loss = NA_real_,
                            building_block_class = NA_character_,
                            building_block_adduct = NA_character_,
                            reference_class = NA_character_,
                            mass_shift = 0,
                            required = FALSE,
                            score_weight = 1) {
  n <- max(
    length(lipid_class),
    length(polarity),
    length(chain_role),
    length(evidence_type),
    length(mz),
    length(neutral_loss),
    length(building_block_class),
    length(building_block_adduct),
    length(reference_class),
    length(mass_shift),
    length(required),
    length(score_weight)
  )

  out <- data.frame(
    lipid_class = rep_len(lipid_class, n),
    polarity = rep_len(polarity, n),
    chain_role = rep_len(chain_role, n),
    evidence_type = rep_len(evidence_type, n),
    mz = as.numeric(rep_len(mz, n)),
    neutral_loss = as.numeric(rep_len(neutral_loss, n)),
    building_block_class = rep_len(building_block_class, n),
    building_block_adduct = rep_len(building_block_adduct, n),
    reference_class = rep_len(reference_class, n),
    mass_shift = as.numeric(rep_len(mass_shift, n)),
    required = rep_len(required, n),
    score_weight = as.numeric(rep_len(score_weight, n)),
    stringsAsFactors = FALSE
  )

  allowed_evidence <- c(
    "building_block", "fragment", "neutral_loss", "composition_difference"
  )
  if (!all(out$evidence_type %in% allowed_evidence)) {
    stop(
      "'evidence_type' must be one of: ",
      paste(allowed_evidence, collapse = ", "),
      call. = FALSE
    )
  }

  invalid_fragment <- out$evidence_type == "fragment" &
    !is.finite(out$mz)
  if (any(invalid_fragment)) {
    stop("Fragment chain rules require a finite 'mz'.", call. = FALSE)
  }

  neutral_loss_rule <- out$evidence_type == "neutral_loss"

  has_fixed_neutral_loss <- is.finite(out$neutral_loss)

  has_dynamic_neutral_loss <- !is.na(out$building_block_class) &
    nzchar(out$building_block_class) &
    !is.na(out$building_block_adduct) &
    nzchar(out$building_block_adduct)

  invalid_neutral_loss <- neutral_loss_rule &
    !has_fixed_neutral_loss &
    !has_dynamic_neutral_loss

  if (any(invalid_neutral_loss)) {
    stop(
      "Neutral-loss chain rules require either a finite 'neutral_loss' or ",
      "'building_block_class' and 'building_block_adduct'.",
      call. = FALSE
    )
  }

  invalid_building_block <- out$evidence_type == "building_block" &
    (
      is.na(out$building_block_class) |
        !nzchar(out$building_block_class) |
        is.na(out$building_block_adduct) |
        !nzchar(out$building_block_adduct)
    )
  if (any(invalid_building_block)) {
    stop(
      "Building-block chain rules require 'building_block_class' and ",
      "'building_block_adduct'.",
      call. = FALSE
    )
  }

  invalid_composition_difference <- out$evidence_type ==
    "composition_difference" &
    (
      is.na(out$reference_class) |
        !nzchar(out$reference_class) |
        is.na(out$building_block_adduct) |
        !nzchar(out$building_block_adduct)
    )
  if (any(invalid_composition_difference)) {
    stop(
      "Composition-difference chain rules require 'reference_class' and ",
      "'building_block_adduct'.",
      call. = FALSE
    )
  }

  if (anyNA(out$mass_shift)) {
    stop("'mass_shift' must be numeric and non-missing.", call. = FALSE)
  }

  if (anyNA(out$score_weight) || any(out$score_weight < 0)) {
    stop(
      "'score_weight' must contain non-negative numeric values.",
      call. = FALSE
    )
  }

  out$rule_id <- .make_chain_rule_id(out)

  out[, c(
    "lipid_class", "polarity", "rule_id",
    "chain_role", "evidence_type",
    "mz", "neutral_loss",
    "building_block_class", "building_block_adduct",
    "reference_class", "mass_shift",
    "required", "score_weight"
  )]
}

#' Make chain fragmentation rule identifiers
#'
#' Build deterministic identifiers for chain-level fragmentation rules.
#'
#' @param rules `data.frame` produced by .chain_rules_df().
#'
#' @return Character vector with rule identifiers.
#'
#' @noRd
#'
#' @keywords internal
.make_chain_rule_id <- function(rules) {
  required <- c(
    "chain_role", "evidence_type", "mz", "neutral_loss",
    "building_block_class", "building_block_adduct",
    "reference_class", "mass_shift"
  )

  missing <- setdiff(required, colnames(rules))
  if (length(missing)) {
    stop(
      "'rules' is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  out <- character(nrow(rules))

  idx <- rules$evidence_type == "fragment"
  out[idx] <- paste(
    "chain", "fragment", rules$chain_role[idx],
    format(rules$mz[idx], digits = 8, trim = TRUE),
    sep = "_"
  )

  idx <- rules$evidence_type == "neutral_loss" &
    is.finite(rules$neutral_loss)
  out[idx] <- paste(
    "chain", "neutral_loss", rules$chain_role[idx],
    format(rules$neutral_loss[idx], digits = 8, trim = TRUE),
    sep = "_"
  )

  idx <- rules$evidence_type == "neutral_loss" &
    !is.finite(rules$neutral_loss)
  out[idx] <- paste(
    "chain", "neutral_loss", rules$chain_role[idx],
    rules$building_block_class[idx],
    rules$building_block_adduct[idx],
    format(rules$mass_shift[idx], digits = 8, trim = TRUE),
    sep = "_"
  )

  idx <- rules$evidence_type == "building_block"
  out[idx] <- paste(
    "chain", "building_block", rules$chain_role[idx],
    rules$building_block_class[idx],
    rules$building_block_adduct[idx],
    sep = "_"
  )

  idx <- rules$evidence_type == "composition_difference"
  out[idx] <- paste(
    "chain", "composition_difference", rules$chain_role[idx],
    rules$reference_class[idx],
    rules$building_block_adduct[idx],
    format(rules$mass_shift[idx], digits = 8, trim = TRUE),
    sep = "_"
  )

  if (any(!nzchar(out))) {
    stop("Could not build rule identifiers for some chain rules.",
         call. = FALSE)
  }

  out
}

#' Return default chain-level MS2 fragmentation rules
#'
#' Return default rules used to detect lipid chains from MS2 evidence.
#'
#' These rules are defined at lipid-class and polarity level, independently of
#' the precursor adduct. Evidence from all precursor adducts within the same
#' `adduct_group_id` should be considered during scoring.
#'
#' @examples
#' rules <- defaultChainFragmentRules()
#'
#' head(rules)
#'
#' @return `data.frame` with chain-level MS2 rules.
#'
#' @export
defaultChainFragmentRules <- function() {
  rules <- list(

    # phospholipids: LPL/PL
    .chain_rules_df(
      lipid_class = c("LPA", "LPE", "LPG", "LPI", "LPS", "LPC"),
      polarity = "negative",
      chain_role = "fa",
      building_block_class = "FA",
      building_block_adduct = "M-H"
    ),

    .chain_rules_df(
      lipid_class = c("LPC", "LPE"),
      polarity = "positive",
      chain_role = "fa",
      building_block_class = "MG",
      building_block_adduct = "M+H-H2O"
    ),

    .chain_rules_df(
      lipid_class = "PC",
      polarity = "positive",
      chain_role = "sn1",
      building_block_class = "LPC",
      building_block_adduct = c("M+H", "M+H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "PC-O",
      polarity = "positive",
      chain_role = "alkyl",
      building_block_class = "LPC-O",
      building_block_adduct = c("M+H", "M+H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "PC-P",
      polarity = "positive",
      chain_role = "vinyl",
      building_block_class = "LPC-P",
      building_block_adduct = c("M+H", "M+H-H2O")
    ),

    .chain_rules_df(
      lipid_class = c("PC", "PC-O", "PC-P"),
      polarity = "positive",
      chain_role = "sn2",
      building_block_class = "LPC",
      building_block_adduct = c("M+H", "M+H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "PC",
      polarity = "positive",
      chain_role = "sn2",
      evidence_type = "composition_difference",
      reference_class = "LPC",
      building_block_adduct = c("M+H", "M+H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "PC-O",
      polarity = "positive",
      chain_role = "sn2",
      evidence_type = "composition_difference",
      reference_class = "LPC-O",
      building_block_adduct = c("M+H", "M+H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "PC-P",
      polarity = "positive",
      chain_role = "sn2",
      evidence_type = "composition_difference",
      reference_class = "LPC-P",
      building_block_adduct = c("M+H", "M+H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "PE",
      polarity = "positive",
      chain_role = "sn1",
      building_block_class = c("LPE", "LPE", "MG"),
      building_block_adduct = c("M+H", "M+H-H2O", "M+H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "PE",
      polarity = "positive",
      chain_role = "sn2",
      building_block_class = c("LPE", "LPE", "MG"),
      building_block_adduct = c("M+H", "M+H-H2O", "M+H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "PE-O",
      polarity = "positive",
      chain_role = "alkyl",
      building_block_class = "LPE-O",
      building_block_adduct = c("M+H", "M+H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "PE-P",
      polarity = "positive",
      chain_role = "vinyl",
      building_block_class = "LPE-P",
      building_block_adduct = c("M+H", "M+H-H2O")
    ),

    .chain_rules_df(
      lipid_class = c("PE-O", "PE-P"),
      polarity = "positive",
      chain_role = "sn2",
      building_block_class = c("MG"),
      building_block_adduct = c("M+H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "PC",
      polarity = "negative",
      chain_role = c("sn1", "sn2"),
      building_block_class = "LPC",
      building_block_adduct = c("M-CH3")
    ),

    .chain_rules_df(
      lipid_class = "PC",
      polarity = "negative",
      chain_role = c("sn1", "sn2"),
      building_block_class = "LPC",
      building_block_adduct = c("M-H")
    ),

    .chain_rules_df(
      lipid_class = "PC-O",
      polarity = "negative",
      chain_role = "alkyl",
      building_block_class = "LPC-O",
      building_block_adduct = c("M-CH3", "M-CH3-H2O")
    ),

    .chain_rules_df(
      lipid_class = "PC-O",
      polarity = "negative",
      chain_role = "sn2",
      building_block_class = "FA",
      building_block_adduct = c("M-H", "M-CO2-H")
    ),

    .chain_rules_df(
      lipid_class = "PC-P",
      polarity = "negative",
      chain_role = "vinyl",
      building_block_class = "LPC-P",
      building_block_adduct = c("M-CH3", "M-CH3-H2O")
    ),

    .chain_rules_df(
      lipid_class = "PC-P",
      polarity = "negative",
      chain_role = "sn2",
      building_block_class = "FA",
      building_block_adduct = c("M-H", "M-CO2-H")
    ),

    .chain_rules_df(
      lipid_class = "PE",
      polarity = "negative",
      chain_role = c("sn1", "sn2"),
      building_block_class = "LPE",
      building_block_adduct = "M-H"
    ),

    .chain_rules_df(
      lipid_class = "PE-O",
      polarity = "negative",
      chain_role = "alkyl",
      building_block_class = "LPE-O",
      building_block_adduct = c("M-H", "M-H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "PE-P",
      polarity = "negative",
      chain_role = "vinyl",
      building_block_class = "LPE-P",
      building_block_adduct = c("M-H", "M-H-H2O")
    ),

    .chain_rules_df(
      lipid_class = c("PE-O", "PE-P"),
      polarity = "negative",
      chain_role = "sn2",
      building_block_class = "FA",
      building_block_adduct = c("M-H")
    ),

    .chain_rules_df(
      lipid_class = "PG",
      polarity = "negative",
      chain_role = c("sn1", "sn2"),
      building_block_class = "LPG",
      building_block_adduct = "M-H"
    ),

    .chain_rules_df(
      lipid_class = "PG",
      polarity = "positive",
      chain_role = c("sn1", "sn2"),
      building_block_class = c("MG"),
      building_block_adduct = "M+H-H2O"
    ),

    .chain_rules_df(
      lipid_class = "PI",
      polarity = "negative",
      chain_role = "sn1",
      building_block_class = c("LPI", "LPA"),
      building_block_adduct = "M-H-H2O"
    ),

    .chain_rules_df(
      lipid_class = "PI",
      polarity = "negative",
      chain_role = "sn2",
      building_block_class = c("LPI", "LPA"),
      building_block_adduct = "M-H"
    ),

    .chain_rules_df(
      lipid_class = "PI",
      polarity = "positive",
      chain_role = c("sn1", "sn2"),
      building_block_class = c("LPI"),
      building_block_adduct = "M+H-H2O"
    ),

    .chain_rules_df(
      lipid_class = "PI",
      polarity = "positive",
      chain_role = c("sn1", "sn2"),
      building_block_class = c("MG"),
      building_block_adduct = "M+H-H2O"
    ),

    .chain_rules_df(
      lipid_class = "PS",
      polarity = "negative",
      chain_role = "sn1",
      building_block_class = c("LPA"),
      building_block_adduct = c("M-H", "M-H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "PS",
      polarity = "negative",
      chain_role = "sn2",
      building_block_class = c("LPA"),
      building_block_adduct = c("M-H", "M-H-H2O")
    ),

    # .chain_rules_df(
    #   lipid_class = "PS",
    #   polarity = "negative",
    #   chain_role = "sn2",
    #   evidence_type = "composition_difference",
    #   reference_class = "LPA",
    #   building_block_adduct = c("M-H", "M-H-H2O")
    # ),

    .chain_rules_df(
      lipid_class = c("PC", "PE", "PG", "PI", "PS"),
      polarity = "negative",
      chain_role = "sn2",
      building_block_class = "FA",
      building_block_adduct = "M-H"
    ),

    .chain_rules_df(
      lipid_class = "CL",
      polarity = "negative",
      chain_role = c("sn1", "sn2", "sn3", "sn4"),
      building_block_class = "LPA",
      building_block_adduct = "M-H-H2O"
    ),

    .chain_rules_df(
      lipid_class = "CL",
      polarity = "negative",
      chain_role = c("sn1", "sn2", "sn3", "sn4"),
      building_block_class = "LPA",
      building_block_adduct = "M-H"
    ),

    .chain_rules_df(
      lipid_class = "CL",
      polarity = "negative",
      chain_role = c("sn1", "sn2", "sn3", "sn4"),
      building_block_class = "FA",
      building_block_adduct = "M-H"
    ),

    # sphingolipids
    .chain_rules_df(
      lipid_class = "Cer_d",
      polarity = "positive",
      chain_role = "sph",
      building_block_class = "Sph_d",
      building_block_adduct = "M+H-2H2O"
    ),

    .chain_rules_df(
      lipid_class = "Cer_d",
      polarity = "positive",
      chain_role = "sn2",
      evidence_type = "composition_difference",
      reference_class = "Sph_d",
      building_block_adduct = "M+H-2H2O"
    ),

    .chain_rules_df(
      lipid_class = "Cer_t",
      polarity = "positive",
      chain_role = "phyto",
      building_block_class = "Sph_t",
      building_block_adduct = "M+H-2H2O"
    ),

    .chain_rules_df(
      lipid_class = "Cer_t",
      polarity = "positive",
      chain_role = "sn2",
      evidence_type = "composition_difference",
      reference_class = "Sph_t",
      building_block_adduct = "M+H-2H2O"
    ),

    .chain_rules_df(
      lipid_class = "Cer_d",
      polarity = "negative",
      chain_role = "sph",
      building_block_class = "Sph_d",
      building_block_adduct = c("M-H-H2O", "M-H-2H2O")
    ),

    .chain_rules_df(
      lipid_class = "Cer_d",
      polarity = "negative",
      chain_role = "sph",
      evidence_type = "neutral_loss",
      building_block_class = "Sph_d",
      building_block_adduct = "M",
      mass_shift = -44.05
    ),

    .chain_rules_df(
      lipid_class = "Cer_d",
      polarity = "negative",
      chain_role = "fa",
      building_block_class = "FA",
      building_block_adduct = c("M", "M-H-H2O"),
      mass_shift = c(-1.9918, 0)
    ),

    .chain_rules_df(
      lipid_class = "Cer_t",
      polarity = "negative",
      chain_role = "phyto",
      building_block_class = "Sph_t",
      building_block_adduct = c("M-H-2H2O")
    ),

    .chain_rules_df(
      lipid_class = "Cer_t",
      polarity = "negative",
      chain_role = "phyto",
      evidence_type = "neutral_loss",
      building_block_class = "Sph_t",
      building_block_adduct = "M",
      mass_shift = -44.05
    ),

    .chain_rules_df(
      lipid_class = "Cer_t",
      polarity = "negative",
      chain_role = "fa",
      building_block_class = "FA",
      building_block_adduct = c("M", "M-H-H2O"),
      mass_shift = c(-1.9918, 0)
    ),

    .chain_rules_df(
      lipid_class = "CerP_d",
      polarity = "positive",
      chain_role = "sph",
      building_block_class = "Sph_d",
      building_block_adduct = "M+H-2H2O"
    ),

    .chain_rules_df(
      lipid_class = "CerP_t",
      polarity = "positive",
      chain_role = "phyto",
      building_block_class = "Sph_t",
      building_block_adduct = "M+H-2H2O"
    ),

    .chain_rules_df(
      lipid_class = "CerP_d",
      polarity = "negative",
      chain_role = "sph",
      building_block_class = "SphP_d",
      building_block_adduct = c("M-H", "M-H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "CerP_d",
      polarity = "negative",
      chain_role = "fa",
      building_block_class = "FA",
      building_block_adduct = "M",
      mass_shift = -1.9918
    ),

    .chain_rules_df(
      lipid_class = "CerP_t",
      polarity = "negative",
      chain_role = "phyto",
      building_block_class = "SphP_t",
      building_block_adduct = c("M-H", "M-H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "CerP_t",
      polarity = "negative",
      chain_role = "fa",
      building_block_class = "FA",
      building_block_adduct = "M",
      mass_shift = -1.9918
    ),

    .chain_rules_df(
      lipid_class = "AcylCer_d",
      polarity = "negative",
      chain_role = "sn1",
      evidence_type = "composition_difference",
      reference_class = "Cer_d",
      building_block_adduct = "M-H"
    ),

    .chain_rules_df(
      lipid_class = "AcylCer_d",
      polarity = "negative",
      chain_role = "sn1",
      building_block_class = "FA",
      building_block_adduct = "M-H-H2O",
    ),

    .chain_rules_df(
      lipid_class = "AcylCer_d",
      polarity = "negative",
      chain_role = "sph",
      building_block_class = "Sph_d",
      building_block_adduct = c("M", "M-H-H2O"),
      mass_shift = c(-62.06001, 0)
    ),

    .chain_rules_df(
      lipid_class = "AcylCer_d",
      polarity = "negative",
      chain_role = "sn3",
      building_block_class = "FA",
      building_block_adduct = "M",
      mass_shift = c(-1.9918, -19.0179)
    ),

    .chain_rules_df(
      lipid_class = "AcylCer_d",
      polarity = "positive",
      chain_role = "sn1",
      evidence_type = "composition_difference",
      reference_class = "Cer_d",
      building_block_adduct = c("M+H", "M+H-H2O", "M+H-2H2O")
    ),

    .chain_rules_df(
      lipid_class = "AcylCer_d",
      polarity = "positive",
      chain_role = "sph",
      building_block_class = "Sph_d",
      building_block_adduct = c("M+H-H2O", "M+H-2H2O")
    ),

    .chain_rules_df(
      lipid_class = "AcylCer_d",
      polarity = "positive",
      chain_role = "sn3",
      building_block_class = "FA",
      building_block_adduct = "M",
      mass_shift = 0.02329
    ),

    .chain_rules_df(
      lipid_class = "AcylCer_t",
      polarity = "negative",
      chain_role = "sn1",
      evidence_type = "composition_difference",
      reference_class = "Cer_t",
      building_block_adduct = "M-H"
    ),

    .chain_rules_df(
      lipid_class = "AcylCer_t",
      polarity = "negative",
      chain_role = "sn1",
      building_block_class = "FA",
      building_block_adduct = "M-H-H2O",
    ),

    .chain_rules_df(
      lipid_class = "AcylCer_t",
      polarity = "negative",
      chain_role = "phyto",
      building_block_class = "Sph_t",
      building_block_adduct = c("M", "M-H-H2O"),
      mass_shift = c(-62.06001, 0)
    ),

    .chain_rules_df(
      lipid_class = "AcylCer_t",
      polarity = "negative",
      chain_role = "sn3",
      building_block_class = "FA",
      building_block_adduct = "M",
      mass_shift = c(-1.9918, -19.0179)
    ),

    .chain_rules_df(
      lipid_class = "AcylCer_t",
      polarity = "positive",
      chain_role = "sn1",
      evidence_type = "composition_difference",
      reference_class = "Cer_t",
      building_block_adduct = c("M+H", "M+H-H2O", "M+H-2H2O")
    ),

    .chain_rules_df(
      lipid_class = "AcylCer_t",
      polarity = "positive",
      chain_role = "phyto",
      building_block_class = "Sph_t",
      building_block_adduct = c("M+H-H2O", "M+H-2H2O")
    ),

    .chain_rules_df(
      lipid_class = "AcylCer_t",
      polarity = "positive",
      chain_role = "sn3",
      building_block_class = "FA",
      building_block_adduct = "M",
      mass_shift = 0.02329
    ),

    .chain_rules_df(
      lipid_class = "AHexCer_d",
      polarity = "positive",
      chain_role = "sn1",
      evidence_type = "composition_difference",
      reference_class = "Cer_d",
      building_block_adduct = c("M+H", "M+H-H2O", "M+H-2H2O")
    ),

    .chain_rules_df(
      lipid_class = "AHexCer_d",
      polarity = "positive",
      chain_role = "sph",
      building_block_class = "Sph_d",
      building_block_adduct = c("M+H-H2O", "M+H-2H2O")
    ),

    .chain_rules_df(
      lipid_class = "AHexCer_d",
      polarity = "positive",
      chain_role = "sn3",
      building_block_class = "FA",
      building_block_adduct = "M",
      mass_shift = 0.02329
    ),

    .chain_rules_df(
      lipid_class = "AHexCer_t",
      polarity = "positive",
      chain_role = "sn1",
      evidence_type = "composition_difference",
      reference_class = "Cer_t",
      building_block_adduct = c("M+H", "M+H-H2O", "M+H-2H2O")
    ),

    .chain_rules_df(
      lipid_class = "AHexCer_t",
      polarity = "positive",
      chain_role = "phyto",
      building_block_class = "Sph_t",
      building_block_adduct = c("M+H-H2O", "M+H-2H2O")
    ),

    .chain_rules_df(
      lipid_class = "AHexCer_t",
      polarity = "positive",
      chain_role = "sn3",
      building_block_class = "FA",
      building_block_adduct = "M",
      mass_shift = 0.02329
    ),

    .chain_rules_df(
      lipid_class = "SM_d",
      polarity = "positive",
      chain_role = "sph",
      building_block_class = "Sph_d",
      building_block_adduct = c("M+H-2H2O", "M+H-H2O")
    ),

    .chain_rules_df(
      lipid_class = "SM_d",
      polarity = "positive",
      chain_role = "sn2",
      evidence_type = "composition_difference",
      reference_class = "Sph_d",
      building_block_adduct = c("M+H-2H2O")
    ),

    .chain_rules_df(
      lipid_class = "SM_d",
      polarity = "negative",
      chain_role = "sph",
      building_block_class = "Sph_d",
      building_block_adduct = "M",
      mass_shift = 150.032
    ),

    .chain_rules_df(
      lipid_class = "SM_d",
      polarity = "negative",
      chain_role = "sph",
      building_block_class = "Sph_d",
      building_block_adduct = "M-H-H2O"
    ),

    .chain_rules_df(
      lipid_class = "SM_d",
      polarity = "negative",
      chain_role = "fa",
      building_block_class = "FA",
      building_block_adduct = "M",
      mass_shift = -1.9918
    ),

    .chain_rules_df(
      lipid_class = "SM_d",
      polarity = "negative",
      chain_role = "fa",
      evidence_type = "composition_difference",
      reference_class = "Sph",
      building_block_adduct = "M",
      mass_shift = 150.032
    ),

    .chain_rules_df(
      lipid_class = "SM_t",
      polarity = "negative",
      chain_role = "phyto",
      building_block_class = "Sph_t",
      building_block_adduct = "M",
      mass_shift = 150.032
    ),

    .chain_rules_df(
      lipid_class = "SM_t",
      polarity = "negative",
      chain_role = "phyto",
      building_block_class = "Sph_t",
      building_block_adduct = "M-H-H2O"
    ),

    .chain_rules_df(
      lipid_class = "SM_t",
      polarity = "negative",
      chain_role = "fa",
      building_block_class = "FA",
      building_block_adduct = "M",
      mass_shift = -1.9918
    ),

    .chain_rules_df(
      lipid_class = "SM_t",
      polarity = "negative",
      chain_role = "fa",
      evidence_type = "composition_difference",
      reference_class = "Sph_t",
      building_block_adduct = "M",
      mass_shift = 150.032
    ),

    # fatty acyl
    .chain_rules_df(
      lipid_class = "Carnitine",
      polarity = "positive",
      chain_role = "fa",
      building_block_class = "FA",
      building_block_adduct = "M+H-H2O"
    ),

    .chain_rules_df(
      lipid_class = "CE",
      polarity = "positive",
      chain_role = "fa",
      building_block_class = "FA",
      building_block_adduct = "M+H-H2O"
    ),

    # glycerolipids
    .chain_rules_df(
      lipid_class = "DG",
      polarity = "positive",
      chain_role = c("sn1", "sn2"),
      building_block_class = "MG",
      building_block_adduct = "M+H-H2O"
    ),

    .chain_rules_df(
      lipid_class = "TG",
      polarity = "positive",
      chain_role = c("sn1", "sn2", "sn3"),
      evidence_type = "composition_difference",
      reference_class = "DG",
      building_block_adduct = "M+H-H2O"
    )
  )

  out <- do.call(rbind, rules)
  rownames(out) <- NULL
  out
}

#' Build chain inference rules
#'
#' Build a standardized data frame describing when missing chains can be inferred
#' from the total lipid composition and the chains directly observed by MS2.
#'
#' Inference is controlled per lipid class and polarity. This is useful for
#' classes where one chain-derived signal is often absent or very weak, but the
#' missing chain can be inferred from the total composition.
#'
#' @param lipid_class Character vector with lipid class names.
#' @param polarity Character vector with ionization polarity.
#' @param chain_roles Character vector describing the expected chain roles for
#' the class, separated by `"|"`, e.g. `"sn1|sn2"` or `"fa|fa|fa"`.
#' @param allow_chain_inference Logical vector indicating whether inference is
#' allowed.
#' @param min_observed_chains Integer vector with the minimum number of directly
#' observed chains required before inference is allowed.
#' @param max_inferred_chains Integer vector with the maximum number of chains
#' that can be inferred.
#' @param inferred_chain_roles Character vector with roles allowed to be inferred,
#' separated by `"|"`. For example, `"sn2"` for PC in positive mode.
#' Use an empty string when inference is disabled.
#' @param inference_reason Character vector documenting why inference is allowed
#' or disabled.
#'
#' @return `data.frame` with chain inference rules.
#'
#' @noRd
#'
#' @keywords internal
.chain_inference_rules_df <- function(lipid_class,
                                      polarity,
                                      chain_roles,
                                      allow_chain_inference,
                                      min_observed_chains,
                                      max_inferred_chains,
                                      inferred_chain_roles = "",
                                      inferred_db_class = "",
                                      inference_reason = "") {
  n <- max(
    length(lipid_class),
    length(polarity),
    length(chain_roles),
    length(allow_chain_inference),
    length(min_observed_chains),
    length(max_inferred_chains),
    length(inferred_chain_roles),
    length(inferred_db_class),
    length(inference_reason)
  )

  out <- data.frame(
    lipid_class = rep_len(lipid_class, n),
    polarity = rep_len(polarity, n),
    chain_roles = rep_len(chain_roles, n),
    allow_chain_inference = rep_len(allow_chain_inference, n),
    min_observed_chains = as.integer(rep_len(min_observed_chains, n)),
    max_inferred_chains = as.integer(rep_len(max_inferred_chains, n)),
    inferred_chain_roles = rep_len(inferred_chain_roles, n),
    inferred_db_class = rep_len(inferred_db_class, n),
    inference_reason = rep_len(inference_reason, n),
    stringsAsFactors = FALSE
  )

  if (anyNA(out$min_observed_chains) ||
      anyNA(out$max_inferred_chains) ||
      any(out$min_observed_chains < 0L) ||
      any(out$max_inferred_chains < 0L)) {
    stop(
      "'min_observed_chains' and 'max_inferred_chains' must be non-negative integers.",
      call. = FALSE
    )
  }

  out
}

#' Return default chain-inference rules
#'
#' Define when missing chain roles can be inferred from the difference between
#' the candidate total composition and the observed chain compositions.
#'
#' @examples
#' rules <- defaultChainInferenceRules()
#'
#' head(rules)
#'
#' @return `data.frame` with chain-inference rules.
#'
#' @export
defaultChainInferenceRules <- function() {
  rules <- list(
    .chain_inference_rules_df(
      lipid_class = "PC",
      polarity = "positive",
      chain_roles = "sn1|sn2",
      allow_chain_inference = TRUE,
      min_observed_chains = 1L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "sn2",
      inferred_db_class = "FA",
      inference_reason = "minor sn2-derived fragments in positive mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "PC-O",
      polarity = "positive",
      chain_roles = "alkyl|sn2",
      allow_chain_inference = TRUE,
      min_observed_chains = 1L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "sn2",
      inferred_db_class = "FA",
      inference_reason = "minor sn2-derived fragments in positive mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "PC-P",
      polarity = "positive",
      chain_roles = "vinyl|sn2",
      allow_chain_inference = TRUE,
      min_observed_chains = 1L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "sn2",
      inferred_db_class = "FA",
      inference_reason = "minor sn2-derived fragments in positive mode"
    ),

    # .chain_inference_rules_df(
    #   lipid_class = "PS",
    #   polarity = "negative",
    #   chain_roles = "sn1|sn2",
    #   allow_chain_inference = TRUE,
    #   min_observed_chains = 1L,
    #   max_inferred_chains = 1L,
    #   inferred_chain_roles = "sn2",
    #   inferred_db_class = "FA",
    #   inference_reason = "minor sn2-derived fragments in negative mode"
    # ),

    .chain_inference_rules_df(
      lipid_class = "SM_d",
      polarity = "positive",
      chain_roles = "sph|fa",
      allow_chain_inference = TRUE,
      min_observed_chains = 1L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "fa",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in positive mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "SM_t",
      polarity = "positive",
      chain_roles = "phyto|fa",
      allow_chain_inference = TRUE,
      min_observed_chains = 1L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "fa",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in positive mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "SM_d",
      polarity = "negative",
      chain_roles = "sph|fa",
      allow_chain_inference = TRUE,
      min_observed_chains = 1L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "fa",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in negative mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "SM_t",
      polarity = "negative",
      chain_roles = "phyto|fa",
      allow_chain_inference = TRUE,
      min_observed_chains = 1L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "fa",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in negative mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "Cer_d",
      polarity = "positive",
      chain_roles = "sph|fa",
      allow_chain_inference = TRUE,
      min_observed_chains = 1L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "fa",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in positive mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "Cer_t",
      polarity = "positive",
      chain_roles = "phyto|fa",
      allow_chain_inference = TRUE,
      min_observed_chains = 1L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "fa",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in positive mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "HexCer_d",
      polarity = "positive",
      chain_roles = "sph|fa",
      allow_chain_inference = TRUE,
      min_observed_chains = 1L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "fa",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in positive mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "HexCer_t",
      polarity = "positive",
      chain_roles = "phyto|fa",
      allow_chain_inference = TRUE,
      min_observed_chains = 1L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "fa",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in positive mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "CerP_d",
      polarity = "positive",
      chain_roles = "sph|fa",
      allow_chain_inference = TRUE,
      min_observed_chains = 1L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "fa",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in positive mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "CerP_t",
      polarity = "positive",
      chain_roles = "phyto|fa",
      allow_chain_inference = TRUE,
      min_observed_chains = 1L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "fa",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in positive mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "AcylCer_d",
      polarity = "positive",
      chain_roles = "sn1|sph|sn3",
      allow_chain_inference = TRUE,
      min_observed_chains = 2L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "sn1",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in positive mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "AHexCer_d",
      polarity = "positive",
      chain_roles = "sn1|sph|sn3",
      allow_chain_inference = TRUE,
      min_observed_chains = 2L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "sn1",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in positive mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "AcylCer_t",
      polarity = "positive",
      chain_roles = "sn1|phyto|sn3",
      allow_chain_inference = TRUE,
      min_observed_chains = 2L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "sn1",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in positive mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "AHexCer_t",
      polarity = "positive",
      chain_roles = "sn1|phyto|sn3",
      allow_chain_inference = TRUE,
      min_observed_chains = 2L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "sn1",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in positive mode"
    ),

    .chain_inference_rules_df(
      lipid_class = "CL",
      polarity = "negative",
      chain_roles = "sn1|sn2|sn3|sn4",
      allow_chain_inference = TRUE,
      min_observed_chains = 3L,
      max_inferred_chains = 1L,
      inferred_chain_roles = "sn1|sn2|sn3|sn4",
      inferred_db_class = "FA",
      inference_reason = "minor fa-derived fragments in negative mode"
    ),

    .chain_inference_rules_df(
      lipid_class = c("TG", "CL"),
      polarity = c("positive", "negative"),
      chain_roles = c("sn1|sn2|sn3", "sn1|sn2|sn3|sn4"),
      allow_chain_inference = FALSE,
      min_observed_chains = c(2L, 3L),
      max_inferred_chains = 1L,
      inferred_chain_roles = c("sn3", "sn4"),
      inferred_db_class = c("FA", "FA"),
      inference_reason = c(
        "disabled by default for triacylglycerols",
        "disabled by default for cardiolipins"
      )
    )
  )

  out <- do.call(rbind, rules)
  rownames(out) <- NULL
  out
}

#' Check whether missing chains can be inferred
#'
#' @param lipid_class `character(1)`. Lipid class.
#' @param polarity `character(1)`. Ionization polarity.
#' @param observed_n `integer(1)`. Number of directly observed chains.
#' @param missing_n `integer(1)`. Number of missing chains.
#' @param inference_rules Chain inference rules returned by
#' [defaultChainInferenceRules()].
#'
#' @return `logical(1)`.
#'
#' @noRd
#'
#' @keywords internal
canInferChains <- function(lipid_class,
                           polarity,
                           observed_n,
                           missing_n,
                           inference_rules = defaultChainInferenceRules()) {
  rule <- inference_rules[
    inference_rules$lipid_class == lipid_class &
      inference_rules$polarity == polarity,
    ,
    drop = FALSE
  ]

  if (!nrow(rule)) {
    return(FALSE)
  }

  if (!isTRUE(rule$allow_chain_inference[[1L]])) {
    return(FALSE)
  }

  observed_n >= rule$min_observed_chains[[1L]] &&
    missing_n > 0L &&
    missing_n <= rule$max_inferred_chains[[1L]]
}

#' Parse chain role strings
#'
#' @param x `character(1)`. Chain roles separated by `"|"`.
#'
#' @return Character vector with chain roles.
#'
#' @noRd
#'
#' @keywords internal
.parse_chain_roles <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("'x' must be a non-empty character(1).", call. = FALSE)
  }

  strsplit(x, "\\|", perl = TRUE)[[1L]]
}

#' Get allowed inferred chain roles
#'
#' @param lipid_class `character(1)`. Lipid class.
#' @param polarity `character(1)`. Ionization polarity.
#' @param inference_rules Chain inference rules.
#'
#' @return Character vector with roles allowed for inference.
#'
#' @noRd
#'
#' @keywords internal
.get_allowed_inferred_roles <- function(lipid_class,
                                        polarity,
                                        inference_rules) {
  rule <- inference_rules[
    inference_rules$lipid_class == lipid_class &
      inference_rules$polarity == polarity,
    ,
    drop = FALSE
  ]

  if (!nrow(rule) || !isTRUE(rule$allow_chain_inference[[1L]])) {
    return(character())
  }

  .parse_chain_roles(rule$inferred_chain_roles[[1L]])
}

#' Infer one missing chain from total composition
#'
#' Infer a single missing chain by subtracting observed chains from the total
#' composition.
#'
#' @param total `data.frame` or named numeric vector with columns/elements
#' `carbon`, `double_bond`, and optionally `oxygen`.
#' @param observed_chains `data.frame` with columns `carbon`, `double_bond`,
#' and optionally `oxygen`.
#'
#' @return One-row `data.frame` with inferred composition.
#'
#' @noRd
#'
#' @keywords internal
.infer_complement_chain <- function(total, observed_chains) {
  total <- .normalize_chain_composition(total)

  if (!nrow(observed_chains)) {
    observed_chains <- data.frame(
      carbon = numeric(),
      double_bond = numeric(),
      oxygen = numeric()
    )
  }

  observed_chains <- .normalize_chain_composition(observed_chains)

  out <- data.frame(
    carbon = total$carbon - sum(observed_chains$carbon),
    double_bond = total$double_bond - sum(observed_chains$double_bond),
    oxygen = total$oxygen - sum(observed_chains$oxygen),
    stringsAsFactors = FALSE
  )

  if (any(out < 0)) {
    return(data.frame(
      carbon = numeric(),
      double_bond = numeric(),
      oxygen = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  out
}

#' Return default chain-combination rules
#'
#' Define the expected chain-role architecture for each lipid class and polarity.
#'
#' @examples
#' rules <- defaultChainCombinationRules()
#'
#' head(rules)
#'
#' @return `data.frame`.
#'
#' @export
defaultChainCombinationRules <- function() {
  data.frame(
    lipid_class = c(
      rep(c("PC", "PC-O", "PC-P",
            "PE", "PE-O", "PE-P",
            "PG", "PI"), each = 2), "PS", "PA",
      "DG", "TG", "CL",
      rep(c("Cer_d", "Cer_t",
            "CerP_d", "CerP_t",
            "SM_d", "SM_t"), each = 2),
      "AcylCer_d", "AcylCer_t", "HexCer_d", "HexCer_t"
    ),
    polarity = c(
      rep(c("positive", "negative"), 8),
      rep("negative", 2),
      rep("positive", 2), "negative",
      rep(c("positive", "negative"), 8)
    ),
    chain_roles = c(
      rep("sn1|sn2", 2), rep("alkyl|sn2", 2), rep("vinyl|sn2", 2),
      rep("sn1|sn2", 2), rep("alkyl|sn2", 2), rep("vinyl|sn2", 2),
      rep("sn1|sn2", 6),
      "sn1|sn2", "sn1|sn2|sn3", "sn1|sn2|sn3|sn4",
      rep("sph|fa", 2), rep("phyto|fa", 2),
      rep("sph|fa", 2), rep("phyto|fa", 2),
      rep("sph|fa", 2), rep("phyto|fa", 2),
      "sn1|sph|sn3", "sn1|phyto|sn3",
      "sph|fa", "phyto|fa"
    ),
    stringsAsFactors = FALSE
  )
}

#' Get chain-combination roles
#'
#' @param lipid_class Lipid class.
#' @param polarity Ionization polarity.
#' @param combination_rules Chain-combination rule table.
#'
#' @return Character vector of chain roles in expected order.
#'
#' @noRd
#'
#' @keywords internal
.get_combination_roles <- function(lipid_class,
                                   polarity,
                                   combination_rules) {
  hit <- combination_rules[
    combination_rules$lipid_class == lipid_class &
      combination_rules$polarity == polarity,
    ,
    drop = FALSE
  ]

  if (!nrow(hit)) {
    hit <- combination_rules[
      combination_rules$lipid_class == lipid_class,
      ,
      drop = FALSE
    ]
  }

  if (!nrow(hit)) {
    return(character())
  }

  strsplit(hit$chain_roles[[1L]], "\\|", perl = TRUE)[[1L]]
}


#' Normalize chain composition table
#'
#' @param x `data.frame` or named numeric vector.
#'
#' @return `data.frame` with columns `carbon`, `double_bond`, and `oxygen`.
#'
#' @noRd
#'
#' @keywords internal
.normalize_chain_composition <- function(x) {
  if (is.numeric(x) && !is.null(names(x))) {
    x <- as.data.frame(as.list(x), stringsAsFactors = FALSE)
  }

  if (!is.data.frame(x)) {
    stop("'x' must be a data.frame or a named numeric vector.", call. = FALSE)
  }

  required <- c("carbon", "double_bond")
  missing <- setdiff(required, colnames(x))

  if (length(missing)) {
    stop(
      "Missing composition columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (!"oxygen" %in% colnames(x)) {
    x$oxygen <- 0
  }

  out <- data.frame(
    carbon = as.numeric(x$carbon),
    double_bond = as.numeric(x$double_bond),
    oxygen = as.numeric(x$oxygen),
    stringsAsFactors = FALSE
  )

  if (anyNA(out)) {
    stop("Composition values must be numeric and non-missing.", call. = FALSE)
  }

  out
}

#' Infer missing chains for one candidate
#'
#' Infer missing chain compositions for a candidate if allowed by chain inference
#' rules.
#'
#' This function currently infers one missing chain at a time. Classes requiring
#' inference of multiple alternative chain sets should be handled by a higher
#' level combinatorial chain resolver.
#'
#' @param lipid_class `character(1)`. Lipid class.
#' @param polarity `character(1)`. Ionization polarity.
#' @param total `data.frame` or named numeric vector with total composition.
#' @param observed_chains `data.frame` with observed chain compositions and a
#' `chain_role` column.
#' @param inference_rules Chain inference rules.
#' @param dbs databases
#'
#' @return `data.frame` with inferred chain compositions.
#'
#' @noRd
#'
#' @keywords internal
inferMissingChains <- function(lipid_class,
                               polarity,
                               total,
                               observed_chains,
                               inference_rules = defaultChainInferenceRules(),
                               dbs = NULL) {
  empty <- function() {
    data.frame(
      chain_role = character(),
      carbon = numeric(),
      double_bond = numeric(),
      oxygen = numeric(),
      evidence_status = character(),
      stringsAsFactors = FALSE
    )
  }

  if (!nrow(observed_chains)) {
    return(empty())
  }

  observed_chains <- .normalize_chain_composition(observed_chains)
  observed_n <- nrow(observed_chains)

  rule <- inference_rules[
    inference_rules$lipid_class == lipid_class &
      inference_rules$polarity == polarity,
    ,
    drop = FALSE
  ]

  if (!nrow(rule) || !isTRUE(rule$allow_chain_inference[[1L]])) {
    return(empty())
  }

  chain_roles <- .parse_chain_roles(rule$chain_roles[[1L]])
  missing_n <- length(chain_roles) - observed_n

  if (!canInferChains(
    lipid_class = lipid_class,
    polarity = polarity,
    observed_n = observed_n,
    missing_n = missing_n,
    inference_rules = inference_rules
  )) {
    return(empty())
  }

  allowed_roles <- .get_allowed_inferred_roles(
    lipid_class = lipid_class,
    polarity = polarity,
    inference_rules = inference_rules
  )

  if (missing_n != 1L || length(allowed_roles) != 1L) {
    return(empty())
  }

  inferred <- .infer_complement_chain(
    total = total,
    observed_chains = observed_chains
  )

  if (!nrow(inferred)) {
    return(empty())
  }

  inferred_role <- allowed_roles[[1L]]

  if (!.is_inferred_chain_in_db(
    inferred = inferred,
    inferred_role = inferred_role,
    lipid_class = lipid_class,
    polarity = polarity,
    inference_rules = inference_rules,
    dbs = dbs
  )) {
    return(empty())
  }

  data.frame(
    chain_role = inferred_role,
    carbon = inferred$carbon,
    double_bond = inferred$double_bond,
    oxygen = inferred$oxygen,
    evidence_status = "inferred",
    stringsAsFactors = FALSE
  )
}

#' Build relative-intensity rules
#'
#' @param lipid_class Character vector.
#' @param polarity Character vector.
#' @param rule_id Character vector.
#' @param numerator Character vector.
#' @param denominator Character vector.
#' @param min_ratio Numeric vector.
#' @param required Logical vector.
#' @param score_weight Numeric vector.
#'
#' @return `data.frame`.
#'
#' @noRd
#'
#' @keywords internal
.intensity_rules_df <- function(lipid_class, polarity, rule_id,
                                numerator, denominator, min_ratio,
                                required = TRUE, score_weight = 1) {
  data.frame(
    lipid_class = lipid_class,
    polarity = polarity,
    rule_id = rule_id,
    numerator = numerator,
    denominator = denominator,
    min_ratio = min_ratio,
    required = required,
    score_weight = score_weight,
    stringsAsFactors = FALSE
  )
}

#' Return default relative-intensity rules
#'
#' @examples
#' rules <- defaultRelativeIntensityRules()
#'
#' head(rules)
#'
#' @return `data.frame`.
#'
#' @export
defaultRelativeIntensityRules <- function() {
  rules <- list(
    .intensity_rules_df(
      "PC", "positive",
      "LPC_sn1_over_LPC_sn2",
      "LPC_sn1", "LPC_sn2", 1.5
    ),
    .intensity_rules_df(
      "PC", "negative",
      "LPC_sn1_over_LPC_sn2",
      "LPC_sn1", "LPC_sn2", 1.5
    ),
    # .intensity_rules_df(
    #   "PC-O", "positive",
    #   "LPCO_sn1_over_LPC_sn2",
    #   "LPC-O_sn1", "LPC_sn2", 2
    # ),
    # .intensity_rules_df(
    #   "PC-O", "negative",
    #   "FA_sn2_over_LPC-O_sn1",
    #   "FA_sn2", "LPC-O_sn1", 2
    # ),
    # .intensity_rules_df(
    #   "PC-P", "positive",
    #   "LPC_sn2_over_LPC-P_sn1",
    #   "LPC_sn2", "LPC-P_sn1", 2
    # ),
    # .intensity_rules_df(
    #   "PC-P", "negative",
    #   "FA_sn2_over_LPC-P_sn1",
    #   "FA_sn2", "LPC-P_sn1", 2
    # ),
    .intensity_rules_df(
      "PE", "positive",
      "LPE_sn1_over_LPE_sn2",
      "LPE_sn1", "LPE_sn2", 1.5
    ),
    .intensity_rules_df(
      "PE", "positive",
      "MG_sn2_over_MG_sn1",
      "MG_sn2", "MG_sn1", 1.5
    ),
    .intensity_rules_df(
      "PE", "negative",
      "LPE_sn1_over_LPE_sn2",
      "LPE_sn1", "LPE_sn2", 1.5
    ),
    # .intensity_rules_df(
    #   "PE-O", "positive",
    #   "LPEO_sn1_over_MG_sn2",
    #   "LPE_sn1", "MG_sn2", 2
    # ),
    # .intensity_rules_df(
    #   "PE-O", "negative",
    #   "LPE-O_sn1_over_MG_sn2",
    #   "LPE-O_sn1", "MG_sn2", 2
    # ),
    # .intensity_rules_df(
    #   "PE-P", "positive",
    #   "MG_sn2_over_LPE-P_sn1",
    #   "MG_sn2", "LPE-P_sn1", 3
    # ),
    # .intensity_rules_df(
    #   "PE-P", "negative",
    #   "FA_sn2_over_LPE-P_sn1",
    #   "FA_sn2", "LPE-P_sn1", 2
    # ),
    .intensity_rules_df(
      "PG", "positive",
      "MG_sn2_over_MG_sn1",
      "MG_sn2", "MG_sn1", 1.5
    ),
    .intensity_rules_df(
      "PG", "negative",
      "LPG_sn1_over_LPG_sn2",
      "LPG_sn1", "LPG_sn2", 1.5
    ),
    .intensity_rules_df(
      "PI", "negative",
      "LPI_sn1_over_LPI_sn2",
      "LPI_sn1", "LPI_sn2", 1.5
    ),
    .intensity_rules_df(
      "PI", "negative",
      "LPA_sn1_over_LPA_sn2",
      "LPA_sn1", "LPA_sn2", 1.5
    ),
    .intensity_rules_df(
      "PI", "positive",
      "MG_sn1_over_MG_sn2",
      "MG_sn1", "MG_sn2", 1.5
    ),
    .intensity_rules_df(
      "PI", "positive",
      "LPI_sn1_over_LPI_sn2",
      "LPI_sn1", "LPI_sn2", 1.5
    ),
    .intensity_rules_df(
      "PS", "negative",
      "LPA_sn1_over_LPA_sn2",
      "LPA_sn1", "LPA_sn2", 1.5
    )#,
    # .intensity_rules_df(
    #   "DG", "positive",
    #   "MG_sn1_over_MG_sn2",
    #   "MG_sn1", "MG_sn2", 1.5
    # ),
    # .intensity_rules_df(
    #   "TG", "positive",
    #   "DG_loss_sn2_over_sn1",
    #   "DG_sn2", "dg_loss_sn1", 1.5
    # )
  )

  do.call(rbind, rules)
}


.is_inferred_chain_in_db <- function(inferred,
                                     inferred_role,
                                     lipid_class,
                                     polarity,
                                     inference_rules,
                                     dbs = NULL) {
  if (is.null(dbs)) {
    return(TRUE)
  }

  if (!nrow(inferred)) {
    return(FALSE)
  }

  db_class <- .get_allowed_inferred_db_class(
    lipid_class = lipid_class,
    polarity = polarity,
    inference_rules = inference_rules
  )

  if (is.na(db_class) || !nzchar(db_class)) {
    db_class <- .get_inferred_role_db_class(
      inferred_role = inferred_role,
      lipid_class = lipid_class
    )
  }

  db_tab <- .get_lipid_db_table(dbs, db_class)

  if (!nrow(db_tab)) {
    return(FALSE)
  }

  inferred_composition <- .chain_composition_to_string(inferred)

  inferred_composition %in% db_tab$sum_composition
}

.get_allowed_inferred_db_class <- function(lipid_class,
                                           polarity,
                                           inference_rules) {
  rule <- inference_rules[
    inference_rules$lipid_class == lipid_class &
      inference_rules$polarity == polarity,
    ,
    drop = FALSE
  ]

  if (!nrow(rule) || !"inferred_db_class" %in% colnames(rule)) {
    return(NA_character_)
  }

  db_class <- rule$inferred_db_class[[1L]]

  if (is.na(db_class) || !nzchar(db_class)) {
    return(NA_character_)
  }

  db_class
}

.get_inferred_role_db_class <- function(inferred_role, lipid_class = NULL) {
  if (inferred_role %in% c("fa", "sn1", "sn2", "sn3", "sn4")) {
    return("FA")
  }

  if (identical(inferred_role, "sph")) {
    return("Sph_d")
  }

  if (identical(inferred_role, "phyto")) {
    return("Sph_t")
  }

  if (identical(inferred_role, "alkyl")) {
    return("FA")
  }

  if (identical(inferred_role, "vinyl")) {
    return("FA")
  }

  if (inferred_role %in% c("dg", "DG")) {
    return("DG")
  }

  lipid_class
}

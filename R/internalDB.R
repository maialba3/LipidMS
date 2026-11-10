# ==============================================================================
# Internal lipid database generation utilities
# ==============================================================================

.atomic_mass <- c(
  C = 12.000000,
  H = 1.007825,
  N = 14.00305,
  O = 15.9949,
  P = 30.97393,
  S = 31.97207
)

# ------------------------------------------------------------------------------
# Parsing helpers
# ------------------------------------------------------------------------------

#' Parse a lipid chain or sum composition
#'
#' Parse a lipid composition string in lipidomics shorthand format:
#' `"n_carbon:n_double_bonds"`, `"n_carbon:n_double_bonds;O"`
#' `"n_carbon:n_double_bonds;O2"`.
#'
#' @param x `character(1)`. Chain or sum composition, e.g. `"18:1"`,
#' `"18:1;O2"` or `"34:1"`.
#'
#' @return Named numeric vector with elements `n_carbon`, `n_double_bonds`
#' and `n_oxygen`.
#'
#' @noRd
#'
#' @keywords internal
.parse_chain <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("'x' must be a non-empty character(1).", call. = FALSE)
  }

  x <- trimws(x)
  pattern <- "^([0-9]+):([0-9]+)(?:;O([0-9]*))?$"
  parts <- regmatches(x, regexec(pattern, x, perl = TRUE))[[1L]]

  if (length(parts) == 0L) {
    stop(
      paste(
        "'x' must be in one of these formats:",
        "'n_carbon:n_double_bonds', 'n_carbon:n_double_bonds;O',",
        "or 'n_carbon:n_double_bonds;O2'."
      ),
      call. = FALSE
    )
  }

  n_carbon <- as.numeric(parts[2L])
  n_double_bonds <- as.numeric(parts[3L])

  n_oxygen <- 0
  if (length(parts) >= 4L && !is.na(parts[4L]) && nzchar(parts[4L])) {
    n_oxygen <- as.numeric(parts[4L])
  } else if (grepl(";O$", x, perl = TRUE)) {
    n_oxygen <- 1
  }

  values <- c(
    n_carbon = n_carbon,
    n_double_bonds = n_double_bonds,
    n_oxygen = n_oxygen
  )

  if (anyNA(values) || any(values < 0)) {
    stop("'x' must contain non-negative numeric values.", call. = FALSE)
  }

  values
}

#' Parse a vector of lipid chain compositions
#'
#' @param x Character vector of chain compositions.
#'
#' @return Numeric matrix with columns `n_carbon`, `n_double_bonds` and
#' `n_oxygen`.
#'
#' @noRd
#'
#' @keywords internal
.parse_chain_vector <- function(x) {
  if (!is.character(x) || !length(x)) {
    stop("'x' must be a non-empty character vector.", call. = FALSE)
  }

  parsed <- lapply(x, .parse_chain)
  out <- do.call(rbind, parsed)
  storage.mode(out) <- "numeric"
  colnames(out) <- c("n_carbon", "n_double_bonds", "n_oxygen")
  out
}

#' Validate chain vector
#'
#' @param x Character vector of lipid chain compositions.
#' @param arg Argument name used in error messages.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @noRd
#'
#' @keywords internal
.validate_chain_vector <- function(x, arg = "x") {
  if (!is.character(x) || !length(x)) {
    stop("'", arg, "' must be a non-empty character vector.", call. = FALSE)
  }
  invisible(.parse_chain_vector(x))
  invisible(TRUE)
}

#' Validate chain pools
#'
#' @param chain_pools Named list of chain pools.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @noRd
#'
#' @keywords internal
.validate_chain_pools <- function(chain_pools) {
  if (!is.list(chain_pools) || !length(chain_pools) ||
      is.null(names(chain_pools))) {
    stop("'chain_pools' must be a named list.", call. = FALSE)
  }

  for (nm in names(chain_pools)) {
    if (!length(chain_pools[[nm]])) {
      next
    }
    .validate_chain_vector(chain_pools[[nm]], arg = nm)
  }

  invisible(TRUE)
}

#' Parse all chain pools into numeric matrices
#'
#' @param chain_pools Named list of chain pools as character vectors.
#'
#' @return Named list of numeric matrices with columns `n_carbon`,
#' `n_double_bonds` and `n_oxygen`.
#'
#' @noRd
#'
#' @keywords internal
.parse_chain_pools <- function(chain_pools) {
  .validate_chain_pools(chain_pools)

  lapply(chain_pools, function(x) {
    if (!length(x)) {
      out <- matrix(numeric(0), nrow = 0L, ncol = 3L)
      colnames(out) <- c("n_carbon", "n_double_bonds", "n_oxygen")
      return(out)
    }
    .parse_chain_vector(x)
  })
}

# ------------------------------------------------------------------------------
# Sum composition helpers
# ------------------------------------------------------------------------------

#' Convert numeric sum compositions to lipid shorthand strings
#'
#' @param x Numeric matrix with columns `n_carbon`, `n_double_bonds` and
#' `n_oxygen`.
#'
#' @return Character vector in the format `"x:y"` or `"x:y;Oz"`.
#'
#' @noRd
#'
#' @keywords internal
.sum_compositions_to_string <- function(x) {
  if (!is.matrix(x) || ncol(x) != 3L) {
    stop("'x' must be a numeric matrix with 3 columns.", call. = FALSE)
  }

  out <- paste(x[, "n_carbon"], x[, "n_double_bonds"], sep = ":")

  has_oxygen <- x[, "n_oxygen"] > 0
  if (any(has_oxygen)) {
    out[has_oxygen] <- paste0(
      out[has_oxygen],
      ";O",
      x[has_oxygen, "n_oxygen"]
    )
  }

  out
}

#' Convert sum compositions to a numeric table
#'
#' @param sum_compositions Character vector of lipid sum compositions.
#'
#' @return `data.frame` with columns `sum_composition`, `n_carbon`,
#' `n_double_bonds` and `n_oxygen`.
#'
#' @noRd
#'
#' @keywords internal
.sum_compositions_to_dataframe <- function(sum_compositions) {
  parsed <- lapply(sum_compositions, .parse_chain)

  data.frame(
    sum_composition = sum_compositions,
    n_carbon = vapply(parsed, `[[`, numeric(1), "n_carbon"),
    n_double_bonds = vapply(parsed, `[[`, numeric(1), "n_double_bonds"),
    n_oxygen = vapply(parsed, `[[`, numeric(1), "n_oxygen"),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Formula and exact mass helpers
# ------------------------------------------------------------------------------

#' Build molecular formula string
#'
#' @param elements Named numeric vector with element counts.
#'
#' @return `character(1)` molecular formula.
#'
#' @noRd
#'
#' @keywords internal
.build_formula <- function(elements) {
  elements <- elements[elements != 0]
  paste0(names(elements), unname(elements), collapse = "")
}

#' Calculate neutral exact mass
#'
#' @param elements Named numeric vector with element counts.
#'
#' @return `numeric(1)` neutral exact mass.
#'
#' @noRd
#'
#' @keywords internal
.calculate_exact_mass <- function(elements) {
  sum(unname(elements) * .atomic_mass[names(elements)])
}

#' Calculate molecular formula and neutral exact mass
#'
#' @param sum_composition `character(1)`. Lipid sum composition.
#' @param carbon_offset `numeric(1)`. Carbon increment contributed by the
#' lipid scaffold.
#' @param hydrogen_intercept `numeric(1)`. Constant term in the hydrogen
#' formula.
#' @param n_oxygen `numeric(1)`. Oxygen atoms contributed by the lipid scaffold.
#' @param n_nitrogen `numeric(1)`. Nitrogen atoms.
#' @param n_phosphorus `numeric(1)`. Phosphorus atoms.
#' @param n_sulfur `numeric(1)`. Sulfur atoms.
#'
#' @return Named list with elements `molecular_formula` and
#' `neutral_exact_mass`.
#'
#' @noRd
#'
#' @keywords internal
.calculate_formula_mass <- function(sum_composition, carbon_offset = 0,
                                    hydrogen_intercept,
                                    n_oxygen = 0, n_nitrogen = 0,
                                    n_phosphorus = 0, n_sulfur = 0) {
  cb <- .parse_chain(sum_composition)

  carbon <- cb[["n_carbon"]] + carbon_offset
  hydrogen <- cb[["n_carbon"]] * 2 + hydrogen_intercept -
    cb[["n_double_bonds"]] * 2
  oxygen_total <- n_oxygen + cb[["n_oxygen"]]

  elements <- c(
    C = carbon,
    H = hydrogen,
    N = n_nitrogen,
    O = oxygen_total,
    P = n_phosphorus,
    S = n_sulfur
  )

  list(
    molecular_formula = .build_formula(elements),
    neutral_exact_mass = .calculate_exact_mass(elements)
  )
}

# ------------------------------------------------------------------------------
# Mass model registry
# ------------------------------------------------------------------------------

#' Return lipid mass models
#'
#' @return Named list of functions.
#'
#' @noRd
#'
#' @keywords internal
lipidMassModels <- function() {
  list(
    Sph_d = function(x) .calculate_formula_mass(
      x, hydrogen_intercept = 3, n_oxygen = 2, n_nitrogen = 1
    ),
    Sph_t = function(x) .calculate_formula_mass(
      x, hydrogen_intercept = 5, n_oxygen = 3, n_nitrogen = 1
    ),
    SphP_d = function(x) .calculate_formula_mass(
      x, hydrogen_intercept = 4, n_oxygen = 5, n_nitrogen = 1,
      n_phosphorus = 1
    ),
    SphP_t = function(x) .calculate_formula_mass(
      x, hydrogen_intercept = 6, n_oxygen = 6, n_nitrogen = 1,
      n_phosphorus = 1
    ),
    Cer_d = function(x) .calculate_formula_mass(
      x, hydrogen_intercept = 1, n_oxygen = 3, n_nitrogen = 1
    ),
    Cer_t = function(x) .calculate_formula_mass(
      x, hydrogen_intercept = 3, n_oxygen = 4, n_nitrogen = 1
    ),
    CerP_d = function(x) .calculate_formula_mass(
      x, hydrogen_intercept = 2, n_oxygen = 6, n_nitrogen = 1,
      n_phosphorus = 1
    ),
    CerP_t = function(x) .calculate_formula_mass(
      x, hydrogen_intercept = 4, n_oxygen = 7, n_nitrogen = 1,
      n_phosphorus = 1
    ),
    HexCer_d = function(x) .calculate_formula_mass(
      x, carbon_offset = 6, hydrogen_intercept = 11,
      n_oxygen = 8, n_nitrogen = 1
    ),
    HexCer_t = function(x) .calculate_formula_mass(
      x, carbon_offset = 6, hydrogen_intercept = 13,
      n_oxygen = 9, n_nitrogen = 1
    ),
    AHexCer_d = function(x) .calculate_formula_mass(
      x, carbon_offset = 6, hydrogen_intercept = 9,
      n_oxygen = 9, n_nitrogen = 1
    ),
    AHexCer_t = function(x) .calculate_formula_mass(
      x, carbon_offset = 6, hydrogen_intercept = 11,
      n_oxygen = 10, n_nitrogen = 1
    ),
    SHexCer_d = function(x) .calculate_formula_mass(
      x, carbon_offset = 6, hydrogen_intercept = 11,
      n_oxygen = 11, n_nitrogen = 1, n_sulfur = 1
    ),
    SHexCer_t = function(x) .calculate_formula_mass(
      x, carbon_offset = 6, hydrogen_intercept = 13,
      n_oxygen = 12, n_nitrogen = 1, n_sulfur = 1
    ),
    AcylCer_d = function(x) .calculate_formula_mass(
      x, hydrogen_intercept = -1, n_oxygen = 4, n_nitrogen = 1
    ),
    AcylCer_t = function(x) .calculate_formula_mass(
      x, hydrogen_intercept = 1, n_oxygen = 5, n_nitrogen = 1
    ),
    SM_d = function(x) .calculate_formula_mass(
      x, carbon_offset = 5, hydrogen_intercept = 13,
      n_oxygen = 6, n_nitrogen = 2, n_phosphorus = 1
    ),
    SM_t = function(x) .calculate_formula_mass(
      x, carbon_offset = 5, hydrogen_intercept = 13,
      n_oxygen = 7, n_nitrogen = 2, n_phosphorus = 1
    ),
    Carnitine = function(x) .calculate_formula_mass(
      x, carbon_offset = 7, hydrogen_intercept = 13,
      n_oxygen = 4, n_nitrogen = 1
    ),
    CE = function(x) .calculate_formula_mass(
      x, carbon_offset = 27, hydrogen_intercept = 44,
      n_oxygen = 2
    ),
    FA = function(x) .calculate_formula_mass(
      x, hydrogen_intercept = 0, n_oxygen = 2
    ),
    HFA = function(x) .calculate_formula_mass(
      x, hydrogen_intercept = 0, n_oxygen = 3
    ),
    FAHFA = function(x) .calculate_formula_mass(
      x, hydrogen_intercept = -2, n_oxygen = 4
    ),
    LPA = function(x) .calculate_formula_mass(
      x, carbon_offset = 3, hydrogen_intercept = 7,
      n_oxygen = 7, n_phosphorus = 1
    ),
    `LPA-O` = function(x) .calculate_formula_mass(
      x, carbon_offset = 3, hydrogen_intercept = 9,
      n_oxygen = 6, n_phosphorus = 1
    ),
    PA = function(x) .calculate_formula_mass(
      x, carbon_offset = 3, hydrogen_intercept = 5,
      n_oxygen = 8, n_phosphorus = 1
    ),
    LPE = function(x) .calculate_formula_mass(
      x, carbon_offset = 5, hydrogen_intercept = 12,
      n_oxygen = 7, n_nitrogen = 1, n_phosphorus = 1
    ),
    `LPE-O` = function(x) .calculate_formula_mass(
      x, carbon_offset = 5, hydrogen_intercept = 14,
      n_oxygen = 6, n_nitrogen = 1, n_phosphorus = 1
    ),
    `LPE-P` = function(x) .calculate_formula_mass(
      x, carbon_offset = 5, hydrogen_intercept = 12,
      n_oxygen = 6, n_nitrogen = 1, n_phosphorus = 1
    ),
    PE = function(x) .calculate_formula_mass(
      x, carbon_offset = 5, hydrogen_intercept = 10,
      n_oxygen = 8, n_nitrogen = 1, n_phosphorus = 1
    ),
    `PE-O` = function(x) .calculate_formula_mass(
      x, carbon_offset = 5, hydrogen_intercept = 12,
      n_oxygen = 7, n_nitrogen = 1, n_phosphorus = 1
    ),
    `PE-P` = function(x) .calculate_formula_mass(
      x, carbon_offset = 5, hydrogen_intercept = 10,
      n_oxygen = 7, n_nitrogen = 1, n_phosphorus = 1
    ),
    LPG = function(x) .calculate_formula_mass(
      x, carbon_offset = 6, hydrogen_intercept = 13,
      n_oxygen = 9, n_phosphorus = 1
    ),
    PG = function(x) .calculate_formula_mass(
      x, carbon_offset = 6, hydrogen_intercept = 11,
      n_oxygen = 10, n_phosphorus = 1
    ),
    LPI = function(x) .calculate_formula_mass(
      x, carbon_offset = 9, hydrogen_intercept = 17,
      n_oxygen = 12, n_phosphorus = 1
    ),
    PI = function(x) .calculate_formula_mass(
      x, carbon_offset = 9, hydrogen_intercept = 15,
      n_oxygen = 13, n_phosphorus = 1
    ),
    LPS = function(x) .calculate_formula_mass(
      x, carbon_offset = 6, hydrogen_intercept = 12,
      n_oxygen = 9, n_nitrogen = 1, n_phosphorus = 1
    ),
    PS = function(x) .calculate_formula_mass(
      x, carbon_offset = 6, hydrogen_intercept = 10,
      n_oxygen = 10, n_nitrogen = 1, n_phosphorus = 1
    ),
    LPC = function(x) .calculate_formula_mass(
      x, carbon_offset = 8, hydrogen_intercept = 18,
      n_oxygen = 7, n_nitrogen = 1, n_phosphorus = 1
    ),
    `LPC-P` = function(x) .calculate_formula_mass(
      x, carbon_offset = 8, hydrogen_intercept = 18,
      n_oxygen = 6, n_nitrogen = 1, n_phosphorus = 1
    ),
    `LPC-O` = function(x) .calculate_formula_mass(
      x, carbon_offset = 8, hydrogen_intercept = 20,
      n_oxygen = 6, n_nitrogen = 1, n_phosphorus = 1
    ),
    PC = function(x) .calculate_formula_mass(
      x, carbon_offset = 8, hydrogen_intercept = 16,
      n_oxygen = 8, n_nitrogen = 1, n_phosphorus = 1
    ),
    `PC-P` = function(x) .calculate_formula_mass(
      x, carbon_offset = 8, hydrogen_intercept = 16,
      n_oxygen = 7, n_nitrogen = 1, n_phosphorus = 1
    ),
    `PC-O` = function(x) .calculate_formula_mass(
      x, carbon_offset = 8, hydrogen_intercept = 18,
      n_oxygen = 7, n_nitrogen = 1, n_phosphorus = 1
    ),
    MG = function(x) .calculate_formula_mass(
      x, carbon_offset = 3, hydrogen_intercept = 6,
      n_oxygen = 4
    ),
    DG = function(x) .calculate_formula_mass(
      x, carbon_offset = 3, hydrogen_intercept = 4,
      n_oxygen = 5
    ),
    TG = function(x) .calculate_formula_mass(
      x, carbon_offset = 3, hydrogen_intercept = 2,
      n_oxygen = 6
    ),
    MGDG = function(x) .calculate_formula_mass(
      x, carbon_offset = 6, hydrogen_intercept = 10,
      n_oxygen = 10
    ),
    DGDG = function(x) .calculate_formula_mass(
      x, carbon_offset = 12, hydrogen_intercept = 20,
      n_oxygen = 15
    ),
    BMP = function(x) .calculate_formula_mass(
      x, carbon_offset = 6, hydrogen_intercept = 8,
      n_oxygen = 10, n_phosphorus = 1
    ),
    CL = function(x) .calculate_formula_mass(
      x, carbon_offset = 9, hydrogen_intercept = 14,
      n_oxygen = 17, n_phosphorus = 2
    ),
    NAGly = function(x) .calculate_formula_mass(
      x, carbon_offset = 2, hydrogen_intercept = 3,
      n_oxygen = 3, n_nitrogen = 1
    ),
    NAGlySer = function(x) .calculate_formula_mass(
      x, carbon_offset = 3, hydrogen_intercept = 5,
      n_oxygen = 4, n_nitrogen = 1
    ),
    NAOrn = function(x) .calculate_formula_mass(
      x, carbon_offset = 5, hydrogen_intercept = 10,
      n_oxygen = 2, n_nitrogen = 2
    ),
    NAE = function(x) .calculate_formula_mass(
      x, carbon_offset = 2, hydrogen_intercept = 5,
      n_oxygen = 2, n_nitrogen = 1
    ),
    LNAPE = function(x) .calculate_formula_mass(
      x, carbon_offset = 7, hydrogen_intercept = 13,
      n_oxygen = 8, n_nitrogen = 1, n_phosphorus = 1
    )
  )
}

#' Get a mass model
#'
#' @param name `character(1)`. Mass model name.
#' @param registry Named list of mass model functions.
#'
#' @return A mass model function.
#'
#' @noRd
#'
#' @keywords internal
.get_mass_model <- function(name, registry = lipidMassModels()) {
  fun <- registry[[name]]
  if (is.null(fun)) {
    stop("Unknown mass model: '", name, "'.", call. = FALSE)
  }
  fun
}

# ------------------------------------------------------------------------------
# Chain pools
# ------------------------------------------------------------------------------

#' Return default chain pools
#'
#' @examples
#' pools <- defaultChainPools()
#' names(pools)
#'
#' head(pools$fa)
#'
#' @return Named list of default chain pools.
#'
#' @export
defaultChainPools <- function() {
  list(
    fa = c(
      "8:0", "10:0", "12:0", "14:0", "14:1", "15:0", "16:0", "16:1", "16:2",
      "17:0", "17:1", "18:0", "18:1", "18:2", "18:3", "18:4", "19:0", "19:1",
      "20:0", "20:1", "20:2", "20:3", "20:4", "20:5", "22:0", "22:1", "22:2",
      "22:3", "22:4", "22:5", "22:6", "24:0", "24:1", "24:2", "24:3", "24:4",
      "24:5", "24:6", "26:0", "26:1", "26:2", "26:3", "26:4",  "26:5", "26:6"
    ),
    sph = c("16:0", "16:1", "18:0", "18:1", "18:2", "20:1", "24:1"),
    phyto = c("18:0"),
    alkyl = c("16:0", "16:1", "18:0", "18:1"),
    vinyl = c("16:0","16:1", "18:0", "18:1"),
    ox_fa = character(0)
  )
}

# ------------------------------------------------------------------------------
# Class specifications
# ------------------------------------------------------------------------------

#' Return default lipid class specifications
#'
#' @examples
#' specs <- defaultLipidClassSpecs()
#'
#' head(specs)
#' unique(specs$class)
#'
#' @return `data.frame` describing how each lipid class database is generated.
#'
#' @export
defaultLipidClassSpecs <- function() {
  classes <- c(
    "Cer_d", "Cer_t", "CerP_d", "CerP_t", "HexCer_d", "HexCer_t",
    "AHexCer_d", "AHexCer_t", "SHexCer_d", "SHexCer_t", "SM_d", "SM_t",
    "AcylCer_d", "AcylCer_t",
    "FA", "HFA", "Carnitine", "CE", "Sph_d", "Sph_t", "SphP_d", "SphP_t",
    "LPA", "LPA-O", "LPE", "LPE-O", "LPE-P", "LPG", "LPI", "LPS",
    "LPC", "LPC-O", "LPC-P", "MG",
    "DG", "PC", "PC-O", "PC-P", "PE", "PE-O", "PE-P", "PG", "PI",
    "PS", "PA", "FAHFA", "MGDG", "DGDG", "BMP", "LNAPE",
    "TG", "CL", "NAGly", "NAGlySer", "NAOrn", "NAE"
  )

  specs <- data.frame(class = classes, stringsAsFactors = FALSE)

  specs$backbone <- NA_character_
  specs$chain_pool_sequence <- NA_character_
  specs$mass_model <- specs$class
  specs$headgroup_code <- ""
  specs$headgroup_name <- "none"
  specs$headgroup_count <- 0L
  specs$headgroup_category <- "neutral"
  specs$headgroup_isomer_group <- ""
  specs$adducts_pos <- ""
  specs$adducts_neg <- ""
  specs$exclude_if_unique_pos <- ""
  specs$exclude_if_unique_neg <- ""

  sph_d <- c(
    "Cer_d", "CerP_d", "HexCer_d", "AHexCer_d", "SHexCer_d",
    "SM_d", "AcylCer_d"
  )
  sph_t <- c(
    "Cer_t", "CerP_t", "HexCer_t", "AHexCer_t", "SHexCer_t",
    "SM_t", "AcylCer_t"
  )

  specs$backbone[
    specs$class %in% c(sph_d, sph_t, "Sph_d", "Sph_t", "SphP_d", "SphP_t")
  ] <- "sphingoid"

  specs$chain_pool_sequence[specs$class %in% c(
    "Cer_d", "CerP_d", "HexCer_d", "AHexCer_d", "SHexCer_d", "SM_d"
  )] <- "sph|fa"

  specs$chain_pool_sequence[specs$class %in% c(
    "Cer_t", "CerP_t", "HexCer_t", "AHexCer_t", "SHexCer_t", "SM_t"
  )] <- "phyto|fa"

  specs$chain_pool_sequence[specs$class == "AcylCer_d"] <- "fa|sph|fa"
  specs$chain_pool_sequence[specs$class == "AcylCer_t"] <- "fa|phyto|fa"
  specs$chain_pool_sequence[specs$class == "Sph_d"] <- "sph"
  specs$chain_pool_sequence[specs$class == "Sph_t"] <- "phyto"
  specs$chain_pool_sequence[specs$class == "SphP_d"] <- "sph"
  specs$chain_pool_sequence[specs$class == "SphP_t"] <- "phyto"

  specs$headgroup_code[specs$class %in% c("CerP_d", "CerP_t", "SphP_d", "SphP_t")] <- "P"
  specs$headgroup_name[specs$class %in% c("CerP_d", "CerP_t", "SphP_d", "SphP_t")] <- "phosphate"
  specs$headgroup_count[specs$class %in% c("CerP_d", "CerP_t", "SphP_d", "SphP_t")] <- 1L
  specs$headgroup_category[specs$class %in% c("CerP_d", "CerP_t", "SphP_d", "SphP_t")] <- "phospho"

  specs$headgroup_code[specs$class %in% c("HexCer_d", "HexCer_t")] <- "Hex"
  specs$headgroup_name[specs$class %in% c("HexCer_d", "HexCer_t")] <- "hexose"
  specs$headgroup_count[specs$class %in% c("HexCer_d", "HexCer_t")] <- 1L
  specs$headgroup_category[specs$class %in% c("HexCer_d", "HexCer_t")] <- "glyco"
  specs$headgroup_isomer_group[specs$class %in% c("HexCer_d", "HexCer_t")] <- "GalCer/GlcCer"

  specs$headgroup_code[specs$class %in% c("AHexCer_d", "AHexCer_t")] <- "AHex"
  specs$headgroup_name[specs$class %in% c("AHexCer_d", "AHexCer_t")] <- "aminohexose"
  specs$headgroup_count[specs$class %in% c("AHexCer_d", "AHexCer_t")] <- 1L
  specs$headgroup_category[specs$class %in% c("AHexCer_d", "AHexCer_t")] <- "glyco"

  specs$headgroup_code[specs$class %in% c("SHexCer_d", "SHexCer_t")] <- "SHex"
  specs$headgroup_name[specs$class %in% c("SHexCer_d", "SHexCer_t")] <- "sulfohexose"
  specs$headgroup_count[specs$class %in% c("SHexCer_d", "SHexCer_t")] <- 1L
  specs$headgroup_category[specs$class %in% c("SHexCer_d", "SHexCer_t")] <- "glyco"

  specs$headgroup_code[specs$class %in% c("SM_d", "SM_t")] <- "PC"
  specs$headgroup_name[specs$class %in% c("SM_d", "SM_t")] <- "phosphocholine"
  specs$headgroup_count[specs$class %in% c("SM_d", "SM_t")] <- 1L
  specs$headgroup_category[specs$class %in% c("SM_d", "SM_t")] <- "phospho"

  glycerol_classes <- c(
    "LPA", "LPA-O", "LPE", "LPE-O", "LPE-P", "LPG", "LPI", "LPS",
    "LPC", "LPC-O", "LPC-P", "MG", "DG", "PC", "PC-O", "PC-P",
    "PE", "PE-O", "PE-P", "PG", "PI", "PS", "PA",
    "MGDG", "DGDG", "BMP", "LNAPE", "TG", "CL"
  )
  specs$backbone[specs$class %in% glycerol_classes] <- "glycerol"

  specs$chain_pool_sequence[specs$class %in% c(
    "LPA", "LPE", "LPG", "LPI", "LPS", "LPC", "MG"
  )] <- "fa"

  specs$chain_pool_sequence[specs$class %in% c("LPA-O", "LPE-O", "LPC-O")] <- "alkyl"
  specs$chain_pool_sequence[specs$class %in% c("LPE-P", "LPC-P")] <- "vinyl"

  specs$chain_pool_sequence[specs$class %in% c(
    "DG", "PA", "PC", "PE", "PG", "PI", "PS", "FAHFA",
    "MGDG", "DGDG", "BMP", "LNAPE"
  )] <- "fa|fa"

  specs$chain_pool_sequence[specs$class %in% c("PC-O", "PE-O")] <- "alkyl|fa"
  specs$chain_pool_sequence[specs$class %in% c("PC-P", "PE-P")] <- "vinyl|fa"
  specs$chain_pool_sequence[specs$class == "TG"] <- "fa|fa|fa"
  specs$chain_pool_sequence[specs$class == "CL"] <- "fa|fa|fa|fa"

  specs$headgroup_code[specs$class %in% c("LPA", "LPA-O", "PA")] <- "PA"
  specs$headgroup_name[specs$class %in% c("LPA", "LPA-O", "PA")] <- "phosphate"
  specs$headgroup_count[specs$class %in% c("LPA", "LPA-O", "PA")] <- 1L
  specs$headgroup_category[specs$class %in% c("LPA", "LPA-O", "PA")] <- "phospho"

  specs$headgroup_code[specs$class %in% c("LPE", "LPE-O", "LPE-P", "PE", "PE-O", "PE-P")] <- "PE"
  specs$headgroup_name[specs$class %in% c("LPE", "LPE-O", "LPE-P", "PE", "PE-O", "PE-P")] <- "phosphoethanolamine"
  specs$headgroup_count[specs$class %in% c("LPE", "LPE-O", "LPE-P", "PE", "PE-O", "PE-P")] <- 1L
  specs$headgroup_category[specs$class %in% c("LPE", "LPE-O", "LPE-P", "PE", "PE-O", "PE-P")] <- "phospho"

  specs$headgroup_code[specs$class %in% c("LPG", "PG")] <- "PG"
  specs$headgroup_name[specs$class %in% c("LPG", "PG")] <- "phosphoglycerol"
  specs$headgroup_count[specs$class %in% c("LPG", "PG")] <- 1L
  specs$headgroup_category[specs$class %in% c("LPG", "PG")] <- "phospho"

  specs$headgroup_code[specs$class %in% c("LPI", "PI")] <- "PI"
  specs$headgroup_name[specs$class %in% c("LPI", "PI")] <- "phosphoinositol"
  specs$headgroup_count[specs$class %in% c("LPI", "PI")] <- 1L
  specs$headgroup_category[specs$class %in% c("LPI", "PI")] <- "phospho"

  specs$headgroup_code[specs$class %in% c("LPS", "PS")] <- "PS"
  specs$headgroup_name[specs$class %in% c("LPS", "PS")] <- "phosphoserine"
  specs$headgroup_count[specs$class %in% c("LPS", "PS")] <- 1L
  specs$headgroup_category[specs$class %in% c("LPS", "PS")] <- "phospho"

  specs$headgroup_code[specs$class %in% c("LPC", "LPC-O", "LPC-P", "PC", "PC-O", "PC-P")] <- "PC"
  specs$headgroup_name[specs$class %in% c("LPC", "LPC-O", "LPC-P", "PC", "PC-O", "PC-P")] <- "phosphocholine"
  specs$headgroup_count[specs$class %in% c("LPC", "LPC-O", "LPC-P", "PC", "PC-O", "PC-P")] <- 1L
  specs$headgroup_category[specs$class %in% c("LPC", "LPC-O", "LPC-P", "PC", "PC-O", "PC-P")] <- "phospho"

  specs$headgroup_code[specs$class == "MGDG"] <- "Hex"
  specs$headgroup_name[specs$class == "MGDG"] <- "hexose"
  specs$headgroup_count[specs$class == "MGDG"] <- 1L
  specs$headgroup_category[specs$class == "MGDG"] <- "glyco"

  specs$headgroup_code[specs$class == "DGDG"] <- "Hex2"
  specs$headgroup_name[specs$class == "DGDG"] <- "dihexose"
  specs$headgroup_count[specs$class == "DGDG"] <- 2L
  specs$headgroup_category[specs$class == "DGDG"] <- "glyco"

  specs$headgroup_code[specs$class == "BMP"] <- "BMP"
  specs$headgroup_name[specs$class == "BMP"] <- "bis(monoacylglycero)phosphate"
  specs$headgroup_count[specs$class == "BMP"] <- 1L
  specs$headgroup_category[specs$class == "BMP"] <- "phospho"

  specs$headgroup_code[specs$class == "LNAPE"] <- "LNAPE"
  specs$headgroup_name[specs$class == "LNAPE"] <- "N-acyl phosphoethanolamine"
  specs$headgroup_count[specs$class == "LNAPE"] <- 1L
  specs$headgroup_category[specs$class == "LNAPE"] <- "phospho"

  specs$backbone[specs$class %in% c("FA", "HFA", "Carnitine", "FAHFA", "NAGly", "NAGlySer", "NAOrn", "NAE")] <- "fatty_acid"
  specs$backbone[specs$class == "CE"] <- "sterol"
  specs$chain_pool_sequence[specs$class %in% c("FA", "HFA", "Carnitine", "CE", "NAGly", "NAGlySer", "NAOrn", "NAE")] <- "fa"

  specs$headgroup_code[specs$class == "HFA"] <- "OH"
  specs$headgroup_name[specs$class == "HFA"] <- "hydroxy"
  specs$headgroup_count[specs$class == "HFA"] <- 1L
  specs$headgroup_category[specs$class == "HFA"] <- "modified_fatty_acid"

  specs$headgroup_code[specs$class == "Carnitine"] <- "Carnitine"
  specs$headgroup_name[specs$class == "Carnitine"] <- "carnitine"
  specs$headgroup_count[specs$class == "Carnitine"] <- 1L
  specs$headgroup_category[specs$class == "Carnitine"] <- "fatty_acyl_derivative"

  specs$headgroup_code[specs$class == "CE"] <- "Sterol"
  specs$headgroup_name[specs$class == "CE"] <- "cholesteryl"
  specs$headgroup_count[specs$class == "CE"] <- 1L
  specs$headgroup_category[specs$class == "CE"] <- "sterol_ester"

  specs$headgroup_code[specs$class == "NAGly"] <- "Gly"
  specs$headgroup_name[specs$class == "NAGly"] <- "glycine"
  specs$headgroup_count[specs$class == "NAGly"] <- 1L
  specs$headgroup_category[specs$class == "NAGly"] <- "amide"

  specs$headgroup_code[specs$class == "NAGlySer"] <- "GlySer"
  specs$headgroup_name[specs$class == "NAGlySer"] <- "glycylserine"
  specs$headgroup_count[specs$class == "NAGlySer"] <- 1L
  specs$headgroup_category[specs$class == "NAGlySer"] <- "amide"

  specs$headgroup_code[specs$class == "NAOrn"] <- "Orn"
  specs$headgroup_name[specs$class == "NAOrn"] <- "ornithine"
  specs$headgroup_count[specs$class == "NAOrn"] <- 1L
  specs$headgroup_category[specs$class == "NAOrn"] <- "amide"

  specs$headgroup_code[specs$class == "NAE"] <- "Ethanolamide"
  specs$headgroup_name[specs$class == "NAE"] <- "ethanolamide"
  specs$headgroup_count[specs$class == "NAE"] <- 1L
  specs$headgroup_category[specs$class == "NAE"] <- "amide"

  specs$adducts_pos[specs$class %in% c("Cer_d", "Cer_t")] <- "M+H-H2O|M+Na|M+H"
  specs$adducts_neg[specs$class %in% c("Cer_d", "Cer_t")] <- "M+CH3COO|M+CH3COO-CH3|M-H"

  specs$adducts_pos[specs$class %in% c("CerP_d", "CerP_t")] <- "M+H"
  specs$adducts_neg[specs$class %in% c("CerP_d", "CerP_t")] <- "M+CH3COO|M+CH3COO-CH3|M-H"

  specs$adducts_pos[specs$class %in% c("HexCer_d", "HexCer_t", "AHexCer_d", "AHexCer_t", "SHexCer_d", "SHexCer_t")] <- "M+H|M+Na"
  specs$adducts_neg[specs$class %in% c("HexCer_d", "HexCer_t", "AHexCer_d", "AHexCer_t")] <- "M+CH3COO|M+CH3COO-CH3|M-H"
  specs$adducts_neg[specs$class %in% c("SHexCer_d", "SHexCer_t")] <- "M-H"

  specs$adducts_pos[specs$class %in% c("SM_d", "SM_t")] <- "M+H|M+Na"
  specs$adducts_neg[specs$class %in% c("SM_d", "SM_t")] <- "M+CH3COO|M+CH3COO-CH3|M-H|M-CH3"

  specs$adducts_pos[specs$class %in% c("AcylCer_d", "AcylCer_t")] <- "M+H|M+H-H2O|M+Na"
  specs$adducts_neg[specs$class %in% c("AcylCer_d", "AcylCer_t")] <- "M+CH3COO|M+CH3COO-CH3|M-H|M-CH3"

  specs$adducts_neg[specs$class == "FA"] <- "M-H|2M-H"
  specs$adducts_neg[specs$class == "HFA"] <- "M-H"
  specs$adducts_pos[specs$class == "Carnitine"] <- "M+H|M+Na"
  specs$adducts_pos[specs$class == "CE"] <- "2M+NH4|2M+Na|M+NH4|M+Na"
  specs$adducts_pos[specs$class %in% c("Sph_d", "Sph_t")] <- "M+H"
  specs$adducts_neg[specs$class %in% c("Sph_d", "Sph_t")] <- "M-H"
  specs$adducts_pos[specs$class %in% c("SphP_d", "SphP_t")] <- "M+H"
  specs$adducts_neg[specs$class %in% c("SphP_d", "SphP_t")] <- "M-H"

  specs$adducts_neg[specs$class %in% c("LPA", "LPA-O", "LPG", "LPI", "PA", "FAHFA", "MGDG", "DGDG", "BMP", "LNAPE")] <- "M-H"
  specs$adducts_pos[specs$class %in% c("LPE", "LPE-O", "LPE-P", "LPC", "LPC-O", "LPC-P")] <- "M+H|M+Na"
  specs$adducts_neg[specs$class == "LPE"] <- "M-H"
  specs$adducts_neg[specs$class %in% c("LPE-O", "LPE-P")] <- "M-H|M+NaCH3COO"
  specs$adducts_neg[specs$class == "LPS"] <- "M-H|M+Na-2H"
  specs$adducts_neg[specs$class %in% c("LPC", "LPC-O", "LPC-P")] <- "M+CH3COO|M+CH3COO-CH3|M-CH3|M-H"

  specs$adducts_pos[specs$class == "MG"] <- "M+H-H2O|M+NH4|M+Na"
  specs$adducts_pos[specs$class == "DG"] <- "M+H-H2O|M+NH4|M+Na"

  specs$adducts_pos[specs$class %in% c("PC", "PC-O", "PC-P")] <- "M+H|M+Na"
  specs$adducts_neg[specs$class %in% c("PC", "PC-O", "PC-P")] <- "M+CH3COO|M+CH3COO-CH3|M-CH3|M-H"

  specs$adducts_pos[specs$class %in% c("PE", "PE-O", "PE-P")] <- "M+H|M+Na"
  specs$adducts_neg[specs$class == "PE"] <- "M-H"
  specs$adducts_neg[specs$class %in% c("PE-O", "PE-P")] <- "M-H|M+NaCH3COO"

  specs$adducts_pos[specs$class %in% c("PG", "PI")] <- "M+H|M+NH4|M+Na"
  specs$adducts_neg[specs$class %in% c("PG", "PI")] <- "M-H"

  specs$adducts_neg[specs$class == "PS"] <- "M-H|M+Na-2H"

  specs$adducts_pos[specs$class %in% c("MGDG", "DGDG", "BMP", "TG", "CL")] <- "M+NH4|M+Na"
  specs$adducts_neg[specs$class == "CL"] <- "M-H|M+Na-2H"

  specs$adducts_pos[specs$class %in% c("NAGly", "NAGlySer", "NAOrn", "NAE")] <- "M+H|M+Na"
  specs$adducts_neg[specs$class %in% c("NAGly", "NAGlySer", "NAOrn", "NAE")] <- "M-H"

  specs$exclude_if_unique_neg[
    specs$class %in% c(
      "Cer_d", "Cer_t",
      "CerP_d", "CerP_t",
      "HexCer_d", "HexCer_t",
      "AHexCer_d", "AHexCer_t",
      "SM_d", "SM_t",
      "AcylCer_d", "AcylCer_t",
      "LPC", "LPC-O", "LPC-P",
      "PC", "PC-O", "PC-P"
    )
  ] <- "M-CH3|M+CH3COO-CH3"

  specs$exclude_if_unique_neg[
    specs$class %in% c(
      "PS"
    )
  ] <- "M+Na-2H"

  specs$exclude_if_unique_pos[
    specs$class %in% c(
      "MG", "DG", "TG",
      "Cer_d", "Cer_t",
      "AcylCer_d", "AcylCer_t"
    )
  ] <- "M+H-H2O"

  stopifnot(nrow(specs) == length(classes))
  specs
}

# ------------------------------------------------------------------------------
# Validation
# ------------------------------------------------------------------------------

#' Validate lipid class specifications
#'
#' @param specs `data.frame` with class generation specifications.
#' @param chain_pools Named list of chain pools.
#' @param mass_registry Named list of mass model functions.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @noRd
#'
#' @keywords internal
.validate_lipid_class_specs <- function(specs, chain_pools,
                                        mass_registry = lipidMassModels()) {
  required_cols <- c(
    "class", "backbone", "chain_pool_sequence", "mass_model",
    "headgroup_code", "headgroup_name", "headgroup_count",
    "headgroup_category", "headgroup_isomer_group",
    "adducts_pos", "adducts_neg"
  )

  missing_cols <- setdiff(required_cols, colnames(specs))
  if (length(missing_cols)) {
    stop(
      "Missing columns in 'specs': ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyDuplicated(specs$class)) {
    stop("Lipid class names in 'specs' must be unique.", call. = FALSE)
  }

  if (!all(specs$mass_model %in% names(mass_registry))) {
    stop("Unknown mass models found in 'specs'.", call. = FALSE)
  }

  if (anyNA(specs$chain_pool_sequence) ||
      any(!nzchar(specs$chain_pool_sequence))) {
    stop("'chain_pool_sequence' must be defined for all classes.", call. = FALSE)
  }

  if (anyNA(specs$backbone) || any(!nzchar(specs$backbone))) {
    stop("'backbone' must be defined for all classes.", call. = FALSE)
  }

  if (anyNA(specs$headgroup_count) || any(specs$headgroup_count < 0)) {
    stop("'headgroup_count' must contain non-negative integers.", call. = FALSE)
  }

  pool_names <- names(chain_pools)
  spec_pools <- unique(unlist(strsplit(
    specs$chain_pool_sequence, "\\|", perl = TRUE
  )))

  if (!all(spec_pools %in% pool_names)) {
    missing_pools <- setdiff(spec_pools, pool_names)
    stop(
      "Unknown chain pools found in 'chain_pool_sequence': ",
      paste(missing_pools, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# Combination generators
# ------------------------------------------------------------------------------

#' Return unique rows from a numeric sum composition matrix
#'
#' @param x Numeric matrix with columns `n_carbon`, `n_double_bonds` and
#' `n_oxygen`.
#'
#' @return Numeric matrix with unique rows.
#'
#' @noRd
#'
#' @keywords internal
.unique_sum_composition_matrix <- function(x) {
  if (!is.matrix(x) || ncol(x) != 3L) {
    stop("'x' must be a numeric matrix with 3 columns.", call. = FALSE)
  }

  x[!duplicated(x), , drop = FALSE]
}

#' Parse chain pool sequence from class specification
#'
#' @param spec One-row `data.frame` with class specification.
#'
#' @return Character vector with chain pool names in sequence order.
#'
#' @noRd
#'
#' @keywords internal
.get_chain_pool_sequence <- function(spec) {
  chain_pool_sequence <- spec$chain_pool_sequence[[1L]]

  if (!is.character(chain_pool_sequence) ||
      length(chain_pool_sequence) != 1L ||
      is.na(chain_pool_sequence) ||
      !nzchar(chain_pool_sequence)) {
    stop(
      "'chain_pool_sequence' must be a non-empty character(1).",
      call. = FALSE
    )
  }

  strsplit(chain_pool_sequence, "\\|", perl = TRUE)[[1L]]
}

#' Generate sum compositions from multiple chain pools
#'
#' @param pools List of numeric matrices with columns `n_carbon`,
#' `n_double_bonds` and `n_oxygen`.
#'
#' @return Character vector of unique sum compositions.
#'
#' @noRd
#'
#' @keywords internal
.generate_multi_pool_sum_compositions <- function(pools) {
  if (!is.list(pools) || !length(pools)) {
    stop("'pools' must be a non-empty list.", call. = FALSE)
  }

  if (any(vapply(pools, nrow, integer(1)) == 0L)) {
    return(character())
  }

  n_rows <- vapply(pools, nrow, integer(1))
  idx_grid <- expand.grid(
    lapply(n_rows, seq_len),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  idx_grid <- as.matrix(idx_grid)

  n_carbon <- numeric(nrow(idx_grid))
  n_double_bonds <- numeric(nrow(idx_grid))
  n_oxygen <- numeric(nrow(idx_grid))

  for (i in seq_along(pools)) {
    n_carbon <- n_carbon + pools[[i]][idx_grid[, i], "n_carbon"]
    n_double_bonds <- n_double_bonds +
      pools[[i]][idx_grid[, i], "n_double_bonds"]
    n_oxygen <- n_oxygen + pools[[i]][idx_grid[, i], "n_oxygen"]
  }

  sum_compositions <- cbind(
    n_carbon = n_carbon,
    n_double_bonds = n_double_bonds,
    n_oxygen = n_oxygen
  )
  sum_compositions <- .unique_sum_composition_matrix(sum_compositions)

  .sum_compositions_to_string(sum_compositions)
}

#' Generate sum compositions for one lipid class
#'
#' @param spec One-row `data.frame` with class specification.
#' @param chain_pools Named list of parsed chain pools as numeric matrices.
#'
#' @return Character vector of sum compositions.
#'
#' @noRd
#'
#' @keywords internal
.generate_sum_compositions_from_spec <- function(spec, chain_pools) {
  pool_names <- .get_chain_pool_sequence(spec)

  missing_pools <- setdiff(pool_names, names(chain_pools))
  if (length(missing_pools)) {
    stop(
      "Unknown chain pools in specification: ",
      paste(missing_pools, collapse = ", "),
      call. = FALSE
    )
  }

  selected_pools <- lapply(pool_names, function(z) chain_pools[[z]])
  .generate_multi_pool_sum_compositions(selected_pools)
}

# ------------------------------------------------------------------------------
# Database builders
# ------------------------------------------------------------------------------

#' Build a lipid database table from sum compositions
#'
#' @param sum_compositions Character vector of lipid sum compositions.
#' @param spec One-row `data.frame` with class specification.
#' @param mass_registry Named list of mass model functions.
#'
#' @return `data.frame`.
#'
#' @noRd
#'
#' @keywords internal
.build_lipid_db <- function(sum_compositions, spec,
                            mass_registry = lipidMassModels()) {
  lipid_class <- spec$class[[1L]]
  mass_model <- spec$mass_model[[1L]]
  fun <- .get_mass_model(mass_model, registry = mass_registry)
  sum_compositions <- unique(sum_compositions)

  if (!length(sum_compositions)) {
    return(data.frame(
      class = character(),
      sum_composition = character(),
      n_carbon = numeric(),
      n_double_bonds = numeric(),
      n_oxygen = numeric(),
      molecular_formula = character(),
      neutral_exact_mass = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  sum_df <- .sum_compositions_to_dataframe(sum_compositions)
  res <- lapply(sum_compositions, fun)

  data.frame(
    class = lipid_class,
    sum_composition = sum_df$sum_composition,
    n_carbon = sum_df$n_carbon,
    n_double_bonds = sum_df$n_double_bonds,
    n_oxygen = sum_df$n_oxygen,
    molecular_formula = vapply(res, `[[`, character(1), "molecular_formula"),
    neutral_exact_mass = vapply(res, `[[`, numeric(1), "neutral_exact_mass"),
    stringsAsFactors = FALSE
  )
}

#' Build one lipid class database from its specification
#'
#' @param spec One-row `data.frame`.
#' @param chain_pools Named list of parsed chain pools as numeric matrices.
#' @param mass_registry Named list of mass model functions.
#'
#' @return `data.frame`.
#'
#' @noRd
#'
#' @keywords internal
.build_lipid_db_from_spec <- function(spec, chain_pools,
                                      mass_registry = lipidMassModels()) {
  sum_compositions <- .generate_sum_compositions_from_spec(spec, chain_pools)
  .build_lipid_db(
    sum_compositions = sum_compositions,
    spec = spec,
    mass_registry = mass_registry
  )
}

#' Build adduct table from lipid class specifications
#'
#' @param class_specs `data.frame` with lipid class specifications.
#'
#' @return `data.frame` with columns `class`, `polarity`, `adduct`,
#' `adduct_priority` and `exclude_if_unique`.
#'
#' @noRd
#'
#' @keywords internal
.build_adduct_table <- function(class_specs) {
  out <- vector("list", nrow(class_specs) * 2L)
  k <- 1L

  for (i in seq_len(nrow(class_specs))) {
    class_i <- class_specs$class[[i]]
    adducts_pos <- class_specs$adducts_pos[[i]]
    adducts_neg <- class_specs$adducts_neg[[i]]

    exclude_pos <- if ("exclude_if_unique_pos" %in% colnames(class_specs)) {
      class_specs$exclude_if_unique_pos[[i]]
    } else {
      ""
    }

    exclude_neg <- if ("exclude_if_unique_neg" %in% colnames(class_specs)) {
      class_specs$exclude_if_unique_neg[[i]]
    } else {
      ""
    }

    if (!is.na(adducts_pos) && nzchar(adducts_pos)) {
      pos_vals <- strsplit(adducts_pos, "\\|", perl = TRUE)[[1L]]

      exclude_vals <- if (!is.na(exclude_pos) && nzchar(exclude_pos)) {
        strsplit(exclude_pos, "\\|", perl = TRUE)[[1L]]
      } else {
        character()
      }

      out[[k]] <- data.frame(
        class = class_i,
        polarity = "positive",
        adduct = pos_vals,
        adduct_priority = seq_along(pos_vals),
        exclude_if_unique = as.integer(pos_vals %in% exclude_vals),
        stringsAsFactors = FALSE
      )

      k <- k + 1L
    }

    if (!is.na(adducts_neg) && nzchar(adducts_neg)) {
      neg_vals <- strsplit(adducts_neg, "\\|", perl = TRUE)[[1L]]

      exclude_vals <- if (!is.na(exclude_neg) && nzchar(exclude_neg)) {
        strsplit(exclude_neg, "\\|", perl = TRUE)[[1L]]
      } else {
        character()
      }

      out[[k]] <- data.frame(
        class = class_i,
        polarity = "negative",
        adduct = neg_vals,
        adduct_priority = seq_along(neg_vals),
        exclude_if_unique = as.integer(neg_vals %in% exclude_vals),
        stringsAsFactors = FALSE
      )

      k <- k + 1L
    }
  }

  if (k == 1L) {
    return(data.frame(
      class = character(),
      polarity = character(),
      adduct = character(),
      adduct_priority = integer(),
      exclude_if_unique = integer(),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, out[seq_len(k - 1L)])
  rownames(out) <- NULL
  out
}

#' Build class metadata table from lipid class specifications
#'
#' @param class_specs `data.frame` with lipid class specifications.
#'
#' @return `data.frame` with class-level metadata.
#'
#' @noRd
#'
#' @keywords internal
.build_class_metadata_table <- function(class_specs) {
  class_specs[, c(
    "class", "backbone", "chain_pool_sequence", "mass_model",
    "headgroup_code", "headgroup_name", "headgroup_count",
    "headgroup_category", "headgroup_isomer_group"
  ), drop = FALSE]
}

# ------------------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------------------

#' Create customizable lipid databases
#'
#' Create lipid annotation databases from lipid class specifications and chain
#' pools.
#'
#' @param lipid_classes Character vector with lipid classes to generate. Use
#' `NULL` to generate all classes defined in `class_specs`.
#' @param chain_pools Named list of chain pools. See [defaultChainPools()].
#' @param class_specs `data.frame` with lipid class specifications. See
#' [defaultLipidClassSpecs()].
#' @param mass_registry Named list of mass model functions.
#'
#' @examples
#' dbs <- createLipidDbs(
#'     lipid_classes = c("PC", "PE")
#' )
#'
#' names(dbs)
#'
#' @return A named list with elements `dbs`, `adducts` and `class_metadata`.
#' `adducts` contains columns `class`, `polarity`, `adduct` and
#' `adduct_priority`, where lower values indicate higher priority within each
#' lipid class and polarity.
#'
#' @export
createLipidDbs <- function(lipid_classes = NULL,
                           chain_pools = defaultChainPools(),
                           class_specs = defaultLipidClassSpecs(),
                           mass_registry = lipidMassModels()) {
  .validate_chain_pools(chain_pools)
  .validate_lipid_class_specs(
    class_specs,
    chain_pools = chain_pools,
    mass_registry = mass_registry
  )

  parsed_chain_pools <- .parse_chain_pools(chain_pools)

  if (is.null(lipid_classes)) {
    lipid_classes <- class_specs$class
  }

  if (!all(lipid_classes %in% class_specs$class)) {
    missing_classes <- setdiff(lipid_classes, class_specs$class)
    stop(
      "Unknown lipid classes requested: ",
      paste(missing_classes, collapse = ", "),
      call. = FALSE
    )
  }

  specs_sub <- class_specs[match(lipid_classes, class_specs$class), ,
                           drop = FALSE]

  dbs <- lapply(seq_len(nrow(specs_sub)), function(i) {
    .build_lipid_db_from_spec(
      spec = specs_sub[i, , drop = FALSE],
      chain_pools = parsed_chain_pools,
      mass_registry = mass_registry
    )
  })
  dbs <- stats::setNames(dbs, specs_sub$class)

  list(
    dbs = dbs,
    adducts = .build_adduct_table(specs_sub),
    class_metadata = .build_class_metadata_table(specs_sub)
  )
}

# ------------------------------------------------------------------------------
# Editable configuration helpers
# ------------------------------------------------------------------------------

#' Convert a character vector to a pipe-separated string
#'
#' @param x Character vector.
#'
#' @return `character(1)`.
#'
#' @noRd
#'
#' @keywords internal
.collapse_pipe <- function(x) {
  paste(x, collapse = "|")
}

#' Create a default editable lipid database configuration
#'
#' Create a template data frame with lipid class metadata, mass models, adducts
#' and chain pool sequences. The chain pool values are stored in a companion
#' table attached as the `pool_table` attribute.
#'
#' @param path Optional `character(1)` path to save the configuration file.
#' @param sep Field separator for text output. Ignored for `.csv` paths.
#' @param class_specs `data.frame` with lipid class specifications.
#' @param chain_pools Named list of chain pools.
#'
#' @examples
#' config <- defaultLipidDbConfig()
#'
#' head(config)
#' attr(config, "pool_table")
#'
#' @return `data.frame` with class-level configuration. The `pool_table`
#' attribute contains the editable chain pool definitions.
#'
#' @export
defaultLipidDbConfig <- function(path = NULL, sep = ",",
                                 class_specs = defaultLipidClassSpecs(),
                                 chain_pools = defaultChainPools()) {
  .validate_chain_pools(chain_pools)
  .validate_lipid_class_specs(
    class_specs,
    chain_pools = chain_pools
  )

  cfg <- class_specs[, c(
    "class", "backbone", "chain_pool_sequence", "mass_model",
    "headgroup_code", "headgroup_name", "headgroup_count",
    "headgroup_category", "headgroup_isomer_group",
    "adducts_pos", "adducts_neg"
  ), drop = FALSE]

  pool_names <- unique(unlist(strsplit(
    cfg$chain_pool_sequence,
    "\\|",
    perl = TRUE
  )))

  pool_values <- vapply(
    pool_names,
    function(z) .collapse_pipe(chain_pools[[z]]),
    character(1)
  )

  pool_table <- data.frame(
    pool_name = pool_names,
    pool_values = pool_values,
    stringsAsFactors = FALSE
  )

  if (!is.null(path)) {
    if (!is.character(path) || length(path) != 1L) {
      stop("'path' must be a character(1) or NULL.", call. = FALSE)
    }

    if (grepl("\\.csv$", path, ignore.case = TRUE)) {
      utils::write.csv(cfg, file = path, row.names = FALSE)

      pool_path <- sub("\\.csv$", "_pools.csv", path, ignore.case = TRUE)
      utils::write.csv(pool_table, file = pool_path, row.names = FALSE)
    } else {
      utils::write.table(
        cfg,
        file = path,
        sep = sep,
        quote = TRUE,
        row.names = FALSE,
        col.names = TRUE
      )
    }
  }

  attr(cfg, "pool_table") <- pool_table
  cfg
}

#' Build chain pools from editable configuration
#'
#' @param config `data.frame` returned by [defaultLipidDbConfig()] or edited by
#' the user.
#' @param pool_table Optional `data.frame` with columns `pool_name` and
#' `pool_values`. If `NULL`, the function tries to read it from the `pool_table`
#' attribute of `config`.
#'
#' @examples
#' config <- defaultLipidDbConfig()
#'
#' pools <- chainPoolsFromConfig(config)
#'
#' names(pools)
#' head(pools$fa)
#'
#' @return Named list of chain pools.
#'
#' @export
chainPoolsFromConfig <- function(config, pool_table = NULL) {
  if (is.null(pool_table)) {
    pool_table <- attr(config, "pool_table")
  }

  if (is.null(pool_table) || !is.data.frame(pool_table)) {
    stop(
      "'pool_table' must be provided or stored as an attribute of 'config'.",
      call. = FALSE
    )
  }

  required <- c("pool_name", "pool_values")
  if (!all(required %in% colnames(pool_table))) {
    stop(
      "'pool_table' must contain columns 'pool_name' and 'pool_values'.",
      call. = FALSE
    )
  }

  pool_map <- stats::setNames(
    vector("list", nrow(pool_table)),
    pool_table$pool_name
  )

  for (i in seq_len(nrow(pool_table))) {
    pool_name <- pool_table$pool_name[[i]]
    pool_values <- pool_table$pool_values[[i]]

    if (is.na(pool_values) || !nzchar(pool_values)) {
      pool_map[[pool_name]] <- character()
      next
    }

    pool_map[[pool_name]] <- strsplit(
      pool_values,
      "\\|",
      perl = TRUE
    )[[1L]]
  }

  .validate_chain_pools(pool_map)
  pool_map
}

#' Build lipid class specifications from editable configuration
#'
#' @param config `data.frame` returned by [defaultLipidDbConfig()] or edited by
#' the user.
#'
#' @examples
#' config <- defaultLipidDbConfig()
#'
#' specs <- classSpecsFromConfig(config)
#'
#' head(specs)
#'
#' @return `data.frame` with lipid class specifications.
#'
#' @export
classSpecsFromConfig <- function(config) {
  required <- c(
    "class", "backbone", "chain_pool_sequence", "mass_model",
    "headgroup_code", "headgroup_name", "headgroup_count",
    "headgroup_category", "headgroup_isomer_group",
    "adducts_pos", "adducts_neg"
  )

  if (!all(required %in% colnames(config))) {
    missing_cols <- setdiff(required, colnames(config))
    stop(
      "Invalid configuration table. Missing columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  config[, required, drop = FALSE]
}

#' Create lipid databases from an editable configuration
#'
#' @param config `data.frame` editable lipid class configuration.
#' @param lipid_classes Optional character vector with lipid classes to generate.
#' @param pool_table Optional `data.frame` with chain pool definitions.
#' @param mass_registry Named list of mass model functions.
#'
#' @examples
#' config <- defaultLipidDbConfig()
#'
#' dbs <- createLipidDbsFromConfig(
#'     config,
#'     lipid_classes = c("PC", "PE")
#' )
#'
#' names(dbs)
#'
#' @return A named list with elements `dbs`, `adducts` and `class_metadata`.
#'
#' @export
createLipidDbsFromConfig <- function(config, lipid_classes = NULL,
                                     pool_table = NULL,
                                     mass_registry = lipidMassModels()) {
  chain_pools <- chainPoolsFromConfig(
    config = config,
    pool_table = pool_table
  )

  class_specs <- classSpecsFromConfig(config)

  createLipidDbs(
    lipid_classes = lipid_classes,
    chain_pools = chain_pools,
    class_specs = class_specs,
    mass_registry = mass_registry
  )
}

#' Create LipidMS databases from an Excel configuration file
#'
#' Read lipid-class and chain-pool configuration tables from an Excel
#' workbook and use them to create the LipidMS internal databases.
#'
#' The workbook must contain:
#'
#' \itemize{
#'   \item A sheet named `"lipid_classes_config"` containing the lipid-class
#'   configuration.
#'   \item A sheet named `"chain_pools"` containing the chain-pool definitions.
#' }
#'
#' @param file `character(1)`. Path to an `.xlsx` or `.xls` configuration
#'   file.
#' @param lipid_classes Optional character vector specifying the lipid classes
#'   to include. Passed to [createLipidDbsFromConfig()].
#' @param mass_registry Named list of lipid mass models. By default,
#'   lipidMassModels().
#' @param config_sheet `character(1)`. Name of the worksheet containing the
#'   lipid-class configuration.
#' @param pool_sheet `character(1)`. Name of the worksheet containing the
#'   chain-pool definitions.
#'
#' @examples
#' \dontrun{
#' dbs <- createLipidDbsFromFile(
#'   "LipidMS_default_db_config.xlsx"
#' )
#'
#' dbs <- createLipidDbsFromFile(
#'   "LipidMS_default_db_config.xlsx",
#'   lipid_classes = c("PC", "PE", "TG")
#' )
#' }
#'
#' @return The object returned by [createLipidDbsFromConfig()].
#'
#' @details
#' This function requires the `readxl` package.
#'
#' @export
createLipidDbsFromFile <- function(
    file,
    lipid_classes = NULL,
    mass_registry = lipidMassModels(),
    config_sheet = "lipid_classes_config",
    pool_sheet = "chain_pools") {

  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop(
      "Package 'readxl' is required to read Excel configuration files.",
      call. = FALSE
    )
  }

  if (!is.character(file) ||
      length(file) != 1L ||
      is.na(file) ||
      !nzchar(file)) {
    stop(
      "'file' must be a non-empty character(1).",
      call. = FALSE
    )
  }

  if (!file.exists(file)) {
    stop(
      "Configuration file does not exist: ",
      file,
      call. = FALSE
    )
  }

  extension <- tolower(tools::file_ext(file))

  if (!extension %in% c("xlsx", "xls")) {
    stop(
      "'file' must be an Excel file with extension '.xlsx' or '.xls'.",
      call. = FALSE
    )
  }

  if (!is.character(config_sheet) ||
      length(config_sheet) != 1L ||
      is.na(config_sheet) ||
      !nzchar(config_sheet)) {
    stop(
      "'config_sheet' must be a non-empty character(1).",
      call. = FALSE
    )
  }

  if (!is.character(pool_sheet) ||
      length(pool_sheet) != 1L ||
      is.na(pool_sheet) ||
      !nzchar(pool_sheet)) {
    stop(
      "'pool_sheet' must be a non-empty character(1).",
      call. = FALSE
    )
  }

  available_sheets <- readxl::excel_sheets(file)

  missing_sheets <- setdiff(
    c(config_sheet, pool_sheet),
    available_sheets
  )

  if (length(missing_sheets)) {
    stop(
      "Missing required worksheet",
      if (length(missing_sheets) > 1L) "s" else "",
      ": ",
      paste(missing_sheets, collapse = ", "),
      ". Available worksheets: ",
      paste(available_sheets, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  config <- as.data.frame(
    readxl::read_excel(
      path = file,
      sheet = config_sheet,
      .name_repair = "check_unique"
    ),
    stringsAsFactors = FALSE
  )

  pool_table <- as.data.frame(
    readxl::read_excel(
      path = file,
      sheet = pool_sheet,
      .name_repair = "check_unique"
    ),
    stringsAsFactors = FALSE
  )

  if (!nrow(config)) {
    stop(
      "Worksheet '",
      config_sheet,
      "' contains no lipid-class definitions.",
      call. = FALSE
    )
  }

  if (!nrow(pool_table)) {
    stop(
      "Worksheet '",
      pool_sheet,
      "' contains no chain-pool definitions.",
      call. = FALSE
    )
  }

  createLipidDbsFromConfig(
    config = config,
    lipid_classes = lipid_classes,
    pool_table = pool_table,
    mass_registry = mass_registry
  )
}

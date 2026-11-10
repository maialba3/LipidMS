.step_time <- function(label, expr, verbose = TRUE) {
  if (isTRUE(verbose)) {
    cat("\n---- START annotateLipids::", label, " ----\n", sep = "")
  }
  t0 <- Sys.time()
  value <- force(expr)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (isTRUE(verbose)) {
    cat(
      "---- END annotateLipids::", label,
      " | elapsed_min = ", round(elapsed, 2), " ----\n",
      sep = ""
    )
  }
  value
}

#' Default lipid annotation parameters
#'
#' Return the default parameters used by [annotateLipids()] for MS1/MS2
#' lipid annotation.
#'
#' @param polarity `character(1)`. Ionization polarity. Must be either
#' `"positive"` or `"negative"`.
#'
#' @examples
#' param <- defaultAnnotationParameters("negative")
#' param
#'
#' param_pos <- defaultAnnotationParameters("positive")
#' param_pos
#'
#' @return Named list.
#'
#' @export
defaultAnnotationParameters <- function(polarity) {

  if (missing(polarity)) {
    stop(
      "'polarity' must be specified as either 'positive' or 'negative'.",
      call. = FALSE
    )
  }

  polarity <- match.arg(
    tolower(polarity),
    c("positive", "negative")
  )

  list(
    polarity = polarity,

    dmzprecursor = 5,
    dmzproducts = 10,
    rttol = 3,

    isotope_mzabs = 0.01,
    isotope_maxcharge = 2L,
    isotope_maxiso = 5L,
    isotope_minfrac = 0.5,
    isotope_rttol = 5,
    intval = "into",

    min_ms2_intensity = 100,
    remove_isotopes = TRUE,
    only_confirmed_monoisotopic = TRUE
  )
}

#' Default DIA reconstruction parameters
#'
#' Return the default parameters used by [annotateLipids()] for DIA/SWATH
#' chromatographic peak detection and spectrum reconstruction.
#'
#' @param findChromPeaksIsolationWindowParam
#' A `findChromPeaksIsolationWindowParam` object. If `NULL`, the user must
#' provide it before calling [annotateLipids()].
#'
#' @param min_ms2_intensity Minimum reconstructed fragment intensity.
#'
#' @examples
#' paramdia <- defaultDIAExtractionParameters()
#' paramdia
#'
#' @return Named list.
#'
#' @export
defaultDIAExtractionParameters <- function(
    findChromPeaksIsolationWindowParam = NULL,
    min_ms2_intensity = 100
) {
  list(
    findChromPeaksIsolationWindowParam =
      findChromPeaksIsolationWindowParam,
    expandRt = 0,
    minCor = 0.7,
    diffRt = 3,
    intensity = "into",
    minint = min_ms2_intensity,
    BPPARAM = BiocParallel::SerialParam()
  )
}

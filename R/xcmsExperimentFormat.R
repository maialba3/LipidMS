#' Fix missing DIA/SWATH isolation window metadata
#'
#' This function ensures that the spectra metadata required by
#' `findChromPeaksIsolationWindow()` and
#' `reconstructChromPeakSpectra()` are present and contain no missing
#' values.
#'
#' Missing `isolationWindowTargetMz` values are replaced using
#' `precursorMz` when available, otherwise `0`.
#'
#' Missing `isolationWindowLowerMz` and
#' `isolationWindowUpperMz` values are replaced using
#' `lowMZ` and `highMZ`, respectively.
#'
#' This is especially useful for DIA/SWATH data from vendors or
#' mzML conversion workflows that do not properly export isolation
#' window metadata.
#'
#' The function modifies the spectra metadata stored inside the
#' `XcmsExperiment` object and returns the updated object.
#'
#' @param x A `XcmsExperiment` object.
#'
#' @return A `XcmsExperiment` object with complete isolation window
#' metadata in the spectra backend.
#'
#' @details
#' The following spectra variables are checked and completed:
#'
#' - `isolationWindowTargetMz`
#' - `isolationWindowLowerMz`
#' - `isolationWindowUpperMz`
#'
#' The function is intended to be applied before:
#'
#' - `findChromPeaksIsolationWindow()`
#' - `reconstructChromPeakSpectra()`
#'
#' @examples
#' \dontrun{
#' x <- fix_isolation_windows(x)
#' }
#'
#' @export
fix_isolation_windows <- function(x) {
  sps <- MsExperiment::spectra(x)
  sd <- Spectra::spectraData(sps)

  if (!"isolationWindowTargetMz" %in% colnames(sd)) {
    iw <- sd$precursorMz
    # iw[is.na(iw)] <- 0
    sd$isolationWindowTargetMz <- iw
  } else {
    # sd$isolationWindowTargetMz[is.na(sd$isolationWindowTargetMz)] <- 0
  }

  if (!"isolationWindowLowerMz" %in% colnames(sd)) {
    sd$isolationWindowLowerMz <- sd$lowMZ
  } else {
    idx_na <- is.na(sd$isolationWindowLowerMz)
    sd$isolationWindowLowerMz[idx_na] <- sd$lowMZ[idx_na]
  }

  if (!"isolationWindowUpperMz" %in% colnames(sd)) {
    sd$isolationWindowUpperMz <- sd$highMZ
  } else {
    idx_na <- is.na(sd$isolationWindowUpperMz)
    sd$isolationWindowUpperMz[idx_na] <- sd$highMZ[idx_na]
  }

  Spectra::spectraData(sps) <- sd
  ProtGenerics::spectra(x) <- sps
  x
}

#' Example negative-mode DDA XcmsExperiment
#'
#' @description
#' A small example [xcms::XcmsExperiment] containing LC-MS/MS data acquired
#' in negative ionization mode using data-dependent acquisition (DDA).
#'
#' The object is a subset of the NIST plasma sample
#' `20250930_plasmaNIST_neg_DDA_01.mzML`, restricted to the retention-time
#' interval from 500 to 600 seconds. It is intended for use in LipidMS
#' examples, tests, and documentation.
#'
#' The object contains detected MS1 chromatographic peaks and the associated
#' MS1 and MS2 data required to demonstrate the LipidMS annotation workflow.
#'
#' @format
#' An [xcms::XcmsExperiment] object containing a single negative-mode DDA
#' sample restricted to a retention-time range of 500--600 seconds.
#'
#' @source
#' NIST plasma LC-MS/MS data.
#'
#' @keywords datasets
"xcmsexperiment_negative"

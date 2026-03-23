# LipidMSapp()
#' LipidMS shiny app
#'
#' Interactive UI for LipidMS
#' 
#' @param max_upload_mb max size to upload
#'
#' @examples
#' \dontrun{
#' # example data files can be download from github.com/maialba3/LipidMSv2.0_exampleFiles
#'
#' library(LipidMS)
#' LipidMSapp()
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
LipidMSapp <- function(max_upload_mb = getOption("LipidMS.shiny.maxRequestSizeMB", Inf)) {
  for (pkg in c("shiny", "shinythemes", "shinyjs", "plotly")){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' must be installed to use this function."),
           call. = FALSE)
    }
  }
  
  # env var opcional para override
  env_mb <- Sys.getenv("LIPIDMS_SHINY_MAX_MB", unset = "")
  if (nzchar(env_mb)) suppressWarnings(max_upload_mb <- as.numeric(env_mb))
  
  # bytes (Inf si quieres ilimitado)
  bytes <- if (is.finite(max_upload_mb)) max_upload_mb * 1024^2 else Inf
  
  old <- getOption("shiny.maxRequestSize")
  on.exit(options(shiny.maxRequestSize = old), add = TRUE)
  options(shiny.maxRequestSize = bytes)
  
  if (interactive()){
    shiny::runApp(system.file('LipidMSapp', package='LipidMS'))
  }else {
    stop("Only available for interactive sessions")
  }
}

# LipidMSapp()
#' LipidMS shiny app
#'
#' Interactive UI for LipidMS
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
LipidMSapp <- function() {
  if (interactive()){
    shiny::runApp(system.file('LipidMSapp', package='LipidMS'), quiet = TRUE)
  }else {
    stop("Only available for interactive sessions")
  }
}

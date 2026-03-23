# .with_future_plan
#' Temporarily set a future plan and globals max size while running a block
#'
#' Internal helper to encapsulate future::plan(...) and future.globals.maxSize
#' for a given expression, restoring the user's state on exit.
#'
#' @param workers ncores
#' @param globals_gb max size for globals
#' @param expr expr
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
.with_future_plan <- function(workers = 1L,
                              globals_gb = 24,
                              expr) {
  old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)
  old_glob <- getOption("future.globals.maxSize"); on.exit(options(future.globals.maxSize = old_glob), add = TRUE)
  old_msg  <- getOption("future.startup.message"); on.exit(options(future.startup.message = old_msg), add = TRUE)
  
  options(future.globals.maxSize = as.numeric(globals_gb) * 1024^3)
  options(future.startup.message = FALSE)
  
  maxw <- suppressWarnings(as.integer(future::availableCores()))
  if (is.na(maxw)) maxw <- 1L
  w <- max(1L, min(as.integer(workers), maxw))
  
  if (w > 1L) {
    cl <- parallelly::makeClusterPSOCK(
      workers      = w,
      rscript_args = c("--vanilla", "--default-packages=methods"),
      outfile      = ""   # keep stdout/stderr visible
    )
    
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterCall(cl, function(pkgs) {
          mute <- function(expr) {
            tf <- tempfile()
            con <- file(tf, open = "w")
            sink(con); sink(con, type = "message")
            on.exit({
              sink(type = "message"); sink()
              try(close(con), silent = TRUE)
              try(unlink(tf), silent = TRUE)
            }, add = TRUE)
            force(expr)
            invisible(NULL)
          }
        })
    
    future::plan(future::cluster, workers = cl)
    
  } else {
    future::plan(future::sequential)
  }
  
  force(expr)
}

#============================================================================#
# UI helpers
#============================================================================#

options(shiny.maxRequestSize = 610000000)  # Per request; job quota is validated separately.

.param_row <- function(title, help, input_ui) {
  fluidRow(
    column(
      6,
      h5(strong(title)),
      h6(em(help), style = "color:grey")
    ),
    column(3, input_ui)
  )
}

.param_row_2 <- function(title, help, input_ui_1, input_ui_2,
                         label_1 = "min", label_2 = "max") {
  fluidRow(
    column(
      6,
      h5(strong(title)),
      h6(em(help), style = "color:grey")
    ),
    column(3, input_ui_1),
    column(3, input_ui_2)
  )
}

.section <- function(title, ...) {
  div(class = "section-box", h4(strong(title)), ...)
}

.app_title <- function() {
  div(
    style = "text-align:left; margin:5px 0 15px 0; padding-left:15px;",
    tags$h4(
      "LipidMS: lipid annotation for untargeted LC-MS/MS lipidomics",
      style = "color:#669999; font-weight:600; margin:0;"
    )
  )
}

.navbar_logo <- function() {
  img(
    src = "iconohorizontal.png",
    style = "height:60px; width:110; margin-top:0px; margin-left:-15px;"
  )
}

prevent_sleep <- function() {

  if (.Platform$OS.type != "windows") {
    return(NULL)
  }

  if (!requireNamespace("processx", quietly = TRUE)) {
    stop("Package 'processx' is required.")
  }

  script <- tempfile(fileext = ".ps1")

  writeLines(c(
    'Add-Type -TypeDefinition @"',
    'using System;',
    'using System.Runtime.InteropServices;',
    'public class Power {',
    '  [DllImport("kernel32.dll")]',
    '  public static extern uint SetThreadExecutionState(uint flags);',
    '}',
    '"@',
    '',
    'while ($true) {',
    '  [Power]::SetThreadExecutionState(0x80000001) | Out-Null',
    '  Start-Sleep -Seconds 30',
    '}'
  ), script)

  proc <- tryCatch(
    processx::process$new(
      command = "powershell.exe",
      args = c(
        "-NoProfile",
        "-NonInteractive",
        "-ExecutionPolicy", "Bypass",
        "-File", script
      ),
      stdout = "|",
      stderr = "|",
      cleanup_tree = TRUE
    ),
    error = function(e) {
      unlink(script)
      stop(e)
    }
  )

  Sys.sleep(1)

  if (!proc$is_alive()) {
    error_msg <- paste(proc$read_error_lines(), collapse = "\n")
    unlink(script)
    stop("Could not prevent sleep: ", error_msg)
  }

  list(process = proc, script = script)
}


allow_sleep <- function(handle) {

  if (is.null(handle)) {
    return(invisible(NULL))
  }

  if (handle$process$is_alive()) {
    handle$process$kill()
  }

  unlink(handle$script)

  invisible(NULL)
}

#===============================================================================
# LipidMS Shiny/server job-preparation helpers
#===============================================================================
# These functions are used by app.R to prepare a self-contained job
# folder and generate a runnable LipidMS workflow script from a template.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

.as_r_literal <- function(x) {
  paste(deparse(x), collapse = "\n")
}

.replace_template_tokens <- function(template, values) {
  for (nm in names(values)) {
    template <- gsub(
      pattern = paste0("#", nm, "#"),
      replacement = .as_r_literal(values[[nm]]),
      x = template,
      fixed = TRUE
    )
  }
  template
}

.sanitize_job_name <- function(x) {
  x <- trimws(as.character(x %||% "LipidMS_job"))
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) x <- "LipidMS_job"
  x
}

.new_job_id <- function(jobname = "LipidMS_job") {
  prefix <- .sanitize_job_name(jobname)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  if (requireNamespace("uuid", quietly = TRUE)) {
    paste(prefix, stamp, uuid::UUIDgenerate(), sep = "_")
  } else {
    rand <- paste(sample(c(letters, LETTERS, 0:9), 8, replace = TRUE), collapse = "")
    paste(prefix, stamp, rand, sep = "_")
  }
}

# .ensure_dir <- function(path) {
#   dir.create(path, recursive = TRUE, showWarnings = FALSE)
#   normalizePath(path, winslash = "/", mustWork = FALSE)
# }

.normalize_sample_file <- function(x) {
  x <- basename(as.character(x))
  x <- sub("\\.mzXML$", ".mzML", x, ignore.case = TRUE)
  x
}

.normalize_acquisition_mode <- function(x, sample = NULL) {
  y <- toupper(trimws(as.character(x)))
  y[is.na(y) | !nzchar(y)] <- "FULLSCAN"
  y[y %in% c("MS", "MS1", "FULLSCAN", "FULLSCAN/MS", "FULLSCAN_MS")] <- "FULLSCAN"
  y[y %in% c("DIA", "SWATH")] <- "DIA"
  y[y %in% c("DDA")] <- "DDA"

  if (!is.null(sample)) {
    s <- basename(as.character(sample))
    y[grepl("DDA", s, ignore.case = TRUE)] <- "DDA"
    y[grepl("DIA|SWATH", s, ignore.case = TRUE)] <- "DIA"
  }

  y[!y %in% c("FULLSCAN", "DDA", "DIA")] <- "FULLSCAN"
  y
}

.read_lipidms_metadata <- function(path, sep = ",") {

  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required to read Excel metadata.", call. = FALSE)
    }
    metadata <- as.data.frame(readxl::read_excel(path))
  } else {
    metadata <- utils::read.table(
      path,
      sep = sep,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "\"",
      comment.char = ""
    )
  }

  colnames(metadata) <- trimws(colnames(metadata))

  cn <- tolower(colnames(metadata))
  cn <- gsub("[ .-]+", "_", cn)

  colnames(metadata) <- cn

  sample_aliases <- c(
    "sample",
    "samples",
    "file",
    "filename",
    "file_name",
    "mzml",
    "mzml_file",
    "raw_file"
  )

  acquisition_aliases <- c(
    "acquisition_mode",
    "acquisition",
    "mode",
    "ms_mode",
    "type",
    "data_acquisition"
  )

  sample_type_aliases <- c(
    "sample_type",
    "sampletype",
    "sample_group",
    "sample_groups",
    "samplegroup",
    "samplegroups",
    "group",
    "groups",
    "group_id",
    "groupid"
  )

  sample_col <- intersect(sample_aliases, colnames(metadata))[1]
  acquisition_col <- intersect(acquisition_aliases, colnames(metadata))[1]
  sample_type_col <- intersect(sample_type_aliases, colnames(metadata))[1]

  if (is.na(sample_col)) {
    stop(
      "Metadata is missing a sample/file column or the delimiter is wrong. Accepted names include: ",
      paste(sample_aliases, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.na(acquisition_col)) {
    stop(
      "Metadata is missing an acquisition column or the delimiter is wrong. Accepted names include: ",
      paste(acquisition_aliases, collapse = ", "),
      call. = FALSE
    )
  }

  metadata$sample <- as.character(metadata[[sample_col]])
  metadata$acquisition_mode <- .normalize_acquisition_mode(
    metadata[[acquisition_col]],
    sample = metadata$sample
  )

  metadata$sample <- .normalize_sample_file(metadata$sample)

  # Normalize sample type
  if (!is.na(sample_type_col)) {

    metadata$sample_type <- trimws(
      as.character(metadata[[sample_type_col]])
    )

    metadata$sample_type[
      is.na(metadata$sample_type) |
        !nzchar(metadata$sample_type)
    ] <- "sample"

  } else {

    metadata$sample_type <- rep("sample", nrow(metadata))

  }

  metadata
}

.match_metadata_and_files <- function(metadata, uploaded_files) {
  if (is.null(uploaded_files) || !nrow(uploaded_files)) {
    stop("No mzML files were provided.", call. = FALSE)
  }

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  uploaded_files <- as.data.frame(uploaded_files, stringsAsFactors = FALSE)

  if (!all(c("datapath", "name") %in% colnames(uploaded_files))) {
    stop("uploaded_files must contain columns 'datapath' and 'name'.", call. = FALSE)
  }

  metadata$sample_original <- as.character(metadata$sample)
  metadata$sample_basename <- .normalize_sample_file(metadata$sample)
  uploaded_files$file_basename <- .normalize_sample_file(uploaded_files$name)

  common <- intersect(metadata$sample_basename, uploaded_files$file_basename)
  common <- common[!duplicated(common)]

  if (!length(common)) {
    stop(
      paste(
        "No common samples were found between metadata$sample and uploaded files.",
        "Check that metadata$sample contains mzML file names."
      ),
      call. = FALSE
    )
  }

  metadata <- metadata[match(common, metadata$sample_basename), , drop = FALSE]
  uploaded_files <- uploaded_files[match(common, uploaded_files$file_basename), , drop = FALSE]

  metadata$sample <- uploaded_files$file_basename
  metadata$acquisition_mode <- .normalize_acquisition_mode(
    metadata$acquisition_mode,
    sample = metadata$sample
  )

  metadata$sample_basename <- NULL
  rownames(metadata) <- NULL

  list(
    metadata = metadata,
    files = uploaded_files,
    missing_in_files = setdiff(unique(.normalize_sample_file(metadata$sample_original)), common),
    missing_in_metadata = setdiff(unique(.normalize_sample_file(uploaded_files$name)), common)
  )
}

.collect_lipidms_params <- function(input) {
  list(
    jobname = .sanitize_job_name(input$jobname %||% "LipidMS_job"),
    polarity = input$to_process_polarity %||% "negative",
    lipid_db_source = input$lipid_db_source %||% "default",

    cwp_ms1_ppm = as.numeric(input$cwp_ms1_ppm %||% 25),
    cwp_ms1_peakwidth_min = as.numeric(input$cwp_ms1_peakwidth_min %||% 5),
    cwp_ms1_peakwidth_max = as.numeric(input$cwp_ms1_peakwidth_max %||% 50),
    cwp_ms1_snthresh = as.numeric(input$cwp_ms1_snthresh %||% 3),
    cwp_ms1_prefilter_n = as.integer(input$cwp_ms1_prefilter_n %||% 5),
    cwp_ms1_prefilter_intensity = as.numeric(input$cwp_ms1_prefilter_intensity %||% 100),
    cwp_ms1_mzCenterFun = input$cwp_ms1_mzCenterFun %||% "wMean",
    cwp_ms1_integrate = as.integer(input$cwp_ms1_integrate %||% 2),

    cwp_ms2_ppm = as.numeric(input$cwp_ms2_ppm %||% 25),
    cwp_ms2_peakwidth_min = as.numeric(input$cwp_ms2_peakwidth_min %||% 3),
    cwp_ms2_peakwidth_max = as.numeric(input$cwp_ms2_peakwidth_max %||% 50),
    cwp_ms2_snthresh = as.numeric(input$cwp_ms2_snthresh %||% 2),
    cwp_ms2_prefilter_n = as.integer(input$cwp_ms2_prefilter_n %||% 3),
    cwp_ms2_prefilter_intensity = as.numeric(input$cwp_ms2_prefilter_intensity %||% 50),
    cwp_ms2_mzCenterFun = input$cwp_ms2_mzCenterFun %||% "wMean",
    cwp_ms2_integrate = as.integer(input$cwp_ms2_integrate %||% 2),
    cwp_ms2_mzdiff = as.numeric(input$cwp_ms2_mzdiff %||% -0.001),

    grouping1_bw = as.numeric(input$grouping1_bw %||% 5),
    grouping1_minFraction = as.numeric(input$grouping1_minFraction %||% 0.9),
    grouping1_binSize = as.numeric(input$grouping1_binSize %||% 0.01),

    paramRT_minFraction = as.numeric(input$paramRT_minFraction %||% 0.9),
    paramRT_extraPeaks = as.integer(input$paramRT_extraPeaks %||% 1),
    paramRT_smooth = input$paramRT_smooth %||% "loess",
    paramRT_span = as.numeric(input$paramRT_span %||% 0.4),
    paramRT_family = input$paramRT_family %||% "gaussian",
    paramRT_subsetAdjust = input$paramRT_subsetAdjust %||% "average",

    grouping2_bw = as.numeric(input$grouping2_bw %||% 5),
    grouping2_minFraction = as.numeric(input$grouping2_minFraction %||% 0.5),
    grouping2_binSize = as.numeric(input$grouping2_binSize %||% 0.01),

    dmzprecursor = as.numeric(input$dmzprecursor %||% 10),
    dmzproducts = as.numeric(input$dmzproducts %||% 15),
    rttol = as.numeric(input$rttol %||% 3),
    min_ms2_intensity = as.numeric(input$min_ms2_intensity %||% 100),

    isotope_mzabs = as.numeric(input$isotope_mzabs %||% 0.005),
    isotope_maxcharge = as.integer(input$isotope_maxcharge %||% 2),
    isotope_maxiso = as.integer(input$isotope_maxiso %||% 5),
    isotope_minfrac = as.numeric(input$isotope_minfrac %||% 0.5),
    isotope_rttol = as.numeric(input$isotope_rttol %||% 3),
    intval = input$intval %||% "into",
    remove_isotopes = isTRUE(input$remove_isotopes %||% TRUE),
    only_confirmed_monoisotopic = isTRUE(input$only_confirmed_monoisotopic %||% FALSE),

    paramdia_expandRt = as.numeric(input$paramdia_expandRt %||% 0),
    paramdia_minCor = as.numeric(input$paramdia_minCor %||% 0.7),
    paramdia_diffRt = as.numeric(input$paramdia_diffRt %||% 3),
    paramdia_intensity = input$paramdia_intensity %||% "into",

    type_lipidclasses = input$type_lipidclasses %||% "lipidclasses_default",
    classes_sep = input$classes_sep %||% ";",
    adducts_sep = input$adducts_sep %||% ",",

    parallel_mode = input$parallel_mode,
    workers = input$parallel_workers
  )
}

.validate_mzml_files <- function(files) {
  if (is.null(files) || !is.data.frame(files) || nrow(files) < 1L) {
    stop("mzML files are required.", call. = FALSE)
  }
  if (nrow(files) > 200L) {
    stop("A maximum of 200 mzML files can be processed per job.", call. = FALSE)
  }
  if (is.null(files$name) || anyNA(files$name) ||
      any(!grepl("\\.mzML$", files$name, ignore.case = TRUE))) {
    stop("Only .mzML files are accepted for LC-MS input.", call. = FALSE)
  }
  invisible(TRUE)
}


.validate_job_uploads <- function(input, require_all = TRUE) {
  groups <- c("mzml_files", "metadata_file", "lipid_db_config_file")
  if (require_all) {
    if (is.null(input$mzml_files) || is.null(input$metadata_file)) {
      stop("mzML files and metadata are required.", call. = FALSE)
    }
  }
  entries <- list()
  for (kind in groups) {
    files <- input[[kind]]
    if (is.null(files)) next
    if (!is.data.frame(files) || !all(c("name", "datapath") %in% names(files)) ||
        !nrow(files) || anyNA(files$name) || anyNA(files$datapath)) {
      stop("Invalid uploaded file information.", call. = FALSE)
    }
    if (kind != "mzml_files" && nrow(files) != 1L) {
      stop("Select only one file for each auxiliary input.", call. = FALSE)
    }
    if (kind == "mzml_files") .validate_mzml_files(files)
    for (i in seq_len(nrow(files))) {
      entries[[length(entries) + 1L]] <- list(
        kind = kind, name = files$name[[i]],
        path = normalizePath(files$datapath[[i]], winslash = "/", mustWork = TRUE)
      )
    }
  }
  invisible(TRUE)
}

.validate_metadata <- function(metadata) {

  if (nrow(metadata) < 2L) {
    stop(
      "Metadata must contain at least 2 samples.",
      call. = FALSE
    )
  }

  acquisition <- toupper(
    trimws(
      as.character(metadata$acquisition_mode)
    )
  )

  has_ms2 <- any(
    acquisition %in% c("DDA", "DIA", "SWATH")
  )

  if (!has_ms2) {
    stop(
      paste(
        "At least one sample must be acquired in",
        "DDA, DIA or SWATH mode."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.available_cores <- function() {
  n <- parallel::detectCores(logical = FALSE)
  if (is.na(n) || n < 1L) n <- 1L
  as.integer(n)
}

make_bpparam <- function(parallel_mode = "auto", workers = NULL) {

  parallel_mode <- match.arg(
    parallel_mode,
    c("auto", "manual", "serial")
  )

  available <- .available_cores()

  if (parallel_mode == "serial") {
    return(BiocParallel::SerialParam())
  }

  if (parallel_mode == "auto") {
    workers <- max(1L, ceiling(available/2))
  } else {
    if (is.null(workers) || length(workers) != 1L ||
        !is.finite(workers) || workers < 1) {
      stop("'workers' must be a positive number in manual mode.")
    }
  }

  workers <- max(1L, min(as.integer(workers), available))

  if (workers == 1L) {
    return(BiocParallel::SerialParam())
  }

  if (.Platform$OS.type == "windows") {
    return(BiocParallel::SnowParam(
      workers = workers,
      type = "SOCK"
    ))
  }

  BiocParallel::MulticoreParam(workers = workers)
}

#-------------------------------------------------------------------------------
# Tab UI helpers
#-------------------------------------------------------------------------------

.tab_data_import <- function() {
  tabPanel(
    "Data import", value = "tab1",
    fluidRow(
      column(
        4,
        .section(
          "Input",
          radioButtons(
            "to_process_polarity", "Polarity",
            choices = c("Positive" = "positive", "Negative" = "negative"),
            selected = "negative"
          ),
          fileInput(
            "mzml_files", "Choose mzML files",
            multiple = TRUE,
            accept = c(".mzML")
          ),
          uiOutput("mzml_validation_message"),
          helpText("Maximum 200 mzML files; 600 MB per file; 25 GB total per job."),
          fileInput(
            "metadata_file", "Metadata file (maximum 1 MB)",
            multiple = FALSE,
            accept = c(".csv", ".txt", ".tsv", ".xlsx", ".xls")
          ),
          uiOutput("metadata_validation_message"),
          helpText(
            "Required columns: sample and acquisition_mode.",
            tags$br(),
            "sample: mzML file name.",
            tags$br(),
            "acquisition_mode: FULLSCAN, DDA, DIA, SWATH."
          ),
          selectInput(
            "metadata_sep", "Metadata delimiter",
            choices = c("Comma (,)" = ",", "Semicolon (;)" = ";", "Tab" = "\t"),
            selected = ","
          )
        ),
        .section(
          "Lipid database configuration",
          radioButtons(
            "lipid_db_source", NULL,
            choices = c(
              "Default LipidMS database" = "default",
              "Upload custom Excel configuration" = "custom_excel"
            ),
            selected = "default"
          ),

          downloadButton(
            "download_default_lipid_db_config",
            "Download default Excel template"
          ),

          conditionalPanel(
            condition = "input.lipid_db_source == 'custom_excel'",
            fileInput(
              "lipid_db_config_file",
              "Upload edited lipid database configuration (.xls/.xlsx, maximum 1 MB)",
              multiple = FALSE,
              accept = c(".xls", ".xlsx")
            ),
            uiOutput("lipid_db_validation_message"),
            helpText(
              "The Excel file must contain the sheets: lipid_classes_config and chain_pools."
            )
          )
        ),
        uiOutput("upload_validation_message"),
        actionButton("JumpTo2", "Next >", width = "120px")
      ),
      column(
        8,
        h3("Import summary"),
        tableOutput("import_summary"),
        h4("Matching preview"),
        tableOutput("matching_preview")
      )
    )
  )
}

.tab_ms1 <- function() {
  tabPanel(
    "MS1 peak-picking", value = "tab2",
    .section(
      "MS1 CentWaveParam",
      .param_row("ppm", "Mass accuracy for MS1 peak detection.", numericInput("cwp_ms1_ppm", NULL, value = 25, min = 0, max = 100, step = 1)),
      .param_row_2("peakwidth", "Expected chromatographic peak width in seconds.", numericInput("cwp_ms1_peakwidth_min", "min", value = 5, min = 0, max = 500), numericInput("cwp_ms1_peakwidth_max", "max", value = 50, min = 0, max = 500)),
      .param_row("snthresh", "Signal-to-noise threshold.", numericInput("cwp_ms1_snthresh", NULL, value = 3, min = 0, max = 100)),
      .param_row_2("prefilter", "Minimum number of scans and intensity before peak detection.", numericInput("cwp_ms1_prefilter_n", "n", value = 5, min = 0, max = 100), numericInput("cwp_ms1_prefilter_intensity", "intensity", value = 100, min = 0, max = 1e8)),
      .param_row("mzCenterFun", "Function used to calculate the m/z center.", selectInput("cwp_ms1_mzCenterFun", NULL, choices = c("wMean", "mean", "apex"), selected = "wMean")),
      .param_row("integrate", "Integration method.", selectInput("cwp_ms1_integrate", NULL, choices = c("1" = 1, "2" = 2), selected = 2))
    ),
    actionButton("GoBackTo1", "< Previous"),
    actionButton("JumpTo3", "Next >")
  )
}

.tab_ms2 <- function() {
  tabPanel(
    "MS2/DIA", value = "tab3",
    .section(
      "MS2 CentWaveParam / DIA extraction",
      .param_row("ppm", "Mass accuracy for MS2 chromatographic peak detection.", numericInput("cwp_ms2_ppm", NULL, value = 25, min = 0, max = 100)),
      .param_row_2("peakwidth", "Expected MS2 chromatographic peak width in seconds.", numericInput("cwp_ms2_peakwidth_min", "min", value = 3, min = 0, max = 500), numericInput("cwp_ms2_peakwidth_max", "max", value = 50, min = 0, max = 500)),
      .param_row("snthresh", "Signal-to-noise threshold.", numericInput("cwp_ms2_snthresh", NULL, value = 2, min = 0, max = 100)),
      .param_row_2("prefilter", "Minimum number of scans and intensity before MS2 peak detection.", numericInput("cwp_ms2_prefilter_n", "n", value = 3, min = 0, max = 100), numericInput("cwp_ms2_prefilter_intensity", "intensity", value = 50, min = 0, max = 1e8)),
      .param_row("mzCenterFun", "Function used to calculate the m/z center.", selectInput("cwp_ms2_mzCenterFun", NULL, choices = c("wMean", "mean", "apex"), selected = "wMean")),
      .param_row("integrate", "Integration method.", selectInput("cwp_ms2_integrate", NULL, choices = c("1" = 1, "2" = 2), selected = 2)),
      .param_row("mzdiff", "Minimum m/z difference for overlapping peaks.", numericInput("cwp_ms2_mzdiff", NULL, value = -0.001, min = -1, max = 1, step = 0.001)),
      tags$hr(),
      .param_row("expandRt", "Retention-time expansion for DIA reconstruction.", numericInput("paramdia_expandRt", NULL, value = 0, min = 0, max = 100)),
      .param_row("minCor", "Minimum correlation for DIA reconstruction.", numericInput("paramdia_minCor", NULL, value = 0.7, min = 0, max = 1, step = 0.05)),
      .param_row("diffRt", "Maximum retention-time difference for DIA reconstruction.", numericInput("paramdia_diffRt", NULL, value = 3, min = 0, max = 100)),
      .param_row("intensity", "Intensity column for DIA reconstruction.", selectInput("paramdia_intensity", NULL, choices = c("into", "maxo"), selected = "into"))
    ),
    actionButton("GoBackTo2", "< Previous"),
    actionButton("JumpTo4", "Next >")
  )
}

.tab_grouping <- function() {
  tabPanel(
    "Grouping / RT", value = "tab4",
    .section(
      "Initial grouping",
      .param_row("bw", "Bandwidth for initial peak-density grouping.", numericInput("grouping1_bw", NULL, value = 5, min = 0, max = 100)),
      .param_row("minFraction", "Minimum sample fraction for initial grouping.", numericInput("grouping1_minFraction", NULL, value = 0.8, min = 0, max = 1, step = 0.05)),
      .param_row("binSize", "m/z bin size for initial grouping.", numericInput("grouping1_binSize", NULL, value = 0.01, min = 0, max = 1, step = 0.001))
    ),
    .section(
      "Retention-time correction",
      .param_row("minFraction", "Minimum fraction for peak groups used in RT correction.", numericInput("paramRT_minFraction", NULL, value = 0.8, min = 0, max = 1, step = 0.05)),
      .param_row("extraPeaks", "Extra peaks allowed in peak groups.", numericInput("paramRT_extraPeaks", NULL, value = 1, min = 0, max = 20)),
      .param_row("smooth", "Smoothing method.", selectInput("paramRT_smooth", NULL, choices = c("loess", "linear"), selected = "loess")),
      .param_row("span", "LOESS span.", numericInput("paramRT_span", NULL, value = 0.4, min = 0, max = 1, step = 0.05)),
      .param_row("family", "LOESS family.", selectInput("paramRT_family", NULL, choices = c("gaussian", "symmetric"), selected = "gaussian")),
      .param_row("subsetAdjust", "How to adjust samples outside the subset.", selectInput("paramRT_subsetAdjust", NULL, choices = c("average", "previous"), selected = "average"))
    ),
    .section(
      "Final grouping",
      .param_row("bw", "Bandwidth for final grouping.", numericInput("grouping2_bw", NULL, value = 5, min = 0, max = 100)),
      .param_row("minFraction", "Minimum sample fraction for final grouping.", numericInput("grouping2_minFraction", NULL, value = 0.15, min = 0, max = 1, step = 0.05)),
      .param_row("binSize", "m/z bin size for final grouping.", numericInput("grouping2_binSize", NULL, value = 0.005, min = 0, max = 1, step = 0.001))
    ),
    actionButton("GoBackTo3", "< Previous"),
    actionButton("JumpTo5", "Next >")
  )
}

.tab_annotation <- function() {
  tabPanel(
    "Annotation", value = "tab5",
    .section(
      "Annotation parameters",
      .param_row("dmzprecursor", "Precursor m/z tolerance in ppm.", numericInput("dmzprecursor", NULL, value = 10, min = 0, max = 100)),
      .param_row("dmzproducts", "Product/fragment m/z tolerance in ppm.", numericInput("dmzproducts", NULL, value = 15, min = 0, max = 100)),
      .param_row("rttol", "Retention-time tolerance for coelution/MS2 matching.", numericInput("rttol", NULL, value = 3, min = 0, max = 100)),
      .param_row("min_ms2_intensity", "Minimum MS2 fragment intensity kept.", numericInput("min_ms2_intensity", NULL, value = 100, min = 0, max = 1e8)),
      tags$hr(),
      .param_row("isotope_mzabs", "Absolute m/z tolerance for isotope annotation.", numericInput("isotope_mzabs", NULL, value = 0.005, min = 0, max = 1, step = 0.001)),
      .param_row("isotope_maxcharge", "Maximum isotope charge.", numericInput("isotope_maxcharge", NULL, value = 2, min = 1, max = 10)),
      .param_row("isotope_maxiso", "Maximum isotope rank.", numericInput("isotope_maxiso", NULL, value = 5, min = 1, max = 20)),
      .param_row("isotope_minfrac", "Minimum isotope fraction.", numericInput("isotope_minfrac", NULL, value = 0.5, min = 0, max = 1, step = 0.05)),
      .param_row("isotope_rttol", "RT tolerance for isotope annotation.", numericInput("isotope_rttol", NULL, value = 3, min = 0, max = 100)),
      .param_row("intval", "Intensity value used by CAMERA/isotope annotation.", selectInput("intval", NULL, choices = c("into", "maxo", "intb"), selected = "into")),
      checkboxInput("remove_isotopes", "Remove isotope peaks before MS1 candidate search", value = TRUE),
      checkboxInput("only_confirmed_monoisotopic", "Export only features with confirmed monoisotopic peak", value = TRUE)
    ),
    actionButton("GoBackTo4", "< Previous"),
    actionButton("JumpTo6", "Next >")
  )
}

.tab_parallelization <- function() {

  available <- .available_cores()

  tabPanel(
    "Parallelization",
    value = "tab6",

    .section(
      "Parallel processing",

      radioButtons(
        "parallel_mode",
        "Execution mode",
        choices = c(
          "Sequential (1 core)" = "serial",
          "Automatic" = "auto",
          "Manual" = "manual"
        ),
        selected = "auto"
      ),

      helpText(
        paste("Available physical cores:", available)
      ),

      conditionalPanel(
        condition = "input.parallel_mode == 'manual'",

        numericInput(
          "parallel_workers",
          "Number of workers",
          value = max(1L, ceiling(available/2)),
          min = 1,
          max = available,
          step = 1
        )
      )
    ),

    actionButton("GoBackTo5", "< Previous"),
    actionButton("JumpTo7", "Next >")
  )
}

.max_job_name_chars <- 29L

.job_name_length_message <- function(count) {
  sprintf(
    "Job name has %d characters. The maximum is %d. Please shorten it before running the job.",
    count,
    .max_job_name_chars
  )
}

.tab_server_job <- function() {
  tabPanel(
    "Server job", value = "tab7",
    fluidRow(
      column(
        5,
        .section(
          "Job settings",
          textInput("jobname", "Job name", value = "LipidMS_job"),
          uiOutput("job_name_length_warning"),
        ),
        tags$hr(),
        actionButton("Run_analysis", "Run job", class = "btn-success",
                     icon = icon("gears"), disabled = "disabled"),
        uiOutput("run_validation_message"),
        actionButton("GoBackTo5", "< Previous")
      ),
      column(
        7,
        h3("Selected parameters"),
        uiOutput("parameters_output"),
        downloadButton("download_parameters_run", "Download parameters")
      )
    )
  )
}

.tab_about <- function() {
  external_link <- function(label, url) {
    tags$a(label, href = url, target = "_blank", rel = "noopener noreferrer")
  }
  tabPanel(
    "About", value = "tab_about",
    fluidRow(
      column(
        10, offset = 0, style = "overflow-wrap:anywhere;",
        .section(
          "Contact",
          tags$ul(
            tags$li(tags$a("maribel_alcoriza@iislafe.es", href = "mailto:maribel_alcoriza@iislafe.es")),
            tags$li(tags$a("juancarlos_garcia@iislafe.es", href = "mailto:juancarlos_garcia@iislafe.es"))
          )
        ),
        .section(
          "Useful links",
          tags$ul(
            tags$li(external_link("LipidMS — GitHub repository", "https://github.com/maialba3/LipidMS")),
            tags$li(external_link("UBMP — Research Group", "https://iislafe.portalinvestigacion.com/grupos/356")),
            tags$li(external_link("Unidad Analítica — IIS La Fe", "https://www.iislafe.es/es/apoyo-a-la-investigacion/plataformas-cientifico-tecnologicas/unidad-analitica/"))
          )
        ),
        .section(
          "References",
          tags$ol(
            tags$li(
              "Alcoriza-Balaguer MI et al. (2019). LipidMS: An R Package for Lipid Annotation in Untargeted Liquid Chromatography-Data Independent Acquisition-Mass Spectrometry Lipidomics. ",
              tags$em("Analytical Chemistry"), ", 91(1), 836–845. ",
              external_link("https://doi.org/10.1021/acs.analchem.8b03409", "https://doi.org/10.1021/acs.analchem.8b03409")
            ),
            tags$li(
              "Alcoriza-Balaguer MI et al. (2022). LipidMS 3.0: an R-package and a web-based tool for LC-MS/MS data processing and lipid annotation. ",
              tags$em("Bioinformatics"), ", 38(20), 4826–4828. ",
              external_link("https://doi.org/10.1093/bioinformatics/btac581", "https://doi.org/10.1093/bioinformatics/btac581")
            )
          )
        )
      )
    )
  )
}

.tab_generated_job <- function() {

  tabPanel(
    "Results", value = "tab8",

    fluidRow(

      # Left column: downloads
      column(
        4,
        style = "overflow-wrap:anywhere;",

        .section(
          "Download results",
          div(
            style = "display: flex;
            flex-direction: column;
            align-items: flex-start;
            gap: 12px;",

            uiOutput("download_featureTable_ui"),
            uiOutput("download_chromPeakTable_ui"),
            uiOutput("download_xcmsObject_ui")
          )
        )
      ),

      # Right column: result tables
      column(
        8,
        tabsetPanel(
          id = "results_tabs",
          type = "tabs",
          selected = "features",

          tabPanel(
            "Feature table",
            value = "features",
            br(),
            uiOutput("featureTable_ui")
          ),

          tabPanel(
            "Chromatographic peaks",
            value = "peaks",
            br(),
            uiOutput("chromPeakTable_ui")
          ),

          tabPanel(
            "MS1 candidates",
            value = "ms1",
            br(),
            h5("Candidate lipids based only on MS1 information"),
            uiOutput("ms1Candidates_ui")
          ),

          tabPanel(
            "MS2 evidence",
            value = "ms2",
            br(),
            h5("All MS2 evidence for annotated and non-annotated lipids"),
            uiOutput("ms2Evidence_ui")
          )
        )
      )
    )
  )
}

#-------------------------------------------------------------------------------
# User interface
#-------------------------------------------------------------------------------

.lipidms_theme <- bslib::bs_theme(version = 5, bootswatch = "flatly")

ui <- fluidPage(
  theme = .lipidms_theme,
  tags$head(
    tags$script(src = "upload_validation.js"),
    tags$style(HTML("
  body {
    padding-top: 10px !important;
    font-family: 'Inter', Arial, sans-serif !important;
  }

  #upload_status {
    position: fixed; bottom: 18px; right: 18px; z-index: 1040;
    max-width: min(420px, calc(100vw - 36px)); max-height: 45vh; overflow: auto;
    background: var(--bs-body-bg, white); border: 1px solid #669999;
    border-radius: 6px; padding: 12px; box-shadow: 0 2px 10px #0002;
    font-size: 14px;
  }
  #upload_status[data-status='error'] { border-color: #b00020; color: #b00020; }
  #upload_status ul { margin-bottom: 0; padding-left: 20px; }

  .section-box {
    border: 1px solid #ddd;
    border-radius: 6px;
    padding: 15px;
    margin: 10px 0;
    background-color: #f9f9f9;
  }

  .shiny-notification {
    background-color: #fc804a;
    border: 1px solid #e74818;
    color: white;
    font-size: 15px;
    padding: 15px;
    border-radius: 8px;
  }

  .shiny-progress-notification .progress-message {
    display: block;
    font-size: 1em;
    font-weight: 600;
  }

  .shiny-progress-notification .progress-detail {
    display: block;
    color: #680d0d;
    font-size: 0.85em;
    margin-top: 4px;
  }

  .shiny-notification-message {
    background-color: #669973;
    border: 1px solid #4e7a5d;
  }

  input[type='checkbox'],
  input[type='radio'] {
    accent-color: #669999;
  }

  input[type='checkbox']:checked,
  input[type='radio']:checked {
    background-color: #669999 !important;
    border-color: #669999 !important;
  }

  input[type='number'] {
    height: 35px;
    font-size: 12px;
  }

  .local-script-preview {
    margin-bottom: 12px;
  }

  .local-script-preview > summary {
    cursor: pointer;
    list-style: none;
    display: flex;
    align-items: center;
    gap: 8px;
    color: #4e7a7a;
    font-weight: 600;
  }

  .local-script-preview > summary::-webkit-details-marker {
    display: none;
  }

  .local-script-preview-chevron {
    display: inline-block;
    transition: transform 0.15s ease;
  }

  .local-script-preview[open] .local-script-preview-chevron {
    transform: rotate(90deg);
  }

  /* BOTONES */

  .btn-primary,
  .btn-success {

    color: white !important;
    background-color: #669999 !important;
    border-color: #669999 !important;
  }

  .btn-primary:hover,
  .btn-success:hover {

    background-color: #5c8c8c !important;
    border-color: #5c8c8c !important;
  }

  .btn-primary:focus,
  .btn-success:focus,
  .btn-primary:active,
  .btn-success:active {

    background-color: #4e7a7a !important;
    border-color: #4e7a7a !important;
  }

  td {
    max-width: 420px;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
  }

  td:hover {
    overflow: visible;
    white-space: normal;
    background-color: #f9f9f9;
  }

  /* NAVBAR */

  .navbar {
    background-color: #669999 !important;
    border-color: #669999 !important;
    border-radius: 0 !important;
    box-shadow: none !important;
  }

  .navbar-brand {
    padding: 0 15px !important;
    display: flex !important;
    align-items: center !important;
  }

  .navbar .navbar-nav .nav-link {
    border-bottom: 0;
    color: white !important;
    font-size: 15px;
    font-weight: 500;
  }

  .navbar .navbar-nav .nav-link:hover {
    background-color: #5c8c8c !important;
    color: white !important;
  }

  .navbar .navbar-nav .nav-link.active,
  .navbar .navbar-nav .nav-link.active:hover,
  .navbar .navbar-nav .nav-link.active:focus {

    background-color: #4e7a7a !important;
    color: white !important;
  }

"))
  ),
  uiOutput("main_navbar"),
  tags$div(id = "upload_status", hidden = "hidden", role = "status", `aria-live` = "polite"),

  div(
    style = "
    width: 100%;
    text-align: center;
    margin-top: 40px;
    margin-bottom: 20px;
  ",

    tags$hr(
      style = "
      border: 0;
      border-top: 1px solid #d0d0d0;
      margin: 0 5% 25px 5%;
    "
    ),

    img(
      src = "logos_institucionales.png",
      alt = "Institutional partners",
      style = "
      width: 40%;
      height: 40%;
      margin-left: 5%;
    "
    )
  )
)


#============================================================================#
# Shiny application server
#============================================================================#
##############################################################################
server <- function(input, output, session) {

  output$main_navbar <- renderUI({
    tabs <- list(
      .tab_data_import(),
      .tab_ms1(),
      .tab_ms2(),
      .tab_grouping(),
      .tab_annotation(),
      .tab_parallelization(),
      .tab_server_job()
    )

    if (!is.null(job_info())) {
      tabs <- c(tabs, list(.tab_generated_job()))
    }
    tabs <- c(tabs, list(.tab_about()))

    do.call(
      navbarPage,
      c(
        list(
          title = .navbar_logo(),
          windowTitle = "LipidMS",
          id = "inTabset",
          selected = "tab1",
          theme = .lipidms_theme,
          header = .app_title()
        ),
        tabs
      )
    )
  })

  # updateTabsetPanel(
  #   session,
  #   "inTabset",
  #   selected = "tab1"
  # )

  observeEvent(input$JumpTo2, updateTabsetPanel(session, "inTabset", selected = "tab2"))
  observeEvent(input$GoBackTo1, updateTabsetPanel(session, "inTabset", selected = "tab1"))
  observeEvent(input$JumpTo3, updateTabsetPanel(session, "inTabset", selected = "tab3"))
  observeEvent(input$GoBackTo2, updateTabsetPanel(session, "inTabset", selected = "tab2"))
  observeEvent(input$JumpTo4, updateTabsetPanel(session, "inTabset", selected = "tab4"))
  observeEvent(input$GoBackTo3, updateTabsetPanel(session, "inTabset", selected = "tab3"))
  observeEvent(input$JumpTo5, updateTabsetPanel(session, "inTabset", selected = "tab5"))
  observeEvent(input$GoBackTo4, updateTabsetPanel(session, "inTabset", selected = "tab4"))
  observeEvent(input$JumpTo6, updateTabsetPanel(session, "inTabset", selected = "tab6"))
  observeEvent(input$GoBackTo5, updateTabsetPanel(session, "inTabset", selected = "tab5"))
  observeEvent(input$JumpTo7, updateTabsetPanel(session, "inTabset", selected = "tab7"))
  observeEvent(input$GoBackTo6, updateTabsetPanel(session, "inTabset", selected = "tab6"))

  # File checks depend only on actual uploads, not progress, navigation or parameters.
  upload_fields <- c("mzml_files", "metadata_file", "lipid_db_config_file")
  file_checks <- setNames(lapply(upload_fields, function(field) reactive({
    files <- input[[field]]
    if (is.null(files)) return(character())
    tryCatch({
      .validate_job_uploads(setNames(list(files), field), require_all = FALSE)
      character()
    }, error = function(e) conditionMessage(e))
  })), upload_fields)
  metadata_result <- reactive({
    if (is.null(input$metadata_file) || length(file_checks$metadata_file())) return(NULL)
    tryCatch({
      md <- .read_lipidms_metadata(input$metadata_file$datapath,
                                   sep = input$metadata_sep %||% ",")
      .validate_metadata(md)
      list(data = md, error = character())
    }, error = function(e) list(data = NULL, error = paste("Metadata validation failed:", conditionMessage(e))))
  })
  parameter_ui <- tagList(.tab_ms1(), .tab_ms2(), .tab_grouping(), .tab_annotation())
  parameter_tags <- c(Filter(function(tag) identical(tag$attribs$type, "number"),
                             htmltools::tagQuery(parameter_ui)$find("input")$selectedTags()),
                      htmltools::tagQuery(parameter_ui)$find("select")$selectedTags())

  validation_state <- reactive({
    client <- input$upload_validation_client %||% list()
    errors <- setNames(rep(list(character()), length(upload_fields) + 1L), c(upload_fields, "job"))
    pending <- missing <- character()
    current <- setNames(rep(FALSE, length(upload_fields)), upload_fields)
    required <- c("mzml_files", "metadata_file",
                  if (identical(input$lipid_db_source, "custom_excel")) "lipid_db_config_file")
    total <- 0
    labels <- c(mzml_files = "mzML files", metadata_file = "metadata file",
                lipid_db_config_file = "custom Excel configuration")
    for (field in upload_fields) {
      files <- input[[field]]
      errors[[field]] <- as.character(unlist(c(client$errors[[field]], client$transfer_errors[[field]]), use.names = FALSE))
      selected <- client$selections[[field]]
      current[[field]] <- !isTRUE(client$pending[[field]]) &&
        (is.null(selected) || (is.data.frame(files) && nrow(files) == length(selected) &&
                                 identical(vapply(selected, function(x) as.character(x$name), ""), as.character(files$name)) &&
                                 identical(vapply(selected, function(x) as.numeric(x$size), 0), as.numeric(files$size))))
      if (!length(errors[[field]]) && !current[[field]] && length(selected)) {
        count <- length(selected)
        pending <- c(pending, if (field == "mzml_files") sprintf("Waiting for %d mzML files to finish uploading.", count)
                     else paste("Waiting for", labels[[field]], "to finish uploading."))
      } else if (current[[field]] && !length(errors[[field]])) {
        errors[[field]] <- file_checks[[field]]()
      }
      if (field %in% required && (is.null(files) || (is.data.frame(files) && !nrow(files)) || (!is.null(selected) && !length(selected))) &&
          !length(selected) && !length(errors[[field]])) {
        missing <- c(missing, paste("Select", labels[[field]], "in Data import."))
      }
      if (current[[field]] && is.data.frame(files)) total <- total + sum(file.info(files$datapath)$size, na.rm = TRUE)
    }
    errors$job <- as.character(unlist(client$errors$job, use.names = FALSE))
    if (!length(errors$job) && total > 25000000000) {
      errors$job <- "The combined mzML and auxiliary files exceed 25 GB (decimal) per job."
    }
    md <- matched <- NULL
    if (current[["metadata_file"]] && !length(errors$metadata_file)) {
      result <- metadata_result()
      if (!is.null(result)) {
        md <- result$data
        errors$metadata_file <- result$error
      }
    }
    if (!is.null(md) && current[["mzml_files"]] && !length(errors$mzml_files) && !is.null(input$mzml_files)) {
      matched <- tryCatch(.match_metadata_and_files(md, input$mzml_files), error = function(e) {
        errors$job <<- c(errors$job, paste("Sample matching failed:", conditionMessage(e)))
        NULL
      })
    }
    file_errors <- unique(unlist(errors, use.names = FALSE))
    files_ready <- !length(c(file_errors, pending, missing)) && !is.null(matched)
    settings <- settings_errors <- character()
    for (tag in parameter_tags) {
      id <- tag$attribs$id
      value <- input[[id]]
      if (is.null(value) || !length(value) || anyNA(value) ||
          (identical(tag$attribs$type, "number") && (!is.numeric(value) || any(!is.finite(value))))) {
        settings_errors <- c(settings_errors, paste("Configure required parameter", id, "in the analysis tabs."))
      }
    }
    errors$settings <- settings_errors
    settings <- c(settings, settings_errors)
    reasons <- unique(c(file_errors, pending, missing, settings))
    list(errors = errors, pending = pending, missing = missing, metadata = md,
         matched = matched, files_ready = files_ready, ready = files_ready && !length(settings),
         reasons = reasons, revision = client$revision %||% 0L,
         status = if (length(c(file_errors, settings_errors))) "error" else if (length(pending)) "uploading"
         else if (files_ready) "complete" else "waiting")
  })
  validation_message <- function(field) {
    errors <- validation_state()$errors[[field]]
    if (!length(errors)) return(NULL)
    tags$div(role = "alert", style = "color:#b00020; font-weight:600;", lapply(errors, tags$p))
  }
  output$mzml_validation_message <- renderUI(validation_message("mzml_files"))
  output$metadata_validation_message <- renderUI(validation_message("metadata_file"))
  output$lipid_db_validation_message <- renderUI(validation_message("lipid_db_config_file"))
  output$upload_validation_message <- renderUI(validation_message("job"))
  output$run_validation_message <- renderUI({
    state <- validation_state()
    if (state$ready) return(NULL)
    tags$div(role = "status", class = "text-danger mt-2",
             tags$p("Before running the job:"), tags$ul(lapply(state$reasons, tags$li)),
             tags$p("Correct file selections and metadata in Data import."))
  })
  observe({
    state <- validation_state()
    session$sendCustomMessage("lipidms-upload-state", list(
      revision = state$revision, errors = state$errors, pending = state$pending,
      missing = state$missing, ready = state$ready, status = state$status))
  })
  metadata_raw <- reactive(validation_state()$metadata)
  matched_preview <- reactive(validation_state()$matched)

  output$import_summary <- renderTable({
    state <- validation_state()
    n_files <- if (is.null(input$mzml_files)) 0L else nrow(input$mzml_files)
    md_status <- if (length(state$pending)) "Uploading / validating selected files" else if (is.null(input$metadata_file)) "Not loaded" else {
      md <- metadata_raw()
      if (is.null(md)) {
        return(data.frame(message = "Metadata could not be read. Check column names and separator."))
      }
      paste0(nrow(md), " rows; ", ncol(md), " columns")
    }
    data.frame(
      item = c("Input mode", "mzML files", "Metadata", "Metadata status", "Polarity"),
      value = c("mzML files + metadata", n_files, input$metadata_file$name %||% NA_character_, md_status, input$to_process_polarity)
    )
  })

  output$matching_preview <- renderTable({
    if (length(validation_state()$pending)) {
      return(data.frame(message = "Waiting for uploaded files to finish validation."))
    }
    if (is.null(input$mzml_files) ||
        is.null(input$metadata_file)) {
      return(data.frame(
        message = "Upload mzML files and metadata to preview matching."
      ))
    }

    m <- matched_preview()

    if (is.null(m)) {
      return(data.frame(
        message = "Metadata/files could not be matched. Check metadata columns, separator and file names."
      ))
    }

    sampletype <- if ("sampletype" %in% names(m$metadata)) {
      as.character(m$metadata$sampletype)
    } else {
      rep(NA_character_, nrow(m$metadata))
    }
    data.frame(
      sample = m$metadata$sample,
      acquisition_mode = m$metadata$acquisition_mode,
      sample_type = m$metadata$sample_type,
      stringsAsFactors = FALSE
    )
  })

  # Use the default template bundled with the app.
  output$download_default_lipid_db_config <- downloadHandler(
    filename = function() {
      "LipidMS_default_db_config.xlsx"
    },
    content = function(file) {

      template <- system.file(
        "extdata",
        "LipidMS_default_db_config.xlsx",
        package = "LipidMS"
      )

      if (!nzchar(template) || !file.exists(template)) {
        stop("Default lipid database configuration template not found.")
      }

      if (!file.copy(template, file, overwrite = TRUE)) {
        stop("Could not copy the default lipid database configuration template.")
      }
    }
  )

  output$parameters_output <- renderUI({
    p <- .collect_lipidms_params(input)
    HTML(paste0(
      "<ul>",
      "<li><b>Input:</b> mzML files + metadata</li>",
      "<li><b>Job name:</b> ", p$jobname, "</li>",
      "<li><b>Polarity:</b> ", p$polarity, "</li>",
      "<li><b>MS1 ppm:</b> ", p$cwp_ms1_ppm, "</li>",
      "<li><b>MS2 ppm:</b> ", p$cwp_ms2_ppm, "</li>",
      "<li><b>dmzprecursor / dmzproducts / rttol:</b> ", p$dmzprecursor, " / ", p$dmzproducts, " / ", p$rttol, "</li>",
      "<li><b>remove_isotopes:</b> ", p$remove_isotopes, "</li>",
      "</ul>"
    ))
  })

  output$download_parameters_run <- downloadHandler(
    filename = function() paste0(.sanitize_job_name(input$jobname), "_selected_parameters.html"),
    content = function(file) {
      p <- .collect_lipidms_params(input)
      html <- paste0("<html><body><h3>LipidMS selected parameters</h3><pre>",
                     paste(capture.output(str(p)), collapse = "\n"),
                     "</pre></body></html>")
      writeLines(html, con = file)
    }
  )

  job_info <- reactiveVal(NULL)

  output$job_name_length_warning <- renderUI({
    req(input$jobname)
    count <- nchar(input$jobname, type = "chars")
    if (count <= .max_job_name_chars) return(NULL)
    tags$p(
      .job_name_length_message(count),
      role = "alert",
      style = "color:#b00020; font-weight:600; margin-top:-8px;"
    )
  })


  #==========================================================================#
  # Lipid DB and classes
  #==========================================================================#

  dbsReactive <- reactive({

    if (identical(input$lipid_db_source, "custom_excel")) {

      req(input$lipid_db_config_file)

      createLipidDbsFromFile(
        input$lipid_db_config_file$datapath
      )

    } else {

      createLipidDbs()

    }
  })

  #==========================================================================#
  # Run analysis
  #==========================================================================#

  xReactive <- reactiveVal(NULL)
  resultsReactive <- reactiveVal(NULL)

  observeEvent(input$Run_analysis, {
    req(input$jobname)
    state <- validation_state()
    if (!isTRUE(state$ready)) {
      showNotification("Check uploaded files and analysis parameters before running.",
                       type = "error")
      return(invisible(NULL))
    }

    run_analysis <- function() {

      tryCatch({

        matched <- req(state$matched)
        dbs <- req(dbsReactive())
        p <- .collect_lipidms_params(input)
        if (nchar(input$jobname) > .max_job_name_chars) {
          stop(.job_name_length_message(nchar(input$jobname)))
        }
        if (p$cwp_ms1_peakwidth_min >= p$cwp_ms1_peakwidth_max ||
            p$cwp_ms2_peakwidth_min >= p$cwp_ms2_peakwidth_max) {
          stop("Peakwidth minimum must be smaller than maximum.")
        }

        # para evitar que window entre en suspension
        sleep_handle <- prevent_sleep()
        on.exit(
          allow_sleep(sleep_handle),
          add = TRUE
        )

        withProgress(message = "Running LipidMS analysis", value = 0, {
          incProgress(0.0, detail = "Preparing input files")
          metadata <- matched$metadata
          files <- as.character(matched$files$datapath)
          rownames(metadata) <- make.unique(basename(files))
          metadata$sample <- basename(files)

          bpparam <- make_bpparam(input$parallel_mode, input$parallel_workers)
          ms1param <- xcms::CentWaveParam(
            ppm = p$cwp_ms1_ppm,
            peakwidth = c(p$cwp_ms1_peakwidth_min, p$cwp_ms1_peakwidth_max),
            snthresh = p$cwp_ms1_snthresh,
            prefilter = c(p$cwp_ms1_prefilter_n, p$cwp_ms1_prefilter_intensity),
            mzCenterFun = p$cwp_ms1_mzCenterFun,
            integrate = p$cwp_ms1_integrate
          )
          ms2param <- xcms::CentWaveParam(
            ppm = p$cwp_ms2_ppm,
            peakwidth = c(p$cwp_ms2_peakwidth_min, p$cwp_ms2_peakwidth_max),
            snthresh = p$cwp_ms2_snthresh,
            prefilter = c(p$cwp_ms2_prefilter_n, p$cwp_ms2_prefilter_intensity),
            mzCenterFun = p$cwp_ms2_mzCenterFun,
            integrate = p$cwp_ms2_integrate,
            mzdiff = p$cwp_ms2_mzdiff
          )
          grouping1 <- xcms::PeakDensityParam(
            sampleGroups = metadata$sample_type,
            bw = p$grouping1_bw,
            minFraction = p$grouping1_minFraction,
            binSize = p$grouping1_binSize
          )
          grouping2 <- xcms::PeakDensityParam(
            sampleGroups = metadata$sample_type,
            bw = p$grouping2_bw,
            minFraction = p$grouping2_minFraction,
            binSize = p$grouping2_binSize
          )
          rtparam <- xcms::PeakGroupsParam(
            minFraction = p$paramRT_minFraction,
            extraPeaks = p$paramRT_extraPeaks,
            smooth = p$paramRT_smooth,
            span = p$paramRT_span,
            family = p$paramRT_family,
            subsetAdjust = p$paramRT_subsetAdjust
          )

          paramannot <- defaultAnnotationParameters(p$polarity)
          for (nm in intersect(names(paramannot), names(p))) {
            paramannot[[nm]] <- p[[nm]]
          }
          paramannot$BPPARAM_isotopes <- bpparam

          paramdia <- defaultDIAExtractionParameters(
            findChromPeaksIsolationWindowParam = ms2param,
            min_ms2_intensity = p$min_ms2_intensity
          )
          paramdia$expandRt <- p$paramdia_expandRt
          paramdia$minCor <- p$paramdia_minCor
          paramdia$diffRt <- p$paramdia_diffRt
          paramdia$intensity <- p$paramdia_intensity
          paramdia$BPPARAM <- bpparam

          mse <- MsExperiment::readMsExperiment(
            spectraFiles = files,
            sampleData = S4Vectors::DataFrame(metadata)
          )
          incProgress(0.05, detail = "Finding MS1 chromatographic peaks")
          x <- xcms::findChromPeaks(mse, param = ms1param, BPPARAM = bpparam)
          incProgress(0.20, detail = "Initial grouping")
          x <- xcms::groupChromPeaks(x, param = grouping1)
          incProgress(0.05, detail = "Retention-time correction")
          x <- xcms::adjustRtime(x, param = rtparam, BPPARAM = bpparam)
          incProgress(0.05, detail = "Final grouping")
          x <- xcms::groupChromPeaks(x, param = grouping2)
          incProgress(0.10, detail = "Filling chromatographic peaks")
          x <- xcms::fillChromPeaks(x, param = xcms::ChromPeakAreaParam(), BPPARAM = bpparam)
          incProgress(0.05, detail = "Annotating lipids")
          x <- annotateLipids(
            x = x, paramannot = paramannot, paramdia = paramdia,
            dbs = dbs, rules = defaultFragmentRules(), BPPARAM = bpparam
          )
          incProgress(0.45, detail = "Preparing results")
          y <- cbind(
            as.data.frame(xcms::featureDefinitions(x)),
            as.data.frame(xcms::featureValues(x, value = "into", filled = TRUE))
          )
          if ("has_monoisotopic_peak" %in% names(y) &&
              isTRUE(p$only_confirmed_monoisotopic)) {
            y <- y[!is.na(y$has_monoisotopic_peak) &
                     y$has_monoisotopic_peak, , drop = FALSE]
          }
          y <- y[, !names(y) %in% c("peakidx", "monoisotopic_peak_ids"),
                 drop = FALSE]

          lipidms_metadata <- S4Vectors::metadata(x)$LipidMS
          xReactive(x)
          resultsReactive(list(
            x = x, feature_table = y,
            ms1_candidates = lipidms_metadata$ms1_candidates,
            ms2_evidence = lipidms_metadata$ms2_annotation$evidence
          ))
          job_info(list(name = p$jobname, completed = Sys.time()))
          incProgress(0.05, detail = "Done")
        })
        updateTabsetPanel(session, "inTabset", selected = "tab8")

      }, error = function(e) {
        message("LipidMS analysis failed: ", conditionMessage(e))
        showNotification(
          paste("Analysis failed:", conditionMessage(e)),
          type = "error",
          duration = NULL
        )
      })

    } # End run_analysis()

    run_analysis()

  }) # End observeEvent()

  #==========================================================================#
  # Results outputs
  #==========================================================================#

  output$featureTable_ui <- renderUI({
    req(resultsReactive())
    DT::dataTableOutput("featureTable")
  })

  output$featureTable <- DT::renderDataTable({
    req(resultsReactive())
    DT::datatable(
      resultsReactive()$feature_table,
      options = list(scrollX = TRUE, pageLength = 25)
    )
  })

  output$chromPeakTable_ui <- renderUI({
    req(xReactive())
    DT::dataTableOutput("chromPeakTable")
  })

  output$chromPeakTable <- DT::renderDataTable({
    req(xReactive())

    x <- xReactive()

    out <- cbind(
      as.data.frame(xcms::chromPeaks(x)),
      as.data.frame(xcms::chromPeakData(x))
    )

    DT::datatable(
      out,
      options = list(scrollX = TRUE, pageLength = 25)
    )
  })

  output$ms1Candidates_ui <- renderUI({
    req(resultsReactive())
    DT::dataTableOutput("ms1Candidates")
  })

  output$ms1Candidates <- DT::renderDataTable({
    req(resultsReactive())
    DT::datatable(
      resultsReactive()$ms1_candidates,
      extensions = "Buttons",
      options = list(dom = "Bfrtip", buttons = c("copy", "csv", "excel"),
                     scrollX = TRUE, pageLength = 25)
    )
  })

  output$ms2Evidence_ui <- renderUI({
    req(resultsReactive())
    DT::dataTableOutput("ms2Evidence")
  })

  output$ms2Evidence <- DT::renderDataTable({
    req(resultsReactive())
    DT::datatable(
      resultsReactive()$ms2_evidence,
      extensions = "Buttons",
      options = list(dom = "Bfrtip", buttons = c("copy", "csv", "excel"),
                     scrollX = TRUE, pageLength = 25)
    )
  })

  #==========================================================================#
  # Downloads
  #==========================================================================#

  output$download_featureTable_ui <- renderUI({
    req(resultsReactive())
    downloadButton("download_featureTable", "Feature table")
  })

  output$download_featureTable <- downloadHandler(
    filename = function() {
      paste0(input$jobname, "_features_", input$to_process_polarity, ".xlsx")
    },
    content = function(file) {
      req(resultsReactive())
      writexl::write_xlsx(resultsReactive()$feature_table, path = file)
    }
  )

  output$download_chromPeakTable_ui <- renderUI({
    req(xReactive())
    downloadButton("download_chromPeakTable", "Chrom peaks")
  })

  output$download_chromPeakTable <- downloadHandler(
    filename = function() {
      paste0(input$jobname, "_chromPeaks_", input$to_process_polarity, ".csv")
    },
    content = function(file) {
      req(xReactive())

      x <- xReactive()

      out <- cbind(
        as.data.frame(xcms::chromPeaks(x)),
        as.data.frame(xcms::chromPeakData(x))
      )

      utils::write.csv(out, file, row.names = TRUE)
    }
  )

  output$download_xcmsObject_ui <- renderUI({
    req(xReactive())
    downloadButton("download_xcmsObject", "XcmsExperiment object")
  })

  output$download_xcmsObject <- downloadHandler(
    filename = function() {
      paste0(input$jobname, "_annotation_", input$to_process_polarity, ".rds")
    },
    content = function(file) {
      req(xReactive())
      saveRDS(xReactive(), file = file)
    }
  )
}

shinyApp(ui, server)


#============================================================================#
# UI helpers
#============================================================================#

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


#============================================================================#
# User interface
#============================================================================#
options(shiny.maxRequestSize = 500 * 1024^3)

ui <- fluidPage(
  tagList(
    tags$head(
      tags$link(
        rel = "stylesheet",
        href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.5.0/css/all.min.css"
      ),
      tags$style(HTML("
        body {
          padding-top: 10px !important; font-family: 'Inter', sans-serif !important;
        }

        .shiny-notification {
          background-color: #fc804a;
          border: 1px solid #e74818;
          color: white;
          font-size: 16px;
          padding: 20px;
          max-width: 400px;
          height: auto;
          border-radius: 8px;
          margin-right: 20px;
          margin-left: auto;
          box-shadow: 0px 4px 10px rgba(0, 0, 0, 0.2);
        }

        .shiny-notification-message {
          background-color: #669973;
          border: 1px solid #4e7a5d;
        }

        input[type='checkbox'] {
          accent-color: #669999;
        }

        .checkbox {
          margin-left: 15px;
        }

        .radio input[type='radio'] {
          accent-color: #669999;
        }

        .radio {
          margin-left: 15px;
        }

        .radio label {
          display: inline-block;
          padding: 10px 20px;
          margin-right: 10px;
          border-radius: 6px;
          background-color: transparent;
          border: 1px solid transparent;
          cursor: pointer;
          transition: background-color 0.3s ease;
        }

        .radio label.selected-radio {
          background-color: #f8f8f8;
          border: none;
          box-shadow: none;
          width: 95%;
          font-weight: bold;
        }

        input[type='number'] {
          height: 35px;
          font-size: 11px;
        }

        .section-box {
          border: 1px solid #ddd;
          border-radius: 6px;
          padding: 15px;
          margin: 10px 10px;
          background-color: #f9f9f9;
        }

        td {
          max-width: 350px;
          white-space: nowrap;
          overflow: hidden;
          text-overflow: ellipsis;
        }

        td:hover {
          overflow: visible;
          white-space: normal;
          background-color: #f9f9f9;
        }
      ")),
      tags$script(HTML("
        $(document).on('shiny:inputchanged', function(event) {
          if (event.name === 'input_mode') {
            var $group = $('input[name=\"input_mode\"]').closest('.form-group');
            $group.find('label').removeClass('selected-radio');
            $group.find('input[name=\"input_mode\"]:checked').closest('label').addClass('selected-radio');
          }

          if (event.name === 'type_lipidclasses') {
            var $group = $('input[name=\"type_lipidclasses\"]').closest('.form-group');
            $group.find('label').removeClass('selected-radio');
            $group.find('input[name=\"type_lipidclasses\"]:checked').closest('label').addClass('selected-radio');
          }
        });

        $(document).ready(function() {
          $('input[name=\"input_mode\"]:checked').closest('label').addClass('selected-radio');
          $('input[name=\"type_lipidclasses\"]:checked').closest('label').addClass('selected-radio');
        });
      "))
    ),

    navbarPage(
      position = "static-top",
      theme = shinythemes::shinytheme("flatly"),
      title = " ",
      windowTitle = "LipidMS",
      id = "inTabset",

      #======================================================================#
      # Tab 1: Data import
      #======================================================================#

      tabPanel(
        title = "Data import",
        value = "tab1",

        div(
          style = "display: flex; align-items: center; margin: 0px 0px 20px 20px;",
          img(
            src = "iconohorizontal.png",
            style = "margin-right: 10px; height: 50px;"
          ),
          tags$h5(
            HTML("<strong style='color:#2C3E50; font-size:18px;'>Lipid annotation for LC-MS/MS data</strong>")
          )
        ),

        sidebarPanel(
          width = 3,

          fluidRow(
            radioButtons(
              "input_mode",
              "Input data",
              choices = c(
                "mzML files + metadata.csv" = "mzml",
                "XcmsExperiment object" = "xcms"
              ),
              selected = "mzml"
            ),
            style = "margin: 20px 0px 0px 0px"
          ),

          fluidRow(
            column(
              12,
              div(
                class = "section-box",
                radioButtons(
                  "to_process_polarity",
                  "Polarity",
                  choices = c(
                    "Positive" = "positive",
                    "Negative" = "negative"
                  ),
                  selected = "negative"
                )
              )
            )
          ),

          conditionalPanel(
            condition = "input.input_mode == 'mzml'",
            fluidRow(
              column(
                12,
                div(
                  class = "section-box",
                  fileInput(
                    "mzml_files",
                    "Choose mzML file/s",
                    multiple = TRUE,
                    accept = c(".mzML", ".mzXML")
                  ),
                  fileInput(
                    "metadata_file",
                    "Metadata file",
                    multiple = FALSE,
                    accept = c(
                      ".csv",
                      ".txt",
                      ".xlsx",
                      ".xls"
                    )
                  ),
                  selectInput(
                    "metadata_sep",
                    "Metadata delimiter (only for csv/txt)",
                    choices = c(
                      "Comma (,)" = ",",
                      "Semicolon (;)" = ";",
                      "Tab" = "\t"
                    ),
                    selected = ","
                  )
                )
              )
            )
          ),

          conditionalPanel(
            condition = "input.input_mode == 'xcms'",
            fluidRow(
              column(
                12,
                div(
                  class = "section-box",
                  fileInput(
                    "xcms_file",
                    "Choose XcmsExperiment file (.rds)",
                    multiple = FALSE,
                    accept = c(".rds")
                  )
                )
              )
            )
          ),

          tags$hr(style = "border-color: white"),

          fluidRow(
            column(
              12,
              tags$label(
                "Choose lipid class source for analysis",
                style = "font-size: 16px; font-weight: bold;"
              ),
              radioButtons(
                "type_lipidclasses",
                label = NULL,
                choices = list(
                  "Default list" = "lipidclasses_default",
                  "Custom list" = "lipidclasses_custom"
                ),
                selected = "lipidclasses_default",
                width = "100%"
              )
            )
          ),

          conditionalPanel(
            condition = "input.type_lipidclasses == 'lipidclasses_custom'",
            fluidRow(
              column(
                12,
                div(
                  class = "section-box",
                  fluidRow(
                    column(
                      6,
                      tags$label("Column delimiter", style = "font-weight: bold;"),
                      selectInput(
                        "classes_sep",
                        NULL,
                        choices = c(
                          "Semicolon (;)" = ";",
                          "Comma (,)" = ",",
                          "Pipe (|)" = "|"
                        ),
                        selected = ";"
                      )
                    ),
                    column(
                      6,
                      tags$label("Adducts delimiter", style = "font-weight: bold;"),
                      selectInput(
                        "adducts_sep",
                        NULL,
                        choices = c(
                          "Semicolon (;)" = ";",
                          "Comma (,)" = ",",
                          "Pipe (|)" = "|"
                        ),
                        selected = ","
                      )
                    )
                  ),
                  fileInput(
                    "classes_file",
                    "Choose custom lipid class file",
                    multiple = FALSE,
                    accept = c(".csv", ".txt")
                  )
                )
              )
            )
          ),

          tags$hr(style = "border-color: white"),

          fluidRow(
            actionButton("JumpTo2", "Next >", width = "100px"),
            style = "margin: 20px 0px 0px 0px"
          )
        ),

        mainPanel(
          tableOutput("import_summary")
        )
      ),


      #======================================================================#
      # Tab 2: MS1 peak picking
      #======================================================================#

      tabPanel(
        title = "MS1 peak-picking",
        value = "tab2",

        mainPanel(
          h4(strong("MS1 CentWaveParam")),

          .param_row(
            "ppm",
            "Mass accuracy for MS1 peak detection.",
            numericInput("cwp_ms1_ppm", NULL, value = 25, min = 0, max = 100, step = 1)
          ),

          .param_row_2(
            "peakwidth",
            "Expected chromatographic peak width range in seconds.",
            numericInput("cwp_ms1_peakwidth_min", "min", value = 5, min = 0, max = 500, step = 1),
            numericInput("cwp_ms1_peakwidth_max", "max", value = 50, min = 0, max = 500, step = 1)
          ),

          .param_row(
            "snthresh",
            "Signal-to-noise threshold.",
            numericInput("cwp_ms1_snthresh", NULL, value = 3, min = 0, max = 100, step = 1)
          ),

          .param_row_2(
            "prefilter",
            "Minimum number of scans and intensity threshold before peak detection.",
            numericInput("cwp_ms1_prefilter_n", "n", value = 5, min = 0, max = 100, step = 1),
            numericInput("cwp_ms1_prefilter_intensity", "intensity", value = 100, min = 0, max = 1e7, step = 10)
          ),

          .param_row(
            "mzCenterFun",
            "Function used to calculate the m/z center.",
            selectInput(
              "cwp_ms1_mzCenterFun",
              NULL,
              choices = c("wMean", "mean", "apex"),
              selected = "wMean"
            )
          ),

          .param_row(
            "integrate",
            "Integration method.",
            selectInput(
              "cwp_ms1_integrate",
              NULL,
              choices = c("1" = 1, "2" = 2),
              selected = 2
            )
          ),

          tags$hr(),

          fluidRow(
            actionButton("GoBackTo1", "< Previous", width = "100px", style = "margin:40px 0px"),
            actionButton("JumpTo3", "Next >", width = "100px", style = "margin:40px 0px")
          )
        )
      ),


      #======================================================================#
      # Tab 3: MS2 / DIA peak picking
      #======================================================================#

      tabPanel(
        title = "MS2 / DIA",
        value = "tab3",

        mainPanel(
          h4(strong("MS2 CentWaveParam")),

          .param_row(
            "ppm",
            "Mass accuracy for MS2/DIA peak detection.",
            numericInput("cwp_ms2_ppm", NULL, value = 25, min = 0, max = 100, step = 1)
          ),

          .param_row_2(
            "peakwidth",
            "Expected MS2 chromatographic peak width range in seconds.",
            numericInput("cwp_ms2_peakwidth_min", "min", value = 3, min = 0, max = 500, step = 1),
            numericInput("cwp_ms2_peakwidth_max", "max", value = 50, min = 0, max = 500, step = 1)
          ),

          .param_row(
            "snthresh",
            "Signal-to-noise threshold.",
            numericInput("cwp_ms2_snthresh", NULL, value = 2, min = 0, max = 100, step = 1)
          ),

          .param_row_2(
            "prefilter",
            "Minimum number of scans and intensity threshold before MS2 peak detection.",
            numericInput("cwp_ms2_prefilter_n", "n", value = 3, min = 0, max = 100, step = 1),
            numericInput("cwp_ms2_prefilter_intensity", "intensity", value = 50, min = 0, max = 1e7, step = 10)
          ),

          .param_row(
            "mzCenterFun",
            "Function used to calculate the m/z center.",
            selectInput(
              "cwp_ms2_mzCenterFun",
              NULL,
              choices = c("wMean", "mean", "apex"),
              selected = "wMean"
            )
          ),

          .param_row(
            "integrate",
            "Integration method.",
            selectInput(
              "cwp_ms2_integrate",
              NULL,
              choices = c("1" = 1, "2" = 2),
              selected = 2
            )
          ),

          .param_row(
            "mzdiff",
            "Minimum m/z difference for overlapping peaks.",
            numericInput("cwp_ms2_mzdiff", NULL, value = -0.001, min = -1, max = 1, step = 0.001)
          ),

          tags$hr(),

          h4(strong("DIA / SWATH reconstruction parameters")),

          .param_row(
            "expandRt",
            "Retention-time expansion around candidate peaks.",
            numericInput("paramdia_expandRt", NULL, value = 0, min = 0, max = 100, step = 1)
          ),

          .param_row(
            "minCor",
            "Minimum chromatographic correlation for DIA fragment association.",
            numericInput("paramdia_minCor", NULL, value = 0.7, min = 0, max = 1, step = 0.05)
          ),

          .param_row(
            "diffRt",
            "Maximum RT difference for DIA fragment association.",
            numericInput("paramdia_diffRt", NULL, value = 3, min = 0, max = 100, step = 1)
          ),

          .param_row(
            "intensity",
            "Intensity column used for DIA reconstruction.",
            selectInput(
              "paramdia_intensity",
              NULL,
              choices = c("into", "maxo"),
              selected = "into"
            )
          ),

          tags$hr(),

          fluidRow(
            actionButton("GoBackTo2", "< Previous", width = "100px", style = "margin:40px 0px"),
            actionButton("JumpTo4", "Next >", width = "100px", style = "margin:40px 0px")
          )
        )
      ),


      #======================================================================#
      # Tab 4: Grouping + RT correction
      #======================================================================#

      tabPanel(
        title = "Grouping / RT",
        value = "tab4",

        mainPanel(
          h4(strong("Initial grouping")),

          .param_row(
            "bw",
            "Bandwidth for chromatographic peak grouping.",
            numericInput("grouping1_bw", NULL, value = 5, min = 0, max = 100, step = 1)
          ),

          .param_row(
            "minFraction",
            "Minimum fraction of samples required in initial grouping.",
            numericInput("grouping1_minFraction", NULL, value = 0.9, min = 0, max = 1, step = 0.05)
          ),

          .param_row(
            "binSize",
            "m/z bin size for initial grouping.",
            numericInput("grouping1_binSize", NULL, value = 0.01, min = 0, max = 1, step = 0.001)
          ),

          tags$hr(),

          h4(strong("Retention-time correction: PeakGroupsParam")),

          .param_row(
            "minFraction",
            "Minimum fraction of samples for peak groups used in RT correction.",
            numericInput("paramRT_minFraction", NULL, value = 0.9, min = 0, max = 1, step = 0.05)
          ),

          .param_row(
            "extraPeaks",
            "Number of additional peaks allowed in peak groups.",
            numericInput("paramRT_extraPeaks", NULL, value = 1, min = 0, max = 20, step = 1)
          ),

          .param_row(
            "smooth",
            "Smoothing method for RT correction.",
            selectInput(
              "paramRT_smooth",
              NULL,
              choices = c("loess", "linear"),
              selected = "loess"
            )
          ),

          .param_row(
            "span",
            "Span parameter for loess smoothing.",
            numericInput("paramRT_span", NULL, value = 0.4, min = 0, max = 1, step = 0.05)
          ),

          .param_row(
            "family",
            "Loess family.",
            selectInput(
              "paramRT_family",
              NULL,
              choices = c("gaussian", "symmetric"),
              selected = "gaussian"
            )
          ),

          .param_row(
            "subsetAdjust",
            "How subset correction is applied.",
            selectInput(
              "paramRT_subsetAdjust",
              NULL,
              choices = c("average", "previous"),
              selected = "average"
            )
          ),

          tags$hr(),

          h4(strong("Final grouping")),

          .param_row(
            "bw",
            "Bandwidth for final chromatographic peak grouping.",
            numericInput("grouping2_bw", NULL, value = 5, min = 0, max = 100, step = 1)
          ),

          .param_row(
            "minFraction",
            "Minimum fraction of samples required in final grouping.",
            numericInput("grouping2_minFraction", NULL, value = 0.5, min = 0, max = 1, step = 0.05)
          ),

          .param_row(
            "binSize",
            "m/z bin size for final grouping.",
            numericInput("grouping2_binSize", NULL, value = 0.01, min = 0, max = 1, step = 0.001)
          ),

          tags$hr(),

          fluidRow(
            actionButton("GoBackTo3", "< Previous", width = "100px", style = "margin:40px 0px"),
            actionButton("JumpTo5", "Next >", width = "100px", style = "margin:40px 0px")
          )
        )
      ),


      #======================================================================#
      # Tab 5: Annotation
      #======================================================================#

      tabPanel(
        title = "Annotation",
        value = "tab5",

        mainPanel(
          h4(strong("Precursor and product matching")),

          .param_row(
            "dmzprecursor",
            "Mass tolerance for precursor ions, in ppm.",
            numericInput("dmzprecursor", NULL, value = 10, min = 0, max = 100, step = 1)
          ),

          .param_row(
            "dmzproducts",
            "Mass tolerance for product ions, in ppm.",
            numericInput("dmzproducts", NULL, value = 15, min = 0, max = 100, step = 1)
          ),

          .param_row(
            "rttol",
            "Retention-time tolerance for linking precursor and product ions.",
            numericInput("rttol", NULL, value = 3, min = 0, max = 100, step = 1)
          ),

          .param_row(
            "min_ms2_intensity",
            "Minimum MS2 fragment intensity used during annotation.",
            numericInput("min_ms2_intensity", NULL, value = 100, min = 0, max = 1e7, step = 10)
          ),

          tags$hr(),

          h4(strong("Isotope annotation")),

          .param_row(
            "isotope_mzabs",
            "Absolute m/z tolerance for isotope annotation.",
            numericInput("isotope_mzabs", NULL, value = 0.005, min = 0, max = 1, step = 0.001)
          ),

          .param_row(
            "isotope_maxcharge",
            "Maximum isotope charge.",
            numericInput("isotope_maxcharge", NULL, value = 2, min = 1, max = 10, step = 1)
          ),

          .param_row(
            "isotope_maxiso",
            "Maximum isotope rank.",
            numericInput("isotope_maxiso", NULL, value = 5, min = 1, max = 20, step = 1)
          ),

          .param_row(
            "isotope_minfrac",
            "Minimum fraction for isotope annotation.",
            numericInput("isotope_minfrac", NULL, value = 0.5, min = 0, max = 1, step = 0.05)
          ),

          .param_row(
            "isotope_rttol",
            "Retention-time tolerance for isotope annotation.",
            numericInput("isotope_rttol", NULL, value = 3, min = 0, max = 100, step = 1)
          ),

          .param_row(
            "intval",
            "Intensity value used for isotope annotation.",
            selectInput(
              "intval",
              NULL,
              choices = c("into", "maxo"),
              selected = "into"
            )
          ),

          tags$hr(),

          h4(strong("Annotation filters")),

          fluidRow(
            column(
              12,
              checkboxInput(
                "only_confirmed_monoisotopic",
                "Only keep confirmed monoisotopic features",
                value = TRUE
              ),
              checkboxInput(
                "remove_isotopes",
                "Remove isotope peaks before MS1 candidate search",
                value = TRUE
              )
            )
          ),

          tags$hr(),

          fluidRow(
            column(
              12,
              checkboxGroupInput(
                "lipidclasses_annotate",
                "Lipid classes to annotate:",
                choices = NULL,
                selected = NULL,
                inline = TRUE
              )
            )
          ),

          tags$hr(),

          fluidRow(
            actionButton("GoBackTo4", "< Previous", width = "100px", style = "margin:40px 0px"),
            actionButton("JumpTo6", "Next >", width = "100px", style = "margin:40px 0px")
          )
        )
      ),


      #======================================================================#
      # Tab 6: Run
      #======================================================================#

      tabPanel(
        title = "Run",
        value = "tab6",

        div(
          style = "display: flex; align-items: center; margin: 0px 0px 20px 20px;",
          img(
            src = "iconohorizontal.png",
            style = "margin-right: 10px; height: 50px;"
          ),
          tags$h5(
            HTML("<strong style='color:#2C3E50; font-size:18px;'>Lipid annotation for LC-MS/MS data</strong>")
          )
        ),

        sidebarPanel(
          width = 3,

          textInput(
            "jobname",
            "Job Name",
            value = paste("Job", as.character(Sys.Date()), sep = "_")
          ),

          fluidRow(
            column(
              9,
              textInput(
                "resdir",
                "Results folder",
                value = "./results"
              )
            ),
            column(
              3,
              br(),
              shinyFiles::shinyFilesButton(
                "resdir_button",
                "Browse",
                "Select results folder",
                multiple = FALSE
              )
            )
          ),

          numericInput(
            "ncores",
            "Number of cores",
            value = 4,
            min = 1,
            max = 32,
            step = 1
          ),

          selectInput(
            "parallel_backend",
            "Parallel backend",
            choices = c(
              "SerialParam" = "SerialParam",
              "SnowParam" = "SnowParam",
              "MulticoreParam" = "MulticoreParam"
            ),
            selected = "SerialParam"
          ),

          tags$hr(style = "border-color: white"),

          actionButton(
            "Run_analysis",
            "Run analysis",
            style = "color: #fff; background-color: #669999; border-color: #000; width: 100%;"
          ),

          tags$hr(style = "border-color: white"),

          fluidRow(
            actionButton("GoBackTo5", "< Previous", style = "margin:10px 10px"),
            actionButton("JumpTo7", "Next >", style = "margin:10px 10px")
          )
        ),

        column(
          9,
          tags$div(
            style = "padding: 0px; height: 75vh; overflow-y: auto; position: relative;",
            downloadButton(
              "download_parameters_run",
              "Selected Parameters",
              style = "position: sticky; top: 10px; right: 10px; z-index: 1000; float: right;"
            ),
            uiOutput("parameters_output")
          )
        )
      ),


      #======================================================================#
      # Tab 7: Results
      #======================================================================#

      tabPanel(
        title = "Results",
        value = "tab7",

        fluidRow(
          column(
            width = 2,
            div(
              style = "display: flex; align-items: center; margin: 0px 0px 20px 20px;",
              img(
                src = "iconohorizontal.png",
                style = "margin-right: 10px; height: 50px;"
              )
            ),

            sidebarPanel(
              width = 12,

              tags$h5(
                HTML("<strong style='color:#2C3E50; font-size:15px;'>Explore Results</strong>")
              ),

              tabsetPanel(
                id = "tabs",
                type = "tabs",
                tabPanel("Feature table", value = "features"),
                tabPanel("Chrom peaks", value = "chrompeaks"),
                tabPanel("MS1 candidates", value = "ms1candidates"),
                tabPanel("MS2 evidence", value = "ms2evidence")
              ),

              tags$hr(style = "border-color: white"),

              tags$h5(
                HTML("<strong style='color:#2C3E50; font-size:15px;'>Download results</strong>")
              ),

              fluidRow(
                tags$div(uiOutput("download_featureTable_ui"), style = "margin-bottom: 3px;"),
                tags$div(uiOutput("download_chromPeakTable_ui"), style = "margin-bottom: 3px;"),
                tags$div(uiOutput("download_xcmsObject_ui"), style = "margin-bottom: 3px;")
              ),

              tags$hr(style = "border-color: white"),

              fluidRow(
                actionButton(
                  "GoBackTo6",
                  "< Previous",
                  width = "100px",
                  style = "margin:40px 0px"
                )
              )
            )
          ),

          column(
            width = 10,

            conditionalPanel(
              condition = "input.tabs == 'features'",
              uiOutput("featureTable_ui")
            ),

            conditionalPanel(
              condition = "input.tabs == 'chrompeaks'",
              uiOutput("chromPeakTable_ui")
            ),

            conditionalPanel(
              condition = "input.tabs == 'ms1candidates'",
              uiOutput("ms1Candidates_ui")
            ),

            conditionalPanel(
              condition = "input.tabs == 'ms2evidence'",
              uiOutput("ms2Evidence_ui")
            )
          )
        )
      )
    )
  )
)

#============================================================================#
# Shiny application server
#============================================================================#
##############################################################################
server <- function(input, output, session) {

  session$onSessionEnded(function() {
    stopApp()
  })

  #==========================================================================#
  # Reactive values
  #==========================================================================#

  xReactive <- reactiveVal(NULL)
  resultsReactive <- reactiveVal(NULL)
  dbsReactive <- reactiveVal(NULL)
  rulesReactive <- reactiveVal(NULL)

  #==========================================================================#
  # Tab visibility depending on input mode
  #==========================================================================#

  observe({
    req(input$input_mode)

    if (input$input_mode == "xcms") {
      hideTab("inTabset", "tab2")
      hideTab("inTabset", "tab3")
      hideTab("inTabset", "tab4")

      showTab("inTabset", "tab5")
      showTab("inTabset", "tab6")
    } else {
      showTab("inTabset", "tab2")
      showTab("inTabset", "tab3")
      showTab("inTabset", "tab4")
      showTab("inTabset", "tab5")
      showTab("inTabset", "tab6")
    }

    if (is.null(resultsReactive())) {
      hideTab("inTabset", "tab7")
    } else {
      showTab("inTabset", "tab7")
    }
  })

  #==========================================================================#
  # Navigation
  #==========================================================================#

  observeEvent(input$JumpTo2, {
    if (input$input_mode == "xcms") {
      updateTabsetPanel(session, "inTabset", selected = "tab5")
    } else {
      updateTabsetPanel(session, "inTabset", selected = "tab2")
    }
  })

  observeEvent(input$GoBackTo1, {
    updateTabsetPanel(session, "inTabset", selected = "tab1")
  })

  observeEvent(input$JumpTo3, {
    updateTabsetPanel(session, "inTabset", selected = "tab3")
  })

  observeEvent(input$GoBackTo2, {
    updateTabsetPanel(session, "inTabset", selected = "tab2")
  })

  observeEvent(input$JumpTo4, {
    updateTabsetPanel(session, "inTabset", selected = "tab4")
  })

  observeEvent(input$GoBackTo3, {
    updateTabsetPanel(session, "inTabset", selected = "tab3")
  })

  observeEvent(input$JumpTo5, {
    updateTabsetPanel(session, "inTabset", selected = "tab5")
  })

  observeEvent(input$GoBackTo4, {
    if (input$input_mode == "xcms") {
      updateTabsetPanel(session, "inTabset", selected = "tab1")
    } else {
      updateTabsetPanel(session, "inTabset", selected = "tab4")
    }
  })

  observeEvent(input$JumpTo6, {
    updateTabsetPanel(session, "inTabset", selected = "tab6")
  })

  observeEvent(input$GoBackTo5, {
    updateTabsetPanel(session, "inTabset", selected = "tab5")
  })

  observeEvent(input$JumpTo7, {
    updateTabsetPanel(session, "inTabset", selected = "tab7")
  })

  observeEvent(input$GoBackTo6, {
    updateTabsetPanel(session, "inTabset", selected = "tab6")
  })

  #==========================================================================#
  # Parallel backend
  #==========================================================================#

  make_bpparam <- reactive({
    ncores <- as.integer(input$ncores)

    if (input$parallel_backend == "SnowParam") {
      BiocParallel::SnowParam(ncores)
    } else if (input$parallel_backend == "MulticoreParam") {
      BiocParallel::MulticoreParam(ncores)
    } else {
      BiocParallel::SerialParam()
    }
  })

  #==========================================================================#
  # Parameter builders
  #==========================================================================#

  cwp_ms1_reactive <- reactive({
    xcms::CentWaveParam(
      ppm = input$cwp_ms1_ppm,
      peakwidth = c(input$cwp_ms1_peakwidth_min,
                    input$cwp_ms1_peakwidth_max),
      snthresh = input$cwp_ms1_snthresh,
      prefilter = c(input$cwp_ms1_prefilter_n,
                    input$cwp_ms1_prefilter_intensity),
      mzCenterFun = input$cwp_ms1_mzCenterFun,
      integrate = as.integer(input$cwp_ms1_integrate)
    )
  })

  cwp_ms2_reactive <- reactive({
    xcms::CentWaveParam(
      ppm = input$cwp_ms2_ppm,
      peakwidth = c(input$cwp_ms2_peakwidth_min,
                    input$cwp_ms2_peakwidth_max),
      snthresh = input$cwp_ms2_snthresh,
      prefilter = c(input$cwp_ms2_prefilter_n,
                    input$cwp_ms2_prefilter_intensity),
      mzCenterFun = input$cwp_ms2_mzCenterFun,
      integrate = as.integer(input$cwp_ms2_integrate),
      mzdiff = input$cwp_ms2_mzdiff
    )
  })

  paramannot_reactive <- reactive({
    list(
      polarity = input$to_process_polarity,

      dmzprecursor = input$dmzprecursor,
      dmzproducts = input$dmzproducts,
      rttol = input$rttol,

      isotope_mzabs = input$isotope_mzabs,
      isotope_maxcharge = as.integer(input$isotope_maxcharge),
      isotope_maxiso = as.integer(input$isotope_maxiso),
      isotope_minfrac = input$isotope_minfrac,
      isotope_rttol = input$isotope_rttol,
      intval = input$intval,

      min_ms2_intensity = input$min_ms2_intensity,
      only_confirmed_monoisotopic = input$only_confirmed_monoisotopic,
      remove_isotopes = input$remove_isotopes
    )
  })

  paramdia_reactive <- reactive({
    list(
      findChromPeaksIsolationWindowParam = cwp_ms2_reactive(),
      expandRt = input$paramdia_expandRt,
      minCor = input$paramdia_minCor,
      diffRt = input$paramdia_diffRt,
      intensity = input$paramdia_intensity,
      minint = input$min_ms2_intensity
    )
  })

  paramRT_reactive <- reactive({
    xcms::PeakGroupsParam(
      minFraction = input$paramRT_minFraction,
      extraPeaks = input$paramRT_extraPeaks,
      smooth = input$paramRT_smooth,
      span = input$paramRT_span,
      family = input$paramRT_family,
      peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
      subset = integer(),
      subsetAdjust = input$paramRT_subsetAdjust
    )
  })

  paramGrouping_initial_reactive <- reactive({
    req(files_reactive())

    xcms::PeakDensityParam(
      sampleGroups = rep(1, length(files_reactive())),
      bw = input$grouping1_bw,
      minFraction = input$grouping1_minFraction,
      binSize = input$grouping1_binSize
    )
  })

  paramGrouping_final_reactive <- reactive({
    req(files_reactive())

    xcms::PeakDensityParam(
      sampleGroups = rep(1, length(files_reactive())),
      bw = input$grouping2_bw,
      minFraction = input$grouping2_minFraction,
      binSize = input$grouping2_binSize
    )
  })

  #==========================================================================#
  # Metadata and files
  #==========================================================================#

  metadata_raw_reactive <- reactive({

    req(input$metadata_file)

    ext <- tolower(tools::file_ext(input$metadata_file$name))

    metadata <- tryCatch({

      switch(
        ext,

        csv = utils::read.csv(
          input$metadata_file$datapath,
          sep = input$metadata_sep,
          stringsAsFactors = FALSE,
          check.names = FALSE
        ),

        txt = utils::read.table(
          input$metadata_file$datapath,
          sep = input$metadata_sep,
          header = TRUE,
          stringsAsFactors = FALSE,
          check.names = FALSE
        ),

        xlsx = as.data.frame(
          readxl::read_excel(input$metadata_file$datapath),
          stringsAsFactors = FALSE
        ),

        xls = as.data.frame(
          readxl::read_excel(input$metadata_file$datapath),
          stringsAsFactors = FALSE
        ),

        stop("Unsupported metadata format. Use csv, txt, xlsx or xls.")
      )

    }, error = function(e) {
      validate(
        need(
          FALSE,
          paste0(
            "Metadata file could not be read. ",
            "Please check the file format and delimiter. Details: ",
            conditionMessage(e)
          )
        )
      )
    })

    validate(
      need(
        is.data.frame(metadata),
        "Metadata could not be read as a table. Please check the delimiter."
      )
    )

    if (ncol(metadata) == 1L && ext %in% c("csv", "txt")) {
      validate(
        need(
          FALSE,
          paste0(
            "Only one column was detected in the metadata file. ",
            "The delimiter is probably incorrect. ",
            "Please check whether the file uses comma (,), semicolon (;) or tab."
          )
        )
      )
    }

    required_cols <- c("sample", "acquisition_mode")

    validate(
      need(
        all(required_cols %in% colnames(metadata)),
        paste0(
          "Metadata column names are not valid. ",
          "The metadata file must contain exactly these required columns: ",
          paste(required_cols, collapse = ", "),
          ". Columns detected: ",
          paste(colnames(metadata), collapse = ", "),
          ". Please rename the columns or check the delimiter."
        )
      )
    )

    metadata
  })

  metadata_reactive <- reactive({

    metadata <- metadata_raw_reactive()

    out <- S4Vectors::DataFrame(
      sample = basename(metadata$sample),
      acquisition_mode = metadata$acquisition_mode,
      row.names = metadata$sample
    )

    acquisition <- toupper(as.character(out$acquisition_mode))

    acquisition[grepl("FULLSCAN|FULLSCAN/MS|MS1|MS$", acquisition)] <- "fullscan"
    acquisition[grepl("DDA", acquisition)] <- "DDA"
    acquisition[grepl("SWATH|DIA", acquisition)] <- "DIA"

    validate(
      need(
        all(acquisition %in% c("FULLSCAN", "DDA", "DIA", "fullscan")),
        "Column 'acquisition_mode' must contain values compatible with fullscan/MS, DDA or DIA/SWATH."
      )
    )

    out$acquisition_mode <- acquisition

    out
  })

  files_reactive <- reactive({
    req(input$mzml_files)
    input$mzml_files$datapath
  })

  volumes <- c(
    Home = normalizePath("~"),
    getVolumes()()
  )

  shinyDirChoose(
    input,
    "resdir_button",
    roots = volumes
  )

  #==========================================================================#
  # Lipid DB and classes
  #==========================================================================#

  observe({
    dbs <- createLipidDbs()
    dbsReactive(dbs)

    if (!is.null(dbs$adducts)) {
      classes <- sort(unique(dbs$adducts$class))

      updateCheckboxGroupInput(
        session,
        "lipidclasses_annotate",
        choices = classes,
        selected = classes
      )
    }
  })

  #==========================================================================#
  # Results folder
  #==========================================================================#

  observe({

    req(input$resdir_button)

    path <- parseDirPath(
      volumes,
      input$resdir_button
    )

    if (length(path) > 0) {

      updateTextInput(
        session,
        "resdir",
        value = path
      )

    }

  })

  #==========================================================================#
  # Import summary
  #==========================================================================#

  output$import_summary <- renderTable({

    if (input$input_mode == "mzml") {

      metadata_status <- "Not loaded"

      if (!is.null(input$metadata_file)) {
        md <- metadata_raw_reactive()
        metadata_status <- paste0(
          nrow(md),
          " rows; columns: ",
          paste(colnames(md), collapse = ", ")
        )
      }

      data.frame(
        item = c(
          "Input mode",
          "mzML files",
          "Metadata file",
          "Metadata status",
          "Polarity"
        ),
        value = c(
          "mzML files + metadata",
          if (!is.null(input$mzml_files)) nrow(input$mzml_files) else 0,
          if (!is.null(input$metadata_file)) input$metadata_file$name else NA,
          metadata_status,
          input$to_process_polarity
        )
      )

    } else {

      data.frame(
        item = c(
          "Input mode",
          "XcmsExperiment file",
          "Polarity"
        ),
        value = c(
          "XcmsExperiment object (.rds)",
          if (!is.null(input$xcms_file)) input$xcms_file$name else NA,
          input$to_process_polarity
        )
      )
    }
  })

  #==========================================================================#
  # Selected parameters summary
  #==========================================================================#

  output$parameters_output <- renderUI({
    HTML(paste0(
      "<h4 style='color:#669999; font-weight:bold;'>",
      "<i class='fa-solid fa-circle-info' style='color:#669999; margin-right:8px;'></i>",
      "Selected Parameters</h4>",

      "<h5><b>Input</b></h5>",
      "<ul>",
      "<li><b>Input mode:</b> ", input$input_mode, "</li>",
      "<li><b>Polarity:</b> ", input$to_process_polarity, "</li>",
      "<li><b>Job name:</b> ", input$jobname, "</li>",
      "<li><b>Results folder:</b> ", input$resdir, "</li>",
      "</ul>",

      if (input$input_mode == "mzml") {
        paste0(
          "<h5><b>MS1 CentWaveParam</b></h5>",
          "<ul>",
          "<li>ppm: ", input$cwp_ms1_ppm, "</li>",
          "<li>peakwidth: c(", input$cwp_ms1_peakwidth_min, ", ",
          input$cwp_ms1_peakwidth_max, ")</li>",
          "<li>snthresh: ", input$cwp_ms1_snthresh, "</li>",
          "<li>prefilter: c(", input$cwp_ms1_prefilter_n, ", ",
          input$cwp_ms1_prefilter_intensity, ")</li>",
          "<li>mzCenterFun: ", input$cwp_ms1_mzCenterFun, "</li>",
          "<li>integrate: ", input$cwp_ms1_integrate, "</li>",
          "</ul>"
        )
      } else {
        ""
      },

      "<h5><b>Annotation</b></h5>",
      "<ul>",
      "<li>dmzprecursor: ", input$dmzprecursor, "</li>",
      "<li>dmzproducts: ", input$dmzproducts, "</li>",
      "<li>rttol: ", input$rttol, "</li>",
      "<li>min_ms2_intensity: ", input$min_ms2_intensity, "</li>",
      "<li>remove_isotopes: ", input$remove_isotopes, "</li>",
      "<li>only_confirmed_monoisotopic: ",
      input$only_confirmed_monoisotopic, "</li>",
      "</ul>"
    ))
  })

  output$download_parameters_run <- downloadHandler(
    filename = function() {
      paste0(input$jobname, "_selected_parameters.html")
    },
    content = function(file) {
      html <- as.character(output$parameters_output())
      writeLines(html, con = file)
    }
  )

  #==========================================================================#
  # Load XcmsExperiment object
  #==========================================================================#

  load_xcms_object <- function(path) {

    x <- readRDS(path)

    validate(
      need(
        inherits(x, "XcmsExperiment"),
        "The loaded object is not an XcmsExperiment."
      )
    )

    x
  }

  #==========================================================================#
  # Run analysis
  #==========================================================================#

  observeEvent(input$Run_analysis, {

    req(input$jobname)

    dir.create(input$resdir, showWarnings = FALSE, recursive = TRUE)

    withProgress(message = "Running LipidMS analysis", value = 0, {

      incProgress(0.05, detail = "Preparing lipid databases and rules")

      dbs <- createLipidDbs()
      rules <- defaultFragmentRules()

      dbsReactive(dbs)
      rulesReactive(rules)

      if (input$input_mode == "mzml") {

        req(input$mzml_files, input$metadata_file)

        incProgress(0.10, detail = "Reading metadata and mzML files")

        metadata <- metadata_reactive()
        files <- files_reactive()

        uploaded_names <- basename(input$mzml_files$name)
        metadata_names <- basename(as.character(metadata$sample))

        missing_in_files <- setdiff(metadata_names, uploaded_names)
        missing_in_metadata <- setdiff(uploaded_names, metadata_names)

        if (length(missing_in_files) > 0) {
          showNotification(
            paste0(
              "Samples present in metadata but not uploaded: ",
              paste(missing_in_files, collapse = ", "),
              ". They will be ignored."
            ),
            type = "warning",
            duration = NULL
          )
        }

        if (length(missing_in_metadata) > 0) {
          showNotification(
            paste0(
              "Uploaded mzML/mzXML files not present in metadata: ",
              paste(missing_in_metadata, collapse = ", "),
              ". They will be ignored."
            ),
            type = "warning",
            duration = NULL
          )
        }

        keep_metadata <- metadata_names %in% uploaded_names

        metadata <- metadata[keep_metadata, , drop = FALSE]
        metadata_names <- basename(as.character(metadata$sample))

        file_idx <- match(metadata_names, uploaded_names)

        validate(
          need(
            nrow(metadata) > 0,
            "No common samples were found between metadata and uploaded mzML/mzXML files."
          ),
          need(
            all(!is.na(file_idx)),
            "Some metadata samples could not be matched to uploaded mzML/mzXML files."
          )
        )

        files <- files[file_idx]

        rownames(metadata) <- basename(files)
        metadata$sample <- basename(files)

        validate(
          need(
            length(files) == nrow(metadata),
            paste0(
              "Internal matching error: ",
              length(files),
              " files but ",
              nrow(metadata),
              " metadata rows."
            )
          )
        )

        BiocParallel::register(make_bpparam())

        mse <- MsExperiment::readMsExperiment(
          spectraFiles = files,
          sampleData = metadata
        )

        incProgress(0.30, detail = "Finding MS1 chromatographic peaks")

        x <- xcms::findChromPeaks(
          mse,
          param = cwp_ms1_reactive()
        )

        incProgress(0.35, detail = "Initial grouping")

        x <- xcms::groupChromPeaks(
          x,
          param = paramGrouping_initial_reactive()
        )

        incProgress(0.40, detail = "Retention-time correction")

        x <- xcms::adjustRtime(
          x,
          param = paramRT_reactive(),
          BPPARAM = make_bpparam()
        )

        incProgress(0.50, detail = "Final grouping")

        x <- xcms::groupChromPeaks(
          x,
          param = paramGrouping_final_reactive()
        )

        incProgress(0.55, detail = "Filling chromatographic peaks")

        x <- xcms::fillChromPeaks(
          x,
          param = xcms::ChromPeakAreaParam()
        )

        # saveRDS(
        #   x,
        #   file = file.path(
        #     input$resdir,
        #     paste0(
        #       input$jobname,
        #       "_xcmsexperiment_",
        #       input$to_process_polarity,
        #       ".rds"
        #     )
        #   )
        # )

      } else {

        req(input$xcms_file)

        incProgress(0.55, detail = "Loading XcmsExperiment")

        x <- load_xcms_object(
          path = input$xcms_file$datapath
        )
      }

      incProgress(0.60, detail = "Annotating lipids")

      paramannot <- paramannot_reactive()
      paramdia <- paramdia_reactive()

      x <- annotateLipids(
        x = x,
        paramannot = paramannot,
        paramdia = paramdia,
        dbs = dbs,
        rules = rules,
        BPPARAM = make_bpparam()
      )

      incProgress(0.90, detail = "Preparing feature table")

      y <- cbind(
        as.data.frame(xcms::featureDefinitions(x)),
        as.data.frame(xcms::featureValues(
          x,
          value = "into",
          filled = TRUE
        ))
      )

      if ("has_monoisotopic_peak" %in% colnames(y) &&
          isTRUE(input$only_confirmed_monoisotopic)) {
        y <- y[y$has_monoisotopic_peak, , drop = FALSE]
      }

      y <- y[, !colnames(y) %in% c("peakidx", "monoisotopic_peak_ids"),
             drop = FALSE]

      saveRDS(
        x,
        file = file.path(
          input$resdir,
          paste0(
            input$jobname,
            "_annotation_",
            input$to_process_polarity,
            ".rds"
          )
        )
      )

      writexl::write_xlsx(
        y,
        path = file.path(
          input$resdir,
          paste0(
            input$jobname,
            "_features_",
            input$to_process_polarity,
            ".xlsx"
          )
        )
      )

      xReactive(x)
      resultsReactive(list(
        x = x,
        feature_table = y,
        ms1_candidates = S4Vectors::metadata(x)$LipidMS$ms1_candidates,
        ms2_evidence = S4Vectors::metadata(x)$LipidMS$ms2_annotation$evidence
      ))

      incProgress(1, detail = "Done")

      showTab("inTabset", "tab7")
      updateTabsetPanel(session, "inTabset", selected = "tab7")
    })
  })

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
      options = list(scrollX = TRUE, pageLength = 25)
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
      options = list(scrollX = TRUE, pageLength = 25)
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

# Run the application
shinyApp(ui = ui, server = server)


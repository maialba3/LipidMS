#============================================================================#
# Data
#============================================================================#
##############################################################################
# Reactive variables
msdataReactive <- reactiveVal(NULL) # msdata
fileClassesReactive <- reactiveVal(NULL) # Contains a file with lipid classes and adducts.

dbChoicesReactive <- reactiveVal(NULL) # The possible DB values for each lipid class for filtering and prediction.
polarityReactive <- reactiveVal(NULL) # The polarity of the analysis.
analysisReactive <- reactiveVal(NULL) # The type of data (msbatch, msobject, external).

parametersReactive <- reactiveVal(NULL) # Text of parameters displayed on screen in Run.
choicesNamesReactive <- reactiveVal(NULL) # Lipid classes selected for filtering and prediction


#============================================================================#
# User interface of the Shiny application
#============================================================================#
##############################################################################

ui <- fluidPage(
  tagList(
    tags$head(
      # Load Font Awesome for icons
      tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.5.0/css/all.min.css"),
      tags$style(HTML("

    body {
      padding-top: 10px !important;
    }

/* General style (in case no type is defined)  */
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

/* GREEN (success - message) */
.shiny-notification-message {
  background-color: #669973;
  border: 1px solid #4e7a5d;
}

/* Spinning circle animated loader */
.loader {
  display: inline-block;
  width: 16px;
  height: 16px;
  border: 3px solid white;
  border-top: 3px solid #e74818; /* top border color */
  border-radius: 50%;
  animation: spin 1s linear infinite;
  margin-left: 10px;
  vertical-align: middle;
}

@keyframes spin {
  0% { transform: rotate(0deg); }
  100% { transform: rotate(360deg); }
}


    #console_output {
      overflow-y: auto; /* Enable vertical scrolling if the content overflows. */
      white-space: pre-wrap; /* Adjust the text to preserve the format */
    }

    .dropstart .dropdown-menu {
      position: fixed !important;
      top: 50% !important;
      left: 50% !important;
      transform: translate(-50%, -50%) !important;
      width: 80vw !important; /* % of the total Shiny window it occupies when opened */
      max-width: 1100px;
      height: 95vh;
      overflow: auto;
      z-index: 1050;
      padding: 20px;
      background-color: white;
      box-shadow: 0px 4px 12px rgba(0, 0, 0, 0.15);
      border-radius: 8px;
      justify-content: center;
    }

        /* Checkbox */
        input[type='checkbox'] {
      accent-color: #669999;
        }
        /* Left margin for checkboxes */
        .checkbox {
         margin-left: 15px;
        }


        /* Radio buttons */
    .radio input[type='radio'] {
      accent-color: #669999;
    }
    /* Left margin for Radio buttons */
    .radio {
     margin-left: 15px;
    }


* Estilo visual del recuadro */
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

    /* Estilo cuando está seleccionado */
    .radio label.selected-radio {
      background-color: #f8f8f8;
      border: none;
      box-shadow: none;
      width: 95%;
      font-weight: bold;
    }



    .progress-bar {
      height: 80px;
      font-size: 10px;
      font-weight: bold;
      background-color: #669999;
      color: white;
    }

        /* Slider bar */
    .irs-bar,
    .irs-bar-edge {
      background: #669999 !important;
      border-color: #669999 !important;
    }


    /* Value label (if displayed in the middle of the slider) */
    .irs-single,
    .irs-from,
    .irs-to {
      background: #669999 !important;
      border-color: #669999 !important;
    }

    /* Disable TabPanel */
      .nav-tabs li.disabled a {
        color: grey !important;
        background-color: #f5f5f5 !important;
        cursor: not-allowed !important;
      }


  /* When the checkbox is disabled, the cursor also changes over the text */
  .disabled-checkbox label {
    cursor: not-allowed !important;
    opacity: 0.6;
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
      Shiny.addCustomMessageHandler('toggleTabs', function(message) {
        message.ids.forEach(function(id) {
          var tab = $('a[data-value=\"' + id + '\"]');
          if (message.disable) {
            tab.addClass('disabled-tab');
            tab.on('click.disabled', function(e) {
              e.preventDefault();
              e.stopImmediatePropagation();
              return false;
            });
            tab.parent().addClass('disabled');
          } else {
            tab.removeClass('disabled-tab');
            tab.off('click.disabled');
            tab.parent().removeClass('disabled');
          }
        });
      });


      // Aplica la clase al label del radio seleccionado
      $(document).on('shiny:inputchanged', function(event) {
        if (event.name === 'to_process_analysis') {
          // Solo afecta al grupo de to_process_analysis
          var $group = $('input[name=\"to_process_analysis\"]').closest('.form-group');
          $group.find('label').removeClass('selected-radio');
          $group.find('input[name=\"to_process_analysis\"]:checked').closest('label').addClass('selected-radio');
        }

        if (event.name === 'type_lipidclasses') {
          // Solo afecta al grupo de type_lipidclasses
          var $group = $('input[name=\"type_lipidclasses\"]').closest('.form-group');
          $group.find('label').removeClass('selected-radio');
          $group.find('input[name=\"type_lipidclasses\"]:checked').closest('label').addClass('selected-radio');
        }
      });


      // Aplica la clase al cargar la app
      $(document).ready(function() {
        $('input[name=\"to_process_analysis\"]:checked').closest('label').addClass('selected-radio');
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
      
      # tab1: Import Data
      tabPanel(
        title = "Data import", value = "tab1",
        
        div(
          style = "display: flex; align-items: center; margin: 0px 0px 20px 20px;",
          img(
            src = 'iconohorizontal.png',
            style = "margin-right: 10px; height: 50px;"
          ),
          tags$h5(
            HTML("<strong style='color:#2C3E50; font-size:18px;'> Lipid annotation for LC-MS/MS data</strong>"),
            style = "margin-right: 30px;"
          ),
          
          tags$div(
            style = "position: relative; display: inline-block;",
            
            # Circle with an icon acting as a button to view the tutorial image
            actionButton(
              inputId = "info_btn",
              label = NULL,
              icon = icon("info-circle", class = "fa-lg"),
              style = "
      background-color: #669999; 
      color: white; 
      border-radius: 50%; 
      border: 2px solid #669999;
      width: 35px; 
      height: 35px; 
      padding: 0;
      font-size: 18px;"
            ),
            
            tags$div(
              id = "dropdown_menu",
              style = "
      display: none; 
      position: fixed; 
      top: 50%; 
      left: 50%; 
      transform: translate(-50%, -50%);
      height: 90vh; 
      width: 80vw;
      overflow-y: auto; 
      background-color: rgba(255, 255, 255, 1); 
      padding: 20px; 
      font-size: 12px; 
      text-align: right; 
      border: none; 
      border-radius: 5px; 
      box-shadow: 0 4px 20px rgba(0,0,0,0.3);
      z-index: 1000;",
              
              # Zoom controls
              tags$div(
                style = "position: sticky; top: 0; background-color: rgba(255,255,255,0.95); z-index: 10; padding-bottom: 5px; text-align: right;",
                actionButton(
                  "zoom_out_dropdown_info", "-", 
                  style = "width: 25px; height: 25px; padding: 0; font-size: 16px; vertical-align: middle;"
                ),
                span(
                  id = "zoom_level_dropdown_info", "100%", 
                  style = "font-weight: bold; margin: 0 8px; vertical-align: middle;"
                ),
                actionButton(
                  "zoom_in_dropdown_info", "+", 
                  style = "width: 25px; height: 25px; padding: 0; font-size: 16px; vertical-align: middle;"
                )
              ),
              
              # Image
              tags$img(
                id = "info_dropdown",
                src = "info_dataImport.png",
                style = "width: 100%; height: auto; display: block; margin: 0 auto;"
              ),
              
              # Script JavaScript for zoom
              tags$script(HTML("
      var zoomDropdown_info = 100;
      var imgDropdown_info = document.getElementById('info_dropdown');
      var zoomLevelSpanDropdown_info = document.getElementById('zoom_level_dropdown_info');

      document.getElementById('zoom_in_dropdown_info').onclick = function() {
        if (zoomDropdown_info < 200) {
          zoomDropdown_info += 10;
          imgDropdown_info.style.width = zoomDropdown_info + '%';
          zoomLevelSpanDropdown_info.textContent = zoomDropdown_info + '%';
        }
      };
      document.getElementById('zoom_out_dropdown_info').onclick = function() {
        if (zoomDropdown_info > 30) {
          zoomDropdown_info -= 10;
          imgDropdown_info.style.width = zoomDropdown_info + '%';
          zoomLevelSpanDropdown_info.textContent = zoomDropdown_info + '%';
        }
      };
    "))
            ),
            
            # Script to show/hide the dropdown when clicking the circle
            tags$script(HTML("
    var infoBtn = document.getElementById('info_btn');
    var dropdownMenu = document.getElementById('dropdown_menu');

    infoBtn.onclick = function(event) {
      event.stopPropagation(); // avoid closing immediately
      if (dropdownMenu.style.display === 'none') {
        dropdownMenu.style.display = 'block';
      } else {
        dropdownMenu.style.display = 'none';
      }
    };

    // Cerrar al hacer click fuera
    document.addEventListener('click', function(event) {
      if (!dropdownMenu.contains(event.target) && !infoBtn.contains(event.target)) {
        dropdownMenu.style.display = 'none';
      }
    });

  "))
          )
          
        ),
        
        sidebarPanel(
          width = 3,
            fluidRow(
              # Select how to process the data
              radioButtons("to_process_analysis", "How do you want to process your data?",
                           choices = c("Batch Processing" = "msbatch",
                                       "Single Sample Processing" = "msobject"),
                           selected = "msbatch"),
              style = "margin: 20px 0px 0px 0px"
            ),
          fluidRow(
            column(12,
                   div(style = "border: 1px solid #ddd; border-radius: 6px; padding: 15px; margin: 0px 10px; background-color: #f9f9f9;",
            fluidRow(
              # Select analysis polarity
              radioButtons("to_process_polarity", "Polarity",
                           choices = c("Positive" = "positive",
                                       "Negative" = "negative"),
                           selected = NULL),
              style = "margin: 10px 0px 0px 0px"
            ),
            fluidRow(
              # Import mzXML files
              fileInput("file1", "Choose mzXML File/s",
                        multiple = TRUE,
                        accept = c("text/plain", ".mzXML")),
              style = "margin: 10px 0px 0px 0px"
            ),
            fluidRow(
              # Import metadata
              fileInput("metadata", "Metadata csv file",
                        multiple = TRUE,
                        accept = c("text/plain", ".csv")),
              style = "margin: 0px 0px 0px 0px"
            )
            ))),

          tags$hr(style = "border-color: white"),
          
          fluidRow(
            column(12,
                   tags$div(
                     tags$label("Choose lipid class source for analysis", style = "font-size: 16px; font-weight: bold;"),
                     radioButtons("type_lipidclasses", label = NULL,
                                  choices = list(
                                    "Default list" = "lipidclasses_default",
                                    "Custom list" = "lipidclasses_custom"),
                                  selected = "lipidclasses_default",
                                  width = "100%")
                   )
            )
          ),
          # If the default list is selected, nothing appears in this tab (the default list is displayed in the annotation and RT prediction tabs, the same in both, synchronized).
          
          conditionalPanel(
            condition = "input.type_lipidclasses == 'lipidclasses_custom'", # If the customized list is selected…
            
            fluidRow(
              column(12,
                     div(style = "border: 1px solid #ddd; border-radius: 6px; padding: 15px; margin: 10px 10px; background-color: #f9f9f9;",
                         # Separadores
                         fluidRow(
                           column(6,
                                  tags$label("Column delimiter", style = "margin-bottom: 0px; font-weight: bold;"),
                                  selectInput("classes_sep", NULL,
                                              choices = c("Semicolon (;)" = ";", "Comma (,)" = ",", "Pipe (|)" = "|"),
                                              selected = ";")
                           ),
                           column(6,
                                  tags$label("Adducts delimiter", style = "margin-bottom: 0px; font-weight: bold;"),
                                  selectInput("adducts_sep", NULL,
                                              choices = c("Semicolon (;)" = ";", "Comma (,)" = ",", "Pipe (|)" = "|"),
                                              selected = ",")
                           )
                         ),
                         
                         fluidRow(
                           column(12,
                                  tags$label("Choose a file to upload", style = "margin-top: 5px; margin-bottom: 0px; font-weight: bold;"),
                                  uiOutput("file_upload_ui")
                           )
                         )
                     )
              )
            )
            
          )
          ,
          
          # Previous and Next buttons
          tags$hr(style = "border-color: white"),
          fluidRow(
            actionButton("JumpTo2", "Next >", width = "100px"),
            style = "margin: 20px 0px 0px 0px"
          )
        ),
        
        mainPanel(
          conditionalPanel(
            condition = "output.import_summary != null"
          ),
          tableOutput("import_summary") # Here the table with the imported data information will be displayed
        )
      ),
      
      
      # Tab 2: Peak Selection
      tabPanel(title = "Peak-picking", value = "tab2",
               tags$head(
                   tags$style(HTML("input[type=\"number\"] {
                               height: 35px; font-size:11px
                               }"))),
                 mainPanel(
                   position="left",
                   fluidRow(column(6, h5(strong("dmzagglom (in ppm)")),
                                   h6(em("m/z tolerance used for partitioning and clustering. 5 by default."), style = "color:grey")),
                            column(3, numericInput("dmzagglom_ms1", "MS1", value = 15, min = 0, max = 100, step = 1)),
                            column(3, numericInput("dmzagglom_ms2", "MS2", value = 15, min = 0, max = 100, step = 1))),
                   fluidRow(column(6, h5(strong("drtagglom (in seconds)")),
                                   h6(em("rt window used for partitioning (in seconds). 25 by default."), style = "color:grey")),
                            column(3, numericInput("drtagglom_ms1", "MS1", value = 500, min = 0, max = 1000, step = 1)),
                            column(3, numericInput("drtagglom_ms2", "MS2", value = 500, min = 0, max = 1000, step = 1))),
                   fluidRow(column(6, h5(strong("drtclust (in seconds)")),
                                   h6(em("rt window used for clustering (in seconds). 25 by default."), style = "color:grey")),
                            column(3, numericInput("drtclust_ms1", "MS1", value = 100, min = 0, max = 1000, step = 1)),
                            column(3, numericInput("drtclust_ms2", "MS2", value = 200, min = 0, max = 1000, step = 1))),
                   fluidRow(column(6, h5(strong("minpeak")),
                                   h6(em("minimum number of measurements required for a peak. By default, 5 for MS1 and 4 for MS2."), style = "color:grey")),
                            column(3, numericInput("minpeak_ms1", "MS1", value = 5, min = 1, max = 25, step = 1)),
                            column(3, numericInput("minpeak_ms2", "MS2", value = 4, min = 1, max = 25, step = 1))),
                   fluidRow(column(6, h5(strong("minint")),
                                   h6(em("minimum intensity of a peak. By default, 1000 for MS1 and 100 for MS2."), style = "color:grey")),
                            column(3, numericInput("minint_ms1", "MS1", value = 1000, min = 0, max = 1e5, step = 1)),
                            column(3, numericInput("minint_ms2", "MS2", value = 100, min = 0, max = 1e5, step = 1))),
                   fluidRow(column(6, h5(strong("drtgap (in seconds)")),
                                   h6(em("maximum rt gap length to be filled. 5 by default."), style = "color:grey")),
                            column(3, numericInput("drtgap_ms1", "MS1", value = 5, min = 0, max = 30, step = 1)),
                            column(3, numericInput("drtgap_ms2", "MS2", value = 5, min = 0, max = 30, step = 1))),
                   fluidRow(column(6, h5(strong("drtminpeak (in seconds)")),
                                   h6(em("minimum rt width of a peak. 15 by default. At least minpeak within the drtminpeak window are required to define a peak."), style = "color:grey")),
                            column(3, numericInput("drtminpeak_ms1", "MS1", value = 15, min = 0, max = 500, step = 1)),
                            column(3, numericInput("drtminpeak_ms2", "MS2", value = 15, min = 0, max = 500, step = 1))),
                   fluidRow(column(6, h5(strong("drtmaxpeak (in seconds)")),
                                   h6(em("maximum rt width of a single peak. 100 by default."), style = "color:grey")),
                            column(3, numericInput("drtmaxpeak_ms1", "MS1", value = 100, min = 0, max = 500, step = 1)),
                            column(3, numericInput("drtmaxpeak_ms2", "MS2", value = 200, min = 0, max = 500, step = 1))),
                   fluidRow(column(6, h5(strong("maxeicpeaks")),
                                   h6(em("maximum number of peaks within one EIC. By default, 5 for 10 and 3 for MS2."), style = "color:grey")),
                            column(3, numericInput("recurs_ms1", "MS1", value = 5, min = 1, max = 50, step = 1)),
                            column(3, numericInput("recurs_ms2", "MS2", value = 10, min = 1, max = 50, step = 1))),
                   fluidRow(column(6, h5(strong("weight")),
                                   h6(em("weight for assigning measurements to a peak. By default, 2 for MS1 and 3 for MS2."), style = "color:grey")),
                            column(3, numericInput("weight_ms1", "MS1", value = 2, min = 1, max = 10, step = 1)),
                            column(3, numericInput("weight_ms2", "MS2", value = 3, min = 1, max = 10, step = 1))),
                   fluidRow(column(6, h5(strong("SN")),
                                   h6(em("signal-to-noise ratio. By default, 3 for MS1 and 2 for MS2."), style = "color:grey")),
                            column(3, numericInput("sn_ms1", "MS1", value = 3, min = 1, max = 10, step = 1)),
                            column(3, numericInput("sn_ms2", "MS2", value = 2, min = 1, max = 10, step = 1))),
                   fluidRow(column(6, h5(strong("SB")),
                                   h6(em("signal-to-base ratio. By default, 3 for MS1 and 2 for MS2."), style = "color:grey")),
                            column(3, numericInput("sb_ms1", "MS1", value = 3, min = 1, max = 10, step = 1)),
                            column(3, numericInput("sb_ms2", "MS2", value = 2, min = 1, max = 10, step = 1))),
                   fluidRow(column(6, h5(strong("dmzIso (in ppm)")),
                                   h6(em("mass tolerance for isotope matching. 5 by default."), style = "color:grey")),
                            column(3, numericInput("dmzIso_ms1", "MS1", value = 5, min = 0, max = 100, step = 1)),
                            column(3, numericInput("dmzIso_ms2", "MS2", value = 5, min = 0, max = 10, step = 1))),
                   fluidRow(column(6, h5(strong("drtIso")),
                                   h6(em("rt window for isotope matching. 5 by default."), style = "color:grey")),
                            column(3, numericInput("drtIso_ms1", "MS1", value = 5, min = 0, max = 100, step = 1)),
                            column(3, numericInput("drtIso_ms2", "MS2", value = 5, min = 0, max = 100, step = 1))),
                   tags$hr(),
                   
                   # Previous and Next buttons
                   fluidRow(actionButton("GoBackTo1", "< Previous", width = "100px", 
                                         style="margin:40px 0px"), 
                            actionButton("JumpTo3", "Next >", width = "100px", 
                                         style="margin:40px 0px"))
                 )
               
      ),
      
      # tab3: Process the data
      tabPanel(title = "Batch processing", value = "tab3",
               tags$head(
                   tags$style(HTML("input[type=\"number\"] { height: 35px; font-size:11px }"))
                 ),
                 mainPanel(
                   fluidRow(column(6, h5(strong("dmzalign")),
                                   h6(em("mass tolerance between peak groups for alignment (in ppm). 5 by default."), style = "color:grey")),
                            column(3, numericInput("dmzalign", "", value = 5, min = 0, max = 100, step = 1))),
                   fluidRow(column(6, h5(strong("drtalign")),
                                   h6(em("maximum rt distance between peaks for alignment (in seconds). 30 by default."), style = "color:grey")),
                            column(3, numericInput("drtalign", "", value = 30, min = 0, max = 200, step = 1))),
                   fluidRow(column(6, h5(strong("span")),
                                   h6(em("span parameter for loess rt smoothing. 0.4 by default."), style = "color:grey")),
                            column(3, numericInput("span", "", value = 0.4, min = 0, max = 1, step = 0.05))),
                   fluidRow(column(6, h5(strong("minsamplesfracalign")),
                                   h6(em("minimum samples fraction represented in each cluster used for alignment. 0.75 by default."), style = "color:grey")),
                            column(3, numericInput("minsamplesfracalign", "", value = 0.75, min = 0, max = 1, step = 0.05))),
                   fluidRow(column(6, h5(strong("dmzgroup")),
                                   h6(em("mass tolerance between peak groups for grouping (in ppm). 5 by default."), style = "color:grey")),
                            column(3, numericInput("dmzgroup", "", value = 5, min = 0, max = 100, step = 1))),
                   fluidRow(column(6, h5(strong("drtagglomgroup")),
                                   h6(em("maximum rt distance in mz partitions for grouping (in seconds). 30 by default. It shouldn't be smaller than drtgroup."), style = "color:grey")),
                            column(3, numericInput("drtagglomgroup", "", value = 30, min = 0, max = 200, step = 1))),
                   fluidRow(column(6, h5(strong("drtgroup")),
                                   h6(em("maximum rt distance between peaks for grouping (in seconds). 30 by default."), style = "color:grey")),
                            column(3, numericInput("drtgroup", "", value = 15, min = 0, max = 200, step = 1))),
                   fluidRow(column(6, h5(strong("minsamplesfracgroup")),
                                   h6(em("minimum samples fraction represented in each cluster used for grouping. 0.25 by default."), style = "color:grey")),
                            column(3, numericInput("minsamplesfracgroup", "", value = 0.25, min = 0, max = 1, step = 0.05))),
                   fluidRow(column(6, h5(strong("parallel")),
                                   h6(em("parallelize processing. FALSE by default."), style = "color:grey")),
                            column(3, selectInput("parallel", "", choices = list("TRUE" = TRUE, "FALSE" = FALSE), selected = TRUE))),
                   fluidRow(column(6, h5(strong("ncores")),
                                   h6(em("number of cores to parallelize"), style = "color:grey")),
                            column(3, numericInput("ncores", "", value = 4, min = 1, max = 16, step = 1))),
                   fluidRow(column(6, h5(strong("global_gb")),
                                   h6(em("memory available for each core when parallelizing"), style = "color:grey")),
                            column(3, numericInput("global_gb", "", value = 1024, min = 8, max = Inf, step = 1))),
                   tags$hr(),
                   # Previous and next buttons
                   fluidRow(actionButton("GoBackTo2", "< Previous", width = "100px", style="margin:40px 0px"),
                            actionButton("JumpTo4", "Next >", width = "100px", style="margin:40px 0px"))
                 )

      ),
      
      # tab4: Annotate
      tabPanel(title = "Annotation", value = "tab4",
               mainPanel(fluidRow(column(6, h5(strong("dmzprecursor")),
                                           h6(em("Mass tolerance for precursor ions. 5 by default."), style = "color:grey")),
                                    column(3, numericInput("dmzprecursor", "", value = 5, min = 0, max = 100, step = 1))),
                           fluidRow(column(6, h5(strong("dmzproducts")),
                                           h6(em("mass tolerance for product ions. 10 by default."), style = "color:grey")),
                                    column(3, numericInput("dmzproducts", "", value = 10, min = 0, max = 100, step = 1))),
                           fluidRow(column(6, h5(strong("rttol")),
                                           h6(em("total rt window for coelution between precursor and product ions. 5 by default."), style = "color:grey")),
                                    column(3, numericInput("rttol", "", value = 5, min = 0, max = 50, step = 1))),
                           fluidRow(column(6, h5(strong("coelcutoff")),
                                           h6(em("coelution score threshold between parent and fragment ions. Only applied if rawData info is supplied. 0.7 by default."), style = "color:grey")),
                                    column(3, numericInput("coelcutoff", "", value = 0.7, min = 0, max = 1, step = 0.05))),
                           tags$hr(),
                           
                           # The list of lipid classes and adducts (default or customized) for annotation appears, where you can select which specific classes you want or do not want (synchronized with the RT prediction list—any change in one updates the other).
                           fluidRow(column(12, checkboxGroupInput("lipidclasses_annotate", "Lipid classes to annotate:",
                                                                  choices = NULL, 
                                                                  selected = NULL,
                                                                  inline = TRUE))),
                           tags$hr(),
                           
                           # Previous and next buttons
                           fluidRow(actionButton("GoBackTo3", "< Previous", width = "100px",
                                                 style="margin:40px 0px"),
                                    actionButton("JumpTo5", "Next >", width = "100px",
                                                 style="margin:40px 0px"))
                 )
      ),
      
      
      # tab7: Function Execution 
      tabPanel(title = "Run", value = "tab5",
               
               # LipidMS text and logo at the top
               div(
                 style = "display: flex; align-items: center; margin: 0px 0px 20px 20px;",
                 img(src = 'iconohorizontal.png',
                     style = "margin-right: 10px; height: 50px;"),
                 tags$h5(
                   HTML("<strong style='color:#2C3E50; font-size:18px;'> Lipid annotation for LC-MS/MS data</strong>"),
                   style = "margin-right: 30px;"
                 )),
               
               # Panel with options to run the analysis
               sidebarPanel(width = 3,
                            fluidRow(
                              # Add the name that will later appear when downloading the files
                              textInput("jobname", "Job Name", 
                                        value = paste("Job", as.character(Sys.Date()), sep = "_")),
                              style = "margin: 20px 0px 0px 0px"
                            ),
                            tags$hr(style = "border-color: white"),

                            
                            # Button to perform the analysis, executes all functions
                            fluidRow(actionButton("Run_analysis", "Run analysis",
                                                  style="color: #fff; background-color: #669999; border-color: #000; width: 100%;"),
                                     style = "margin: 0px 0px 0px 0px"),
                            tags$hr(style = "border-color: white"),
                            
                            # Previous and Next buttons
                            uiOutput("navButtonsRun")
                          ),
               
               column(9, 
                      tags$div(
                        style = "padding: 0px; height: 75vh; overflow-y: auto; position: relative;",
                        
                        # Button to download an .html of the parameters selected so far for the analysis
                        downloadButton(
                          "download_parameters_run", 
                          "Selected Parameters", 
                          style = "position: sticky; top: 10px; right: 10px; z-index: 1000; float: right;"),
                        
                        # Informative text showing the parameters selected so far for the analysis, updating as changes are made
                        uiOutput("parameters_output")
                      )
               )
               
      ),
      
      # tab8: Results
      tabPanel(title = "Results", value = "tab6", # If there are no results, this tab is not displayed
               tags$head(
                 tags$style(HTML("
    .nav-tabs li a {
      display: block;
      width: 100%;
      text-align: center;
      color: #577b9e;
      font-weight: bold;
      border: 1px solid #f9f9f9 !important;
      margin-bottom: 3px;
    }

        .nav-tabs > li.active > a,
    .nav-tabs > li.active > a:focus,
    .nav-tabs > li.active > a:hover {
      color: #2c3e50 !important;
      background-color: #f9f9f9;
      border: 2px solid #435f7a !important;
      margin-bottom: 3px;
    }

    .nav-tabs > li > a:hover {
      color: #2c3e50;
    }

    .nav-tabs li {
      width: 100%;
    }

    .no-width-change-tabs .nav-tabs li a {
      display: block;
      width: 50%;
      text-align: center;
    }
    .no-width-change-tabs .nav-tabs li {
      width: 50%;
    }
  ")),
                 
               ),
               
               # Sidebar panel with the results tabs
               fluidRow(
                 column(width = 2,
                        div(
                          style = "display: flex; align-items: center; margin: 0px 0px 20px 20px;",
                          img(src = 'iconohorizontal.png',
                              style = "margin-right: 10px; height: 50px;")),
                        # Only tabs with data to display will be shown; otherwise, they will not appear
                        sidebarPanel(width = 12,
                                     
                                     # SECTION 1: Interactive result exploration
                                     tags$h5(HTML("<strong style='color:#2C3E50; font-size:15px;'> Explore Results</strong>")),
                                     tags$p("Use the tabs above to view summary tables, annotated peaklists and batch results.",
                                            style = "color:#7F8C8D; font-size:13px; margin-bottom:10px;"),
                                     # Tabs
                                     tabsetPanel(id = "tabs", type = "tabs",
                                                 tabPanel("Summary Tables Annotation", value = "summary"),
                                                 tabPanel("Annotated Peaklists", value = "peaklists"),
                                                 tabPanel("Batch Results", value = "batch")
                                     ),
                                     
                                     # SECTION 2: Download individual result files
                                     tags$hr(style = "border-color: white"),
                                     tags$h5(HTML("<strong style='color:#2C3E50; font-size:15px;'> Download individual files</strong>")),
                                     fluidRow(
                                       tags$div(uiOutput("download_summaryTable_ui"), style = "margin-right: 10px; margin-bottom: 3px;"),
                                       tags$div(uiOutput("download_annotatedPeaklist_ui"), style = "margin-right: 10px; margin-bottom: 3px;"),
                                       tags$div(uiOutput("download_batchResults_ui"), style = "margin-right: 10px;"),
                                       style = "margin: 5px 0px 0px 0px"
                                     ),
                                     
                                     # SECTION 3: Download full rtdata object
                                     tags$hr(style = "border-color: white"),
                                     tags$h5(HTML("<strong style='color:#2C3E50; font-size:15px;'> Download full rtdata object</strong>")),
                                     tags$p("Export the entire analysis session, including all results.",
                                            style = "color:#7F8C8D; font-size:13px; margin-bottom:10px;"),
                                     fluidRow(
                                       tags$div(uiOutput("download_msdata_ui"), style = "margin-bottom: 5px;"),
                                       tags$div(uiOutput("export_msdata_ui")),
                                       style = "margin: 5px 0px 0px 0px"
                                     ),
                                     
                                     tags$hr(style = "border-color: white"),
                                     
                                     # Previous and Next buttons
                                     fluidRow(actionButton("GoBackTo5", "< Previous",
                                                           width = "100px", style = "margin:40px 0px"),
                                              style = "margin: 20px 0px 0px 0px")
                        )
                        
                 ),
                 
                 # Remaining space for dynamic content that displays the results of the selected tab
                 column(width = 10, 
                        conditionalPanel(
                          condition = "input.tabs == 'summary'", # If 'Summary Tables Annotation' is selected...
                          uiOutput("summaryTable_ui") # displays the summary results
                        ),
                        
                        conditionalPanel(
                          condition = "input.tabs == 'peaklists'", # If 'Annotated Peaklists' is selected...
                          uiOutput("annotatedPeaklist_ui") # displays the peaklist results
                        ),
                        
                        conditionalPanel(
                          condition = "input.tabs == 'batch'", # If 'Batch Results' is selected...
                          uiOutput("features_ui") # displays the feature results
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
  
  #============================================================================#
  # Button functionality: Next, Previous
  #============================================================================#
  observeEvent(input$JumpTo2, {updateTabsetPanel(session, "inTabset", selected = "tab2")}) # Data Import -> Peak-picking
  observeEvent(input$GoBackTo1, {updateTabsetPanel(session, "inTabset", selected = "tab1")}) # Peak-picking -> Data Import
  observeEvent(input$JumpTo3, {updateTabsetPanel(session, "inTabset", selected = "tab3")}) # Peak-picking -> Batch processing
  observeEvent(input$GoBackTo2, {updateTabsetPanel(session, "inTabset", selected = "tab2")}) # Batch processing -> Peak-picking
  observeEvent(input$JumpTo4, {updateTabsetPanel(session, "inTabset", selected = "tab4")}) # Batch processing -> Annotation
  observeEvent(input$GoBackTo3, {updateTabsetPanel(session, "inTabset", selected = "tab3")}) # Annotation -> Batch processing
  observeEvent(input$JumpTo5, {updateTabsetPanel(session, "inTabset", selected = "tab5")}) # Annotation ->  Run
  observeEvent(input$GoBackTo4, {updateTabsetPanel(session, "inTabset", selected = "tab4")}) # Run -> Annotation
  observeEvent(input$JumpTo6, {updateTabsetPanel(session, "inTabset", selected = "tab6")}) # Run -> Results
  observeEvent(input$GoBackTo5, {updateTabsetPanel(session, "inTabset", selected = "tab5")}) # Results <- Run
  
  observe({
    msdata <- msdataReactive()
    analysis_type <- analysisReactive()
    
    if (is.null(msdata)) {
      hideTab("inTabset", "tab6")
      output$navButtonsRun <- renderUI({
        fluidRow(
          actionButton("GoBackTo4", "< Previous", style = "margin:10px 10px"),
          style = "margin: 0px 0px 0px 0px"
        )
      })
    } else {
      showTab("inTabset", "tab6")
      if (analysis_type == "msobject") {
        hideTab("tabs", "batch")
      }
      
      output$navButtonsRun <- renderUI({
        fluidRow(
          actionButton("GoBackTo4", "< Previous", style = "margin:10px 10px"),
          actionButton("JumpTo6", "Next >", style = "margin:10px 10px"),
          style = "margin: 0px 0px 0px 0px"
        )
      })
    }
  })
  
  output$file_upload_ui <- renderUI({
    fileInput("classes_file", NULL,
              multiple = FALSE,
              accept = c("text/plain", ".csv"),
              width = "100%",
              placeholder = "No file selected")
  })

  #============================================================================#
  # Dynamic text: Selected parameters from Run
  #============================================================================#
  output$parameters_output <- renderUI({

    # tittle: Selected parameters
    text_output <- paste0(
      "<h4 style='color:#669999; font-weight:bold; margin-bottom:15px'>",
      "<i class='fa-solid fa-circle-info' style='color:#669999; margin-right:8px;'></i>",
      "Selected Parameters</h4>",
      
      "<p style='color:#666666; font-style:italic; margin-top:25px; margin-bottom:25px;'>",
      "If you want to save the selected parameters as an <b>.html</b> file, make sure to do so <b>before clicking the 'Run'</b> button.<br>",
      "Otherwise, the displayed values may change dynamically during the analysis.",
      "</p>"
    )
    
    # Choices: How do you want to process your data
    choices_process <- c(
      "Batch Processing" = "msbatch",
      "Single Sample Processing" = "msobject"
    )
    # Choices: Polarity
    choices_polarity <- c(
      "Positive" = "positive",
      "Negative" = "negative"
    )
    
    # Text displayed on screen
    # Job name
    text_output <- paste0(
      text_output,
      "<table style='border-collapse: collapse; margin-left: 10px; margin-bottom: 10px;'>",
      
      "<tr>",
      "<td style='font-weight:bold; color:#2c3e50; padding-right:30px;'>●&nbsp;&nbsp; Job name</td>",
      "<td style='color:#2c3e50;'>",
      input$jobname,
      "</td>",
      "</tr>",

      
      # Only for raw data: How do you want to process your data
        paste0(
          "<tr>",
          "<td style='font-weight:bold; color:#2c3e50; padding-right:30px;'>●&nbsp;&nbsp; How do you want to process your data?</td>",
          "<td style='color:#2c3e50;'>",
          names(choices_process)[choices_process == input$to_process_analysis],
          "</td>",
          "</tr>"
        )
      ,
      
      # Polarity
      "<tr>",
      "<td style='font-weight:bold; color:#2c3e50; padding-right:30px;'>●&nbsp;&nbsp; Polarity</td>",
      "<td style='color:#2c3e50;'>",
      names(choices_polarity)[choices_polarity == polarityReactive()],
      "</td>",
      "</tr>",
      
      "</table>",
      "<hr style='border:none; height:1px; background-color:#ffffff;'>"
    )
    
    # Choices: type list lipid classes (default or custom)
    if (input$type_lipidclasses == "lipidclasses_default") {
      text_type_lipidclasses <- "&nbsp;&nbsp; (default list)"
    } else if (input$type_lipidclasses == "lipidclasses_custom") {
      text_type_lipidclasses <- "&nbsp;&nbsp; (custom list)"
    }
    
    
    # --- PEAK-PICKING: Annotation ---
      text_output <- paste0(text_output,
                            "<p><i class='fa-solid fa-square-check' style='color:#95a5a6;'></i> <b> Peak-picking</b></p>",
                            
                            # Table with header
                            "<table style='border-collapse: collapse; margin-left: 50px;'>",
                            
                            # Header row
                            "<tr>",
                            "<th style='padding-right: 60px;'></th>",
                            "<th style='color:#2c3e50; text-align:center;padding-right:50px;'>MS1</th>",
                            "<th style='color:#2c3e50; text-align:center;padding-right:30px;'>MS2</th>",
                            "</tr>",
                            
                            # Row 1
                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right:50px;'>●&nbsp;&nbsp; dmzagglom</td>",
                            "<td>", input$dmzagglom_ms1, " ppm</td>",
                            "<td>", input$dmzagglom_ms2, " ppm</td>",
                            "</tr>",
                            
                            # Row 2
                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50;'>●&nbsp;&nbsp; drtagglom</td>",
                            "<td>", input$drtagglom_ms1, " s</td>",
                            "<td>", input$drtagglom_ms2, " s</td>",
                            "</tr>",
                            
                            # Row 3 ...
                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50;'>●&nbsp;&nbsp; drtclust</td>",
                            "<td>", input$drtclust_ms1, " s</td>",
                            "<td>", input$drtclust_ms2, " s</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50;'>●&nbsp;&nbsp; minpeak</td>",
                            "<td>", input$minpeak_ms1, "</td>",
                            "<td>", input$minpeak_ms2, "</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50;'>●&nbsp;&nbsp; minint</td>",
                            "<td>", input$minint_ms1, "</td>",
                            "<td>", input$minint_ms2, "</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50;'>●&nbsp;&nbsp; drtgap</td>",
                            "<td>", input$drtgap_ms1, " s</td>",
                            "<td>", input$drtgap_ms2, " s</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50;'>●&nbsp;&nbsp; drtminpeak</td>",
                            "<td>", input$drtminpeak_ms1, " s</td>",
                            "<td>", input$drtminpeak_ms2, " s</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50;'>●&nbsp;&nbsp; drtmaxpeak</td>",
                            "<td>", input$drtmaxpeak_ms1, " s</td>",
                            "<td>", input$drtmaxpeak_ms2, " s</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50;'>●&nbsp;&nbsp; maxeicpeaks</td>",
                            "<td>", input$recurs_ms1, "</td>",
                            "<td>", input$recurs_ms2, "</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50;'>●&nbsp;&nbsp; weight</td>",
                            "<td>", input$weight_ms1, "</td>",
                            "<td>", input$weight_ms2, "</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50;'>●&nbsp;&nbsp; SN</td>",
                            "<td>", input$sn_ms1, "</td>",
                            "<td>", input$sn_ms2, "</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50;'>●&nbsp;&nbsp; SB</td>",
                            "<td>", input$sb_ms1, "</td>",
                            "<td>", input$sb_ms2, "</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50;'>●&nbsp;&nbsp; dmzIso</td>",
                            "<td>", input$dmzIso_ms1, " ppm</td>",
                            "<td>", input$dmzIso_ms2, " ppm</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50;'>●&nbsp;&nbsp; drtIso</td>",
                            "<td>", input$drtIso_ms1, " s</td>",
                            "<td>", input$drtIso_ms2, " s</td>",
                            "</tr>",
                            
                            "</table>"
      )
      
      
      # BATCH PROCESSING - ANNOTATION (TRUE) 
      text_output <- paste0(text_output,
                            "<hr>",
                            "<p><i class='fa-solid fa-square-check' style='color:#95a5a6;'></i> <b> Batch processing</b></p>",
                            
                            # Table with header
                            "<table style='border-collapse: collapse; margin-left: 50px;'>",
                            
                            # Header row
                            "<tr>",
                            "<th style='padding-right: 50px;'></th>", 
                            "<th'></th>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; dmzalign</td>",
                            "<td>", input$dmzalign, " ppm</td>",
                            "</tr>",   

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; drtalign</td>",
                            "<td>", input$drtalign, " s</td>",
                            "</tr>",  

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; span</td>",
                            "<td>", input$span, " </td>",
                            "</tr>",  

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; minsamplesfracalign</td>",
                            "<td>", input$minsamplesfracalign, "</td>",
                            "</tr>", 

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; dmzgroup</td>",
                            "<td>", input$dmzgroup, " ppm</td>",
                            "</tr>", 

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; drtagglomgroup</td>",
                            "<td>", input$drtagglomgroup, " s</td>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; drtgroup</td>",
                            "<td>", input$drtgroup, " s</td>",
                            "</tr>", 

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; minsamplesfracgroup</td>",
                            "<td>", input$minsamplesfracgroup, "</td>",
                            "</tr>", 

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; parallel</td>",
                            "<td>", input$parallel, "</td>",
                            "</tr>", 

                            
                            "</table>"
                            
      )
      
      # ANNOTATION - ANNOTATION (TRUE) 
      text_output <- paste0(text_output,
                            "<hr>",
                            "<p><i class='fa-solid fa-square-check' style='color:#95a5a6;'></i> <b> Annotation</b></p>",

                            "<table style='border-collapse: collapse; margin-left: 50px;'>",

                            "<tr>",
                            "<th style='padding-right: 50px;'></th>",
                            "<th'></th>",
                            "</tr>",

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; dmzprecursor</td>",
                            "<td>", input$dmzprecursor, "</td>",
                            "</tr>",   

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; dmzproducts</td>",
                            "<td>", input$dmzproducts, "</td>",
                            "</tr>",  

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; rttol</td>",
                            "<td>", input$rttol, " s</td>",
                            "</tr>",  

                            "<tr>",
                            "<td style='font-weight:bold; color:#2c3e50; padding-right: 30px;'>●&nbsp;&nbsp; coelcutoff</td>",
                            "<td>", input$coelcutoff, "</td>",
                            "</tr>", 
                            
                            "</table>"
                            
      )
      
      # dropdown row with the lipid classes to annotate selected by the user
      if (!is.null(choicesNamesReactive())) {
        lipid_adducts <- names(choicesNamesReactive())
        # selected lipid classes by the user
        selected_lipids <- lipid_adducts[lipid_adducts %in% input$lipidclasses_annotate]
        
        unique_id <- paste0("lipid_details_", as.integer(runif(1, 1000, 9999)))  # unique ID
        
        if (length(selected_lipids) == 0) {
          lipid_text <- "<div style='margin-left: 30px; color:gray; font-style: italic;'>No lipid classes selected to annotate.</div>"
        } else {
          # Table without header
          lipid_text <- "<table style='border-collapse: collapse; margin-left: 30px; margin-top:10px;'>"
          
          for (lipid in selected_lipids) {
            html_entry <- choicesNamesReactive()[[lipid]]
            
            # Extract lipid class
            lipid_class <- sub(".*<b[^>]*>([^<]+)</b>.*", "\\1", html_entry)
            lipid_class <- gsub("&nbsp;", "", lipid_class)
            
            # Extract adducts and clean
            adducts_raw <- sub(".*<small>\\(([^<]+)\\)</small>.*", "\\1", html_entry)
            adducts_clean <- gsub("&nbsp;", "", adducts_raw)
            adduct_list <- unlist(strsplit(adducts_clean, ";"))
            adducts_text <- paste(trimws(adduct_list), collapse = paste0(" ", input$adducts_sep, " "))
            
            # Compact row with lipid in bold
            lipid_text <- paste0(
              lipid_text,
              "<tr>",
              "<td style='padding: 1px 15px 1px 0; color:#2c3e50; font-size: 85%;'><b>", lipid_class, "</b></td>",
              "<td style='padding: 1px 0; color:#444; font-size: 85%;'>", adducts_text, "</td>",
              "</tr>"
            )
          }
          
          lipid_text <- paste0(lipid_text, "</table>")
        }
        
        
        
        # Expandable
        expandable_row <- paste0(
          "<div style='margin-left: 50px; margin-top: 20px;'>",
          "<div style='color:#2c3e50; font-weight:bold; cursor:pointer;' onclick='toggleDetails(\"", unique_id, "\")'>",
          "<span id='", unique_id, "_toggle'>△</span> Show lipid classes to annotate",
          "<span style='font-weight:normal;'>", text_type_lipidclasses, "</span></div>",
          "<div id='", unique_id, "' style='display:block; margin-top: 10px; color:#444;'>",
          lipid_text,
          "</div>",
          "</div>"
        )
        
        # Output final lipid classes list
        text_output <- paste0(
          text_output,
          expandable_row,
          "<script>
      function toggleDetails(id) {
        var content = document.getElementById(id);
        var toggleIcon = document.getElementById(id + '_toggle');
        if (content.style.display === 'none') {
          content.style.display = 'block';
          toggleIcon.innerHTML = '△';
        } else {
          content.style.display = 'none';
          toggleIcon.innerHTML = '▽';
        }
      }
    </script>"
        )
      }
      
    
    parametersReactive(text_output)
    HTML(text_output)
  })
  

  #============================================================================#
  # Download .html with the selected parameters
  #============================================================================#
  output$download_parameters_run <- downloadHandler(
    filename = function(){paste(input$jobname, "parameters.html", sep="_")},
    content = function(file) {
      raw_html <- parametersReactive()
      
      full_html <- paste0(
        "<!DOCTYPE html>\n<html>\n<head>\n",
        "<meta charset='UTF-8'>\n",
        "<link rel='stylesheet' href='https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css'>\n",
        "<link rel='stylesheet' href='https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.5.0/css/all.min.css'>\n",
        "<style>\n",
        "body { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; font-size: 14px; color: #2c3e50; padding: 30px; margin-left: 50px; }\n",
        "table { width: auto; margin-bottom: 20px; }\n",
        "th, td { padding: 6px 10px; border: none; }\n",
        "</style>\n",
        "</head>\n<body>\n",
        raw_html,
        "\n</body>\n</html>"
      )
      
      writeLines(full_html, file)
    }
  )
  
  
  
  #============================================================================#
  # Import Data text
  #============================================================================#
  output$import_summary <- renderUI({
    
    # --- Main title ---
    text_output <- paste0(
      "<h4 style='color:#669999; font-weight:bold; margin-bottom:15px'>",
      "<i class='fa-solid fa-circle-info' style='color:#669999; margin-right:8px;'></i>",
      "Data Import</h4>"
    )
    
      
      # --- FILES ---
      # If .mzXML files are uploaded
      if (!is.null(input$file1) && length(input$file1$name) > 0) {
        files_checkbox <- "<p><input type='checkbox' checked disabled><b> .mzXML files</b></p>"
        
        unique_id_files <- paste0("files_details_", as.integer(runif(1, 1000, 9999)))  
        files_text <- "<ul style='margin-left: 10px;'>"
        for (f in input$file1$name) {
          files_text <- paste0(files_text, "<li>", f, "</li>")
        }
        files_text <- paste0(files_text, "</ul>")
        
        files_expandable <- paste0(
          "<hr>",
          files_checkbox,
          "<div style='margin-left:50px; margin-right:80px; border:1px solid #ddd; border-radius:6px; padding:5px; background-color:#f9f9f9; font-size:14px; color:#2C3E50;'>",
          "<div style='color:#2c3e50; font-weight:bold; cursor:pointer;' onclick='toggleDetails(\"", unique_id_files, "\")'>",
          "<span id='", unique_id_files, "_toggle'>△</span> Uploaded files</div>",
          "<div id='", unique_id_files, "' style='display:block; margin-top: 10px; color:#444;'>",
          files_text,
          "</div></div>"
        )
        
        text_output <- paste0(text_output, files_expandable)
        
      } else {
        # No files uploaded yet
        text_output <- paste0(
          text_output,
          "<hr>",
          "<p><input type='checkbox' disabled><b> .mzXML files</b></p>",
          "<div style='margin-left:50px; margin-right:80px; border:1px solid #ddd; border-radius:6px; padding:5px; background-color:#f9f9f9; font-size:14px; color:#2C3E50;'>
                  <p style='margin-left:20px; margin-top:10px; color:gray; font-style:italic; font-size:14px;'>
          No data files uploaded yet. Please upload a valid <b><b>.mzXML</b></b> files.</p></div>"
        )
      }
      
      # --- METADATA ---
      # If metadata CSV is uploaded
      if (!is.null(input$metadata) && file.exists(input$metadata$datapath)) {
        metadata_checkbox <- "<p><input type='checkbox' checked disabled><b> Metadata file</b></p>"
        
        # Read CSV
        metadata2 <- utils::read.csv(input$metadata$datapath, sep = ",")
        if(!all(grepl(".mzXML", metadata2$sample))) {
          metadata2$sample <- paste0(metadata2$sample, ".mzXML")
        }
        
        # HTML table with first 3 columns
        metadata_table <- "<table style='border-collapse: collapse; margin-left: 30px;'>"
        metadata_table <- paste0(
          metadata_table,
          "<tr>",
          paste0("<th style='padding: 5px 10px;'>", colnames(metadata2)[1:3], "</th>", collapse = ""),
          "</tr>"
        )
        for(i in 1:nrow(metadata2)) {
          metadata_table <- paste0(
            metadata_table,
            "<tr>",
            "<td style='padding: 2px 30px 2px 10px;'>", metadata2[i,1], "</td>",
            "<td style='padding: 2px 10px;'>", metadata2[i,2], "</td>",
            "<td style='padding: 2px 10px;'>", metadata2[i,3], "</td>",
            "</tr>"
          )
        }
        metadata_table <- paste0(metadata_table, "</table>")
        
        unique_id_metadata <- paste0("metadata_details_", as.integer(runif(1,1000,9999)))
        
        metadata_expandable <- paste0(
          "<hr>",
          metadata_checkbox,
          "<div style='margin-left:50px; margin-right:80px; border:1px solid #ddd; border-radius:6px; padding:5px; background-color:#f9f9f9; font-size:14px; color:#2C3E50;'>",
          "<div style='color:#2c3e50; font-weight:bold; cursor:pointer;' onclick='toggleDetails(\"", unique_id_metadata, "\")'>",
          "<span id='", unique_id_metadata, "_toggle'>△</span> Metadata details</div>",
          "<div id='", unique_id_metadata, "' style='display:block; margin-top: 10px; color:#444;'>",
          metadata_table,
          "</div></div>"
        )
        
        text_output <- paste0(text_output, metadata_expandable)
        
      } else {
        # No metadata uploaded yet
        
        text_expandible <- "<p style='margin-left:20px; margin-top:10px; color:gray; font-style:italic; font-size:14px;'>
           No metadata file uploaded yet. Please upload a valid <b>.csv</b> file.</p>"
        
        text_output <- paste0(
          text_output,
          "<hr>",
          "<p><input type='checkbox' disabled><b> Metadata file</b></p>",
          
          "<div style='margin-left:50px; margin-right:80px; border:1px solid #ddd; border-radius:6px; padding:5px; background-color:#f9f9f9; font-size:14px; color:#2C3E50;'>",
          text_expandible,
          
          "<p style='margin-left:20px; color:gray; font-style:italic; font-size:14px;'>
                          The metada file must contain the following columns, separated by <b>','</b> :</p>

                          <table style='margin-left:50px; margin-bottom:10px; border-collapse: collapse; font-size:14px; color:gray; font-style:italic;'>
                  <thead>
                    <tr style='background-color:#f2f2f2;'>
                      <th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap;'>sample</th>
                      <th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap;'>acquisitionmode</th>
                      <th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap;'>sampletype</th>
                    </tr>
                  </thead>
                  <tbody>
                    <tr>
                      <td style='border: 1px solid #ddd; padding: 6px;'>Sample name</td>
                      <td style='border: 1px solid #ddd; padding: 6px;'>acquisitionmode<br>[DIA, DDA or MS]</td>
                      <td style='border: 1px solid #ddd; padding: 6px;'>sample type <br>[eg. Blank, Pool...]</td>
                    </tr>
                  </tbody>
                </table>

              </div>"
        )
      }
      
      # --- JavaScript toggle function ---
      text_output <- paste0(
        text_output,
        "<script>
      function toggleDetails(id) {
        var content = document.getElementById(id);
        var toggleIcon = document.getElementById(id + '_toggle');
        if (content.style.display === 'none') {
          content.style.display = 'block';
          toggleIcon.innerHTML = '△';
        } else {
          content.style.display = 'none';
          toggleIcon.innerHTML = '▽';
        }
      }
    </script>"
      )
    
    
    # Lipid classes (default or custom) ---
    if (!is.null(input$type_lipidclasses)) {
      
      if (input$type_lipidclasses == "lipidclasses_default") {
        text_type_lipidclasses <- "&nbsp;&nbsp; (default list)"
        lipid_adducts <- names(choicesNamesReactive())  # lipids default
        lipid_available <- length(lipid_adducts) > 0
        
      } else if (input$type_lipidclasses == "lipidclasses_custom") {
        text_type_lipidclasses <- "&nbsp;&nbsp; (custom list)"
        
        # If a custom file is loaded (e.g., classes_file), we display its classes
        if (!is.null(fileClassesReactive()) && length(choicesNamesReactive()) > 0) {
          lipid_adducts <- names(choicesNamesReactive())
          lipid_available <- TRUE
        } else {
          lipid_adducts <- character(0)
          lipid_available <- FALSE
        }
      }
      
      # If a list is available...
      if (lipid_available) {
        lipid_checkbox <- "<p><input type='checkbox' checked disabled><b> Lipid classes list for analysis</b></p>"
        
        lipid_text <- "<table style='border-collapse: collapse; margin-left: 30px; margin-top:10px;'>"
        for (lipid in lipid_adducts) {
          html_entry <- choicesNamesReactive()[[lipid]]
          
          # Extract name and adducts from the original HTML
          lipid_class <- sub(".*<b[^>]*>([^<]+)</b>.*", "\\1", html_entry)
          lipid_class <- gsub("&nbsp;", "", lipid_class)
          
          adducts_raw <- sub(".*<small>\\(([^<]+)\\)</small>.*", "\\1", html_entry)
          adducts_clean <- gsub("&nbsp;", "", adducts_raw)
          adduct_list <- unlist(strsplit(adducts_clean, ";"))
          adducts_text <- paste(trimws(adduct_list), collapse = paste0(" ", input$adducts_sep, " "))
          
          lipid_text <- paste0(
            lipid_text,
            "<tr>",
            "<td style='padding: 1px 15px 1px 0; color:#2c3e50; font-size: 85%; font-weight:bold;'>", lipid_class, "</td>",
            "<td style='padding: 1px 0; color:#666; font-size: 85%;'>", adducts_text, "</td>",
            "</tr>"
          )
        }
        lipid_text <- paste0(lipid_text, "</table>")
        
        unique_id <- paste0("lipid_details_", as.integer(runif(1, 1000, 9999)))

        # Expandable
        lipid_expandable <- paste0(
          "<hr>",
          lipid_checkbox,
          "<div style='margin-left:50px; margin-right:80px; border:1px solid #ddd; border-radius:6px; padding:5px; background-color:#f9f9f9; font-size:14px; color:#2C3E50;'>",
          "<div style='color:#2c3e50; font-weight:bold; cursor:pointer;' onclick='toggleDetails(\"", unique_id, "\")'>",
          "<span id='", unique_id, "_toggle'>△</span> Show lipid classes",
          "<span style='font-weight:normal;'>", text_type_lipidclasses, "</span></div>",
          "<div id='", unique_id, "' style='display:block; margin-top: 10px; color:#666;'>",
          lipid_text,
          "</div></div>"
        )
        
        text_output <- paste0(text_output, lipid_expandable)
        
        # Script JS toggle
        text_output <- paste0(
          text_output,
          "<script>
      function toggleDetails(id) {
        var content = document.getElementById(id);
        var toggleIcon = document.getElementById(id + '_toggle');
        if (content.style.display === 'none') {
          content.style.display = 'block';
          toggleIcon.innerHTML = '△';
        } else {
          content.style.display = 'none';
          toggleIcon.innerHTML = '▽';
        }
      }
    </script>"
        )
        
      } else {
        # If there is no list, display a message
        lipid_checkbox <- "<p><input type='checkbox' disabled><b> Lipid classes list for analysis</b></p>"
        
        if (input$type_lipidclasses == "lipidclasses_custom") {
          
            text_expandible <- "<p style='margin-left:20px; margin-top:10px; color:gray; font-style:italic; font-size:14px;'>
            No lipid classes available for prediction.<br>
            Please upload a valid <b>.csv</b> file with your custom lipid classes and adducts.</p>"

          text_output <- paste0(
            text_output,
            "<hr>",
            "<p><input type='checkbox' disabled><b> Lipid classes list for analysis</b></p>",
            
            "<div style='margin-left:50px; margin-right:80px; border:1px solid #ddd; border-radius:6px; padding:5px; background-color:#f9f9f9; font-size:14px; color:#2C3E50;'>",
            text_expandible,
            
            "<p style='margin-left:20px; color:gray; font-style:italic; font-size:14px;'>
                          The custom lipid classes file must contain the following columns:</p>

                          <table style='margin-left:50px; margin-bottom:10px; border-collapse: collapse; font-size:14px; color:gray; font-style:italic;'>
                  <thead>
                    <tr style='background-color:#f2f2f2;'>
                      <th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap;'>Classes</th>
                      <th style='border: 1px solid #ddd; padding: 6px; white-space: nowrap;'>Adducts</th>
                    </tr>
                  </thead>
                  <tbody>
                    <tr>
                      <td style='border: 1px solid #ddd; padding: 6px;'>Name of the<br>lipid</td>
                      <td style='border: 1px solid #ddd; padding: 6px;'>Ion type<br>(e.g., [M+H]+)</td>
                    </tr>
                  </tbody>
                </table>

              </div>"
          )
        } else {
          text_expandible <- "<div style='margin-left: 30px; color:gray; font-style: italic;'>No lipid classes available.</div>"
          text_output <- paste0(
            text_output,
            "<hr>",
            "<p><input type='checkbox' disabled><b> Lipid classes list for analysis</b></p>",
            text_expandible,
            "<div style='margin-left:50px; margin-right:80px; border:1px solid #ddd; border-radius:6px; padding:5px; background-color:#f9f9f9; font-size:14px; color:#2C3E50;'>
                  <p style='margin-left:20px; margin-top:10px; color:gray; font-style:italic; font-size:14px;'>
                  No lipid classes available because the polarity of the analysis has not been defined due to missing data.</p>
                  <p style='margin-left:20px; color:gray; font-style:italic; font-size:14px;'>
                  Please upload a valid file (<b>.csv</b>, <b>.xlsx</b>, or <b>.txt</b>) or import the data directly from the R environment, depending on your preferred method of data input.</p>


              </div>"
          )
        }
      }
    }
    
    HTML(text_output)
  })
  
  
  
  #============================================================================#
  # Change reactive values when there are changes in...
  #============================================================================#
  # data_type, msdata, how to process data, polarity
  observe({
      polarityReactive(input$to_process_polarity) 
      analysisReactive(input$to_process_analysis)
  })
  
  
  #============================================================================#
  # Change lipid class list according to default or custom                    
  #============================================================================#
  
  observeEvent(
    list(input$type_lipidclasses, fileClassesReactive(), polarityReactive()), {
    if (input$type_lipidclasses == "lipidclasses_default") {
      # Forzar el reinicio del fileInput cambiando el label
      # Forzar re-render del fileInput
      output$file_upload_ui <- renderUI({
        fileInput("classes_file", NULL,
                  multiple = FALSE,
                  accept = c("text/plain", ".csv"),
                  width = "100%",
                  placeholder = "No file selected")
      })
      
      
      # También puedes actualizar otros reactivos si es necesario
      polarity <- polarityReactive()
      dbs <- LipidMS:::dbsAdducts(polarity = polarity)
      fileClassesReactive(NULL)
      
      select_lipidclasses(analysis_step = "lipidclasses_annotate", dbs)
      
      # si es lista custom...
    } else if (input$type_lipidclasses == "lipidclasses_custom") {
        if (is.null(fileClassesReactive())) { # si no se ha subido archivo lista vacia
          select_lipidclasses(analysis_step = "lipidclasses_annotate", dbs = NULL)
        } else { # si se ha subido leer archivo
          file_path <- fileClassesReactive()
          dbs <- LipidMS:::dbsAdducts(
            file = file_path,
            sep = input$classes_sep,
            sep_adducts = input$adducts_sep,
            polarity = NULL,
            classes = NULL
          )
          select_lipidclasses(analysis_step = "lipidclasses_annotate", dbs = dbs)
        }
    }
  })

  
  #============================================================================#
  # Import a .RData file.
  #============================================================================#
  observeEvent(input$RDatafile, {
    
    # Load notification.
    id_loading <- showNotification(
      HTML(paste0("Importing .RData file.<br>
      It may take a few minutes...<span class='loader'></span>")),
      type = "default",
      duration = NULL
    )

    # Create a temporary environment and load the file.
    temp_env <- new.env()
    
    tryCatch({
      load(input$RDatafile$datapath, envir = temp_env)  # Load the file into the temporary environment.
      
      msdata_name <- ls(temp_env)  # Get the name of the loaded object.
      msdata <- temp_env[[msdata_name]]  # Extract the object.
      
      # Check if it has the correct structure.
      msdata_structure <- LipidMS:::checkmsdatastructure(msdata)
      analysis_type <- msdata_structure$structure  # Save the type of analysis.
      
      # If it has the correct structure...
      if (analysis_type %in% c("msbatch", "msobject")) {
        analysisReactive(analysis_type) # Save the type of analysis.
        
        msdataReactive(msdata)  # Save the object in reactive.
        
        # Get the polarity of the loaded object.
        polarity <- LipidMS:::determinepolarity(msdata, analysis_type)
        polarityReactive(polarity)
        
        if (is.null(polarity) || !(analysis_type %in% c("msbatch", "msobject"))) {
          showNotification(
            paste0("The object '", msdata_name, "' does not have the required format."),
            type = "error"
          )
          removeNotification(id_loading)
          return()
        }
        
        # Update the interface with the results.
        update_summaryTable_ui(session)
        update_annotatedPeaklist_ui(session)
        update_features_ui(session)
        
        # Notification that it has been imported successfully.
        showNotification(
          paste0("The object '", msdata_name, "' was imported successfully."),
          type = "message",
          duration = NULL
        ) 
      }
      
    }, error = function(e) {
      showNotification(
        paste("Error importing the file:", e$message),
        type = "error"
      )
    }, finally = {
      removeNotification(id_loading)
    })
  })
  
  
  #============================================================================#
  # Import msbatch from the global environment into the app.
  #============================================================================#
  observeEvent(input$import_msdata, {

    msdata_name <- input$import_name_msdata  # Get the name of the loaded object.
    
    id_loading <- showNotification(
      HTML(paste0("Importing the object ", msdata_name, " from the global environment... <span class='loader'></span>")),
      type = "default",
      duration = NULL
    )
    
    # Check if the object exists in the global environment.
    if (!exists(msdata_name, envir = .GlobalEnv)) {
      showNotification(
        paste0("The object '", msdata_name, "' does not exist in the global environment."),
        type = "error",
        duration = NULL
      )
      removeNotification(id_loading)
      return()
    }
    
    # Get the object from the global environment.
    msdata <- get(msdata_name, envir = .GlobalEnv)
    
    # Check if it has the correct structure.
    msdata_structure <- LipidMS:::checkmsdatastructure(msdata)
    analysis_type <- msdata_structure$structure  # Save the type of analysis.
    
    # If it has the correct structure...
    if (analysis_type %in% c("msbatch", "msobject")) {
      msdataReactive(msdata)
      analysisReactive(analysis_type)
            
      # Get the polarity of the loaded object.
      polarity <- LipidMS:::determinepolarity(msdata, analysis_type)
      polarityReactive(polarity)
      
      if (is.null(polarity) || !(analysis_type %in% c("msbatch", "msobject"))) {
        showNotification(
          paste0("The object '", msdata_name, "' does not have the required format."),
          type = "error"
        )
        removeNotification(id_loading)
        return()
      }
      
      # Update the interface with the results.
      update_summaryTable_ui(session)
      update_annotatedPeaklist_ui(session)
      update_features_ui(session)
      
      # Notification that it has been imported successfully.
      showNotification(
        paste0("The object '", msdata_name, "' was imported successfully."),
        type = "message",
        duration = NULL
      )
      
      # Remove the initial load notification.
      removeNotification(id_loading)
      
    } else {
      # If it doesn’t have the correct structure, it notifies you and stops.
      showNotification(
        paste0("The object '", msdata_name, "' does not have the required format."),
        type = "error"
      )
      removeNotification(id_loading)
      return()
    }

  })
  
  
  #============================================================================#
  # Load file with the classes and adducts to predict.
  #============================================================================#
  observeEvent(input$classes_file, {
    # Show loading message
    loading_id <- showNotification(HTML("Importing lipid class file... <span class='loader'></span>"),
                                   type = "default", duration = NULL)
    
    # Validate file input
    if (is.null(input$classes_file) || input$classes_file$datapath == "") {
      removeNotification(loading_id)
      showNotification("No file uploaded", type = "error")
    } else {
      # Safely attempt to run select_lipidclasses
      tryCatch({
        # Leer el archivo directamente desde input
        file_path <- input$classes_file$datapath
        
        fileClassesReactive(file_path)
        
        removeNotification(loading_id)
        showNotification("Lipid class file imported successfully", type = "message", duration = NULL)
        
      }, error = function(e) {
        removeNotification(loading_id)
        showNotification(paste("Error importing lipid class file:", e$message), type = "error", duration = NULL)
      })
    }
  })
  

  #============================================================================#
  # Export msbatch from the global environment to the app.
  #============================================================================#
  observeEvent(input$export_msdata, {
    
    msdata <- msdataReactive()
    analysis_type <- analysisReactive()
    
    # Keep the global environment.
    envir <- ".GlobalEnv"
    if (envir == ".GlobalEnv") {
      envir <- tryCatch({
        get(envir, envir = .GlobalEnv)
      }, error = function(e) {
        NULL
      })
    } else {
      stop("Wrong envir. It must be '.GlobalEnv'.")
    }
    
    # Send the msdata to the global environment.
    if (!is.null(msdata)) {
      if (analysis_type == "msobject") { # If it is an msobject…
        id_loading <- showNotification(HTML("Exporting the object 'msobject' to the global environment... <span class='loader'></span>"), 
                                       type = "default", duration = NULL)
        assign("msobject", msdata, envir = envir) # Default name "msobject".
        
        removeNotification(id_loading)
        showNotification("The object 'msobject' has been successfully exported.", type = "message", duration = NULL) 
        
      } else if (analysis_type == "msbatch") { # If it is an msbatch…
        id_loading <- showNotification(HTML("Exporting the object 'msbatch' to the global environment... <span class='loader'></span>"),
                                       type = "default", duration = NULL)
        assign("msbatch", msdata, envir = envir) # Default name "msbatch".
        
        removeNotification(id_loading)
        showNotification("The object 'msbatch' has been successfully exported.", type = "message", duration = NULL) 
      } 
      
    } else {
      # SIf the object doesn’t exist, nothing will be exported.
      if (analysis_type == "msobject") {
        showNotification("The object 'msobject' does not exist.", type = "error", duration = NULL)
      } else if (analysis_type == "msbatch") {
        showNotification("The object 'msbatch' does not exist.", type = "error", duration = NULL) 
      } 
    }

  })


  # Make the download button appear only if data is loaded.
  output$download_summaryTable_ui <- renderUI({
    analysis_type <- analysisReactive()
    msdata <- msdataReactive()
    
    # Recalcular objetos válidos
    valid_objects <- list()
    if (analysis_type == "msobject" && !is.null(msdata)) {
      valid_objects <- Filter(function(obj) !is.null(obj$annotation$results), msdata)
    } 
    
    has_valid_data <- length(valid_objects) > 0
    
    # Mostrar el botón solo si hay datos válidos
    if (has_valid_data) {
      downloadButton("downloadSummary", "Summary Tables", 
                     class = "btn-primary",
                     style = "color: #fff; background-color: #435f7a ;
                              border-color: #2c3e50 ; width: 100%;
                              text-align: center; white-space: normal;
                              margin: 10px 0px 5px 0px;")
    } else {
      NULL
    }
  })
  
  
  output$download_annotatedPeaklist_ui <- renderUI({
    analysis_type <- analysisReactive()
    msdata <- msdataReactive()
    
    # Recalcular objetos válidos
    valid_objects <- list()
    if (analysis_type == "msobject" && !is.null(msdata)) {
      valid_objects <- Filter(function(obj) !is.null(obj$annotation$annotatedPeaklist), msdata)
    } 
    
    has_valid_data <- length(valid_objects) > 0
    
    if (has_valid_data) {
      downloadButton("downloadPeaklist", "Annotated Peaklists", 
                     class = "btn-primary",
                     style = "color: #fff; background-color: #435f7a ;
                              border-color: #2c3e50 ; width: 100%;
                              text-align: center; white-space: normal;
                              margin: 10px 0px 5px 0px;")
    } else {
      NULL
    }
  })
  
  
  output$download_batchResults_ui <- renderUI({
    analysis_type <- analysisReactive()
    msdata <- msdataReactive()
    
    # Solo mostrar si hay datos válidos
    has_valid_data <- !is.null(msdata$features) && nrow(msdata$features) > 0
    
    if (analysis_type == "msbatch" && has_valid_data) {
      tags$div(
        style = "display: flex; justify-content: flex-end; gap: 10px;",
        downloadButton("downloadBatchResults", "Batch Results", class = "btn btn-primary"),
        downloadButton("downloadPlots", "Batch Plots", 
                       class = "btn-primary",
                       style = "color: #fff; background-color: #435f7a ;
                              border-color: #2c3e50 ; width: 100%;
                              text-align: center; white-space: normal;
                              margin: 10px 0px 5px 0px;")
      )
    } else {
      NULL
    }
  })
  
  
  # Make the download button appear only if data is loaded.
  output$download_msdata_ui <- renderUI({
    if (!is.null(msdataReactive())) {
      downloadButton("download_msdata", ".RData", class = "btn-primary", 
                     style = "color: #fff; background-color: #2c3e50 ;
                              border-color: #2c3e50 ; width: 100%;
                              text-align: center; white-space: normal;")
    }
  })
  
  # Make the export button appear only if data is loaded.
  output$export_msdata_ui <- renderUI({
    if (!is.null(msdataReactive())) {
      actionButton("export_msdata",  "Export in R environment", 
                   style = "color: #fff; background-color: #2c3e50 ;
                              border-color: #2c3e50 ; width: 100%;
                              text-align: center; white-space: normal;")
    }
  })
  
  
  #============================================================================#
  # Run Button
  #============================================================================#
  observeEvent(input$Run_analysis, {
    analysis_type <- analysisReactive()
    
    if (analysis_type == "msobject" | analysis_type == "single") { # MSOBJECTS ("single") ## ## ## ## ## ## 
      
      # apply annotation (single) ............................................
          req(input$file1, input$metadata)
          
          id_loading <- showNotification(HTML("Processing data... <span class='loader'></span>"),
                                         type = "default", duration = NULL)
          
          # Metadata
          metadata <- utils::read.csv(input$metadata$datapath, sep=",")
          if (!all(grepl(".mzXML", metadata$sample))){
            metadata$sample <- paste(metadata$sample, ".mzXML", sep = "")
          }
          
          # .mzXML files
          files <- data.frame(sample = input$file1$name,
                              path = input$file1$datapath)
          
          metadata <- as.data.frame(merge(metadata, files, by = "sample"))
          
          # Keep only the lipid classes to annotate.
          if (input$to_process_polarity == "positive") {
            lipidClassesPos <- input$lipidclasses_annotate
            lipidClassesNeg <- c() 
          } else if (input$to_process_polarity == "negative") {
            lipidClassesPos <- c()
            lipidClassesNeg <- input$lipidclasses_annotate 
          }
          
          print("Run")
          print(lipidClassesNeg)
          # Annotation
          msdata <- singleProcessing(files = metadata$path, 
                                     filesname= metadata$sample,
                                     acquisitionmode = metadata$acquisitionmode, 
                                     polarity = input$to_process_polarity,
                                     input$dmzagglom_ms1,input$dmzagglom_ms2,
                                     input$drtagglom_ms1, input$drtagglom_ms2,
                                     input$drtclust_ms1, input$drtclust_ms2,
                                     input$minpeak_ms1, input$minpeak_ms2,
                                     input$drtgap_ms1, input$drtgap_ms2,
                                     input$drtminpeak_ms1, input$drtminpeak_ms2,
                                     input$drtmaxpeak_ms1, input$drtmaxpeak_ms2,
                                     input$recurs_ms1, input$recurs_ms2,
                                     input$sb_ms1, input$sb_ms2,
                                     input$sn_ms1, input$sn_ms2,
                                     input$minint_ms1, input$minint_ms2,
                                     input$weight_ms1, input$weight_ms2,
                                     input$dmzIso_ms1, input$dmzIso_ms2,
                                     input$drtIso_ms1, input$drtIso_ms2,
                                     input$dmzprecursor, input$dmzproducts,
                                     input$rttol, input$coelcutoff,
                                     input$jobname,
                                     lipidClassesPos,
                                     lipidClassesNeg)
          msdataReactive(msdata) 
            
          # Update the user interface.
          update_summaryTable_ui(session)
          update_annotatedPeaklist_ui(session)
          update_features_ui(session)
          
          # Notification that the annotation has been completed.
          removeNotification(id_loading)
          showNotification("Processing has been successfully completed.", type = "message", duration = NULL)
    } 
    
    else if (analysis_type == "msbatch" | analysis_type == "batch") { # MSBATCH ## ## ## ## ## ## ## #
      
      # apply annotation (msbatch) .............................................

          req(input$file1, input$metadata)
  
          id_loading <- showNotification(HTML("Processing data... <span class='loader'></span>"),
                                         type = "default", duration = NULL)
          
          # Metadata
          metadata <- utils::read.csv(input$metadata$datapath, sep=",")
          
          if (!all(grepl(".mzXML", metadata$sample))){
            metadata$sample <- paste(metadata$sample, ".mzXML", sep = "")
          }
          
          # .mzXML file
          files <- data.frame(sample = input$file1$name,
                              path = input$file1$datapath)
          
          metadata <- merge(metadata, files, by = "sample")
          
          
          # Keep only the lipid classes to annotate.
          if (input$to_process_polarity == "positive") {
            lipidClassesPos <- input$lipidclasses_annotate
            lipidClassesNeg <- c()
          } else if (input$to_process_polarity == "negative") {
            lipidClassesPos <- c()
            lipidClassesNeg <- input$lipidclasses_annotate 
          }
          
          # Annotation
          msdata <- batchProcessing(metadata = metadata,
                                    input$to_process_polarity,
                                    input$dmzagglom_ms1,input$dmzagglom_ms2,
                                    input$drtagglom_ms1, input$drtagglom_ms2,
                                    input$drtclust_ms1, input$drtclust_ms2,
                                    input$minpeak_ms1, input$minpeak_ms2,
                                    input$drtgap_ms1, input$drtgap_ms2,
                                    input$drtminpeak_ms1, input$drtminpeak_ms2,
                                    input$drtmaxpeak_ms1, input$drtmaxpeak_ms2,
                                    input$recurs_ms1, input$recurs_ms2,
                                    input$sb_ms1, input$sb_ms2,
                                    input$sn_ms1, input$sn_ms2,
                                    input$minint_ms1, input$minint_ms2,
                                    input$weight_ms1, input$weight_ms2,
                                    input$dmzIso_ms1, input$dmzIso_ms2,
                                    input$drtIso_ms1, input$drtIso_ms2,
                                    input$dmzalign, input$drtalign,
                                    input$span, input$minsamplesfracalign,
                                    input$dmzgroup, input$drtagglomgroup,
                                    input$drtgroup, input$minsamplesfracgroup,
                                    input$parallel, 
                                    input$ncores, input$global_gb,
                                    input$dmzprecursor, input$dmzproducts,
                                    input$rttol, input$coelcutoff,
                                    input$jobname,
                                    lipidClassesPos,
                                    lipidClassesNeg)
          
          msdataReactive(msdata)
          
          # Update the user interface.
          update_summaryTable_ui(session)
          update_annotatedPeaklist_ui(session)
          update_features_ui(session)
          
          # Notification that the prediction has been completed.
          removeNotification(id_loading)
          showNotification("Processing has been successfully completed.", type = "message", duration = NULL)
    }
    
  })
  
  

  #============================================================================#
  # -- "SummaryTable Results" interface.
  #============================================================================#
  update_summaryTable_ui <- function(session) {

    analysis_type <- analysisReactive()
    msdata <- msdataReactive()

    # Filter only the objects that have information in annotation$results.
    valid_objects <- list()

    if (analysis_type == "msobject" && !is.null(msdata)) {
      valid_objects <- Filter(function(obj) !is.null(obj$annotation$results), msdata)
    } else if (analysis_type == "msbatch" && !is.null(msdata$msobjects)) {
      valid_objects <- Filter(function(obj) !is.null(obj$annotation$results), msdata$msobjects)
    }

    has_valid_data <- length(valid_objects) > 0

    # Render dynamic UI.
    output$summaryTable_ui <- renderUI({
      if (!has_valid_data) {
        # If there is no valid information, display a message.
        tags$div(
          style = "margin-top: 20px; font-size: 18px; font-weight: bold; text-align: center;",
          "There is no information available to display the table.")
      } else {
        # If there is valid information
        tags$div(
          style = "margin-top: 20px;",
          tags$div(
            # Object selector to view its summary.
            style = "width: 100%;",
            selectInput(
              inputId = "selected_summary",
              label = "Select an object to view its summary table:",
              choices = sapply(valid_objects, function(obj) obj$metaData$generalMetadata$file),
              selected = valid_objects[[1]]$metaData$generalMetadata$file,
              width = "100%"
              )
          ),
          DT::DTOutput("summaryTable") # Table
        )
      }
    })

    # Render the interactive table only if there is valid data.
    if (has_valid_data) {
      output$summaryTable <- DT::renderDT({
        selected_object <- valid_objects[[which(sapply(valid_objects,
                                                       function(obj) obj$metaData$generalMetadata$file) == input$selected_summary)]]
        # Table of the selected object.
        DT::datatable(
          selected_object$annotation$results,
          options = list(
            pageLength = 10,
            autoWidth = TRUE,
            scrollX = TRUE,
            dom = 'frtip',
            paging = TRUE,
            pageLength = 10,
            pagingType = "full_numbers"
          )
        )
      })
    }
  }

  
  #============================================================================#
  # -- "AnnotatedPeaklist Results" interface.
  #============================================================================#
  update_annotatedPeaklist_ui <- function(session) {
    
    analysis_type <- analysisReactive()
    msdata <- msdataReactive()
    
    # Filter only the objects that have information in annotation$annotatedPeaklist.
    valid_objects <- list()
    
    if (analysis_type == "msobject" && !is.null(msdata)) {
      valid_objects <- Filter(function(obj) !is.null(obj$annotation$annotatedPeaklist), msdata)
    } else if (analysis_type == "msbatch" && !is.null(msdata$msobjects)) {
      valid_objects <- Filter(function(obj) !is.null(obj$annotation$annotatedPeaklist), msdata$msobjects)
    }
    
    has_valid_data <- length(valid_objects) > 0
    
    # Render dynamic UI.
    output$annotatedPeaklist_ui <- renderUI({
      if (!has_valid_data) {
        # If there is no valid information, display a message.
        tags$div(
          style = "margin-top: 20px; font-size: 18px; font-weight: bold; text-align: center;",
          "There is no information available to display the annotated peak table."
        )
      } else {
        # If there is valid information
        if (analysis_type == "msobject") { # If it is a msobject...
          tags$div(
            style = "margin-top: 20px;",
            # Object selector to view its summary.
            tags$div(
              style = "width: 100%;",
              selectInput(
                inputId = "selected_msobject", 
                label = "Select a file to view its table:",
                choices = sapply(valid_objects, function(obj) obj$metaData$generalMetadata$file),
                selected = valid_objects[[1]]$metaData$generalMetadata$file,
                width = "100%"
              )
            ),
            DT::DTOutput("annotatedPeaklist_table") # Table
          )
        } else if (analysis_type == "msbatch") { # If it is a msobject...
          tags$div(
            style = "margin-top: 20px;",
            # Object selector to view its summary.
            tags$div(
              style = "width: 100%;",
              selectInput(
                inputId = "selected_file", 
                label = "Select a file to view its table:",
                choices = sapply(valid_objects, function(obj) obj$metaData$generalMetadata$file),
                selected = valid_objects[[1]]$metaData$generalMetadata$file,
                width = "100%"
              )
            ),
            DT::DTOutput("annotatedPeaklist_msbatch_table") # Table
          )
        }
      }
    })
    
    # Render the interactive table only if there is valid data (msobject).
    if (analysis_type == "msobject" && has_valid_data) {
      output$annotatedPeaklist_table <- DT::renderDT({
        selected_object <- valid_objects[[which(sapply(valid_objects, 
                                                       function(obj) obj$metaData$generalMetadata$file) == input$selected_msobject)]]
        
        DT::datatable(
          selected_object$annotation$annotatedPeaklist,
          options = list(
            pageLength = 10,
            autoWidth = TRUE,
            scrollX = TRUE,
            dom = 'frtip',
            paging = TRUE,
            pageLength = 10,
            pagingType = "full_numbers"
          )
        )
      })
    }
    
    # Render the interactive table only if there is valid data (msbatch).
    if (analysis_type == "msbatch" && has_valid_data) {
      output$annotatedPeaklist_msbatch_table <- DT::renderDT({
        selected_object <- valid_objects[[which(sapply(valid_objects, 
                                                       function(obj) obj$metaData$generalMetadata$file) == input$selected_file)]]
        
        DT::datatable(
          selected_object$annotation$annotatedPeaklist,
          options = list(
            pageLength = 10,
            autoWidth = TRUE,
            scrollX = TRUE,
            dom = 'frtip',
            paging = TRUE,
            pageLength = 10,
            pagingType = "full_numbers"
          )
        )
      })
    }
  }
  
  
  #============================================================================#
  # -- "Batch Results" interface.
  #============================================================================#
  update_features_ui <- function(session) {
    
    analysis_type <- analysisReactive()
    msdata <- msdataReactive()
    
    if (analysis_type == "msobject") { # If it is an msobject, it should be empty.
      output$features_ui <- renderUI({
        tagList(
          tags$div(
            style = "margin-top: 20px;")
        )
      })
    } else if (analysis_type %in% c("msbatch", "external")) { # If it is a msbatch o external...
      # Render dynamic content in the UI.
      if (analysis_type == "msbatch") { # If it is an msbatch…
        output$features_ui <- renderUI({
          tagList(
            tags$div(
              style = "margin-top: 20px;",
              DT::DTOutput("features_table")) # Table
          )
        }) 
      } 
      # Render the interactive table with the data from msdata$features
      output$features_table <- DT::renderDT({
        
        DT::datatable(
          msdata$features,  
          options = list(
            pageLength = 10,
            autoWidth = TRUE,  # Automatic column width adjustment
            scrollX = TRUE,  # Horizontal scroll
            dom = 'frtip'  # Table elements (buttons, search, etc.)
          ),
          extensions = 'Buttons',  # Enable export buttons
          class = "display nowrap"  # Prevent the table from overflowing outside the container
        )
      })

    }
  }
  
  

  #============================================================================#
  # Downloads
  #============================================================================#
  # Download SUMMARY TABLES RESULTS ...........................................
  output$downloadSummary <- downloadHandler(
    filename = function(){paste(input$jobname, "SummaryTables.zip", sep="_")},
    content = function(file){
      
      notif_id <- showNotification(HTML("Downloading 'Summary'... <span class='loader'></span>"),
                                   type = "default", duration = NULL)
      
      # go to a temp dir to avoid permission issues
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      files <- NULL;
      
      msdata <- msdataReactive()
      analysis_type <- analysisReactive()
      
      if (analysis_type == "msobject") {
        # loop through the sheets
        for (i in 1:length(msdata)){
          #write each sheet to a csv file, save the name
          fileName <- paste(gsub(".mzXML", "" , input$file1$name[i]), "_summaryResults.csv", sep="")
          write.csv(msdata[[i]]$annotation$results, fileName, row.names = FALSE)
          files <- c(fileName, files)
        }
      } else if (analysis_type == "msbatch") {
        # loop through the sheets
        for (i in 1:length(msdata$msobjects)){
          if (msdata$msobjects[[i]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
            # write each sheet to a csv file, save the name
            fileName <- gsub(".mzXML", "_summaryResults.csv" , msdata$metaData$sample[i])
            write.csv(msdata$msobjects[[i]]$annotation$results, fileName, row.names = FALSE)
            files <- c(fileName, files)
          }
        }
      }
      
      # create the zip file
      zip(file, files)
      
      removeNotification(notif_id)
      showNotification("'Summary' has been successfully downloaded.", type = "message", duration = NULL)
    }
  )
  
  # Download ANNOTATED PEAKLIST RESULTS .......................................
  output$downloadPeaklist <- downloadHandler(
    filename = function(){paste(input$jobname, "Peaklists.zip", sep="_")},
    content = function(file){
      notif_id <- showNotification(HTML("Downloading 'Peaklists'... <span class='loader'></span>"),
                                   type = "default", duration = NULL)
      # go to a temp dir to avoid permission issues
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      files <- NULL;
      
      msdata <- msdataReactive()
      analysis_type <- analysisReactive()
      
      if (analysis_type == "msobject") {
        # loop through the sheets
        for (i in 1:length(msdata)){
          # write each sheet to a csv file, save the name
          fileName <- paste(gsub(".mzXML", "" , input$file1$name[i]), "_annotatedPeaklist.csv", sep="")
          write.csv(msdata[[i]]$annotation$annotatedPeaklist, fileName, row.names = FALSE)
          files <- c(fileName, files)
        }
      } else if (analysis_type == "msbatch") {
        # loop through the sheets
        for (i in 1:length(msdata$msobjects)){
          if (msdata$msobjects[[i]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
            # write each sheet to a csv file, save the name
            fileName <- gsub(".mzXML", "_annotatedPeaklist.csv" , msdata$metaData$sample[i])
            write.csv(msdata$msobjects[[i]]$results$annotatedPeaklist, fileName, row.names = FALSE)
            files <- c(fileName, files)
          }
        }
      }
      
      # create the zip file
      zip(file, files)
      
      removeNotification(notif_id)
      showNotification("'Peaklists' has been successfully downloaded.", type = "message", duration = NULL)
    }
  )
  
  # Download PLOTS RESULTS .....................................................
  output$downloadPlots <- downloadHandler(
    filename = function(){paste(input$jobname, "Plots.zip", sep="_")},
    content = function(file){
      
      analysis_type <- analysisReactive()
      
      if (analysis_type == "msbatch") {
        notif_id <- showNotification(HTML("Downloading 'Plots'... <span class='loader'></span>"),
                                     type = "default", duration = NULL)
        
        msdata <- msdataReactive()
        
        # go to a temp dir to avoid permission issues
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        files <- NULL;
        
        # loop through the msobjects
        for (i in 1:length(msdata$msobjects)){
          if (msdata$msobjects[[i]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
            fileName <- gsub(".mzXML", "_plots.pdf" , msdata$metaData$sample[i])
            if (msdata$msobjects[[i]]$metaData$generalMetadata$acquisitionmode == "DIA"){
              height <- 7
            } else {
              height <- 8
            }
            grDevices::pdf(file = fileName, width = 8, height = height)
            if (length(msdata$msobjects[[i]]$annotation$plots) > 0){
              for (pl in 1:length(msdata$msobjects[[i]]$annotation$plots)){
                print(msdata$msobjects[[i]]$annotation$plots[[pl]])
              }
            }
            grDevices::dev.off()
            files <- c(fileName, files)
          }
        }
        # create the zip file
        zip(file, files) 
        
        removeNotification(notif_id)
        showNotification("'Plots' has been successfully downloaded.", type = "message", duration = NULL)
      }
    }
  )
  
  # Download BATCH RESULTS .....................................................
  output$downloadBatchResults <- downloadHandler(
    filename = function(){paste(input$jobname, "BatchResults.zip", sep="_")},
    content = function(file){
      
      analysis_type <- analysisReactive()
      
      if (analysis_type == "msbatch") {
        notif_id <- showNotification(HTML("Downloading 'Batch Results'... <span class='loader'></span>"),
                                     type = "default", duration = NULL) 
        
        msdata <- msdataReactive()
        
        # go to a temp dir to avoid permission issues
        owd <- setwd(tempdir())
        on.exit(setwd(owd));
        
        # extract data
        peaklist <- msdata$features
        peaklistNoIso <- peaklist[peaklist$isotope %in% c("", "[M+0]"),]
        
        # write files
        write.csv(peaklist, "FeaturesMatrix.csv", row.names = FALSE)
        write.csv(peaklistNoIso, "FeaturesMatrixIsotopesRemoved.csv", row.names = FALSE)
        pdf("RTdevplot.pdf")
        LipidMS::rtdevplot(msdata)
        LipidMS::rtdevplot(msdata, colorbygroup = FALSE)
        dev.off()
        pdf("TIC.pdf", height = 7, width = 10)
        LipidMS::plotticmsbatch(msdata)
        LipidMS::plotticmsbatch(msdata, colorbygroup = FALSE)
        dev.off()
        
        
        files <- c("FeaturesMatrix.csv", "FeaturesMatrixIsotopesRemoved.csv",
                   "RTdevplot.pdf", "TIC.pdf")
        
        #create the zip file
        zip(file, files) 
        
        removeNotification(notif_id)
        showNotification("'Batch Results' has been successfully downloaded.", type = "message", duration = NULL)
        
      } 
    
    }
  )
  

  
  # Download .RDATA ...........................................................
  output$download_msdata <- downloadHandler(
    filename = function() {
      msdata <- msdataReactive()
      analysis_type <- analysisReactive()
      
      if (is.null(msdata)) {
        return(NULL)  # If there is no data, skip the download.
      } else {
        paste(input$jobname, paste(analysis_type,"RData", sep = "."), sep = "_") 
      }
    },
    content = function(file) {
      msdata <- msdataReactive()
      analysis_type <- analysisReactive()
      
      if (is.null(msdata)) {
        return()  # If there is no data, do nothing.
      } else {
        if (analysis_type == "msobject") { # name 'msobject'
          msobject <- msdata
          save(msobject, file = file)
        } else if (analysis_type == "msbatch") { # name 'msbatch'
          msbatch <- msdata
          save(msbatch, file = file)
        } 
      }
    }
  )

  
  
  #============================================================================#
  # Annotation Functions
  #============================================================================#
  # BATCH PROCESSING ...........................................................
  batchProcessing <- function(metadata, polarity,
                              dmzagglom_ms1, dmzagglom_ms2, drtagglom_ms1, drtagglom_ms2,
                              drtclust_ms1, drtclust_ms2, minpeak_ms1, minpeak_ms2,
                              drtgap_ms1, drtgap_ms2, drtminpeak_ms1, drtminpeak_ms2,
                              drtmaxpeak_ms1, drtmaxpeak_ms2, recurs_ms1, recurs_ms2,
                              sb_ms1, sb_ms2, sn_ms1, sn_ms2, minint_ms1, minint_ms2,
                              weight_ms1, weight_ms2, dmzIso_ms1, dmzIso_ms2, drtIso_ms1,
                              drtIso_ms2, dmzalign, drtalign, span, minsamplesfracalign, 
                              dmzgroup, drtagglomgroup, drtgroup, minsamplesfracgroup,
                              parallel, 
                              ncores, global_gb, 
                              dmzprecursor, dmzproducts, rttol, coelcutoff,
                              jobname, lipidClassesPos, lipidClassesNeg){
    
    samplenames <- metadata$sample
    metadata$sample <- metadata$path
    
    # Convierte explícitamente el valor de input$parallel a lógico
    if (is.null(parallel) || !parallel %in% c("TRUE", "FALSE")) {
      parallel <- FALSE
    } else {
      parallel <- as.logical(parallel)
    }
    
    # Detectar ncores
    if (parallel) {
      ncores <- ncores
    } else {
      ncores <- 1
    }
    
    
    msbatch <- batchdataProcessing(files = metadata$sample,
                                   metadata = metadata,
                                   polarity = polarity,
                                   dmzagglom = c(dmzagglom_ms1, dmzagglom_ms2),
                                   drtagglom = c(drtagglom_ms1, drtagglom_ms2),
                                   drtclust = c(drtclust_ms1, drtclust_ms2),
                                   minpeak = c(minpeak_ms1, minpeak_ms2),
                                   drtgap = c(drtgap_ms1, drtgap_ms2),
                                   drtminpeak = c(drtminpeak_ms1, drtminpeak_ms2),
                                   drtmaxpeak = c(drtmaxpeak_ms1, drtmaxpeak_ms2),
                                   recurs = c(recurs_ms1, recurs_ms2),
                                   sb = c(sb_ms1, sb_ms2),
                                   sn = c(sn_ms1, sn_ms2),
                                   minint = c(minint_ms1, minint_ms2),
                                   weight = c(weight_ms1, weight_ms2),
                                   dmzIso = c(dmzIso_ms1, dmzIso_ms2),
                                   drtIso = c(drtIso_ms1, drtIso_ms2),
                                   parallel = parallel,
                                   ncores = ncores,
                                   global_gb = getOption("LipidMS.future.globals.maxSizeGB", global_gb))
    
    msbatch$metaData$sample <- samplenames
    
    # Alignment
    msbatch <- alignmsbatch(msbatch, dmz = dmzalign, drt = drtalign, span = span, 
                            minsamplesfrac = minsamplesfracalign, 
                            parallel = parallel, ncores = ncores,
                            global_gb = getOption("LipidMS.future.globals.maxSizeGB", global_gb))
    
    # Grouping
    msbatch <- groupmsbatch(msbatch, dmz = dmzgroup, drtagglom = drtagglomgroup,
                            drt = drtgroup, minsamplesfrac = minsamplesfracgroup, 
                            parallel = parallel, ncores = ncores,
                            global_gb = getOption("LipidMS.future.globals.maxSizeGB", global_gb))
    
    # Fill missing peaks
    msbatch <- fillpeaksmsbatch(msbatch)
    
    # Lipid Annotation
    msbatch <- annotatemsbatch(msbatch,
                               ppm_precursor = dmzprecursor,
                               ppm_products = dmzproducts,
                               rttol = rttol,
                               coelCutoff = coelcutoff,
                               lipidClassesPos = lipidClassesPos,
                               lipidClassesNeg = lipidClassesNeg,
                               parallel = FALSE)
    
    for (m in 1:length(msbatch$msobjects)){
      if (msbatch$msobjects[[m]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
        msbatch$msobjects[[m]] <- plotLipids(msbatch$msobjects[[m]])
      }
    }
    
    return(msbatch)
  }
  
  # SINGLE PROCESSING ..........................................................
  singleProcessing <- function(files, filesname, acquisitionmode, polarity,
                               dmzagglom_ms1, dmzagglom_ms2, drtagglom_ms1, drtagglom_ms2,
                               drtclust_ms1, drtclust_ms2, minpeak_ms1, minpeak_ms2,
                               drtgap_ms1, drtgap_ms2, drtminpeak_ms1, drtminpeak_ms2,
                               drtmaxpeak_ms1, drtmaxpeak_ms2, recurs_ms1, recurs_ms2,
                               sb_ms1, sb_ms2, sn_ms1, sn_ms2, minint_ms1, minint_ms2,
                               weight_ms1, weight_ms2, dmzIso_ms1, dmzIso_ms2, drtIso_ms1,
                               drtIso_ms2, dmzprecursor, dmzproducts, rttol, coelcutoff,
                               jobname, lipidClassesPos, lipidClassesNeg){
    
    msobjects <- list()
    
    for (f in 1:length(files)){
      msobjects[[f]] <- dataProcessing(file = files[f],
                                       polarity = polarity,
                                       acquisitionmode = acquisitionmode[f],
                                       dmzagglom = c(dmzagglom_ms1, dmzagglom_ms2),
                                       drtagglom = c(drtagglom_ms1, drtagglom_ms2),
                                       drtclust = c(drtclust_ms1, drtclust_ms2),
                                       minpeak = c(minpeak_ms1, minpeak_ms2),
                                       drtgap = c(drtgap_ms1, drtgap_ms2),
                                       drtminpeak = c(drtminpeak_ms1, drtminpeak_ms2),
                                       drtmaxpeak = c(drtmaxpeak_ms1, drtmaxpeak_ms2),
                                       recurs = c(recurs_ms1, recurs_ms2),
                                       sb = c(sb_ms1, sb_ms2),
                                       sn = c(sn_ms1, sn_ms2),
                                       minint = c(minint_ms1, minint_ms2),
                                       weight = c(weight_ms1, weight_ms2),
                                       dmzIso = c(dmzIso_ms1, dmzIso_ms2),
                                       drtIso = c(drtIso_ms1, drtIso_ms2))
    }
    
    # If polarity is positive
    if (polarity == "positive"){
      for (m in 1:length(msobjects)){
        msobjects[[m]] <- idPOS(msobjects[[m]],
                                ppm_precursor = dmzprecursor,
                                ppm_products = dmzproducts,
                                rttol = rttol,
                                coelCutoff = coelcutoff,
                                lipidClasses = lipidClassesPos)
      }
    }
    
    # If polarity is negative
    if (polarity == "negative"){
      for (m in 1:length(msobjects)){
        msobjects[[m]] <- idNEG(msobjects[[m]],
                                ppm_precursor = dmzprecursor,
                                ppm_products = dmzproducts,
                                rttol = rttol,
                                coelCutoff = coelcutoff,
                                lipidClasses = lipidClassesNeg)
      }
    }
    
    for (m in 1:length(msobjects)){
      msobjects[[m]] <- plotLipids(msobjects[[m]])
    }
    
    return(msobjects)
  }
  
  
  
  # Function to dynamically change the selection of lipid classes based on the selected polarity.
  select_lipidclasses <- function(analysis_step, dbs) {
    polarity <- polarityReactive()
    
    if (analysis_step == "lipidclasses_predict") {
      label <- "Lipid classes to predict"
    } else if (analysis_step == "lipidclasses_annotate") {
      label <- "Lipid classes to annotate"
    }
    
    # Custom labels based on list type.
    label_custom <- switch(
      analysis_step,
      "lipidclasses_annotate" = "Lipid classes to annotate (csv list)",
      "lipidclasses_predict" = "Lipid classes to predict (csv list)"
    )
    
    label_nopolarity <- "No data imported"
    label_nofile <- "No file selected"
    
    # label if polarity is NULL
    if (is.null(polarity)) {
      label_html <- HTML(paste0(
        "<b>", label_custom, "</b><br>",
        "<span style='color:gray;'>&nbsp;&nbsp;&nbsp;&nbsp;", label_nopolarity, "</span>"
      ))
      
      updateCheckboxGroupInput(session, analysis_step,
                               label = label_html,
                               choices = list(),
                               selected = NULL)
      return()
    }
    
    # --- we only get here if polarity is not NULL ---

    # names of lipid classes and adducts
    if(!is.null(dbs)) {
      choicesNames <- list()
      for (i in 1:length(dbs)) {
        lipid <- names(dbs)[i]
        choicesNames[[i]] <- HTML(paste0(
          "<b style='color:black;'>", lipid, "  &nbsp;&nbsp;&nbsp;&nbsp; </b>
      <span style='color:gray;'><small>(",
          paste(colnames(dbs[[lipid]])[-c(1:3)], collapse = ";  &nbsp;"),
          ")</small></span>"
        ))
      }
      
      lipid_adducts <- choicesNames
      names(lipid_adducts) <- names(dbs)
      choicesNamesReactive(lipid_adducts)
      
      choicesValues <- as.list(names(dbs)) 
    }
    
    # label based on polarity
    if (polarity == "positive") {
      label <- paste(label, "for ESI + :")
    } else if (polarity == "negative") {
      label <- paste(label, "for ESI - :")
    }
    
    # Display checkbox group based on class type (default or custom)
    if (input$type_lipidclasses == "lipidclasses_default") {
      updateCheckboxGroupInput(session, analysis_step,
                               label = label,
                               choiceNames = choicesNames,
                               choiceValues = choicesValues,
                               selected = choicesValues)
      
    } else if (input$type_lipidclasses == "lipidclasses_custom") {
      if (!is.null(dbs)) {
        updateCheckboxGroupInput(session, analysis_step,
                                 label = label_custom,
                                 choiceNames = choicesNames,
                                 choiceValues = choicesValues,
                                 selected = choicesValues)
      } else {
        label_html <- HTML(paste0(
          "<b>", label_custom, "</b><br>",
          "<span style='color:gray;'>&nbsp;&nbsp;&nbsp;&nbsp;", label_nofile, "</span>"
        ))
        
        updateCheckboxGroupInput(session, analysis_step,
                                 label = label_html,
                                 choices = list(),
                                 selected = NULL)
      }
    }
  }
  
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)


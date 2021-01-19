options(shiny.maxRequestSize=500*1024^2)

shinyServer(function(input, output) {
  observeEvent(input$do,
               {output$txt_result <- renderText({"Job completed"})})
  observeEvent(input$do, {
    annotated_msobjects <- startWork(input$file1$datapath, input$file1$name,
                                     input$sI_acquisitionmode, input$sI_polarity,
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
                                     input$lipidClassesPos,
                                     input$lipidClassesNeg)

    output$summaryTable <- renderUI({
      lapply(1:length(annotated_msobjects), function(i) {
        output[[annotated_msobjects[[i]]$metaData$generalMetadata$file]] <-
          renderTable({
          annotated_msobjects[[i]]$results
        })
      })
    })

    output$annotatedPeaklist <- renderUI({
      lapply(1:length(annotated_msobjects), function(i) {
        output[[annotated_msobjects[[i]]$metaData$generalMetadata$file]] <- renderTable({
          annotated_msobjects[[i]]$annotatedPeaklist
        })
      })
    })

    output$downloadSummary <- downloadHandler(
      filename = function(){paste(input$jobname, "SummaryTables.zip", sep="_")},
      content = function(file){
        #go to a temp dir to avoid permission issues
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        files <- NULL;

        #loop through the sheets
        for (i in 1:length(annotated_msobjects)){
          #write each sheet to a csv file, save the name
          fileName <- paste(gsub(".mzXML", "" , input$file1$name[i]), "_summaryResults.csv", sep="")
          write.csv(annotated_msobjects[[i]]$results, fileName, row.names = FALSE)
          files <- c(fileName, files)
        }
        #create the zip file
        zip(file, files)
      }
    )

    output$downloadPeaklist <- downloadHandler(
      filename = function(){paste(input$jobname, "AnnotatedPeaklists.zip", sep="_")},
      content = function(file){
        #go to a temp dir to avoid permission issues
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        files <- NULL;

        #loop through the sheets
        for (i in 1:length(annotated_msobjects)){
          #write each sheet to a csv file, save the name
          fileName <- paste(gsub(".mzXML", "" , input$file1$name[i]), "_annotatedPeaklist.csv", sep="")
          write.csv(annotated_msobjects[[i]]$annotatedPeaklist, fileName, row.names = FALSE)
          files <- c(fileName, files)
        }
        #create the zip file
        zip(file, files)
      }
    )

    output$downloadPlots <- downloadHandler(
      filename = function(){paste(input$jobname, "Plots.zip", sep="_")},
      content = function(file){
        #go to a temp dir to avoid permission issues
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        files <- NULL;

        # loop through the msobjects
        for (i in 1:length(annotated_msobjects)){
          fileName <- paste(gsub(".mzXML", "" , input$file1$name[i]), "_plots.pdf", sep="")
          if (annotated_msobjects[[i]]$metaData$acquisitionmode == "DIA"){
            height <- 7
          } else {
            height <- 8
          }
          grDevices::pdf(file = fileName, width = 8, height = height)
          for (pl in 1:length(annotated_msobjects[[i]]$plots)){
            print(annotated_msobjects[[i]]$plots[[pl]])
          }
          grDevices::dev.off()
          files <- c(fileName, files)
        }
        #create the zip file
        zip(file, files)
      }
    )
  })
})

startWork <- function(files, filesname, acquisitionmode, polarity,
                     dmzagglom_ms1, dmzagglom_ms2, drtagglom_ms1, drtagglom_ms2,
                     drtclust_ms1, drtclust_ms2, minpeak_ms1, minpeak_ms2,
                     drtgap_ms1, drtgap_ms2, drtminpeak_ms1, drtminpeak_ms2,
                     drtmaxpeak_ms1, drtmaxpeak_ms2, recurs_ms1, recurs_ms2,
                     sb_ms1, sb_ms2, sn_ms1, sn_ms2, minint_ms1, minint_ms2,
                     weight_ms1, weight_ms2, dmzIso_ms1, dmzIso_ms2, drtIso_ms1,
                     drtIso_ms2, dmzprecursor, dmzproducts, rttol, coelcutoff,
                     jobname, lipidClassesPos, lipidClassesNeg, plotresults){

  ################################################################################
  # dataProcessing
  msobjects <- list()

  for (f in 1:length(files)){
    msobjects[[f]] <- dataProcessing(file = files[f],
                                     polarity = polarity,
                                     acquisitionmode = acquisitionmode,
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

  ################################################################################
  # annotation
  annotated_msobjects <- list()

  # If polarity is positive
  if (polarity == "positive"){
    for (m in 1:length(msobjects)){
      annotated_msobjects[[m]] <- idPOS(msobjects[[m]],
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
      annotated_msobjects[[m]] <- idNEG(msobjects[[m]],
                                        ppm_precursor = dmzprecursor,
                                        ppm_products = dmzproducts,
                                        rttol = rttol,
                                        coelCutoff = coelcutoff,
                                        lipidClasses = lipidClassesNeg)
    }
  }

  for (m in 1:length(annotated_msobjects)){
    annotated_msobjects[[m]] <- plotLipids(annotated_msobjects[[m]])
  }

  return(annotated_msobjects)
}

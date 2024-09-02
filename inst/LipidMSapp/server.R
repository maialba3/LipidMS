options(shiny.maxRequestSize=500*1024^2)

shinyServer(function(input, output, session) {
  
  observeEvent(input$JumpTo2, {
    updateTabsetPanel(session, "inTabset",
                      selected = "tab2")
  })
  
  observeEvent(input$GoBackTo1, {
    updateTabsetPanel(session, "inTabset",
                      selected = "tab1")
  })
  
  observeEvent(input$JumpTo3, {
    updateTabsetPanel(session, "inTabset",
                      selected = "tab3")
  })
  
  observeEvent(input$GoBackTo2, {
    updateTabsetPanel(session, "inTabset",
                      selected = "tab2")
  })
  
  observeEvent(input$JumpTo4, {
    updateTabsetPanel(session, "inTabset",
                      selected = "tab4")
  })
  
  observeEvent(input$GoBackTo3, {
    updateTabsetPanel(session, "inTabset",
                      selected = "tab3")
  })
  
  observeEvent(input$JumpTo5, {
    updateTabsetPanel(session, "inTabset",
                      selected = "tab5")
  })
  
  observeEvent(input$GoBackTo4, {
    updateTabsetPanel(session, "inTabset",
                      selected = "tab4")
  })
  
  observe({
    req(input$file1, input$metadata)
    metadata2 <- read.csv(input$metadata$datapath, sep=",")
    if (!all(grepl(".mzXML", metadata2$sample))){
      metadata2$sample <- paste(metadata2$sample, ".mzXML", sep = "")
    }
    files <- data.frame(sample = input$file1$name,
                        path = input$file1$datapath)
    metadata2 <- merge(metadata2, files, by = "sample")
    output$metadata <- renderTable({metadata2[,1:3]})
  })
  
  
  # observeEvent(input$do,
  #              {output$txt_result <- renderText({"Job completed"})})
  observeEvent(input$do, {
    req(input$file1, input$metadata)
    if (input$analysis == "single"){
      metadata <- read.csv(input$metadata$datapath, sep=",")
      if (!all(grepl(".mzXML", metadata$sample))){
        metadata$sample <- paste(metadata$sample, ".mzXML", sep = "")
      }
      files <- data.frame(sample = input$file1$name,
                          path = input$file1$datapath)
      metadata <- merge(metadata, files, by = "sample")
      msobjects <- singleProcessing(metadata$path, metadata$sample,
                                              metadata$acquisitionmode, input$sI_polarity,
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
                                              input$lipidClassesPos,
                                              input$lipidClassesNeg)
      
      output$summaryTable <- renderUI({
        lapply(1:length(msobjects), function(i) {
          output[[msobjects[[i]]$metaData$generalMetadata$file]] <-
            renderTable({
              msobjects[[i]]$annotation$results
            })
        })
      })
      
      output$annotatedPeaklist <- renderUI({
        lapply(1:length(msobjects), function(i) {
          output[[msobjects[[i]]$metaData$generalMetadata$file]] <- renderTable({
            msobjects[[i]]$annotation$annotatedPeaklist
          })
        })
      })
      
      output$features <- renderTable({data.frame("No results for batch processing")})
      
      output$downloadSummary <- downloadHandler(
        filename = function(){paste(input$jobname, "SummaryTables.zip", sep="_")},
        content = function(file){
          #go to a temp dir to avoid permission issues
          owd <- setwd(tempdir())
          on.exit(setwd(owd))
          files <- NULL;
          
          #loop through the sheets
          for (i in 1:length(msobjects)){
            #write each sheet to a csv file, save the name
            fileName <- paste(gsub(".mzXML", "" , input$file1$name[i]), "_summaryResults.csv", sep="")
            write.csv(msobjects[[i]]$annotation$results, fileName, row.names = FALSE)
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
          for (i in 1:length(msobjects)){
            #write each sheet to a csv file, save the name
            fileName <- paste(gsub(".mzXML", "" , input$file1$name[i]), "_annotatedPeaklist.csv", sep="")
            write.csv(msobjects[[i]]$annotation$annotatedPeaklist, fileName, row.names = FALSE)
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
          for (i in 1:length(msobjects)){
            fileName <- paste(gsub(".mzXML", "" , input$file1$name[i]), "_plots.pdf", sep="")
            if (msobjects[[i]]$metaData$generalMetadata$acquisitionmode == "DIA"){
              height <- 7
            } else {
              height <- 8
            }
            grDevices::pdf(file = fileName, width = 8, height = height)
            if (length(msobjects[[i]]$annotation$plots) > 0){
              for (pl in 1:length(msobjects[[i]]$annotation$plots)){
                print(msobjects[[i]]$annotation$plots[[pl]])
              }
            }
            grDevices::dev.off()
            files <- c(fileName, files)
          }
          #create the zip file
          zip(file, files)
        }
      )
      
    } else if (input$analysis == "batch"){
      metadata <- read.csv(input$metadata$datapath, sep=",")
      if (!all(grepl(".mzXML", metadata$sample))){
        metadata$sample <- paste(metadata$sample, ".mzXML", sep = "")
      }
      files <- data.frame(sample = input$file1$name,
                          path = input$file1$datapath)
      metadata <- merge(metadata, files, by = "sample")
      msbatch <- batchProcessing(metadata = metadata,
                                 input$sI_polarity,
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
                                 input$parallel, input$ncores,
                                 input$dmzprecursor, input$dmzproducts,
                                 input$rttol, input$coelcutoff,
                                 input$jobname,
                                 input$lipidClassesPos,
                                 input$lipidClassesNeg)
      
      output$summaryTable <- renderUI({
        lapply(1:length(msbatch$msobjects), function(i) {
          output[[msbatch$msobjects[[i]]$metaData$generalMetadata$file]] <-
            renderTable({
              msbatch$msobjects[[i]]$annotation$results
            })
        })
      })
      
      output$annotatedPeaklist <- renderUI({
        lapply(1:length(msbatch$msobjects), function(i) {
          output[[msbatch$msobjects[[i]]$metaData$generalMetadata$file]] <- renderTable({
            msbatch$msobjects[[i]]$annotation$annotatedPeaklist
          })
        })
      })
      
      output$features <- renderTable({msbatch$features})
      
      output$downloadSummary <- downloadHandler(
        filename = function(){paste(input$jobname, "SummaryTables.zip", sep="_")},
        content = function(file){
          # go to a temp dir to avoid permission issues
          owd <- setwd(tempdir())
          on.exit(setwd(owd))
          files <- NULL;
          
          # loop through the sheets
          for (i in 1:length(msbatch$msobjects)){
            if (msbatch$msobjects[[i]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
              # write each sheet to a csv file, save the name
              fileName <- gsub(".mzXML", "_summaryResults.csv" , msbatch$metaData$sample[i])
              write.csv(msbatch$msobjects[[i]]$annotation$results, fileName, row.names = FALSE)
              files <- c(fileName, files)
            }
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
          for (i in 1:length(msbatch$msobjects)){
            if (msbatch$msobjects[[i]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
              #write each sheet to a csv file, save the name
              fileName <- gsub(".mzXML", "_annotatedPeaklist.csv" , msbatch$metaData$sample[i])
              write.csv(msbatch$msobjects[[i]]$results$annotatedPeaklist, fileName, row.names = FALSE)
              files <- c(fileName, files)
            }
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
          for (i in 1:length(msbatch$msobjects)){
            if (msbatch$msobjects[[i]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
              fileName <- gsub(".mzXML", "_plots.pdf" , msbatch$metaData$sample[i])
              if (msbatch$msobjects[[i]]$metaData$generalMetadata$acquisitionmode == "DIA"){
                height <- 7
              } else {
                height <- 8
              }
              grDevices::pdf(file = fileName, width = 8, height = height)
              if (length(msbatch$msobjects[[i]]$annotation$plots) > 0){
                for (pl in 1:length(msbatch$msobjects[[i]]$annotation$plots)){
                  print(msbatch$msobjects[[i]]$annotation$plots[[pl]])
                }
              }
              grDevices::dev.off()
              files <- c(fileName, files)
            }
          }
          #create the zip file
          zip(file, files)
        }
      )
      
      output$downloadBatchResults <- downloadHandler(
        filename = function(){paste(input$jobname, "BatchResults.zip", sep="_")},
        content = function(file){
          #go to a temp dir to avoid permission issues
          owd <- setwd(tempdir())
          on.exit(setwd(owd));
          
          # extract data
          peaklist <- msbatch$features
          peaklistNoIso <- peaklist[peaklist$isotope %in% c("", "[M+0]"),]
          
          # write files
          write.csv(peaklist, "FeaturesMatrix.csv", row.names = FALSE)
          write.csv(peaklistNoIso, "FeaturesMatrixIsotopesRemoved.csv", row.names = FALSE)
          pdf("RTdevplot.pdf")
          rtdevplot(msbatch)
          rtdevplot(msbatch, colorbygroup = FALSE)
          dev.off()
          pdf("TIC.pdf", height = 7, width = 10)
          plotticmsbatch(msbatch)
          plotticmsbatch(msbatch, colorbygroup = FALSE)
          dev.off()
          
          
          files <- c("FeaturesMatrix.csv", "FeaturesMatrixIsotopesRemoved.csv",
                     "RTdevplot.pdf", "TIC.pdf")
          
          #create the zip file
          zip(file, files)
        }
      )
    }
  })
})

batchProcessing <- function(metadata, polarity,
                     dmzagglom_ms1, dmzagglom_ms2, drtagglom_ms1, drtagglom_ms2,
                     drtclust_ms1, drtclust_ms2, minpeak_ms1, minpeak_ms2,
                     drtgap_ms1, drtgap_ms2, drtminpeak_ms1, drtminpeak_ms2,
                     drtmaxpeak_ms1, drtmaxpeak_ms2, recurs_ms1, recurs_ms2,
                     sb_ms1, sb_ms2, sn_ms1, sn_ms2, minint_ms1, minint_ms2,
                     weight_ms1, weight_ms2, dmzIso_ms1, dmzIso_ms2, drtIso_ms1,
                     drtIso_ms2, dmzalign, drtalign, span, minsamplesfracalign, 
                     dmzgroup, drtagglomgroup, drtgroup, minsamplesfracgroup,
                     parallel, ncores, dmzprecursor, dmzproducts, rttol, coelcutoff,
                     jobname, lipidClassesPos, lipidClassesNeg){

  #==============================================================================#
  # Procesamiento
  #==============================================================================#
  
  ###################
  # Peak-picking
  samplenames <- metadata$sample
  metadata$sample <- metadata$path
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
                                 ncores = ncores)
  
  msbatch$metaData$sample <- samplenames
  
  ###################
  # Alignment
  msbatch <- alignmsbatch(msbatch, dmz = dmzalign, drt = drtalign, span = span, 
                          minsamplesfrac = minsamplesfracalign, 
                          parallel = parallel, ncores = ncores)
  
  ###################
  # Grouping
  msbatch <- groupmsbatch(msbatch, dmz = dmzgroup, drtagglom = drtagglomgroup,
                          drt = drtgroup, minsamplesfrac = minsamplesfracgroup, 
                          parallel = parallel, ncores = ncores)
  
  
  #####################
  # Fill missing peaks
  msbatch <- fillpeaksmsbatch(msbatch)
  
  
  ###################
  # Lipid Annotation
  msbatch <- annotatemsbatch(msbatch,
                             ppm_precursor = dmzprecursor,
                             ppm_products = dmzproducts,
                             rttol = rttol,
                             coelCutoff = coelcutoff,
                             lipidClassesPos = lipidClassesPos,
                             lipidClassesNeg = lipidClassesNeg)
  
  for (m in 1:length(msbatch$msobjects)){
    if (msbatch$msobjects[[m]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
      msbatch$msobjects[[m]] <- plotLipids(msbatch$msobjects[[m]])
    }
  }

  return(msbatch)
}


singleProcessing <- function(files, filesname, acquisitionmode, polarity,
                             dmzagglom_ms1, dmzagglom_ms2, drtagglom_ms1, drtagglom_ms2,
                             drtclust_ms1, drtclust_ms2, minpeak_ms1, minpeak_ms2,
                             drtgap_ms1, drtgap_ms2, drtminpeak_ms1, drtminpeak_ms2,
                             drtmaxpeak_ms1, drtmaxpeak_ms2, recurs_ms1, recurs_ms2,
                             sb_ms1, sb_ms2, sn_ms1, sn_ms2, minint_ms1, minint_ms2,
                             weight_ms1, weight_ms2, dmzIso_ms1, dmzIso_ms2, drtIso_ms1,
                             drtIso_ms2, dmzprecursor, dmzproducts, rttol, coelcutoff,
                             jobname, lipidClassesPos, lipidClassesNeg){
  
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
  msobjects <- list()
  
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

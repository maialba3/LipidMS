# plotLipids
#' Plot informative peaks for lipids annotation
#'
#' Plot informative peaks for each lipid annotated using idPOS and idNEG (or
#' similar functions).
#'
#' @param msobject annotated msobject.
#' @param spar smoothing parameter. Numeric value between 0 and 1.
#'
#' @return msobject with a plots element which contains a list of plots.
#' Plots on the left side represent raw values while plots on the left side are
#' smoothed or clean scans (MS2 in DDA).
#'
#' @details Peak intensities are relative to the maximum intensity of each peak
#' to ease visualization.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#' 
#' library(LipidMS)
#' library(LipidMSdata2)
#'
#' msobject <- idPOS(LipidMSdata2::msobjectDIApos)
#' msobject <- plotLipids(msobject)
#' 
#' # display the first plot
#' msobject$plots[[1]]
#' msobject$plots[["yourpeakIDofinterest"]]
#' 
#' # save all plot to a pdf file
#' pdf("plotresults.pdf")
#' for (p in 1:length(msobject$plots)){
#'   print(msobject$plots[[p]])
#' }
#' dev.off()
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
plotLipids <- function(msobject, spar = 0.4){

  ##############################################################################
  # Check arguments and inputs
  if (!"results" %in% names(msobject)){
    stop("No results to be plotted")
  }
  if (spar < 0 | spar > 1){
    spar <- 0.4
    warning("spar parameter set to 0.4")
  }
  if ("plots" %in% names(msobject)){
    cat("\n Removing previous plots...")
    msobject$plots <- NULL
  }
  results <- msobject$results
  msobject$plots <- list()
  ##############################################################################
  # For each lipid in the results data frame, extract peaks
  while (nrow(results) > 0){
    r <- results$peakID[1]
    rclass <- results$Class[1]
    ############################################################################
    # MS1 peaks
    toremove <- which(results$peakID == r & results$Class == rclass)
    parent <- results[results$peakID == r & results$Class == rclass,]
    if (nrow(parent) > 1){
      parent$ID[1] <- paste(parent$ID, collapse = "|")
      parent$peakID[1] <- paste(unique(parent$peakID), collapse = "|")
      parent <- parent[1,]
    }
    adducts <- unlist(strsplit(parent$Adducts, ";"))
    ms1 <- c()
    for (a in adducts){
      ss <- msobject$detailsAnnotation[[parent$Class]]$candidates
      ss <- ss[grepl(gsub("+", "\\+", as.character(paste("^", a, sep="")), fixed = TRUE),
                     msobject$detailsAnnotation[[parent$Class]]$candidates$adducts),]
      ss <- ss[ss$cb == parent$CDB,]
      ss <- ss[order(abs(ss$RT - parent$RT), decreasing = FALSE),]
      ms1 <- rbind(ms1, ss[1,])
    }
    ms1$adducts <- as.character(sapply(ms1$adducts, function(x) unlist(strsplit(x, ";"))[[1]]))
    peaksMS1 <- ms1$peakID
    mzpeaksMS1 <- ms1$m.z # to use in case data is DDA
    namesMS1 <- paste(as.character(round(ms1$m.z, 3)), ms1$adducts, sep="_")
    
    # eics <- list()
    # for (m in 1:length(mzpeaksMS1)){
    #   eic
    # }
    
    ############################################################################
    # MS2 peaks
    peaksMS2 <- c()
    scansMS2 <- c() # to use in case data is DDA
    namesMS2 <- c()
    class <- c()
    chains <- c()

    # get index for all adducts
    c <- which(msobject$detailsAnnotation[[parent$Class]]$candidates$peakID %in% ms1$peakID)

    # extract class fragments
    if ("classfragments" %in% names(msobject$detailsAnnotation[[parent$Class]])){
      class <- do.call(rbind, msobject$detailsAnnotation[[parent$Class]]$classfragments[c])
      if (length(class) > 0){
        if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
          peaksMS2 <- c(peaksMS2, class$peakID)
        } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
          peaksMS2 <- c(peaksMS2, class$m.z)
          scansMS2 <- c(scansMS2, class$peakID)
        }
        namesMS2 <- c(namesMS2, paste(as.character(round(class$m.z, 3)), "class fragment", sep="_"))
      }
    }

    # extract chain fragments
    if ("chainfragments" %in% names(msobject$detailsAnnotation[[parent$Class]])){
      if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
        for (i in c){
          ch <- c()
          if (length(msobject$detailsAnnotation[[parent$Class]]$chainfragments[[i]]) > 0){
            ch <- do.call(rbind, msobject$detailsAnnotation[[parent$Class]]$chainfragments[[i]])
            if (length(ch) > 0){
              chains <- rbind(chains, ch)
              ch <- c()
            }
          }
        }
        if (length(chains) > 0){
          peaksMS2 <- c(peaksMS2, chains$peakID)
          namesMS2 <- c(namesMS2, paste(as.character(round(chains$m.z, 3)), chains$db, chains$cb, chains$adduct, sep="_"))
          namesMS2 <- namesMS2[peaksMS2 != ""]
          peaksMS2 <- peaksMS2[peaksMS2 != ""]
          namesMS2 <- unique(namesMS2)
          peaksMS2 <- unique(peaksMS2)
        }
      } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
        for (i in c){
          ch <- c()
          if (length(msobject$detailsAnnotation[[parent$Class]]$chainfragments[[i]]) > 0){
            ch <- do.call(rbind, msobject$detailsAnnotation[[parent$Class]]$chainfragments[[i]])
            if (length(ch) > 0){
              chains <- rbind(chains, ch)
              ch <- c()
            }
          }
        }
        if (length(chains) > 0){
          chains <- unique(chains)
          chains <- chains[chains$m.z != 0,]
          peaksMS2 <- c(peaksMS2, chains$m.z)
          scansMS2 <- c(scansMS2, chains$peakID)
          namesMS2 <- c(namesMS2, paste(as.character(round(chains$m.z, 3)), chains$db, chains$cb, chains$adduct, sep="_"))
        }
      }
    }

    # id data is DDA, extract spectrum data
    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      rawMS <- c()
      for (i in c){
        f <- c()
        if (length(msobject$detailsAnnotation[[parent$Class]]$coelfrags[[i]]) > 0){
          f <- msobject$detailsAnnotation[[parent$Class]]$coelfrags[[i]]
          if (length(f) > 0){
            rawMS <- rbind(rawMS, f)
            f <- c()
          }
        }
      }
    }

    ############################################################################
    # Plot
    colorsMS1 <- c("#42858C", "#FE9300", "#870E75", "#3E71A8", "#FE6900")
    colorsMS2 <- c("#7F8E39", "#5F3659", "#E5C616", "#16A08CFF", "#628395",
                   "#C5D86D", "#969696FF", "#358359FF", "#9F4D23FF", "#D86C4FFF",
                   "#170C2EFF", "#473B75FF", "#F19C1FFF")

    if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
      nplots <- 1 + length(unique(rawMS$peakID))
    } else {
      nplots <- 2
    }
    grDevices::pdf(NULL) # use a pdf NULL device to save plots to an object
    grDevices::dev.control(displaylist = "enable")
    graphics::par(mfrow=c(nplots, 2), mar = c(3,4,4,1), mgp=c(2,1,0), bg = "white")
    # plot MS1 info
    if (length(peaksMS1) > 0){
      ssms1 <- msobject$MS1[msobject$MS1$peakID %in% peaksMS1,]
      minrt1 <- min(ssms1$RT)
      maxrt1 <- max(ssms1$RT)
      ints <- c()
      for (p in 1:length(peaksMS1)){
        toplot <- msobject$MS1[msobject$MS1$peakID == peaksMS1[p],]
        toplot <- toplot[order(toplot$RT, decreasing = FALSE),]
        ints <- append(ints, max(toplot$int))
        toplot$int <- toplot$int*100/max(toplot$int)
        if (p == 1){
          plot(toplot$RT, toplot$int, type = "l", col = scales::alpha(colorsMS1[p], 0.8),
               xlim = c(minrt1-5, maxrt1+5), ylim = c(0, 110),
               lwd = 2.5, ylab = "Rel. Intensity", xlab = "RT (sec)",
               main = paste("MS1: ", paste(parent$ID, as.character(round(parent$m.z, 2)),
                                           as.character(round(parent$RT, 1)), sep="_"), sep = ""),
               las = 1, cex.axis = 0.7, cex.lab = 1, cex.main = 1)
        } else {
          graphics::lines(toplot$RT, toplot$int, col = scales::alpha(colorsMS1[p], 0.8), lwd = 2.5)
        }
      }
      graphics::legend("topright", legend=namesMS1,
             col=colorsMS1[1:length(peaksMS1)], lty=1, lwd = 2, cex=0.6)
      graphics::legend("bottomright", title = "Max. intensity",
             legend=formatC(ints, format = "e", digits = 2),
             col=colorsMS1[1:length(peaksMS1)], lty=1, lwd = 2, cex=0.6)

      # smoothed
      for (p in 1:length(peaksMS1)){
        toplot <- msobject$MS1[msobject$MS1$peakID == peaksMS1[p],]
        toplot <- toplot[order(toplot$RT, decreasing = FALSE),]
        pred <- tryCatch({stats::predict(stats::smooth.spline(toplot$RT, toplot$int, spar = spar),
                                  x = toplot$RT)},
                         error = function(e) {return(list(x = toplot$RT,
                                                          y = toplot$int))})
        toplot$RT <- pred$x
        toplot$int <- pred$y
        toplot$int <- toplot$int*100/max(toplot$int)
        if (p == 1){
          plot(toplot$RT, toplot$int, type = "l", col = scales::alpha(colorsMS1[p], 0.8),
               xlim = c(minrt1-5, maxrt1+5), ylim = c(0, 110),
               lwd = 2.5, ylab = "Rel. Intensity", xlab = "RT (sec)",
               main = paste("MS1: ", paste(parent$ID, as.character(round(parent$m.z, 2)),
                                           as.character(round(parent$RT, 1)), sep="_"), sep = ""),
               las = 1, cex.axis = 0.7, cex.lab = 1, cex.main = 1, lty = 5)
        } else {
          graphics::lines(toplot$RT, toplot$int, col = scales::alpha(colorsMS1[p], 0.8), lwd = 2.5,
                lty = 5)
        }
      }
      graphics::legend("topright", legend=namesMS1,
             col=colorsMS1[1:length(peaksMS1)], lty = 5, lwd = 2, cex=0.6)
      graphics::legend("bottomright", title = "Max. intensity",
             legend=formatC(ints, format = "e", digits = 2),
             col=colorsMS1[1:length(peaksMS1)], lty = 5, lwd = 2, cex=0.6)
    }

    # plot MS2 info
    if (length(peaksMS2) > 0){
      # if data is DIA
      if (msobject$metaData$generalMetadata$acquisitionmode == "DIA"){
        ssms2 <- msobject$MS2[msobject$MS2$peakID %in% peaksMS2,]
        minrt2 <- min(ssms2$RT)
        maxrt2 <- max(ssms2$RT)
        maxint2 <- max(ssms2$int)
        ints2 <- c()
        for (p in 1:length(peaksMS2)){
          toplot <- msobject$MS2[msobject$MS2$peakID == peaksMS2[p],]
          toplot <- toplot[order(toplot$RT, decreasing = FALSE),]
          ints2 <- append(ints2, max(toplot$int))
          toplot$int <- toplot$int/max(toplot$int)
          toplot$int <- toplot$int*100/max(toplot$int)
          if (p == 1){
            plot(toplot$RT, toplot$int, type = "l", col = scales::alpha(colorsMS2[p], 0.8),
                 xlim = c(minrt1-5, maxrt1+5), ylim = c(0, 110),
                 lwd = 2.5, ylab = "Rel. Intensity", xlab = "RT (sec)",
                 main = paste("MS2: ", paste(parent$ID, as.character(round(parent$m.z, 2)),
                                             as.character(round(parent$RT, 1)), sep="_"), sep = ""),
                 las = 1, cex.axis = 0.7, cex.lab = 1, cex.main = 1)
          } else {
            graphics::lines(toplot$RT, toplot$int, col = scales::alpha(colorsMS2[p], 0.8), lwd = 2.5)
          }
        }
        graphics::legend("topright", legend=namesMS2,
               col=colorsMS2[1:length(peaksMS2)], lty=1, lwd = 2, cex=0.6)
        graphics::legend("bottomright", title = "Max. intensity",
               legend=formatC(ints2, format = "e", digits = 2),
               col=colorsMS2[1:length(peaksMS2)], lty=1, lwd = 2, cex=0.6)

        # smoothed
        for (p in 1:length(peaksMS2)){
          toplot <- msobject$MS2[msobject$MS2$peakID == peaksMS2[p],]
          toplot <- toplot[order(toplot$RT, decreasing = FALSE),]
          toplot$int <- toplot$int/max(toplot$int)
          pred <- tryCatch({stats::predict(stats::smooth.spline(toplot$RT, toplot$int, spar = spar),
                                    x = toplot$RT)},
                           error = function(e) {return(list(x = toplot$RT,
                                                            y = toplot$int))})
          toplot$RT <- pred$x
          toplot$int <- pred$y
          toplot$int <- toplot$int*100/max(toplot$int)
          if (p == 1){
            plot(toplot$RT, toplot$int, type = "l", col = scales::alpha(colorsMS2[p], 0.8),
                 xlim = c(minrt1-5, maxrt1+5), ylim = c(0, 110),
                 lwd = 2.5, ylab = "Rel. Intensity", xlab = "RT (sec)",
                 main = paste("MS2: ", paste(parent$ID, as.character(round(parent$m.z, 2)),
                                             as.character(round(parent$RT, 1)), sep="_"), sep = ""),
                 las = 1, cex.axis = 0.7, cex.lab = 1, cex.main = 1, lty = 5)
          } else {
            graphics::lines(toplot$RT, toplot$int, col = scales::alpha(colorsMS2[p], 0.8), lwd = 2.5,
                  lty = 5)
          }
        }
        graphics::legend("topright", legend=namesMS2,
               col=colorsMS2[1:length(peaksMS2)], lty = 5, lwd = 2, cex=0.6)
        graphics::legend("bottomright", title = "Max. intensity",
               legend=formatC(ints2, format = "e", digits = 2),
               col=colorsMS2[1:length(peaksMS2)], lty = 5, lwd = 2, cex=0.6)
        # if data is DDA
      } else if (msobject$metaData$generalMetadata$acquisitionmode == "DDA"){
        # for each scan
        for (s in unique(rawMS$peakID)){
          # subset raw data
          ssrawMS <- rawMS[rawMS$peakID == s,]
          ssmaxrawMS <- max(ssrawMS$int)
          ssrawMS$int <- ssrawMS$int*100/max(ssrawMS$int)
          ssrawMS$int[ssrawMS$int < 2] <- ssrawMS$int[ssrawMS$int < 2] + 2
          mz2 <- peaksMS2[scansMS2 == s]
          namesmz2 <- namesMS2[scansMS2 == s]
          namesmz2 <- namesmz2[order(mz2, decreasing = FALSE)]
          mz2 <- mz2[order(mz2, decreasing = FALSE)]
          
          # assign colors
          ssrawMS$color <- "black"
          ssrawMS$color[ssrawMS$m.z %in% mz2] <- colorsMS2[1:sum(ssrawMS$m.z %in% mz2)]

          # Find precursor in the MS/MS spectrum
          scanprec <- unlist(strsplit(s, "_"))
          collisionenergy <- as.numeric(scanprec[2])
          scanprec <- as.numeric(scanprec[3])
          precursor <- msobject$metaData$scansMetadata$precursor[
            which(msobject$metaData$scansMetadata$msLevel == 2 &
                    msobject$metaData$scansMetadata$collisionEnergy == collisionenergy)[scanprec]]
          prec <- as.numeric(unlist(sapply(precursor, mzMatch, ssrawMS$m.z, ppm = 10)))
          if (length(prec) > 0){
            prec <- prec[seq(1, length(prec), 2)]
            mzprec <- ssrawMS$m.z[mz1p]
            nameprec <- paste(round(mzprec, 3), "_precursor", sep="")
            if (ssrawMS$color[prec] == "black"){
              ssrawMS$color[prec] <- colorsMS2[length(namesmz2)+1]
              namesmz2 <- c(namesmz2, nameprec)
            }
          }

          #plot
          plot(ssrawMS$m.z, ssrawMS$int, type = "h", col = scales::alpha(ssrawMS$color, 1),
               xlim = c(0, max(ssrawMS$m.z)+20), ylim = c(0, 132),
               lwd = 1, ylab = "Rel. Intensity", xlab = "m/z",
               main = paste("MS2: ", paste(parent$ID, as.character(round(parent$m.z, 2)),
                                           as.character(round(parent$RT, 1)), sep="_"),
                            paste("\nPrecursor: ", round(precursor, 3), sep = ""),
                            sep = ""),
               las = 1, cex.axis = 0.7, cex.lab = 1, cex.main = 1, lty = 1, yaxt = "n" )
          graphics::axis(2,at=seq(2, 122, 20), labels = seq(0, 120, 20))
          graphics::legend("topright", legend=namesmz2,
                 col=colorsMS2[1:length(namesmz2)], lty = 1, lwd = 2, cex=0.6)

          # clean
          ssrawMSclean <- ssrawMS[ssrawMS$color != "black",]
          plot(ssrawMSclean$m.z, ssrawMSclean$int, type = "h",
               col = scales::alpha(ssrawMSclean$color, 1),
               xlim = c(0, max(ssrawMS$m.z)+20), ylim = c(0, 132),
               lwd = 1.5, ylab = "Rel. Intensity", xlab = "m/z",
               main = paste("MS2: ", paste(parent$ID, as.character(round(parent$m.z, 2)),
                                           as.character(round(parent$RT, 1)), sep="_"),
                            paste("\nPrecursor: ", round(precursor, 3), sep = ""),
                            sep = ""),
               las = 1, cex.axis = 0.7, cex.lab = 1, cex.main = 1, lty = 1, yaxt = "n" )
          graphics::axis(2,at=seq(2, 102, 20), labels = seq(0, 100, 20))
          graphics::legend("topright", legend=namesmz2,
                 col=colorsMS2[1:length(namesmz2)], lty = 1, lwd = 2, cex=0.6)
        }
      }
    }
    msobject$plots[[r]] <- grDevices::recordPlot() # save plot
    invisible(grDevices::dev.off()) # close pdf NULL device
    results <- results[-toremove,]
  }
  return(msobject)
}

# createLipidDB
#' Customizable lipid DBs creator
#'
#' It allows to create easy-customizable lipid DBs for annotation with LipidMS
#' package.
#'
#' @param lipid character value indicating the class of lipid. See Details.
#' @param chains character vector indicating the FA chains to be employed
#' @param chains2 character vector containing the sphingoid bases to be employed
#' if required.
#'
#' @return List with the requested dbs (data frames)
#'
#' @details \code{lipidClass} argument needs to be one of the following
#' character values: "Cer", "CerP", "GlcCer", "SM", "Carnitine", "CE", "FA",
#' "HFA", "Sph" (sphingoid bases), "SphP", "MG", "LPA", , "LPC",
#' "LPE", "LPG", "LPI", "LPS", "FAHFA", "DG", "PC", "PE", "PG", "PI", "PS",
#' "PA", "TG", "CL" or "all".
#'
#' @examples
#' fas <- c("8:0", "10:0", "12:0", "14:0", "14:1", "15:0", "16:0", "16:1",
#' "17:0", "18:0", "18:1", "18:2", "18:3", "18:4", "20:0", "20:1", "20:2",
#' "20:3", "20:4", "20:5", "22:0", "22:1", "22:2", "22:3", "22:4", "22:5",
#' "22:6", "24:0", "24:1", "26:0")
#' sph <- c("16:0", "16:1", "18:0", "18:1")
#' newdb <- createLipidDB(lipid = "PC", chains = fas, chains2 = sph)
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
createLipidDB <- function(lipid, chains, chains2){
  customizedDataSets <- list()
  if (sum(lipid == "Cer" || lipid == "CerP" || lipid == "GlcCer" ||
          lipid == "SM") == 1){
    db <- dbSphingolipids(chains = chains, chains2 = chains2, lipid = lipid)
    db <- data.frame(formula=db$formula, total=db$total,
                     Mass=as.numeric(db$Mass), ID = paste(lipid, "(", db$total, ")", sep=""),
                     stringsAsFactors = F)
    if (lipid == "Cer"){
      customizedDataSets[["cerdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                  Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "CerP"){
      customizedDataSets[["cerPdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "GlcCer"){
      customizedDataSets[["glccerdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                     Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "SM"){
      customizedDataSets[["smdb"]] <- data.frame(formula=db$formula,
                                                 total=db$total, Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
  } else if (sum(lipid == "FA" || lipid == "HFA" || lipid == "Carnitine" ||
                 lipid == "LPA" || lipid == "LPE" || lipid == "LPG" ||
                 lipid == "LPI" || lipid == "LPS" || lipid == "LPC" ||
                 lipid == "MG" || lipid == "CE" || lipid == "Sph" ||
                 lipid == "SphP" || lipid == "LPEo" || lipid == "LPAo" ||
                 lipid == "LPCp" || lipid == "LPCo" || lipid == "LPEp") == 1){
    if (lipid %in% c("Sph", "SphP")){
      db <- dbOneChain(chains = chains2, lipid = lipid)
      db <- data.frame(formula=db$formula, total=db$total,
                       Mass=as.numeric(db$Mass), ID = paste(lipid, "(", db$total, ")", sep=""),
                       stringsAsFactors = F)
    } else {
      db <- dbOneChain(chains = chains, lipid = lipid)
      db <- data.frame(formula=db$formula, total=db$total,
                       Mass=as.numeric(db$Mass), ID = paste(lipid, "(", db$total, ")", sep=""),
                       stringsAsFactors = F)
    }
    if (lipid == "FA"){
      customizedDataSets[["fadb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "HFA"){
      customizedDataSets[["hfadb"]] <- data.frame(formula=db$formula, total=db$total,
                                                  Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "Carnitine"){
      customizedDataSets[["carnitinedb"]] <- data.frame(formula=db$formula, total=db$total,
                                                        Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPA"){
      customizedDataSets[["lysopadb"]] <- data.frame(formula=db$formula, total=db$total,
                                                     Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPE"){
      customizedDataSets[["lysopedb"]] <- data.frame(formula=db$formula, total=db$total,
                                                     Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPG"){
      customizedDataSets[["lysopgdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                     Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPI"){
      customizedDataSets[["lysopidb"]] <- data.frame(formula=db$formula, total=db$total,
                                                     Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPS"){
      customizedDataSets[["lysopsdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                     Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPC"){
      customizedDataSets[["lysopcdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                     Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "MG"){
      customizedDataSets[["mgdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "CE"){
      customizedDataSets[["CEdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "Sph"){
      customizedDataSets[["sphdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                  Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "SphP"){
      customizedDataSets[["sphPdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPEo"){
      customizedDataSets[["lysopeodb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPAo"){
      customizedDataSets[["lysopaodb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPCp"){
      customizedDataSets[["lysopcpdb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPCo"){
      customizedDataSets[["lysopcodb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPEp"){
      customizedDataSets[["lysopepdb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
  } else if (sum(lipid == "DG" || lipid == "PC" || lipid == "PE" ||
                 lipid == "PG" || lipid == "PI" || lipid == "PS" || lipid == "PIP" ||
                 lipid == "PIP2" ||lipid == "PIP3" || lipid == "FAHFA" ||
                 lipid == "PA") == 1 || lipid == "PEo" || lipid == "PCp" ||
             lipid == "PCo" || lipid == "PEp"){
    db <- dbTwoChains(chains = chains, lipid = lipid)
    db <- data.frame(formula=db$formula, total=db$total,
                     Mass=as.numeric(db$Mass), ID = paste(lipid, "(", db$total, ")", sep=""),
                     stringsAsFactors = F)
    if (lipid == "FAHFA"){
      customizedDataSets[["fahfadb"]] <- data.frame(formula=db$formula, total=db$total,
                                                    Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "DG"){
      customizedDataSets[["dgdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PE"){
      customizedDataSets[["pedb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PG"){
      customizedDataSets[["pgdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PI"){
      customizedDataSets[["pidb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PIP"){
      customizedDataSets[["pipdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                  Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PIP2"){
      customizedDataSets[["pip2db"]] <- data.frame(formula=db$formula, total=db$total,
                                                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PIP3"){
      customizedDataSets[["pip3db"]] <- data.frame(formula=db$formula, total=db$total,
                                                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PS"){
      customizedDataSets[["psdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PC"){
      customizedDataSets[["pcdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PA"){
      customizedDataSets[["padb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PEo"){
      customizedDataSets[["peodb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PCp"){
      customizedDataSets[["pcpdb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PCo"){
      customizedDataSets[["pcodb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PEp"){
      customizedDataSets[["pepdb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
  } else if (lipid == "TG") {
    db <- dbThreeChains(chains = chains, lipid = lipid)
    db <- data.frame(formula=db$formula, total=db$total,
                     Mass=as.numeric(db$Mass), ID = paste(lipid, "(", db$total, ")", sep=""),
                     stringsAsFactors = F)
    customizedDataSets[["tgdb"]] <- data.frame(formula=db$formula, total=db$total,
                                               Mass=as.numeric(db$Mass), stringsAsFactors = F)
  } else if (lipid == "CL") {
    db <- dbFourChains(chains = chains, lipid = lipid)
    db <- data.frame(formula=db$formula, total=db$total,
                     Mass=as.numeric(db$Mass), ID = paste(lipid, "(", db$total, ")", sep=""),
                     stringsAsFactors = F)
    customizedDataSets[["cldb"]] <- data.frame(formula=db$formula, total=db$total,
                                               Mass=as.numeric(db$Mass), stringsAsFactors = F)
  } else if (lipid == "all"){
    ceramides <- dbSphingolipids(chains = chains, chains2 = chains2,
                                 lipid = "Cer")
    ceramides <- data.frame(formula=ceramides$formula, total=ceramides$total,
                            Mass=as.numeric(ceramides$Mass), ID = paste("Cer(", ceramides$total,
                                                                        ")", sep=""), stringsAsFactors = F)
    ceramidesP <- dbSphingolipids(chains = chains, chains2 = chains2,
                                  lipid = "CerP")
    ceramidesP <- data.frame(formula=ceramidesP$formula, total=ceramidesP$total,
                             Mass=as.numeric(ceramidesP$Mass), ID = paste("CerP(", ceramidesP$total,
                                                                          ")", sep=""), stringsAsFactors = F)
    glccer <- dbSphingolipids(chains = chains, chains2 = chains2,
                              lipid = "GlcCer")
    glccer <- data.frame(formula=glccer$formula, total=glccer$total,
                         Mass=as.numeric(glccer$Mass), ID = paste("GlcCer(", glccer$total,
                                                                  ")", sep=""), stringsAsFactors = F)
    sm <- dbSphingolipids(chains = chains, chains2 = chains2, lipid = "SM")
    sm <- data.frame(formula=sm$formula, total=sm$total,
                     Mass=as.numeric(sm$Mass), ID = paste("SM(", sm$total,
                                                          ")", sep=""), stringsAsFactors = F)
    fa <- dbOneChain(chains = chains, lipid = "FA")
    fa <- data.frame(formula=fa$formula, total=fa$total,
                     Mass=as.numeric(fa$Mass), ID = paste("FA(", fa$total,
                                                          ")", sep=""), stringsAsFactors = F)
    hfa <- dbOneChain(chains = chains, lipid = "HFA")
    hfa <- data.frame(formula=hfa$formula, total=hfa$total,
                      Mass=as.numeric(hfa$Mass), ID = paste("HFA(", hfa$total,
                                                            ")", sep=""), stringsAsFactors = F)
    carnitine <- dbOneChain(chains = chains, lipid = "Carnitine")
    carnitine <- data.frame(formula=carnitine$formula, total=carnitine$total,
                            Mass=as.numeric(carnitine$Mass), ID = paste("Carnitine(", carnitine$total,
                                                                        ")", sep=""), stringsAsFactors = F)
    CE <- dbOneChain(chains = chains, lipid = "CE")
    CE <- data.frame(formula=CE$formula, total=CE$total,
                     Mass=as.numeric(CE$Mass), ID = paste("CE(", CE$total,
                                                          ")", sep=""), stringsAsFactors = F)
    mg <- dbOneChain(chains = chains, lipid = "MG")
    mg <- data.frame(formula=mg$formula, total=mg$total,
                     Mass=as.numeric(mg$Mass), ID = paste("MG(", mg$total,
                                                          ")", sep=""), stringsAsFactors = F)
    sph <- dbOneChain(chains = chains2, lipid = "Sph")
    sph <- data.frame(formula=sph$formula, total=sph$total,
                      Mass=as.numeric(sph$Mass), ID = paste("Sph(", sph$total,
                                                            ")", sep=""), stringsAsFactors = F)
    sphP <- dbOneChain(chains = chains2, lipid = "SphP")
    sphP <- data.frame(formula=sphP$formula, total=sphP$total,
                       Mass=as.numeric(sphP$Mass), ID = paste("SphP(", sphP$total,
                                                              ")", sep=""), stringsAsFactors = F)
    lysopc <- dbOneChain(chains = chains, lipid = "LPC")
    lysopc <- data.frame(formula=lysopc$formula, total=lysopc$total,
                         Mass=as.numeric(lysopc$Mass), ID = paste("LPC(", lysopc$total,
                                                                  ")", sep=""), stringsAsFactors = F)
    lysope <- dbOneChain(chains = chains, lipid = "LPE")
    lysope <- data.frame(formula=lysope$formula, total=lysope$total,
                         Mass=as.numeric(lysope$Mass), ID = paste("LPE(", lysope$total,
                                                                  ")", sep=""), stringsAsFactors = F)
    lysopg <- dbOneChain(chains = chains, lipid = "LPG")
    lysopg <- data.frame(formula=lysopg$formula, total=lysopg$total,
                         Mass=as.numeric(lysopg$Mass), ID = paste("LPG(", lysopg$total,
                                                                  ")", sep=""), stringsAsFactors = F)
    lysopi <- dbOneChain(chains = chains, lipid = "LPI")
    lysopi <- data.frame(formula=lysopi$formula, total=lysopi$total,
                         Mass=as.numeric(lysopi$Mass), ID = paste("LPI(", lysopi$total,
                                                                  ")", sep=""), stringsAsFactors = F)
    lysops <- dbOneChain(chains = chains, lipid = "LPS")
    lysops <- data.frame(formula=lysops$formula, total=lysops$total,
                         Mass=as.numeric(lysops$Mass), ID = paste("LPS(", lysops$total,
                                                                  ")", sep=""), stringsAsFactors = F)
    lysopa <- dbOneChain(chains = chains, lipid = "LPA")
    lysopa <- data.frame(formula=lysopa$formula, total=lysopa$total,
                         Mass=as.numeric(lysopa$Mass), ID = paste("LPA(", lysopa$total,
                                                                  ")", sep=""), stringsAsFactors = F)
    pc <- dbTwoChains(chains = chains, lipid = "PC")
    pc <- data.frame(formula=pc$formula, total=pc$total,
                     Mass=as.numeric(pc$Mass), ID = paste("PC(", pc$total,
                                                          ")", sep=""), stringsAsFactors = F)
    pe <- dbTwoChains(chains = chains, lipid = "PE")
    pe <- data.frame(formula=pe$formula, total=pe$total,
                     Mass=as.numeric(pe$Mass), ID = paste("PE(", pe$total,
                                                          ")", sep=""), stringsAsFactors = F)
    pg <- dbTwoChains(chains = chains, lipid = "PG")
    pg <- data.frame(formula=pg$formula, total=pg$total,
                     Mass=as.numeric(pg$Mass), ID = paste("PG(", pg$total,
                                                          ")", sep=""), stringsAsFactors = F)
    pi <- dbTwoChains(chains = chains, lipid = "PI")
    pi <- data.frame(formula=pi$formula, total=pi$total,
                     Mass=as.numeric(pi$Mass), ID = paste("PI(", pi$total,
                                                          ")", sep=""), stringsAsFactors = F)
    pip <- dbTwoChains(chains = chains, lipid = "PIP")
    pip <- data.frame(formula=pip$formula, total=pip$total,
                      Mass=as.numeric(pip$Mass), ID = paste("PIP(", pip$total,
                                                            ")", sep=""), stringsAsFactors = F)
    pip2 <- dbTwoChains(chains = chains, lipid = "PIP2")
    pip2 <- data.frame(formula=pip2$formula, total=pip2$total,
                       Mass=as.numeric(pip2$Mass), ID = paste("PIP2(", pip2$total,
                                                              ")", sep=""), stringsAsFactors = F)
    pip3 <- dbTwoChains(chains = chains, lipid = "PIP3")
    pip3 <- data.frame(formula=pip3$formula, total=pip3$total,
                       Mass=as.numeric(pip3$Mass), ID = paste("PIP3(", pip3$total,
                                                              ")", sep=""), stringsAsFactors = F)
    ps <- dbTwoChains(chains = chains, lipid = "PS")
    ps <- data.frame(formula=ps$formula, total=ps$total,
                     Mass=as.numeric(ps$Mass), ID = paste("PS(", ps$total,
                                                          ")", sep=""), stringsAsFactors = F)
    pa <- dbTwoChains(chains = chains, lipid = "PA")
    pa <- data.frame(formula=pa$formula, total=pa$total,
                     Mass=as.numeric(pa$Mass), ID = paste("PA(", pa$total,
                                                          ")", sep=""), stringsAsFactors = F)
    fahfa <- dbTwoChains(chains = chains, lipid = "FAHFA")
    fahfa <- data.frame(formula=fahfa$formula, total=fahfa$total,
                        Mass=as.numeric(fahfa$Mass), ID = paste("FAHFA(", fahfa$total,
                                                                ")", sep=""), stringsAsFactors = F)
    dg <- dbTwoChains(chains = chains, lipid = "DG")
    dg <- data.frame(formula=dg$formula, total=dg$total,
                     Mass=as.numeric(dg$Mass), ID = paste("DG(", dg$total,
                                                          ")", sep=""), stringsAsFactors = F)
    tg <- dbThreeChains(chains = chains, lipid = "TG")
    tg <- data.frame(formula=tg$formula, total=tg$total,
                     Mass=as.numeric(tg$Mass), ID = paste("TG(", tg$total,
                                                          ")", sep=""), stringsAsFactors = F)
    cl <- dbFourChains(chains = chains, lipid = "CL")
    cl <- data.frame(formula=cl$formula, total=cl$total,
                     Mass=as.numeric(cl$Mass), ID = paste("CL(", cl$total,
                                                          ")", sep=""), stringsAsFactors = F)
    peo <- dbTwoChains(chains = chains, lipid = "PEo")
    peo <- data.frame(formula=peo$formula, total=peo$total,
                      Mass=as.numeric(peo$Mass),
                      ID = paste("PEo(", peo$total, ")", sep=""),
                      stringsAsFactors = F)
    lysopeo <- dbOneChain(chains = chains, lipid = "LPEo")
    lysopeo <- data.frame(formula=lysopeo$formula, total=lysopeo$total,
                          Mass=as.numeric(lysopeo$Mass),
                          ID = paste("LPEo(", pe$total, ")", sep=""),
                          stringsAsFactors = F)
    lysopao <- dbOneChain(chains = chains, lipid = "LPAo")
    lysopao <- data.frame(formula=lysopao$formula, total=lysopao$total,
                          Mass=as.numeric(lysopao$Mass),
                          ID = paste("LPAo(", pe$total, ")", sep=""),
                          stringsAsFactors = F)
    pcp <- dbTwoChains(chains = chains, lipid = "PCp")
    pcp <- data.frame(formula=pcp$formula, total=pcp$total,
                      Mass=as.numeric(pcp$Mass),
                      ID = paste("PCp(", pcp$total, ")", sep=""),
                      stringsAsFactors = F)
    lysopcp <- dbOneChain(chains = chains, lipid = "LPCp")
    lysopcp <- data.frame(formula=lysopcp$formula, total=lysopcp$total,
                          Mass=as.numeric(lysopcp$Mass),
                          ID = paste("LPCp(", lysopcp$total, ")", sep=""),
                          stringsAsFactors = F)
    pco <- dbTwoChains(chains = chains, lipid = "PCo")
    pco <- data.frame(formula=pco$formula, total=pco$total,
                      Mass=as.numeric(pco$Mass),
                      ID = paste("PCo(", pco$total, ")", sep=""),
                      stringsAsFactors = F)
    lysopco <- dbOneChain(chains = chains, lipid = "LPCo")
    lysopco <- data.frame(formula=lysopco$formula, total=lysopco$total,
                          Mass=as.numeric(lysopco$Mass),
                          ID = paste("LPCo(", lysopco$total, ")", sep=""),
                          stringsAsFactors = F)
    pep <- dbTwoChains(chains = chains, lipid = "PEp")
    pep <- data.frame(formula=pep$formula, total=pep$total,
                      Mass=as.numeric(pep$Mass),
                      ID = paste("PEp(", pco$total, ")", sep=""),
                      stringsAsFactors = F)
    lysopep <- dbOneChain(chains = chains, lipid = "LPEp")
    lysopep <- data.frame(formula=lysopep$formula, total=lysopep$total,
                          Mass=as.numeric(lysopep$Mass),
                          ID = paste("LPEp(", lysopep$total, ")", sep=""),
                          stringsAsFactors = F)
    customizedDataSets[["cerdb"]] <- ceramides
    customizedDataSets[["cerPdb"]] <- ceramidesP
    customizedDataSets[["glccerdb"]] <- glccer
    customizedDataSets[["smdb"]] <- sm
    customizedDataSets[["fadb"]] <- fa
    customizedDataSets[["hfadb"]] <- hfa
    customizedDataSets[["carnitinedb"]] <- carnitine
    customizedDataSets[["lysopadb"]] <- lysopa
    customizedDataSets[["lysopedb"]] <- lysope
    customizedDataSets[["lysopgdb"]] <- lysopg
    customizedDataSets[["lysopidb"]] <- lysopi
    customizedDataSets[["lysopsdb"]] <- lysops
    customizedDataSets[["lysopcdb"]] <- lysopc
    customizedDataSets[["mgdb"]] <- mg
    customizedDataSets[["CEdb"]] <- CE
    customizedDataSets[["sphdb"]] <- sph
    customizedDataSets[["sphPdb"]] <- sphP
    customizedDataSets[["fahfadb"]] <- fahfa
    customizedDataSets[["pedb"]] <- pe
    customizedDataSets[["pgdb"]] <- pg
    customizedDataSets[["pidb"]] <- pi
    customizedDataSets[["pipdb"]] <- pip
    customizedDataSets[["pip2db"]] <- pip2
    customizedDataSets[["pip3db"]] <- pip3
    customizedDataSets[["psdb"]] <- ps
    customizedDataSets[["padb"]] <- pa
    customizedDataSets[["pcdb"]] <- pc
    customizedDataSets[["dgdb"]] <- dg
    customizedDataSets[["tgdb"]] <- tg
    customizedDataSets[["cldb"]] <- cl
    customizedDataSets[["peodb"]] <- peo
    customizedDataSets[["lysopeodb"]] <- lysopeo
    customizedDataSets[["lysopaodb"]] <- lysopao
    customizedDataSets[["pcpdb"]] <- pcp
    customizedDataSets[["lysopcpdb"]] <- lysopcp
    customizedDataSets[["pcodb"]] <- pco
    customizedDataSets[["lysopcodb"]] <- lysopco
    customizedDataSets[["pepdb"]] <- pep
    customizedDataSets[["lysopepdb"]] <- lysopep
  }
  return(customizedDataSets)
}

# assignDB
#' Load LipidMS default data bases
#'
#' load all LipidMS default data bases required to run identification functions.
#'
#' @return list of data frames
#'
#' @examples
#' \dontrun{
#' dbs <- assignDB()
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
assignDB <- function(){
  dbs <- list()
  dbs[["cerdb"]] <- LipidMS::cerdb
  dbs[["smdb"]] <- LipidMS::smdb
  dbs[["fadb"]] <- LipidMS::fadb
  dbs[["hfadb"]] <- LipidMS::hfadb
  dbs[["carnitinedb"]] <- LipidMS::carnitinesdb
  dbs[["lysopadb"]] <- LipidMS::lysopadb
  dbs[["lysopedb"]] <- LipidMS::lysopedb
  dbs[["lysopgdb"]] <- LipidMS::lysopgdb
  dbs[["lysopidb"]] <- LipidMS::lysopidb
  dbs[["lysopsdb"]] <- LipidMS::lysopsdb
  dbs[["lysopcdb"]] <- LipidMS::lysopcdb
  dbs[["mgdb"]] <- LipidMS::mgdb
  dbs[["CEdb"]] <- LipidMS::CEdb
  dbs[["sphdb"]] <- LipidMS::sphdb
  dbs[["sphPdb"]] <- LipidMS::sphPdb
  dbs[["fahfadb"]] <- LipidMS::fahfadb
  dbs[["pcdb"]] <- LipidMS::pcdb
  dbs[["pedb"]] <- LipidMS::pedb
  dbs[["pgdb"]] <- LipidMS::pgdb
  dbs[["pidb"]] <- LipidMS::pidb
  dbs[["psdb"]] <- LipidMS::psdb
  dbs[["padb"]] <- LipidMS::padb
  dbs[["dgdb"]] <- LipidMS::dgdb
  dbs[["tgdb"]] <- LipidMS::tgdb
  dbs[["cldb"]] <- LipidMS::cldb
  dbs[["badb"]] <- LipidMS::badb
  dbs[["baconjdb"]] <- LipidMS::baconjdb
  dbs[["nlsphdb"]] <- LipidMS::nlsphdb
  dbs[["adductsTable"]] <- LipidMS::adductsTable
  return(dbs)
}

# getInclusionList
#' Obtain an inclusion list from the annotation results
#'
#' Obtain an inclusion list from the annotation results.
#'
#' @param results data frame. Output of identification functions.
#' @param adductsTable data frame with the adducts allowed and their mass
#' difference.
#'
#' @return Data frame with 6 columns: formula, RT, neutral mass, m/z, adduct
#' and the compound name.
#'
#' @examples
#' \dontrun{
#' devtools::install_github("maialba3/LipidMSdata2")
#'
#' library(LipidMS)
#' msobject <- idPOS(LipidMSdata2::msobjectDIApos)
#' getInclusionList(msobject)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
getInclusionList <- function(results, adductsTable = LipidMS::adductsTable){
  Form_Mn <- apply(results, 1, getFormula)
  if (class(Form_Mn) == "matrix"){
    new <- list()
    for (i in 1:ncol(Form_Mn)){
      new[[i]] <- c(Form_Mn[1,i], Form_Mn[2,i])
    }
    Form_Mn <- new
  }
  na <- which(unlist(lapply(Form_Mn, length)) == 0)
  if (length(na) > 0){
    Form_Mn[[na]] <- c(NA, NA)
  }
  Formula <- unlist(lapply(Form_Mn, "[[", 1))
  RT <- results$RT
  Mn <- as.numeric(unlist(lapply(Form_Mn, "[[", 2)))
  adducts <- sapply(as.vector(results$Adduct), strsplit, ";")
  mzs <- rep(list(vector()), nrow(results))
  for (i in 1:nrow(results)){
    ad <- adducts[[i]]
    for (a in 1:length(ad)){
      adinfo <- adductsTable[adductsTable$adduct == ad[a],]
      mz <- (adinfo$n*Mn[i]+adinfo$mdif)/abs(adinfo$charge)
      mzs[[i]] <- append(mzs[[i]], mz)
    }
  }
  Name <- results$ID
  inclusionList <- vector()
  for (i in 1:nrow(results)){
    for (a in 1:length(adducts[[i]])){
      inclusionList <- rbind(inclusionList,
                             data.frame(Formula[i], RT[i], Mn[i],
                                        mzs[[i]][a], adducts[[i]][a],
                                        Name[i], stringsAsFactors = F))
    }
  }
  colnames(inclusionList) <- c("Formula", "RT", "Mn", "m.z", "Adduct", "Name")
  inclusionList <- unique(inclusionList)
  return(inclusionList)
}

# searchIsotopes
#' Target isotopes search
#'
#' This function uses annotation results of an unlabelled sample to search
#' for labelled compounds in a labelled sample.
#'
#' @param msobject a msobject generated by any of the identification functions.
#' @param label isotope employed for the experiment. It can be "13C" or "D".
#' @param adductsTable adducts table employed for lipids annotation.
#' @param rttol rt window in seconds.
#' @param ppm mass error tolerance.
#' @param coelCutoff coelution score threshold between isotopes. By default, 0.8.
#'
#' @return List with the isotopes for each compound in the results data frame.
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
searchIsotopes <- function(msobject,
                           label,
                           adductsTable = LipidMS::adductsTable,
                           rttol = 10,
                           ppm = 10,
                           coelCutoff = 0.8){

  results <- msobject$results
  MS1 <- msobject$peaklist$MS1

  Form_Mn <- apply(results, 1, getFormula)
  formula <- unlist(lapply(Form_Mn, "[[", "Formula"))
  results$Formula <- formula
  comp <- do.call(rbind, sapply(results$Formula, function(x) {
    c <- CHNOSZ::makeup(x)
    c <- c[c("C", "H", "N", "O", "P")]
    names(c) <- c("C", "H", "N", "O", "P")
    return(as.data.frame(c))
  }))
  colnames(comp) <- c("C", "H", "N", "O", "P")
  results$C <- comp[,"C"]
  results$H <- comp[,"H"]
  if (label == "13C"){
    massdiff <- 1.0033548
    label <- "C"
  } else if (label == "D"){
    massdiff <- 1.006277
    label <- "H"
  }
  isotopes <- apply(results, 1, function(x){
    mz <- as.numeric(x["m.z"])
    rt <- as.numeric(x["RT"])
    top <- as.numeric(x[label])
    intensities <- vector()
    mzs <- vector()
    rts <- vector()
    for (i in 0:top){
      if (nrow(MS1[abs(MS1$RT - rt) < rttol,]) > 0){
        subsetMS1 <- MS1[abs(MS1$RT - rt) < rttol,]
        m <- mzMatch(mz+i*massdiff, subsetMS1[abs(subsetMS1$RT - rt) <
                                                rttol, "m.z"], ppm)
        if (length(m) > 0){
          m <- m[seq(1, length(m), 2)]
          sel <- subsetMS1[m,]
          sel$coelScore <- coelutionScore(as.character(x["peakID"]), sel$peakID, msobject$MS1)
          sel <- sel[sel$coelScore >= coelCutoff,]
          if (nrow(sel) > 0){
            int <- sel[order(abs(sel$RT - rt)),"int"][1]
            RT <- sel[order(abs(sel$RT - rt)),"RT"][1]
            MZ <- sel[order(abs(sel$RT - rt)),"m.z"][1]
          } else {
            int <- 0
            RT <- rt
            MZ <- mz+i*massdiff
          }
        } else {
          int <- 0
          RT <- rt
          MZ <- mz+i*massdiff
        }
      } else {
        int <- 0
        RT <- rt
        MZ <- mz+i*massdiff
      }
      intensities <- append(intensities, int)
      mzs <- append(mzs, MZ)
      rts <- append(rts, RT)
    }
    adduct <- unlist(strsplit(as.character(x["Adducts"]), ";"))[1]
    return(data.frame(Name = paste(x["ID"], " [M+", 0:top, "]",
                                   sep = ""),
                      Adduct = rep(adduct, length(mzs)),
                      Formula = rep(x["Formula"], length(mzs)),
                      m.z = mzs, RT = rts,  int = intensities,
                      stringsAsFactors = F))
  })
  return(isotopes)
}

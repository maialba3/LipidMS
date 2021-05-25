
################################################################################
# LipidMS 3.0.0:
################################################################################

library(LipidMS)

#==============================================================================#
# Input
#==============================================================================#
# Example data files can be downloaded from 
# https://drive.google.com/drive/folders/1hSYrQBkh-rAA-oiaKqGkrNL7uWQraV75?usp=sharing

setwd("yourWorkingDirectory")

# csv file with 3 columns: sample (mzXML file names), acquisitionmode 
# (MS, DIA or DDA) and sampletype (QC, group1, group2, etc.)
metadata <- read.csv("Metadata.csv", sep=",", dec = ".") 

#==============================================================================#
# Processing parameters
#==============================================================================#

###################
# Peak-picking
polarity <- "positive" # 6550 abrir hasta 50 ppm, con el orbi dejar en 15-20 ppm
dmzagglom <- 15
drtagglom <- 500
drtclust <- 100
minpeak <- c(5, 3)
drtgap <- 5
drtminpeak <- 15
drtmaxpeak <- c(100, 200)
recurs <- 5
sb <- 2
sn <- 2
minint <- c(1000, 500)
weight <- c(2, 3)
dmzIso <- 5
drtIso <- 5

###################
# Batch processing
dmzalign <- 5
drtalign <- 30
span <- 0.4
minsamplesfracalign <- 0.75
dmzgroup <- 5
drtagglomgroup <- 30
drtgroup <- 15
minsamplesfracgroup <- 0.25
parallel <- TRUE
ncores <- 2


#==============================================================================#
# Processing
#==============================================================================#

###################
# Peak-picking
msbatch <- batchdataProcessing(metadata = metadata,
                               polarity = polarity,
                               dmzagglom = dmzagglom,
                               drtagglom = drtagglom,
                               drtclust = drtclust,
                               minpeak = minpeak,
                               drtgap = drtgap,
                               drtminpeak = drtminpeak,
                               drtmaxpeak = drtmaxpeak,
                               recurs = recurs,
                               sb = sb,
                               sn = sn,
                               minint = minint,
                               weight = weight,
                               dmzIso = dmzIso,
                               drtIso = drtIso,
                               parallel = parallel,
                               ncores = ncores)

save(msbatch, file="msbatch.rda.gz", compress = TRUE)

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
msbatch <- annotatemsbatch(msbatch)


#==============================================================================#
# Output
#==============================================================================#

###################
# features
peaklist <- msbatch$features
peaklistNoIso <- peaklist[peaklist$isotope %in% c("", "[M+0]"),]

View(peaklistNoIso)

write.csv(peaklist, file="peaklist.csv")
write.csv(peaklistNoIso, file="peaklistNoIso.csv")

###################
# annotations

# results
for (i in 1:length(msbatch$msobjects)){
  if (msbatch$msobjects[[i]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    fileName <- gsub(".mzXML", "_summaryResults.csv" , 
                     msbatch$msobjects[[i]]$metaData$generalMetadata$file[i])
    write.csv(msbatch$msobjects[[i]]$annotation$results, fileName, row.names = FALSE)
  }
}

# Annotated Peaklists
for (i in 1:length(msbatch$msobjects)){
  if (msbatch$msobjects[[i]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    fileName <- gsub(".mzXML", "_annotatedPeaklist.csv" , 
                     msbatch$msobjects[[i]]$metaData$generalMetadata$file[i])
    write.csv(msbatch$msobjects[[i]]$annotation$annotatedPeaklist, fileName, row.names = FALSE)
  }
}

###################
# plots

pdf("RTdevplot.pdf")
rtdevplot(msbatch)
rtdevplot(msbatch, colorbygroup = FALSE)
dev.off()

pdf("TIC.pdf", height = 7, width = 10)
plotticmsbatch(msbatch)
plotticmsbatch(msbatch, colorbygroup = FALSE)
dev.off()

# identified lipids plots
for (m in 1:length(msbatch$msobjects)){
  if (msbatch$msobjects[[m]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    msbatch$msobjects[[m]] <- plotLipids(msbatch$msobjects[[m]])
  }
}
# print them to pdf
for (s in 1:length(msbatch$msobjects)){
  if (msbatch$msobjects[[s]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
    print(s)
    if (msbatch$msobjects[[s]]$metaData$generalMetadata$acquisitionmode == "DIA"){
      height <- 7
    } else {
      height <- 9
    }
    pdf(file = gsub(".mzXML", "_plots.pdf", msbatch$msobjects[[s]]$metaData$generalMetadata$file), 
        width = 8, height = height)
    for ( pl in 1:length(msbatch$msobjects[[s]]$annotation$plots)){
      print(msbatch$msobjects[[s]]$annotation$plots[[pl]])
    }
    dev.off()
  }
}



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
files <- dir()[grepl("mzXML", dir())]

#==============================================================================#
# Processing parameters
#==============================================================================#

acquisitionmode <- rep("DIA", length(files)) # it assumes you add MS, DDA or DIA to the file names
acquisitionmode[grep("DDA", files)] <- "DDA" 
acquisitionmode[grep("MS", files)] <- "MS" 
# otherwise, create a vector such as: acquisitionmode <- c("DIA", "MS", "MS", "MS", "DDA",...)
polarity <- "positive" # change to positive if required
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

# Annotation parameters
dmzprecursor <- 5
dmzproducts <- 10
rttol <- 6
coelcutoff <- 0.7


#==============================================================================#
# Processing
#==============================================================================#

# for single file processing, only samples acquired in DIA or DDA will be processed
files <- files[acquisitionmode %in% c("DIA", "DDA")]
acquisitionmode <- acquisitionmode[acquisitionmode %in% c("DIA", "DDA")]

#################
# Peak-picking
msobjects <- list()

for (f in 1:length(files)){
  msobjects[[f]] <- dataProcessing(file = files[f],
                                   acquisitionmode = acquisitionmode[f],
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
                                   drtIso = drtIso)
}
save(msobjects, file = "msobjects.rda")

#################
# Annotation
msobjects <- list() # we could also overwrite the previous msobjects list

# If polarity is positive
if (polarity == "positive"){
  for (m in 1:length(msobjects)){
    msobjects[[m]] <- idPOS(msobjects[[m]],
                                      ppm_precursor = dmzprecursor,
                                      ppm_products = dmzproducts,
                                      rttol = rttol,
                                      coelCutoff = coelcutoff,
                                      lipidClasses = c("MG", "LPC", "LPE", "PC", 
                                                       "PE", "Sph", "SphP", "Cer", 
                                                       "SM", "Carnitines", "CE", 
                                                       "DG", "TG"))
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
                                      lipidClasses = c("FA", "FAHFA", "LPC", "LPE", 
                                                       "LPG", "LPI", "LPS", "PC", 
                                                       "PE", "PG", "PI", "PS", "Sph", 
                                                       "SphP", "Cer", "CL", "BA"))
  }
}
save(msobjects, file = "msobjects.rda")

#==============================================================================#
# Output
#==============================================================================#

load("msobjects.rda")

View(msobjects[[1]]$results)

# write tables to csv files
for (f in 1:length(msobjects)){
  fname <- gsub(".mzXML", "", msobjects[[f]]$metaData$generalMetadata$file)
  write.csv(msobjects[[f]]$annotation$results, file = paste0(fname, "_summaryResults.csv"),
            row.names = FALSE)
  write.csv(msobjects[[f]]$annotation$annotatedPeaklist, 
            file = paste0(fname, "_annotatedPeaklist.csv"),
            row.names = FALSE)
}

# plots
for (f in 1:length(msobjects)){
  msobjects[[f]] <- plotLipids(msobjects[[f]])
  
  fname <- gsub(".mzXML", "", msobjects[[f]]$metaData$generalMetadata$file)
  pdf(paste0(fname, ".pdf"))
  for (p in 1:length(msobjects[[f]]$annotation$plots)){
    print(msobjects[[f]]$annotation$plots[[p]])
  }
  dev.off()
}

save(msobjects, file = "msobjects.rda")

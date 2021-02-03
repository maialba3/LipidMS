
################################################################################
# LipidMS 2.0.0:
################################################################################

library(LipidMS)

################################################################################
# Input
setwd("yourWorkingDirectory")
files <- dir()[grepl("mzXML", dir())]

# Processing parameters
acquisitionmode <- rep("DIA", length(files)) # it assumes you add DDA or DIA to the file names
acquisitionmode[grep("DDA", files)] <- "DDA" 
# otherwise, create a vector such as: acquisitionmode <- c("DIA", "DIA", "DIA", "DDA", "DDA",...)
polarity <- "negative" # change to positive if required
dmzagglom <- 5
drtagglom <- 25
drtclust <- 25
minpeak <- c(4, 3)
drtgap <- 5
drtminpeak <- 20
drtmaxpeak <- 200
recurs <- 5
sb <- 2
sn <- 2
minint <- c(500, 100)
weight <- c(2, 3)
dmzIso <- 5
drtIso <- 5

# Annotation parameters
dmzprecursor <- 5
dmzproducts <- 10
rttol <- 6
coelcutoff <- 0.7


################################################################################
# Processing
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

################################################################################
# Annotation
annotated_msobjects <- list() # we could also overwrite the previous msobjects list

# If polarity is positive
if (polarity == "positive"){
  for (m in 1:length(msobjects)){
    annotated_msobjects[[m]] <- idPOS(msobjects[[m]],
                                      ppm_precursor = dmzprecursor,
                                      ppm_products = dmzproducts,
                                      rttol = rttol,
                                      coelCutoff = coelcutoff,
                                      lipidClasses = c("MG", "LPC", "LPE", "PC", "PE", "Sph", 
                                                       "SphP", "Cer", "SM", "Car", "CE", "DG", "TG"))
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
                                      lipidClasses = c("FA", "FAHFA", "LPC", "LPE", "LPG", "LPI", 
                                                       "LPS", "PC", "PE", "PG", "PI", "PS", "Sph", 
                                                       "SphP", "Cer", "CL", "BA"))
  }
}
save(annotated_msobjects, file = "annotated_msobjects.rda")

################################################################################
# output
load("annotated_msobjects.rda")

View(annotated_msobjects[[1]]$results)

# write tables to csv files
for (f in 1:length(annotated_msobjects)){
  fname <- gsub(".mzXML", "", annotated_msobjects[[f]]$metaData$generalMetadata$file)
  write.csv(annotated_msobjects[[f]]$results, file = paste0(fname, "_summaryResults.csv"),
            row.names = FALSE)
  write.csv(annotated_msobjects[[f]]$annotatedPeaklist, 
            file = paste0(fname, "_annotatedPeaklist.csv"),
            row.names = FALSE)
}

# plots
for (f in 1:length(annotated_msobjects)){
  annotated_msobjects[[f]] <- plotLipids(annotated_msobjects[[f]])
  
  fname <- gsub(".mzXML", "", annotated_msobjects[[f]]$metaData$generalMetadata$file)
  pdf(paste0(fname, ".pdf"))
  for (p in 1:length(annotated_msobjects[[f]]$plots)){
    print(annotated_msobjects[[f]]$plots[[p]])
  }
  dev.off()
}

save(annotated_msobjects, file = "annotated_msobjects.rda")

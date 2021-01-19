# sumChains
#' Calculate total number of carbons and double bounds of lipid chains
#'
#' Given the structure of a lipid specie, it sums up the chains.
#'
#' @param chains character value with the configuration of the chains separated
#' by a white space
#' @param n number of chains
#'
#' @return Character value indicating the total number of carbons and double
#' bounds
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
sumChains <- function(chains, n){
  cb <- unlist(strsplit(chains, "[: ]"))
  if (n == 2){
    sum_c <- as.numeric(cb[1])+as.numeric(cb[3])
    sum_b <- as.numeric(cb[2])+as.numeric(cb[4])
  } else if (n == 3){
    sum_c <- as.numeric(cb[1])+as.numeric(cb[3])+as.numeric(cb[5])
    sum_b <- as.numeric(cb[2])+as.numeric(cb[4])+as.numeric(cb[6])
  } else if (n == 4){
    sum_c <- as.numeric(cb[1])+as.numeric(cb[3])+as.numeric(cb[5])+as.numeric(cb[7])
    sum_b <- as.numeric(cb[2])+as.numeric(cb[4])+as.numeric(cb[6])+as.numeric(cb[8])
  }
  total <- paste(c(sum_c, sum_b), collapse=":")
  return(total)
}

# mzMatch
#' mz match withing a vector of mz values
#'
#' This function searches marches between a given mz and a vector of mz values
#' with certain mass  tolerance and returns the index of the matched values. It
#' is used by identification functions to find candidates of each class of lipid
#' based on full MS information.
#'
#' @param mz mz value to be matched
#' @param mzvector vector of mz values
#' @param ppm mass error tolerance
#'
#' @return Numeric vector indicating the index of matched mz values and ppms for
#' each one of those matches (match1, ppm1, match2, ppm2, etc.)
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
mzMatch <- function(mz, mzvector, ppm){
  matches_ppm <- vector()
  for (i in 1:length(mzvector)){
    ppm_observed <- abs(((mz-mzvector[i])/mz)*1000000)
    if (ppm_observed <= ppm){
      matches_ppm <- append(matches_ppm, c(i, ppm_observed))
    }
  }
  return(matches_ppm)
}

# cbs
#' Total number of carbons and double bounds
#'
#' This function matches mz values with neutral masses from a dataframe which
#' links masses and structures (carbons and double bounds) and extracts the
#' structural information. It is used by identification functions to look for
#' the structure of the previously chosen candidates.
#'
#' @param mz mz value to be matched
#' @param ppm mass error tolerance
#' @param db database
#' @param charge numeric value indicating the charge of the ion
#'
#' @return Character value or vector indicating structural information
#' (carbons:bounds)
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
cbs <- function(mz, ppm, db, charge=0){
  cb <- as.vector(db[abs(db["Mass"]+charge-mz) <
      ppm*(db["Mass"]+charge)/1000000, "total"])
  if (length(cb) > 1){
    cb <- cb[1]
  }
  if (length(cb) == 0){
    cb <- ""
  }
  return(cb)
}

# filtermsms
#' Presence or absence of an mz value withing a vector of mz values
#'
#' This function indicates the presence or absence of a fragment within a set
#' of mz values with certain tolerance. It is used by identification functions
#' to look for the generic fragments of each class of lipid.
#'
#' @param fragments vector of mz values
#' @param frag mz to be matched
#' @param ppm mass tolerance
#'
#' @return Logical value
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
filtermsms <- function(fragments, frag, ppm){
  sel <- fragments[which(abs((fragments$m.z - frag)/frag)*1000000 < ppm),]
  sel <- sel[which.max(sel$coelScore),]
  return(sel)
}

# frags
#' Search for fragments of interest withing a list of coeluting fragments
#'
#' Given a set of coeluting fragments, this function searches for matches within
#' a database. It is used by identification functions to extract fragments of
#' interest based on the fragmentation patterns of each class of lipid.
#'
#' @param df data frame containing coeluting fragments
#' @param ppm mass tolerance
#' @param db database (data frame with two columns) where to look into
#' @param charge mdiff
#'
#' @return Data frame containing matched ions information
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
frags <- function(df, ppm, db,  mdiff, charge, n){
  if (nrow(df) > 0){
    cb <- data.frame(0, 0, 0, 0, 0, 0)
    colnames(cb) <- c("cb", "m.z", "RT", "int", "peakID", "coelScore")
    found <- FALSE
    if (nrow(df) > 0){
      for (x in 1:nrow(df)){
        y <- abs((abs((n*db$Mass+mdiff)/charge)-df[x,"m.z"])*1000000/
                   abs((n*db$Mass+mdiff)/charge)) < ppm
        if (sum(y) > 0){
          cb <- rbind(cb, data.frame(c(cb = db[which(y == TRUE), "total"],
                                       df[x,]), stringsAsFactors = F))
          found <- TRUE
        }
      }
      if (found == FALSE){
        cb <- data.frame()
      } else {
        cb <- cb[2:nrow(cb),]
        if (sum(duplicated(cb$cb)) > 0){
          cb <- joinfrags(cb)
        }
      }
    }
  } else {
    cb <- data.frame()
  }
  return(cb)
}

# joinfrags
#' Join fragments information when several peaks of the same fragment are
#' coeluting with a unique candidate
#'
#' Function employed by \link{frags}.
#'
#' @param df data frame containing coeluting fragments
#'
#' @return Data frame
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
joinfrags <- function(df){
  new <- vector()
  for (f in unique(df$cb)){
    subset <- df[df$cb == f,]
    subset$int <- rep(sum(subset$int), nrow(subset))
    new <- rbind(new, subset[which.max(subset$coelScore),])
  }
  return(new)
}

# diffcb
#' Difference between two carbon:bounds structures
#'
#' This function calculates the number of carbon and double bounds that differ
#' between two structures.
#'
#' @param total character value indicating the precursor structure
#' @param frag character value indicating the fragment structure
#' @param db db of chains to be considered
#'
#' @return Character value
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
diffcb <- function(total, frag, db){
  fas <- db$total
  total <- unlist(strsplit(total, ":"))
  fr <- unlist(sapply(frag, strsplit, ":"))
  x <- as.numeric(total) - as.numeric(fr)
  frags <- paste(x[seq(1, length(x), 2)], x[seq(2, length(x), 2)], sep=":")
  del <- setdiff(frags, fas)
  frags[frags %in% del] <- ""
  if (length(frags) > 0){
    return(as.vector(frags))
  } else {
    return("")
  }
}

# select
#' Check matches between chains composition and precursor structures
#'
#' This function checks if the sum up of the chains structure match the
#' precursor structure. It is used by \link{combineChains}.
#'
#' @param chains character value containing chains structure separated by a
#' white space.
#' @param parent precursor ion structure
#' @param n number of chains
#'
#' @return Logical value
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
select <- function (chains, parent, n){
  chains <- unlist(strsplit(chains, "[: ]"))
  if (n == 3) {
    sum_c <- as.numeric(chains[1]) + as.numeric(chains[3]) +
      as.numeric(chains[5])
    sum_e <- as.numeric(chains[2]) + as.numeric(chains[4]) +
      as.numeric(chains[6])
  } else if (n == 2) {
    sum_c <- as.numeric(chains[1]) + as.numeric(chains[3])
    sum_e <- as.numeric(chains[2]) + as.numeric(chains[4])
  } else if (n == 4){
    sum_c <- as.numeric(chains[1])+as.numeric(chains[3])+as.numeric(chains[5])+
      as.numeric(chains[7])
    sum_e <- as.numeric(chains[2])+as.numeric(chains[4])+as.numeric(chains[6])+
      as.numeric(chains[8])
  }
  equal <- paste(c(sum_c, sum_e), collapse = ":") == parent
  return(equal)
}

# findPrecursor
#' Find candidate precursor from fullMS function
#'
#' This function is employed by all identification function in this package to
#' find possible precursors in the fullMS function.
#'
#' @param MS1 data frame containing m/z, RT and intensity of ions from the first
#' function
#' @param db database to be searched for matches
#' @param massdif mass difference between neutral mass and the adduct expected
#' @param rt rt window
#' @param n numeric value indicating whether to look for monomers (1), dimers
#' (2), etc.
#' @param charge numeric value indicating the expected charge of the ions
#'
#'
#' @return Subset of the original data frame without adducts
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
findPrecursor <- function(MS1, db, ppm, massdif, rt, n=1, charge=1){
  precursors <- unlist(lapply(abs(n*db$Mass+massdif)/abs(charge), mzMatch,
                              MS1$m.z, ppm))
  if (length(precursors) > 0){
    matches <- precursors[seq(1, length(precursors), 2)]
    ppms <- precursors[seq(2, length(precursors), 2)]
    prec <- MS1[matches,]
    if (class(prec) == "numeric"){
      prec <- as.data.frame(t(prec))
    }
    ppms <- ppms[prec$RT >= rt[1] & prec$RT <= rt[2]]
    prec <- prec[prec$RT >= rt[1] & prec$RT <= rt[2],]
    if (nrow(prec) > 0){
      if (n == 1){
        cb <- sapply(prec$m.z*abs(charge), cbs, ppm, db, massdif)
      } else if (n == 2){
        cb <- vector()
        for (i in 1:nrow(prec)){
          cb <- append(cb, cbs((prec$m.z[i]*abs(charge)-massdif)/2, ppm, db))
        }
      }
      # joining info
      allprec <- cbind(prec, ppms, as.data.frame(cb, stringsAsFactors = F))
      return(allprec)
    } else {
      return(data.frame())
    }
  } else {
    return(data.frame())
  }
}


# crossAdducts
#' Cross different candidates tables to remove false positives.
#'
#' This function crosses tables of precursor candidates identified using different
#' adducts.
#'
#' @param df1 data frame containing identification results using the main adduct
#' @param df2 data frame containing identification results using a secondary
#' adduct
#' @param rttol retention time tolerance in seconds
#' @param rawData raw scans data. Output of \link{dataProcessing} function
#' (MS1$rawData).
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied.
#'
#' @return Subset of the original data frame without adducts
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
crossAdducts <- function(df1, df2, rttol, rawData, coelCutoff){
  if (nrow(df2) > 0 & nrow(df1) == 0){
    df <- df2
    df2 <- data.frame()
  }
  if (nrow(df1) > 0 & nrow(df2) > 0){
    toremove <- vector()
    tokeep2 <- rep(TRUE, nrow(df2))
    for (m in 1:nrow(df1)){
      if (df1$peakID[m] %in% df2$peakID){
        sel <- df2[df2$peakID == df1$peakID[m],]
        matched <- which(abs(sel$RT - df1$RT) < rttol & sel$cb == df1$cb)
        if (length(matched) > 0){
          scores <- vector()
          for (s in 1:length(matched)){
            score <- coelutionScore(df1$peakID[matched[s]],
                                    peak2 = sel$peakID, rawData = rawData)
            scores <- append(scores, score)
          }
          if (any(scores >= coelCutoff)){
            toremove <- append(toremove, TRUE)
          } else {
            toremove <- append(toremove, FALSE)
          }
        } else {
          tokeep2[which(df2$peakID == df1$peakID[m])] <- FALSE
          toremove <- append(toremove, FALSE)
        }
      } else {
        toremove <- append(toremove, FALSE)
      }
    }
    df1 <- df1[!toremove,]
    df2 <- df2[tokeep2,]
    if (nrow(df1) > 0 & nrow(df2) > 0){
      df1$adducts <- as.vector(df1$adducts)
      ad1 <- df1$adducts
      df2$adducts <- as.vector(df2$adducts)
      ad2 <- df2$adducts[1]
      for (c in 1:nrow(df1)){
        common  <- which(df2$cb == df1$cb[c] & abs(df2$RT - df1$RT[c]) < rttol)
        if (length(common) > 0){
          df1$adducts[c] <- paste(c(ad1[c], ad2), collapse = ";")
          for (i in 1:length(common)){
            df2$adducts[common[i]] <- paste(c(ad2, ad1[c]), collapse = ";")
          }

        }
      }
    }
    df <- rbind(df1, df2)
    df <- filtrateAdducts(df)
  }
  return(df)
}

# filtrateAdducts
#' Remove low adduct supported candidates to avoid false positives.
#'
#' In case some feature has been annotated to different candidate species,
#' this function removes the one with less adducts assigned.
#'
#' @param df data frame containing candidates
#'
#' @return Subset of the original data frame
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
filtrateAdducts <- function(df){
  if (nrow(df) > 0){
    toremove <- c()
    for (c in 1:(nrow(df)-1)){
      same <- which(df[,"m.z"] == df[c, "m.z"] &
                      df[,"RT"] == df[c, "RT"])
      same <- same[same != c]
      if (length(same) > 0){
        same <- c(c, same)
        nadducts <- sapply(same, function(x){
          length(unlist(strsplit(df[x,"adducts"], ";")))-1
        })
        if (any(same[1] < same[2:length(same)])){
          toremove <- append(c, toremove)
        }
      }
    }
    if (length(toremove) > 0){
      df <- df[-toremove,]
    }
  }
  return(df)
}

# coelutionScore
#' calculate coelution score between two peaks
#'
#' Calculate coelution score between two peaks.
#'
#' @param peak1 character vector specifying the peakID of the first peak.
#' @param peak2 character vector specifying the peakID of the second peak.
#' @param rawData data frame with raw data for each scan. it need to have at
#' least 5 columns: m.z, RT, int, Scan (ordinal number for a given MS function)
#' and peakID (peakID to which it has been assigned).
#'
#' #' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
coelutionScore <- function(peak1, peak2, rawData){
  if (nrow(rawData) > 0){
    chrom1 <- rawData[rawData$peakID == peak1, c("int", "Scan")]
    chrom1 <- chrom1[order(chrom1$Scan),]
    pred1 <- tryCatch({stats::predict(stats::smooth.spline(chrom1$Scan, chrom1$int, spar = 0.4),
                               x = chrom1$Scan)},
                      error = function(e) {return(list(x = chrom1$Scan, y = chrom1$int))})
    if(length(pred1) > 0){
      chrom1 <- data.frame(int = pred1$y, Scan = floor(chrom1$Scan))
    }
    scores <- sapply(peak2, function(x) {
      chrom2 <- rawData[rawData$peakID == x, c("int", "Scan")]
      chrom2 <- chrom2[order(chrom2$Scan),]
      pred2 <- tryCatch({stats::predict(stats::smooth.spline(chrom2$Scan, chrom2$int, spar = 0.4),
                                 x = chrom2$Scan)},
                        error = function(e) {return(list(x = chrom2$Scan,
                                                         y = chrom2$int))})
      if (length(pred2) > 0 & length(pred1) > 0){
        chrom2 <- data.frame(int = pred2$y, Scan = floor(chrom2$Scan))
        merged <- merge(chrom1, chrom2, by="Scan")
      } else {
        merged <- data.frame()
      }
      if(nrow(merged) > 0){
        score <- cor(merged[,"int.x"], merged[,"int.y"])
        if (is.na(score)){
          score <- 0
        }
      } else {
        score <- 0
      }
      return(score)
    })
  } else {
    scores <- rep(0, length(peak2))
  }
  return(scores)
}

# checkIntRules
#' Check intensity rules
#'
#' Check intensity rules to confirm chains structure.
#'
#' @param intrules character vector specifying the fragments to compare. See
#' details.
#' @param rates character vector with the expected rates given as a string
#' (i.e. "3/1"). See details.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param nchains number of chains of the targeted lipid class.
#' @param combinations output of \link{combineChains}
#' @param sn list of chain fragments identified. Object fragments of the
#' \link{combineChains} output.
#'
#' @details This function will be employed when the targeted lipid class has
#' more than one chain.
#'
#' Taking PG subclass as an example, intensities of lysoPG fragments
#' (informative for sn1) can be employed to confirm the chains structure
#' (intrules = c("lysopg_sn1/lysopg_sn1")).
#' In this case, the intensity of lysoPG resulting from the loss of the FA chain
#' in sn2 is at least 3 times higher (rates = c("3/1")) than the lysoPG
#' resulting from the loss of the FA chain in sn1.
#'
#' For the intrules argument, "/" will be use to separate the fragments to
#' compare, and "_" will be use to indicate in which list of fragments we
#' need to look for their intensities. This will depend on the chain fragments
#' rules defined previiously.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
checkIntRules <- function(intrules, rates, intrequired, nchains, combinations,
                          sn){
  if (nchains == 1){
    passed <- FALSE
  } else if (length(intrules) == 0 | "Unknown" %in% intrules){
    if (nrow(combinations) > 0){
      passed <- rep(FALSE, nrow(combinations))
      for (c in 1:nrow(combinations)){
        if (nchains == 2){
          if (combinations[c,1] == combinations[c,2]){
            passed[c] <- TRUE
          }
        } else if (nchains == 3){
          if (combinations[c,1] == combinations[c,2] &
              combinations[c,1] == combinations[c,3]){
            passed[c] <- TRUE
          }
        } else if (nchains == 4){
          if (combinations[c,1] == combinations[c,2] &
              combinations[c,1] == combinations[c,3] &
              combinations[c,1] == combinations[c,4]){
            passed[c] <- TRUE
          }
        }
      }
    } else {
      return(F)
    }
  } else if (nchains == 2){
    if (nrow(combinations) == 0){
      return(F)
    } else {
      verified <- matrix(NA,  ncol=length(intrules), nrow=nrow(combinations))
      for (c in 1:nrow(combinations)){
        if (combinations[c,1] == combinations[c,2]){
          verified[c,] <-T
        } else {
          for (i in 1:length(intrules)){
            comp <- unlist(strsplit(intrules[i], "[_/]"))
            list <- sn[[c]]
            if (comp[2] == "sn1"){
              f1 <- 1
            } else {
              f1 <- 2
            }
            if (comp[4] == "sn1"){
              f2 <- 1
            } else {
              f2 <- 2
            }
            int1 <- sum(as.numeric(list$int[list$cb == combinations[c,f1] &
                                           list$db == comp[1]]))
            if (length(int1) == 0){int1 <- 0} else {
              int1 <- int1/sum(combinations[c,] == combinations[c,f1])
            }
            int2 <- sum(as.numeric(list$int[list$cb == combinations[c,f2] &
                                          list$db == comp[3]]))
            if (length(int2) == 0){int2 <- 0} else {
              int2 <- int2/sum(combinations[c,] == combinations[c,f2])
            }
            if (length(int1) == 1 & length(int2) == 1){
              if (int1 == 0 & int2 == 0){
                verified[c,i] <- F
              } else {
                if (eval(parse(text=rates[i])) >= 1){
                  check <- int1/int2 >= eval(parse(text=rates[i]))
                } else {
                  check <- int1/int2 <= eval(parse(text=rates[i]))
                }
                verified[c,i] <- check
              }
            } else {
              if (eval(parse(text=rates[i])) >= 1){
                check <- any(int1[1]/int2 >= eval(parse(text=rates[i])))
              } else {
                check <- any(int1[1]/int2 < eval(parse(text=rates[i])))
              }
              verified[c,i] <- check
            }
          }
        }
      }
      if (any(intrequired)){ # when there are required fragments, we check those
        passed <- unlist((apply(verified, 1, function(x){
          sum(x)
        }))) >= sum(intrequired)
      } else { # if there isnt any required fragment, any of them will be enough
        passed <- apply(verified, 1, sum) > 0
      }
    }
  } else if (nchains == 3){
    if (nrow(combinations) == 0){
      return(F)
    } else {
      verified <- matrix(NA,  ncol=length(intrules), nrow=nrow(combinations))
      for (c in 1:nrow(combinations)){
        if (combinations[c,1] == combinations[c,2] &
            combinations[c,1] == combinations[c,3]){
          verified[c,] <- T
        } else {
          for (i in 1:length(intrules)){
            comp <- unlist(strsplit(intrules[i], "[_/]"))
            list <- sn[[c]]
            if (comp[2] == "sn1"){
              f1 <- 1
            } else if (comp[2] == "sn2"){
              f1 <- 2
            } else {
              f1 <- 3
            }
            if (comp[4] == "sn1"){
              f2 <- 1
            } else if (comp[4] == "sn2"){
              f2 <- 2
            } else {
              f2 <- 3
            }
            int1 <- sum(as.numeric(list$int[list$cb == combinations[c,f1] &
                                              list$db == comp[1]]))
            if (length(int1) == 0){int1 <- 0} else {
              int1 <- int1/sum(combinations[c,] == combinations[c, f1])
            }
            int2 <- sum(as.numeric(list$int[list$cb == combinations[c,f2] &
                                              list$db == comp[3]]))
            if (length(int2) == 0){int2 <- 0} else {
              int2 <- int2/sum(combinations[c,] == combinations[c, f2])
            }
            if (length(int1) == 1 & length(int2) == 1){
              if (int1 == 0 & int2 == 0){
                verified[c,i] <- F
              } else {
                if (eval(parse(text=rates[i])) >= 1){
                  check <- int1/int2 >= eval(parse(text=rates[i]))
                } else {
                  check <- int1/int2 <= eval(parse(text=rates[i]))
                }
                verified[c,i] <- check
              }
            } else {
              if (eval(parse(text=rates[i])) >= 1){
                check <- any(int1[1]/int2 >= eval(parse(text=rates[i])))
              } else {
                check <- any(int1[1]/int2 < eval(parse(text=rates[i])))
              }
              verified[c,i] <- check
            }
          }
        }
      }
      if (any(intrequired)){ # when there are required fragments, we check those
        passed <- unlist((apply(verified, 1, function(x){
          sum(x)
        }))) >= sum(intrequired)
      } else { # if there isnt any required fragment, any of them will be enough
        passed <- apply(verified, 1, sum) > 0
      }
    }
  }
  return(passed)
}

# findMS2precursor
#' find lisnks between MS1 peaks and precursors selected for MS2 in DDA
#'
#' Find lisnks between MS1 peaks and precursors selected for MS2 in DDA.
#'
#' @param mz mz
#' @param minint minimum intensity
#' @param maxint maximum intensity
#' @param precursors data frame with all precursors
#' @param ppm mass tolerance
#'
#' @returns peak-pick based on previous EIC clusters generated by \link{clustering}
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
findMS2precursor <- function(mz, minrt, maxrt, precursors, ppm){
  fprec <- precursors[which(precursors$retentionTime >= minrt &
                              precursors$retentionTime <= maxrt),]
  if (nrow(fprec) > 0){
    matches <- mzMatch(mz, fprec$precursor, ppm)
    if (length(matches) > 0){
      scans <- fprec$Scan[matches[seq(1, length(matches), 2)]]
    } else {
      scans <- c()
    }
  } else {
    scans <- c()
  }
  return(scans)
}


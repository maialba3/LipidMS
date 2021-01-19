# MassSph
#' Calculate formula and mass of sphingoid bases
#'
#' Calculate formula and mass of sphingoid bases
#'
#' @param Sph character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassSph <- function(Sph){
  cb <- unlist(strsplit(Sph, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2+3-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO2"),
    collapse="")
  mass <- C*12+H*1.007825+14.00305+2*15.9949
  return(c(formula, mass))
}

# MassSphP
#' Calculate formula and mass of sphingoid phosphate bases
#'
#' Calculate formula and mass of sphingoid phosphate bases
#'
#' @param SphP character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassSphP <- function(SphP){
  cb <- unlist(strsplit(SphP, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2+4-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO5P"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+5*15.9949+30.97393
  return(c(formula, mass))
}

# MassCer
#' Calculate formula and mass of ceramides
#'
#' Calculate formula and mass of ceramides
#'
#' @param cer character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassCer <- function(cer){
  cb <- unlist(strsplit(cer, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2+1-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO3"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+3*15.9949
  return(c(formula, mass))
}

# MassCerP
#' Calculate formula and mass of ceramides phosphate
#'
#' Calculate formula and mass of ceramides phosphate
#'
#' @param cerP character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassCerP <- function(cerP){
  cb <- unlist(strsplit(cerP, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2+2-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO6P"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+6*15.9949+30.97393
  return(c(formula, mass))
}

# MassGlcCer
#' Calculate formula and mass of glucoceramides
#'
#' Calculate formula and mass of glucoceramides
#'
#' @param glccer character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassGlcCer <- function(glccer){
  cb <- unlist(strsplit(glccer, "[: ]"))
  C <- as.numeric(cb[1])+6
  H <- as.numeric(cb[1])*2+11-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO8"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+8*15.9949
  return(c(formula, mass))
}

# MassSM
#' Calculate formula and mass of sphingomyelines
#'
#' Calculate formula and mass of sphingomyelines
#'
#' @param SM character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassSM <- function(SM){
  cb <- unlist(strsplit(SM, "[: ]"))
  C <- as.numeric(cb[1])+5
  H <- as.numeric(cb[1])*2+13-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "N2O6P"),
                   collapse="")
  mass <- C*12+H*1.007825+2*14.00305+6*15.9949+30.97393
  return(c(formula, mass))
}

# MassCarnitine
#' Calculate formula and mass of carnitines
#'
#' Calculate formula and mass of carnitines
#'
#' @param carnitine character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassCarnitines <- function(carnitine){
  cb <- unlist(strsplit(carnitine, "[: ]"))
  C <- as.numeric(cb[1])+7
  H <- as.numeric(cb[1])*2+13-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO4"),
                   collapse="")
  mass <- C*12+H*1.007825++14.00305+4*15.9949
  return(c(formula, mass))
}

# MassCE
#' Calculate formula and mass of cholesterol esthers
#'
#' Calculate formula and mass of cholesterol esthers
#'
#' @param CE character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassCE <- function(CE){
  cb <- unlist(strsplit(CE, "[: ]"))
  C <- as.numeric(cb[1])+27
  H <- as.numeric(cb[1])*2+44-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O2"),
                   collapse="")
  mass <- C*12+H*1.007825+2*15.9949
  return(c(formula, mass))
}

# MassFA
#' Calculate formula and mass of fatty acids
#'
#' Calculate formula and mass of fatty acids
#'
#' @param FA character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassFA <- function(FA){
  cb <- unlist(strsplit(FA, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O2"),
                   collapse="")
  mass <- C*12+H*1.007825+2*15.9949
  return(c(formula, mass))
}

# MassHFA
#' Calculate formula and mass of hydroxi fatty acids
#'
#' Calculate formula and mass of hydroxi fatty acids
#'
#' @param HFA character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassHFA <- function(HFA){
  cb <- unlist(strsplit(HFA, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O3"),
                   collapse="")
  mass <- C*12+H*1.007825+3*15.9949
  return(c(formula, mass))
}

# MassFAHFA
#' Calculate formula and mass of FAHFA
#'
#' Calculate formula and mass of FAHFA
#'
#' @param FAHFA character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassFAHFA <- function(FAHFA){
  cb <- unlist(strsplit(FAHFA, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2-2-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O4"),
                   collapse="")
  mass <- C*12+H*1.007825+4*15.9949
  return(c(formula, mass))
}

# MassLysoPA
#' Calculate formula and mass of LPA
#'
#' Calculate formula and mass of LPA
#'
#' @param LPA character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPA <- function(LPA){
  cb <- unlist(strsplit(LPA, "[: ]"))
  C <- as.numeric(cb[1])+3
  H <- as.numeric(cb[1])*2+7-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H),
                     "O7P"), collapse="")
  mass <- C*12+H*1.007825+7*15.9949+30.97393
  return(c(formula, mass))
}

# MassLysoPAo
#' Calculate formula and mass of LPAo
#'
#' Calculate formula and mass of LPAo
#'
#' @param LPAo character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPAo <- function(LPAo){
  cb <- unlist(strsplit(LPAo, "[: ]"))
  C <- as.numeric(cb[1])+3
  H <- as.numeric(cb[1])*2+9-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H),
                     "O6P"), collapse="")
  mass <- C*12+H*1.007825+6*15.9949+30.97393
  return(c(formula, mass))
}

# MassPA
#' Calculate formula and mass of PA
#'
#' Calculate formula and mass of PA
#'
#' @param PA character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPA <- function(PA){
  cb <- unlist(strsplit(PA, "[: ]"))
  C <- as.numeric(cb[1])+3
  H <- as.numeric(cb[1])*2+5-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O8P"),
                   collapse="")
  mass <- C*12+H*1.007825+8*15.9949+30.97393
  return(c(formula, mass))
}

# MassLysoPE
#' Calculate formula and mass of LPE
#'
#' Calculate formula and mass of LPE
#'
#' @param LPE character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPE <- function(LPE){
  cb <- unlist(strsplit(LPE, "[: ]"))
  C <- as.numeric(cb[1])+5
  H <- as.numeric(cb[1])*2+12-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO7P"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+7*15.9949+30.97393
  return(c(formula, mass))
}

# MassLysoPEo
#' Calculate formula and mass of LPEo
#'
#' Calculate formula and mass of LPEo
#'
#' @param LPEo character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPEo <- function(LPEo){
  cb <- unlist(strsplit(LPEo, "[: ]"))
  C <- as.numeric(cb[1])+5
  H <- as.numeric(cb[1])*2+14-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO6P"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+6*15.9949+30.97393
  return(c(formula, mass))
}

# MassLysoPEp
#' Calculate formula and mass of LPEp
#'
#' Calculate formula and mass of LPEp
#'
#' @param LPEp character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPEp <- function(LPEp){
  cb <- unlist(strsplit(LPEp, "[: ]"))
  C <- as.numeric(cb[1])+5
  H <- as.numeric(cb[1])*2+12-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO6P"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+6*15.9949+30.97393
  return(c(formula, mass))
}

# MassPE
#' Calculate formula and mass of PE
#'
#' Calculate formula and mass of PE
#'
#' @param PE character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPE <- function(PE){
  cb <- unlist(strsplit(PE, "[: ]"))
  C <- as.numeric(cb[1])+5
  H <- as.numeric(cb[1])*2+10-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO8P"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+8*15.9949+30.97393
  return(c(formula, mass))
}

# MassPEo
#' Calculate formula and mass of plasmanyl PE
#'
#' Calculate formula and mass of plasmanyl PE
#'
#' @param PEp character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPEo <- function(PEo){
  cb <- unlist(strsplit(PEo, "[: ]"))
  C <- as.numeric(cb[1])+5
  H <- as.numeric(cb[1])*2+12-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO7P"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+7*15.9949+30.97393
  return(c(formula, mass))
}

# MassPEp
#' Calculate formula and mass of plasmenyl PE
#'
#' Calculate formula and mass of plasmenyl PE
#'
#' @param PEp character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPEp <- function(PEp){
  cb <- unlist(strsplit(PEp, "[: ]"))
  C <- as.numeric(cb[1])+5
  H <- as.numeric(cb[1])*2+10-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO7P"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+7*15.9949+30.97393
  return(c(formula, mass))
}

# MassLysoPG
#' Calculate formula and mass of LPG
#'
#' Calculate formula and mass of LPG
#'
#' @param LPG character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPG <- function(LPG){
  cb <- unlist(strsplit(LPG, "[: ]"))
  C <- as.numeric(cb[1])+6
  H <- as.numeric(cb[1])*2+13-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H),
                     "O9P"), collapse="")
  mass <- C*12+H*1.007825+9*15.9949+30.97393
  return(c(formula, mass))
}

# MassPG
#' Calculate formula and mass of PG
#'
#' Calculate formula and mass of PG
#'
#' @param PG character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPG <- function(PG){
  cb <- unlist(strsplit(PG, "[: ]"))
  C <- as.numeric(cb[1])+6
  H <- as.numeric(cb[1])*2+11-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H),
                     "O10P"), collapse="")
  mass <- C*12+H*1.007825+10*15.9949+30.97393
  return(c(formula, mass))
}

# MassLysoPI
#' Calculate formula and mass of LPI
#'
#' Calculate formula and mass of LPI
#'
#' @param LPI character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPI <- function(LPI){
  cb <- unlist(strsplit(LPI, "[: ]"))
  C <- as.numeric(cb[1])+9
  H <- as.numeric(cb[1])*2+17-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O12P"),
                   collapse="")
  mass <- C*12+H*1.007825+12*15.9949+30.97393
  return(c(formula, mass))
}

# MassPI
#' Calculate formula and mass of PI
#'
#' Calculate formula and mass of PI
#'
#' @param PI character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPI <- function(PI){
  cb <- unlist(strsplit(PI, "[: ]"))
  C <- as.numeric(cb[1])+9
  H <- as.numeric(cb[1])*2+15-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O13P"),
                   collapse="")
  mass <- C*12+H*1.007825+13*15.9949+30.97393
  return(c(formula, mass))
}

# MassPIP
#' Calculate formula and mass of PIP
#'
#' Calculate formula and mass of PIP
#'
#' @param PIP character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPIP <- function(PIP){
  cb <- unlist(strsplit(PIP, "[: ]"))
  C <- as.numeric(cb[1])+9
  H <- as.numeric(cb[1])*2+16-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O16P2"),
                   collapse="")
  mass <- C*12+H*1.007825+16*15.9949+2*30.97393
  return(c(formula, mass))
}

# MassPIP2
#' Calculate formula and mass of PIP2
#'
#' Calculate formula and mass of PIP2
#'
#' @param PIP2 character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPIP2 <- function(PIP2){
  cb <- unlist(strsplit(PIP2, "[: ]"))
  C <- as.numeric(cb[1])+9
  H <- as.numeric(cb[1])*2+17-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O19P3"),
                   collapse="")
  mass <- C*12+H*1.007825+19*15.9949+3*30.97393
  return(c(formula, mass))
}

# MassPIP3
#' Calculate formula and mass of PIP3
#'
#' Calculate formula and mass of PIP3
#'
#' @param PIP3 character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPIP3 <- function(PIP3){
  cb <- unlist(strsplit(PIP3, "[: ]"))
  C <- as.numeric(cb[1])+9
  H <- as.numeric(cb[1])*2+18-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O22P4"),
                   collapse="")
  mass <- C*12+H*1.007825+22*15.9949+4*30.97393
  return(c(formula, mass))
}

# MassLysoPS
#' Calculate formula and mass of LysoPS
#'
#' Calculate formula and mass of LysoPS
#'
#' @param LPS character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPS <- function(LPS){
  cb <- unlist(strsplit(LPS, "[: ]"))
  C <- as.numeric(cb[1])+6
  H <- as.numeric(cb[1])*2+12-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO9P"),
                   collapse="")
  mass <- C*12+H*1.007825+9*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassPS
#' Calculate formula and mass of PS
#'
#' Calculate formula and mass of PS
#'
#' @param PS character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPS <- function(PS){
  cb <- unlist(strsplit(PS, "[: ]"))
  C <- as.numeric(cb[1])+6
  H <- as.numeric(cb[1])*2+10-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO10P"),
                   collapse="")
  mass <- C*12+H*1.007825+10*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassLysoPC
#' Calculate formula and mass of LysoPC
#'
#' Calculate formula and mass of LysoPC
#'
#' @param LPC character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPC <- function(LPC){
  cb <- unlist(strsplit(LPC, "[: ]"))
  C <- as.numeric(cb[1])+8
  H <- as.numeric(cb[1])*2+18-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO7P"),
                   collapse="")
  mass <- C*12+H*1.007825+7*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassLysoPCp
#' Calculate formula and mass of LysoPCp
#'
#' Calculate formula and mass of LysoPCp
#'
#' @param LPCp character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPCp <- function(LPCp){
  cb <- unlist(strsplit(LPCp, "[: ]"))
  C <- as.numeric(cb[1])+8
  H <- as.numeric(cb[1])*2+18-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO6P"),
                   collapse="")
  mass <- C*12+H*1.007825+6*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassLysoPCo
#' Calculate formula and mass of LysoPCo
#'
#' Calculate formula and mass of LysoPCo
#'
#' @param LPCo character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPCo <- function(LPCo){
  cb <- unlist(strsplit(LPCo, "[: ]"))
  C <- as.numeric(cb[1])+8
  H <- as.numeric(cb[1])*2+20-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO6P"),
                   collapse="")
  mass <- C*12+H*1.007825+6*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassPC
#' Calculate formula and mass of PC
#'
#' Calculate formula and mass of PC
#'
#' @param PC character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPC <- function(PC){
  cb <- unlist(strsplit(PC, "[: ]"))
  C <- as.numeric(cb[1])+8
  H <- as.numeric(cb[1])*2+16-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO8P"),
                   collapse="")
  mass <- C*12+H*1.007825+8*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassPCp
#' Calculate formula and mass of PCp
#'
#' Calculate formula and mass of PCp
#'
#' @param PCp character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPCp <- function(PCp){
  cb <- unlist(strsplit(PCp, "[: ]"))
  C <- as.numeric(cb[1])+8
  H <- as.numeric(cb[1])*2+16-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO7P"),
                   collapse="")
  mass <- C*12+H*1.007825+7*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassPCo
#' Calculate formula and mass of PCo
#'
#' Calculate formula and mass of PCo
#'
#' @param PCo character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPCo <- function(PCo){
  cb <- unlist(strsplit(PCo, "[: ]"))
  C <- as.numeric(cb[1])+8
  H <- as.numeric(cb[1])*2+18-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO7P"),
                   collapse="")
  mass <- C*12+H*1.007825+7*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassMG
#' Calculate formula and mass of MG
#'
#' Calculate formula and mass of MG
#'
#' @param MG character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassMG <- function(MG){
  cb <- unlist(strsplit(MG, "[: ]"))
  C <- as.numeric(cb[1])+3
  H <- as.numeric(cb[1])*2+6-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O4"),
                   collapse="")
  mass <- C*12+H*1.007825+4*15.9949
  return(c(formula, mass))
}

# MassDG
#' Calculate formula and mass of DG
#'
#' Calculate formula and mass of DG
#'
#' @param DG character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassDG <- function(DG){
  cb <- unlist(strsplit(DG, "[: ]"))
  C <- as.numeric(cb[1])+3
  H <- as.numeric(cb[1])*2+4-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O5"),
                   collapse="")
  mass <- C*12+H*1.007825+5*15.9949
  return(c(formula, mass))
}

# MassTG
#' Calculate formula and mass of TG
#'
#' Calculate formula and mass of TG
#'
#' @param TG character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassTG <- function(TG){
  cb <- unlist(strsplit(TG, "[: ]"))
  C <- as.numeric(cb[1])+3
  H <- as.numeric(cb[1])*2+2-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O6"),
                   collapse="")
  mass <- C*12+H*1.007825+6*15.9949
  return(c(formula, mass))
}

# MassCL
#' Calculate formula and mass of CL
#'
#' Calculate formula and mass of CL
#'
#' @param CL character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassCL <- function(CL){
  cb <- unlist(strsplit(CL, "[: ]"))
  C <- as.numeric(cb[1])+9
  H <- as.numeric(cb[1])*2+14-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O17P2"),
                   collapse="")
  mass <- C*12+H*1.007825+17*15.9949+2*30.97393
  return(c(formula, mass))
}

# dbSphingolipids
#' Creation of a database for Cer, CerP, GlcCer and SM
#'
#' Creation of a database for Cer, CerP, GlcCer and SM
#'
#' @param chains character vector indicating the FAs to be employed
#' @param chains2 character vector indicating the sphingoid bases to be employed
#' @param lipid character value indication the class of lipid.
#'
#' @return data frame containing formula, mass and total number of carbons and
#' insaturations.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
dbSphingolipids <- function(chains, chains2, lipid){
  comb <- vector()
  for (i in chains2){
    comb <- append(comb, sapply(chains, paste, i, collapse=" "))
  }
  comb <- unique(comb)
  total  <- unique(unlist(lapply(comb, sumChains, n = 2)))
  if (lipid == "Cer"){
    fm <- lapply(total, MassCer)
  } else if (lipid == "GlcCer"){
    fm <- lapply(total, MassGlcCer)
  } else if (lipid == "CerP"){
    fm <- lapply(total, MassCerP)
  } else if (lipid == "SM"){
    fm <- lapply(total, MassSM)
  }
  db <- data.frame(formula=unlist(lapply(fm, "[[", 1)), total=total,
                   Mass=as.numeric(unlist(lapply(fm, "[[", 2))), stringsAsFactors = F)
  return(db)
}

# dbOneChain
#' Creation of a database for Carnitines, CE, FA, HFA, LPL, MG, sphingoid
#' bases and sphingoid bases phosphate.
#'
#' Creation of a database for Carnitines, CE, FA, HFA, LPL, MG, sphingoid
#' bases and sphingoid bases phosphate.
#'
#' @param chains character vector indicating the FAs to be employed
#' @param lipid character value indication the class of lipid.
#'
#' @return data frame containing formula, mass and total number of carbons and
#' insaturations
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
dbOneChain <- function(chains, lipid){
  fas <- unique(chains)
  if (lipid == "Carnitine"){
    fm <- lapply(fas, MassCarnitines)
  } else if (lipid == "CE"){
    fm <- lapply(fas, MassCE)
  } else if (lipid == "FA"){
    fm <- lapply(fas, MassFA)
  } else if (lipid == "HFA"){
    fm <- lapply(fas, MassHFA)
  } else if(lipid == "LPA"){
    fm <- lapply(fas, MassLysoPA)
  } else if(lipid == "LPE"){
    fm <- lapply(fas, MassLysoPE)
  } else if(lipid == "LPG"){
    fm <- lapply(fas, MassLysoPG)
  } else if(lipid == "LPI"){
    fm <- lapply(fas, MassLysoPI)
  } else if(lipid == "LPS"){
    fm <- lapply(fas, MassLysoPS)
  } else if(lipid == "LPC"){
    fm <- lapply(fas, MassLysoPC)
  } else if(lipid == "MG"){
    fm <- lapply(fas, MassMG)
  } else if(lipid == "Sph"){
    fm <- lapply(fas, MassSph)
  } else if(lipid == "SphP"){
    fm <- lapply(fas, MassSphP)
  }
  db <- data.frame(formula=unlist(lapply(fm, "[[", 1)), total=fas,
                   Mass=as.numeric(as.vector(unlist(lapply(fm, "[[", 2)))),
                   stringsAsFactors = F)
  return(db)
}

# dbTwoChains
#' Creation of a database for FAHFA, DG and PL.
#'
#' Creation of a database for FAHFA, DG and PL.
#'
#' @param chains character vector indicating the FAs to be employed
#' @param lipid character value indication the class of lipid.
#'
#' @return data frame containing formula, mass and total number of carbons and
#' insaturations
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
dbTwoChains <- function(chains, lipid){
  fas <- rep(chains, 2)
  comb <- unique(utils::combn(fas, 2, simplify=F))
  total  <- unique(unlist(lapply(comb, sumChains, n = 2)))
  if (lipid == "PE"){
    fm <- lapply(total, MassPE)
  } else if (lipid == "PG"){
    fm <- lapply(total, MassPG)
  } else if (lipid == "PI"){
    fm <- lapply(total, MassPI)
  } else if (lipid == "PIP"){
    fm <- lapply(total, MassPIP)
  } else if (lipid == "PIP2"){
    fm <- lapply(total, MassPIP2)
  } else if (lipid == "PIP3"){
    fm <- lapply(total, MassPIP3)
  } else if (lipid == "PC"){
    fm <- lapply(total, MassPC)
  } else if (lipid == "DG"){
    fm <- lapply(total, MassDG)
  } else if (lipid == "FAHFA"){
    fm <- lapply(total, MassFAHFA)
  } else if (lipid == "PS"){
    fm <- lapply(total, MassPS)
  } else if (lipid == "PA"){
    fm <- lapply(total, MassPA)
  }
  db <- data.frame(formula=unlist(lapply(fm, "[[", 1)), total=total,
                   Mass=as.numeric(as.vector(unlist(lapply(fm, "[[", 2)))),
                   stringsAsFactors = F)
  return(db)
}

# dbThreeChains
#' Creation of a database for TG.
#'
#' Creation of a database for TG.
#'
#' @param chains character vector indicating the FAs to be employed
#' @param lipid character value indication the class of lipid.
#'
#' @return data frame containing formula, mass and total number of carbons and
#' insaturations
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
dbThreeChains <- function(chains, lipid){
  fas <- rep(chains, 3)
  comb <- unique(utils::combn(fas, 3, simplify=F))
  total  <- unique(unlist(lapply(comb, sumChains, n = 3)))
  if (lipid == "TG"){
    fm <- lapply(total, MassTG)
  }
  db <- data.frame(formula=unlist(lapply(fm, "[[", 1)), total=total,
                   Mass=as.numeric(as.vector(unlist(lapply(fm, "[[", 2)))),
                   stringsAsFactors = F)
  return(db)
}

# dbFourChains
#' Creation of a database for C.
#'
#' Creation of a database CL.
#'
#' @param chains character vector indicating the FAs to be employed
#' @param lipid character value indication the class of lipid.
#'
#' @return data frame containing formula, mass and total number of carbons and
#' insaturations
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
dbFourChains <- function(chains, lipid){
  fas <- rep(chains, 4)
  comb <- unique(utils::combn(fas, 4, simplify=F))
  total  <- unique(unlist(lapply(comb, sumChains, n = 4)))
  if (lipid == "CL"){
    fm <- lapply(total, MassCL)
  }
  db <- data.frame(formula=unlist(lapply(fm, "[[", 1)), total=total,
                   Mass=as.numeric(as.vector(unlist(lapply(fm, "[[", 2)))),
                   stringsAsFactors = F)
  return(db)
}

# getFormula
#' Get formula and neutral mass for annotated compounds
#'
#' Get formula and neutral mass for annotated compounds.
#'
#' @param df data frame with the input results
#'
#' @return Data frame
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
getFormula <- function(df){
  lipidClass <- df["Class"]
  cdb <- df["CDB"]
  if (lipidClass == "BA"){
    db <- LipidMS::badb
  } else if (lipidClass == "Carnitine"){
    db <- LipidMS::carnitinesdb
  } else if (lipidClass == "Cer"){
    db <- LipidMS::cerdb
  } else if (lipidClass == "CL"){
    db <- LipidMS::cldb
  } else if (lipidClass == "DG"){
    db <- LipidMS::dgdb
  } else if (lipidClass == "CE"){
    db <- LipidMS::CEdb
  } else if (lipidClass == "FA"){
    db <- LipidMS::fadb
  } else if (lipidClass == "FAHFA"){
    db <- LipidMS::fahfadb
  } else if (lipidClass == "HFA"){
    db <- LipidMS::hfadb
  } else if (lipidClass == "LPC"){
    db <- LipidMS::lysopcdb
  } else if (lipidClass == "LPE"){
    db <- LipidMS::lysopedb
  } else if (lipidClass == "LPG"){
    db <- LipidMS::lysopgdb
  } else if (lipidClass == "LPI"){
    db <- LipidMS::lysopidb
  } else if (lipidClass == "LPS"){
    db <- LipidMS::lysopsdb
  } else if (lipidClass == "MG"){
    db <- LipidMS::mgdb
  } else if (lipidClass == "PC"){
    db <- LipidMS::pcdb
  } else if (lipidClass == "PE"){
    db <- LipidMS::pedb
  } else if (lipidClass == "PG"){
    db <- LipidMS::pgdb
  } else if (lipidClass == "PI"){
    db <- LipidMS::pidb
  } else if (lipidClass == "PS"){
    db <- LipidMS::psdb
  } else if (lipidClass == "SM"){
    db <- LipidMS::smdb
  } else if (lipidClass == "Sph"){
    db <- LipidMS::sphdb
  } else if (lipidClass == "SphP"){
    db <- LipidMS::sphPdb
  } else if (lipidClass == "TG"){
    db <- LipidMS::tgdb
  } else if (lipidClass == "TG"){
    db <- LipidMS::tgdb
  }
  index <- which(db$total == as.character(cdb))
  if(length(index) > 0){
    formula <- db$formula[index]
    Mn <- db$Mass[index]
  } else {
    tempdb <- createLipidDB(lipidClass,
                            chains = c(as.character(cdb),
                                       "0:0", "0:0", "0:0"),
                            chains2 = "0:0")[[1]]
    formula <- tempdb[tempdb$total == as.character(cdb), 1]
    Mn <- tempdb[tempdb$total == as.character(cdb), 3]
  }
  return(data.frame(Formula = as.character(formula),
                    Mn = as.numeric(Mn), stringsAsFactors = FALSE))
}



#include <list>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <vector>
#include <algorithm>


using namespace std;

extern "C"{
  
  
  //****************************************************************************
  //** Max distance between clusters
  //****************************************************************************
  SEXP clustdist(SEXP mins, 
                 SEXP maxs
  ){
    
    PROTECT(mins = AS_NUMERIC(mins));
    PROTECT(maxs = AS_NUMERIC(maxs));
    double *minsv;
    minsv = NUMERIC_POINTER(mins);
    double *maxsv;
    maxsv = NUMERIC_POINTER(maxs);
    int lengMins = LENGTH(mins);
    int lengMaxs = LENGTH(maxs);
    int leng = lengMins*lengMaxs;
    SEXP cdiff;
    PROTECT(cdiff = NEW_NUMERIC(leng));
    SETLENGTH(cdiff, leng);
    double *cdiffv;
    cdiffv = NUMERIC_POINTER(cdiff);
    for(int n = 0; n < leng; n++){
      *(cdiffv + n) = NA_REAL;
    }
    
    for (int x = 0; x < lengMins - 1; x++){
      for(int y = x+1; y < lengMaxs; y++){
        *(cdiffv + (x*lengMins+y)) = max(fabs(*(minsv + x) - *(maxsv + y)), fabs(*(minsv + y) - *(maxsv + x)));
      }
    }
    
    SETLENGTH(cdiff, leng);
    UNPROTECT(3);
    
    return cdiff;
  }
  
}




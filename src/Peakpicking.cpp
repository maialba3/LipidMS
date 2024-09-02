/*
Peakpicking functions imported from enviPick R package:
Partitioning, clustering & peak detection for LC-MS centroided data
author: Martin Loos, Martin.Loos@eawag.ch
*/

#include <list>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <vector>
#include "auxPeakPicking.h"

using namespace std;

extern "C"{

      /************************************************************************/
      /* agglomerative partitioning 2D ****************************************/
      /************************************************************************/
      SEXP agglom(     SEXP mz, /* must be sorted */
                       SEXP rt,
                       SEXP ppm,
                       SEXP dmz,
                       SEXP drt
                      ){


            PROTECT(mz = AS_NUMERIC(mz));
            PROTECT(rt = AS_NUMERIC(rt));
            PROTECT(ppm = AS_INTEGER(ppm));
            PROTECT(dmz = AS_NUMERIC(dmz));
            PROTECT(drt = AS_NUMERIC(drt));
            double *mass;
            mass = NUMERIC_POINTER(mz);
            double *ret;
            ret = NUMERIC_POINTER(rt);
            int ppm2 = INTEGER_VALUE(ppm);
            double dmass = NUMERIC_VALUE(dmz);
            double dret = NUMERIC_VALUE(drt);
            int n, m, p, k = 0, found = 0;
            double lowmass, highmass, lowret, highret;
            int leng = LENGTH(rt);
            SEXP outit;
            PROTECT(outit = NEW_INTEGER(leng));
            int *at;
            at = INTEGER_POINTER(outit);
            for(n = 0; n < leng; n++){*(at + n) = 0;}
            SETLENGTH(outit, leng);
            int *these;
            these = new int[leng];
            int *those;
            those = new int[leng];
            int untilthese = 0, untilthose = 0, atthese;

            for(n = 0; n < leng; n++){
                if(*(at + n) == 0){ /* unassigned entry ? */
                    k++; // increase partition ID
                    *(at + n) = k;
                    found = 1;
                    these[0] = n;
                    untilthese = 1;
                    atthese = 1;
                    while(found == 1){
                        found = 0;
                        if(atthese == 1){
                            untilthose = 0;
                            for(m = 0; m < untilthese; m++){
                                if(ppm2 == 1){
                                   lowmass = (*(mass + these[m]) - ((*(mass + these[m]) * dmass) / 1E6));
                                   highmass = (*(mass + these[m]) + ((*(mass + these[m]) * dmass) / 1E6));
                                }else{
                                   lowmass = (*(mass + these[m]) - (dmass));
                                   highmass = (*(mass + these[m]) + (dmass));
                                }
                                lowret = (*(ret + these[m]) - dret);
                                highret = (*(ret + these[m]) + dret);
                                if(these[m] > 0){
                                    for(p = (these[m] - 1); p >= 0; p--){ /* towards lower mass */
                                        if((*(mass + p)) <= (lowmass)){
                                            break;
                                        }else{
                                            if(*(at + p) == 0){
                                                if((*(ret + p) >= lowret) && (*(ret + p) <= highret)){
                                                    those[untilthose] = p;
                                                    untilthose++;
                                                    *(at + p) = k;
                                                    found = 1;
                                                }
                                            }
                                        }
                                    }
                                }
                                if(these[m] < (leng - 1)){
                                    for(p = (these[m] + 1); p < leng; p++){ /* towards higher mass */
                                        if((*(mass + p)) >= (highmass)){
                                            break;
                                        }else{
                                            if(*(at + p) == 0){
                                                if((*(ret + p) >= lowret) && (*(ret + p) <= highret)){
                                                    those[untilthose] = p;
                                                    untilthose++;
                                                    *(at + p) = k;
                                                    found = 1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            atthese = 0;
                        }else{
                            untilthese = 0;
                            for(m = 0; m < untilthose; m++){
                                if(ppm2 == 1){
                                   lowmass = (*(mass + those[m]) - ((*(mass + those[m]) * dmass) / 1E6));
                                   highmass = (*(mass + those[m]) + ((*(mass + those[m]) * dmass) / 1E6));
                                }else{
                                   lowmass = (*(mass + those[m]) - (dmass));
                                   highmass = (*(mass + those[m]) + (dmass));
                                }
                                lowret = (*(ret + those[m]) - dret);
                                highret = (*(ret + those[m]) + dret);
                                if(those[m] > 0){
                                    for(p = (those[m] - 1); p >= 0; p--){
                                        if((*(mass + p)) <= (lowmass)){
                                            break;
                                        }else{
                                            if(*(at + p) == 0){
                                                if((*(ret + p) >= lowret) & (*(ret + p) <= highret)){
                                                    these[untilthese] = p;
                                                    untilthese++;
                                                    *(at + p) = k;
                                                    found = 1;
                                                }
                                            }
                                        }
                                    }
                                }
                                if(those[m] < (leng-1)){
                                    for(p = (those[m] + 1); p < leng; p++){
                                       if((*(mass + p)) >= (highmass)){
                                           break;
                                        }else{
                                             if(*(at + p) == 0){
                                                 if((*(ret+p) >= lowret) & (*(ret + p) <= highret)){
                                                     these[untilthese] = p;
                                                     untilthese++;
                                                     *(at + p) = k;
                                                     found = 1;
                                                 }
                                             }
                                        }
                                    }
                                }
                            }
                            atthese = 1;
                        }
                    }
                }
            }

            delete[] these;
            delete[] those;
            SETLENGTH(outit, leng);
            UNPROTECT(6);
            return outit;

      }

	
      /************************************************************************/
      /* assemble indices for return value of agglom **************************/
      /************************************************************************/
      SEXP indexed(SEXP index, /* must be sorted */
				           SEXP intensity,
                   SEXP minpeak,
                   SEXP maxint,
                   SEXP maxindex
                   ){

            PROTECT(index = AS_NUMERIC(index));
            PROTECT(intensity = AS_NUMERIC(intensity));
            PROTECT(minpeak = AS_INTEGER(minpeak));
            PROTECT(maxint = AS_NUMERIC(maxint));
            PROTECT(maxindex = AS_NUMERIC(maxindex));
            double *ind;
            ind = NUMERIC_POINTER(index);
            double *inte;
            inte = NUMERIC_POINTER(intensity);
            int minpeak2 = INTEGER_VALUE(minpeak);
            double maxint2 = NUMERIC_VALUE(maxint);
			double tempmax = 0;
            int maxind = INTEGER_VALUE(maxindex);
			int leng = LENGTH(index);
            int n, from, to, counted, atind;
            SEXP outit;
            PROTECT(outit = allocMatrix(INTSXP, maxind, 3));
            int *at;
            at = INTEGER_POINTER(outit);
            for(n=0; n<(maxind*3); n++){
                *(at + n) = 0;
            }

            from = 1;
            to = 1;
            counted = 1;
            tempmax = *inte;
            atind = 0;
            for(n = 1; n < leng; n++){
                if(*(ind + n) != *(ind + n - 1)){
                	if( ((tempmax >= maxint2) || (counted >= minpeak2)) && (*(ind + n - 1) != 0) ){
                   		*(at + atind) = from;
                        *(at + maxind + atind) = to;
                   		*(at + (maxind*2) + atind) = counted;
                   		atind++;
                    };
                	from = (n + 1);
                	to = from;
                    counted = 1;
                    tempmax = *(inte + n);
                }else{
                	if(*(inte+n) > tempmax){
               	        tempmax = *(inte + n);
                	};
                    to++;
                    counted++;
                }
            }
            n--;
            if(  ((tempmax >= maxint2)||(counted >= minpeak2)) && (*(ind + n) != 0)  ){
                *(at + atind) = from;
                *(at + maxind + atind) = to;
                *(at + (maxind*2) + atind) = counted;
            }

            SETLENGTH(outit, maxind * 3);
            UNPROTECT(6);
            return outit;

      }


      /************************************************************************/
      /* assemble agglom part ID vector for measurements **********************/
      /************************************************************************/
      SEXP partID(       SEXP index,
                         SEXP leng
                         ){

            PROTECT(index = AS_INTEGER(index));
            PROTECT(leng = AS_INTEGER(leng));
            int lengit = LENGTH(index);
            lengit=lengit/3;
            int *indexed;
            indexed = INTEGER_POINTER(index);
            int leng2 = INTEGER_VALUE(leng);
            SEXP outit;
            PROTECT(outit = NEW_INTEGER(leng2));
            int *at;
            int n,m,id=1;
            at = INTEGER_POINTER(outit);
            for(n=0;n<leng2;n++){
                *(at+n) = 0;
            }

            for(n=0;n<lengit;n++){
                for((m=*(indexed+n));(m<=*(indexed+lengit+n));m++){
                    *(at+(m-1))=id;
                }
                id++;
            }

            UNPROTECT(3);
            return outit;
      }


		/************************************************************************/
		/* EIC clustering 2D ****************************************************/
		/************************************************************************/
		SEXP getEIC(   	SEXP mz,
                        SEXP RT,
                        SEXP intens,
                        SEXP orderedint,
                        SEXP orderedret,
                        SEXP dmzdens,
                        SEXP ppm2,
                        SEXP drtdens,
                        SEXP merged2
                        ){

			PROTECT(mz = AS_NUMERIC(mz));
			PROTECT(RT = AS_NUMERIC(RT));
			PROTECT(intens = AS_NUMERIC(intens));
			PROTECT(orderedint = AS_INTEGER(orderedint));
			PROTECT(orderedret = AS_INTEGER(orderedret));
			PROTECT(dmzdens = AS_NUMERIC(dmzdens));
			PROTECT(ppm2 = AS_INTEGER(ppm2));
			PROTECT(drtdens = AS_NUMERIC(drtdens));
			PROTECT(merged2 = AS_INTEGER(merged2));
			double *ret, *mass, *intensity;
			mass = NUMERIC_POINTER(mz);
			ret = NUMERIC_POINTER(RT);
			intensity = NUMERIC_POINTER(intens);
			int *ordint, *ordret;
			ordint = INTEGER_POINTER(orderedint);
			ordret = INTEGER_POINTER(orderedret);
			double dmzdens2 = NUMERIC_VALUE(dmzdens);
			int ppm3 = INTEGER_VALUE(ppm2);
			double drtdens_1 = NUMERIC_VALUE(drtdens);
			int merged3 = INTEGER_VALUE(merged2);
			int leng = LENGTH(RT);
			int m, n, i, k = 0, clustnumb, maxat = 0, maxit = 0;
			double delmz;
			SEXP clusters;
			PROTECT(clusters = allocMatrix(REALSXP, leng, 13));
			double *clus;
			clus = REAL(clusters);
			for(m = 0;m < 13;m++){
				for(n = 0; n < leng; n++){
					clus[(m * leng) + n] = 0;
				}
			}
			int *at; // vector to store cluster candidates
			at = new int[leng];
	
			/* initialize with most intense measurement ************************/
			clustnumb = 1;
			if(ppm3 == 1){delmz = ((dmzdens2**(mass + (*(ordint) - 1))) / 1e6);}else{delmz = dmzdens2;}
			clus[(0 * leng)] = (*(mass + (*(ordint) - 1)) - (2 * delmz)); 	// low mass boundary based on mass tolerance
			clus[(1 * leng)] = (*(mass + (*(ordint) - 1)) + (2 * delmz));	// high mass boundary based on mass tolerance
			clus[(2 * leng)] = (*(ret + (*(ordint) - 1)) - drtdens_1);  	// low RT boundary
			clus[(3 * leng)] = (*(ret + (*(ordint) - 1)) + drtdens_1);  	// high RT boundary
			clus[(4 * leng)] = 1;                                			// number of measurements
			clus[(5 * leng)] = *(mass + (*(ordint) - 1));            		// mass sum, later set to mean during merging
			clus[(6 * leng)] = clustnumb;                        			// cluster ID
			clus[(7 * leng)] = 0;                                			// merged (1) or not (0)?
			clus[(8 * leng)] = *(intensity + (*(ordint) - 1));       		// maximum intensity in a cluster
			clus[(9 * leng) + (*(ordint) - 1)] = clustnumb;          		// cluster ID for measurement
			clus[(10 * leng)] = 0;                               			// variance, set later when merged
			clus[(11 * leng)] = *(mass + (*(ordint) - 1));           		// lowest centroid mass in cluster
			clus[(12 * leng)] = *(mass + (*(ordint) - 1));           		// highest centroid mass in cluster

			/* assign all other peaks ***********************************************/
			for(n = 1; n < leng; n++){ // n = 0 assigned above
			   
				/* (a) check for possible fit to existing clusters *****************/
				maxat = 0;
				for(m = 0; m < clustnumb; m++){
					if(*(mass + (*(ordint + n) - 1)) <= clus[(1 * leng) + m]){ // below high mass boundary?
						if(*(mass + (*(ordint + n) - 1)) >= clus[(0 * leng) + m]){ // above low mass boundary?
							if(*(ret + (*(ordint + n) - 1)) >= clus[(2 * leng) + m]){ // above low RT boundary?
								if(*(ret + (*(ordint + n) - 1)) <= clus[(3 *leng) + m]){ // below high RT boundary?
									at[maxat] = m; // store candidate cluster
									maxat++;
								}
							}
						}
					}
				}
			   
				/* (b) not assignable - create new cluster ***************************/
				if(maxat == 0){
					clustnumb++;
					clus[(9 * leng) + (*(ordint + n) - 1)] = clustnumb;
					if(ppm3 == 1){delmz = ((dmzdens2**(mass + (*(ordint + n) - 1))) / 1e6);}else{delmz = dmzdens2;}
					clus[0 + (clustnumb - 1)] = (*(mass + (*(ordint + n) - 1)) - (2 * delmz));
					clus[(1 * leng) + (clustnumb - 1)] = (*(mass + (*(ordint + n) - 1)) + (2 * delmz));
					clus[(2 * leng) + (clustnumb - 1)] = (*(ret + (*(ordint + n) - 1)) - drtdens_1);
					clus[(3 * leng) + (clustnumb - 1)] = (*(ret + (*(ordint + n) - 1)) + drtdens_1);
					clus[(4 * leng) + (clustnumb - 1)] = 1;
					clus[(5 * leng) + (clustnumb - 1)] = *(mass + (*(ordint + n) - 1));
					clus[(6 * leng) + (clustnumb - 1)] = clustnumb;
					clus[(7 * leng) + (clustnumb - 1)] = 0;
					clus[(8 * leng) + (clustnumb - 1)] = *(intensity + (*(ordint + n) - 1));
					clus[(11 * leng) + (clustnumb - 1)] = *(mass + (*(ordint + n) - 1));
					clus[(12 * leng) + (clustnumb - 1)] = *(mass + (*(ordint + n) - 1));
					continue; // go to next centroid
				}
			   
				/* (c) assignable - but check for other centroids at same RT first **/
				for(i = 0; i < leng; i++){ // first, find same index for RT order as for intensity order for this centroid
					if(*(ordint + n) == *(ordret + i)){
						k = i;
						break;
					}
				}

				maxit = maxat;
				if(k > 0){ 
					for(i = (k - 1);i >= 0; i--){ // backward over RT-order
						if(  *(ret + (*(ordret + k) - 1)) == *(ret + (*(ordret + i) - 1))  ){ // same RT?
							if((clus[(9 * leng) + (*(ordret + i) - 1)]) == 0){ // not in cluster yet
								continue;
							}else{
								for(m = 0; m < maxat; m++){
									if(clus[(9 * leng) + (*(ordret + i) - 1)] == (at[m] + 1)){
										at[m] = -9999;
										maxit--;
										break;
									}
								}
							}
						}else{
							break;
						}
					}
				}

				if(k < (leng - 1)){ // forward over RT-order
					for(i = (k + 1); i < leng; i++){
						if(*(ret + (*(ordret + k) - 1)) == *(ret + (*(ordret + i) - 1))){ // same RT?
							if((clus[(9 * leng) + (*(ordret + i) - 1)]) == 0){ // not in cluster yet
								continue;
							}else{
								for(m = 0; m < maxat; m++){
									if(clus[(9 * leng) + (*(ordret + i) - 1)] == (at[m] + 1)){
										at[m] = -9999;
										maxit--;
										break;
									}
								}
							}
						}else{
                           break;
						}
					}
				}

				if(maxit < maxat){ // rewrite maxat to omit -9999 entries
					i = maxat;
					k = 0;
					for(m = 0; m < i; m++){
						if(at[m] == -9999){
							k++;
							maxat--;
						}else{
							at[m - k] = at[m];
						}
					}
				}
				if(maxat != maxit) Rprintf("\n ERROR: maxat != maxit");

				/* (d) not assignable because of same RT centroids - create new cluster */
				if(maxat == 0){
					clustnumb++;
					clus[(9 * leng) + (*(ordint + n) - 1)] = clustnumb;
					if(ppm3 == 1){delmz = ((dmzdens2**(mass + (*(ordint + n) - 1))) / 1e6);}else{delmz = dmzdens2;}
					clus[0 + (clustnumb - 1)] = (*(mass + (*(ordint + n) - 1)) - (2 * delmz));
					clus[(1 * leng) + (clustnumb - 1)] = (*(mass + (*(ordint + n) - 1)) + (2 * delmz));
					clus[(2 * leng) + (clustnumb - 1)] = (*(ret + (*(ordint + n) - 1)) - drtdens_1);
					clus[(3 * leng) + (clustnumb - 1)] = (*(ret + (*(ordint + n) - 1)) + drtdens_1);
					clus[(4 * leng) + (clustnumb - 1)] = 1;
					clus[(5 * leng) + (clustnumb - 1)] = *(mass + (*(ordint + n) - 1));
					clus[(6 * leng) + (clustnumb - 1)] = clustnumb;
					clus[(7 * leng) + (clustnumb - 1)] = 0;
					clus[(8 * leng) + (clustnumb - 1)] = *(intensity + (*(ordint + n) - 1));
					clus[(11 * leng) + (clustnumb - 1)] = *(mass + (*(ordint + n) - 1));
					clus[(12 * leng) + (clustnumb - 1)] = *(mass + (*(ordint + n) - 1));
					continue;
				}
			   
				/* (e) fits exactly to one cluster *************************************/
				if(maxat == 1){
					clus[(9 * leng) + (*(ordint + n) - 1)] = (at[0] + 1);
					if(ppm3 == 1){delmz = ((dmzdens2**(mass + (*(ordint + n) - 1))) / 1e6);}else{delmz = dmzdens2;}
					/* shrink lower & upper mass bounds */
					if(clus[(0 * leng) + (at[0])] < (*(mass + (*(ordint + n) - 1)) - (2 * delmz))){clus[(0 * leng) + (at[0])] = (*(mass + (*(ordint + n) - 1)) - (2 * delmz));}
					if(clus[(1 * leng) + (at[0])] > (*(mass + (*(ordint + n) - 1)) + (2 * delmz))){clus[(1 * leng) + (at[0])] = (*(mass + (*(ordint + n) - 1)) + (2 * delmz));}
					/* adapt lower & upper RT bounds */
					if(clus[(2 * leng) + (at[0])] > (*(ret + (*(ordint + n) - 1)) - drtdens_1)){clus[(2 * leng) + (at[0])] = (*(ret + (*(ordint + n) - 1)) - drtdens_1);}
					if(clus[(3 * leng) + (at[0])] < (*(ret + (*(ordint + n) - 1)) + drtdens_1)){clus[(3 * leng) + (at[0])] = (*(ret + (*(ordint + n) - 1)) + drtdens_1);}
					clus[(4 * leng) + (at[0])] = clus[(4 * leng) + (at[0])] + 1;
					clus[(5 * leng) + (at[0])] = clus[(5 * leng) + (at[0])] + *(mass + (*(ordint + n) - 1));
					if(*(mass + (*(ordint + n) - 1)) < clus[(11 * leng) + (at[0])]){clus[(11 * leng) + (at[0])] = *(mass + (*(ordint + n) - 1));};
					if(*(mass + (*(ordint + n) - 1)) > clus[(12 * leng) + (at[0])]){clus[(12 * leng) + (at[0])] = *(mass + (*(ordint + n) - 1));};
					continue;
				}
			   
				/* (f) fits to several clusters - m/z density clustering *******/
				if(maxat > 1){
					delmz = fabs(*(mass + (*(ordint + n) - 1)) - (clus[(5 * leng) + (at[0])] / clus[(4 * leng) + (at[0])] )); // set to mass difference to first cluster
					for(m = 1; m < maxat; m++){
						if(
							fabs(*(mass + (*(ordint + n) - 1)) - (clus[(5 * leng) + (at[m])] / clus[(4 * leng) + (at[m])] )) < delmz
						){
							at[0] = at[m];
							delmz = fabs(*(mass + (*(ordint + n) - 1)) - (clus[(5 * leng) + (at[0])] / clus[(4 * leng) + (at[0])] ));
						}
					}
					clus[(9 * leng) + (*(ordint + n) - 1)] = (at[0] + 1);
					if(ppm3 == 1){delmz = ((dmzdens2**(mass + (*(ordint + n) - 1))) / 1e6);}else{delmz = dmzdens2;}
					/* shrink lower & upper mass bounds */
					if(clus[(0 * leng) + (at[0])] < (*(mass + (*(ordint + n) - 1)) - (2 * delmz))){clus[(0 * leng) + (at[0])] = (*(mass + (*(ordint + n) - 1)) - (2 * delmz));}
					if(clus[(1 * leng) + (at[0])] > (*(mass + (*(ordint + n) - 1)) + (2 * delmz))){clus[(1 * leng) + (at[0])] = (*(mass + (*(ordint + n) - 1)) + (2 * delmz));}
					/* adapt lower & upper RT bounds */
					if(clus[(2 * leng) + (at[0])] > (*(ret + (*(ordint + n) - 1)) - drtdens_1)){clus[(2 * leng) + (at[0])] = (*(ret + (*(ordint + n) - 1)) - drtdens_1);}
					if(clus[(3 * leng) + (at[0])] < (*(ret + (*(ordint + n) - 1)) + drtdens_1)){clus[(3 * leng) + (at[0])] = (*(ret + (*(ordint + n) - 1)) + drtdens_1);}
					clus[(4 * leng) + (at[0])] = clus[(4 * leng) + (at[0])] + 1;
					clus[(5 * leng) + (at[0])] = clus[(5 * leng) + (at[0])] + *(mass + (*(ordint + n) - 1));
					if(*(mass + (*(ordint + n) - 1)) < clus[(11 * leng) + (at[0])]){clus[(11 * leng) + (at[0])] = *(mass + (*(ordint + n) - 1));};
					if(*(mass + (*(ordint + n) - 1)) > clus[(12 * leng) + (at[0])]){clus[(12 * leng) + (at[0])] = *(mass + (*(ordint + n) - 1));};
					continue;
				}
			}
			
			if((merged3 != 1) || (clustnumb == 1)){
               delete[] at;
               UNPROTECT(10);
               return(clusters);		   
			}
		   
		   
		   
		   
			/* *************************************************************** */
			/* merge isobaric mass cluster *************************************/

            /* (a) set mean mass ***********************************************/
            for(n = 0; n < clustnumb; n++){
                clus[(5 * leng) + (n)] = (clus[(5 * leng) + (n)] / clus[(4 * leng) + (n)]);
            }
			
            /* (b) define exclusion-list of cluster with same RT ***************/
            std::vector<int> blackone;
            std::vector<int> blacktwo;
            for(n = 1; n < leng; n++){ // over centroids alias ordret positions
                for(k = (n - 1); k >= 0; k--){
                    if(*(ret + (*(ordret + k) - 1)) == *(ret + (*(ordret + n) - 1))){ // identical RT
                        if( /* precheck: confine to overlapping mass + tolerance */
                            (clus[(0 * leng) + (int(clus[(9 * leng) + (*(ordret + n) - 1)]) - 1)] < clus[(1 * leng) + (int(clus[(9 * leng) + (*(ordret + k) - 1)]) - 1)]) &&
                            (clus[(1 * leng) + (int(clus[(9 * leng) + (*(ordret + n) - 1)]) - 1)] > clus[(0 * leng) + (int(clus[(9 * leng) + (*(ordret + k) - 1)]) - 1)])
                        ){   /* smallest cluster index always first */
                            if((clus[(9 * leng) + (*(ordret + k) - 1)]) < (clus[(9 * leng) + (*(ordret + n) - 1)])){
                                blackone.push_back(int(clus[(9 * leng) + (*(ordret + k) - 1)]));
                                blacktwo.push_back(int(clus[(9 * leng) + (*(ordret + n) - 1)]));
                            }else{
                                blackone.push_back(int(clus[(9 * leng) + (*(ordret + n) - 1)]));
                                blacktwo.push_back(int(clus[(9 * leng) + (*(ordret + k) - 1)]));
                            }
                        }
                    }else{
                        break; // breaks k - different RTs
                    }
                }
            }
            int blacksize = blackone.size(); // to check whether any exclusions were found at all
			
			
			
            /* (c) find mergeable clusters: indices & mass differences *********/
            /*     cross-check for overlaps in RT between cluster **************/
            std::vector<int> clusterone;		// store index of mergeable cluster 1
            std::vector<int> clustertwo;		// store index of mergeable cluster 2
            std::vector<double> clusterdiff; 	// store cluster difference in mean mass
            int doit;
            for(n = 0; n < (clustnumb - 1); n++){
                for(m = (n + 1); m < clustnumb; m++){
                    if( /* precheck: confine to overlapping mass + tolerance ... */
                        (clus[(0 * leng) + (n)] < clus[(1 * leng) + (m)]) &&
                        (clus[(1*leng) + (n)] > clus[(0 * leng) + (m)])
                    ){ 
                        if( /* ... and check for overlapping exact centroid masses */
                            (
                                (clus[(0 * leng) + (n)] <= clus[(11 * leng) + (m)]) &&
                                (clus[(1 * leng)+(n)] >= clus[(12 * leng) + (m)])
                            )||(
                                (clus[(0 * leng) + (m)] <= clus[(11 * leng) + (n)]) &&
                                (clus[(1 * leng) + (m)] >= clus[(12 * leng) + (n)])
                            )
                        ){
                            if(blacksize == 0){
                                clusterone.push_back(n + 1);
                                clustertwo.push_back(m + 1);
                                clusterdiff.push_back(fabs(clus[(5 * leng) + (n)] - clus[(5 * leng) + (m)]));
                            }else{
                                doit = 0;
                                for(k = 0; k < blacksize; k++){
                                    if(blackone[k] == (n + 1)){
                                        if(blacktwo[k] == (m + 1)){
                                            doit = 1;
                                            break;
                                        }
                                    }
                                }
                                if(doit == 0){
                                    clusterone.push_back(n + 1);
                                    clustertwo.push_back(m + 1);
                                    clusterdiff.push_back(fabs(clus[(5 * leng) + (n)] - clus[(5 * leng) + (m)]));
                                }
                            }
                        }
                    }
                }
            }
			
            /* (d) merge cluster ***********************************************/
            int mergesize = clusterone.size(), stay, gone;
            std::vector<int> clusterase;
            if(mergesize > 0){
                while(mergesize > 0){
                    
					/* find cluster pair closest in mean mass difference */
                    delmz = clusterdiff[0];
                    m = 0;
					if(mergesize > 1){
                        for(k = 1; k < mergesize; k++){
                            if(clusterdiff[k] < delmz){
                                m = k;
                                delmz = clusterdiff[k];
                            }
                        }
                  	}
                    
					/* reassign cluster IDs of centroids */
					for(n = 0; n < leng; n++){
                        if(clus[(9 * leng) + n] == clustertwo[m]){
                            clus[(9 * leng) + n] = clusterone[m];
                        }
                    }
					
					/* update cluster entries */
                    clus[(4 * leng) + (clusterone[m] - 1)] = (clus[(4 * leng) + (clusterone[m] - 1)] + clus[(4 * leng) + (clustertwo[m] - 1)]); /* number of measurements */
					clus[(6 * leng) + (clustertwo[m] - 1)] = clusterone[m]; /* cluster ID */
                    clus[(7 * leng) + (clustertwo[m] - 1)] = 1;            /* merged? */
                   
					/* delete that pair & all links to second (=merged) cluster & check blacklist  */
					stay = clusterone[m];					
					gone = clustertwo[m];
                    mergesize = clusterone.size();
                    n = 0;
                    while(n < mergesize){ /* remove merges */
                        if( (clusterone[n] == gone) || (clustertwo[n] == gone) ){
                            clusterone.erase(clusterone.begin() + n);
                            clustertwo.erase(clustertwo.begin() + n);
                            clusterdiff.erase(clusterdiff.begin() + n);
                            mergesize = clusterone.size();
                        }else{
                            n++;
                        }
                    }
                    if((blacksize > 0) && (mergesize > 0)){
                        clusterase.clear(); /* check blacklist-entries */
                        for(k = 0; k < blacksize; k++){
                            if(blackone[k] == gone){
                                clusterase.push_back(blacktwo[k]);
                            }
                            if(blacktwo[k] == gone){
                                clusterase.push_back(blackone[k]);
                            }
                        }
                        if(clusterase.size() > 0){ /* remove indirect blacklisted links */
                            for(m = 0; (unsigned) m < clusterase.size(); m++){
                                mergesize = clusterone.size();
                                n = 0;
                                while(n < mergesize){
                                    if( (clustertwo[n] == stay) || (clusterone[n] == clusterase[m])){
                                        clusterone.erase(clusterone.begin() + n);
                                        clustertwo.erase(clustertwo.begin() + n);
                                        clusterdiff.erase(clusterdiff.begin() + n);
                                        mergesize = clusterone.size();
                                    }else{
                                        n++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
			
			/* output ******************************************************/
			delete[] at;
			UNPROTECT(10);
			return(clusters);

		}


		/************************************************************************/
		/* gap filling, max for same RT *****************************************/
		/************************************************************************/
		SEXP gapfill(  	SEXP RT,
						SEXP intens,
						SEXP RTorder,  	/* increasing RT order */
						SEXP mz,
						SEXP index,		/* some centroid index */
						SEXP scans,		/* RTs, order increasing */
						SEXP drtsmall
					){

			PROTECT(RT = AS_NUMERIC(RT));
			PROTECT(intens = AS_NUMERIC(intens));
			PROTECT(RTorder = AS_INTEGER(RTorder));
			PROTECT(mz = AS_NUMERIC(mz));
			PROTECT(index = AS_NUMERIC(index));
			PROTECT(scans = AS_NUMERIC(scans));
			PROTECT(drtsmall = AS_NUMERIC(drtsmall));

			int leng2 = LENGTH(RT);
			int leng3 = LENGTH(scans);
			int n, m, l, k;
			double drt = NUMERIC_VALUE(drtsmall);
			double *ret, *inte, *mass, *scanned, *ind;
			ret = NUMERIC_POINTER(RT);
			inte = NUMERIC_POINTER(intens);
			mass = NUMERIC_POINTER(mz);
			scanned = NUMERIC_POINTER(scans);
			ind = NUMERIC_POINTER(index);
			int *ordret;
			ordret = INTEGER_POINTER(RTorder);

			SEXP ans;
			PROTECT(ans = allocMatrix(REALSXP, leng3, 10));
			double *rans;
			rans = REAL(ans);
			for(m = 0; m < leng3; m++){
               rans[m] = 0; 							// mass
               rans[m + (1 * leng3)] = 0; 				// intensity
               rans[m + (2 * leng3)] = *(scanned + m); 	// RT
               rans[m + (3 * leng3)] = 0; 				// index
               rans[m + (4 * leng3)] = 0; 				// Filter
               rans[m + (5 * leng3)] = 0; 				// 1st peak-pick
               rans[m + (6 * leng3)] = 0; 				// with peak criteria
               rans[m + (7 * leng3)] = 0; 				// baseline #1
               rans[m + (8 * leng3)] = 0; 				// baseline #2
               rans[m + (9 * leng3)] = 0; 				// 2nd peak-pick
			}

		   
			/* for identical RT, fill with intensity-closest *******************/
			l = 0; // start at first RT scan
			for(n = 0; n < leng2; n++){
				//for(m = 0; m < leng3; m++){
				for(m = l; m < leng3; m++){
					if(*(ret + (*(ordret + n) - 1)) == rans[m + (2 * leng3)]){
						if(rans[m + (1 * leng3)] == 0){ // no entry yet
							rans[m] = *(mass + (*(ordret + n) - 1));
							rans[m + (1 * leng3)] = *(inte + (*(ordret + n) - 1));
							rans[m + (3 * leng3)] = *(ind + (*(ordret + n) - 1));
						}else{				
							if( (rans[m + (1 * leng3)]) < *(inte + (*(ordret + n) - 1))){ // only fill in if intensity is higher
								rans[m] = *(mass + (*(ordret + n) - 1));
								rans[m + (1 * leng3)] = *(inte + (*(ordret + n) - 1));
								rans[m + (3 * leng3)] = *(ind + (*(ordret + n) - 1));					    
							}
						}
						l = m;
						break;
					};
				};
			};
		    
           /* gap-filling - interpolation *************************************/
           k = 0;
           for(m = k; m < (leng3 - 2); m++){
               if(rans[m] != 0){
                   for(n = (m + 1); n < leng3; n++){
                       if(rans[n] != 0){
                           break;
                       }
                   }
                   if((n - m) > 1){ /* at least one measurement in between */
                       if(fabs(rans[n + (2 * leng3)] - rans[m + (2 * leng3)]) <= drt){
                           for(l = (m + 1); l < n; l++){
                               rans[l + (1 * leng3)] = rans[m + (1 * leng3)] + (
                                   ((rans[n + (1 * leng3)] - rans[m + (1 * leng3)])) /
                                   (fabs(rans[n + (2 * leng3)] - rans[m + (2 * leng3)])) *
                                   (fabs(rans[l + (2 * leng3)] - rans[m + (2 * leng3)]))
                               );
                           }
                       }
                       k = n;
                   }else{
                       k = n;
                   }
               }
           }

           UNPROTECT(8);
           return ans;
      }


      /************************************************************************/
      /* peakpicking **********************************************************/
      /************************************************************************/
       SEXP pickpeak(   SEXP out1, /* generated in gapfill */
                        SEXP drtlarge,
                        SEXP drttotal,
                        SEXP minpeak,
                        SEXP recurs,
                        SEXP weight,
                        SEXP SB,
                        SEXP SN,
                        SEXP minint,
                        SEXP upint,
                        SEXP ended,
                        SEXP win
            ){

           PROTECT(out1 = AS_NUMERIC(out1));
           PROTECT(drtlarge = AS_NUMERIC(drtlarge));
           PROTECT(drttotal = AS_NUMERIC(drttotal));
           PROTECT(minpeak = AS_INTEGER(minpeak));
           PROTECT(recurs = AS_INTEGER(recurs));
           PROTECT(weight = AS_INTEGER(weight));
           PROTECT(SB = AS_NUMERIC(SB));
           PROTECT(SN = AS_NUMERIC(SN));
           PROTECT(minint = AS_NUMERIC(minint));
           PROTECT(upint = AS_NUMERIC(upint));
           PROTECT(ended = AS_INTEGER(ended));
           PROTECT(win = AS_INTEGER(win));

           double drt2 = NUMERIC_VALUE(drtlarge);
           double drt4 = NUMERIC_VALUE(drttotal);
           double weight2 = NUMERIC_VALUE(weight);
           double minint2 = NUMERIC_VALUE(minint);
           double upint2 = NUMERIC_VALUE(upint);
           double SB2 = NUMERIC_VALUE(SB);
           double SN2 = NUMERIC_VALUE(SN);
           int win2 = INTEGER_VALUE(win);
           int ended2 = INTEGER_VALUE(ended);
           int minpeak2 = INTEGER_VALUE(minpeak);
           int recurs2 = INTEGER_VALUE(recurs);
           int leng3 = LENGTH(out1);
           leng3 = leng3 / 10;

           double *rans;
           rans = REAL(out1);

           double often=0;

           /* 1st recursive peak detection ************************************/
           peakdetect(leng3, &often, recurs2, rans, weight2, drt4);
           if(often>0){
               /* check if peak-criteria fulfilled & filter *******************/
               peakcrit1(rans, leng3, minpeak2, SB2, minint2, upint2, ended2, drt2, &often);
           };
           if(often>0){
               /* subtract + interpolate + get baseline ***********************/
               /* 2nd peak criteria check *************************************/
               peakcrit2(rans, leng3, minpeak2, minint2, upint2, win2, often, SN2);
           }
           UNPROTECT(12);
           return out1;

      }
	  
}
/*
Auxilliary peakpicking functions imported from enviPick R package:
Partitioning, clustering & peak detection for LC-MS centroided data
author: Martin Loos, Martin.Loos@eawag.ch
*/

#ifndef auxPeakPicking_h
#define auxPeakPicking_h

inline void peakdetect(
     int leng3,
     double *often,
     int recurs2,
     double *rans,
     double weight2,
     double drt4){

     double *delint, *score1, *score2, *score3;
     delint= new double[leng3-1];
     score1= new double[leng3-1];
     score2= new double[leng3-1];
     score3= new double[leng3-1];
     double maxint;
     int p=0,m,n,from,to,optimfrom,optimto;
     often[0]=0;

     while((p!=-1)&(often[0]<recurs2)){
         p=-1; /* find most intensive peak */
         maxint=0;
         for(n=0;n<leng3;n++){
             if(rans[n+(5*leng3)]==0){
                 if(rans[n+(1*leng3)]>=maxint){
                     maxint=rans[n+(1*leng3)];
                     p=n;
                 }
             }
         }
         if(p!=-1){    /* any peaks left? */
             for(n=0;n<(leng3-1);n++){ /* reset */
                 score1[n]=0;
                 score2[n]=0;
                 score3[n]=0;
             }
             often[0]++;
             /* to the left ***************************************************/
             if(p>0){
                 /* get del_int************************************************/
                 from=p;
                 for(n=(p-1);n>=0;n--){
                     if(
                         (rans[n+(5*leng3)]==0)&
                         (fabs(rans[p+(2*leng3)]-rans[n+(2*leng3)])<=drt4)&
                         (rans[n+(1*leng3)]!=0)
                     ){
                         delint[n]=( (rans[n+1+(4*leng3)]-rans[n+(4*leng3)]) / fabs(rans[n+(2*leng3)]-rans[n+1+(2*leng3)]) );
                         from=n;
                     }else{
                         break;
                     }
                 }
                 if(from<p){
                     /* get scores ********************************************/
                     for(n=(p-1);n>=from;n--){
                         maxint=rans[n+(4*leng3)];
                         for(m=n;m<p;m++){
                             /* on score 1 ************************************/
                             score1[n]=score1[n]+delint[m];
                             /* on score 2 ************************************/
                             if(rans[m+(4*leng3)]>maxint){
                                 maxint=rans[m+(4*leng3)];
                             }
                             score2[n]=(score2[n]+(fabs(rans[m+(4*leng3)]-maxint)));
                         }
                         /* on score 3 ****************************************/
                         score3[n]=(score1[n]-(score2[n]*weight2));
                     }
                     optimfrom=from;
                     for(n=(p-1);n>=from;n--){
                         if(score3[n]>score3[optimfrom]){
                             optimfrom=n;
                         }
                     }
                 }else{
                     optimfrom=p;
                 }
             }else{
                 optimfrom=p;
             };
             for(n=optimfrom;n<=p;n++){
                 rans[n+(5*leng3)]=often[0];
             }
             /* to the right **************************************************/
             if(p<(leng3-1)){
                 /* get del_int************************************************/
                 to=p;
                 for(n=(p+1);n<leng3;n++){
                     if(
                         (rans[n+(5*leng3)]==0)&
                         (fabs(rans[p+(2*leng3)]-rans[n+(2*leng3)])<=drt4)&
                         (rans[n+(1*leng3)]!=0)
                     ){
                         delint[n-1]=( (rans[n-1+(4*leng3)]-rans[n+(4*leng3)]) / fabs(rans[n+(2*leng3)]-rans[n-1+(2*leng3)]) );
                         to=n;
                     }else{
                         break;
                     }
                 }
                 if(to>p){
                     /* get scores ********************************************/
                     for(n=(p+1);n<=to;n++){
                         maxint=rans[n+(4*leng3)];
                         for(m=n;m>p;m--){
                             /* on score 1 ************************************/
                             score1[n-1]=score1[n-1]+delint[m-1];
                             /* on score 2 ************************************/
                             if(rans[m+(4*leng3)]>maxint){
                                 maxint=rans[m+(4*leng3)];
                             }
                             score2[n-1]=(score2[n-1]+(fabs(rans[m+(4*leng3)]-maxint)));
                         }
                         /* on score 3 ****************************************/
                         score3[n-1]=(score1[n-1]-(score2[n-1]*weight2));
                     }
                     optimto=to;
                     for(n=(p+1);n<=to;n++){
                         if(score3[n-1]>score3[optimto-1]){
                             optimto=n;
                         }
                     }
                 }else{
                     optimto=p;
                 }
             }else{
                 optimto=p;
             };
             for(n=p;n<=optimto;n++){
                 rans[n+(5*leng3)]=often[0];
             }
         }
     }
     delete[] delint;
     delete[] score1;
     delete[] score2;
     delete[] score3;

}


inline void peakcrit1(
       double *rans,
       int leng3,
       int minpeak2,
       double SB2,
       double minint2,
       double upint2,
       double ended2,
       double drt2,
       double *often
       ){

		int n, m, to, p, dont, a, b;
		double maxint;
		bool found_end, found_minpeaks, evaluate_lower, evaluate_upper;
		int *traced1, *traced2, *traced3, *traced4, *peaked;
		traced1 = new int[int(often[0])];
		traced2 = new int[int(often[0])];
		traced3 = new int[int(often[0])];
		traced4 = new int[int(often[0])];
		peaked = new int[int(often[0])];
		for(n = 0; n < int(often[0]); n++){
			traced1[n] = 0;
			traced2[n] = 0;
			traced3[n] = 0;
			traced4[n] = 0;
			peaked[n] = 0;
		};

		/* extract criteria from peaks - new version ***************************/
		n = 0;
		while(n < leng3){
			if(rans[n + (5 * leng3)] != 0){ /* 1st pick - got start of peak = n */
				found_end = false;
				for(m = n; m < leng3; m++){ /* 1st pick - get end of peak = m */
					if(rans[n + (5 * leng3)] != rans[m + (5 * leng3)]){
						found_end = true;
						break;
					} /* else reached end < leng3 */
				}
				if(found_end){m--;}
				/*  -> evaluate minint2 and upint2 */
				maxint = 0;
				for(a = n; a <= m; a++){
					if(rans[a + (4 * leng3)] > maxint) maxint = rans[a + (4 * leng3)];
				}
				if(maxint >= minint2){
					traced3[int(rans[n + (5 * leng3)] - 1)] = 1;
				}
				if(maxint >= upint2){
					traced4[int(rans[n + (5 * leng3)] - 1)] = 1;
				}
				/*  -> evaluate minpeak within drt2 */
				found_minpeaks = false;
				if((m - n + 1) >= minpeak2){
					for(a = n; a <= (m - 1); a++){
						if(rans[a] != 0){
							p = 1;
							for(b = (a + 1); b <= m; b++){
								if(fabs(rans[a + (2 * leng3)] - rans[b + (2 * leng3)]) > drt2){
									break;
								}
								if(rans[b] != 0){
									p++;
								}
								if(p >= minpeak2){
									found_minpeaks = true;
									break;
								}
							}
						}
						if(found_minpeaks){break;}
					}
				}
				if(found_minpeaks){
					traced1[int(rans[n + (5 * leng3)] - 1)] = 1;
				}
				/* -> evaluate two-sided S/B, if next centroids on sides are not peak candidates themselves */
				evaluate_lower = true;
				if(n == 0){ /* do not evaluate lower S/B if rigth at RT start ... */
					evaluate_lower = false;
				}else{ /* ... or right adjacent to another peak candidate*/
					if(rans[(n - 1) + (5 * leng3)] != 0){evaluate_lower = false;}
				}

				evaluate_upper = true;
				if(m == (leng3 - 1)){ /* do not evaluate lower S/B if rigth at RT end ... */
					evaluate_upper = false;
				}else{ /* ... or right adjacent to another peak candidate*/
					if(rans[(m + 1) + (5 * leng3)] != 0){evaluate_upper = false;}
				}
				if(!evaluate_lower && !evaluate_upper){
					traced2[int(rans[n + (5 * leng3)] - 1)] = 1;
				}else{
					if(evaluate_lower && !evaluate_upper){
						if((maxint / rans[n + (4 * leng3)]) >= SB2){
							traced2[int(rans[n + (5 * leng3)] - 1)] = 1;
						}
					}else{
						if(!evaluate_lower && evaluate_upper){
							if((maxint / rans[m + (4 * leng3)]) >= SB2){
								traced2[int(rans[n + (5 * leng3)] - 1)] = 1;
							}
						}else{
							if(
								((maxint / rans[n + (4 * leng3)]) >= SB2) &&
								((maxint / rans[m + (4 * leng3)]) >= SB2)
							){
								traced2[int(rans[n + (5 * leng3)] - 1)] = 1;
							}

						}
					}

				}
				n = (m + 1); /* continue n after 1st pick */
			}else{
				n++;
			}
		}

       /* evaluate criteria - write peaks *************************************/
       to = 1;
       dont = 0;
       for(n = 0; n < int(often[0]); n++){
                   if(
                       ((traced1[n] == 1) && (traced2[n] == 1) && (traced3[n] == 1) && (dont < ended2)) ||
                       ((traced4[n] == 1) && (dont < ended2))
                   ){
                       peaked[n] = to;
                       to++;
                   }else{
                     peaked[n] = 0;
                     dont++;
                   }
       }
       for(n = 0; n < leng3; n++){
           if(rans[n + (5 * leng3)] != 0){
               if(peaked[int(rans[n + (5 * leng3)] - 1)] != 0){
                   rans[n + (6 * leng3)] = peaked[int(rans[n + (5 * leng3)] - 1)];
               }
           }
       }
       often[0] = to;

       delete[] traced1;
       delete[] traced2;
       delete[] traced3;
       delete[] traced4;
       delete[] peaked;

}


inline void peakcrit2(
       double *rans,
       int leng3,
       int minpeak2,
       double minint2,
       double upint2,
       int win2,
       double often,
       double SN2
       ){

       int n, m, to, p;
       double maxint,delint,delRT,meanint;
       int *traced1, *traced2, *traced3, *traced4, *peaked;
       traced1 = new int[int(often)];
       traced2 = new int[int(often)];
       traced3 = new int[int(often)];
       traced4 = new int[int(often)];
       peaked = new int[int(often)];
       for(n=0;n<int(often);n++){
               traced1[n]=0;
               traced2[n]=0;
               traced3[n]=0;
               traced4[n]=0;
               peaked[n]=0;
       };

       /* subtract peaks & fill in ********************************************/
       for(n=0;n<leng3;n++){
           if(rans[n+(6*leng3)]==0){
               rans[n+(7*leng3)]=rans[n+(4*leng3)];
           }else{
               for(m=n;m<leng3;m++){
                   if(rans[m+(6*leng3)]==0){
                       m=m-1;
                       break;
                   };
               };
               if((n>0)&(m<(leng3-1))){
                   delint=(rans[(m+1)+(4*leng3)]-rans[(n-1)+(4*leng3)]);
                   delRT=fabs(rans[(m+1)+(2*leng3)]-rans[(n-1)+(2*leng3)]);
                   delint=delint/delRT;
                   for(to=n;to<=m;to++){
                       rans[to+(7*leng3)]=( rans[(n-1)+(4*leng3)] +
                           ( delint * fabs(rans[(to)+(2*leng3)]-rans[(n-1)+(2*leng3)]) )
                       );
                   }
               }else{ /* special cases - peak at beginning or end */
                   if(n==0){
                       delint=(rans[(m+1)+(4*leng3)]);
                       delRT=fabs(rans[(m+1)+(2*leng3)]);
                       delint=delint/delRT;
                       for(to=n;to<=m;to++){
                           rans[to+(7*leng3)]=( 0 +
                               ( delint * rans[(to)+(2*leng3)] )
                           );
                       }
                   }
                   if(m==leng3){
                       delint=(-1*rans[(n-1)+(4*leng3)]);
                       delRT=fabs(rans[(leng3-1)+(2*leng3)]-rans[(n-1)+(2*leng3)]);
                       delint=delint/delRT;
                       for(to=n;to<=m;to++){
                           rans[to+(7*leng3)]=( rans[(n-1)+(4*leng3)] +
                               ( delint * fabs(rans[(to)+(2*leng3)]-rans[(n-1)+(2*leng3)]) )
                           );
                       }
                   } /* if both m,n -> 0 to 0 baseline*/
               }
               n=m;
           }
       }
       /* smooth baseline - window mean ***************************************/
       for(n=0;n<leng3;n++){
           maxint=0;
           p=0;
           for((m=(n-win2));(m<=(n+win2));m++){
               if((m>=0)&(m<leng3)){
                   maxint=maxint+rans[m+(7*leng3)];
                   p++;
               }
           }
           rans[n+(8*leng3)]=(maxint/double(p));
       }
       for(n=0;n<leng3;n++){
           rans[n+(7*leng3)]=rans[n+(8*leng3)];
           rans[n+(8*leng3)]=0;
       }
       /* subtract baseline ***************************************************/
       for(n=0;n<leng3;n++){
        rans[n+(8*leng3)]=rans[n+(4*leng3)]-rans[n+(7*leng3)];
       }
       for(n=0;n<leng3;n++){
           if(rans[n+(8*leng3)]<0){
               rans[n+(8*leng3)]=0;
           }
       }
       /* 2nd peak picking ****************************************************/
       /* define noise: mean deviation from baseline **************************/
       meanint = 0;
       p = 0;
       for(n = 0; n < leng3; n++){
			/* must have an intensity, must not be part of a peak candidate, baseline intensity must be present */
           if((rans[n + (4 * leng3)] != 0) && (rans[n + (6 * leng3)] == 0) && (rans[n + (7 * leng3)] != 0) ){
               meanint = meanint + rans[n + (8 * leng3)]; /* = baseline-corrected intensity */
               p++;
           }
       }
       if(p != 0){
           meanint = meanint / double(p);
       };
       /* extract criteria from peaks *****************************************/
       for(n=0;n<leng3;n++){
                   if(rans[n+(6*leng3)]!=0){
                       p=1;
                       maxint=0;
                       for(m=n;m<leng3;m++){
                           if(rans[m+(6*leng3)]==rans[n+(6*leng3)]){
                               if(rans[m+(8*leng3)]!=0){
                                   p++;
                                   if(rans[m+(8*leng3)]>maxint){
                                       maxint=rans[m+(8*leng3)];
                                   }
                               }
                           }else{
                               n=m;
                               break;
                           }
                       }
                       if(p>=minpeak2){
                           traced1[int(rans[(m-1)+(6*leng3)]-1)]=1;
                       }
                       if(meanint != 0){ /* any noise catched? */
                           if((maxint / meanint) >= SN2){
                               traced2[int(rans[(m-1)+(6*leng3)]-1)]=1;
                           }
                       }else{
                           traced2[int(rans[(m-1)+(6*leng3)]-1)]=1;
                       }
                       if(maxint>=minint2){
                           traced3[int(rans[(m-1)+(6*leng3)]-1)]=1;
                       }
                       if(maxint>=upint2){
                           traced4[int(rans[(m-1)+(6*leng3)]-1)]=1;
                       }

                   }
       }
       /* evaluate criteria - write peaks *************************************/
       to=1;
       for(n=0;n<int(often);n++){
                   if(
                       ((traced1[n]==1)&&(traced2[n]==1)&&(traced3[n]==1)) ||
                       (traced4[n]==1)
                   ){
                       peaked[n]=to;
                       to++;
                   }else{
                       peaked[n]=0;
                   }
       }
       for(n=0;n<leng3;n++){
           if(rans[n+(6*leng3)]!=0){
               if(peaked[int(rans[n+(6*leng3)]-1)]!=0){
                   rans[n+(9*leng3)]=peaked[int(rans[n+(6*leng3)]-1)];
               }
           }
       }

       delete[] traced1;
       delete[] traced2;
       delete[] traced3;
       delete[] traced4;
       delete[] peaked;

}


#endif


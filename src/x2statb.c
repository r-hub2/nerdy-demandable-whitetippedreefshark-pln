#ifndef R
#include <malloc.h>
#else
#include <R.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* version when not all categ vectors present */
/* X2 = sum_{a: O(a)>0} (Oa-Ea)^2/Ea  + [m^n-sum_{a:O(a)>0} E(a)]
*/
double x2statb(int n,int m, int nrec, double **dat, double *fr, double **alp, double *b,
  int nn, double *x, double *w, int nq)
{ int i,*ii,*kk,ic;
  double x2,epr,tem,sumexp;
  double plgngh(int, int, double **, double *, int, int *, int *,
    double *, double *, int);
  /*extern int nn;*/

#ifndef R
  ii=(int *) malloc((n+1) * sizeof(int));
  kk=(int *) malloc((n+1) * sizeof(int));
#else
  ii=(int *) Calloc((n+1), int);
  kk=(int *) Calloc((n+1), int);    
#endif
  for(i=1;i<=n;i++) ii[i]=i;
  for(ic=0,x2=0.,sumexp=0.;ic<nrec;ic++)
  { if(fr[ic]==0) continue;
    for(i=1;i<=n;i++) kk[i]=dat[ic][i-1];
    epr=plgngh(n,m,alp,b,n,ii,kk,x,w,nq);
    epr*=nn;
    sumexp+=epr;
    /*for(i=1;i<=n;i++) printf("%d ", kk[i]);*/
    /*printf("O=%d, E=%.1f\n", 1.,epr);*/
    tem=fr[ic]-epr;
    x2+=tem*tem/epr; 
  }
  x2+= nn - sumexp;
  /*printf("x2=%f, df(x2)=%d\n", x2,n2-np-1); */
#ifndef R
  free(ii); free(kk);
#else
  Free(ii); Free(kk);
#endif
  return x2;
}


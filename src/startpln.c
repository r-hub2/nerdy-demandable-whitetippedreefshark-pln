#ifndef R
#include <malloc.h>
#else
#include <R.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* summarize data into univariate margins in order to get starting
   point estimate for polytomous logit-normit model,
   based on assumption of independent item,
   input can include optional frequencies
   categories assumed to be 0 1 ... m-1
*/
#ifdef MAIN
/* gcc -DMAIN -o startpln startpln.c -lm */
/* startpln 5 3 67 < sim3fr.dat */
main(int argc, char *argv[])
{ int i,n,nn,m,j,nrec;
  /*int dat[NN][N],fr[NN];*/
  double **dat,*fr;
  /*void summfr(int, int, int, int, int [][N], int [ ], double [M*N]);*/
  void startpln(int, int, int, int, double **, double *, double *);
  /*double start[M*N];*/
  double *start;
  double **dmatrix(int,int);

  
  if(argc<4) 
  { printf("Usage: %s #items #categ #records [-fr]\n", argv[0]); exit(0); }
  n=atoi(argv[1]); m=atoi(argv[2]); nrec=atoi(argv[3]);
  dat=dmatrix(nrec,n);  /* reverse order in interface to matlab ? */
  fr=(double *) malloc((nrec) * sizeof(double));
  start=(double *) malloc((m-1)*n * sizeof(double));

  /* use index 0 for row and column */
  for(i=0,nn=0;i<nrec;i++)
  { for(j=0;j<n;j++) scanf("%lf", &dat[i][j]); 
    scanf("%lf", &fr[i]); nn+=fr[i]; 
  }

  printf("nn=%d\n", nn); startpln(n,m,nn,nrec,dat,fr,start); 
  

  /* also find starting points for betas based on correlation? */
  printf("\nstarting point for alphas for pln model\n");
  for(i=0;i<n*(m-1);i++) 
  { printf("%f ", start[i]); if((i+1)%5==0) printf("\n"); }
  printf("\n");
  free(start);
  free(fr);
  free(dat[0]); free(dat);
  exit(0);
}
#endif

/* get 1D margins,
 * nothing simpler than loops if #categories can be m>2
 */
void startpln(int n, int m, int nn, int nrec, double **dat, double *fr, 
   double *start)
{ int i,j,k,ii;
  /*int s[M],tot;*/
  double *s,tot;
  double alp;

#ifndef R
  s=(double *) malloc(m * sizeof(double));
#else
  s=(double *) Calloc(m, double);
#endif
  /* univariate margins */
  ii=0;
  for(j=0;j<n;j++)
  { for(k=0;k<m;k++) s[k]=0;
    for(i=0;i<nrec;i++) 
    { k=dat[i][j]; s[k]+=fr[i]; }
#ifdef MAIN
    printf("margin %d:\n", j+1);
#endif
    for(k=0,tot=0;k<m;k++) 
    { tot+=s[k]; 
      alp=tot/nn;
      alp=log((1.-alp)/alp);
#ifdef MAIN
      printf("%d %.0f %f\n", k, s[k], alp);
#endif
      if(k<m-1) { start[ii]=alp; ii++; }
    }
  }
#ifndef R
  free(s);
#else
  Free(s);   
#endif
}

#ifndef ALL
/* allocate matrices such that rows occupy contiguous space */
double **dmatrix(int nrow, int ncol)
{ int i; double **amat,*avec;
#ifndef R
  avec=(double *)malloc((unsigned) (nrow*ncol)*sizeof(double));
  amat=(double **)malloc((unsigned) nrow*sizeof(double*));
#else
  avec=(double *)Calloc((unsigned) (nrow*ncol), double);
  amat=(double **)Calloc((unsigned) nrow, double*);    
#endif
  for(i=0;i<nrow;i++) amat[i]=avec+i*ncol;
  return amat;
}

int **imatrix(int nrow, int ncol)
{ int i; int **amat,*avec;
#ifndef R
  avec=(int *)malloc((unsigned) (nrow*ncol)*sizeof(int));
  amat=(int **)malloc((unsigned) nrow*sizeof(int*));
#else
  avec=(int *)Calloc((unsigned) (nrow*ncol),int);
  amat=(int **)Calloc((unsigned) nrow, int*);    
#endif
  for(i=0;i<nrow;i++) amat[i]=avec+i*ncol;
  return amat;
}
#endif

#ifdef R
void Rstartpln( int *nitem, int *ncateg, int *nrec, double *dataset, double *testout)  

{
 double **dat,*fr; 
  double **dmatrix(int,int); 
  int i,j,nn; 
  void startpln(int, int, int, int, double **, double *, double *);

  /*printf("%d %d %d %f\n", nitem, ncateg, nrec,*aa); */
  /* convert to matrices in C (i.e. row/column transpose */ 
  dat=dmatrix(*nrec,*nitem); 
  fr=(double *) Calloc(*nrec, double);


  /* nn = total of fr[]  */
  for(i=0,nn=0;i<*nrec;i++)
  { for(j=0;j<*nitem;j++) dat[i][j] = *(dataset+(j*(*nrec)+i));  
    fr[i]=*(dataset+(i+*nitem*(*nrec))); nn+=fr[i]; 
    /*if(i<10) printf("%f %f\n", dat[i][nitem-1], fr[i]);*/ 
  }
  startpln(*nitem,*ncateg,nn,*nrec,dat,fr,testout);
  Free(fr); 
  Free(dat[0]); Free(dat); 
} 

# endif

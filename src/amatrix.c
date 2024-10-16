#ifndef R
#include <malloc.h>
#else
#include <R.h>
#endif

#include <stdlib.h>
/* allocate matrices such that rows occupy contiguous space */
double **dmatrix(int nrow, int ncol)
{ int i; double **amat,*avec;

#ifndef R
  avec=(double *)malloc((unsigned) (nrow*ncol)*sizeof(double));
  amat=(double **)malloc((unsigned) nrow*sizeof(double*));
#else
  avec=(double *) calloc((unsigned) (nrow*ncol), sizeof(double));
  amat=(double **) calloc((unsigned) nrow,sizeof(double*));   
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
  avec=(int *) calloc((unsigned) (nrow*ncol), sizeof(int));
  amat=(int **) calloc((unsigned) nrow, sizeof(int*));    
#endif
  for(i=0;i<nrow;i++) amat[i]=avec+i*ncol;
  return amat;
}

double ***gmatrix(int r, int c, int d)
{ int i; 
  double ***acube,**amat,*avec;
#ifndef R
  avec=(double *)malloc((unsigned) (r*c*d)*sizeof(double));
  amat=(double **)malloc((unsigned) (r*c)*sizeof(double*));
  acube=(double ***)malloc((unsigned) r*sizeof(double**));
#else
  avec=(double *) calloc((unsigned) (r*c*d), sizeof(double));
  amat=(double **) calloc((unsigned) (r*c), sizeof(double*));
  acube=(double ***) calloc((unsigned) r, sizeof(double**));
#endif
  for(i=0;i<(r*c);i++) amat[i]=avec+i*d;
  for(i=0;i<r;i++) acube[i]=amat+i*c;
  return acube;
}




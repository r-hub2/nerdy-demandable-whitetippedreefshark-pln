#ifndef R
#include <malloc.h>
#else
#include <R.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef MAXDOUBLE
#define MAXDOUBLE 1.797693134862315708e+308
#endif
/* version for nrbcpln.c which has one more argument in func than
   nrmlepln.c */

void nrminbcl(int np, double *val, 
    void (*func)(int, double *, double *, double **, int, 
      int, int, int, int, double *, double *, double **, double *,
      double ***, double ***, double ***, double *pp),
    double *lb, double *ub, int mxiter, double eps, int iprint, double bdd,
    int *iconv, double *fnval, double **invhes,
    int n, int m, int nq, int nrec, double *x, double *w, double **dat, double *fr,
    double ***g, double ***g1, double ***g2, double *pp)
/* np=dimension of problem, 
   val[]=starting point on entry, solution on completion
   func = function to be minimized with gradient and Hessian
   lb = lower bounds
   ub = upper bounds
   mxiter = max number of iterations
   eps = tolerance for Gaussian elimination, linear system solver
   iprint = control on amount of printing
             0 for no printing of iterations
             1 for printing x^{(k)} on each iteration
             2 for printing x^{(k)} and linear system on each iteration
   bdd is bound on difference of 2 consecutive iterations
   (useful is starting point is far from solution and function
    is far from convex)
   Return parameters:
   iconv = 1 if converged, 0 otherwise
   fnval = function value at minimum
   invhes[][] = inverse Hessian on return, allocate in calling routine
*/
{ double out;
  int j,k,iter,outb;
  double **dmatrix(int, int);
  void lsolve(double **, int n, int nk, double *det, double tol);
  double det,**dd,mxdif,tem;

  dd=dmatrix(np,np+1);
  iter=0; 
#ifndef R
  if(iprint>=2) printf("eps= %f\n", eps);
//#else
//  if(iprint==1) Rprintf("eps= %f\n", eps);
#endif
    
  do
  { func(np,val,&out,dd,iprint,
      n,m,nq,nrec,x,w,dat,fr,g,g1,g2,pp);
    lsolve(dd,np,np+1,&det,1.e-6);
    outb=0;
    for(j=0,mxdif=0.;j<np;j++) 
    { tem=dd[j][np]; val[j]-=tem;
      if(val[j]<=lb[j] || val[j]>=ub[j]) outb=1;
      if(fabs(tem)>mxdif) mxdif=fabs(tem);
      /* if(out==DBL_MAX) outb=1;  // Solaris only has MAXDOUBLE
         Linux has both */
      if(out==MAXDOUBLE) outb=1;
    }
    /* modified NR to limit size of step to within bounds */

    while(outb==1 || mxdif>bdd)
    { mxdif/=2.;
      for(j=0,outb=0;j<np;j++) 
      { dd[j][np]/=2; tem=dd[j][np]; val[j]+=tem; 
        if(val[j]<=lb[j] || val[j]>=ub[j]) { outb=1; }
      }
    }

    iter++;
#ifndef R
    if(iprint>=1)
    { printf("iter. %d :", iter);
      printf(" fnval = %f; ", out);
      printf("mxdif =%f : ",  mxdif);
      for(j=0;j<np;j++) printf("%f ", val[j]); printf("\n");
    }
#else
    if(iprint>=1)
    { Rprintf("iter. %d :", iter);
      Rprintf(" fnval = %f; ", out);
      Rprintf("mxdif =%f : ",  mxdif);
      for(j=0;j<np;j++) Rprintf("%f ", val[j]); Rprintf("\n");
    }
#endif
  }

  while(mxdif>eps && iter<mxiter);

  if(iter>=mxiter)
  {
#ifndef R
    printf("Did not converge\n");
#else
    Rprintf("Did not converge\n");
#endif
    *iconv=0;
  }
  else { *iconv=1; }
#ifndef R
  free(dd[0]); free(dd);
#else
  Free(dd[0]); Free(dd);
#endif
  *fnval=out;

  /* inverse Hessian */
  dd=dmatrix(np,2*np);
  func(np,val,&out,dd,iprint,
    n,m,nq,nrec,x,w,dat,fr,g,g1,g2,pp);
  for(j=0;j<np;j++)
  { for(k=np;k<2*np;k++) dd[j][k]=0;
    dd[j][j+np]=1;
#ifndef R
    if(iprint>=2)
    { for(k=0;k<np;k++) printf("%f ", dd[j][k]); printf("\n"); }
/*#else
    if(iprint==1)
    { for(k=0;k<np;k++) Rprintf("%f ", dd[j][k]); Rprintf("\n"); }*/
#endif
  }
  lsolve(dd,np,2*np,&det,1.e-6);
  /* inverse now in right partition of matrix */
  for(j=0;j<np;j++)
  { for(k=0;k<np;k++) invhes[j][k]=dd[j][k+np]; }
#ifndef R
  free(dd[0]); free(dd);
#else
  Free(dd[0]); Free(dd);
#endif
}

#ifndef ALL
void lsolve(double **A, int n, int nk, double *det, double tol)
/* Routine to solve the linear system  a*z=b */
/* Inputs:
   A = matrix of the linear system (A=[a|b])
        A must be allocated as an (n*nk) matrix 
       (the 0th row and 0th column are used here)
   n = order of the matrix 'a'
   nk = number of columns of 'A' to use, nk-n is the number of
        right hand side vectors, or the column dimension of 'b'
   tol= tolerance for the diagonal elements in partial pivoting

   Outputs:
   det = determinant
   A in columns n,...,nk-1 : the solution vectors z corresponding to
        different right hand side vectors in b
   A in columns 0,...,n-1:  these will be changed to a "reduced" matrix
        which is triangular.
*/

{ int parity,i,j,k,l,mrow;
  double mx,t,sum;
  parity=1;
  
/*  loop for Gaussian elimination with partial pivoting */
  for(k=0;k<n-1;k++)
  { mx=fabs(A[k][k]); mrow=k;
    for(l=k+1;l<n;l++)
    { if(fabs(A[l][k])>mx) { mx=fabs(A[l][k]); mrow=l;} }
    if(mx<=tol) { *det=0.; return;}
    if(mrow>k)
    { parity*= -1;
      for(i=0;i<nk;i++)
      { t=A[mrow][i]; A[mrow][i]=A[k][i]; A[k][i]=t;}
    }
    for(i=k+1;i<n;i++)
    { A[i][k]/=A[k][k];
      for(j=k+1;j<nk;j++) A[i][j]-=A[i][k]*A[k][j];
    }
  }
  if(fabs(A[n-1][n-1])<=tol) { *det=0.; return;}
  for(i=0,*det=parity;i<n;i++) *det*=A[i][i];
  if(n==nk) return;
  
/*  loop for backsolving to find solution for right hand sides */
  for(l=n;l<nk;l++)
  { A[n-1][l]/=A[n-1][n-1];
    for(i=n-2;i>=0;i--)
    { for(j=i+1,sum=0.;j<n;j++) sum+=A[i][j]*A[j][l];
      A[i][l]=(A[i][l]-sum)/A[i][i];
    }
  }
}
#endif

#ifdef ALONE
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
#endif

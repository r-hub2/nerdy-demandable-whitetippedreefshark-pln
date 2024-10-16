#include <stdio.h>
#include <math.h>
//#define NP 101
//void geppldet(double A[][2*NP], int n, int nk, double *ldet, double tol, 
void geppldet(double **A, int n, int nk, double *ldet, double tol, 
  int *par)
/* for large matrices, return log(det) */
/* routine for solving a linear system of equations */
{ int parity,i,j,k,l,mrow;
  double mx,t,sum;
  parity=1;
  for(k=1;k<n;k++)
  { mx=fabs(A[k][k]); mrow=k;
    for(l=k+1;l<=n;l++)
    { if(fabs(A[l][k])>mx) { mx=fabs(A[l][k]); mrow=l; } }
    if(mx<=tol) { *ldet=-1.; *par=0; return; }
    if(mrow>k)
    { parity*= -1;
      for(i=1;i<=nk;i++)
      { t=A[mrow][i]; A[mrow][i]=A[k][i]; A[k][i]=t; }
    }
    for(i=k+1;i<=n;i++)
    { A[i][k]/=A[k][k];
      for(j=k+1;j<=nk;j++) A[i][j]-=A[i][k]*A[k][j];
    }
  }
  if(fabs(A[n][n])<=tol) { *ldet=-1.; *par=0; return; }
  for(i=1,*ldet=0,*par=parity;i<=n;i++) 
  { *ldet+=log(fabs(A[i][i]));
    if(A[i][i]<0) *par= -(*par);
  }
  if(n==nk) return;
  for(l=n+1;l<=nk;l++)
  { A[n][l]/=A[n][n];
    for(i=n-1;i>=1;i--)
    { for(j=i+1,sum=0.;j<=n;j++) sum+=A[i][j]*A[j][l];
      A[i][l]=(A[i][l]-sum)/A[i][i];
    }
  }
}

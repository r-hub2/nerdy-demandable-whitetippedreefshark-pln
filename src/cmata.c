#ifndef R
#include <malloc.h>
#else
#include <R.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* input Xi, Delta, output Cmat=Deltac*(Deltac'*Xi*Deltac)^{-1} *Deltac'
 where Deltac is the orthogonal complement of Delta
  dimension of Xi is sxs, dimension of Delta is sxq (s>q),
  dimension of Deltac is sx(s-q) and Deltac'*Delta=0
*/
/*#define NM 1500 */
/*#define NP 51 */
#ifdef MAIN
/* this main program has not been converted to pointers so won't work now */
main(int argc, char *argv[])
{ int s,q,i,j;
  double xi[NM][NM],del[NM][NP],cmat[NM][NM];
  void qfmat(int, int, double [][NM], double[][NP], double [][NM]);
  
  scanf("%d %d", &s,&q);
  for(i=1;i<=s;i++)
  { for(j=1;j<=s;j++) scanf("%lf", &xi[i][j]); }
  for(i=1;i<=s;i++)
  { for(j=1;j<=q;j++) scanf("%lf", &del[i][j]); }
  qfmat(s,q,xi,del,cmat);
  /*cmat(ss,np,Xir,Delr,Sinv);*/
  for(i=1;i<=s;i++)
  { for(j=1;j<=s;j++) printf("%.5f ", cmat[i][j]); printf("\n"); }
  printf("===============================================\n");
  exit(0);
}
#endif

void qfmat(int ss, int q, double **Xir, double **Delr, double **Cmat)
/*void qfmat(int ss,int q,double Xir[][NM],double Delr[][NP], double Cmat[][NM])*/
{ int ipr,i,j,k,ssq,parity;
  /*void nullsp(double [][NP], int, int, double [][NM], int);*/
  void nullsp(double **, int, int, double **, int);
  /*double Delc[NM][NM],A[NM][2*NM],b[NM][NM],sum,ldet;*/
  double **Delc,**A,**b,sum,ldet;
  /*void gepp3(double A[][2*NM], int, int nk, double *, double, int *);*/
  void gepp3(double **A, int, int nk, double *, double, int *);
  double **dmatrix(int,int);

  ssq=ss-q;
  Delc=dmatrix(ss+1,ssq+1);
  ipr=0;
  nullsp(Delr, ss, q, Delc, ipr);
  A=dmatrix(ssq+1,ssq+ss+1);
  b=dmatrix(ss+1,ssq+1);
  /*printf("dimension of Delc should be %d x %d\n", ss,ss-q);*/
  /*for(i=1;i<=ss;i++)*/
  /*{ for(j=1;j<=ssq;j++) printf("%.4f ", Delc[i][j]);*/
  /*  printf("\n");*/
  /*}*/

  /* Cmat = Delc (Delc' * Xir * Delc)^{-1} Delc' */
  for(i=1;i<=ss;i++)
  { for(j=1;j<=ssq;j++)
    { for(k=1,sum=0.;k<=ss;k++) sum+=Xir[i][k]*Delc[k][j];
      b[i][j]=sum;
    }
  }
  for(i=1;i<=ssq;i++)
  { for(j=1;j<=ssq;j++)
    { for(k=1,sum=0.;k<=ss;k++) sum+=Delc[k][i]*b[k][j];
      A[i][j]=sum;
    }
    for(j=1;j<=ss;j++) A[i][ssq+j]=Delc[j][i];
  }
#ifndef R
  free(b[0]); free(b); 
#else
  Free(b[0]); Free(b); 
#endif
  gepp3(A,ssq,ssq+ss,&ldet,1.e-12,&parity);
#ifndef R
  if(parity==0) { printf("zero determinant for A\n"); }
#endif
  /* soln in columns ssq+1 ... ssq+ss */
  for(i=1;i<=ss;i++)
  { for(j=1;j<=ss;j++)
    { for(k=1,sum=0.;k<=ssq;k++) sum+=Delc[i][k]*A[k][ssq+j];
      Cmat[i][j]=sum;
    }
  }

  /* print Cmat */
  /*for(i=1;i<=ss;i++)*/
  /*{ for(j=1;j<=ss;j++) printf("%.5f ", Cmat[i][j]); printf("\n"); }*/
#ifndef R  
  free(A[0]); free(A); free(Delc[0]); free(Delc);
#else
  Free(A[0]); Free(A); Free(Delc[0]); Free(Delc);
#endif
}

void gepp3(double **A, int n, int nk, double *ldet, double tol, int *par)
/*void gepp3(double A[][2*NM], int n, int nk, double *ldet, double tol, int *par)*/
/* for large matrices, return log(det) */
/* return par=0 for zero determinant (singular) */
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


/* check on getting a basis for the null space of an m*n matrix
   or rank r, r<=m<=n */
/* null space or orthogonal complement */
void nullsp(double **a, int n, int m, double **acomp, int ipr)
/*void nullsp(double a[][NP], int n, int m, double acomp[][NM], int ipr)*/
{ /*double a0[NP][NM],v[NM][NM],z[NM];*/
  double **a0,**v,*z;
  double eps,azero;
  int i,j,r,k;
  /*void a1svd2(int ,int , double [][NM], double, double [][NM], double z[], int);*/
  void a1svd2(int ,int , double **, double, double **, double *z, int);
  double **dmatrix(int,int);

  eps=1.e-10;
  /*ipr=1;*/
  azero=1.e-10;
  a0=dmatrix(m+1,n+1);
  for(i=1;i<=m;i++) { for(j=1;j<=n;j++) a0[i][j]=a[j][i]; }

  /* need original copy of a which gets changes in a1svd*/
  v=dmatrix(n+1,n+1);
  #ifndef R
  z=(double *) malloc((n+1) * sizeof(double));
  #else
  z=(double *) Calloc((n+1), double);
  #endif
  a1svd2(m,n,a0,eps,v,z,ipr);

#ifndef R
  free(a0[0]); free(a0);
#else
  Free(a0[0]); Free(a0);    
#endif
  /* extract last columns of v */
  for(i=1;i<=m;i++) { if(z[i]<azero) break; } 
  r=i-1;
#ifndef R
  if(r!=m) printf("nullsp: r=%d != m=%d\n", r,m);
#endif
  /*r=m; // this should be the case if full column rank*/
  /* check that a*v2=0 where v2 is the last n-r columns of v */
  /*printf("dimension of V2=%d x %d\n", n,n-r);*/

  /* return orthogonal complement of a[][] in acomp[][] */
  for(k=1;k<=n;k++) 
  { for(j=1;j<=n-r;j++) acomp[k][j]=v[k][j+r]; }
#ifndef R
  free(v[0]); free(v); free(z);
#else
  Free(v[0]); Free(v); Free(z);
#endif
}


/* this routine does not use the index of 0 for rows and
   columns, and matrices are not passed as pointers
*/
/*void a1svd2(int m, int n, double a[][NM], double eps, double v[][NM],*/
void a1svd2(int m, int n, double **a, double eps, double **v,
  double *z, int ipr)
/* algorithm 1  singular value decomposition by column orthogona- 
   lisation via plane rotations */
/* Nash, J.C. (1990). Compact Numerical Methods for Computers: Linear
   Algebra and Function Minimisation, second edition. Hilger, New York.*/
/* J.C. Nash   July 1978, February 1980, April 1989 */
/* translated from Fortran to C by H. Joe */
/* m by n  matrix a  is decomposed to  u*z*v^t 
   a    =   array containing a [input],  u [output] 
   na   =   first dimension of a (not needed for C) 
   eps  =   machine precision 
   v    =   array in which orthogonal matrix v is accumulated 
   nv   =   first dimension of v  (not needed for C)
   z    =   vector of singular values 
   ipr   =  if ipr>0  then print some extra output
*/

{ int i,j,k,j1,n1,count;
  double tol,p,q,r,vv,c,s;
  /* underflow avoidance strategy */
  double small;
  small=1.0e-36;
  /* above is value for IBM */
  tol=n*n*eps*eps;
  for(i=1;i<=n;i++)
  { for(j=1;j<=n;j++) v[i][j]=0.0;
    v[i][i]=1.0;
  }
  n1=n-1;
  count=n*(n-1)/2;
  while(count>0)
  { for(j=1;j<=n1;j++)
    { j1=j+1;
      for(k=j1;k<=n;k++)
      { p=0.0; q=0.0; r=0.0;
        for(i=1;i<=m;i++)
        { if(fabs(a[i][j])>small&&fabs(a[i][k])>small) p+=a[i][j]*a[i][k];
          if(fabs(a[i][j])>small) q+=a[i][j]*a[i][j];
          if(fabs(a[i][k])>small) r+=a[i][k]*a[i][k];
        }
        if(q<r) { c=0.0; s=1.0;}
        else
        { if(r<=tol)
          { count--; continue;}
          if((p*p)/(q*r)<tol) { count--; continue;}
          q-=r; vv=sqrt(4.0*p*p+q*q);
          c=sqrt((vv+q)/(2.0*vv)); s=p/(vv*c);
        }
        for(i=1;i<=m;i++)
        { r=a[i][j]; a[i][j]=r*c+a[i][k]*s;
          a[i][k]=-r*s+a[i][k]*c;
        }
        for(i=1;i<=n;i++)
        { r=v[i][j]; v[i][j]=r*c+v[i][k]*s;
          v[i][k]=-r*s+v[i][k]*c;
        }
      }
    }
#ifndef R    
    if(ipr>0) printf("%4d rotations\n", count);
#endif
  }
  
  for(j=1;j<=n;j++)
  { q=0.0;
    for(i=1;i<=m;i++) q+=a[i][j]*a[i][j];
    q=sqrt(q);
    z[j]=q;
#ifndef R
    if(ipr>0) printf("sv(%3d)= %16.8e\n", j,q);
#endif
    if(q<tol) continue;
    for(i=1;i<=m;i++) a[i][j]/=q;
  }
}


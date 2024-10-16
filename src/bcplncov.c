#ifndef R
#include <malloc.h>
#else
#include <R.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// M2 = (M-1)^2+1
//#define M2 26
//  JK = N choose 2
//#define JK 200
/* computation of Godambe matrix for bivariate composite likelihood method,
   polytomous logit-normit model, */
// parameters a11 ... a_{1,m-1} a21 ... a_{n,m-1} b1 ... bn
#define SQRTPI 1.772453850905516
#ifdef MAINC
main(int argc, char *argv[])
{ int i,m,j,n,nobs,ip,jp,np,iprint;
  //double V[NP][NP],param[NP],alp[N][M],b[N];
  double **V,*param,**alp,*b;
  void bclcov(int *, int *, int *, double *param, double **V, int *);
  double **dmatrix(int,int);

  setbuf(stdout,NULL);
  scanf("%d %d", &n,&m); // number of items and categories
  nobs=1;
  np=m*n; 
  //V=dmatrix(np+1,np+1);
  V=dmatrix(np,np);

  param=(double *) malloc((np+1) * sizeof(double));
  alp=dmatrix(n+1,m);
  b=(double *) malloc((n+1) * sizeof(double));

  // alp[i][1],...,alp[i][m-1] should be in decreasing order
  ip=0; // change to ip++ after param 
  for(i=1;i<=n;i++)
  { for(j=1;j<m;j++) 
    { scanf("%lf", &alp[i][j]); param[ip]=alp[i][j]; ip++; }
  }
  for(i=1;i<=n;i++) 
  { scanf("%lf", &b[i]); param[ip]=b[i]; ip++; }
  scanf("%d", &nobs);  // sample size
  printf("#items=%d, #categ=%d, sample size=%d\n", n,m,nobs);
  printf("Parameter values alp[i][j], j=1..m-1, i=1..n; b[i]\n");
  for(i=1;i<=n;i++)
  { for(j=1;j<m;j++) printf("%f ", alp[i][j]);  printf(" : "); }
  printf("\n");
  for(i=1;i<=n;i++) printf("%f ", b[i]); printf("\n");
  iprint=1;
  //bclcov(&n,&m,&nobs,alp,b,V,&iprint);
  bclcov(&n,&m,&nobs,param,V,&iprint);
  // print V and SE from sqrt(diag(V))
  printf("\nmain : estimated covariance matrix of BCL estimate\n");
  for(ip=1;ip<=np;ip++)
  { for(jp=1;jp<=np;jp++) printf("%f ", V[ip-1][jp-1]);  printf("\n"); }
  printf("\nBCLestimates and SEs:\n");
  for(ip=1;ip<=np;ip++)
  { printf("%f %f\n", param[ip-1], sqrt(V[ip-1][ip-1]) ); }
  printf("===============================================\n");
  
  free(param);
  free(V[0]); free(V);
  free(b); free(alp[0]); free(alp);
  exit(0);
}
#endif

/* estimating equation approach
   solve psi(theta,data)=0   where psi involves a sum of length N=sample size,
   and theta is a parameter vector of length q.
   Here psi is the vector of derivatives of the bivariate composite likelihood 
   and for the polytomous logit-normit model, q=m*n,
   where n=#items, m=#categories.
   The D=H hessian matrix comes from N^{-1} {partial psi/ partial theta'}
   so  D is qxq. [H is now the common notation in the composite likelihood
   literature, e.g., Varin, Reid & Firth (2011), Statistica Sinica, 21, 5-42.]
   The M=J matrix comes from N^{-1} Var(psi) or the covariance matrix
   of one term in the sum psi. [J is now the common notation for this
   matrix in the composite likelihood literature.]
   The V matrix for the covariance matrix of the BCL estimate is
     D^{-1} M (D')^{-1}
*/
/* n0=#items
   m=#categories
   param: vector of parameters which is converted to 
       alp= matrix of alphas (cutpoints)
       b = vector of slopes
   V = covariance matrix of BCL estimator, returned by this function
   iprint (set =1 to get print of intermediate calculations in Unix)
*/
#ifndef R
void bclcov(int *n00, int *m0, int *nobs0, double *param, double **V,  
   int *iprint)
{
  int i,j,rr,np,n2;
  int jk,n20,j1p,j2p,jk2,parity;
  double sum,tem; 
  void invjk(int jk, int n, int *, int *, int *cb);
  //double derjk[JK][M2][NP],D[NP][NP],MM[NP][NP],A[NP][2*NP],ldet,pr,der[NP];
  double ***derjk,**D,**MM,**A,ldet,pr,*der;
  //void geppldet(double [][2*NP], int n, int nk, double *ldet, double tol,
  void geppldet(double **,int n,int nk, double *ldet, double tol, int *parity);
  void plgnderghi(int n, int m, double **alp, double *b, int r, 
   int *ii, int *kk, double *pr, double *der, int ideriv,
   double *x, double *w, int nq);
  void d2v(int d, int k, int ii, int *jj);
  int ii[5],kk[5],*cb;
  int ip,jp,j1,j2,k1,k2,iv,r,ic,m2,ijk,ic1,ic2;
  double **dmatrix(int,int);
  double ***gmatrix(int,int,int);
  void gauher(double *, double *, int);
  int nq,n0,m,nobs;
  double *x,*w;
  double **alp, *b;

  n0= *n00; m=*m0; nobs=*nobs0;
  nq=48;
  printf("\nnq=%d\n",nq);

  x=(double *) malloc((nq+1) * sizeof(double));
  w=(double *) malloc((nq+1) * sizeof(double));

  gauher(x,w,nq);
  for (j=1;j<=nq;j++) x[j]*=M_SQRT2;
  for (j=1;j<=nq;j++) w[j]/=SQRTPI;
  //for (j=1;j<=nq;j++) printf("%3d %14.6e %14.6e\n",j,x[j],w[j]);

  for(i=1,n20=1;i<=n0;i++) n20*=m;
  if(*iprint==1) printf("\n#items=%d, #cells=%d^%d=%d\n", n0,m,n0,n20);
  np=m*n0; m2=m*m;
  A=dmatrix(np+1,2*np+1);
  MM=dmatrix(np+1,np+1);
  D=dmatrix(np+1,np+1);
  derjk=gmatrix(n0*(n0-1)/2+1,(m*m)+1,np+1);
  der=(double *) malloc((np+1) * sizeof(double));
  cb=(int *) malloc((n0+2) * sizeof(int));

  for(i=1;i<=n0;i++) cb[i]= (i*(i-1))/2;
  ip=0;
  alp=dmatrix(n0+1,m);
  b=(double *) malloc((n0+1) * sizeof(double));

  ip=0; // ip++ after param later
  for(i=1;i<=n0;i++)
  { for(j=1;j<m;j++) 
    { alp[i][j]=param[ip]; ip++; }
  }
  for(i=1;i<=n0;i++) { b[i]=param[ip]; ip++; }

  // initialize the D matrix
  for(ip=1;ip<=np;ip++)
  { for(jp=1;jp<=np;jp++) D[ip][jp]=0.; }
  // step 1 bivariate quantities
  // step 2, the D matrix plus store terms used in M matrix
  ijk=0;
  for(j2=2;j2<=n0;j2++)
  { ii[2]=j2;
    for(j1=1;j1<j2;j1++)
    { ii[1]=j1; ic=0; ijk++;
      for(k1=0;k1<m;k1++) 
      { kk[1]=k1; 
        for(k2=0;k2<m;k2++)
        { kk[2]=k2;
          plgnderghi(n0, m, alp, b, 2, ii, kk, &pr, der,1, x,w,nq);
          // derjk is needed for computation of M matrix
          for(ip=1;ip<=np;ip++) derjk[ijk][ic][ip]=0.;
          if(pr>0.) 
          { for(ip=1;ip<=np;ip++)
            { for(jp=1;jp<=np;jp++) D[ip][jp]+=der[ip]*der[jp]/pr; }
            //for(ip=1;ip<=np;ip++) printf("%f ", der[ip]);
            //printf(": %f\n", pr);
            for(ip=1;ip<=np;ip++) derjk[ijk][ic][ip]=der[ip]/pr;
          }
          ic++;
        }
      }
    }
  }

  if(*iprint==1)
  { printf("D matrix (symmetric)\n");
    for(ip=1;ip<=np;ip++)
    { for(jp=1;jp<=np;jp++) printf("%f ", D[ip][jp]);  printf("\n"); }
  }
    
  // step 3, the M matrix, 
  for(ip=1;ip<=np;ip++)
  { for(jp=1;jp<=np;jp++) MM[ip][jp]=0.; }
  // handle cases of 3 and 4 distinct indices in j1,j2,j1',j2'
  // printf("now executing a loop of %dx%d\n", cb[n0],cb[n0]); 
  for(jk=1;jk<cb[n0];jk++)
  { invjk(jk,n0,&j1,&j2,cb);
    for(jk2=jk+1;jk2<=cb[n0];jk2++)
    { invjk(jk2,n0,&j1p,&j2p,cb);
      // determine number of distinct elements among j1,j2,j1p,j2p
      if(j1==j1p || j2==j2p || j2==j1p) 
      { //printf(": 3 distinct\n");
        r=3; 
        n2=m*m*m;
        if(j1==j1p)
        { ii[1]=j1; ii[2]=j2; ii[3]=j2p;
          // sum over cases of kk[]
          for(iv=0;iv<n2;iv++)
          { d2v(r,m,iv,kk);
            // actually don't need ders here so use ideriv=0
            plgnderghi(n0, m, alp, b, r, ii, kk, &pr, der, 0, x,w,nq);
            ic1=kk[1]*m+kk[2];
            ic2=kk[1]*m+kk[3];
            for(ip=1;ip<=np;ip++)
            { for(jp=1;jp<=np;jp++)
              { tem=pr*derjk[jk][ic1][ip]*derjk[jk2][ic2][jp]; 
                MM[ip][jp]+=tem; MM[jp][ip]+=tem;
              }
            }
          }
        }
        else if(j2==j2p)
        { ii[1]=j1; ii[2]=j1p; ii[3]=j2p;
          for(iv=0;iv<n2;iv++)
          { d2v(r,m,iv,kk);
            plgnderghi(n0, m, alp, b, r, ii, kk, &pr, der, 0, x,w,nq);
            ic1=kk[1]*m+kk[3];
            ic2=kk[2]*m+kk[3];
            for(ip=1;ip<=np;ip++)
            { for(jp=1;jp<=np;jp++)
              { tem=pr*derjk[jk][ic1][ip]*derjk[jk2][ic2][jp]; 
                MM[ip][jp]+=tem; MM[jp][ip]+=tem;
              }
            }
          }
        }
        else /* j2==j1p */
        { ii[1]=j1; ii[2]=j1p; ii[3]=j2p;
          for(iv=0;iv<n2;iv++)
          { d2v(r,m,iv,kk);
            plgnderghi(n0, m, alp, b, r, ii, kk, &pr, der, 0, x,w,nq);
            ic1=kk[1]*m+kk[2];
            ic2=kk[2]*m+kk[3];
            for(ip=1;ip<=np;ip++)
            { for(jp=1;jp<=np;jp++)
              { tem=pr*derjk[jk][ic1][ip]*derjk[jk2][ic2][jp]; 
                MM[ip][jp]+=tem; MM[jp][ip]+=tem;
              }
            }
          }
        }
      }

      else 
      { //printf(": 4 distinct\n");
        r=4; n2=m*m; n2=n2*n2;
        ii[1]=j1; ii[2]=j2; ii[3]=j1p; ii[4]=j2p;
        for(iv=0;iv<n2;iv++)
        { d2v(r,m,iv,kk);
          plgnderghi(n0, m, alp, b, r, ii, kk, &pr, der, 0, x,w,nq);
          ic1=kk[1]*m+kk[2];
          ic2=kk[3]*m+kk[4];
          for(ip=1;ip<=np;ip++)
          { for(jp=1;jp<=np;jp++)
            { tem=pr*derjk[jk][ic1][ip]*derjk[jk2][ic2][jp]; 
              MM[ip][jp]+=tem; MM[jp][ip]+=tem;
            }
          }
        }
      }
    }
  }
  if(*iprint==1)
  { printf("M matrix\n");
    for(ip=1;ip<=np;ip++)
    { for(jp=1;jp<=np;jp++) printf("%f ", MM[ip][jp]);  printf("\n"); }
  }
    
  // add D matrix to MM
  for(ip=1;ip<=np;ip++)
  { for(jp=1;jp<=np;jp++) MM[ip][jp]+=D[ip][jp]; }
  printf("M matrix\n");
  for(ip=1;ip<=np;ip++)
  { for(jp=1;jp<=np;jp++) printf("%f ", MM[ip][jp]);  printf("\n"); }
    
  // invert D matrix, then D^{-1}*M*D^{-1} (since D is symmetric)
  for(jp=1;jp<=np;jp++)
  { for(ip=1;ip<=np;ip++)
    { A[jp][ip]=D[jp][ip]; A[jp][ip+np]=0.; }
    A[jp][jp+np]=1.;
  }
  geppldet(A,np,2*np,&ldet,1.e-12,&parity);
  if(parity==0 && *iprint==1) { printf("zero determinant for D\n"); }
  //printf("D inverse \n");
  //for(ip=1;ip<=np;ip++)
  //{ for(jp=1;jp<=np;jp++) printf("%f ", A[ip][jp+np]);  printf("\n"); }

  for(jp=1;jp<=np;jp++)
  { for(ip=1;ip<=np;ip++)
    { for(rr=1,sum=0.;rr<=np;rr++) sum+=A[jp][rr+np]*MM[rr][ip];
      A[jp][ip]=sum;
    }
  }
  for(jp=1;jp<=np;jp++)
  { for(ip=1;ip<=np;ip++)
    { for(rr=1,sum=0.;rr<=np;rr++) sum+=A[jp][rr]*A[rr][ip+np];
      V[jp-1][ip-1]=sum;
    }
  }
  if(*iprint==1)
  { printf("V matrix, asymptotic version scaled by sample size \n");
    for(ip=1;ip<=np;ip++)
    { for(jp=1;jp<=np;jp++) printf("%f ", V[ip-1][jp-1]);  printf("\n"); }
  }
  /* for sample of size nobs, V <- V/nobs to get estimated covariance matrix
     of BCL estimator */
  for(ip=1;ip<=np;ip++)
  { for(jp=1;jp<=np;jp++) V[ip-1][jp-1]/=nobs; }
  free(A[0]); free(A);
  free(MM[0]); free(MM);
  free(D[0]); free(D);
  free(derjk[0][0]); free(derjk[0]); free(derjk);
  free(der); free(cb);
  free(w); free(x); 
  free(b); free(alp[0]); free(alp);

}
#endif

/* invert from jk to (j,k) */
void invjk(int jk, int n, int *j, int *k, int *cb)
{ int i;
  for(i=1;i<=n;i++) 
  { if(jk<=cb[i]) break; }
  *k=i; *j=jk-cb[i-1]; 
}

/* polytomous logit-normit probability for margin with derivative */
// ideriv =1 if derivatives needed, otherwise ideriv=0
void plgnderghi(int n, int m, double **alp, double *b, int r, 
   int *ii, int *kk, double *pr, double *der, int ideriv,
   double *x, double *w, int nq)
{ int j,np,id,nd;
  void pghderi(double, double **alp, double *b, 
     int n, int m, int r, int *ii, int *kk, double *, int);
  //double fn[NP];
  double *fn;

  np=m*n;
#ifndef R
  fn=(double *) malloc((np+1) * sizeof(double));
#else
  fn=(double *) Calloc((np+1),double);
#endif
  if(ideriv==1) nd=np; else nd=0;
  // loop for probabilities and partial derivatives
  for(id=0;id<=nd;id++) der[id]=0.;
  for(j=1;j<=nq;j++)
  { pghderi(x[j],alp,b,n,m,r,ii,kk,fn,ideriv);
    // return a vector of fn[], [0] element is the function
    for(id=0;id<=nd;id++) der[id]+=w[j]*fn[id];
  }
  *pr=der[0];
#ifndef R
  free(fn);
#else
  Free(fn);
#endif
}

/* for faster speed, avoid the derivatives for M matrix, 
   speeds up a little for n=20, #categs=4 */ 
// ii[] is subset of 1...n, does not have to be ordered
// ideriv =1 if derivatives needed, otherwise ideriv=0
void pghderi(double xx, double **alp, double *b, int n, int m, int r, 
   int *ii, int *kk, double *fn, int ideriv)
{ double lg1,lg2,diff,tem,da1,da2,db,dtem;
  int i,j,k,np,id,np1,nd;

  np=m*n; np1=np-n;
  if(ideriv==1) nd=np; else nd=0;
  fn[0]=1.;
  for(id=1;id<=np;id++) fn[id]=0.;
  for(j=1,tem=1.;j<=r;j++)
  { i=ii[j]; k=kk[j];
    if(k==0) 
    { lg2=1.; lg1=1./(1.+exp(-alp[i][1]-b[i]*xx)); 
      da2=0.; da1=lg1*(1.-lg1);
      diff=lg2-lg1;
    // skip next lines for M matrix
      if(ideriv>0)
      { if(diff<=0.) dtem=0.; else dtem=-da1/diff;
        fn[(i-1)*(m-1)+1]= dtem;
      }
    }
    else if(k==m-1) 
    { lg2=1./(1.+exp(-alp[i][m-1]-b[i]*xx)); lg1=0.; 
      da1=0.; da2=lg2*(1.-lg2);
      diff=lg2-lg1;
    // skip next lines for M matrix
      if(ideriv>0)
      { if(diff<=0.) dtem=0.; else dtem=da2/diff;
        fn[(i-1)*(m-1)+m-1]= dtem;
      }
    }
    else 
    { lg2=1./(1+exp(-alp[i][k]-b[i]*xx)); lg1=1./(1+exp(-alp[i][k+1]-b[i]*xx)); 
      da1=lg1*(1.-lg1); da2=lg2*(1.-lg2);
      diff=lg2-lg1;
    // skip next lines for M matrix
      if(ideriv>0)
      { if(diff<=0.) dtem=0.; else dtem=-da1/diff;
        fn[(i-1)*(m-1)+k+1]= dtem;
        if(diff<=0.) dtem=0.; else dtem=da2/diff;
        fn[(i-1)*(m-1)+k]= dtem;
      }
    }
    // skip next lines for M matrix
    if(ideriv>0)
    { db=xx*(da2-da1);
      if(diff<=0.) dtem=0.; else dtem=db/diff;
      fn[np1+i]=dtem;
    }
    tem*=diff;
  }
  for(id=0;id<=nd;id++) fn[id]*=tem;
}

#ifdef R
void Rbclcov(int *n00, int *m0, int *nobs0, double *param, double *Vout, int *nq,  
   int *iprint)
{
  double **V;
  int i,j,rr,np,n2;
  int jk,n20,j1p,j2p,jk2,parity;
  double sum,tem; 
  void invjk(int jk, int n, int *, int *, int *cb);
  //double derjk[JK][M2][NP],D[NP][NP],MM[NP][NP],A[NP][2*NP],ldet,pr,der[NP];
  double ***derjk,**D,**MM,**A,ldet,pr,*der;
  //void geppldet(double [][2*NP], int n, int nk, double *ldet, double tol,
  void geppldet(double **,int n,int nk, double *ldet, double tol, int *parity);
  void plgnderghi(int n, int m, double **alp, double *b, int r, 
   int *ii, int *kk, double *pr, double *der, int ideriv,
   double *x, double *w, int nq);
  void d2v(int d, int k, int ii, int *jj);
  int ii[5],kk[5],*cb;
  int ip,jp,j1,j2,k1,k2,iv,r,ic,m2,ijk,ic1,ic2;
  double **dmatrix(int,int);
  double ***gmatrix(int,int,int);
  void gauher(double *, double *, int);
  int n0,m,nobs;
  double *x,*w;
  double **alp, *b;

  n0= *n00; m=*m0; nobs=*nobs0;
  //if(*iprint==1) Rprintf("\nnq=%d\n",*nq);
  x=(double *) Calloc((*nq+1), double);
  w=(double *) Calloc((*nq+1), double);
  gauher(x,w,*nq);
  for (j=1;j<=*nq;j++) x[j]*=M_SQRT2;
  for (j=1;j<=*nq;j++) w[j]/=SQRTPI;
  //for (j=1;j<=nq;j++) printf("%3d %14.6e %14.6e\n",j,x[j],w[j]);

  for(i=1,n20=1;i<=n0;i++) n20*=m;
  //if(*iprint==1) Rprintf("\n#items=%d, #cells=%d^%d=%d\n", n0,m,n0,n20);
  np=m*n0; m2=m*m;
  V=dmatrix(np,np);  
  A=dmatrix(np+1,2*np+1);
  MM=dmatrix(np+1,np+1);
  D=dmatrix(np+1,np+1);
  derjk=gmatrix(n0*(n0-1)/2+1,(m*m)+1,np+1);
  der=(double *) Calloc((np+1), double);
  cb=(int *) Calloc((n0+2), int);
  for(i=1;i<=n0;i++) cb[i]= (i*(i-1))/2;
  ip=0;
  alp=dmatrix(n0+1,m);
  b=(double *) Calloc((n0+1), double);
  ip=0; // ip++ after param later
  for(i=1;i<=n0;i++)
  { for(j=1;j<m;j++) 
    { alp[i][j]=param[ip]; ip++; }
  }
  for(i=1;i<=n0;i++) { b[i]=param[ip]; ip++; }

  // initialize the D matrix
  for(ip=1;ip<=np;ip++)
  { for(jp=1;jp<=np;jp++) D[ip][jp]=0.; }
  // step 1 bivariate quantities
  // step 2, the D matrix plus store terms used in M matrix
  ijk=0;
  for(j2=2;j2<=n0;j2++)
  { ii[2]=j2;
    for(j1=1;j1<j2;j1++)
    { ii[1]=j1; ic=0; ijk++;
      for(k1=0;k1<m;k1++) 
      { kk[1]=k1; 
        for(k2=0;k2<m;k2++)
        { kk[2]=k2;
          plgnderghi(n0, m, alp, b, 2, ii, kk, &pr, der,1, x,w,*nq);
          // derjk is needed for computation of M matrix
          for(ip=1;ip<=np;ip++) derjk[ijk][ic][ip]=0.;
          if(pr>0.) 
          { for(ip=1;ip<=np;ip++)
            { for(jp=1;jp<=np;jp++) D[ip][jp]+=der[ip]*der[jp]/pr; }
            //for(ip=1;ip<=np;ip++) printf("%f ", der[ip]);
            //printf(": %f\n", pr);
            for(ip=1;ip<=np;ip++) derjk[ijk][ic][ip]=der[ip]/pr;
          }
          ic++;
        }
      }
    }
  }

/*  if(*iprint==1)
  { Rprintf("D matrix (symmetric)\n");
    for(ip=1;ip<=np;ip++)
    { for(jp=1;jp<=np;jp++) Rprintf("%f ", D[ip][jp]);  Rprintf("\n"); }
  }*/

  // step 3, the M matrix, 
  for(ip=1;ip<=np;ip++)
  { for(jp=1;jp<=np;jp++) MM[ip][jp]=0.; }
  // handle cases of 3 and 4 distinct indices in j1,j2,j1',j2'
  // printf("now executing a loop of %dx%d\n", cb[n0],cb[n0]); 
  for(jk=1;jk<cb[n0];jk++)
  { invjk(jk,n0,&j1,&j2,cb);
    for(jk2=jk+1;jk2<=cb[n0];jk2++)
    { invjk(jk2,n0,&j1p,&j2p,cb);
      // determine number of distinct elements among j1,j2,j1p,j2p
      if(j1==j1p || j2==j2p || j2==j1p) 
      { //printf(": 3 distinct\n");
        r=3; 
        n2=m*m*m;
        if(j1==j1p)
        { ii[1]=j1; ii[2]=j2; ii[3]=j2p;
          // sum over cases of kk[]
          for(iv=0;iv<n2;iv++)
          { d2v(r,m,iv,kk);
            // actually don't need ders here so use ideriv=0
            plgnderghi(n0, m, alp, b, r, ii, kk, &pr, der, 0, x,w,*nq);
            ic1=kk[1]*m+kk[2];
            ic2=kk[1]*m+kk[3];
            for(ip=1;ip<=np;ip++)
            { for(jp=1;jp<=np;jp++)
              { tem=pr*derjk[jk][ic1][ip]*derjk[jk2][ic2][jp]; 
                MM[ip][jp]+=tem; MM[jp][ip]+=tem;
              }
            }
          }
        }
        else if(j2==j2p)
        { ii[1]=j1; ii[2]=j1p; ii[3]=j2p;
          for(iv=0;iv<n2;iv++)
          { d2v(r,m,iv,kk);
            plgnderghi(n0, m, alp, b, r, ii, kk, &pr, der, 0, x,w,*nq);
            ic1=kk[1]*m+kk[3];
            ic2=kk[2]*m+kk[3];
            for(ip=1;ip<=np;ip++)
            { for(jp=1;jp<=np;jp++)
              { tem=pr*derjk[jk][ic1][ip]*derjk[jk2][ic2][jp]; 
                MM[ip][jp]+=tem; MM[jp][ip]+=tem;
              }
            }
          }
        }
        else /* j2==j1p */
        { ii[1]=j1; ii[2]=j1p; ii[3]=j2p;
          for(iv=0;iv<n2;iv++)
          { d2v(r,m,iv,kk);
            plgnderghi(n0, m, alp, b, r, ii, kk, &pr, der, 0, x,w,*nq);
            ic1=kk[1]*m+kk[2];
            ic2=kk[2]*m+kk[3];
            for(ip=1;ip<=np;ip++)
            { for(jp=1;jp<=np;jp++)
              { tem=pr*derjk[jk][ic1][ip]*derjk[jk2][ic2][jp]; 
                MM[ip][jp]+=tem; MM[jp][ip]+=tem;
              }
            }
          }
        }
      }

      else 
      { //printf(": 4 distinct\n");
        r=4; n2=m*m; n2=n2*n2;
        ii[1]=j1; ii[2]=j2; ii[3]=j1p; ii[4]=j2p;
        for(iv=0;iv<n2;iv++)
        { d2v(r,m,iv,kk);
          plgnderghi(n0, m, alp, b, r, ii, kk, &pr, der, 0, x,w,*nq);
          ic1=kk[1]*m+kk[2];
          ic2=kk[3]*m+kk[4];
          for(ip=1;ip<=np;ip++)
          { for(jp=1;jp<=np;jp++)
            { tem=pr*derjk[jk][ic1][ip]*derjk[jk2][ic2][jp]; 
              MM[ip][jp]+=tem; MM[jp][ip]+=tem;
            }
          }
        }
      }
    }
  }
  /*if(*iprint==1)
  { Rprintf("M matrix\n");
    for(ip=1;ip<=np;ip++)
    { for(jp=1;jp<=np;jp++) Rprintf("%f ", MM[ip][jp]);  Rprintf("\n"); }
  }*/
  // add D matrix to MM
  for(ip=1;ip<=np;ip++)
  { for(jp=1;jp<=np;jp++) MM[ip][jp]+=D[ip][jp]; }

  /*if(*iprint==1)
  {  Rprintf("M matrix\n");
    for(ip=1;ip<=np;ip++)
    { for(jp=1;jp<=np;jp++) Rprintf("%f ", MM[ip][jp]);  Rprintf("\n"); }
  }*/

  // invert D matrix, then D^{-1}*M*D^{-1} (since D is symmetric)
  for(jp=1;jp<=np;jp++)
  { for(ip=1;ip<=np;ip++)
    { A[jp][ip]=D[jp][ip]; A[jp][ip+np]=0.; }
    A[jp][jp+np]=1.;
  }

  geppldet(A,np,2*np,&ldet,1.e-12,&parity);
  //if(parity==0 && *iprint==1) { Rprintf("zero determinant for D\n"); }
  //printf("D inverse \n");
  //for(ip=1;ip<=np;ip++)
  //{ for(jp=1;jp<=np;jp++) printf("%f ", A[ip][jp+np]);  printf("\n"); }

  for(jp=1;jp<=np;jp++)
  { for(ip=1;ip<=np;ip++)
    { for(rr=1,sum=0.;rr<=np;rr++) sum+=A[jp][rr+np]*MM[rr][ip];
      A[jp][ip]=sum;
    }
  }

  for(jp=1;jp<=np;jp++)
  { for(ip=1;ip<=np;ip++)
    { for(rr=1,sum=0.;rr<=np;rr++) sum+=A[jp][rr]*A[rr][ip+np];
      V[jp-1][ip-1]=sum;
    }
  }

  /*if(*iprint==1)
  { Rprintf("V matrix, asymptotic version scaled by sample size \n");
    for(ip=1;ip<=np;ip++)
    { for(jp=1;jp<=np;jp++) Rprintf("%f ", V[ip-1][jp-1]);  Rprintf("\n"); }
  }*/

  /* for sample of size nobs, V <- V/nobs to get estimated covariance matrix
     of BCL estimator */
  for(ip=0;ip<np;ip++)
  { for(jp=0;jp<np;jp++)
    {
      V[ip][jp]/=nobs;
      *(Vout + (jp+ip*np)) = V[ip][jp];
    }
  }
  Free(A[0]); Free(A);
  Free(MM[0]); Free(MM);
  Free(D[0]); Free(D);
  Free(derjk[0][0]); Free(derjk[0]); Free(derjk);
  Free(der); Free(cb);
  Free(w); Free(x); 
  Free(b); Free(alp[0]); Free(alp);
}
#endif

#ifndef R
#include <malloc.h>
#else
#include <R.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* computation of M2 statistic, common beta (Rasch model) */
/* derivatives wrt beta parameters done differently */
#define SQRTPI 1.772453850905516

#ifdef MAINM
main(int argc, char *argv[])
{ int i,j,k,np,r,n,n2,m,nm;
  void gauher(double *, double *, int);
  double **dmatrix(int, int);
  /*int **imatrix(int, int);*/
  double **alp,*b,b0;
  double *pp;
  int nrec;
  double **dat,*fr;
  void summ2fr(int, int, int, int, double **, double *, double *, int *);
  int g,ng,z;
  FILE *in;
  double m2;
  double m2rasch(int, int, int, double *, double **alp, double b, int nn,
    double *x, double *w, int nq, int *c2);
  double *x,*w;
  int *c2;
  int nq,nn;

  setbuf(stdout,NULL);
  if(argc==1) { printf("%s parfile [nq] < datafile", argv[0]); exit(0); }
  if(argc>=3) nq=atoi(argv[2]); else nq=48;
  printf("\nnq=%d\n",nq);
  x=(double *) malloc((nq+1) * sizeof(double));
  w=(double *) malloc((nq+1) * sizeof(double));
  gauher(x,w,nq);
  for (j=1;j<=nq;j++) x[j]*=M_SQRT2;
  for (j=1;j<=nq;j++) w[j]/=SQRTPI;
  /*for (j=1;j<=nq;j++) printf("%3d %14.6e %14.6e\n",j,x[j],w[j]);*/

  in=fopen(argv[1],"r");
  fscanf(in, "%d %d %d", &n,&m,&nrec); /* number of items and categories */
  for(i=1,n2=1;i<=n;i++) n2*=m;
  printf("\nn=%d, #categ=%d, nrec=%d\n", n,n2,nrec);

  dat=dmatrix(nrec,n);
  fr=(double *) malloc((nrec) * sizeof(double));
  for(i=0,nn=0;i<nrec;i++)
  { for(j=0;j<n;j++) scanf("%lf", &dat[i][j]); 
    scanf("%lf", &fr[i]); nn+=fr[i];
    /* later add option for all fr=1 */
  }
  /* print out last row as a check*/
  /*for(j=0;j<n;j++) printf("%.0f ", dat[nrec-1][j]); 
  printf("%.0f", fr[nrec-1]);
  printf("\n");*/
  printf("sample size=%d\n", nn);

  nm=n*(m-1) + (n*(n-1)*(m-1)*(m-1))/2;
  pp=(double *) malloc((nm+1) * sizeof(double));
  summ2fr(n,m,nn,nrec,dat,fr,pp,&ng); 
  /* print out column vector of sample moments */
  printf("sample moments as a vector, dim=%d\n", ng);
  for(g=1;g<=ng;g++)
  { printf("%.3f ", pp[g]); if(g%10==0) printf("\n"); }
  printf("\n");

  /* alp[i][1],...,alp[i][m-1] should be in decreasing order */
  alp=dmatrix(n+1,m);
  /*b=(double *) malloc((n+1) * sizeof(double));*/

  for(i=1;i<=n;i++)
  { for(j=1;j<m;j++) fscanf(in, "%lf", &alp[i][j]); }
  fscanf(in, "%lf", &b0);
  printf("\nparameters\n");
  for(i=1;i<=n;i++)
  { for(j=1;j<m;j++) printf("%f ", alp[i][j]); }
  printf("\n");
  printf("%f ", b0); printf("\n");
  fclose(in);
  /*for(i=1;i<=n;i++) b[i]=b0;*/

  c2=(int *) malloc((n+2) * sizeof(int));
  for(i=2,c2[1]=0;i<=n+1;i++) c2[i]=(i*(i-1))/2;
  /* make call to m2rasch */
  m2=m2rasch(n,m,ng,pp,alp,b0,nn,x,w,nq,c2);
  /* continue here */
  np=n*m-n+1;
  printf("M2=%f, df=%d\n", m2,ng-np);
  
  free(w); free(x);
  free(alp[0]); free(alp); /*free(b);*/
  free(dat[0]); free(dat); free(fr);
  free(pp);
  free(c2); 
  printf("===============================================\n");
  exit(0);
}
#endif

/* M2 statistic given sample moments pp[] in lexico order */
double m2rasch(int n, int m, int ng, double *pp, double **alp, double b0, int nn,
  double *x, double *w, int nq, int *c2)
{ void inv2(int, int, int *, int *, int *);
  int z,g,i,i1,i2,k,k1,k2,r,g1,g2,r1,r2,ir,r3,np,nm;
  int *ii,*kk,*cc,*dd,*ss,*tt;
  void raschdergh(int, int, double **, double, int, int *, int *,
    double *, double *, int, double *, double *, int, int *);
  double plgngh(int, int, double **, double *, int, int *, int *, 
    double *, double *, int);  // for ifndef ALL same as m2.c 
  //double plgngh(int, int, double **, double *, int, int *, int *, int,
  // double *, double *, int, int *);  
  double x2,tem,tem12,m2; /*pr[NM],e[NM],del[NM][NP];*/
  double **del,*pr,*e;
  /*int iimat[NM][2],kkmat[NM][2],rvec[NM];*/
  int **iimat,**kkmat,*rvec;
  void merg(int, int, int *, int *, int *, int *, int *, int *, int *);
  double **xi,**cmat;
  void qfmat(int, int, double **, double **, double **);
  double *b;
  double **dmatrix(int,int);
  int **imatrix(int,int);
  /*extern int *c2,nq,nn;*/

  nm=n*(m-1) + (n*(n-1)*(m-1)*(m-1))/2;
  np=n*m-n+1;

#ifndef R
  ii=(int *) malloc((n+1) * sizeof(int));
  kk=(int *) malloc((n+1) * sizeof(int));
  b=(double *) malloc((n+1) * sizeof(double));
  del=dmatrix(nm+1,np+1);
  pr=(double *) malloc((nm+1) * sizeof(double));
  e=(double *) malloc((nm+1) * sizeof(double));
  rvec=(int *) malloc((nm+1) * sizeof(int));
#else
  ii=(int *) Calloc((n+1), int);
  kk=(int *) Calloc((n+1), int);
  b=(double *) Calloc((n+1), double);
  del=dmatrix(nm+1,np+1);
  pr=(double *) Calloc((nm+1), double);
  e=(double *) Calloc((nm+1), double);
  rvec=(int *) Calloc((nm+1), int);
#endif
  iimat=imatrix(nm+1,2);
  kkmat=imatrix(nm+1,2);
  for(i=1;i<=n;i++) b[i]=b0; /* for plgngh */
  g=0;
  for(z=1;z<=c2[n+1];z++)
  { inv2(z,n,&i1,&i2,c2);
    if(z<=n) 
    { r=1;
      ii[1]=i1;
      for(k=1;k<m;k++)
      { kk[1]=k; g++;
        raschdergh(n,m,alp,b0,r,ii,kk,&pr[g],del[g],nn, x,w,nq,c2);
        iimat[g][0]=i1; kkmat[g][0]=k; rvec[g]=r;
      }
    }
    else
    { r=2;
      ii[1]=i1; ii[2]=i2;
      for(k1=1;k1<m;k1++)
      { for(k2=1;k2<m;k2++)
        { kk[1]=k1; kk[2]=k2; g++;
          raschdergh(n,m,alp,b0,r,ii,kk,&pr[g],del[g],nn, x,w,nq,c2);
          iimat[g][0]=i1; iimat[g][1]=i2; 
          kkmat[g][0]=k1; kkmat[g][1]=k2; rvec[g]=r;
        }
      }
    }
  }
  /*printf("g=%d\n", g);*/
#ifdef MAINM
  printf("\nmodel moments as a vector, dim=%d\n", ng);
  for(g=1;g<=ng;g++)
  { printf("%.3f ", pr[g]); if(g%10==0) printf("\n"); }
  printf("\n");
#endif

  /* statistic 1, X2 ignoring dependence */
  /* also set up for solving linear system */
  xi=dmatrix(ng+1,ng+1);
  for(g=1,x2=0.;g<=ng;g++)
  { tem=pp[g]-pr[g]; e[g]=tem; /* xi[g][ng+1]=tem;*/
    tem=tem*tem/pr[g]; x2+=tem; 
  }
  x2*=nn;
  /*printf("x2=%f\n", x2);*/

#ifndef R
  cc=(int *) malloc((n+1) * sizeof(int));
  dd=(int *) malloc((n+1) * sizeof(int));
  ss=(int *) malloc((n+1) * sizeof(int));
  tt=(int *) malloc((n+1) * sizeof(int));
#else
  cc=(int *) Calloc((n+1), int);
  dd=(int *) Calloc((n+1), int);
  ss=(int *) Calloc((n+1), int);
  tt=(int *) Calloc((n+1), int);
#endif
  /* M2: get the covariance matrix */
  for(g1=1;g1<=ng;g1++)
  { /* diagonal element */
    tem=pr[g1]; xi[g1][g1]=tem*(1.-tem);
    r1=rvec[g1];
    for(ir=1;ir<=r1;ir++) 
    { ii[ir]=iimat[g1][ir-1]; kk[ir]=kkmat[g1][ir-1]; }
    for(g2=g1+1;g2<=ng;g2++)
    { r2=rvec[g2];
      for(ir=1;ir<=r2;ir++) 
      { ss[ir]=iimat[g2][ir-1]; tt[ir]=kkmat[g2][ir-1]; }
      merg(r1,r2,ii,kk,ss,tt,&r3,cc,dd);
      //if(r3>0) tem12=plgngh(n,m,alp,b,r3,cc,dd,nn, x,w,nq,c2);
      // for ifndef ALL, same as m2.c
      if(r3>0) tem12=plgngh(n,m,alp,b,r3,cc,dd,x,w,nq); 
      else tem12=0.;
      tem12-=tem*pr[g2];
      xi[g1][g2]=tem12; xi[g2][g1]=tem12;
    }
  }

  /* Then M2= nn* e'*cmat*e */
  cmat=dmatrix(ng+1,ng+1);
  qfmat(ng,np,xi,del,cmat);
  /* print out cmat matrix for checking */

  for(g1=1,m2=0.;g1<=ng;g1++) 
  { tem=e[g1];
    for(g2=1;g2<=ng;g2++) m2+=tem*e[g2]*cmat[g1][g2];
  }
  m2*=nn;
#ifdef MAINM
  printf("\nM2=%f, df=%d\n", m2,ng-np);
#endif
    
#ifndef R
  free(ii); free(kk); free(cc); free(dd); free(ss); free(tt);
  free(b);
  free(del[0]); free(del);
  free(xi[0]); free(xi);
  free(cmat[0]); free(cmat);
  free(pr); free(e);
  free(rvec);
  free(iimat[0]); free(iimat);
  free(kkmat[0]); free(kkmat);
#else
  Free(ii); Free(kk); Free(cc); Free(dd); Free(ss); Free(tt);
  Free(b);
  Free(del[0]); Free(del);
  Free(xi[0]); Free(xi);
  Free(cmat[0]); Free(cmat);
  Free(pr); Free(e);
  Free(rvec);
  Free(iimat[0]); Free(iimat);
  Free(kkmat[0]); Free(kkmat);    
#endif
  return m2;
}

#ifndef ALL
void inv2(int z, int n, int *i1, int *i2, int *c2)
{ int zz,j;
  /*extern int *c2;*/
  if(z<=n) { *i1=z; *i2=0; }
  else 
  { zz=z-n;
    for(j=2;j<=n;j++)
    { if(zz>c2[j]) continue;
      break;
    }
    *i2=j; *i1=zz-c2[j-1];
  }
}
#endif

#ifndef ALL
/* polytomous logit-normit probability for margin */
/*double plgngh(int n, int m, double **alp, double *b, int r, 
   int *ii, int *kk, int nn,
   double *x, double *w, int nq, int *c2)*/
double plgngh(int n, int m, double **alp, double *b, int r, 
   int *ii, int *kk,
   double *x, double *w, int nq)
// for ifndef ALL, last argument *c2 is not needed
{ int j;
  double pgh(double, double **alp, double *b, 
     int m, int r, int *ii,int *kk);
  double fn,sum;
  /*extern double *x,*w;*/
  /*extern int nq;*/

  /* loop for probabilities and partial derivatives */
  for(j=1,sum=0.;j<=nq;j++)
  { fn=pgh(x[j],alp,b,m,r,ii,kk);
    sum+=w[j]*fn;
  }
  /*sum/=SQRTPI;*/
  return sum;
}
#endif

#ifndef ALL
double pgh(double x, double **alp, double *b, int m, int r, 
   int *ii, int *kk)
{ double lg,tem;
  int i,j,k;
  for(j=1,tem=1.;j<=r;j++)
  { i=ii[j]; k=kk[j];
    if(k==0) { lg=1.-1./(1.+exp(-alp[i][1]-b[i]*x)); }
    else if(k==m-1) { lg=1./(1.+exp(-alp[i][m-1]-b[i]*x)); }
    else { lg=1./(1+exp(-alp[i][k]-b[i]*x))-1./(1+exp(-alp[i][k+1]-b[i]*x)); }
    tem*=lg;
  }
  return tem;
}
#endif

/* polytomous logit-normit probability for margin with derivative */
/* derivatives wrt beta parameters done differently */
void raschdergh(int n, int m, double **alp, double b0, int r, 
   int *ii, int *kk, double *pr, double *der, int nn,
  double *x, double *w, int nq, int *c2)
{ int j,np,id;
  void rpghder(double, double **alp, double b0, 
     int n, int m, int r, int *ii,int *kk, double *);
  double *fn;
  /*extern double *x,*w;*/
  /*extern int nq;*/

  np=n*m-n+1;
#ifndef R
  fn=(double *) malloc((np+1) * sizeof(double));
#else
  fn=(double *) Calloc((np+1), double);
#endif

  /* loop for probabilities and partial derivatives */
  for(id=0;id<=np;id++) der[id]=0.;
  for(j=1;j<=nq;j++)
  { rpghder(x[j],alp,b0,n,m,r,ii,kk,fn);
    /* return a vector of fn[], [0] element is the function */
    for(id=0;id<=np;id++) der[id]+=w[j]*fn[id];
  }
  /*for(id=0;id<=np;id++) der[id]/=SQRTPI;*/
  *pr=der[0];
#ifndef R
  free(fn);
#else
  Free(fn);
#endif
}

void rpghder(double x, double **alp, double b0, int n, int m, int r, 
   int *ii, int *kk, double *fn)
{ double lg1,lg2,diff,tem,da1,da2,db,dtem;
  int i,j,k,np,id,np1;

  /*np=m*n; np1=np-n;*/
  np=m*n-n+1; np1=m*n;
  fn[0]=1.;
  for(id=1;id<=np;id++) fn[id]=0.;
  for(j=1,tem=1.;j<=r;j++)
  { i=ii[j]; k=kk[j];
    if(k==0) 
    { lg2=1.; lg1=1./(1.+exp(-alp[i][1]-b0*x)); 
      da2=0.; da1=lg1*(1.-lg1);
      diff=lg2-lg1;
      if(diff<=0.) dtem=0.; else dtem=-da1/diff;
      fn[(i-1)*(m-1)+1]= dtem;
    }
    else if(k==m-1) 
    { lg2=1./(1.+exp(-alp[i][m-1]-b0*x)); lg1=0.; 
      da1=0.; da2=lg2*(1.-lg2);
      diff=lg2-lg1;
      if(diff<=0.) dtem=0.; else dtem=da2/diff;
      fn[(i-1)*(m-1)+m-1]= dtem;
    }
    else 
    { lg2=1./(1+exp(-alp[i][k]-b0*x)); lg1=1./(1+exp(-alp[i][k+1]-b0*x)); 
      da1=lg1*(1.-lg1); da2=lg2*(1.-lg2);
      diff=lg2-lg1;
      if(diff<=0.) dtem=0.; else dtem=-da1/diff;
      fn[(i-1)*(m-1)+k+1]= dtem;
      if(diff<=0.) dtem=0.; else dtem=da2/diff;
      fn[(i-1)*(m-1)+k]= dtem;
    }
    /* beta parameters */
    db=x*(da2-da1);
    if(diff<=0.) dtem=0.; else dtem=db/diff;
    /*fn[np1+i]=dtem;*/
    fn[np]+=dtem;  /* sum over the beta parameters */
    tem*=diff;
  }
  for(id=0;id<=np;id++) fn[id]*=tem;
  /* fn[0]=tem; */
}

#ifndef ALL
/* summarize data into univariate, bivariate margins 
   input can include optional frequencies
   categories assumed to be 0 1 ... m-1
*/
/* get 1D, 2D margins */
void summ2fr(int n, int m, int nn, int nrec, double **dat, double *fr, 
   double *pp, int *ng)
{ int i,j,k,j1,j2,k1,k2,g;
  int *s,**s2;
  int **imatrix(int, int);

#ifndef R
  s=(int *) malloc((m+1) * sizeof(int));
#else
  s=(int *) Calloc((m+1), int);
#endif
    
  s2=imatrix(m+1,m+1);

  g=0;
  /* univariate margins */
  for(j=0;j<n;j++)
  { for(k=0;k<m;k++) s[k]=0;
    for(i=0;i<nrec;i++) 
    { k=dat[i][j]; s[k]+=fr[i]; }

    for(k=1;k<m;k++) { g++; pp[g]=((double)s[k])/nn; }
  }

  /* bivariate margins */
  for(j2=1;j2<n;j2++)
  { for(j1=0;j1<j2;j1++)
    { for(k1=0;k1<m;k1++) 
      { for(k2=0;k2<m;k2++) s2[k1][k2]=0; }
      for(i=0;i<nrec;i++) 
      { k1=dat[i][j1]; k2=dat[i][j2]; s2[k1][k2]+=fr[i]; }

      for(k1=1;k1<m;k1++) 
      { for(k2=1;k2<m;k2++) 
        { g++; pp[g]=((double)s2[k1][k2])/nn; }
      }
    }
  }
  *ng=g;
#ifndef R
  free(s); free(s2[0]); free(s2);
#else
  Free(s); Free(s2[0]); Free(s2);
#endif
}
#endif

#ifndef ALL
/* input I(Y[i1]=j1,...,Y[im_1]=jm_1), I(Y[s1]=t1,...,Y[sm_2]=tm_2);
   merger of common indices
   */
void merg(int m1, int m2, int *ii,int *j, int *ss,int *t, 
   int *m3, int *cc,int *d)
{ int i,i1,i2;
  i1=1; i2=1; i=1; ii[m1+1]=10000; ss[m2+1]=10000;
  while(i1<=m1 || i2<=m2)
  { if(ii[i1]<ss[i2])
    { cc[i]=ii[i1]; d[i]=j[i1]; i++; i1++; }
    else if (ii[i1]>ss[i2])
    { cc[i]=ss[i2]; d[i]=t[i2]; i++; i2++; }
    else
    { cc[i]=ss[i2]; 
      if(j[i1]==t[i2]) d[i]=t[i2];
      else 
      { d[i]=-1; 
        /* in this case expected value of prodict is zero, */
        /* return a negative value for m3 */
        *m3=-1; return;
      }
      i++; i2++; i1++; 
    }
  }
  *m3=i-1; 
  /*printf("m3=%d\n", *m3);*/
  /*for(i=1;i<=*m3;i++) printf("%d ", cc[i]); printf("\n");*/
}
#endif

#ifdef R
void Rm2rasch(int *nitem, int *ncateg, int *nrec, double *dataset, double *alphas,
              double *bvec, double *samplemout, double *m2statrout, double *dfraschout,
              int *nq)
{
   int i,j,np,nm; 
   void gauher(double *, double *, int); 
   double **dmatrix(int, int); 
   double **alp,*b; 
   double *pp; 
   double **dat,*fr; 
   void summ2fr(int, int, int, int, double **, double *, double *, int *); 
   int ng; 
   double m2; 
   double m2rasch(int, int, int, double *, double **alp, double b, int nn, 
     double *x, double *w, int nq, int *c2); 
   int *c2; 
   double *x,*w; 
   int nn; 

   dat=dmatrix(*nrec,*nitem); 
   fr=(double *) Calloc((*nrec), double);  

   for(i=0,nn=0;i<(*nrec);i++)  
   { for(j=0;j<(*nitem);j++) dat[i][j] = *(dataset + (j*(*nrec)+i));   
     fr[i]=*(dataset + (i+(*nitem)*(*nrec))); nn+=fr[i];  
     /*if(i<10) printf("%f %f\n", dat[i][n-1], fr[i]);*/
   }  
  /* print out last row as a check*/
  /*  for(j=1;j<=n;j++) printf("%d ", dat[nrec-1][j]); 
  printf("%d", fr[nrec-1]);
  printf("\n");
  printf("sample size=%d\n", nn);*/

   nm=(*nitem)*((*ncateg)-1) + ((*nitem)*((*nitem)-1)*((*ncateg)-1)*((*ncateg)-1))/2;  
   pp=(double *) Calloc((nm+1), double);  
   summ2fr(*nitem,*ncateg,nn,*nrec,dat,fr,pp,&ng);   
   /* print out column vector of sample moments */
  /*
  printf("sample moments as a vector, dim=%d\n", ng);
  for(g=1;g<=ng;g++)
  { printf("%.3f ", pp[g]); if(g%10==0) printf("\n"); }
  printf("\n");  */

   x=(double *) Calloc((*nq+1), double);  
   w=(double *) Calloc((*nq+1), double);  
   gauher(x,w,*nq);  
   for (j=1;j<=*nq;j++) x[j]*=M_SQRT2;  
   for (j=1;j<=*nq;j++) w[j]/=SQRTPI;  
   /*for (j=1;j<=nq;j++) printf("%3d %14.6e %14.6e\n",j,x[j],w[j]);*/

   /* alp[i][1],...,alp[i][m-1] should be in decreasing order */

   alp=dmatrix(*nitem+1,*ncateg);  
   for(i=0;i<(*nitem);i++)  
   { for(j=0;j<(*ncateg-1);j++)   
     { alp[i+1][j+1]=*(alphas + ((*ncateg-1)*i+j));  
       /*printf("alpha(%d,%d)=%f\n",i,j,alp[i+1][j+1]); */
     }  
   }  

   /* not sure if this is necessary as b is never used 
	see below comment - CF 20100717 */
   b=(double *) Calloc((*nitem+1), double);  
   for(i=0;i<(*nitem);i++)   
   { b[i+1]=*bvec;  
     /*printf("beta(%d)=%f\n",i,b[i]);*/
   }  
    /*printf("\nparameters\n");
  for(i=1;i<=n;i++)
  { for(j=1;j<m;j++) printf("%f ", alp[i][j]); }
  printf("\n");
  for(i=1;i<=n;i++) printf("%f ", b[i]); printf("\n");*/

   c2=(int *) Calloc((*nitem+2), int);  
   for(i=2,c2[1]=0;i<=(*nitem)+1;i++) c2[i]=(i*(i-1))/2;  
  /* make call to m2rasch */
   m2=m2rasch(*nitem,*ncateg,ng,pp,alp,*bvec,nn,x,w,*nq,c2); 

  /* original call to m2rasch used b0, which corresponds
	 to *bvec in current code. So, b is never passed to 
	the m2rasch function and is never used */
/*   m2=m2rasch(n,m,ng,pp,alp,b0,nn,x,w,nq,c2); */

  /* continue here */
   np=(*nitem)*(*ncateg)-(*nitem)+1; 
   /*printf("M2=%f, df=%d\n", m2,ng-np);*/
  
  /* Creation of output */

   for(j=0;j<ng;j++) *(samplemout + j)=pp[j+1]; 

   *m2statrout = m2; 

   *dfraschout = ng-np; 

   Free(w); Free(x); 
   Free(pp);  
   Free(alp[0]); Free(alp); Free(b);  
   Free(dat[0]); Free(dat); Free(fr);  
   Free(c2);   
}

#endif

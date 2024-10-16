#ifndef R
#include <malloc.h>
#else
#include <R.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* gauss-hermite integ for ordinal logit-normal model */
/* computation of M2 statistic */
#define SQRTPI 1.772453850905516

#ifdef MAINM
main(int argc, char *argv[])
{ int i,j,k,np,r,n,m,nm;
  void gauher(double *, double *, int);
  double **alp,*b;
  double *pp;
  int nrec;
  double **dat,*fr;
  double **dmatrix(int, int);
  /*int **imatrix(int, int);*/
  void summ2fr(int, int, int, int, double **, double *, double *, int *);
  int g,ng,z;
  FILE *in;
  double m2,x2,n2;
  double m2stat(int, int, int, double *, double **, double *,
     int *, int, double *, double *, int);
  /*int *ii,*kk;*/
  double x2statb(int, int, int, double **, double *, double **, double *,
    int, double *, double *, int);
  int *c2;
  double *x,*w;
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
  printf("\nn=%d, #categ=%d, nrec=%d\n", n,m,nrec);
  np=m*n;

  dat=dmatrix(nrec,n);
  fr=(double *) malloc((nrec) * sizeof(double));
  for(i=0,nn=0;i<nrec;i++)
  { for(j=0;j<n;j++) scanf("%lf", &dat[i][j]); 
    scanf("%lf", &fr[i]); nn+=fr[i];
    /* later add option for all fr=1 */
  }
  /* print out last row as a check */
  /*for(j=0;j<n;j++) printf("%0.f ", dat[nrec-1][j]); 
    printf("%0.f", fr[nrec-1]); 
    printf("\n"); */
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
  b=(double *) malloc((n+1) * sizeof(double));

  for(i=1;i<=n;i++)
  { for(j=1;j<m;j++) fscanf(in, "%lf", &alp[i][j]); }
  for(i=1;i<=n;i++) fscanf(in, "%lf", &b[i]);
  printf("\nparameters\n");
  for(i=1;i<=n;i++)
  { for(j=1;j<m;j++) printf("%f ", alp[i][j]); }
  printf("\n");
  for(i=1;i<=n;i++) printf("%f ", b[i]); printf("\n");
  np=n*m;

  fclose(in);

  c2=(int *) malloc((n+2) * sizeof(int));
  x2=x2statb(n,m,nrec,dat,fr,alp,b,nn,x,w,nq);
  printf("\nX2=%f, df=%.0f\n", x2,n2-np-1);
  for(i=2,c2[1]=0;i<=n+1;i++) c2[i]=(i*(i-1))/2;
  /* make call to m2stat */
  m2=m2stat(n,m,ng,pp,alp,b,c2,nn,x,w,nq);
  printf("M2=%f, df=%d\n", m2,ng-np);
  
  free(w); free(x);
  free(alp[0]); free(alp); free(b);
  free(dat[0]); free(dat); free(fr);
  free(pp);
  free(c2); 
  printf("===============================================\n");
  exit(0);
}
#endif


/* M2 statistic given sample moments pp[] in lexico order */
/* also compute residuals (standardized) */
double m2stat(int n, int m, int ng, double *pp, double **alp, double *b,
  int *c2, int nn, double *x, double *w, int nq)
{ void inv2(int, int, int *, int *, int *);
  int z,g,i1,i2,k,k1,k2,r,g1,g2,r1,r2,ir,r3,np,nm;
  int *ii,*kk,*cc,*dd,*ss,*tt;
  void plgndergh(int, int, double **, double *, int, int *, int *,
    double *, double *, double *, double *, int);
  double plgngh(int, int, double **, double *, int, int *, int *,
    double *, double *, int);  
  double x2,tem,tem12,m2; /*pr[NM],e[NM],del[NM][NP];*/
  double **del,*pr,*e;
  /*int iimat[NM][2],kkmat[NM][2],rvec[NM];*/
  int **iimat,**kkmat,*rvec;
  void merg(int, int, int *, int *, int *, int *, int *, int *, int *);
  double **xi,**cmat; 
  /*void qfmat(int, int, double [][NM], double[][NP], double [][NM]);*/
  void qfmat(int, int, double **, double **, double **);
  double **dmatrix(int,int);
  int **imatrix(int,int);
  /*extern int *c2,nn;*/

  nm=n*(m-1) + (n*(n-1)*(m-1)*(m-1))/2;
  np=n*m;
    
#ifndef R
  ii=(int *) malloc((n+1) * sizeof(int));
  kk=(int *) malloc((n+1) * sizeof(int));
  del=dmatrix(nm+1,np+1);
  pr=(double *) malloc((nm+1) * sizeof(double));
  e=(double *) malloc((nm+1) * sizeof(double));
  rvec=(int *) malloc((nm+1) * sizeof(int));
#else
  ii=(int *) calloc((n+1), sizeof(int));
  kk=(int *) calloc((n+1), sizeof(int));
  del=dmatrix(nm+1,np+1);
  pr=(double *) calloc((nm+1), sizeof(double));
  e=(double *) calloc((nm+1), sizeof(double));
  rvec=(int *) calloc((nm+1), sizeof(int));   
#endif

  iimat=imatrix(nm+1,2);
  kkmat=imatrix(nm+1,2);
  g=0;
  for(z=1;z<=c2[n+1];z++)
  { inv2(z,n,&i1,&i2,c2);
    if(z<=n) 
    { r=1;
      ii[1]=i1;
      for(k=1;k<m;k++)
      { kk[1]=k; g++;
        plgndergh(n,m,alp,b,r,ii,kk,&pr[g],del[g],x,w,nq);
        /* printf("g=%d\n", g); */
        /*for(i=1;i<=m*n;i++) printf("%f ", del[g][i]); printf("\n");*/
        iimat[g][0]=i1; kkmat[g][0]=k; rvec[g]=r;
      }
    }
    else
    { r=2;
      ii[1]=i1; ii[2]=i2;
      for(k1=1;k1<m;k1++)
      { for(k2=1;k2<m;k2++)
        { kk[1]=k1; kk[2]=k2; g++;
          plgndergh(n,m,alp,b,r,ii,kk,&pr[g],del[g],x,w,nq);
          /*pr[g]=plgngh(n,m,alp,b,r,ii,kk,x,w,nq);*/
          iimat[g][0]=i1; iimat[g][1]=i2; 
          kkmat[g][0]=k1; kkmat[g][1]=k2; rvec[g]=r;
        }
      }
    }
  }

#ifdef MAINM
  /*printf("g=%d\n", g);*/
  printf("\nmodel moments as a vector, dim=%d\n", ng);
  for(g=1;g<=ng;g++)
  { printf("%.3f ", pr[g]); if(g%10==0) printf("\n"); }
  printf("\n");
#endif

  /* statistic 1, X2 ignoring dependence */
  /* also set up for solving linear system */
  xi=dmatrix(ng+1,ng+1);
  for(g=1,x2=0.;g<=ng;g++)
  { tem=pp[g]-pr[g]; e[g]=tem; /* xi[g][ng+1]=tem; */
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
  cc=(int *) calloc((n+1), sizeof(int));
  dd=(int *) calloc((n+1), sizeof(int));
  ss=(int *) calloc((n+1), sizeof(int));
  tt=(int *) calloc((n+1), sizeof(int));
#endif

  /* M2: get the covariance matrix */
  /*xi=dmatrix(ng+1,ng+1);*/
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
  /*printf("\nM2=%f, df=%d\n", m2,ng-np);*/
#ifndef R
  free(ii); free(kk); free(cc); free(dd); free(ss); free(tt);
  free(del[0]); free(del);
  free(xi[0]); free(xi);
  free(cmat[0]); free(cmat);
  free(pr); free(e);
  free(rvec);
  free(iimat[0]); free(iimat);
  free(kkmat[0]); free(kkmat);
#else
  free(ii); free(kk); free(cc); free(dd); free(ss); free(tt);
  free(del[0]); free(del);
  free(xi[0]); free(xi);
  free(cmat[0]); free(cmat);
  free(pr); free(e);
  free(rvec);
  free(iimat[0]); free(iimat);
  free(kkmat[0]); free(kkmat);   
#endif
  return m2;
}

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

/* polytomous logit-normit probability for margin */
double plgngh(int n, int m, double **alp, double *b, int r, 
   int *ii, int *kk, double *x, double *w, int nq)
{ int j;
  double pgh(double, double **, double *, 
     int, int, int *,int *);
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

/* polytomous logit-normit probability for margin with derivative */
void plgndergh(int n, int m, double **alp, double *b, int r, 
   int *ii, int *kk, double *pr, double *der, double *x, double *w, int nq)
{ int j,np,id;
  void pghder(double, double **, double *, 
     int, int, int, int *,int *, double *);
  double *fn;
  /*extern double *x,*w;*/
  /*extern int nq;*/

  np=m*n;
#ifndef R
  fn=(double *) malloc((np+1) * sizeof(double));
#else
  fn=(double *) calloc((np+1), sizeof(double));
#endif
  /* loop for probabilities and partial derivatives */
  for(id=0;id<=np;id++) der[id]=0.;
  for(j=1;j<=nq;j++)
  { pghder(x[j],alp,b,n,m,r,ii,kk,fn);
    /* return a vector of fn[], [0] element is the function */
    for(id=0;id<=np;id++) der[id]+=w[j]*fn[id];
  }
  /*for(id=0;id<=np;id++) der[id]/=SQRTPI;*/
  *pr=der[0];
#ifndef R
  free(fn);
#else
  free(fn); 
#endif
}

void pghder(double x, double **alp, double *b, int n, int m, int r, 
   int *ii, int *kk, double *fn)
{ double lg1,lg2,diff,tem,da1,da2,db,dtem;
  int i,j,k,np,id,np1;

  np=m*n; np1=np-n;
  fn[0]=1.;
  for(id=1;id<=np;id++) fn[id]=0.;
  for(j=1,tem=1.;j<=r;j++)
  { i=ii[j]; k=kk[j];
    if(k==0) 
    { lg2=1.; lg1=1./(1.+exp(-alp[i][1]-b[i]*x)); 
      da2=0.; da1=lg1*(1.-lg1);
      diff=lg2-lg1;
      /* todo : check that this is OK */
      if(diff<=0.) dtem=0.; else dtem=-da1/diff;
      fn[(i-1)*(m-1)+1]= dtem;
    }
    else if(k==m-1) 
    { lg2=1./(1.+exp(-alp[i][m-1]-b[i]*x)); lg1=0.; 
      da1=0.; da2=lg2*(1.-lg2);
      diff=lg2-lg1;
      if(diff<=0.) dtem=0.; else dtem=da2/diff;
      fn[(i-1)*(m-1)+m-1]= dtem;
    }
    else 
    { lg2=1./(1+exp(-alp[i][k]-b[i]*x)); lg1=1./(1+exp(-alp[i][k+1]-b[i]*x)); 
      da1=lg1*(1.-lg1); da2=lg2*(1.-lg2);
      diff=lg2-lg1;
      if(diff<=0.) dtem=0.; else dtem=-da1/diff;
      fn[(i-1)*(m-1)+k+1]= dtem;
      if(diff<=0.) dtem=0.; else dtem=da2/diff;
      fn[(i-1)*(m-1)+k]= dtem;
    }
    db=x*(da2-da1);
    if(diff<=0.) dtem=0.; else dtem=db/diff;
    fn[np1+i]=dtem;
    tem*=diff;
  }
  for(id=0;id<=np;id++) fn[id]*=tem;
  /* fn[0]=tem; */
}


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

  g=0;
#ifndef R
  s=(int *) malloc(m * sizeof(int));
#else
  s=(int *) calloc(m, sizeof(int));
#endif
  s2=imatrix(m,m);
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
  free(s);
  free(s2[0]); free(s2);
#else
  free(s);
  free(s2[0]); free(s2);
#endif
}

/* input I(Y[i1]=j1,...,Y[im_1]=jm_1), I(Y[s1]=t1,...,Y[sm_2]=tm_2);
   merger of common indices
   */
void merg(int m1, int m2, int *ii, int *j, int *ss, int *t, 
   int *m3, int *cc, int *d)
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

#ifdef R

void Rm2 (int *nitem, int *ncateg, int *nrec, double *dataset, double *alphas,
          double *bvec, double *samplemout, double *m2statout, double *dfout,
          int *nq)
 {  
   int i,j,np,nm; 
   void gauher(double *, double *, int); 
   double **alp,*b; 
   double *pp; 
   double **dat,*fr; 
   double **dmatrix(int, int); 
   void summ2fr(int, int, int, int, double **, double *, double *,  int *); 
   int ng; 
   double m2,x2,n2; 
   double m2stat(int, int, int, double *, double **, double *, 
      int *, int, double *, double *, int); 
   double x2statb(int, int, int, double **, double *, double **, double	*,
                  int, double *, double *, int); 
   int *c2; 
   double *x,*w; 
   int nn; 

   for(i=1,n2=1;i<=*nitem;i++) n2*=(*ncateg); 
   /*printf("\nn=%d, #categ=%d, nrec=%d\n", n,m,nrec);*/
   np=(*ncateg)*(*nitem); 

   dat=dmatrix(*nrec,*nitem); 
   fr=(double *) calloc((*nrec), sizeof(double)); 

   for(i=0,nn=0;i<(*nrec);i++) 
   { for(j=0;j<(*nitem);j++) dat[i][j] = *(dataset + (j*(*nrec)+i));  
     fr[i]=*(dataset + (i+(*nitem)*(*nrec))); nn+=fr[i]; 
     /*if(i<10) printf("%f %f\n", dat[i][n-1], fr[i]);*/
   } 
   /* print out last row as a check*/
   /*  for(j=0;j<n;j++) printf("%d ", dat[nrec-1][j]); 
     printf("%d", fr[nrec-1]); 
   printf("\n"); */
   /* printf("sample size=%d\n", nn); */

   nm=(*nitem)*((*ncateg)-1) + ((*nitem)*((*nitem)-1)*((*ncateg)-1)*	((*ncateg)-1))/2; 
   pp=(double *) calloc((nm+1), sizeof(double)); 
   summ2fr(*nitem,*ncateg,nn,*nrec,dat,fr,pp,&ng);  
   /* print out column vector of sample moments */
   /*
   printf("sample moments as a vector, dim=%d\n", ng);
   for(g=1;g<=ng;g++)
   { printf("%.3f ", pp[g]); if(g%10==0) printf("\n"); }
   printf("\n"); */

   /*printf("\nnq=%d\n",nq);*/
   x=(double *) calloc((*nq+1), sizeof(double)); 
   w=(double *) calloc((*nq+1), sizeof(double)); 
   gauher(x,w,*nq); 
   for (j=1;j<=*nq;j++) x[j]*=M_SQRT2; 
   for (j=1;j<=*nq;j++) w[j]/=SQRTPI; 
   /*for (j=1;j<=nq;j++) printf("%3d %14.6e %14.6e\n",j,x[j],w[j]);*/

   /* alp[i][1],...,alp[i][m-1] should be in decreasing order */

   alp=dmatrix((*nitem)+1,(*ncateg)); 
   for(i=0;i<(*nitem);i++) 
   { for(j=0;j<((*ncateg)-1);j++)  
     { alp[i+1][j+1]=*(alphas + (((*ncateg)-1)*i+j)); 
       /*printf("alpha(%d,%d)=%f\n",i,j,alp[i+1][j+1]);*/
     } 
   } 

   b=(double *) calloc(((*nitem)+1), sizeof(double)); 
   for(i=0;i<(*nitem);i++)  
   { b[i+1]=*(bvec + i); 
     /*printf("beta(%d)=%f\n",i,b[i]);*/
   } 
   /*printf("\nparameters\n");
   for(i=1;i<=n;i++)
   { for(j=1;j<m;j++) printf("%f ", alp[i][j]); }
   printf("\n");
   for(i=1;i<=n;i++) printf("%f ", b[i]); printf("\n");*/
   np=(*nitem)*(*ncateg); 

   c2=(int *) calloc(((*nitem)+2), sizeof(int)); 
   /* x2=x2statb(n,m,dat,fr,alp,b,nn,x,w,nq); */
   x2=x2statb(*nitem,*ncateg,*nrec,dat,fr,alp,b,nn,x,w,*nq); 
   /*printf("\nx2=%f, df=%e\n", x2,n2-np-1);*/
   for(i=2,c2[1]=0;i<=(*nitem)+1;i++) c2[i]=(i*(i-1))/2; 
   /* make call to m2stat */
   m2=m2stat(*nitem,*ncateg,ng,pp,alp,b,c2,nn,x,w,*nq); 

   /*x2=x2statb(n,m,dat,fr,alp,b);*/
   /*printf("\nx2=%f, df=%e\n", x2,n2-np-1);*/


   /* Creation of output */

   for(j=0;j<ng;j++) *(samplemout + j)=pp[j+1]; 
   *m2statout = m2; 
   *dfout = ng-np; 

   free(w); free(x); 
   free(alp[0]); free(alp); free(b); 
   free(dat[0]); free(dat); free(fr); 
   free(pp); 
   free(c2);  
}

#endif

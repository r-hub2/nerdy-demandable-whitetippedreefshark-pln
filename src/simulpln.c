#ifndef R
#include <malloc.h>
#else
#include <Rmath.h>
#include <R.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*#include <malloc.h>*/
/* simulation of data for polytomous logit-normit model 
   with discretized normal latent variable that corresponds
   to Gauss-Hermite integration */
/* gcc -DCALONE -o simulpln simulpln.c ssgauher.c amatrix.c -lm */
#define SQRTPI 1.772453850905516
#ifdef CALONE
int main(int argc, char *argv[])
{ double **dmatrix(int, int);
  void gauher(double *, double *, int);
  int j,k,n2,i,ir,jr,neq,iq;
  double **alp,*b;
  int seed,nsim;
  double z,u,tem;
  /*double rnorm(void);*/
  double *cdfw;
  double *x,*w;
  int nq,n,nrec,m; 
  double **dat,*fr;
  int nn;

  nq=48;
  /*printf("\nnq=%d\n",nq);*/
  x=(double *) malloc((nq+1) * sizeof(double));
  w=(double *) malloc((nq+1) * sizeof(double));
  cdfw=(double *) malloc((nq+1) * sizeof(double));
  gauher(x,w,nq);
  for (j=1;j<=nq;j++) x[j]*=M_SQRT2;
  for (j=1;j<=nq;j++) w[j]/=SQRTPI;
  /*for (j=1;j<=nq;j++) printf("%3d %14.6e %14.6e\n",j,x[j],w[j]); */
  /* GH is like probability mass w[j] on x[j] as approx */
  for(i=1,cdfw[0]=0.;i<=nq;i++) cdfw[i]=cdfw[i-1]+w[i];
  
  /* number of items and categories */
  printf("Enter n=#items, m=#categories, nn=sample size, seed\n");
  scanf("%d %d %d %d", &n,&m,&nn,&seed); 
  for(i=1,n2=1;i<=n;i++) n2*=m;
  printf("\nn=%d, m=%d, m^n=%d, nn=%d, seed=%d\n", n,m,n2,nn,seed);
  dat=dmatrix(nn,n);
  fr=(double *) malloc((nn) * sizeof(double));
  srand(seed);
  /* alp[i][1],...,alp[i][m-1] should be in decreasing order */
  alp=dmatrix(n+1,m);
  b=(double *) malloc((n+1) * sizeof(double));
  for(i=1;i<=n;i++)
  { for(j=1;j<m;j++) scanf("%lf", &alp[i][j]); }
  for(i=1;i<=n;i++) { scanf("%lf", &b[i]); }
  
  printf("Parameter values alp[i][j], j=1..m-1, i=1..n; b[i]\n");
  for(i=1;i<=n;i++)
  { for(j=1;j<m;j++) printf("%f ", alp[i][j]); printf("\n"); }
  for(i=1;i<=n;i++) printf("%f ", b[i]); printf("\n");

  /* simulate data in dat[][], nn=total sample size */
  for(i=0;i<nn;i++)
  { 
#ifdef NO
    /* use the next line for actual normal latent variables */
    /*z=rnorm(); */
#endif
    u=rand()/2147483648.;
    /* replace by GH approx, sample from x[] with pr w[] */
    for(iq=1;iq<=nq;iq++) { if(u<=cdfw[iq]) break; }
    z=x[iq];
    for(j=1;j<=n;j++)
    { u=rand()/2147483648.;
      for(k=0;k<m-1;k++)
      { tem=alp[j][k+1]+b[j]*z; tem=1./(1.+exp(-tem));
        if(u>tem) break;
      }
      dat[i][j-1]=k;
    }
    fr[i]=1;
  }
  nrec=nn;
  /*
  for(i=1;i<=nrec;i++)
  { for(j=1;j<=n;j++) printf("%d ", dat[i][j]);
    printf(" : %d\n", fr[i]);
  } */

  /* condense dat[], check for duplicates, this could be inefficient */
  for(ir=1;ir<nn;ir++)
  { for(jr=0;jr<ir;jr++)
    { for(i=0,neq=0;i<n;i++) neq+=(dat[ir][i]==dat[jr][i]);
      if(neq==n) { fr[jr]++; fr[ir]=0; break; }
    }
  }
  /* printf("converted data with frequencies\n"); */
  for(i=0;i<nrec;i++)
  { if(fr[i]>0)
    { for(j=0;j<n;j++) printf("%.0f ", dat[i][j]);
      printf("   %.0f\n", fr[i]);
    }
  }
  free(w); free(x); free(cdfw);
  free(alp[0]); free(alp); free(b);
  free(dat[0]); free(dat); free(fr);
  exit(0);
}
#endif

#ifdef NO
double rnorm(void)  
{ double x,y,z;
  do
  { x=(random()/2147483648.0)*2.0-1.0;
    y=(random()/2147483648.0)*2.0-1.0;
    z=x*x+y*y;
  }
  while(z>=1);
  z=x*sqrt(-2.0*log(z)/z);
  return(z);
}
#endif

#ifdef R
/* Without being in an R package, compilation is something like the following,
   compilation with -DDEBUG allows printing in the C code when called from R
    gcc -DR -I/usr/include/R -fpic -c simulpln.c ssgauher.c amatrix.c
or
    gcc -DR -DDEBUG -I/usr/include/R -fpic -c simulpln.c ssgauher.c amatrix.c
    gcc -I/usr/include/R -shared -o simulpln.so simulpln.o ssgauher.o amatrix.o -lRmath -lm
*/

void Rsimulpln( int *nitems0, int *ncateg0, int *nsampsize0, double *alpvec, double *beta0, double *Rdat) 
{ double **dmatrix(int, int);
  void gauher(double *, double *, int);
  int j,k,i,ir,jr,neq,iq,ii;
  double **alp,*b;
  double z,u,tem;
  double *cdfw;
  double *x,*w;
  int nq,nitems,nrec,ncateg,nsampsize;
  double **dat,*fr,**datfr;

  nq=48;
  x=(double *) Calloc((nq+1), double);
  w=(double *) Calloc((nq+1), double);
  cdfw=(double *) Calloc((nq+1), double);
  gauher(x,w,nq);
  for (j=1;j<=nq;j++) x[j]*=M_SQRT2;
  for (j=1;j<=nq;j++) w[j]/=SQRTPI;
  /*for (j=1;j<=nq;j++) printf("%3d %14.6e %14.6e\n",j,x[j],w[j]); */
  /* GH is like probability mass w[j] on x[j] as approx */
  for(i=1,cdfw[0]=0.;i<=nq;i++) cdfw[i]=cdfw[i-1]+w[i];
  
  /* number of items , categories and sample size */
  nitems= *nitems0;
  ncateg= *ncateg0;
  nsampsize= *nsampsize0;

  dat=dmatrix(nsampsize,nitems);
  fr=(double *) Calloc((nsampsize), double);

  /* alp[i][1],...,alp[i][ncateg-1] should be in decreasing order */
  alp=dmatrix(nitems+1,ncateg);
  for(i=0;i<nitems;i++)
  { for(j=0;j<(ncateg-1);j++) alp[i+1][j+1]=alpvec[(ncateg-1)*i+j]; }

  b=(double *) Calloc((nitems+1), double);
  for(i=0;i<nitems;i++) b[i+1]=beta0[i];

  GetRNGstate();
  /* simulate data in dat[][], nsampsize=total sample size */
  for(i=0;i<nsampsize;i++)
  { /* GH approx, sample from x[] with pr w[] */
    u=unif_rand();
    for(iq=1;iq<=nq;iq++) { if(u<=cdfw[iq]) break; }
    z=x[iq];
    for(j=1;j<=nitems;j++)
    { 
      u=unif_rand();  // use rng in R
      for(k=0;k<ncateg-1;k++)
      { tem=alp[j][k+1]+b[j]*z; tem=1./(1.+exp(-tem));
        if(u>tem) break;
      }
      dat[i][j-1]=k;
    }
    fr[i]=1;
  }
  PutRNGstate();
  nrec=nsampsize;
  /* condense dat[], check for duplicates, this could be inefficient */
  for(ir=1;ir<nsampsize;ir++)
  { for(jr=0;jr<ir;jr++)
    { for(i=0,neq=0;i<nitems;i++) neq+=(dat[ir][i]==dat[jr][i]);
      if(neq==nitems) { fr[jr]++; fr[ir]=0; nrec--; break; }
    }
  }
  /* printf("converted data with frequencies\n"); */
  datfr=dmatrix(nrec,nitems+1);
  ii=0;
  for(i=0;i<nsampsize;i++)
  { if(fr[i]>0)
    { for(j=0;j<nitems;j++) datfr[ii][j]=dat[i][j];
      datfr[ii][nitems]=fr[i];
      ii++;
    }
  }

  /* C stores by rows, R stores by columns */
  /*for(i=0;i<nrec;i++)
  { for(j=0;j<=nitems;j++) Rdat[i+j*nrec]=datfr[i][j]; }*/

  /* above changed to below, CF, 2013-01-24*/
  /* this makes trimming 0's from Rdat easier in R
   by instead importing the data into a matrix by rows*/
  ii=0;
  for(i=0;i<nrec;i++)
  { for(j=0;j<=nitems;j++)
    {
      Rdat[ii]=datfr[i][j];
      ii++;
    }
  }
    
    
  Free(w); Free(x); Free(cdfw);
  Free(alp[0]); Free(alp); Free(b);
  Free(dat[0]); Free(dat); Free(fr);
  Free(datfr[0]); Free(datfr);  
}
#endif

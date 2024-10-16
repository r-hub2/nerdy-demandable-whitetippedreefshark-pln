#ifndef R
#include <malloc.h>
#else
#include <R.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* log-likelihood for polytomous logit-normit model with common beta */
/* (Rasch model) */
/* read from parameter file in argv[1] including starting point */

#define SQRTPI 1.772453850905516

#ifdef DATA

int main(int argc, char *argv[])
{ double *param,*lb,*ub,fnval,**invhes,bdd;
  double **dmatrix(int, int);
  double ***gmatrix(int, int, int);
  void gauher(double *, double *, int);
  void raschnllk(int, double *, double *, double **, int,
      int, int, int, int, double *, double *, double **, double *,
      double ***, double ***, double ***);
  int j,k,np,r,n2,i,ip,nn,jp;
  int mxiter,iprint,iconv;
  double **alp,b;
  FILE *in;
  void nrmin(int,double *, void (*)(int, double *, double *, double **, int,
      int, int, int, int, double *, double *, double **, double *,
      double ***, double ***, double ***), 
      double *lb, double *ub, int, double,int, double bdd, 
      int *, double *, double **,
      int, int, int, int, double *, double *, double **, double *,
      double ***, double ***, double ***);
  double *x,*w;
  int nq,n,nrec,m;
  double **dat,*fr;
  double ***g,***g1,***g2;

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
  mxiter=20;
  fscanf(in, "%d %d %d", &n,&m,&nrec); /* number of items and categories */
  for(i=1,n2=1;i<=n;i++) n2*=m;
  printf("\nn=%d, #categ=%d, nrec=%d\n", n,n2,nrec);
  /*np=m*n;*/
  np=m*n-n+1;  /* common beta */
  lb=(double *) malloc(np * sizeof(double));
  ub=(double *) malloc(np * sizeof(double));
  for(ip=0;ip<np;ip++) { lb[ip]=-10; ub[ip]=10; }
  ip=0;

  /* alp[i][1],...,alp[i][m-1] should be in decreasing order */
  alp=dmatrix(n+1,m+1);
  param=(double *) malloc((np+1) * sizeof(double));
  for(i=1;i<=n;i++)
  { for(j=1;j<m;j++) 
    { fscanf(in, "%lf", &alp[i][j]); param[ip]=alp[i][j]; ip++; }
  }
  fscanf(in, "%lf", &b); param[ip]=b; ip++; 
  fclose(in);

  //dat=dmatrix(pow(m,n)+1,n+1);
  dat=dmatrix(nrec,n);
  //fr=(double *) malloc((pow(m,n)+1) * sizeof(double));
  fr=(double *) malloc((nrec) * sizeof(double));
  for(i=0,nn=0;i<nrec;i++)
  { for(j=0;j<n;j++) scanf("%lf", &dat[i][j]); 
    scanf("%lf", &fr[i]); nn+=fr[i]; 
    /*for(j=0;j<n;j++) printf("%d ", dat[i][j]); */
    /*printf(": %d\n", fr[i]);*/
  }
  /*for(j=0;j<n;j++) printf("%f ", dat[nrec-1][j]);*/
  printf("sample size=%d\n", nn);

  invhes=dmatrix(np,np);
  iprint=1;
  bdd=5.;
  g=gmatrix(n+1,m+1,nq+1);
  g1=gmatrix(n+1,m+1,nq+1);
  g2=gmatrix(n+1,m+1,nq+1);

  nrmin(np,param,raschnllk,lb,ub,mxiter,1.e-4,iprint,bdd,&iconv,&fnval,invhes,
    n, m, nq, nrec, x, w, dat, fr, g, g1, g2);
  printf("negative log-likelihood = %f\n", fnval);
  printf("Parameter estimates alp[i][j], j=1..m-1, i=1..n; b\n");
  for(ip=0;ip<np;ip++) 
  { printf("%f ", param[ip]); if(ip%8==7) printf("\n"); }
  printf("\n");
  printf("SEs alp[i][j], j=1..m-1, i=1..n; b\n");
  for(ip=0;ip<np;ip++) 
  { printf("%f ", sqrt(invhes[ip][ip])); if(ip%8==7) printf("\n"); }
  printf("\n");
  /* comment out for testing */
  printf("\ninverse Hessian\n");
  for(ip=0;ip<np;ip++) 
  { for(jp=0;jp<np;jp++) printf("%.10f ", invhes[ip][jp]); printf("\n"); }
  
  free(invhes[0]); free(invhes);
  free(w); free(x);
  free(dat[0]); free(dat); free(fr);
  free(param); free(lb); free(ub);
  free(alp[0]); free(alp);
  free(g[0][0]); free(g[0]); free(g);
  free(g1[0][0]); free(g1[0]); free(g1);
  free(g2[0][0]); free(g2[0]); free(g2); 
  return 0;
}
#endif

void raschnllk(int np, double *param, double *nllk, double **dd, int iprint,
    int n, int m, int nq, int nrec, double *x, double *w, double **dat, double *fr,
    double ***g, double ***g1, double ***g2)
{ int i,j,ir,ip,iq,jp,nfr,npful,ipp,jpp;
  double **alp,*bvec;
  double **dmatrix(int, int);
  void plgngh2(int, int, int[], double *, double *, double **,
    int, double *, double *,
    double ***, double ***, double ***);
  double tem,tem2,xx,bb,pr;
  int *kk;
  double *der, **hes;
  /*extern double *x,*w;*/
  /*extern int nq,n,nrec,m;*/
  /*extern double **dat,*fr;*/
  /*extern double ***g,***g1,***g2;*/
 
  alp=dmatrix(n+1,m);
#ifndef R
  bvec=(double *) malloc((n+1) * sizeof(double));
#else
  bvec=(double *) Calloc((n+1), double);
#endif
  npful=m*n;
  ip=0;
  for(i=1;i<=n;i++)
  { for(j=1;j<m;j++) 
    { alp[i][j]=param[ip]; ip++; }
  }
  for(i=1;i<=n;i++) { bvec[i]=param[ip]; }   /* common beta */
  ip++;

  for(ip=0;ip<np;ip++)
  { for(jp=0;jp<=np;jp++) dd[ip][jp]=0.; }

  /* set up global arrays when parameters are fixed in an iteration */
  for(i=1;i<=n;i++)
  { bb=bvec[i];
    for(iq=1;iq<=nq;iq++)
    { g[i][0][iq]=1.; g[i][m][iq]=0.;
      g1[i][0][iq]=0.; g1[i][m][iq]=0.;
      g2[i][0][iq]=0.; g2[i][m][iq]=0.;
      xx=x[iq];
      for(j=1;j<m;j++) 
      { tem=1./(1.+exp(-alp[i][j]-bb*xx)); 
        g[i][j][iq]=tem;
        tem2=tem*(1.-tem);
        g1[i][j][iq]=tem2;
        g2[i][j][iq]=tem2*(1.-2.*tem);
      }
    }
  }
#ifndef R
  kk=(int *) malloc((n+1) * sizeof(int));
  der=(double *) malloc((npful+1) * sizeof(double));
#else
  kk=(int *) Calloc((n+1), int);
  der=(double *) Calloc((npful+1), double);
#endif
  hes=dmatrix(npful+1,npful+1);

  /* ii[] is vector of indices and kk[] is vector of categories */
  /*for(i=1;i<=n;i++) ii[i]=[i];*/
  for(ir=0,*nllk=0.;ir<nrec;ir++)
  { for(i=1;i<=n;i++) kk[i]=dat[ir][i-1];
    nfr=fr[ir];
    if(nfr==0) continue;  /* allow for zero counts in input data */
    plgngh2(n,m,kk,&pr,der,hes,
    nq,x,w,g,g1,g2);
    /* don't update hessian if pr<=0 */
    if(pr<=0.) { pr=1.e-10; /*printf("ir=%d, pr=%f\n", ir,pr);*/ }
    else
    { /* Hessian matrix plus gradient, set up for nrmin in matrix dd[][] */
      /* to handle common b: indices np..npful are for b */
      for(ip=1;ip<=npful;ip++)
      { ipp=ip; if(ip>np) ipp=np;
        dd[ipp-1][np]-= nfr*der[ip]/pr;
        for(jp=1;jp<=npful;jp++)
        { jpp=jp; if(jp>np) jpp=np;
          tem=hes[ip][jp]-der[ip]*der[jp]/pr;
          dd[ipp-1][jpp-1]-= nfr*tem/pr;
        }
      }
    }
    *nllk -= nfr*log(pr); 
    /*printf("ir=%d, pr=%e, nllk=%f\n", ir,pr,*nllk);*/
  }

  /* for testing print out gradient */
#ifndef R
  if(iprint==1)
  { printf("\ngradient\n");
    for(ip=0;ip<np;ip++) printf("%f ", dd[ip][np]); 
    printf("\n\n");
  }
#else
  if(iprint==1)
  { Rprintf("\ngradient\n");
    for(ip=0;ip<np;ip++) Rprintf("%f ", dd[ip][np]); 
    Rprintf("\n\n");
  }   
#endif

#ifndef R
  free(alp[0]); free(alp); free(bvec);
  free(der); free(hes[0]); free(hes);
  free(kk);
#else
  Free(alp[0]); Free(alp); Free(bvec);
  Free(der); Free(hes[0]); Free(hes);
  Free(kk);
#endif
}

#ifdef R

void Rnrmlerasch( int *nitem, int *ncateg, int *nrec, double *dataset, double *alphas,
                 double *bvec, double *abound, double *bbound, double *nllkout,
                 double *mleraschout, double *sevecout, double *invhesout, int *nq,
                 int *mxiter, int *iconv, int *iprint)
{ 
  double *param,*lb,*ub,fnval,**invhes,bdd;
  double **dmatrix(int, int);
  double ***gmatrix(int, int, int);
  void gauher(double *, double *, int);
  void raschnllk(int, double *, double *, double **, int,
      int, int, int, int, double *, double *, double **, double *,
      double ***, double ***, double ***);
  int j,np,n2,i,ip,nn,jp;
  double **alp,b;
  void nrmin(int,double *, void (*)(int, double *, double *, double **, int,
      int, int, int, int, double *, double *, double **, double *,
      double ***, double ***, double ***), 
      double *lb, double *ub, int, double,int, double bdd, 
      int *, double *, double **,
      int, int, int, int, double *, double *, double **, double *,
      double ***, double ***, double ***);
  double *x,*w;
  double **dat,*fr;
  double ***g,***g1,***g2;

  x=(double *) Calloc((*nq+1), double);
  w=(double *) Calloc((*nq+1), double);
  gauher(x,w,*nq);
  for (j=1;j<=*nq;j++) x[j]*=M_SQRT2;
  for (j=1;j<=*nq;j++) w[j]/=SQRTPI; 
  /*for (j=1;j<=nq;j++) printf("%3d %14.6e %14.6e\n",j,x[j],w[j]);*/
  
  /* convert to matrices in C (i.e. row/column transpose) */
  dat=dmatrix(*nrec,*nitem);
  fr=(double *) Calloc((*nrec), double);
  /* nn = total of fr[] */
  for(i=0,nn=0;i<(*nrec);i++)
  { for(j=0;j<(*nitem);j++) dat[i][j] = *(dataset + (j*(*nrec)+i)); 
    fr[i]=*(dataset + (i+(*nitem)*(*nrec))); nn+=fr[i];
    /*if(i<10) printf("%f %f\n", dat[i][n-1], fr[i]);*/
  }

  for(i=1,n2=1;i<=(*nitem);i++) n2*=(*ncateg);
  if(*iprint==1)
  {
    Rprintf("\nn=%d, #categ=%d, nrec=%d\n", *nitem,*ncateg,*nrec);
  }
  np=(*ncateg)*(*nitem)-(*nitem)+1;  /* common beta */

  lb=(double *) Calloc(np, double);
  ub=(double *) Calloc(np, double);
  for(ip=0;ip<np;ip++) { lb[ip]=*(abound+0); ub[ip]=*(abound+1); } 
  /* add boundary to common slope ? */
  lb[np-1]=*(bbound+0);
  ub[np-1]=*(bbound+1);
  /*lb[np-1]=0.;*/

  /* alp[i][1],...,alp[i][m-1] should be in decreasing order */
  ip=0;
  alp=dmatrix((*nitem)+1,(*ncateg));
  param=(double *) Calloc(np, double);
  for(i=0;i<(*nitem);i++)
  { for(j=0;j<((*ncateg)-1);j++)  
    { alp[i+1][j+1]=*(alphas + (((*ncateg)-1)*i+j)); param[ip]=alp[i+1][j+1]; ip++;
      /*if(*iprint==1)
      {
        Rprintf("alpha(%d,%d)=%f\n",i,j,alp[i+1][j+1]);
      }*/
    }
  }
  b=*bvec;
  param[ip]=b; ip++;

  invhes=dmatrix(np,np); 
  bdd=5.;  /* bound the difference between iterations */

  g=gmatrix((*nitem)+1,(*ncateg)+1,*nq+1);
  g1=gmatrix((*nitem)+1,(*ncateg)+1,*nq+1);
  g2=gmatrix((*nitem)+1,(*ncateg)+1,*nq+1);

  nrmin(np,param,raschnllk,lb,ub,*mxiter,1.e-4,*iprint,bdd,iconv,&fnval,invhes,
    *nitem, *ncateg, *nq, *nrec, x, w, dat, fr, g, g1, g2);

  /* below used for testing only */
  /*printf("negative log-likelihood = %f\n", fnval);
  printf("Parameter estimates alp[i][j], j=1..m-1, i=1..n; b\n");
  for(ip=0;ip<np;ip++) 
  { printf("%f ", param[ip]); if(ip%8==7) printf("\n"); }
  printf("\n");
  printf("SEs alp[i][j], j=1..m-1, i=1..n; b\n");
  for(ip=0;ip<np;ip++) 
  { printf("%f ", sqrt(invhes[ip][ip])); if(ip%8==7) printf("\n"); }
  printf("\n");
  */

  /*printf("inverse Hessian\n");
  for(ip=0;ip<np;ip++) 
  { for(jp=0;jp<np;jp++) printf("%.10f ", invhes[ip][jp]); printf("\n"); }*/


  /* Creation of output */

  *nllkout = fnval;

  for(ip=0;ip<np;ip++) *(mleraschout + ip)=param[ip];

  for(ip=0;ip<np;ip++) *(sevecout + ip)=sqrt(invhes[ip][ip]); 

  for(ip=0;ip<np;ip++) 
  { for(jp=0;jp<np;jp++) *(invhesout + (jp+ip*np))=invhes[ip][jp];}


  Free(invhes[0]); Free(invhes);
  Free(w); Free(x);
  Free(dat[0]); Free(dat); Free(fr);
  Free(param); Free(lb); Free(ub);
  Free(alp[0]);
  Free(alp);
  Free(g[0][0]);
  Free(g1[0][0]);
  Free(g2[0][0]);
  Free(g[0]); Free(g);
  Free(g1[0]); Free(g1);
  Free(g2[0]); Free(g2);
}
#endif

#ifndef R
#include <malloc.h>
#else
#include <R.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* log-likelihood for polytomous logit-normit model
   using Newton Raphson with computed Hessian */

#define SQRTPI 1.772453850905516

#ifdef DATA
/* gcc -DDATA -o nrmlepln nrmlepln.c nrmin.o amatrix.o ssgauher.c polyder2.c -lm */

int main(int argc, char *argv[])
{ double *param,*lb,*ub,fnval,**invhes,bdd;
  double **dmatrix(int, int);
  double ***gmatrix(int, int, int);
  void gauher(double *, double *, int);
  void plnnllk(int, double *, double *, double **, int,
      int, int, int, int, double *, double *, double **, double *,
      double ***, double ***, double ***);
  int j,np,n2,i,ip,jp,nn;
  int mxiter,iprint,iconv;
  double **alp,*b,tmp;
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
  double ***g, ***g1, ***g2;

  setbuf(stdout,NULL);
  if(argc==1) { printf("%s parfile [nq] < datafile", argv[0]); exit(0); }
  if(argc>=3) nq=atoi(argv[2]); else nq=48;
  if(nq>48) nq=48;
  printf("\nnq=%d\n",nq);
  x=(double *) malloc((nq+1) * sizeof(double));
  w=(double *) malloc((nq+1) * sizeof(double));
  gauher(x,w,nq);
  for (j=1;j<=nq;j++) x[j]*=M_SQRT2;
  for (j=1;j<=nq;j++) w[j]/=SQRTPI;
  /*for (j=1;j<=nq;j++) printf("%3d %14.6e %14.6e\n",j,x[j],w[j]);*/

  in=fopen(argv[1],"r");
  mxiter=20;
  /*mxiter=1;*/
  fscanf(in, "%d %d %d", &n,&m,&nrec); /* number of items and categories */
  for(i=1,n2=1;i<=n;i++) n2*=m;
  printf("\nn=%d, #categ=%d, nrec=%d\n", n,m,nrec);
  np=m*n;

  ub=(double *) malloc(np * sizeof(double));
  lb=(double *) malloc(np * sizeof(double));
  for(ip=0;ip<np;ip++) { lb[ip]=-10; ub[ip]=10; }
  /* cannot put lower bound of 0 on betas unless
     all items oriented positively */
  for(ip=np-n;ip<np;ip++) { lb[ip]=-1; }
    /* TO DO: set bounds as an argument passed from R*/  
    /*for(ip=np-n;ip<np;ip++) { lb[ip]=0; }*/
  ip=0;

  /* alp[i][1],...,alp[i][m-1] should be in decreasing order */
  alp=dmatrix(n+1,m);
  b=(double *) malloc((n+1) * sizeof(double));

  param=(double *) malloc(np * sizeof(double));
  for(i=1;i<=n;i++)
  { for(j=1;j<m;j++) 
    { fscanf(in, "%lf", &tmp); alp[i][j]=tmp;
      param[ip]=tmp; ip++; }
  }

  for(i=1;i<=n;i++) 
  { fscanf(in, "%lf", &tmp); b[i]=tmp; param[ip]=tmp; ip++; }
  fclose(in);
 
  dat=dmatrix(nrec,n);
  fr=(double *) malloc((nrec) * sizeof(double));
  for(i=0,nn=0;i<nrec;i++)
  { for(j=0;j<n;j++) scanf("%lf", &dat[i][j]);
    scanf("%lf", &fr[i]); nn+=fr[i];
    /*for(j=0;j<n;j++) printf("%d ", dat[i][j]); 
      printf(": %d\n", fr[i]); */
  }
  /*for(j=0;j<n;j++) printf("%f ", dat[nrec-1][j]); */
  printf("sample size=%d\n", nn);

  invhes=dmatrix(np,np);
  iprint=1;
  bdd=5.;  /* bound the difference between iterations */
  g=gmatrix(n+1,m+1,nq+1);
  g1=gmatrix(n+1,m+1,nq+1);
  g2=gmatrix(n+1,m+1,nq+1);

  nrmin(np,param,plnnllk,lb,ub,mxiter,1.e-4,iprint,bdd,&iconv,&fnval,invhes,
    n, m, nq, nrec, x, w, dat, fr, g, g1, g2);
  printf("negative log-likelihood = %f\n", fnval);
  printf("Parameter estimates alp[i][j], j=1..m-1, i=1..n; b[i]\n");
  for(ip=0;ip<np;ip++) 
  { printf("%f ", param[ip]); if(ip%8==7) printf("\n"); }
  printf("\n");
  printf("SEs alp[i][j], j=1..m-1, i=1..n; b[i]\n");
  for(ip=0;ip<np;ip++) 
  { printf("%f ", sqrt(invhes[ip][ip])); if(ip%8==7) printf("\n"); }
  printf("\n");
  /* can comment out for testing */
  printf("\ninverse Hessian\n");
  for(ip=0;ip<np;ip++) 
  { for(jp=0;jp<np;jp++) printf("%.10f ", invhes[ip][jp]); printf("\n"); }
  

  free(invhes[0]); free(invhes);
  free(alp[0]); free(alp);
  free(b);
  free(dat[0]); free(dat);
  free(fr);
  free(param);
  free(g[0][0]); free(g[0]); free(g);
  free(g1[0][0]); free(g1[0]); free(g1);
  free(g2[0][0]); free(g2[0]); free(g2); 
  free(w); free(x);
  free(ub); free(lb);
  exit(0);
}
#endif

void plnnllk(int np, double *param, double *nllk, double **dd, int iprint,
    int n, int m, int nq, int nrec, double *x, double *w, double **dat, double *fr,
    double ***g, double ***g1, double ***g2)
{ int i,j,ir,ip,iq,jp,nfr;
  void plgngh2(int, int, int[], double *, double *, double **,
    int, double *, double *,
    double ***, double ***, double ***);
  double **dmatrix(int, int);
  double tem,tem2,xx,bb,pr;
  int *kk;
  double *der,**hes;
  double **alp,*b;
  /*extern double *x,*w,**dat,*fr;*/
  /*extern int nq,n,nrec,m;*/

  alp=dmatrix(n+1,m);
#ifndef R
  b=(double *) malloc((n+1) * sizeof(double));
#else
  b=(double *) calloc((n+1), sizeof(double));
#endif
  ip=0;
  for(i=1;i<=n;i++)
  { for(j=1;j<m;j++) 
    { alp[i][j]=param[ip]; ip++; }
  }
  for(i=1;i<=n;i++) { b[i]=param[ip]; ip++; }

  for(ip=0;ip<np;ip++)
  { for(jp=0;jp<=np;jp++) dd[ip][jp]=0.; }

  /* set up global arrays when parameters are fixed in an iteration */
  for(i=1;i<=n;i++)
  { bb=b[i];
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
  der=(double *) malloc((np+1) * sizeof(double));
#else
  kk=(int *) calloc((n+1), sizeof(int));
  der=(double *) calloc((np+1), sizeof(double));
#endif
  hes=dmatrix(np+1,np+1);

  for(ir=0,*nllk=0.;ir<nrec;ir++)
  { for(i=1;i<=n;i++) kk[i]=dat[ir][i-1];
    nfr=fr[ir];
    if(nfr==0) continue;
    /*pr=plgngh(n,m,alp,b,n,ii,kk,nq);*/
    /* simplify as ii is fixed here and r=n */
    /*pr=plgngh(n,m,alp,b,kk);*/
    plgngh2(n,m,kk,&pr,der,hes,
      nq,x,w,g,g1,g2);
    /* don't update hessian if pr<=0 */
    if(pr<=0.) { pr=1.e-10; /*printf("ir=%d, pr=%f\n", ir,pr);*/ }
    else
    { /* Hessian matrix plus gradient, set up for nrmin in matrix dd[][] */
      for(ip=1;ip<=np;ip++)
      { dd[ip-1][np]-= nfr*der[ip]/pr;
        for(jp=1;jp<=np;jp++)
        { tem=hes[ip][jp]-der[ip]*der[jp]/pr;
          dd[ip-1][jp-1]-= nfr*tem/pr;
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
  free(kk);
  free(der); free(hes[0]); free(hes);
  free(alp[0]); free(alp); free(b);
#else
  free(kk);
  free(der); free(hes[0]); free(hes);
  free(alp[0]); free(alp); free(b);
#endif
}

#ifdef R

void Rnrmlepln( int *nitem, int *ncateg, int *nrec, double *dataset,
               double *alphas, double *bvec, double *abound, double *bbound,
               double *nllkout, double *mleplnout, double *sevecout,
               double *invhesout, int *nq, int *mxiter,int *iconv, int *iprint) 
{
   double *param,*lb,*ub,fnval,**invhes,bdd; 
   double **dmatrix(int, int); 
   double ***gmatrix(int, int, int); 
   void gauher(double *, double *, int); 
   void plnnllk(int, double *, double *, double **, int, 
       int, int, int, int, double *, double *, double **, double *, 
       double ***, double ***, double ***); 
   int j,np,n2,i,ip,nn,jp; 
   double **alp,*b; 
   void nrmin(int,double *, void (*)(int, double *, double *, double   **, int, 
       int, int, int, int, double *, double *, double **, double *, 
       double ***, double ***, double ***),  
       double *lb, double *ub, int, double,int, double bdd,  
       int *, double *, double **, 
       int, int, int, int, double *, double *, double **, double *, 
       double ***, double ***, double ***); 
   double *x,*w; 
   int n;
   double **dat,*fr; 
   double ***g, ***g1, ***g2; 

   n=*nitem;

   /*nq=48; */
   x=(double *) calloc((*nq+1), sizeof(double)); 
   w=(double *) calloc((*nq+1), sizeof(double)); 
   gauher(x,w,*nq); 

   for (j=1;j<=*nq;j++) x[j]*=M_SQRT2; 
   for (j=1;j<=*nq;j++) w[j]/=SQRTPI;
    
   /*for (j=1;j<=nq;j++) Rprintf("%3d %14.6e %14.6e\n",j,x[j],w[j]);*/
  
   /*printf("%d %d %d %f\n", n, m, nrec,(*mxGetPr(prhs[1])));*/
 
   /* convert to matrices in C (i.e. row/column transpose) */
   dat=dmatrix(*nrec,*nitem); 
   fr=(double *) calloc((*nrec), sizeof(double)); 
   /* nn = total of fr[] */
   for(i=0,nn=0;i<*nrec;i++) 
   { for(j=0;j<*nitem;j++) dat[i][j] = *(dataset + (j*(*nrec)+i));  
     fr[i]=*(dataset + (i+(*nitem)*(*nrec))); nn+=fr[i];
     /*if(i<10) printf("%f %f\n", dat[i][n-1], fr[i]);*/
   } 

   for(i=1,n2=1;i<=*nitem;i++) n2*=(*ncateg);
   if (*iprint==1)
   {
     Rprintf("\nn=%d, #categ=%d, nrec=%d\n", *nitem,*ncateg,*nrec);
   }
   np=(*ncateg)*(*nitem); 

   lb=(double *) calloc(np, sizeof(double)); 
   ub=(double *) calloc(np, sizeof(double)); 
   for(ip=0;ip<np;ip++) { lb[ip]=*(abound+0); ub[ip]=*(abound+1); } 

  /* cannot put lower bound of 0 on betas unless
     all items oriented positively */
   for(ip=np-n;ip<np;ip++) { lb[ip]=*(bbound+0); ub[ip]=*(bbound+1);}
   ip=0; 

   /* alp[i][1],...,alp[i][m-1] should be in decreasing order */

   alp=dmatrix((*nitem)+1,(*ncateg)); 
   param=(double *) calloc(np, sizeof(double)); 
   for(i=0;i<(*nitem);i++)
   { for(j=0;j<((*ncateg)-1);j++)  
     { alp[i+1][j+1]=*(alphas + ((*ncateg)-1)*i+j); param[ip]=alp[i+1][j+1]; ip++;
       /*if(*iprint==1)
       {
         Rprintf("alpha(%d,%d)=%f\n",i,j,alp[i+1][j+1]);
       }*/
     } 
   } 

   b=(double *) calloc(((*nitem)+1), sizeof(double));  
   for(i=0;i<(*nitem);i++)   
   { param[ip]=*(bvec + i); b[i]=*(bvec + i); ip++;  
     /*if(*iprint==1)
     {
       Rprintf("beta(%d)=%f\n",i,b[i]);
     }*/
   }  

   /*for(j=1;j<=n;j++) printf("%d ", dat[nrec][j]);*/

  /*printf("First row: %d %d %d %d %d %d\n", dat[0][0], dat[0][1],
      dat[0][2], dat[0][3], dat[0][4], fr[0]);
    printf("\nn=%d, #categ=%d, nrec=%d\n", n,m,nrec);
    printf("nn=%d\n", nn);
  */

   invhes=dmatrix(np,np);  
   bdd=5.;  /* bound the difference between iterations */

   g=gmatrix((*nitem)+1,(*ncateg)+1,*nq+1);  
   g1=gmatrix((*nitem)+1,(*ncateg)+1,*nq+1);  
   g2=gmatrix((*nitem)+1,(*ncateg)+1,*nq+1);  

   nrmin(np,param,plnnllk,lb,ub,*mxiter,1.e-4,*iprint,bdd,iconv,&fnval,invhes, 
     *nitem, *ncateg, *nq, *nrec, x, w, dat, fr, g, g1, g2);  

  /* below used for testing only */
  /* printf("negative log-likelihood = %f\n", fnval);
  printf("Parameter estimates alp[i][j], j=1..m-1, i=1..n; b[i]\n");
  for(ip=0;ip<np;ip++) 
  { printf("%f ", param[ip]); if(ip%8==7) printf("\n"); }
  printf("\n");
  printf("SEs alp[i][j], j=1..m-1, i=1..n; b[i]\n");

  for(ip=0;ip<np;ip++) 
  { printf("%f ", sqrt(invhes[ip][ip])); if(ip%8==7) printf("\n"); }
  printf("\n"); */
 
  /*printf("inverse Hessian\n");
  for(ip=0;ip<np;ip++) 
  { for(jp=0;jp<np;jp++) printf("%.10f ", invhes[ip][jp]); printf("\n"); }
  */

  /* Creation of output */

   *nllkout = fnval;

   for(ip=0;ip<np;ip++) *(mleplnout + ip)=param[ip];  

   for(ip=0;ip<np;ip++) *(sevecout + ip)=sqrt(invhes[ip][ip]);  

   for(ip=0;ip<np;ip++)   
   { for(jp=0;jp<np;jp++) *(invhesout + (jp+ip*np))=invhes[ip][jp];}  

   free(invhes[0]); free(invhes);  
   free(alp[0]); free(alp);  
   free(b);  
   free(dat[0]); free(dat);  
   free(fr); free(param);  
   free(w); free(x);  
   free(ub); free(lb);  
   free(g[0][0]);  
   free(g[0]); free(g);  
   free(g1[0][0]);  
   free(g1[0]); free(g1);  
   free(g2[0][0]);  
   free(g2[0]); free(g2);  
}

#endif

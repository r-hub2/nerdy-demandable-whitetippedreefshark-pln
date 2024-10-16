/* June 2008, change to pointers and matlab interface, like nrmlepln.c */
#ifndef R
#include <malloc.h>
#else
#include <R.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* bivariate composite log-likelihood for polytomous logit-normit model
   (graded logistic)  using Newton Raphson 
   gcc -DDATA -o nrbcpln nrbcpln.c nrminbcl.c amatrix.c ssgauher.c -lm
*/

#define SQRTPI 1.772453850905516
#ifdef DATA
int main(int argc, char *argv[])
{ double *param,*lb,*ub,fnval,**invhes,bdd;
  double **dmatrix(int, int);
  double ***gmatrix(int, int, int);
  void gauher(double *, double *, int);
  void plnbclk(int, double *, double *, double **, int,
      int, int, int, int, double *, double *, double **, double *,
      double ***, double ***, double ***, double *pp);
  int j,np,i,ip,nn,jp,nm;
  int mxiter,iprint,iconv;
  double **alp,*b,tmp;
  FILE *in;
  void nrminbcl(int,double *, 
      void (*)(int, double *, double *, double **, int,
        int, int, int, int, double *, double *, double **, double *,
        double ***, double ***, double ***, double *pp), 
      double *lb, double *ub, int, double,int, double bdd, 
      int *, double *, double **,
      int, int, int, int, double *, double *, double **, double *,
      double ***, double ***, double ***, double *pp);
  double *x,*w;
  int nq,n,nrec,m;
  double **dat,*fr;
  double ***g, ***g1, ***g2;
  int ig,nb;
  void summ2frbv(int, int, int, int, double **dat, double *fr, double *, int *);
  double *pp;

  setbuf(stdout,NULL);
  if(argc==1) { printf("%s parfile [nq] < datafile", argv[0]); exit(0); }
  if(argc>=3) nq=atoi(argv[2]); else nq=48;
  if(nq>48) nq=48;
  printf("\nnq=%d\n",nq);
  x=(double *) malloc((nq+1) * sizeof(double));
  w=(double *) malloc((nq+1) * sizeof(double));
  gauher(x,w,nq);m
  for (j=1;j<=nq;j++) x[j]*=M_SQRT2;
  for (j=1;j<=nq;j++) w[j]/=SQRTPI;
  /*for (j=1;j<=nq;j++) printf("%3d %14.6e %14.6e\n",j,x[j],w[j]);*/
  
  in=fopen(argv[1],"r");
  mxiter=20;
  /*mxiter=1;*/ /* for testing */
  fscanf(in, "%d %d %d", &n,&m,&nrec); /* number of items and categories */
  printf("\nn=%d, #categ=%d, nrec=%d\n", n,m,nrec);
  np=m*n;

  ub=(double *) malloc(np * sizeof(double));
  lb=(double *) malloc(np * sizeof(double));
  for(ip=0;ip<np;ip++) { lb[ip]=-10; ub[ip]=10; }
  /*for(ip=0;ip<np;ip++) { lb[ip]=-7; ub[ip]=7; }*/
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

  /* bivariate composite: number of moments is n*(n-1)*m*m/2 */
  nm=(n*(n-1)*m*m)/2;
  printf("n=%d, m=%d, nm=%d\n", n,m,nm);
  pp=(double *) malloc((nm+2) * sizeof(double));

  dat=dmatrix(nrec,n);
  fr=(double *) malloc((nrec) * sizeof(double));
  for(i=0,nn=0;i<nrec;i++)
  { for(j=0;j<n;j++) scanf("%lf", &dat[i][j]); 
    scanf("%lf", &fr[i]); nn+=fr[i]; 
    /*for(j=0;j<n;j++) printf("%d ", dat[i][j]); 
      printf(": %d\n", fr[i]); */
  }
  /*for(j=0;j<n;j++) printf("%d ", dat[nrec][j]); */
  printf("nn=%d\n", nn);

  /* convert data to bivariate margins */
  /* order of bivariate margins is 12 13 23 14 24 34 etc */
  summ2frbv(n,m,nn,nrec,dat,fr,pp,&nb); 
  /* print out column vector of sample moments */
  printf("biv props as a vector, dim=%d\n", nb);
  for(ig=1;ig<=nb;ig++)
  { printf("%.3f ", pp[ig]); if(ig%10==0) printf("\n"); }
  printf("\n");

  invhes=dmatrix(np,np);
  iprint=1;
  bdd=5.;  /* bound the difference between iterations */
  g=gmatrix(n+1,m+1,nq+1);
  g1=gmatrix(n+1,m+1,nq+1);
  g2=gmatrix(n+1,m+1,nq+1);

  nrminbcl(np,param,plnbclk,lb,ub,mxiter,1.e-4,iprint,bdd,&iconv,&fnval,invhes,
    n, m, nq, nrec, x, w, dat, fr, g, g1, g2, pp);
  printf("bivariate composite negative log-likelihood = %f\n", fnval);
  printf("Parameter estimates alp[i][j], j=1..m-1, i=1..n; b[i]\n");
  for(ip=0;ip<np;ip++) 
  { printf("%f ", param[ip]); if(ip%8==7) printf("\n"); }
  printf("\n");
  /* replace by Godambe matrix later ? */

  free(invhes[0]); free(invhes);
  free(alp[0]); free(alp);
  free(b);
  free(dat[0]); free(dat);
  free(fr);
  free(param);
  free(g[0][0]); free(g[0]); free(g);
  free(g1[0][0]); free(g1[0]); free(g1);
  free(g2[0][0]); free(g2[0]); free(g2);
  free(w); free(x); free(pp);
  free(ub); free(lb);
  exit(0);
}
#endif

/* bivariate composite negative log-likelihood */
/* extra argument pp  that is not in plnnllk, so version of nrmin is
  different? */
void plnbclk(int np, double *param, double *obj, double **dd, int iprint,
  int n,int m, int nq, int nrec, double *x, double *w, double **dat, double *fr,
    double ***g, double ***g1, double ***g2, double *pp)
{ int i,j,ip,iq,jp;
  /*void plgnghbc(int [], int [], double *, double [], double [][NP]);*/
  void plgnghbc(int [], int[], double *, double *, double **,
    int,int, int, double *, double *,
    double ***, double ***, double ***);
  double **dmatrix(int, int);
  double tem,tem2,xx,bb,pr;
  int kk[3],ii[3];
  double *der,**hes;
  double **alp,*b;
  int j1,j2,k1,k2,ig;
  double nfr;

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
  der=(double *) malloc((np+1) * sizeof(double));
#else
  der=(double *) calloc((np+1), sizeof(double));
#endif
  hes=dmatrix(np+1,np+1);


  /* replace nllk by bivariate composite negative log-likelihood */
  ig=0;
  *obj=0;
  for(j2=2;j2<=n;j2++)
  { ii[2]=j2;
    for(j1=1;j1<j2;j1++)
    { /* margin (j1,j2) */
      ii[1]=j1;
      for(k1=0;k1<m;k1++) 
      { for(k2=0;k2<m;k2++) 
        { ig++; 
          /* from bivariate marginal probability */
          kk[1]=k1; kk[2]=k2; 
          nfr=pp[ig];
          if(nfr>0)
          { plgnghbc(ii,kk,&pr,der,hes,  n,m,nq,x,w,g,g1,g2);
            /* don't update hessian if pr<=0 */
            if(pr<=0.) pr=1.e-10;
            else
            { /* Hessian plus gradient, set up for nrmin in matrix dd[][] */
              for(ip=1;ip<=np;ip++)
              { dd[ip-1][np]-= nfr*der[ip]/pr;
                for(jp=1;jp<=np;jp++)
                { tem=hes[ip][jp]-der[ip]*der[jp]/pr;
                  dd[ip-1][jp-1]-= nfr*tem/pr;
                }
              }
            }
            *obj-=nfr*log(pr);
          }
        }
      }
    }
  }
  /* end of loop *obj is the negative log biv comp likelihood */
  /* for testing print out gradient */
#ifndef R
  if(iprint==1)
  { printf("\ngradient\n");
    for(ip=0;ip<np;ip++) printf("%f ", dd[ip][np]); printf("\n");
    printf("\n");
  }
#else
  if(iprint==1)
  { Rprintf("\ngradient\n");
    for(ip=0;ip<np;ip++) Rprintf("%f ", dd[ip][np]); Rprintf("\n");
    Rprintf("\n");
  }  
#endif
    
#ifndef R
  free(der); free(hes[0]); free(hes);
  free(alp[0]); free(alp); free(b);
#else
  free(der); free(hes[0]); free(hes);
  free(alp[0]); free(alp); free(b);   
#endif
}

/* summarize data into bivariate margins 
   categories assumed to be 0 1 ... m-1,
   get all 2D margins including category 0 */
void summ2frbv(int n, int m, int nn, int nrec, double **dat, double *fr, 
   double *pp, int *nb)
{ int i,j1,j2,k1,k2,ig;
  int **imatrix(int, int);
  int **s2;
  /*int s2[M][M];*/

  s2=imatrix(m,m);
  ig=0;

  /* bivariate margins */
  for(j2=1;j2<n;j2++)
  { for(j1=0;j1<j2;j1++)
    { for(k1=0;k1<m;k1++) 
      { for(k2=0;k2<m;k2++) s2[k1][k2]=0; }
      for(i=0;i<nrec;i++) 
      { k1=dat[i][j1]; k2=dat[i][j2]; s2[k1][k2]+=fr[i]; }

      for(k1=0;k1<m;k1++) 
      { for(k2=0;k2<m;k2++) 
        { ig++; pp[ig]=((double)s2[k1][k2]); }
      }
    }
  }
  *nb=ig;
#ifndef R
  free(s2[0]); free(s2);
#else
  free(s2[0]); free(s2);
#endif
}

/* set up calculation of first and second order derivatives
  for a term in bivariate marginal 
  margin ii[]=(i1,i2), categories kk[]=(k1,k2)
  this replaces plgngh2 of nllk
*/
void plgnghbc(int ii[], int kk[], double *pr, double *der, double **hes,
    int n, int m, int nq, double *x, double *w,
    double ***g, double ***g1, double ***g2)
{ int k,np,id,jd,ii1,ib,ik,ik1,ic,np1;
  int jj1,jk,jk1,k1,k2,jb,jc;
  int iz,ii0,jj0;
  void pghderbc(int ii[], int kk[], double [], double [][7],
      int nq, double *x, double *w,
      double ***g, double ***g1, double ***g2);
  double fn[7],fn2[7][7];

  np=m*n; np1=np-n;
  /* loop for reindexing partial derivatives */
  for(id=0;id<=np;id++) der[id]=0.;
  for(id=0;id<=np;id++) { for(jd=0;jd<=np;jd++) hes[id][jd]=0.; }
  /* do this sum in pghderbc */
  /*for(iq=1;iq<=nq;iq++)*/
  pghderbc(ii,kk,fn,fn2,  nq,x,w,g,g1,g2); 
  *pr=fn[0];
  /* need reindexing here */
  for(iz=1;iz<=2;iz++)
  { k=kk[iz]; ii0=2*iz-1; ii1=ii0+1; ib=4+iz; 
    ik=(ii[iz]-1)*(m-1)+k; ik1=ik+1; ic=np1+ii[iz];
    if(k==0)
    { der[ik1]=fn[ii1]; 
      hes[ik1][ik1]=fn2[ii1][ii1];
      hes[ik1][ic]=fn2[ii1][ib];
      /* add the transpose  */
      hes[ic][ik1]=hes[ik1][ic];
    }
    else if(k==m-1)
    { der[ik]=fn[ii0]; 
      hes[ik][ik]=fn2[ii0][ii0];
      hes[ik][ic]=fn2[ii0][ib];
      hes[ic][ik]=hes[ik][ic];
    }
    else
    { der[ik]=fn[ii0]; der[ik1]=fn[ii1]; 
      hes[ik][ik]=fn2[ii0][ii0]; hes[ik1][ik1]=fn2[ii1][ii1];
      hes[ik][ic]=fn2[ii0][ib]; hes[ik1][ic]=fn2[ii1][ib];
      hes[ic][ik]=hes[ik][ic]; hes[ic][ik1]=hes[ik1][ic];
    }
    der[ic]=fn[ib];
    hes[ic][ic]=fn2[ib][ib];
  }
  /* second order mixed derivatives  only iz=1, jz=2 (no loop needed) */
  k1=kk[1]; ii0=1; ii1=ii0+1; ib=5; 
  ik=(ii[1]-1)*(m-1)+k1; ik1=ik+1; ic=np1+ii[1];
  if(k1==0) ik=0;
  else if(k1==m-1) ik1=0;
  k2=kk[2]; jj0=3; jj1=jj0+1; jb=6; 
  jk=(ii[2]-1)*(m-1)+k2; jk1=jk+1; jc=np1+ii[2];
  if(k2==0) jk=0;
  else if(k2==m-1) jk1=0;
  /*printf("i=%d, j=%d, ik=%d, ik1=%d, jk=%d, jk1=%d, ii=%d, jj=%d\n",
    i,j,ik,ik1,jk,jk1,ii,jj); */
  /* cases k1=0,k2=0 etc 9 cases in all : instead
     if(k1=0) set ik=0 
     if(k1=m) set ik1=0  */
  hes[ik][jk] = fn2[ii0][jj0];
  hes[ik][jk1] = fn2[ii0][jj1];
  hes[ik1][jk] = fn2[ii1][jj0];
  hes[ik1][jk1] = fn2[ii1][jj1];
  hes[ik][jc] = fn2[ii0][jb];
  hes[ik1][jc] = fn2[ii1][jb];
  hes[jk][ic] = fn2[jj0][ib];
  hes[jk1][ic] = fn2[jj1][ib];
  hes[ic][jc] = fn2[ib][jb];
  /* transposes */
  hes[jk][ik]=hes[ik][jk];
  hes[jk1][ik]=hes[ik][jk1];
  hes[jk][ik1]=hes[ik1][jk];
  hes[jk1][ik1]=hes[ik1][jk1];
  hes[jc][ik]=hes[ik][jc];
  hes[jc][ik1]=hes[ik1][jc];
  hes[ic][jk]=hes[jk][ic];
  hes[ic][jk1]=hes[jk1][ic];
  hes[jc][ic]=hes[ic][jc];
}

void pghderbc(int ii[], int kk[], double sfn[], double sfn2[][7],
    int nq, double *x, double *w,
  double ***g, double ***g1, double ***g2)
{ double diff,tem,diff1,diff2;
  int iq,i,j,k,np,id,jd,np1,k1,k2,iiz,jjz,ib,jb,iz;
  double xx,der1,der2;
  double fn[7],fn2[7][7];

  /*np=3*n; np1=2*n;*/
  np=6; np1=4;
  for(id=0;id<=np;id++) sfn[id]=0.;
  for(id=1;id<=np;id++) 
  { for(jd=1;jd<=np;jd++) sfn2[id][jd]=0.; }

  for(iq=1;iq<=nq;iq++)
  { xx=x[iq];
    fn[0]=1.;
    for(id=1;id<=np;id++) fn[id]=0.;
    for(id=1;id<=np;id++) 
    { for(jd=1;jd<=np;jd++) fn2[id][jd]=0.; }
    /* first order derivatives and second order same i */
    for(iz=1,tem=1.;iz<=2;iz++)
    { k=kk[iz]; i=ii[iz]; /* i is index of margin */
      diff=g[i][k][iq]-g[i][k+1][iq];
      tem*=diff;
      iiz=iz*2; ib=np1+iz;
      if(diff>0.)
      { fn[iiz-1]= g1[i][k][iq]/diff;
        fn[iiz]= -g1[i][k+1][iq]/diff;
        fn2[iiz-1][iiz-1]= g2[i][k][iq]/diff;
        fn2[iiz][iiz]= -g2[i][k+1][iq]/diff;
        fn[ib]=xx*(g1[i][k][iq]-g1[i][k+1][iq])/diff;
        fn2[ib][ib]=xx*xx*(g2[i][k][iq]-g2[i][k+1][iq])/diff;
        fn2[iiz-1][ib]=xx*g2[i][k][iq]/diff;
        fn2[iiz][ib]=-xx*g2[i][k+1][iq]/diff;
      }
    }

    /* second order mixed derivatives  only iz=1, jz=2 (no loop needed) */
    iiz=2; ib=np1+1;
    k1=kk[1]; i=ii[1];
    diff1=g[i][k1][iq]-g[i][k1+1][iq];
    der1=g1[i][k1][iq]-g1[i][k1+1][iq];
    if(diff1>0.) 
    { jjz=4; jb=np1+2;
      k2=kk[2]; j=ii[2];
      diff2=g[j][k2][iq]-g[j][k2+1][iq];
      if(diff2>0.)
      { diff=diff1*diff2;
        der2=g1[j][k2][iq]-g1[j][k2+1][iq];
        fn2[iiz-1][jjz-1]= g1[i][k1][iq]*g1[j][k2][iq]/diff;
        fn2[iiz-1][jjz]= -g1[i][k1][iq]*g1[j][k2+1][iq]/diff;
        fn2[iiz][jjz-1]= -g1[i][k1+1][iq]*g1[j][k2][iq]/diff;
        fn2[iiz][jjz]= g1[i][k1+1][iq]*g1[j][k2+1][iq]/diff;
        fn2[iiz-1][jb]= g1[i][k1][iq]*xx*der2/diff;
        fn2[iiz][jb]= -g1[i][k1+1][iq]*xx*der2/diff;
        fn2[jjz-1][ib]= g1[j][k2][iq]*xx*der1/diff;
        fn2[jjz][ib]= -g1[j][k2+1][iq]*xx*der1/diff;
        fn2[ib][jb]=xx*xx*der1*der2/diff;
      }
    }
    tem*=w[iq];
    for(id=0;id<=np;id++) { fn[id]*=tem; sfn[id]+=fn[id]; }
    for(id=1;id<=np;id++) 
    { for(jd=id;jd<=np;jd++) 
      { fn2[id][jd]*=tem; sfn2[id][jd]+=fn2[id][jd]; }
    }
  }
  /* w[i] now divides by SQRTPI */
  /* for(id=0;id<=np;id++) sfn[id]/=SQRTPI;
  for(id=1;id<=np;id++) 
  { for(jd=id;jd<=np;jd++) sfn2[id][jd]/=SQRTPI; } */
  /* loop and print out only for testing */
  /*
  printf("\n");
  for(id=1;id<np;id++) 
  { for(jd=id+1;jd<=np;jd++) sfn2[jd][id]=sfn2[id][jd]; }
  for(id=1;id<=np;id++) 
  { for(jd=1;jd<=np;jd++) printf("%f ", sfn2[id][jd]); printf("\n"); }
  */
}

#ifdef R
void Rnrbcpln (int *nitem, int *ncateg, int *nrec, double *dataset, double *alphas,
               double *bvec, double *abound, double *bbound, double *nbclpln,
               double *bclpln, int *nq, int *mxiter, int *iconv, int *iprint)
{
   double *param,*lb,*ub,fnval,**invhes,bdd;
   double **dmatrix(int, int);
   double ***gmatrix(int, int, int);
   void gauher(double *, double *, int);
   void plnbclk(int, double *, double *, double **, int,
      int, int, int, int, double *, double *, double **, double *,
      double ***, double ***, double ***, double *pp);
   int j,np,i,ip,nn,nm;

  double **alp,*b;
   void nrminbcl(int,double *, void (*)(int, double *, double *, double **, int,
      int, int, int, int, double *, double *, double **, double *,
      double ***, double ***, double ***, double *), 
      double *lb, double *ub, int, double,int, double bdd, 
      int *, double *, double **,
      int, int, int, int, double *, double *, double **, double *,
      double ***, double ***, double ***, double *);
   double *x,*w;
   double **dat,*fr;
   double ***g, ***g1, ***g2;
   int nb;
   void summ2frbv(int, int, int, int, double **dat, double *fr, double *, int *);
   double *pp;
   
   x=(double *) calloc((*nq+1), sizeof(double));
   w=(double *) calloc((*nq+1), sizeof(double));
   gauher(x,w,*nq); 
   for (j=1;j<=*nq;j++) x[j]*=M_SQRT2; 
   for (j=1;j<=*nq;j++) w[j]/=SQRTPI; 

    /*for (j=1;j<=nq;j++) printf("%3d %14.6e %14.6e\n",j,x[j],w[j]);*/
 
  /* printf("%d %d %d %f\n", n, m, nrec,(*mxGetPr(prhs[1]))); */

  /* convert to matrices in C (i.e. row/column transpose) */
   dat=dmatrix(*nrec,*nitem);
   fr=(double *) calloc((*nrec), sizeof(double));
  /* nn = total of fr[] */
   for(i=0,nn=0;i<(*nrec);i++) 
   { for(j=0;j<(*nitem);j++) dat[i][j] = *(dataset + (j*(*nrec)+i)); 
     fr[i]= *(dataset + (i+(*nitem)*(*nrec))); nn+=fr[i];
     /* if(i<10) printf("%f %f\n", dat[i][n-1], fr[i]); */
   }

   /*for(i=1,n2=1;i<=n;i++) n2*=m;*/
   np=(*ncateg)*(*nitem); 
   nm=((*nitem)*((*nitem)-1)*(*ncateg)*(*ncateg))/2; 
   pp=(double *) calloc((nm+2), sizeof(double));

   lb=(double *) calloc(np, sizeof(double));
   ub=(double *) calloc(np, sizeof(double));
   for(ip=0;ip<np;ip++) { lb[ip]=*(abound+0); ub[ip]=*(abound+1); }
    
   if (*iprint==1)
   {
     Rprintf("\nn=%d, #categ=%d, nrec=%d\n", *nitem,*ncateg,*nrec);
   }
   
   /* cannot put lower bound of 0 on betas unless
     all items oriented positively */
   for(ip=np-(*nitem);ip<np;ip++) { lb[ip]=*(bbound+0); ub[ip]=*(bbound+1); }
    /*for(ip=np-(*nitem);ip<np;ip++) { lb[ip]=0; }*/

  ip=0;

  /* alp[i][1],...,alp[i][m-1] should be in decreasing order */
  /*alphas */

   alp=dmatrix((*nitem)+1,(*ncateg));
   param=(double *) calloc(np, sizeof(double));
   for(i=0;i<(*nitem);i++)
   { for(j=0;j<((*ncateg)-1);j++) 
     { alp[i+1][j+1]=*(alphas + ((*ncateg-1)*i+j)); param[ip]=alp[i+1][j+1]; ip++;
       /*printf("alpha(%d,%d)=%f\n",i,j,alp[i+1][j+1]);*/
     }
   }

  /*bvec*/
   b=(double *) calloc(((*nitem)+1), sizeof(double)); 
   for(i=0;i<(*nitem);i++)  
   { param[ip]=*(bvec + i); b[i+1]=*(bvec + i); ip++; 
    /*printf("beta(%d)=%f\n",i,b[i]);*/
   }

  /*for(j=0;j<n;j++) printf("%d ", dat[nrec][j]);*/

  /* convert data to bivariate margins */
  /* order of bivariate margins is 12 13 23 14 24 34 etc */
   summ2frbv(*nitem,*ncateg,nn,*nrec,dat,fr,pp,&nb); 
  /* testing: print out column vector of sample moments */
  /*printf("biv props as a vector, dim=%d\n", nb);
  for(ig=1;ig<=nb;ig++)
  { printf("%.3f ", pp[ig]); if(ig%10==0) printf("\n"); }
  printf("\n"); */

  invhes=dmatrix(np,np);

  bdd=5.;  /* bound the difference between iterations */

  g=gmatrix((*nitem)+1,(*ncateg)+1,*nq+1);
  g1=gmatrix((*nitem)+1,(*ncateg)+1,*nq+1);
  g2=gmatrix((*nitem)+1,(*ncateg)+1,*nq+1);

  nrminbcl(np,param,plnbclk,lb,ub,*mxiter,1.e-4,*iprint,bdd,iconv,&fnval,invhes,
    *nitem, *ncateg, *nq, *nrec, x, w, dat, fr, g, g1, g2, pp);

  /* below used for testing only */
  /* printf("negative BCL = %f\n", fnval);
  printf("Parameter estimates alp[i][j], j=1..m-1, i=1..n; b[i]\n");
  for(ip=0;ip<np;ip++) 
  { printf("%f ", param[ip]); if(ip%8==7) printf("\n"); }
  printf("\n");
  */

	*nbclpln = fnval;

  for(ip=0;ip<np;ip++) *(bclpln + ip)=param[ip];

/* could be replaced with Godambe matrix later */

  free(invhes[0]); free(invhes);
  free(alp[0]); free(alp);
  free(b);
  free(dat[0]); free(dat);
  free(fr); free(param);
  free(w); free(x); free(pp);
  free(ub); free(lb);
  free(g[0][0]);
  free(g[0]); free(g);
  free(g1[0][0]);
  free(g1[0]); free(g1);
  free(g2[0][0]);
  free(g2[0]); free(g2);
}

#endif



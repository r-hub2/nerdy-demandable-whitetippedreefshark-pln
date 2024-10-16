#ifndef R
#include <malloc.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 21
#define M 6
/* gauss-hermite integ for ordinal logit-normal model */
/* compute derivatives and Hessian matrix of pr wrt parameters */

/* gcc -DMAIN -o polyder2 polyder2.c ssgauher.c -lm */
#define SQRTPI 1.772453850905516
#define NP 101
#define NQ 49
#ifdef MAIN

double *x,*w;
int nq;
double ***g,***g1,***g2;

main(int argc, char *argv[])
{ int i,j,np,n,n2,m,id,jd;
  double **dmatrix(int, int);
  double ***gmatrix(int, int, int);
  void gauher(double *, double *, int);
  double alp[N][M],b[N];
  void plgngh2(int, int, int[], double *, double *, double **,
    int, double *, double *, double ***, double ***, double ***);
  double pr,*der,**hes,tem,xx,bb,tem2;
  int kk[N],iq;

  setbuf(stdout,NULL);
  if(argc>=2) nq=atoi(argv[1]); else nq=48;  /* nq <64? */
  if(nq>48) nq=48;
  printf("\nnq=%d\n",nq);
  x=(double *) malloc((nq+1) * sizeof(double));
  w=(double *) malloc((nq+1) * sizeof(double));
  gauher(x,w,nq);
  for (j=1;j<=nq;j++) x[j]*=M_SQRT2;
  for (j=1;j<=nq;j++) w[j]/=SQRTPI;
  /*for (j=1;j<=nq;j++) printf("%3d %14.6e %14.6e\n",j,x[j],w[j]);*/

  scanf("%d %d", &n,&m); /* number of items and categories */
  for(i=1,n2=1;i<=n;i++) n2*=m;
  printf("\nn=%d, m^n=%d\n", n,n2);
  np=m*n;
  /* alp[i][1],...,alp[i][m-1] should be in decreasing order */
  for(i=1;i<=n;i++)
  { for(j=1;j<m;j++) scanf("%lf", &alp[i][j]); }
  for(i=1;i<=n;i++) scanf("%lf", &b[i]);
  /* globally computes terms of integrands (and their derivatives)
     for each x in GH quad pts */

  g=gmatrix(N,M,NQ);
  g1=gmatrix(N,M,NQ);
  g2=gmatrix(N,M,NQ);
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

  
  /* kk[] is vector of categories */
  der=(double *) malloc(NP * sizeof(double));
  hes=dmatrix(NP,NP);
  for(i=1;i<=n;i++) scanf("%d", &kk[i]);
  while(kk[1]>=0)
  { for(i=1;i<=n;i++) printf("k[%d]=%d; ", i,kk[i]); printf("\n");
    plgngh2(n,m,kk,&pr,der,hes,nq,x,w,g,g1,g2);
    printf(" pr=%f\n", pr);
    for(id=1;id<=np;id++) printf("%10.6f ", der[id]); printf("\n"); 
    for(id=1;id<=np;id++) 
    { for(jd=1;jd<=np;jd++) printf("%f ", hes[id][jd]); printf("\n"); }
    for(i=1;i<=n;i++) scanf("%d", &kk[i]);
  }
  free(w); free(x);
  free(der);
  free(hes[0]); free(hes);
  free(g[0][0]); free(g[0]), free(g);
  free(g1[0][0]); free(g1[0]), free(g1);
  free(g2[0][0]); free(g2[0]), free(g2);
  printf("===============================================\n");
  exit(0);
}
#endif

/* more general version with 2nd order derivs in hes[][] */
void plgngh2(int n, int m, int kk[], double *pr, double *der, double **hes,
    int nq, double *x, double *w,
    double ***g, double ***g1, double ***g2)
{ int i,j,k,np,id,jd,ii,ii1,ib,ik,ik1,ic,np1;
  int jj,jj1,jk,jk1,k1,k2,jb,jc;
  void pghder2(int n, int kk[], double [], double [][3*N],
    int, double *, double *,
    double ***, double ***, double ***);
  double fn[3*N],fn2[3*N][3*N];

  np=m*n; np1=np-n;
  /* loop for reindexing partial derivatives  */
  for(id=0;id<=np;id++) der[id]=0;
  for(id=0;id<=np;id++) { for(jd=0;jd<=np;jd++) hes[id][jd]=0; }
  /* do this sum in pghder2 */
  /*for(iq=1;iq<=nq;iq++)*/
  pghder2(n,kk,fn,fn2,nq,x,w,g,g1,g2); 
  *pr=fn[0];
  /* need reindexing here 1 -> k1,  2 -> k1+1, 3 -> k2 etc */
  for(i=1;i<=n;i++)
  { k=kk[i]; ii=2*i-1; ii1=ii+1; ib=2*n+i; 
    ik=(i-1)*(m-1)+k; ik1=ik+1; ic=np1+i;
    if(k==0)
    { der[ik1]=fn[ii1];
      hes[ik1][ik1]=fn2[ii1][ii1];
      hes[ik1][ic]=fn2[ii1][ib];
      /* add the transpose */
      hes[ic][ik1]=hes[ik1][ic];
    }
    else if(k==m-1)
    { der[ik]=fn[ii]; 
      hes[ik][ik]=fn2[ii][ii];
      hes[ik][ic]=fn2[ii][ib];
      hes[ic][ik]=hes[ik][ic];
    }
    else
    { der[ik]=fn[ii]; der[ik1]=fn[ii1]; 
      hes[ik][ik]=fn2[ii][ii]; hes[ik1][ik1]=fn2[ii1][ii1];
      hes[ik][ic]=fn2[ii][ib]; hes[ik1][ic]=fn2[ii1][ib];
      hes[ic][ik]=hes[ik][ic]; hes[ic][ik1]=hes[ik1][ic];
    }
    der[ic]=fn[ib];
    hes[ic][ic]=fn2[ib][ib];
  }
  /* mixed derivatives */
  for(i=1;i<n;i++)
  { k1=kk[i]; ii=2*i-1; ii1=ii+1; ib=2*n+i; 
    ik=(i-1)*(m-1)+k1; ik1=ik+1; ic=np1+i;
    if(k1==0) ik=0;
    else if(k1==m-1) ik1=0;
    for(j=i+1;j<=n;j++)
    { k2=kk[j]; jj=2*j-1; jj1=jj+1; jb=2*n+j; 
      jk=(j-1)*(m-1)+k2; jk1=jk+1; jc=np1+j;
      if(k2==0) jk=0;
      else if(k2==m-1) jk1=0;
      /*printf("i=%d, j=%d, ik=%d, ik1=%d, jk=%d, jk1=%d, ii=%d, jj=%d\n",
        i,j,ik,ik1,jk,jk1,ii,jj); */
      /* cases k1=0,k2=0 etc 9 cases in all : instead 
         if(k1=0) set ik=0 
         if(k1=m) set ik1=0 */
      hes[ik][jk] = fn2[ii][jj];
      hes[ik][jk1] = fn2[ii][jj1];
      hes[ik1][jk] = fn2[ii1][jj];
      hes[ik1][jk1] = fn2[ii1][jj1];
      hes[ik][jc] = fn2[ii][jb];
      hes[ik1][jc] = fn2[ii1][jb];
      hes[jk][ic] = fn2[jj][ib];
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
  }
}

void pghder2(int n, int kk[], double sfn[], double sfn2[][3*N],
  int nq, double *x, double *w,
  double ***g, double ***g1, double ***g2)
{ double diff,tem,diff1,diff2;
  int iq,i,j,k,np,id,jd,np1,k1,k2,ii,jj,ib,jb;
  double xx,der1,der2;
  double fn[3*N],fn2[3*N][3*N];
  /*extern double *x,*w;*/
  /*extern int nq;*/
  /*extern double ***g,***g1,***g2;*/

  np=3*n; np1=2*n;
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
    for(i=1,tem=1.;i<=n;i++)
    { k=kk[i];
      diff=g[i][k][iq]-g[i][k+1][iq];
      tem*=diff;
      ii=i*2; ib=np1+i;
      if(diff>0.)
      { fn[ii-1]= g1[i][k][iq]/diff;
        fn[ii]= -g1[i][k+1][iq]/diff;
        fn2[ii-1][ii-1]= g2[i][k][iq]/diff;
        fn2[ii][ii]= -g2[i][k+1][iq]/diff;
        fn[ib]=xx*(g1[i][k][iq]-g1[i][k+1][iq])/diff;
        fn2[ib][ib]=xx*xx*(g2[i][k][iq]-g2[i][k+1][iq])/diff;
        fn2[ii-1][ib]=xx*g2[i][k][iq]/diff;
        fn2[ii][ib]=-xx*g2[i][k+1][iq]/diff;
      }
    }

    /* second order mixed derivatives */
    for(i=1;i<n;i++)
    { ii=i*2; ib=np1+i;
      k1=kk[i];
      diff1=g[i][k1][iq]-g[i][k1+1][iq];
      der1=g1[i][k1][iq]-g1[i][k1+1][iq];
      if(diff1<=0.) continue;
      for(j=i+1;j<=n;j++)
      { jj=j*2; jb=np1+j;
        k2=kk[j];
        diff2=g[j][k2][iq]-g[j][k2+1][iq];
        if(diff2>0.)
        { diff=diff1*diff2;
          der2=g1[j][k2][iq]-g1[j][k2+1][iq];
          fn2[ii-1][jj-1]= g1[i][k1][iq]*g1[j][k2][iq]/diff;
          fn2[ii-1][jj]= -g1[i][k1][iq]*g1[j][k2+1][iq]/diff;
          fn2[ii][jj-1]= -g1[i][k1+1][iq]*g1[j][k2][iq]/diff;
          fn2[ii][jj]= g1[i][k1+1][iq]*g1[j][k2+1][iq]/diff;
          fn2[ii-1][jb]= g1[i][k1][iq]*xx*der2/diff;
          fn2[ii][jb]= -g1[i][k1+1][iq]*xx*der2/diff;
          fn2[jj-1][ib]= g1[j][k2][iq]*xx*der1/diff;
          fn2[jj][ib]= -g1[j][k2+1][iq]*xx*der1/diff;
          fn2[ib][jb]=xx*xx*der1*der2/diff;
        }
      }
    }
    /*fn[0]=tem;*/
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
  /* set up lower triangle?*/
  /*
  for(id=1;id<np;id++) 
  { for(jd=id+1;jd<=np;jd++) sfn2[jd][id]=sfn2[id][jd]; }
  for(id=1;id<=np;id++) 
  { for(jd=1;jd<=np;jd++) printf("%f ", sfn2[id][jd]); printf("\n"); }
  */
}


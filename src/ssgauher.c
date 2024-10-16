#ifndef R
#include <malloc.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* gauss-hermite quadrature points and weights 
   from p 35-36 of Stroud and Secrest, gaussian quadrature formulas 
   QA 299.4 G3 S7 1966 */
/* converted to orthonormal polys to get more accuracy for n>30 */
/* translated from Fortran to C by H Joe in April 2004 */
/* gcc -DMAIN -o ssgauher ssgauher.c -lm */

#ifdef MAINQ
/* sample main program */
main(int argc, char *argv[])
{ int n,i,np;
  double *x,*w,eps;
  double sum;
  /*void hermit(int, double *, double *, double);*/
  void gauher(double *, double *, int);

  np=50;
  x=(double *) malloc(np * sizeof(double));
  w=(double *) malloc(np * sizeof(double));
  
  eps=3.e-14;
  scanf("%d", &n);
  while(n>0)
  { /*hermit(n,x,w,eps);*/
    gauher(x,w,n);
    for(i=1;i<=n;i++)
    { printf("%5d %25.14f %25.14e\n", i,x[i],w[i]); }
    
    /* test function :  cos */
    for(i=1,sum=0.;i<=n;i++) sum+=w[i]*cos(x[i]);
    printf("gh:   %.10f\n",sum);
    printf("exact %.10f\n", 1.772453850905516*exp(-.25));
    scanf("%d", &n);
  }
}
#endif

#define EPS 3.0e-14
#define MAXIT 10
/*void hermit(int nn, double *x, double *a, double eps)*/
void gauher(double *x, double *a, int nn)
{ /* calculates the zeros x(i) of the nn-th order hermite polynomial 
     the largest zero will be stored in x(1). also calculates the 
     corresponding coefficients a(i) of the nn-th order gauss-hermite 
     quadrature formula of degree 2*nn-1 */
  int n1,i,n2,ni;
  double fn,s,xt,dpn,pn1;
  void hroot(double *, int, double *, double *, double);
  
  fn=(double)nn;
  n1=nn-1;
  n2=(nn+1)/2;
  s=pow(2.*fn+1.,.16667);
  for(i=1;i<=n2;i++)
  { if(i==1) xt=pow(s,3.)-1.85575/s; /* largest zero */
    else if(i==2) xt=xt-1.14*pow(fn,.426)/xt; /* second zero */
    else if(i==3) xt=1.86*xt-.86*x[1]; /* third zero */
    else if(i==4) xt=1.91*xt-.91*x[2]; /* fourth zero */
    else xt=2.*xt-x[i-2]; /* all other zeros */
    
    /*hroot(&xt,nn,&dpn,&pn1,eps);*/
    hroot(&xt,nn,&dpn,&pn1,EPS);
    x[i]=xt;
    a[i]=2./dpn/dpn;
    ni=nn-i+1;
    x[ni]=-xt;
    a[ni]=a[i];
  }
}

void hroot(double *x, int nn, double *dpn, double *pn1, double eps)
{ /* improves the approximate root x 
     in addition also obtain 
     dpn=derivative of h_n at x 
     pn1=value of h_{n-1} at x */
  int iter;
  double d,p,dp;
  void hrecur(double *,double *,double *,double,int);
  iter=0;
  while(iter<MAXIT)
  { iter++;
    hrecur(&p,&dp,pn1,*x,nn);
    d=p/dp;
    (*x)-=d;
    if(fabs(d)<=eps) { *dpn=dp; return; }
  }
#ifndef R
  printf("too many iterations in hroot\n");
#endif
}

/* recursion with orthonormal Hermite */
void hrecur(double *pn, double *dpn, double *pn1, double x, int nn)
{ int j;
  double p1,p,dp1,dp,fj,fj2,q,dq,fj1;
  p1=0.7511255444649425; /* pi^{-.25} */
  dp1=0.;
  p=x*1.062251932027197;  /* pi^{-.25}*sqrt(2) */
  dp=1.062251932027197;
  for(j=2;j<=nn;j++)
  { fj=j;
    fj1=sqrt(2./fj);
    fj2=sqrt((fj-1.)/fj);
    q=x*fj1*p-fj2*p1;
    dq=(x*dp+p)*fj1-fj2*dp1;
    p1=p; p=q; dp1=dp; dp=dq;
  }
  *pn=p; *dpn=dp; *pn1=p1;
}


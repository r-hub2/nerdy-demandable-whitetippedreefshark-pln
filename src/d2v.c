#ifndef R
#include <malloc.h>
#endif

#include <stdio.h>
#include <math.h>

/* decimal to category vector and vice versa */
/* gcc -DMAIN -o d2v d2v.c -lm */
#ifdef MAIN
main()
{ int d,dd,i,j,k,*jj,*a,ii;
  //void d2brev(int, int, int []);
  //int b2drev(int, int []);
  void d2v(int, int, int, int *);
  int v2d(int, int, int *);
  int m;

  m=10;
  jj=(int *) malloc(m * sizeof(int));
  a=(int *) malloc(m * sizeof(int));

  /* d = dimension, dd=k^d */
  scanf("%d %d %d", &k,&d,&i);
  while(k>0)
  { d2v(d,k,i,jj);
    printf("i=%d: ", i);
    for(j=1;j<=d;j++) printf(" %d", jj[j]); printf("\n");
    ii=v2d(d,k,jj);
    printf("ii=%d\n", ii);
    scanf("%d %d %d", &k,&d,&i);
  }
  free(jj); free(a);
}
#endif

/* returns k-category vector jj of dimension d  corr. to ii */
void d2v(int d, int k, int ii, int *jj)
{ int i,t;
  t=ii;
  for(i=d;i>=1;i--)
  { jj[i]=t%k; t=t/k; }
}

/* k-category vector jj to decimal, length d */
int v2d(int d, int k, int jj[])
{ int i,pw,s;
  pw=1;
  for(i=d,s=0;i>=1;i--)
  { s+=pw*jj[i]; pw*=k; }
  return s;
}



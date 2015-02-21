/*
 A re-implementation of CLINDE in C (originally in R: infer.R and pcor.R)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include "parse_option.h"
#include "tsv.h"
/* we include gsl-1.16 in CLINDE locally for convenience.
   It should be properly configured, built and installed.
 */
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cdf.h>

/************************************************************************/
/* mainly for debug */
/*
static void pr_ints(int n, int x[]) {
  int i=0;
  for(i=0; i<n; i++)
    printf(" %d", x[i]);
  printf("\n");
}
static void pr_int_arr(int n, int ns[], int* xs[]) {
  int i=0;
  printf("\n");
  for(i=0; i<n; i++) {
    printf("%d:", i);
    pr_ints(ns[i],xs[i]);
  }
}
static void pr_doubles(int n, double x[]) {
  int i=0;
  for(i=0; i<n; i++)
    printf(" %f", x[i]);
  printf("\n");
}
*/

/***********************************************************************/
static clock_t stop_watch_start = 0;
static void start_stop_watch() {
  printf("Start stop watch.\n");
  stop_watch_start = clock();
}
static void report_stop_watch() {
  printf("Stop watch elapsed time: %f sec\n",
	 (float)(clock() - stop_watch_start)/CLOCKS_PER_SEC);
}

/***********************************************************************/
static void* Malloc_buf(size_t s, char* purpose) {
  void* r = NULL;
  r = malloc(s);
  if(r == NULL) {
    fprintf(stderr, "Error: Cannot allocate %d bytes for %s.\n", s, purpose);
    exit(1);
  }
  return r;
}
#define Malloc(n,t,p) ((t *)Malloc_buf(sizeof(t)*(n), p))

/***********************************************************************/
typedef struct s_array2d {
  int rows,cols;
  double* v;
} array2d;

/* assume row major format */
static double aref(array2d* d, int i, int j) {return d->v[j+i*(d->cols)];}
static void aset(array2d* d, int i, int j, double val) {d->v[j+i*(d->cols)] = val;}

/* assume column major format */
static double caref(array2d* d, int i, int j) {return d->v[i+j*(d->rows)];}
static void caset(array2d* d, int i, int j, double val) {d->v[i+j*(d->rows)] = val;}
static double* get_col(array2d* d, int c) {return d->v+c*(d->rows);}
/****/
static array2d* new_array2d(int n_rows, int n_cols) {
  array2d* r=NULL;
  r = Malloc(1, array2d, "new_array2d()");
  r->rows = n_rows;
  r->cols = n_cols;
  r->v = Malloc(n_rows*n_cols, double, "new_array2d()");
  return r;
}
static void free_array2d(array2d* d) {
  if(!d) return;
  if(d->v) free(d->v);
  free(d);
}

static void print_array2d(FILE* f, array2d* d) {
  int i=0, j=0;
  for(i=0; i<(d->rows); i++) {
    for(j=0; j<(d->cols); j++) {
      fprintf(f, "%f\t", aref(d,i,j));
    }
    fprintf(f,"\n");
  }
}

static void print_carray2d(FILE* f, array2d* d) {
  int i=0, j=0;
  for(i=0; i<(d->rows); i++) {
    for(j=0; j<(d->cols); j++) {
      fprintf(f, "%f\t", caref(d,i,j));
    }
    fprintf(f,"\n");
  }
}

/***************************************************************/
typedef struct s_link {
  int from, to; /* from --> to */
  int delay;
  double score;
  double test_value;
} link;

/***************************************************************/
static double var(int n, double x[]) {
  /* variance of x, n-1 as denominator. */
  int i=0;
  double sx2=0, mx=0, r=0;
  for(i=0; i<n; i++) mx += x[i];
  mx /= n;
  for(i=0; i<n; i++) {
    r = (x[i] - mx);
    sx2 += r*r;
  }
  return sx2/(n-1);
}

static double cov(int n, double x[], double y[]) {
  /* covariance of x and y, n-1 as denominator. */
  int i=0;
  double sxy=0, sx=0, sy=0;
  for(i=0; i<n; i++) {
    sx += x[i];
    sy += y[i];
    sxy += x[i]*y[i];
  }
  return (sxy - sx*sy/n)/(n-1);
}

static double cor(int n, double x[], double y[]) {
  /* correlation of x and y. */
  int i=0;
  double sxy=0, sx2=0, sy2=0, mx=0, my=0, rx=0,ry=0;
  /* get the mean first, use the straight-forward formula */
  for(i=0; i<n; i++) {
    mx += x[i];
    my += y[i];
  }
  mx /= n;
  my /= n;

  for(i=0; i<n; i++) {
    rx = x[i] - mx;
    ry = y[i] - my;
    sx2 += rx*rx;
    sy2 += ry*ry;
    sxy += rx*ry;
  }
  return sxy/sqrt(sx2*sy2);
}

/***** Partial Correlation *****/
void test_func_pcor(int L, double x[], double y[], double* score, double* test_value) {
  double r=0, test_stat=0, p_value=0;
  r = cor(L, x,y);
  test_stat = r*sqrt((L-2)/(1-r*r));
  /* t distribution, df=L-2, two-sided test */
  p_value = 2.0*gsl_cdf_tdist_Q(fabs(test_stat), L-2);

  if(score) *score = -log10(p_value);
  if(test_value) *test_value = r;
}

void test_c_func_pcor(int L, double x[], double y[], int nZ, double Z[],
		      double* score,double* test_value) {
  /* get the partial correlation */
  double Sxx[2][2], r=0, test_stat=0, p_value=0;
  int sign=0, i=0,j=0;
  gsl_matrix *Szz=NULL;
  gsl_permutation *p=NULL;
  gsl_vector *Sxz1=NULL, *Sxz2=NULL, *v1=NULL, *v2=NULL;

  /* Sxx = cov([x,y], [x,y]), which is 2 by 2 */
  Sxx[0][0] = var(L, x);
  Sxx[1][0] = Sxx[0][1] = cov(L, x,y);
  Sxx[1][1] = var(L, y);
  /* Sxz = cov([x,y], Z), which is 2 by nZ */
  Sxz1 = gsl_vector_alloc(nZ);
  Sxz2 = gsl_vector_alloc(nZ);
  for(i=0; i<nZ; i++) {
    gsl_vector_set(Sxz1, i, cov(L, x, Z+i*L));
    gsl_vector_set(Sxz2, i, cov(L, y, Z+i*L));
  }
  /* Szz = cov(Z,Z), which is nZ by nZ */
  Szz = gsl_matrix_alloc(nZ,nZ);
  for(i=0; i<nZ; i++) {
    gsl_matrix_set(Szz,i,i, var(L, Z+i*L));
    for(j=i+1; j<nZ; j++) {
      r = cov(L, Z+i*L, Z+j*L);
      gsl_matrix_set(Szz,i,j, r);
      gsl_matrix_set(Szz,j,i, r);
    }
  }
  /* S.r = Sxx - Sxz * (Szz^-1) * (Sxz'), which is 2 by 2 */
  p = gsl_permutation_calloc(nZ);
  v1 = gsl_vector_alloc(nZ);
  v2 = gsl_vector_alloc(nZ);

  gsl_linalg_LU_decomp(Szz, p, &sign);
  /* [v1 v2] = (Szz^-1) * (Sxz') */
  gsl_linalg_LU_solve(Szz, p, Sxz1, v1);
  gsl_linalg_LU_solve(Szz, p, Sxz2, v2);

  gsl_blas_ddot(Sxz1,v1, &r);
  Sxx[0][0] -= r;
  gsl_blas_ddot(Sxz1,v2, &r); /* should be the same as Sxz2 . v1 */
  Sxx[1][0] -= r;
  Sxx[0][1] -= r;
  gsl_blas_ddot(Sxz2,v2, &r);
  Sxx[1][1] -= r;

  /* partial correlation */
  r = Sxx[1][0]/sqrt(fabs(Sxx[0][0]*Sxx[1][1]));

  /* p-value, t distribution, df = L-2-nZ */
  test_stat = r*sqrt((L-2-nZ)/(1-r*r));
  /*
  p_value = 2.0*gsl_cdf_tdist_Q(fabs(test_stat), L-2-nZ);
  */
  p_value = 2.0*gsl_cdf_ugaussian_Q(fabs(test_stat));

  /* debug */
  /*
  printf("=== test_c_func_pcor: L: %d r: %f test_stat: %f p_value: %f\n", L,r,test_stat,p_value);
  */

  if(score) *score = -log10(p_value);
  if(test_value) *test_value = r;
  /* clean up */
  gsl_permutation_free(p);
  gsl_matrix_free(Szz);
  gsl_vector_free(Sxz1);
  gsl_vector_free(Sxz2);
  gsl_vector_free(v1);
  gsl_vector_free(v2);
}
/***** Mutual Information *****/
/*
  Gaussian kernel probability density to estimate Mutual Information and Conditional Mutual Information
  Formula from Zhang et. al. 2011. "Inferring gene regulatory networks from gene expression data by PC-algorithm based on conditional mutual information". Bioinformatics, 2011.
*/

void test_func_mi(int L, double x[], double y[], double* score, double* test_value) {
  /* the unconditional one
     x and y are vectors of the same length
     Formula: I(X,Y)=(1/2)log(|C(X)|.|C(Y)|/|C(X,Y)|, where C(..) is the covariance matrix, |..| is determinant
     Here we assume x and y are only vectors, so covariance is easy
  */
  double vx=0, vy=0, cxy=0, r=0;
  vx = var(L, x);
  vy = var(L, y);
  cxy = cov(L, x,y);
  /* since here is the special case, and C(X,Y) is 2 by 2, we use the
     formula for its determinant directly */
  r = 0.5*log(vx*vy/(vx*vy-cxy*cxy));
  if(score) *score = r;
  if(test_value) *test_value = r;
}

static double det_cov(int L, double* x, double* y, int nZ, double Z[]) {
  /* det(cov(cbind(x,y,Z))).
     If either x or y is NULL, it is ignored.
   */
  int i=0, j=0, k=0, n=0;
  int sign=0;
  double det=0.0;
  gsl_permutation* p=NULL;
  gsl_matrix* A=NULL;

  /* fill in the matrix */
  n = nZ;
  if(x != NULL) n++;
  if(y != NULL) n++;
  A = gsl_matrix_calloc(n, n);
  /* fill the upper triangular part first */
  i=0;
  if(x != NULL) {
    gsl_matrix_set(A, i,i, var(L, x));
    j=i+1;
    if(y != NULL) {
      gsl_matrix_set(A, i,j, cov(L, x,y));
      j++;
    }
    for(k=0; k<nZ; k++)
      gsl_matrix_set(A, i,j+k, cov(L, x, Z+k*L));
    i++;
  }

  if(y != NULL) {
    gsl_matrix_set(A, i,i, var(L, y));
    for(k=0; k<nZ; k++)
      gsl_matrix_set(A, i,i+1+k, cov(L, y, Z+k*L));
    i++;
  }

  for(j=0; j<nZ; j++) {
    gsl_matrix_set(A, i+j,i+j, var(L, Z+j*L));
    for(k=j+1; k<nZ; k++)
      gsl_matrix_set(A, i+j,i+k, cov(L, Z+j*L, Z+k*L));
  }
  /* copy the upper triangular part to lower part */
  for(i=0; i<n; i++) {
    for(j=i+1; j<n; j++)
      gsl_matrix_set(A, j,i, gsl_matrix_get(A,i,j));
  }

  /* get determinant */
  p = gsl_permutation_calloc(n);
  gsl_linalg_LU_decomp(A, p, &sign);
  det = gsl_linalg_LU_det(A, sign);

  /* clean up */
  gsl_permutation_free(p);
  gsl_matrix_free(A);

  return det;
}
void test_c_func_mi(int L, double x[], double y[], int nZ, double Z[],
		    double* score,double* test_value) {
  /* the conditional one
     x and y should be vectors
     Z has nZ length L vectors in it.
     Formula: (1/2)log((|C(X,Z)|.|C(Y,Z)|)/(|C(Z)|.|C(X,Y,Z)|))
  */
  double cxz=0, cyz=0, cz=0, cxyz=0, r=0;
  cxz = det_cov(L, x,NULL, nZ,Z);
  cyz = det_cov(L, y,NULL, nZ,Z);
  cz = nZ==1 ? var(L, Z) : det_cov(L, NULL,NULL, nZ,Z);
  cxyz = det_cov(L, x,y, nZ,Z);

  r = 0.5*log((cxz*cyz)/(cz*cxyz));

  /* debug */
  /*
  printf("==mi L: %d nZ: %d cxz: %f cyz: %f cz: %f cxyz: %f r: %f\n",
	 L,nZ, cxz,cyz,cz,cxyz,r);
  printf("x:"); pr_doubles(L,x);
  printf("y:"); pr_doubles(L,y);
  {
    int i=0;
    printf("Z:\n");
    for(i=0; i<nZ; i++) {
      printf("%d:", i); pr_doubles(L,Z+i*L);
    }
    printf("\n");
  }
  */

  if(score) *score = r;
  if(test_value) *test_value = r;
}

/***************************************************************/
enum methods {METHOD_PCOR=0, METHOD_MI,  N_METHODS};
char* method_names[N_METHODS] = {"pcor", "mi"};
char* method_long_names[N_METHODS] = {"Partial Correlation", "Mutual Information"};

enum prunings {PRUNING_ALL=0, PRUNING_COMMON,  N_PRUNINGS};
char* pruning_names[N_PRUNINGS] = {"all", "common"};
/***********/

static int min_int_of(int n, int x[]) {
  int i=0;
  int m=x[0];
  for(i=1; i<n; i++)
    if(x[i] < m) m = x[i];
  return m;
}

static int max_int_of(int n, int x[]) {
  int i=0;
  int m=x[0];
  for(i=1; i<n; i++)
    if(x[i] > m) m = x[i];
  return m;
}

static int shift_by(double x[], int n_segs, int lens[],
		    int lag, int nL, double outx[]) {
  /* each segment of x is shifted by lag, and each segment's
     effective length is subtracted by nL.
     If a segment is too short, it is not included.

     Returns the length of the shifted series.
  */
  int i=0, j=0, k=0, n=0;
  for(j=0, k=0; k<n_segs; j+=lens[k++]) {
    n = lens[k] - nL;
    if(n <= 0) continue; /* segment too short */
    memcpy(outx+i, x+j+lag, sizeof(double)*n);
    i += n;
  }
  return i;
}

/* test whether x and y are significantly associated (e.g. by correlation),
   and return a score (larger means more significant) and the test_value
 */
typedef void (*test_func)(int L, double x[], double y[], double* score,double* test_value);

int test_link(int L, double x[], double y[], int max_delay, double st,
	      test_func f, int n_segs, int lens[],
	      link out_links[], double buf[]) {
  /* x and y have length L, representing two time-series.
     st is the score threshold, below which the link is not significant.
     f is a function that tests whether two vectors have significant association.
     lens contains n_segs lengths of the segments in x and y.
     out_links should have room for up to max_delay links.
     buf is for temporary use, and should have length 2L.

     test whether x --> y.
     return the number of links added.
   */
  int i=0, sL=0, n_links=0;
  double score=0, test_value=0;

  /* debug */
  /*
  printf("==== test_link\n");
  printf("lens"); pr_ints(n_segs,lens);
  printf("x"); pr_doubles(L,x);
  printf("y"); pr_doubles(L,y);
  */

  for(i=0; i<=max_delay; i++) {
    sL = shift_by(x,n_segs,lens, 0, i, buf);
    if(sL <= 5) continue; /* not long enough */
    shift_by(y,n_segs,lens, i, i, buf+L);

    /* debug */
    /*
    printf("at delay %d, x", i); pr_doubles(sL,buf);
    printf("at delay %d, y", i); pr_doubles(sL,buf+L);
    */

    f(sL, buf,buf+L, &score,&test_value);

    if(score >= st) { /* add a link */
      out_links[n_links].delay = i;
      out_links[n_links].score = score;
      out_links[n_links].test_value = test_value;
      n_links++;
    }
  }
  return n_links;
}

int infer_grn1(array2d* d, double st, int max_delay, test_func f,
	       int n_segs, int lens[],
	       int out_nL[], link* out_Ls[],
	       link buf_links[], double buf[]) {
  /* Stage 1 of CLINDE, infers the time lag and pairwise association.
     d is n by g array, the n_segs segments of expression data in column major format,
     where rows are time points, and columns are genes.
     The first lens[0] rows are segment 1,
     the next lens[1] rows are segment 2, and so on.
     st is the score threshold, below which the link is not significant.
     f is a function that tests whether two vectors have significant association.

     out_nL should have length g*g, which will be treated in row major
     format, so that out_nL[j+i*g] is the number of links of i --> j.
     out_Ls is similar to out_nL, but holds the pointer to the links,
     i.e. out_Ls[j+i*g] is the array of links of i --> j.

     buf_links should have room for up to g*g*max_delay links,
     and the pointers in out_Ls will point to here.

     buf should have room for up to 2n.

     Returns the number of links added.
   */
  int n=0, g=0, i=0, j=0, k=0, nL=0, n_links=0;
  n = d->rows;
  g = d->cols;

  for(i=0; i<g; i++) {
    for(j=0; j<g; j++) {
      /* test i --> j */
      /* currently no self loops considered */
      if(i == j) {
	out_nL[j+i*g] = 0;
	continue;
      }
      nL = test_link(n, get_col(d,i),get_col(d,j), max_delay,st,
		     f, n_segs,lens, buf_links+n_links,buf);
      out_nL[j+i*g] = nL;
      out_Ls[j+i*g] = buf_links+n_links;
      /* still need to fill in the from and to */
      for(k=0; k<nL; k++) {
	buf_links[n_links + k].from = i;
	buf_links[n_links + k].to = j;
      }
      n_links += nL;
    }
  }
  return n_links;
}

static int best_link(int n, link Ls[]) {
  double s = Ls[0].score;
  int i=1, idx=0;
  for(i=1; i<n; i++) {
    if(Ls[i].score > s) {
      s = Ls[i].score;
      idx = i;
    }
  }
  return idx;
}
int remove_dup_links(int g, int nL[], link* Ls[]) {
  /* both nL and Ls are g by g in row major format.
     nL[j+i*g] stores the number of links of i --> j.
     Ls[j+i*g] stores the array of links of i --> j, with different delays.
     modify Ls such that for each i --> j, there are at most one delay,
     by choosing the one with the highest score.

     Returns the number of links removed.
   */
  int i=0,j=0, k=0, n=0, n_removed=0;
  link* p=NULL;

  for(i=0; i<g; i++) {
    for(j=0; j<g; j++) {
      n = nL[j+i*g];
      if(n <= 1) continue;
      p = Ls[j+i*g];
      k = best_link(n, p);
      if(k > 0) p[0] = p[k];
      n_removed += n-1;
      nL[j+i*g] = 1;
    }
  }
  return n_removed;
}

int remove_dup_links2(int g, int nL[], link* Ls[], int n_links, link* p_links[]) {
  /* both nL and Ls are g by g in row major format.
     nL[j+i*g] stores the number of links of i --> j.
     Ls[j+i*g] stores the array of links of i --> j, with different delays.

     p_links are the n_links pointers to links.
     modify Ls such that for each i --> j, there are at most one delay,
     by choosing the one with the highest score.

     Returns the number of links after removal.
   */
  int i=0,j=0, k=0, n=0;
  int from=0, to=0;
  link* p=NULL;

  /* compact p_links, mark and remove dup links, adjust nL */
  for(i=0, j=0; i<n_links; i++) {

    /* debug */
    /*
    printf("compacting i: %d\tp_links[i]:%p\tfrom: %d\tto:%d\tdelay:%d\t%g\n", i,p_links[i],p_links[i]->from,p_links[i]->to,p_links[i]->delay, p_links[i]->score);
    */

    if(p_links[i]->delay < 0) { /* marked to be removed */
      p_links[i] = NULL;
      continue;
    }

    /* keep this link */
    from = p_links[i]->from;
    to = p_links[i]->to;
    n = nL[to + from*g];

    if(i!=j) p_links[j] = p_links[i];
    j++; /* to keep the link, increment no matter whether i==j */

    if(n > 1) { /* mark other delays to be removed */
      p = Ls[to + from*g];
      for(k=0; k<n; k++) {
	if(p+k != p_links[i])
	  p[k].delay = -1; /* mark it as not used */
      }
      nL[to + from*g] = 1;
    }
  }
  /* adjust Ls */

  /* debug */
  /*
  printf("number of links after no dup: %d\n", j);
  */

  for(i=0; i<j; i++) {
    from = p_links[i]->from;
    to = p_links[i]->to;

    /* debug */
    /*
    printf("i: %d\tp_links[i]:%p\tfrom: %d\tto:%d\tdelay:%d\n", i,p_links[i],from,to,p_links[i]->delay);
    printf("\tLs[to+from*g]: %p\tnL[to+from*g]:%d\n", Ls[to + from*g], nL[to + from*g]);
    fflush(stdout);
    */

    if(p_links[i] != Ls[to + from*g])
      Ls[to + from*g][0] = *(p_links[i]);
  }

  return j;
}

void print_links(FILE* f, int n_links, link* p_links[]) {
  int i=0;
  fprintf(f, "number of links: %d\n", n_links);
  for(i=0; i<n_links; i++) {
    fprintf(f, "To: %d From: %d Delay: %d Coef: %f Score: %f\n",
	    p_links[i]->to, p_links[i]->from,
	    p_links[i]->delay, p_links[i]->test_value,
	    p_links[i]->score);
  }
}

/***************************************************************/
int shift_c_vector(array2d* d, int n_segs, int lens[],
		   int ix, int iy, int lag,
		   int nZ, int iZ[], int dZ[],
		   double out_x[], double out_y[], double out_z[]) {
  /* To shift the x, y and Z by lags as appropriate.
     d is the data in column major format, n by g, with n_segs segments,
     each with lengths in lens.
     ix and iy are the 0-based column index into d.
     lag is the delay of x --> y.

     iZ contains nZ 0-based column indices into d,
     dZ contains the corresponding delays with respect to x.

     output the shifted parts into out_x and out_y for x and y.
     output the shifted parts into out_z for Z, column by column,
     so the total length of out_z will be length(out_x)*nZ.

     Assume the segments are long enough such that shifting does not
     make the segments crossing each other.

     Returns the length of the extracted out_x and out_y.
   */
  int md=0, k=0, L=0, sL=0;
  /* md = min(0, lag, dZ), each delay is to be subtracted from md */
  md = min_int_of(nZ, dZ);
  if(md > 0) md = 0;
  if(md > lag) md = lag;

  /* L = max(0 - md, lag - md, dZ - md).
     L is the length threshold, below which a segment will be shifted to empty
   */
  L = max_int_of(nZ,dZ) - md;
  if((0-md) > L) L = 0 - md;
  if((lag-md) > L) L = lag - md;

  /* take x, original delay is 0 */
  sL = shift_by(get_col(d, ix), n_segs,lens,  0 - md, L, out_x);
  if(sL <= 0) return 0;
  /* take y, original delay is lag */
  shift_by(get_col(d, iy), n_segs,lens,  lag - md, L, out_y);
  /* take Z, oringinal delays in dZ */
  for(k=0; k<nZ; k++) {
    shift_by(get_col(d, iZ[k]), n_segs,lens,  dZ[k] - md, L, out_z+k*sL);
  }
  return sL;
}

/* test whether x and y are significantly associated (e.g. by
   correlation), given Z, and return a score (larger means more
   significant) and the test_value.
   Z has nZ columns, each of length L.
 */
typedef void (*test_c_func)(int L, double x[], double y[], int nZ, double Z[], double* score,double* test_value);

int test_c_link(array2d* d, int n_segs, int lens[],
		double st, test_c_func f,
		int ix, int iy, int lag,
		int n_neis, int neis[],
		int n_d_neis[], int* d_neis[],
		int n_z, int z_set[],
		int ibuf[], double buf[]) {
  /* to do condition test of x --> y | Z.
     most parameters are same as in shift_c_vector().
     neis contains the n_neis indices of the neighbors of x --> y.
     n_d_neis contains the number of delays for each neighbor.
     d_neis contains the array of delays for each neighbor.
     z_set contains n_z indices into neis, z_set is the set of
     Z we are conditioning on.

     ibuf should have length 3g, where g=ncol(d),
     buf should have length (n_z+2)*nrow(d).

     st is score threshold, below which the test is insignificant,
     f is the testing function.

     Returns 0 if any of the test is insignificant.
     Returns 1 otherwise.
  */
  int i=0, n=0, g=0, L;
  int *iZ=NULL, *ic=NULL, *dZ=NULL;
  double *x=NULL, *y=NULL, *Z=NULL;
  double score=0, test_value=0;

  n = d->rows;
  g = d->cols;

  iZ = ibuf;
  ic = iZ+g;
  dZ = ic+g;

  x = buf;
  y = x + n;
  Z = y + n;

  for(i=0; i<n_z; i++) /* become indices of columns into d */
    iZ[i] = neis[z_set[i]];

  /* debug */
  /*
  printf("test_c_link, iZ:"); pr_ints(n_z,iZ); fflush(stdout);
  */

  /* try out the delays of Z */
  for(i=0; i<n_z; i++) ic[i] = 0; /* first set of delays, as counter */
  while(1) {
    /* get the delays */
    for(i=0; i<n_z; i++) dZ[i] = d_neis[z_set[i]][ic[i]];
    /* debug */
    /*
    printf("test_c_link, dZ:"); pr_ints(n_z,dZ); fflush(stdout);
    */

    L = shift_c_vector(d,n_segs,lens, ix,iy,lag, n_z,iZ,dZ, x,y,Z);
    if(L > 5) { /* ignore if too short */
      score = 0;
      f(L, x,y, n_z,Z, &score,&test_value);
      if(score < st) {
      /* test fail for whatever reason, or score not significant */
	return 0;
      }
    }
    /* increment the counter of the delays */
    for(i=n_z-1; i>=0; i--) {
      if(++ic[i] < n_d_neis[z_set[i]]) break; /* incremented, no overflow */
      /* one overflow, reset it, increment the next digit */
      ic[i] = 0;
    }
    if(i < 0) break; /* tried through the delays */
  }

  /* seems good */
  return 1;
}

/***************************************************************/
static int get_neighbors(int g, int ix, int iy, int nL[],
			 int* out_n_x, int* out_n_y,
			 int out_x[], int out_y[]) {
  /* neighbors of x and y except x and y themselves,
     nL[j+i*g] is the number of i --> j links.
     output the indices of the neighbors in out_x and out_y respectively,
     first the parents, then children, also ouput the counts.

     return the total number of neighbors.
   */
  int i=0, n_x=0,n_y=0;

  for(i=0; i<g; i++) {
    if(i==ix || i==iy) continue;
    if(nL[ix + i*g] + nL[i + ix*g] > 0) out_x[n_x++] = i;
    if(nL[iy + i*g] + nL[i + iy*g] > 0) out_y[n_y++] = i;
  }
  /* done */
  if(out_n_x) *out_n_x = n_x;
  if(out_n_y) *out_n_y = n_y;

  return n_x + n_y;
}

static int get_common_neighbors(int g, int ix, int iy, int nL[],
				int out[]) {
  /* common neighbors of x and y except x and y themselves,
     nL[j+i*g] is the number of i --> j links.
     output the indices of the common neighbors in out

     return the total number of common neighbors.
   */
  int i=0, n=0;

  /* into */
  for(i=0; i<g; i++) {
    if(i==ix || i==iy) continue;
    /* (i --> ix or ix --> i) && (i --> iy or iy --> i) */
    if((nL[ix + i*g] + nL[i + ix*g] > 0)
       && (nL[iy + i*g] + nL[i + iy*g] > 0))
      out[n++] = i;
  }

  return n;
}

static void get_the_delays_to_x(int g, int ix, int iy, int lag,
				int nL[], link* Ls[],
				int include_x,int include_y,
				int nZ, int iZ[],
				int out_ndZ[], int* out_dZ[], int ibuf[]) {
  /* get the delays of Z relative to x.
     g, ix, iy and nL have the same meaning as in get_common_neighbors().
     lag is the delay of x --> y.
     iZ contains nZ indices of the Z's.

     include_x indicates whether to include the delays of neighbors of x.
     include_y is similar, but for y.

     output the number of delays for each z in out_ndZ,
     output the delays for each z in out_dZ, and setup the pointer that
     point into ibuf.
     ibuf should be large enough to hold the needed delays.
   */
  int i=0, j=0, k=0, z=0, n=0;
  link* p=NULL;

  for(i=0; i<nZ; i++) {
    z = iZ[i];

    k = 0;
    if(include_x) {
      /* x --> z */
      n = nL[z + ix*g];
      p = Ls[z + ix*g];
      for(j=0; j<n; j++) ibuf[k++] = p[j].delay;
      /* z --> x */
      n = nL[ix + z*g];
      p = Ls[ix + z*g];
      for(j=0; j<n; j++) ibuf[k++] = -(p[j].delay);
    }

    if(include_y) {
      /* y --> z */
      n = nL[z + iy*g];
      p = Ls[z + iy*g];
      for(j=0; j<n; j++) ibuf[k++] = lag + (p[j].delay);
      /* z --> y */
      n = nL[iy + z*g];
      p = Ls[iy + z*g];
      for(j=0; j<n; j++) ibuf[k++] = lag - (p[j].delay);
    }
    /***/
    out_ndZ[i] = k;
    out_dZ[i] = ibuf;
    ibuf += k;
  }
}

static int enumerate_next_subset(int n, int size, int set[]) {
  /* set contains size numbers in 0 .. n-1, in strictly increasing order,
     representing a subset of size numbers in 0 .. n-1.
     Try to get to the next subset by increasing the numbers.
     Return 1 if OK. Return 0 if cannot enumerate the next subset anymore.
     set is modified in place
   */
  int i=0, j=0, k=0, idx;
  for(i=0; i<size; i++) { /* to start from the end */
    idx = size-1-i;
    if(++set[idx] < n-i) { /* the last one can hold a max of n-1, second last can hold a max of n-2 */
      /* OK, reinitialize the rest */
      k = set[idx];
      for(j=idx+1; j<size; j++) set[j] = ++k;
      return 1;
    }
  }
  /* tried through */
  return 0;
}

int filter_links(array2d* d, int n_segs, int lens[],
		 int nL[], link* Ls[],
		 int n_links, link* p_links[],
		 double st, int max_n, test_c_func f, int pruning) {
  /* Stage 2 of CLINDE.
     d is n by g array, the n_segs segments of expression data in column major format,
     where rows are time points, and columns are genes.
     The first lens[0] rows are segment 1,
     the next lens[1] rows are segment 2, and so on.

     st is the score threshold, below which a link is removed on conditional test.
     max_n is the maximum number of neighbors to condition on.
     f is the function for doing conditional test.

     nL[j+i*g] is the number of i --> j links,
     Ls[j+i*g] is the pointer to array of i --> j links.
     p_links contains the pointers to the n_links.

     pruning is either PRUNING_ALL (choose from all neighbors) or
     PRUNING_COMMON (choose only from common neighbors) for condition
     on neighbors.

     nL, Ls and p_links will be modified by removing links.
     Returns the number of links after removal.
   */
  int cn=0, i=0, j=0, is_end=0, n=0,g=0;
  int from=0, to=0, lag=0;
  int n_neis=0, n_x=0,n_y=0, include_x=0,include_y=0;
  int *neis=NULL, *nei_x=NULL,*nei_y=NULL,
    *idx=NULL, *iZ=NULL, *ibuf=NULL,
    *n_d_neis=NULL, *d_neis_buf=NULL;
  int** d_neis = NULL;
  double *buf=NULL;

  if(n_links <= 0) return 0; /* no links, no need to do anything */

  n = d->rows;
  g = d->cols;

  i = 8*g + n_links;

  neis = Malloc(i, int, "filter_links()");
  memset(neis, 0, sizeof(int)*i);
  nei_x = neis;            /* g */
  nei_y = nei_x+g;         /* g */
  idx = nei_y+g;           /* g */
  iZ = idx+g;              /* g */
  ibuf = iZ+g;             /* 3g */
  n_d_neis = ibuf+3*g;       /* g */
  d_neis_buf = n_d_neis+g; /* len n_links */

  d_neis = Malloc(g, int*, "filter_links()");
  memset(d_neis, 0, sizeof(int*)*g);

  buf = Malloc((max_n+2)*n, double, "filter_links()");
  memset(buf, 0, sizeof(double)*(max_n+2)*n);

  /* condition on more and more neighbors */
  for(cn=1; !is_end && cn<=max_n; cn++) {
    is_end = 1;

    /* debug */
    /*
    printf("====== iteration %d\n", cn);
    */

    /* n_links might be updated after the following loop */
    for(i=n_links-1; i>=0; i--) {
      if(p_links[i] == NULL) continue;

      /* debug */
      /*
      printf("p_link[%d] (%p) before adjust: to: %d from: %d lag: %d\n", i, p_links[i], p_links[i]->to, p_links[i]->from, p_links[i]->delay);
      */

      /* test links with lower scores first */
      from = p_links[i]->from;
      while(from < 0) { /* marked to update p_links */
	/* may need to update a few times */
	p_links[i] += from;
	from = p_links[i]->from;
      }
      to = p_links[i]->to;
      lag = p_links[i]->delay;

      /* debug */
      /*
      printf("p_link[%d] (%p) after adjust: to: %d from: %d lag: %d\n", i, p_links[i], p_links[i]->to, p_links[i]->from, p_links[i]->delay);
      printf("== p_links[%d]: to: %d from: %d lag: %d\n", i, to,from,lag);
      */

      /* determine neighbors to condition on */
      include_x=1;
      include_y=1;

      if(pruning == PRUNING_COMMON) {
	n_neis = get_common_neighbors(g, from,to, nL, neis);
      } else {
	get_neighbors(g, from,to, nL, &n_x,&n_y, nei_x,nei_y);
	/* debug */
	/*
	printf("= prune all\n");
	printf("nei_x"); pr_ints(n_x,nei_x);
	printf("nei_y"); pr_ints(n_y,nei_y);
	*/

	if(n_x < n_y) {
	  n_neis = n_x;
	  include_y = 0;
	} else {
	  n_neis = n_y;
	  memcpy(neis, nei_y, n_neis*sizeof(int));
	  include_x = 0;
	}
      }

      /* debug */
      /*
      printf("== neis");
      pr_ints(n_neis, neis);
      */

      if(n_neis < cn) continue;
      /* need to check */
      is_end = 0;

      get_the_delays_to_x(g, from,to,lag, nL,Ls, include_x,include_y,
			  n_neis,neis,
			  n_d_neis,d_neis, d_neis_buf);

      /* debug */
      /*
      printf("n_d_neis"); pr_ints(n_neis, n_d_neis);
      printf("d_neis"); pr_int_arr(n_neis, n_d_neis,d_neis);
      */

      /* generate combinations of neighbors */
      for(j=0; j<cn; j++) idx[j] = j; /* first subset, 0-based index */
      do {
	/* debug */
	/*
	printf("subset"); pr_ints(cn,idx);
	*/

	/* do conditional test */
	if(!test_c_link(d,n_segs,lens, st,f, from,to,lag,
			n_neis,neis, n_d_neis,d_neis,
			cn,idx, ibuf,buf)) {
	  /* remove it */
	  {
	    /* debug */
	    /*
	    printf("=== remove p_links[%d]\n", i);
	    */

	    int m = nL[to + from*g];
	    link* pL = Ls[to + from*g] + m-1;

	    /* debug */
	    /*
	    printf("pL (%p): to: %d from: %d lag: %d\n", pL, pL->to, pL->from, pL->delay);
	    printf("pL - p_links[i]: %d\n", (pL - p_links[i]));
	    */

	    if(pL != p_links[i]) {
	      *(p_links[i]) = *pL;
	      /* to mark the last link to update p_links when encountered */
	      pL->from = -(pL - p_links[i]);
	    }

	    /* debug */
	    /*
	    printf("p_link[%d] (%p) after: to: %d from: %d lag: %d\n", i, p_links[i], p_links[i]->to, p_links[i]->from, p_links[i]->delay);
	    printf("pL (%p) after: to: %d from: %d lag: %d\n", pL, pL->to, pL->from, pL->delay);
	    */

	    nL[to + from*g]--;
	    p_links[i] = NULL;
	  }
	  break;
	}
      } while(enumerate_next_subset(n_neis, cn, idx));
    }
    /* compact p_links if needed */
    for(i=0, j=0; i<n_links; i++) {
      if(p_links[i] != NULL) {
	/* correct the possibly -1 from */
	while(p_links[i]->from < 0) p_links[i] += p_links[i]->from;

	if(i != j) p_links[j] = p_links[i];
	j++; /* to keep the link, increment no matter whether i==j */
      }
    }
    n_links = j;

    /* debug */
    /*
    printf("=== links after iteration %d\n", cn);
    print_links(stdout, n_links, p_links);
    */
  }
  /* clean up */
  free(neis);
  free(d_neis);
  free(buf);
  return n_links;
}

/***************************************************************/
int decreasing_link_score(const void* a, const void* b) {
  link *A,*B;
  A = *(link**)a;
  B = *(link**)b;

  if(A->score > B->score) return -1;
  if(A->score < B->score) return 1;

  return A->delay - B->delay;
}
static int sort_links(int g, int nL[], link* Ls[], link* out_links[]) {
  /* nL[j+i*g] is the number of i --> j links,
     Ls[j+i*g] is the pointer to array of i --> j links.
     sort the links by decreasing score and put the pointers in out_links.

     Returns the number of links.
   */
  int i=0, j=0, k=0, n=0, n_links=0;
  link* p=NULL;

  for(i=0; i<g; i++) {
    for(j=0; j<g; j++) {
      n = nL[j+i*g];
      if(n <= 0) continue;
      p = Ls[j+i*g];
      for(k=0; k<n; k++) {
	out_links[n_links++] = p+k;
      }
    }
  }
  /* sort it */
  qsort(out_links, n_links, sizeof(link*), decreasing_link_score);
  /**/
  return n_links;
}

void clinde(array2d* d, int n_segs, int lens[],
	    double st1, double st2, int max_delay, int max_n,
	    int method, int pruning, 
	    int is_one_delay, int is_no_dup) {
  /* main function of CLINDE.
     d is n by g array, the n_segs segments of expression data in column major format,
     where rows are time points, and columns are genes.
     The first lens[0] rows are segment 1,
     the next lens[1] rows are segment 2, and so on.

     st1 and st2 are the score thresholds for stage 1 and 2 respectively.
     max_n is the maximum number of neighbors to condition on.

     method is either METHOD_PCOR (partial correlation) or METHOD_MI (mutual information).
     pruning is either PRUNING_ALL (choose from all neighbors) or PRUNING_COMMON (choose only from common neighbors) for condition on neighbors in stage 2.

     if is_one_delay is true, remove duplicate links after stage 1.
     if is_no_dup is true, remove duplicate links after stage 2.

   */
  int n_links=0, g=0,n=0;
  int* nL=NULL;
  link** Ls=NULL;
  link** p_links=NULL;
  link* buf_links=NULL;
  double* buf=NULL;

  n = d->rows;
  g = d->cols;

  /* prepare buffers */
  nL = Malloc(g*g, int, "infer_grn()");
  memset(nL, 0, sizeof(int)*g*g);

  Ls = Malloc(g*g, link*, "infer_grn()");
  memset(Ls, 0, sizeof(link*)*g*g);

  p_links = Malloc(g*g*(1+max_delay), link*, "infer_grn()");
  memset(p_links, 0, sizeof(link*)*g*g*(1+max_delay));

  buf_links = Malloc(g*g*(1+max_delay), link, "infer_grn()");
  memset(buf_links, 0, sizeof(link)*g*g*(1+max_delay));

  buf = Malloc(2*n, double, "infer_grn()");
  memset(buf, 0, sizeof(double)*2*n);

  /* stage 1 */
  report_stop_watch();
  printf("==== Stage 1\n");
  infer_grn1(d, st1, max_delay,
	     method==METHOD_MI ? test_func_mi : test_func_pcor,
	     n_segs, lens, nL,Ls, buf_links, buf);
  if(is_one_delay) remove_dup_links(g, nL, Ls);
  n_links = sort_links(g, nL, Ls, p_links);

  printf("==== Initial GRN after stage 1:\n");
  print_links(stdout, n_links, p_links);
  printf("==== End stage 1 ==============\n");

  report_stop_watch();
  fflush(stdout);
  /* stage 2 */
  printf("==== Stage 2\n");
  n_links = filter_links(d, n_segs, lens, nL, Ls, n_links, p_links,
			 st2, max_n,
			 method==METHOD_MI ? test_c_func_mi : test_c_func_pcor,
			 pruning);
  if(is_no_dup) n_links = remove_dup_links2(g, nL, Ls, n_links, p_links);

  printf("==== Final GRN after stage 2:\n");
  print_links(stdout, n_links, p_links);
  printf("==== End stage 2 ============\n");
  fflush(stdout);
  /* clean up */
  free(nL);
  free(Ls);
  free(p_links);
  free(buf_links);
  free(buf);
}

/***************************************************************/
array2d* read_expression_data(char* filename) {
  /* return an 2d array (row major) if can read properly. return NULL on error. */
  array2d* r=NULL;
  int rows=0, cols=0;
  double* res=read_TSV(filename, &rows, &cols);
  if(res == NULL) {
    fprintf(stderr, "Error: Cannot properly read %s for expression data.\n", filename);
    return NULL;
  }
  r = Malloc(1, array2d, "read_expression_data()");
  r->rows = rows;
  r->cols = cols;
  r->v = res;
  return r;
}

array2d* read_n_data(int n, char* fn[], int out_lens[]) {
  /* read n expression data, fn are the filenames.
     The expression data should have the same number of columns.
     Put them one by one in the same 2D array, and output the number
     of rows of each data in out_lens.
     Return the array in column major format, so that it is easier to extract a column.
     Return the 2D array if OK. Return NULL if error.
   */
  int i=0,j=0,k=0, ncols=0, nrows=0, nr=0;
  array2d* r = NULL;
  array2d** buf = NULL;
  buf = Malloc(n,array2d*, "read_n_data()");
  memset(buf, 0, sizeof(array2d*)*n);

  printf("About to read %d expression data\n", n);

  for(i=0; i<n; i++) {
    printf("Reading the %dth expression data %s ... ", i+1, fn[i]);
    buf[i] = read_expression_data(fn[i]);
    if(buf[i] == NULL) goto done;
    if(i==0) {
      ncols = buf[i]->cols;
    } else if(ncols != buf[i]->cols) {
      fprintf(stderr, "Error: the %dth expression data %s has %d columns, but previous ones have %d columns.\n",
	      (i+1),fn[i], buf[i]->cols, ncols);
      goto done;
    }
    printf("OK.\n  Read %d rows, %d cols.\n", buf[i]->rows, buf[i]->cols);
    nrows += buf[i]->rows;
  }
  /* copy them to one array, in column major format */
  r = new_array2d(nrows, ncols);
  for(i=0, nr=0; i<n; i++) {
    out_lens[i] = buf[i]->rows;;
    for(j=0; j<buf[i]->rows; j++, nr++)
      for(k=0; k<ncols; k++)
	caset(r, nr,k, aref(buf[i], j,k));
  }
  printf("Done\nRead %d time points, %d genes.\n", nrows, ncols);
  /***/
 done:
  for(i=0; i<n; i++)
    if(buf[i]) free_array2d(buf[i]);
  free(buf);
  return r;
}

/***************************************************************/
#define DEFAULT_ST1        2
#define DEFAULT_ST2        2
#define DEFAULT_MAX_DELAY  4
#define DEFAULT_MAX_N      4
#define DEFAULT_METHOD     "pcor"
#define DEFAULT_PRUNING    "all"

#define xstr(s) str(s)
#define str(s) #s

/* The option structure of this program */
option all_options[] = {
  /* str,type(STRING),is_required(0),description,index(0),end_index(0), val_int(0),val_real(0),val_str(NULL) */
  {"-?",         FLAG,   0,  "Showing the usage.\n", 0,0,  0,  0.0,    NULL},

  {"-data",      LIST,   1,  "File name(s) of the input expression data (tab/space separated), each row is a time point, each column is a gene. Each file should have the same number of columns, but may have different number of rows.\n", 0,-1,  -1,  0.0,    NULL},

  {"-st1",       REAL,   0,  "Score threshold for stage 1. For method=pcor, score is -log10(p), where p is the intended p-value. For method=mi, score is the mutual information. Default " xstr(DEFAULT_ST1) ".\n", 0,0,  0,  DEFAULT_ST1,    NULL},
  {"-st2",       REAL,   0,  "Score threshold for stage 2. Similar to st1. Default " xstr(DEFAULT_ST2) ".\n", 0,0,  0,  DEFAULT_ST2,    NULL},
  {"-max.delay", INT,    0,  "Maximum delay (lags) to use in the inference, default " xstr(DEFAULT_MAX_DELAY) ".\n", 0,0,  DEFAULT_MAX_DELAY,  0.0,    NULL},
  {"-max.n",     INT,    0,  "Maximum number of parents to condition on in stage 2, default " xstr(DEFAULT_MAX_N) ".\n", 0,0,  DEFAULT_MAX_N,  0.0,    NULL},

  {"-method",    STRING, 0,  "Method of testing links in stage 1 and 2, can be pcor (partial correlation) or mi (mutual information). Default " xstr(DEFAULT_METHOD) "\n", 0,0,  0,  0.0,    DEFAULT_METHOD},
  {"-pruning",   STRING, 0,  "Pruning strategy in stage 2, can be all (consider all neighbors of the two vertices when doing conditional test of the link between the two) or common (consider only common neighbors of the two vertices). Default " xstr(DEFAULT_PRUNING) "\n", 0,0,  0,  0.0,    DEFAULT_PRUNING},
  /*
  {"-local.delay", FLAG,   0,  "In stage 1, in addition to score threshold, also use local maximum of absolute correlation to filter possible delays. Default false.\n", 0,0,  0,  0.0,    NULL},
  */
  {"-one.delay", FLAG,   0,  "To keep only one delay (the one with best score, smallest delay) for each link after stage 1. Default false.\n", 0,0,  0,  0.0,    NULL},
  {"-no.dup",    FLAG,   0,  "Remove duplicate links. To keep only one delay (the one with best score, smallest delay) for each link after stage 2. Default false.\n", 0,0,  0,  0.0,    NULL},

  {NULL,      STRING, 0,  NULL, 0,0,  0,  0.0,    NULL} /* end marker */
};

/***********/
int our_stricmp(char* a, char* b) {
  int r, i=0;
  do {
    r = tolower(a[i]) - tolower(b[i]);
    if(r > 0) return 1;
    if(r < 0) return -1;
    i++;
  } while(a[i]!='\0' || b[i]!='\0');
  return 0;
}

int match_option(char* str, int n, char* choices[],
		 int warning_if_not_found, FILE* outs, char* option_name) {
  /* return -1 if not found. return the index if found */
  /* option_name is for the warning */
  int i=0;
  for(i=0; i<n; i++) {
    if(our_stricmp(str, choices[i])==0) return i;
  }
  if(!warning_if_not_found) return -1;

  /* not found, indicate the choices */
  fprintf(outs, "Warning: Unknown option \"%s\", possible choices are:\n", str);
  for(i=0; i<n; i++)
    fprintf(outs, "\t%s\n", choices[i]);
  fprintf(outs, "\n");
  return -1;
}

/***************************************************************/

int main(int argc, char* argv[]) {
  double st1=DEFAULT_ST1, st2=DEFAULT_ST2;
  int max_delay=DEFAULT_MAX_DELAY, max_n=DEFAULT_MAX_N;
  char* s_method=DEFAULT_METHOD;
  char* s_pruning=DEFAULT_PRUNING;
  int method=0, pruning=0;
  int is_one_delay=0, is_no_dup=0;

  int s_idx=0, e_idx=-1;
  int n_data=0, *data_lens=NULL;
  array2d* data=NULL;

  /* do some initialization */
  if(parse_options(argc,argv,all_options)) {
    /* error parsing options */
    usage(stderr, argv[0], all_options);
    exit(-1);
  }
  if(get_flag_option_value("-?",all_options,0))
    usage(stdout, argv[0], all_options);

  if(!get_list_option_range("-data", all_options, &s_idx, &e_idx) || (e_idx < s_idx)) {
    fprintf(stderr, "Input expression file names not given.\n");
    usage(stderr, argv[0], all_options);
    exit(-1);
  }
  n_data = e_idx - s_idx; /* both inclusive, but start includes the -data */

  st1 = get_real_option_value("-st1", all_options, st1);
  st2 = get_real_option_value("-st2", all_options, st2);
  max_delay = get_int_option_value("-max.delay", all_options, max_delay);
  max_n = get_int_option_value("-max.n", all_options, max_n);
  s_method = get_str_option_value("-method", all_options, s_method);
  s_pruning = get_str_option_value("-pruning", all_options, s_pruning);
  is_one_delay = get_flag_option_value("-one.delay", all_options, is_one_delay);
  is_no_dup = get_flag_option_value("-no.dup", all_options, is_no_dup);

  if((method = match_option(s_method, N_METHODS, method_names, 1,stdout,"-method")) < 0) {
    printf("Use default for method instead.\n");
    method = 0;
  }

  if((pruning = match_option(s_pruning, N_PRUNINGS, pruning_names, 1,stdout,"-pruning")) < 0) {
    printf("Use default for pruning instead.\n");
    pruning = 0;
  }

  if(max_delay <= 0) {
    printf("max.delay should be positive integer, but got %d\n", max_delay);
    printf("Use default instead.\n");
    max_delay = DEFAULT_MAX_DELAY;
  }

  /*********/
  printf("====== CLINDE Grn Inference ======\n");
  start_stop_watch();

  printf("Stage 1 score threshold: %f\n", st1);
  printf("Stage 2 score threshold: %f\n", st2);
  printf("Upper bound (exclusive) of delays to test: %d\n", max_delay);
  printf("Maximum Parents to condition on in stage 2: %d\n", max_n);
  printf("Method: %s (%s)\n", method_names[method], method_long_names[method]);
  printf("Pruning: %s\n", pruning_names[pruning]);
  if(is_one_delay)
    printf("Retain only one delay for each link after stage 1.\n");
  if(is_no_dup)
    printf("Retain only one delay for ech link, remove duplicate links after stage 2.\n");
  /*********/

  data_lens = Malloc(n_data, int, "data_lens");
  data = read_n_data(n_data, argv+s_idx+1, data_lens);
  if(data == NULL) return 1;
  report_stop_watch();

  /* debug */
  /*
  printf("data_lens"); pr_ints(n_data, data_lens);
  */

  /*********/
  clinde(data, n_data,data_lens, st1,st2,max_delay,max_n,
	 method,pruning, is_one_delay, is_no_dup);

  report_stop_watch();
  /*********/
  free_array2d(data);
  free(data_lens);

  return 0;
}

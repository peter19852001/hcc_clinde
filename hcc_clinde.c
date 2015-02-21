/*
  Use modified CLINDE and clustering to find hidden common cause frm
  time series data.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <assert.h>
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
#include <gsl/gsl_multifit.h>

/************************************************************************/
/* mainly for debug */

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
static double mean(int n, double x[]) {
  /* arithmetic mean */
  int i=0;
  double s=0;
  for(i=0; i<n; i++) s+=x[i];
  return s/n;
}

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

static void normalize(int n, double x[], double* out_mean, double* out_sd) {
  /* normalize to zero mean and s.d. 1 (unless s.d. is too small) */
  int i=0;
  double sx2=0, mx=0;
  for(i=0; i<n; i++) mx += x[i];
  mx /= n;
  for(i=0; i<n; i++) {
    x[i] -= mx;
    sx2 += x[i]*x[i];
  }
  sx2 = sqrt(sx2/(n-1));
  if(sx2 > 1.0e-10) {
    for(i=0; i<n; i++)
      x[i] /= sx2;
  }

  if(out_mean) *out_mean = mx;
  if(out_sd) *out_sd = sx2;
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

static double seg_cor(int n, double x[], double y[],
		      double* Sxy, double* Sx2, double* Sy2) {
  /* correlation of x and y.
     But also accumulate the following:
       Sxy: sum of (x[i] - mx)*(y[i] - my)
       Sx2: sum of (x[i] - mx)^2
       Sy2: sum of (y[i] - my)^2
  */
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

  if(Sxy) *Sxy += sxy;
  if(Sx2) *Sx2 += sx2;
  if(Sy2) *Sy2 += sy2;

  return sxy/sqrt(sx2*sy2);
}

static void simple_linear_regression(int n, double x[], double y[], double* intercept, double* slope) {
  /* fit y = intercept + slope*x,
     use a straightforward approach.
   */
  *slope = cov(n,x,y)/var(n,x);
  *intercept = mean(n,y) - (*slope)*mean(n,x);
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

  /* debug
  printf("test_c_func_pcor(), r: %f\n", r);
  */

  /* p-value, t distribution, df = L-2-nZ */
  test_stat = r*sqrt((L-2-nZ)/(1.0000001-r*r));

  /* debug
  printf("test_c_func_pcor(), test_stat: %f\n", test_stat);
  */

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
	       int n_segs, int lens[], char* forbidden,
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

     if forbidden is non-NULL, then if forbidden[i*g + j] is 1, then the link i->j is forbidden.

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
      if(i == j || (forbidden!=NULL && forbidden[i*g+j]==1)) {
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

int remove_zero_delay_links(int n_links, link* p_links[]) {
  /* p_links are the n_links pointers to links.
     modify p_links to remove those with zero delays.
     returns the number of links after removal.
   */
  int i=0, j=0;
  for(i=0; i<n_links; i++) {
    if(p_links[i]->delay == 0) continue;
    p_links[j] = p_links[i];
    j++;
  }
  return j;
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

    if(p_links[i] != Ls[to + from*g]) {
      Ls[to + from*g][0] = *(p_links[i]);
      p_links[i] = Ls[to + from*g];
    }
  }

  return j;
}

void print_link(FILE* f, link* L) {
  fprintf(f, "To: %d From: %d Delay: %d Coef: %f Score: %f\n",
	  L->to, L->from, L->delay, L->test_value, L->score);
}
void print_links(FILE* f, int n_links, link* p_links[]) {
  int i=0;
  fprintf(f, "number of links: %d\n", n_links);
  for(i=0; i<n_links; i++) {
    print_link(f, p_links[i]);
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

  /* debug 
  printf("test_c_link, ix:%d, iy:%d, iZ:", ix,iy); pr_ints(n_z,iZ); fflush(stdout);
  */


  /* try out the delays of Z */
  for(i=0; i<n_z; i++) ic[i] = 0; /* first set of delays, as counter */
  while(1) {
    /* get the delays */
    for(i=0; i<n_z; i++) dZ[i] = d_neis[z_set[i]][ic[i]];
    /* debug 
    printf("test_c_link, dZ:"); pr_ints(n_z,dZ); fflush(stdout);
    */

    L = shift_c_vector(d,n_segs,lens, ix,iy,lag, n_z,iZ,dZ, x,y,Z);
    if(L > 5) { /* ignore if too short */
      score = 0;
      f(L, x,y, n_z,Z, &score,&test_value);

      /* debug 
      printf("test_c_link, score: %f\n", score); fflush(stdout);
      */

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

int filter_hcc_links(int ng, int g, int nL[], link* Ls[], int n_links, link* p_links[]) {
  /* ng is the number of observed genes.
     g is the number of observed plus introduced hidden common nodes.
     The index of introduced hidden common nodes are >= ng but < g.

     nL[j+i*g] is the number of i --> j links,
     Ls[j+i*g] is the pointer to array of i --> j links.
     p_links contains the pointers to the n_links.

     nL, Ls and p_links will be modified by removing links.
     Returns the number of links after removal.
  */
  int h=0, x=0, y=0, z=0, i=0;
  link* pL=NULL;
  /* by brute force, for all h->x and h->y, where h is hidden, x and y
     are observed, if there is an observed z such that either
          x->z-> y
       or y->z-> x
       or (z->x and z->y)
     then remove h->x and h->y
   */
  for(h=ng; h<g; h++) {
    for(x=0; x<ng; x++) {
      if(nL[x + h*g] <= 0) continue;
      /* h->x */
      for(y=x+1; y<ng; y++) {
	if(nL[y + h*g] <= 0) continue;
	/* h->y */
	for(z=0; z<ng; z++) {
	  if((nL[z + x*g] > 0 && nL[y + z*g] > 0) /* x->z->y */
	     || (nL[x + z*g] > 0 && (nL[y + z*g] > 0 || nL[z + y*g] > 0))) { /* z->x and (z->y or y->z) */
	    /* OK, mark h->x and h->y for removal by marking the delay as -1 */
	    for(i=nL[x + h*g], pL=Ls[x + h*g]; i>0; i--, pL++) {pL->delay = -1;}
	    for(i=nL[y + h*g], pL=Ls[y + h*g]; i>0; i--, pL++) {pL->delay = -1;}
	    break;
	  }
	}
      }
    }
    /* check if h still have any children */
    for(i=0; i<g; i++) {
      if(nL[i+h*g] > 0 && Ls[i+h*g]->delay>=0) break;
    }
    if(i >= g) { /* h has no children, also mark its parent links for removal */
      for(x=0; x<g; x++) {
	for(i=nL[h + x*g], pL=Ls[h + x*g]; i>0; i--, pL++) {pL->delay = -1;}
      }
    }
  }
  /* now remove the links marked and compact p_links */
  h = 0;
  for(i=0; i<n_links; i++) {
    pL = p_links[i];
    if(pL == NULL) continue;
    if(pL->delay < 0) {
      nL[(pL->to)+(pL->from)*g] = 0;
      p_links[i] = NULL;
    } else {
      p_links[h] = pL;
      h++;
    }
  }
  return h;
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

int clinde(array2d* d, int n_segs, int lens[],
	   double st1, double st2, int max_delay, int max_n,
	   int method, int pruning, 
	   int is_one_delay, int is_no_dup, int is_keep_zero_delay,
	   char* forbidden, int ng, link** out_links) {
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
     if is_keep_zero_delay is true, keep links with zero delay after stage 2 and before removing dups.
     if is_no_dup is true, remove duplicate links after stage 2.

     if forbidden is non-NULL, then if forbidden[i*g + j] is 1, then the link i->j is forbidden.
     ng is the number of observed genes, and the first ng columns of d correspond to the observed genes.
     if out_links is non-NULL and there are links, *out_links will be written with a pointer to the final links, and the buffer should be freed after use.

     returns the number of final links.
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
	     n_segs, lens, forbidden, nL,Ls, buf_links, buf);
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
  /* add a part to remove links from introduced hidden common node */
  if(g > ng) {
    n_links = filter_hcc_links(ng,g, nL,Ls, n_links,p_links);
  }
  /***/
  if(!is_keep_zero_delay) n_links = remove_zero_delay_links(n_links, p_links);
  /* note that nL and Ls may have more links than in p_links,
     so the main reference are n_links and p_links.
   */
  if(is_no_dup) n_links = remove_dup_links2(g, nL, Ls, n_links, p_links);

  printf("==== \n");
  print_links(stdout, n_links, p_links);
  printf("==== End stage 2 ============\n");
  fflush(stdout);

  /* copy the links to out_links if applicable */
  if(out_links != NULL && n_links>0) {
    {
      link* outL = NULL;
      int i=0;
      outL = Malloc(n_links, link, "clinde()");
      for(i=0; i<n_links; i++) {
	outL[i] = *(p_links[i]);
      }
      *out_links = outL;
    }
  }

  /* clean up */
  free(nL);
  free(Ls);
  free(p_links);
  free(buf_links);
  free(buf);

  return n_links;
}

/***************************************************************/
/* A simple heuristic sequential clustering */

/* for recording the series included in each cluster */
typedef struct s_series {
  int id;
  int idx; /* column index */
  int offset;
  double variance; /* variance of the column */
  struct s_series* next;
} series;

typedef struct s_cluster {
  int id;
  int n_series; /* number of series merged to here */
  int capacity; /* length of buffer of mean and n, divided into segments */
  int s,e; /* 0-based start (inclusive) and end index (exclusive) of the section used */
  double* mean; /* the average of each time point, own this */
  double* n; /* the weight of values in the average, for each time point, own this */
  series* merged; /* not own this */
  struct s_cluster* next;
} cluster;

cluster* new_cluster(int id, int n_series, int capacity, int L,
		     double* init_mean, int init_offset,
		     double* init_n, series* se,
		     cluster* next) {
  /* allocate a new cluster with the provided initial information (deep copy).
     If init_mean is NULL, initialize to 0, and init_n is ignored and n set to 0.
     If init_mean is non-NULL, copy its content as the mean; and if init_n is
     provided, it is copied as n, otherwise 1 is used.
  */
  cluster* c = NULL;
  int i=0;

  c = Malloc(1,cluster, "new_cluster()");
  c->id = id;
  c->n_series = n_series;
  c->capacity = capacity;
  c->mean = Malloc(capacity, double, "new_cluster()");
  memset(c->mean, 0, sizeof(double)*capacity);
  c->n = Malloc(capacity, double, "new_cluster()");
  memset(c->n, 0, sizeof(double)*capacity);

  assert(L <= capacity);
  c->s = init_offset;
  c->e = c->s + L;

  if(init_mean != NULL) {
    memcpy(c->mean + c->s, init_mean, sizeof(double)*L);
    if(init_n != NULL) {
      memcpy(c->n + c->s, init_n, sizeof(double)*L);
    } else {
      for(i=(c->s); i<(c->e); i++) c->n[i] = 1;
    }
  }

  c->merged = se; /* assume se is properly null-terminated */
  c->next = next;
  return c;
}

cluster* single_cluster(int id, int capacity, int L, double* init_mean, int init_offset, series* se, cluster* next) {
  cluster* r=NULL;
  r = new_cluster(id, 1, capacity, L, init_mean, init_offset, NULL, se, next);
  /* se->id is assumed to have been set properly */
  se->offset = r->s;
  se->next = NULL;
  return r;
}

void free_clusters(cluster* cs) {
  cluster* p=NULL;
  while(cs != NULL) {
    p = cs;
    cs = cs->next;
    if(p->mean) free(p->mean);
    if(p->n) free(p->n);
    free(p);
  }
}

void merge_to_cluster(cluster* c, series* se, int m, double es[], int offset) {
  /* add es (of length m) to c, where the offset and corr is relative to c.
   */
  int s=0, e=0, i=0;
  double f=1, intercept=0, slope=0;

  if(offset+m > (c->capacity)) {m = c->capacity - offset;}

  /* determine the overlapping part */
  s = offset > (c->s) ? offset : (c->s);
  e = (offset + m) < (c->e) ? (offset + m) : (c->e);

  assert(e > s);

  simple_linear_regression(e-s, c->mean+s, es, &intercept,&slope);
  if(slope < 0) f=-1;

  /* merge, assuming c->mean is the current estimate of the hidden
     common cause.
     fit es = intercept + slope*(hidden), then add es to the estimation
     and re-normalize.
  */
  for(i=0; i<m; i++) {
    c->mean[offset+i] = ((c->mean[offset+i])*(c->n[offset+i])
			 + f*(es[i] - intercept))/(c->n[offset+i] + f*slope);
    c->n[offset+i] += f*slope;
  }
  /* adjust s and e */
  if(offset < (c->s)) c->s = offset;
  if(offset+m > (c->e)) c->e = offset+m; /* offset+m is <= capacity, due to the adjustment above */

  c->n_series++;
  se->offset = offset;
  se->next = c->merged;
  c->merged = se;
}

double correlation_to_cluster(cluster* c, int n_segs, int lens[],
			      double es[], int* offset) {
  /* calculate the max. abs. correlation of es (of length m) to cluster c,
     return the correlation, and output the corresponding delay.
     Assume es is wholly witin the capacity of c.
     The offset is relative to c.
   */
  int s=0, best_offset=0, dx=0,dy=0, nL=0;
  int i=0, md=0, sL=0;
  double best_corr=0, tmp_c=0;
  double sxy=0, sx2=0, sy2=0;

  /* assume that the initial space c->s is the allowed maximum delay. */
  md = c->s;

  assert(c != NULL);

  for(s=0; s<=2*md; s++) {
    /* consider the segments, which are packed together starting at c->s */
    if(s < (c->s)) {
      dx = c->s;
      dy = c->s - s;
    } else {
      dx = s;
      dy = 0;
    }

    sxy=0; sx2=1e-10; sy2=1e-10;
    for(i=0, sL=0; i<n_segs; i++, sL += lens[i]) {
      /* overlapping length within one segment */
      nL = s < (c->s) ? (lens[i] - dy) : (lens[i] + (c->s) - s);

      /* accumulate the correlation parts, assume the means are
	 calculated well for each segment.  This may give slightly
	 different results to concatenating the overlapping segments
	 and calculate the correlation, but allows each segment to
	 have different means.
      */
      seg_cor(nL, c->mean+sL+dx, es+sL+dy, &sxy,&sx2,&sy2);
    }
    tmp_c = sxy/sqrt(sx2*sy2);

    if(fabs(tmp_c) > fabs(best_corr)) {
      best_corr = tmp_c;
      best_offset = s;
    }
  }

  if(offset) *offset = best_offset;
  return best_corr;
}

cluster* closest_cluster(cluster* cs, int n_segs, int lens[], double es[], int* offset, double* corr) {
  /* cs is a list of clusters, go through each one, get the one with highest abs correlation.
     return the cluster if found, return NULL if none satisfy the criteria.
     Also output the corresponding offset to *offset and correlation to *corr.
   */
  cluster* r=NULL;
  double best_corr=0, tmp_c=0;
  int best_offset=0, tmp_d=0;

  assert(cs != NULL);

  r = NULL;
  for(; cs!=NULL; cs=cs->next) {
    tmp_c = correlation_to_cluster(cs, n_segs,lens, es, &tmp_d);
    printf("Correlation to cluster %d is %f, offset is %d\n",
	   cs->id, tmp_c, tmp_d);
    if(fabs(tmp_c) > fabs(best_corr)) {
      r = cs;
      best_corr = tmp_c;
      best_offset = tmp_d;
    }
  }

  if(offset) *offset = best_offset;
  if(corr) *corr = best_corr;

  return r;
}

cluster* add_one_series(cluster* cs, series* se, int n_segs, int lens[], double es[], int max_delay, double eV_upper) {
  /* add es (of n_segs segments, with lengths in lens) to the list of
     clusters in cs, either to an existing and cloest cluster, or
     being a new cluster itself.

     eV_upper is the tolerated upper bound of error variance, used to determine if two series are close enough.
     If the series is close enough, es will become a new cluster.
     Return the modified linked list.
   */
  cluster* p=NULL;
  double corr=0, abs_corr=0, corr_threshold=0;
  int i=0, m=0, offset=0; /* m is total length */
  for(i=0, m=0; i<n_segs; i++) m += lens[i];

  /* at this stage, the segments are packed together with an offset.
     But the buffer of each cluster has room for each segment for shifting.
   */
  if(cs == NULL) {
    printf("New cluster 0\n");
    return single_cluster(0, m+n_segs*2*max_delay, m, es, max_delay, se, NULL);
  }

  p = closest_cluster(cs, n_segs,lens,es, &offset,&corr);

  if(p != NULL) {
    printf("Closest to cluster %d, corr: %f, offset: %d\n", p->id, corr, offset);
    /* determine the correlation threshold */
    corr_threshold = sqrt((1-eV_upper/(se->variance))*(1-eV_upper/(p->merged->variance)));
  }

  abs_corr = fabs(corr);
  if(p!=NULL && abs_corr >= corr_threshold) { /* existing one */
    printf("Add to cluster %d\n", p->id);
    /* merge_to_cluster(p, se, m, es, offset); */ /* ***** may need to change */
    /* add only the series, but the mean remains the first added series */
    p->n_series++;
    se->offset = offset;
    se->next = p->merged;
    p->merged = se;
    return cs;
  } else { /* new one */
    printf("New cluster %d\n", cs->id+1);
    return single_cluster(cs->id+1, m+n_segs*2*max_delay,m, es, max_delay, se, cs);
  }
}

cluster* find_clusters(array2d* d, int n_segs, int lens[],
		       int candidates[], int n_candidates,
		       int max_delay, double eV_upper, series s_buf[]) {
  /* assume d is in column major form, with n_segs segments, and the
     lengths of the segments are in lens.

     Consider the first n_candidates column indices in candidates for clustering.
     eV_upper is the tolerated upper bound of error variance, used to determine if two series are close enough.

     s_buf is a buffer of series for use here, and will be referred to
     in the clusters.

     Return a list of clusters that have at least two series.
   */
  cluster *cs=NULL, *p=NULL, *r_cs=NULL;
  int m=0, i=0;
  m = d->rows;

  for(i=0; i<n_candidates; i++) {
    s_buf[i].id = i;
    s_buf[i].idx = candidates[i];
    s_buf[i].offset = 0;
    s_buf[i].variance = var(m, get_col(d,candidates[i]));
    s_buf[i].next = NULL;
  }

  /* sequentially go through the columns to find clusters */
  for(i=0; i<n_candidates; i++) {
    printf("\n** Consider series %d, column %d\n", i, candidates[i]);
    cs = add_one_series(cs, s_buf+i, n_segs,lens, get_col(d,candidates[i]),
			max_delay, eV_upper);
  }
  /* print the clusters */
  printf("===== Resulting Clusters =====\n");
  for(p=cs; p!=NULL; p=p->next) {
    printf("Cluster %d: length %d, %d series\n", p->id, p->e - p->s, p->n_series);
    /*
    for(i=0; i<(p->e); i++) {
      printf("\t%d: %f\t%f\n", i, p->n[i], p->mean[i]);
    }
    printf("\n");
    */
  }
  /* filter out those with only one series */
  r_cs=NULL;
  while(cs != NULL) {
    p = cs;
    cs = cs->next;
    if(p->n_series <= 1) {
      p->next = NULL;
      free_clusters(p);
    } else {
      p->next = r_cs;
      r_cs = p;
    }
  }

  return r_cs;
}

int number_of_clusters(cluster* cs) {
  int n=0;
  for(n=0; cs!=NULL; n++, cs=cs->next)
    ;
  return n;
}

void re_estimate_mean(cluster* c, array2d* data, int n_segs, int lens[], int md, double buf[], double buf2[]) {
  /* use SVD to do rank-1 approximation to the overlapping parts of
     the series, to get the coefficients for the series.
     buf and buf2 are for tempoary use, and should have length c->cols.

     md is maximum delay.  output the minimum and maximum offset to
     buf[0] and buf2[0] respectively.
  */
  int i=0, j=0, k=0, ns=0, m=0;
  int minO=0, maxO=0, OL=0, L=0, sL=0, cL=0, mL=0, seg=0;
  int nr=0, nc=0, flip=0;
  series* se=NULL;
  gsl_matrix *A=NULL, *V=NULL;
  gsl_vector *singular_vals=NULL, *work=NULL;
  double *p=NULL, tmp=0;

  m = data->rows; /* should be equal to sum of lens */
  ns = c->n_series;

  /* assume c->s is the maximum delay allowed. And the original data
     in c is ignored. Now we re-fill it with mean, and each segment
     has length 2*max_delay + lens[i].

     Assume c has enough capacity, which is true if it is created by
     our other routines.
  */

  printf("Re-estimate mean for cluster %d\n", c->id);
  
  /* minO is min offset, maxO is max offset */
  minO = 2*md;
  maxO = 0;

  for(i=0, se=c->merged; se!=NULL; se=se->next, i++) {
    printf("series %d at offset %d\n", se->idx, se->offset);
    if(se->offset > maxO) maxO = se->offset;
    if(se->offset < minO) minO = se->offset;
  }

  /* also reset c->s and c->e based on the list of merged series */
  c->s = minO; c->e = maxO + m;
  printf("\n");

  /* now (maxO - minO) is the length of non-overlapping part at the
     front for a segment.

     lens[i] - (maxO - minO) is the length of the overlapping part, so
     the sum of overlapping parts of the segments is OL = m -
     n_segs*(maxO - minO), assuming each segment is long enough.

     buf[] is used for the mean of the series.
   */
  for(OL=0, seg=0; seg<n_segs; seg++) {
    L = lens[seg] - (maxO - minO);
    if(L > 0) OL += L;
  }

  /* SVD in GSL assumes the number of rows >= number of columns */
  if(OL >= ns) {
    nr = OL;
    nc = ns;
    flip = 0;
  } else {
    nr = ns;
    nc = OL;
    flip = 1;
  }

  A = gsl_matrix_alloc(nr, nc);
  V = gsl_matrix_alloc(nc, nc);
  singular_vals = gsl_vector_alloc(nc);
  work = gsl_vector_alloc(nc);

  /* fill in A */
  for(i=0, se=c->merged; se!=NULL; se=se->next, i++) {
    /* segment by segment */
    p = get_col(data, se->idx);
    cL = 0;
    for(seg=0; seg<n_segs; seg++) {
      L = lens[seg] - (maxO - minO); /* overlapping length of this segment */
      if(L > 0) { /* ignored if segment too short */
	buf[i] = mean(lens[seg], p);
	for(j=0, k=maxO-(se->offset); j<L; j++, k++) {
	  if(flip) {
	    gsl_matrix_set(A, i,cL+j, p[k] - buf[i]);
	  } else {
	    gsl_matrix_set(A, cL+j,i, p[k] - buf[i]);
	  }
	}
	cL += L; /* cumulative overlapping length */
      }
      p += lens[seg];
    }
  }

  gsl_linalg_SV_decomp(A, V, singular_vals, work);

  /* only need rank-1 approximation, associate the singular value to
     the coefficients */
  memset(c->mean, 0, sizeof(double)*c->capacity);
  memset(c->n, 0, sizeof(double)*c->capacity);
  /* get the overlapping part and coefficients in buf2 */
  tmp = gsl_vector_get(singular_vals, 0);

  /* the coefficients */
  if(flip) {
    for(j=0; j<ns; j++) {
      buf2[j] = tmp*gsl_matrix_get(A, j,0);
    }
  } else {
    for(j=0; j<ns; j++) {
      buf2[j] = tmp*gsl_matrix_get(V, j,0);
    }
  }
  /* for completeness, write back the sum of square coefficients for the
     overlapping part */
  for(tmp=0, j=0; j<ns; j++) tmp += buf2[j]*buf2[j];

  /* mean for overlapping part */
  cL = 0;
  sL = 0;
  mL = 0;
  for(seg=0; seg<n_segs; seg++) {
    L = lens[seg] - (maxO - minO);
    if(L > 0) {
      if(flip) { /* the estimate is in V */
	for(j=0, k=maxO; j<L; j++, k++) {
	  c->mean[mL+k] = gsl_matrix_get(V, cL+j,0);
	  c->n[mL+k] = tmp;
	}
      } else { /* the estimate is in U, which is A */
	for(j=0, k=maxO; j<L; j++, k++) {
	  c->mean[mL+k] = gsl_matrix_get(A, cL+j,0);
	  c->n[mL+k] = tmp;
	}
      }
      cL += L;
    }

    /* try our best to estimate the non-overlapping part */
    for(i=0, se=c->merged; se!=NULL; se=se->next, i++) {
      p = get_col(data, se->idx) + sL;
      buf[i] = mean(lens[seg], p);
      /* accumulate prefix */
      for(j=se->offset, k=0; j<maxO; j++, k++) {
	c->mean[mL+j] += buf2[i]*(p[k] - buf[i]);
	c->n[mL+j] += buf2[i]*buf2[i];
      }
      /* accumulate suffix */
      for(j=se->offset+lens[seg]-1, k=lens[seg]-1; j>=(maxO+L); j--, k--) {
	c->mean[mL+j] += buf2[i]*(p[k] - buf[i]);
	c->n[mL+j] += buf2[i]*buf2[i];
      }
    }

    /* estimate prefix and suffix */
    for(j=minO; j<maxO; j++) {
      c->mean[mL+j] /= c->n[mL+j];
    }
    /* estimate suffix */
    for(j=minO+lens[seg], k=minO; k<maxO; j++, k++) {
      c->mean[mL+j] /= c->n[mL+j];
    }

    /* after re-estimation */
    for(j=minO; j<(maxO+lens[seg]); j++) {
      printf("\t%d: %f\t%f\n", mL+j, c->n[mL+j], c->mean[mL+j]);
    }

    /* next segment */
    sL += lens[seg];
    mL += lens[seg]+2*md; /* next part for c->mean */
  }
  printf("\n");

  buf[0] = minO;
  buf2[0] = maxO;

  /* clean up */
  gsl_matrix_free(A);
  gsl_matrix_free(V);
  gsl_vector_free(singular_vals);
  gsl_vector_free(work);
}

array2d* augment_data_with_clusters(array2d* data, int n_segs, int lens[],
				    int candidates[], int n_candidates, int n_candidates_parents,
				    int max_delay, double eV_upper, char** out_forbidden) {
  /* data is assumed to be in column major form, with n_segs segments,
     and the lengths of the segments are in lens.

     Consider the first n_candidates column indices in candidates for clustering.
     return a possibly augmented data (still in column major form),
     that adds the cluster centers as columns.
     When augmenting the data, the n_candidates+n_candidates_parents column indices in candidates are also included.

     eV_upper is the tolerated upper bound of error variance, used to determine if two series are close enough.

     If out_forbidden is non-NULL, *out_forbidden will be written a newly allocated array
       where out_forbidden[i*ng + j] is 1 if i->j is forbidden,
       where ng is the number of columns in augmented data.
       The array in *out_forbidden should be freed after use.
  */
  int n_cs=0, m=0, n=0, i=0, s=0;
  array2d* a_data=NULL;
  cluster *cs=NULL, *p=NULL;
  series* s_buf=NULL;
  double *buf=NULL;
  int seg=0, cL=0;
  double *u=NULL, *v=NULL;

  m = data->rows;
  n = n_candidates + n_candidates_parents;

  s_buf = Malloc(n, series, "augment_data_with_clusters()");
  buf = Malloc(2*n, double, "augment_data_with_clusters()");

  cs = find_clusters(data, n_segs,lens, candidates, n_candidates,
		     max_delay, eV_upper, s_buf);
  n_cs = number_of_clusters(cs);

  a_data = new_array2d(m, n+n_cs); /* n_cs new columns */

  /* copy the base candidate and parents columns */
  for(i=0; i<n; i++) {
    memcpy(get_col(a_data,i), get_col(data,candidates[i]), sizeof(double)*m);
  }

  /* the new columns */
  for(p=cs, i=n; p!=NULL; p=p->next, i++) {
    printf("Add cluster %d as column %d\n", p->id, i);
    memset(get_col(a_data, i), 0, sizeof(double)*m);

    /* re-estimate the cluster mean */
    re_estimate_mean(p, data, n_segs,lens, max_delay, buf, buf+n);
    /* the minimum and maximum offset are in buf[0] and buf[n] respectively */

    /* for the time being, simply use the last part of the merged
       series for each segment, so that it is before the merged
       series, but with at least one more delay.
     */
    u = get_col(a_data, i);
    v = p->mean;
    for(seg=0, cL=0; seg<n_segs; seg++) {
      s = buf[n] + 1; /* shift at least one */
      printf("Taking from %d to %d, length %d.\n", cL+s, cL+s+lens[seg]-1, lens[seg]-1);
      memcpy(u, v + s, sizeof(double)*(lens[seg]-1));

      /* next segment */
      cL += lens[seg] + 2*max_delay;
      u += lens[seg];
      v += lens[seg] + 2*max_delay;
    }
  }

  /* record the forbidden links, because they are likely children of
     an introduced hidden common node */
  if(out_forbidden) {
    {
      char* f=NULL;
      int ng = n + n_cs;
      series *a=NULL, *b=NULL;
      f = Malloc(ng*ng, char, "augment_data_with_clusters() out_forbidden");
      memset(f, 0, ng*ng*sizeof(char));
      for(p=cs; p!=NULL; p=p->next) {
	/* brute force */
	for(a=p->merged; a!=NULL; a=a->next) {
	  for(b=a->next; b!=NULL; b=b->next) {
	    f[ng*(a->id) + (b->id)] = 1; /* a -> b forbidden */
	    f[ng*(b->id) + (a->id)] = 1; /* b -> a forbidden */
	  }
	}
      }
      *out_forbidden = f;
    }
  }

  /* done, cleanup */
  free_clusters(cs);
  free(s_buf);

  return a_data;
}

/***************************************************************/
double fit_to_parents(array2d* d, int n_segs, int lens[],
		      int ig, int n_parents, link* parents[]) {
  /* d is n by g array, the n_segs segments of expression data in column major format,
     where rows are time points, and columns are genes.
     The first lens[0] rows are segment 1,
     the next lens[1] rows are segment 2, and so on.

     calculate and return the error variance for gene ig, after
     regressing its expression to its parents.

     also update the coefficient in the parent links.
  */
  int i=0, j=0, k=0, u=0, md=0;
  int m=0, n=0, L=0;
  double chisq = 0;
  gsl_matrix *X=NULL, *cov=NULL;
  gsl_vector *y=NULL, *c=NULL;
  gsl_multifit_linear_workspace *work=NULL;

  /* determine the maximum delay */
  for(i=0; i<n_parents; i++) {
    if(parents[i]->delay > md) md = parents[i]->delay;
  }
  /* determine the length after shifting, each segment is deducted by md */
  L = 0;
  for(i=0; i<n_segs; i++) {
    m = lens[i] - md;
    if(m > 0) L += m;
  }
  /* collect the shifted expression, prepare y = Xc */
  n = n_parents + 1; /* fit also the offset */

  /* debug 
  printf("fit_to_parents(): ig: %d\tmd: %d\tL: %d\tn: %d\n", ig,md, L, n);
   */

  X = gsl_matrix_alloc(L, n);
  y = gsl_vector_alloc(L);

  c = gsl_vector_alloc(n);
  cov = gsl_matrix_alloc(n, n);

  j=0;
  k=0;
  for(i=0; i<n_segs; i++) {
    if(lens[i] > md) {
      for(m=md; m<lens[i]; m++, k++) {
	gsl_matrix_set(X, k, 0, 1.0);
	for(u=0; u<n_parents; u++) {
	  gsl_matrix_set(X, k, u+1, caref(d, j+m - parents[u]->delay,
					  parents[u]->from));
	}
	gsl_vector_set(y, k, caref(d, j+m, ig));
      }
    }
    j += lens[i];
  }

  /* regression */
  work = gsl_multifit_linear_alloc (L, n);
  gsl_multifit_linear (X, y, c, cov, &chisq, work);
  gsl_multifit_linear_free (work);
  /* chisq contains the sum of square of residual error */
  /* get the coefficients, ignore the offset */
  for(u=0; u<n_parents; u++) {
    parents[u]->test_value = gsl_vector_get(c,u+1);
    /* debug 
    printf("** Coef for parent %d: %g\n", parents[u]->from, parents[u]->test_value);
     */
  }
  /* debug 
  printf("** Intercept term: %g\n", gsl_vector_get(c,0));
  printf("** chisq: %g\n", chisq);
   */

  /* clean up */
  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  /* done, return the estimate of variance of residual error */
  return chisq/(L-1);
}

int decreasing_double(const void* a, const void* b) {
  double ax=0, ay=0;
  ax = *(double*)a;
  ay = *(double*)b;
  if(ax > ay) return -1;
  if(ax < ay) return 1;
  return 0;
}
static int translate_index(int idx, int n, int candidates[], int ng) {
  /* idx is to be translated. candidates has length n.
     ng is the number of observed genes.
   */
  if(idx < n) return candidates[idx];
  return ng + idx - n;
}
int hcc_clinde(array2d* d, int n_segs, int lens[],
	       double st1, double st2, int max_delay, int max_n,
	       int method, int pruning, 
	       int is_one_delay, int is_no_dup, int is_keep_zero_delay,
	       double eV, double eV_tolerance, link** out_links) {
  /* If eV > 0, it is the expected variance of the error terms, and if
     the error term of a gene is above eV*(1+eV_tolerance), it will be
     candidate of having hidden common node.

     If eV <= 0, the expected variance of the error terms will be
     estimated as the median.

     Other parameters are similar to clinde() above.

     To try to handle hidden common nodes in GRN.
     Main idea:
     First use CLINDE to infer an GRN, then check the error levels of each gene.

     Those genes with error levels very different from expected, and
     their parents which have no parents are considered candidates for
     having hidden common nodes and the links with them are cut.

     Then cluster the series of the candidates to estimate the hidden
     common node(s). And then run CLINDE on these candidates and
     estimated hidden common nodes, with some prunings.

     Finally merge the two GRNs to give the final GRN.
   */
  int ng=0;
  int n_links=0, i=0, j=0, k=0;
  link* init_links=NULL;
  int* n_to=NULL;
  link*** to=NULL;
  link** buf_p_links=NULL;
  double *err_var=NULL, err_var_threshold=0;
  int n_candidates=0, n_candidates_parents=0, *candidates=NULL, *is_in=NULL;
  int ns_links=0; /* for the sub part involving hidden common node */
  link* s_links=NULL;
  int n_final_links=0;

  ng = d->cols;

  /* get the inititial network */
  printf("==== Initial GRN\n");
  n_links = clinde(d, n_segs, lens, st1, st2, max_delay, max_n,
		   method, pruning, is_one_delay, is_no_dup, is_keep_zero_delay,
		   NULL, ng, &init_links);
  printf("==== End Initial GRN\n");

  /* allocate n_links+1 in case n_links is 0 */
  /* organize the links in a different way */
  n_to = Malloc(ng, int, "hcc_clinde()");
  buf_p_links = Malloc(n_links+1, link*, "hcc_clinde()");
  to = Malloc(ng, link**, "hcc_clinde()");

  memset(n_to, 0, ng*sizeof(int));
  for(i=0; i<n_links; i++) { /* count the links into each gene */
    k = init_links[i].to;
    n_to[k]++;
  }
  for(i=0, j=0; i<ng; i++) { /* assign buffer to each gene */
    to[i] = buf_p_links+j;
    j += n_to[i];
    n_to[i] = 0; /* to use as index temporarily */
  }
  for(i=0; i<n_links; i++) { /* fill in the link* */
    k = init_links[i].to;
    to[k][n_to[k]] = init_links+i;
    n_to[k]++;
  }

  /* debug 
  {
    link* tmp=NULL;
    printf("**** Parents of genes\n");
    for(i=0; i<ng; i++) {
      printf("** %d Parents of gene %d:\n", n_to[i], i);
      for(j=0; j<n_to[i]; j++) {
	tmp = to[i][j];
	printf("\tFrom: %d\tTo: %d\tDelay: %d\tCoef: %g\n",
	       tmp->from, tmp->to, tmp->delay, tmp->test_value);
      }
    }
    printf("**** End Parents of genes\n\n");
  }
    */

  /* get error level of each gene */
  err_var = Malloc(2*ng, double, "hcc_clinde()");
  printf("=== Estimated Error Variance\n");
  for(i=0; i<ng; i++) {
    err_var[i] = fit_to_parents(d, n_segs,lens, i,n_to[i],to[i]);
    err_var[i+ng] = err_var[i];
    printf("Gene %d: %g\n", i, err_var[i]);
  }
  printf("=== End Estimated Error Variance\n");

  /* threshold for error variance */
  if(eV <= 0) { /* use the median as estimate */
    qsort(err_var+ng, ng, sizeof(double), decreasing_double);
    eV = err_var[ng+(ng/2)];
  }
  err_var_threshold = eV*(1+eV_tolerance);
  printf("Error variance threshold: %g * ( 1 + %g ): %g\n",
	 eV, eV_tolerance, err_var_threshold);

  /* determine candidate genes as having hidden common node.
     First mark them, then record their index.
   */
  n_candidates = 0;
  candidates = Malloc(2*ng, int, "hcc_clinde()");
  memset(candidates, 0, 2*ng*sizeof(int));
  is_in = candidates + ng;
  for(i=0; i<ng; i++) {
    if(err_var[i] > err_var_threshold) {
      is_in[i] = 1;
      n_candidates++;
    }
  }
  n_candidates_parents = 0;
  /* also the parents of the candidates (which are not candidates) are
     also included, but restricted */
  for(i=0; i<ng; i++) {
    if(is_in[i] == 1) {
      for(j=n_to[i]-1; j>=0; j--) {
	k = to[i][j]->from;
	if(is_in[k]==0) {
	  is_in[k] = 2;
	  n_candidates_parents++;
	}
      }
    }
  }

  j=0;
  for(i=0; i<ng; i++) {if(is_in[i]==1) candidates[j++] = i;}
  for(i=0; i<ng; i++) {if(is_in[i]==2) candidates[j++] = i;}

  printf("==== Unexpected error variance in %d gene(s):", n_candidates);
  for(i=0; i<n_candidates; i++) printf(" %d", candidates[i]);
  printf("\n");
  printf("==== Their parents, %d gene(s):", n_candidates_parents);
  for(i=0; i<n_candidates_parents; i++)
    printf(" %d", candidates[i+n_candidates]);
  printf("\n");

  report_stop_watch();
  /* do the hidden common node part, if applicable */
  if(n_candidates > 0) {
    {
      char* forbidden=NULL;
      array2d* a_data=NULL;
      int n=0, n_ng=0;

      n = n_candidates + n_candidates_parents;
      printf("====== Clustering and Augmenting the Data ======\n");
      a_data = augment_data_with_clusters(d, n_segs,lens,
					  candidates, n_candidates,
					  n_candidates_parents, max_delay,
					  eV*(1+eV_tolerance), &forbidden);
      report_stop_watch();
      n_ng = a_data->cols;
      if(n_ng > n) { /* introduced hidden common node, re-learn the subpart */
	/* also restrict the parents to be only parents.
	   note that forbidden is relative to candidates.
	*/
	for(i=n-1; i>=n_candidates; i--) {
	  for(j=0; j<n_ng; j++) {
	    forbidden[j*n_ng + i] = 1; /* j -> i forbidden */
	  }
	}

	printf("==== GRN for potential hidden common causes\n");
	ns_links = clinde(a_data, n_segs, lens, st1, st2, max_delay, max_n,
			  method, pruning, is_one_delay, is_no_dup,
			  is_keep_zero_delay, forbidden, n, &s_links);
	printf("==== End GRN for potential hidden common causes\n");
	report_stop_watch();

	/* remove the parents of the candidates */
	for(i=0; i<n_candidates; i++) {
	  n_to[candidates[i]] = 0;
	}
	/* translate the candidates subpart back to the index to d,
	   and introduced hidden nodes to new index.
	 */
	for(i=0; i<ns_links; i++) {
	  s_links[i].from = translate_index(s_links[i].from, n, candidates, ng);
	  s_links[i].to = translate_index(s_links[i].to, n, candidates, ng);
	}
      } else {
	printf("==== No candidate hidden common node identified.\n");
      }
      /* small clean up */
      free_array2d(a_data);
      free(forbidden);
    }
  }

  /* output */
  printf("==== Final GRN\n");
  n_final_links = ns_links;
  for(i=0; i<ng; i++) {
    for(j=0; j<n_to[i]; j++)
      print_link(stdout, to[i][j]);
    n_final_links += n_to[i];
  }
  printf("=== possibly with hidden common node\n");
  for(i=0; i<ns_links; i++) {
    print_link(stdout, s_links+i);
  }
  printf("==== End Final GRN\n");

  if(out_links!=NULL && n_final_links>0) {
    {
      link* outL = NULL;
      int k=0;
      outL = Malloc(n_final_links, link, "hcc_clinde()");
      for(i=0; i<ng; i++) {
	for(j=0; j<n_to[i]; j++)
	  outL[k++] = *(to[i][j]);
      }
      for(i=0; i<ns_links; i++) {
	outL[k++] = s_links[i];
      }
      *out_links = outL;
    }
  }

  /* clean up */
  free(n_to);
  free(buf_p_links);
  free(to);
  free(init_links);
  free(err_var);
  free(candidates);

  /* done */
  return n_final_links;
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

void normalize_data(array2d* d, int n_segs, int lens[]) {
  /* Assume d is column major, each row is time point, each column is gene.
     Also, the rows are divided into n_segs segments, each with lengths in lens.

     Normalize each segment of each gene in place, such that the mean
     is 0, and s.d. is 1 (unless the s.d. is too close to 0).
   */
  int ng=0, i=0, seg=0, sL=0;
  double *p=NULL, g_mean=0, g_sd=0;

  ng = d->cols;

  printf("==== Start Normalization ====\n");
  for(seg=0, sL=0; seg<n_segs; seg++) {
    printf("Segment %d:\n", seg);
    for(i=0; i<ng; i++) {
      p = get_col(d, i);
      g_mean = 0;
      g_sd = 0;
      normalize(lens[seg], p + sL, &g_mean,&g_sd);
      printf("\tgene %d: mean: %g\t\ts.d.: %g\n", i, g_mean,g_sd);
    }
    sL += lens[seg];
  }
  printf("==== End Normalization ====\n");
}

/***************************************************************/
#define DEFAULT_ST1        2
#define DEFAULT_ST2        2
#define DEFAULT_MAX_DELAY  4
#define DEFAULT_MAX_N      4
#define DEFAULT_METHOD     "pcor"
#define DEFAULT_PRUNING    "all"
#define DEFAULT_CTL        0.75
#define DEFAULT_CTU        0.86
#define DEFAULT_EV         1
#define DEFAULT_EV_TOLERANCE 0.1

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
  {"-keep.zero.delay", FLAG,  0,  "To keep links with zero delay after stage 2 and before remove dups. Default false.\n", 0,0,  0,  0.0,    NULL},
  {"-one.delay", FLAG,   0,  "To keep only one delay (the one with best score, smallest delay) for each link after stage 1. Default false.\n", 0,0,  0,  0.0,    NULL},
  {"-no.dup",    FLAG,   0,  "Remove duplicate links. To keep only one delay (the one with best score, smallest delay) for each link after stage 2. Default false.\n", 0,0,  0,  0.0,    NULL},
  {"-normalize", FLAG,   0,  "Center and normalize the expressions such that for each gene, the mean is 0, and s.d. is 1. Default false.\n", 0,0,  0,  0.0,    NULL},
  {"-ctL",       REAL,   0,  "Lower Correlation threshold for initial clustering. If ctL <= abs correlation < ctU, consider the two nodes for having hidden common node. Default " xstr(DEFAULT_CTL) ".\n", 0,0,  0,  DEFAULT_CTL,    NULL},
  {"-ctU",       REAL,   0,  "Upper Correlation threshold for initial clustering. If ctU <= abs correlation, assume a direct link but NOT hidden common node. Default " xstr(DEFAULT_CTU) ".\n", 0,0,  0,  DEFAULT_CTU,    NULL},
  {"-eV",        REAL,   0,  "Expected Error Variance if > 0. If <= 0, the expected error variance will be estimated to be the median of the error variances of the initial GRN. If error variance of a gene is > eV*(1+eV.tolerance), that gene will be considered for having hidden common node. Default " xstr(DEFAULT_EV) ".\n", 0,0,  0,  DEFAULT_EV,    NULL},
  {"-eV.tolerance",        REAL,   0,  "If error variance of a gene is > eV*(1+eV.tolerance), that gene will be considered for having hidden common node. Default " xstr(DEFAULT_EV_TOLERANCE) ".\n", 0,0,  0,  DEFAULT_EV_TOLERANCE,    NULL},

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
  double ctL=DEFAULT_CTL, ctU=DEFAULT_CTU;
  double eV=DEFAULT_EV, eV_tolerance=DEFAULT_EV_TOLERANCE;
  int max_delay=DEFAULT_MAX_DELAY, max_n=DEFAULT_MAX_N;
  char* s_method=DEFAULT_METHOD;
  char* s_pruning=DEFAULT_PRUNING;
  int method=0, pruning=0;
  int is_one_delay=0, is_no_dup=0, is_keep_zero_delay=0, is_normalize=0;

  int s_idx=0, e_idx=-1;
  int n_data=0, *data_lens=NULL;
  array2d *data=NULL;

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
  ctL = get_real_option_value("-ctL", all_options, ctL);
  ctU = get_real_option_value("-ctU", all_options, ctU);
  eV = get_real_option_value("-eV", all_options, eV);
  eV_tolerance = get_real_option_value("-eV.tolerance", all_options, eV_tolerance);
  max_delay = get_int_option_value("-max.delay", all_options, max_delay);
  max_n = get_int_option_value("-max.n", all_options, max_n);
  s_method = get_str_option_value("-method", all_options, s_method);
  s_pruning = get_str_option_value("-pruning", all_options, s_pruning);
  is_one_delay = get_flag_option_value("-one.delay", all_options, is_one_delay);
  is_no_dup = get_flag_option_value("-no.dup", all_options, is_no_dup);
  is_keep_zero_delay = get_flag_option_value("-keep.zero.delay", all_options, is_keep_zero_delay);
  is_normalize = get_flag_option_value("-normalize", all_options, is_normalize);

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
  start_stop_watch();
  printf("====== CLINDE Grn Inference ======\n");

  printf("Stage 1 score threshold: %f\n", st1);
  printf("Stage 2 score threshold: %f\n", st2);
  printf("Lower Correlation threshold for clustering: %f\n", ctL);
  printf("Upper Correlation threshold for clustering: %f\n", ctU);
  if(eV > 0)
    printf("Expected Error Variance: %f\n", eV);
  else
    printf("Expected Error Variance: Estimated\n");
  printf("Error Variance Tolerance: %f\n", eV_tolerance);
  printf("Upper bound (exclusive) of delays to test: %d\n", max_delay);
  printf("Maximum Parents to condition on in stage 2: %d\n", max_n);
  printf("Method: %s (%s)\n", method_names[method], method_long_names[method]);
  printf("Pruning: %s\n", pruning_names[pruning]);
  if(is_one_delay)
    printf("Retain only one delay for each link after stage 1.\n");
  if(is_keep_zero_delay)
    printf("To keep links with zero delay after stage 2 and before removing dups.\n");
  if(is_no_dup)
    printf("Retain only one delay for each link, remove duplicate links after stage 2.\n");

  /*********/

  data_lens = Malloc(n_data, int, "data_lens");
  data = read_n_data(n_data, argv+s_idx+1, data_lens);
  if(data == NULL) return 1;
  report_stop_watch();

  /* debug 
  printf("data_lens"); pr_ints(n_data, data_lens);
  */
  if(is_normalize) {
    printf("To normalize expression of each gene to mean 0, s.d. 1.\n");
    normalize_data(data, n_data,data_lens);
  }

  /*********/
  hcc_clinde(data, n_data,data_lens, st1,st2,max_delay,max_n,
	     method,pruning, is_one_delay, is_no_dup, is_keep_zero_delay,
	     eV, eV_tolerance, NULL);

  report_stop_watch();
  /*********/
  free_array2d(data);
  free(data_lens);

  return 0;
}

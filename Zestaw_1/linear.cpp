#include <iostream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <sys/resource.h>
#include <string.h>
#include "linear.h"

double getTime(struct rusage *ru0, struct rusage *ru1){
  double utime = 0, stime = 0, ttime = 0;
  /* Get time in seconds */
  utime = (double) ru1->ru_utime.tv_sec
	+ 1.e-3 * (double) ru1->ru_utime.tv_usec
	- ru0->ru_utime.tv_sec
	- 1.e-6 * (double) ru0->ru_utime.tv_usec;
  stime = (double) ru1->ru_stime.tv_sec
	+ 1.e-6 * (double) ru1->ru_stime.tv_usec
	- ru0->ru_stime.tv_sec
	- 1.e-6 * (double) ru0->ru_stime.tv_usec;
  ttime = stime + utime;
  return ttime;
}

gsl_matrix* generate_random_matrix(int n, int range_max, int range_min){
  /* Generate A(nxn) matrix and vector of absolute terms. Allocate the memory for
  m_check - original matrix to check the output after calculations */
  gsl_matrix * m = gsl_matrix_alloc (n, n);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      gsl_matrix_set (m, i, j, rand()%range_max+range_min);
    }
  }
  return m;
}

gsl_vector* generate_random_vector(int n, int range_max, int range_min){
  gsl_vector * m0 = gsl_vector_alloc(n);
  for (int i = 0; i < n; i++){
    gsl_vector_set (m0, i, rand()%range_max+range_min);
  }
  return m0;
}

gsl_matrix* create_copy(gsl_matrix* matrix, int n){
  gsl_matrix * m_copy = gsl_matrix_alloc (n, n);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      gsl_matrix_set (m_copy, i, j, gsl_matrix_get (matrix, i, j));
    }
  }
  return m_copy;
}

bool check_correctitude(int n, gsl_vector* x, gsl_vector* free_column, gsl_matrix* matrix){
  /* Check if the result is correct */
  /* matrix * x - free_column ?<= eps */
  double eps = 1.e-6;
  gsl_matrix * m1 = rewrite_vector_to_matrix(x, n);
  gsl_matrix * b_vector = rewrite_vector_to_matrix(free_column, n);
  gsl_matrix * output = gsl_matrix_calloc (n, 1);
  /* Mul m_check with m and save it into output */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, matrix,  m1,
                  0.0, output);
  /* Substract b_vector from output */
  gsl_matrix_sub (output, b_vector);
  for (int i = 0; i < n; i++){
       if(gsl_matrix_get (output, i, 0) > eps)
          return false;
  }
  return true;
}

gsl_matrix* rewrite_vector_to_matrix(gsl_vector* vector, int n){
  gsl_matrix * m = gsl_matrix_alloc (n, 1);
  for (int i = 0; i < n; i++){
    gsl_matrix_set (m, i, 0, gsl_vector_get(vector, i));
  }
  return m;
}

void matrix_decomposite(gsl_matrix* matrix, int &s, gsl_permutation *&p, int n){
  p = gsl_permutation_alloc (n);
  gsl_linalg_LU_decomp (matrix, p, &s);
}

void matrix_lu_solve(gsl_matrix* matrix, gsl_vector *&x, gsl_permutation *&p, gsl_vector* free_column, int n){
  x = gsl_vector_alloc (n);
  gsl_linalg_LU_solve (matrix, p, free_column, x);
}

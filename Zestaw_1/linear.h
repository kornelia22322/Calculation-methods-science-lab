#ifndef _LINEAR
#define _LINEAR

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <sys/resource.h>
#include <ctime>

double getTime(struct rusage *ru0, struct rusage *ru1);

gsl_matrix* generate_random_matrix(int n, int range_max, int range_min);

gsl_vector* generate_random_vector(int n, int range_max, int range_min);

gsl_matrix* create_copy(gsl_matrix* matrix, int n);

bool check_correctitude(int n, gsl_vector* x, gsl_vector* free_column, gsl_matrix* matrix);

gsl_matrix* rewrite_vector_to_matrix(gsl_vector* vector, int n);

void matrix_decomposite(gsl_matrix* matrix, int &s, gsl_permutation *&p, int n);

void matrix_lu_solve(gsl_matrix* matrix, gsl_vector *&x, gsl_permutation *&p, gsl_vector* free_column, int n);

#endif

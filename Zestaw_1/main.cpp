#include "linear.h"
#include <iostream>
#include <sys/resource.h>

using namespace std;

int check_values[] = {10, 50, 100, 200, 300, 400, 500, 600, 700, 800, 1000};

int main (void) {
  int n;
  struct rusage t0, t1, t2;
  FILE *file, *file2;
  file = fopen("plotdata.dat", "w");
  file2 = fopen("plotdata2.dat", "w");

  for(int j=0;j<10;j++){
    n = check_values[j];

    /* Generate A(nxn) matrix and vector of absolute terms. Allocate the memory for
    m_check - original matrix to check the output after calculations */
    gsl_matrix * m = generate_random_matrix(n,1,100);
    gsl_matrix * m_check = create_copy(m,n);
    gsl_vector * m0 = generate_random_vector(n,1,100);

    gsl_vector *x;
    int s;
    gsl_permutation * p;

    getrusage(RUSAGE_SELF, &t0);
    matrix_decomposite (m, s, p, n);

    getrusage(RUSAGE_SELF, &t1);
    matrix_lu_solve (m, x, p, m0, n);

    getrusage(RUSAGE_SELF, &t2);

    gsl_vector_fprintf (stdout, x, "%g");
    fprintf (file,  "%d %lf\n", n, getTime(&t0, &t1));
    fprintf (file2,  "%d %lf\n", n, getTime(&t1, &t2));

    gsl_matrix_free (m);
    gsl_matrix_free (m_check);
    gsl_permutation_free (p);
    gsl_vector_free (x);
  }

  return 0;
}

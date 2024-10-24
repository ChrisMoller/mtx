#include "../mtx_config.h"

#include<cmath>
#include<complex>
#include<iostream>
#include<fstream>
#include<string>

#include <alloca.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_URL
#undef PACKAGE_VERSION
#undef VERSION

#include "Shape.hh"
#include "Matrix.hh"

/***
    a←4 4 ⍴ ¯1 1 ¯1 1 ¯8 4 ¯2 1 27 9 3 1 64 16 4 1
 ***/

Matrix
getEigenvectors (Matrix *mtx)
{
  Matrix res (mtx->rows (), mtx->cols ());

  double *data =
    (double *)alloca (mtx->rows () * mtx->cols () * sizeof(double));
  int i = 0;
  for (int j = 0; j < mtx->rows (); j++) {
    for (int k = 0; k < mtx->cols (); k++, i++) {
      data[i] = mtx->val (j, k).real ();
    }
  }
  
  gsl_matrix_view m = gsl_matrix_view_array (data, mtx->rows (), mtx->cols ());

  gsl_vector_complex *eval = gsl_vector_complex_alloc (mtx->rows ());
  gsl_matrix_complex *evec =
    gsl_matrix_complex_alloc (mtx->rows (), mtx->cols ());

  gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc (mtx->rows ());

  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);

  gsl_eigen_nonsymmv_free (w);

  gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

  int k = 0;
  
  for (int i = 0; i < mtx->cols (); i++) {
    gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
    
    for (int j = 0; j < mtx->rows (); j++) {
      gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
      complex<double> v (GSL_REAL(z), GSL_IMAG(z));
      res.val (i, j, v);
    }
  }

  return res;
}

vector<complex<double>>
getEigenvalues (Matrix *mtx)
{
  vector<complex<double>> res (mtx->cols ());

  double *data =
    (double *)alloca (mtx->rows () * mtx->cols () * sizeof(double));
  int i = 0;
  for (int j = 0; j < mtx->rows (); j++) {
    for (int k = 0; k < mtx->cols (); k++, i++) {
      data[i] = mtx->val (j, k).real ();
    }
  }
  
  gsl_matrix_view m = gsl_matrix_view_array (data, mtx->rows (), mtx->cols ());

  gsl_vector_complex *eval = gsl_vector_complex_alloc (mtx->rows ());
  gsl_matrix_complex *evec =
    gsl_matrix_complex_alloc (mtx->rows (), mtx->cols ());

  gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc (mtx->rows ());

  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);

  gsl_eigen_nonsymmv_free (w);

  gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

  int k = 0;
  
  for (int i = 0; i < mtx->cols (); i++) {
    gsl_complex eval_i = gsl_vector_complex_get (eval, i);

    complex<double> v (GSL_REAL(eval_i), GSL_IMAG(eval_i));
    res[i] = v;
  }
  return res;
}



/*
    This file is part of GNU APL, a free implementation of the
    ISO/IEC Standard 13751, "Programming Language APL, Extended"

    mtx Copyright (C) 2024  Dr. C. H. L. Moller

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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


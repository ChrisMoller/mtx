#include "../mtx_config.h"

#include<cmath>
#include<complex>
#include<iostream>
#include<fstream>
#include<string>

#ifdef HAVE_EIGEN
#include <eigen3/Eigen/Dense>
#endif

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

#ifdef HAVE_EIGEN
Matrix
getEigenvectors (Matrix *mtx)
{
  Matrix res (mtx->rows (), mtx->cols ());
  Eigen::MatrixXcd mtxE (mtx->rows (), mtx->cols ());
  for (int j = 0; j < mtx->rows (); j++) {
    for (int k = 0; k < mtx->cols (); k++) {
      mtxE (j, k) = mtx->val (j, k);
    }
  }

  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces(mtxE);

  if (ces.info () != Eigen::Success)
    DOMAIN_ERROR;

  for (int l = 0; l < mtx->cols (); l++) {
    for (int m = 0; m < mtx->rows (); m++) {
      res.val (l, m, (ces.eigenvectors ().col (l))(m));
    }
  }
  return res;
}
#endif


#ifdef HAVE_EIGEN
vector<complex<double>>
getEigenvalues (Matrix *mtx)
{
  vector<complex<double>> res (mtx->cols ());
  Eigen::MatrixXcd mtxE (mtx->rows (), mtx->cols ());
  for (int j = 0; j < mtx->rows (); j++) {
    for (int k = 0; k < mtx->cols (); k++) {
      mtxE (j, k) = mtx->val (j, k);
    }
  }

  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces(mtxE, false);

  if (ces.info () != Eigen::Success)
    DOMAIN_ERROR;

  for (int l = 0; l < mtx->rows (); l++)
    res[l] = ces.eigenvalues ()(l);

  return res;
}
#endif


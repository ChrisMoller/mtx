#pragma once

#include<complex>
#include "Matrix.hh"
#include "eigens.hh"

vector<complex<double>> getEigenvalues (Matrix *mtx);
Matrix getEigenvectors (Matrix *mtx);

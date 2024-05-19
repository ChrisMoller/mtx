/*
    This file is part of GNU APL, a free implementation of the
    ISO/IEC Standard 13751, "Programming Language APL, Extended"

    Copyright (C) 2008-2013  Dr. Jürgen Sauermann
    edif Copyright (C) 2024  Dr. C. H. L. Moller

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

/***
    det←{mtx['d']⍵}
    cross←{⍺ mtx ⍵}
 ***/

#include<cmath>
#include<complex>
#include<iostream>
#include<fstream>
#include<string>

#include <dlfcn.h>

#include <eigen3/Eigen/Dense>

#include "Native_interface.hh"
#include "APL_types.hh"
#include "Shape.hh"
#include "Value.hh"

#ifdef HAVE_CONFIG_H
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_URL
#undef PACKAGE_VERSION
#undef VERSION
#include "../config.h"
#endif

using namespace std;

class NativeFunction;

extern "C" void * get_function_mux(const char * function_name);

class Matrix
{
public:
  Matrix (int r, int c);
  ~Matrix ();
  complex<double> val (int r, int c);
  void val (int r, int c, complex<double> v);
  int  rows ();
  int  cols ();
  void show ();   
  
private:
  int krows;
  int kcols;
  vector<complex<double>> *vals;
};

Matrix::Matrix (int r, int c)
{
  krows = r;
  kcols = c;
  vals = new vector<complex<double>>(r * c);
}

Matrix::~Matrix ()
{
  delete vals;
}

complex<double>
Matrix::val (int r, int c)
{
  complex<double>rc (0.0, 0.0);
  rc = (*vals)[c + r * kcols];
  return rc;
}

void
Matrix::val (int r, int c, complex<double> v)
{
  (*vals)[c + r * kcols] = v;
}

int
Matrix::rows ()
{
  return krows;
}

int
Matrix::cols ()
{
  return kcols;
}

void
Matrix::show ()
{
  for (int r = 0; r < krows; r++) {
    for (int c = 0; c < kcols; c++) {
      fprintf (stderr, "%gj%g ",
	       this->val (r, c).real (),
	       this->val (r, c).imag ());
    }
    fprintf (stderr, "\n");
  }
}

enum {
  OP_UNKNOWN,
  OP_DETERMINANT,
  OP_CROSS_PRODUCT,
  OP_VECTOR_ANGLE,
  OP_EIGENVECTORS,
  OP_EIGENVALUES
};

static bool
close_fun(Cause cause, const NativeFunction * caller)
{
   return true;
}

// prevent compiler warning
bool (*close_fun_is_unused)(Cause, const NativeFunction *) = &close_fun;


Fun_signature
get_signature()
{
  return SIG_Z_A_F2_B;
}


static Token
eval_fill_B(Value_P B, const NativeFunction * caller)
{
  fprintf (stderr, "eval_fill_B\n");
  return Token(TOK_APL_VALUE1, Str0(LOC));
}

static Token
eval_fill_AB(Value_P A, Value_P B, const NativeFunction * caller)
{
  fprintf (stderr, "eval_fill_AB\n");
  return Token(TOK_APL_VALUE1, Str0(LOC));
}

static Token
eval_ident_Bx(Value_P B, sAxis x, const NativeFunction * caller)
{
  fprintf (stderr, "eval_ident_B\n");
  return Token(TOK_APL_VALUE1, Str0(LOC));
}

static Matrix *
genCofactor (Matrix *mtx, int r, int c)
{
  Matrix *cf = new Matrix (mtx->rows ()-1, mtx->cols ()-1);

  int rx = 0;
  for (int rr = 0; rr < mtx->rows (); rr++) {
    if (rr != r) {
      int cx = 0;
      for (int cc = 0; cc < mtx->cols (); cc++) {
	if (cc != c) {
	  cf->val (rx, cx++, mtx->val (rr, cc));
	}
      }
      rx++;
    }
  }

  return cf;
}

static Matrix
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

  for (int l = 0; l < mtx->cols (); l++) {
    for (int m = 0; m < mtx->rows (); m++) {
      res.val (l, m, (ces.eigenvectors ().col (l))(m));
    }
  }
  return res;
}

static vector<complex<double>>
getEigenvalues (Matrix *mtx)
{
  vector<complex<double>> res (mtx->cols ());
  Eigen::MatrixXcd mtxE (mtx->rows (), mtx->cols ());
  for (int j = 0; j < mtx->rows (); j++) {
    for (int k = 0; k < mtx->cols (); k++) {
      mtxE (j, k) = mtx->val (j, k);
    }
  }

  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces(mtxE);

  for (int l = 0; l < mtx->rows (); l++)
    res[l] = ces.eigenvalues ()(l);

  return res;
}

static complex<double>
getDet (Matrix *mtx)
{	
  complex<double> det (0.0, 0.0);
  if (mtx->rows () == 2)
    det = (mtx->val (0, 0) * mtx->val (1,1)) -
	   (mtx->val (0, 1) * mtx->val (1,0));
  else {
    for (int k = 0; k < mtx->rows (); k++) {
      double sign = pow (-1.0, (double)(k));
      Matrix *cf = genCofactor (mtx, 0, k);
      complex<double> id = getDet (cf);
      delete cf;
      det += sign * mtx->val (0, k) * id;
    }
  }

  return det;
}  

static vector<complex<double>>
getCross (Matrix *mtx)
{
      double data[] = { -1.0, 1.0, -1.0, 1.0,
                    -8.0, 4.0, -2.0, 1.0,
                    27.0, 9.0, 3.0, 1.0,
                    64.0, 16.0, 4.0, 1.0 };

  complex<double> det (0.0, 0.0);
  vector<complex<double>> rc(mtx->cols ());
  if (mtx->rows () == 2)
    det = (mtx->val (0, 0) * mtx->val (1,1)) -
	   (mtx->val (0, 1) * mtx->val (1,0));
  else {
    for (int k = 0; k < mtx->rows (); k++) {
      double sign = pow (-1.0, (double)(k));
      Matrix *cf = genCofactor (mtx, 0, k);
      complex<double> id = getDet (cf);
      delete cf;
      rc[k] = sign * id;
    }
  }

  return rc;
}  

// https://www.microapl.com/apl_help/ch_020_030_030.htm

static Token
eval_XB(Value_P X, Value_P B, const NativeFunction * caller)
{
  /***
      det n x n
      cp  n x n+1
   ***/

  Value_P rc = Str0(LOC);
  
  int op = OP_UNKNOWN;
  
  if (X->is_char_string ()) {
    const UCS_string  ustr = X->get_UCS_ravel();
    UTF8_string which (ustr);
    switch(*which.c_str ()) {
    case 'd':
    case 'D': op = OP_DETERMINANT; break;
    case 'c':
    case 'C': op = OP_CROSS_PRODUCT; break;
    }
    if (op == OP_UNKNOWN) {
      if (!strcasecmp (which.c_str (), "eigenvector"))
	op = OP_EIGENVECTORS;
      else if (!strcasecmp (which.c_str (), "eigenvalue"))
	op = OP_EIGENVALUES;
    }
    if (op == OP_UNKNOWN) {
      SYNTAX_ERROR;
    }
  }
  else if (X->is_numeric_scalar()) 
    op = X->get_sole_integer ();

  const CellType celltype = B->deep_cell_types();

  if (celltype & CT_NUMERIC) {
    const ShapeItem count   = B->element_count();
    const auto      rank    = B->get_rank();
    switch((int)rank) {
    case 0:		// scalar
      {
	switch(op) {
	case OP_CROSS_PRODUCT:
	  UERR << "Scalar argument.  No cross product posible." << endl;
	  RANK_ERROR;
	  break;
	case OP_EIGENVECTORS:
	case OP_EIGENVALUES:
	  UERR << "Scalar argument.  No eigen ops posible." << endl;
	  RANK_ERROR;
	  break;
	}
      }
      break;
    case 1:		// vector --only count == 1 vectors will work for det
      {
	switch(op) {
	case OP_DETERMINANT:
	  if (B->is_empty ()) {
	    LENGTH_ERROR;
	    UERR << "Null argument." << endl;
	  }
	  else {
	    if (count == 1) {
	      const Cell & Bx = B->get_cfirst ();
	      APL_Float xvr = Bx.get_real_value ();
	      if (Bx.is_complex_cell ()) {
		APL_Float xvi = Bx.get_imag_value ();
		rc = ComplexScalar(xvr, xvi, LOC);
	      }
	      else 
		rc = FloatScalar(xvr, LOC);
	    }
	    else {
	      UERR << "Vector argument.  No determinant posible." << endl;
	      LENGTH_ERROR;
	    }
	  }
	  break;
	case OP_CROSS_PRODUCT:
	  UERR << "Vector argument.  No cross product posible." << endl;
	  RANK_ERROR;
	  break;
	case OP_EIGENVECTORS:
	case OP_EIGENVALUES:
	  UERR << "Vector argument.  No eigen ops posible." << endl;
	  RANK_ERROR;
	  break;
	}
      }
      break;
    case 2:		// matrix
      {
	if (B->is_empty ()) {
	  LENGTH_ERROR;	
	  UERR << "Null argument." << endl;
	}
	
	ShapeItem rows = B->get_shape_item(0);
	ShapeItem cols = B->get_shape_item(1);

	switch(op) {
	case OP_DETERMINANT:
	  rc = B;
	  break;
	case OP_CROSS_PRODUCT:
	  rows++;
	  if (rows != cols) {
	    UERR <<
	"For cross product, the shape of the argument must be [n-1 n]" << endl;
	    RANK_ERROR;
	  }
	}
	if (op != OP_CROSS_PRODUCT && rows != cols) {
	  UERR << "Not a square matrix." << endl;
	  RANK_ERROR;
	}
	
	Matrix *mtx = new Matrix (rows, cols);

	{
	  int p = 0;

	  loop (r, rows) {
	    loop (c, cols) {
	      if (op == OP_CROSS_PRODUCT && r == 0)
		mtx->val (0, c, complex (1.0, 0.0));
	      else {
		const Cell & Bv = B->get_cravel (p++);
		APL_Float xvr = Bv.get_real_value ();
		APL_Float xvi = Bv.is_complex_cell ()
		  ? Bv.get_imag_value () : 0.0;
		mtx->val (r, c, complex (xvr, xvi));
	      }
	    }
	  }
	}

	switch(op) {
	case OP_EIGENVECTORS:
	  {
	    Matrix res = getEigenvectors (mtx);
	    Shape shape_Z;
	    shape_Z.add_shape_item(mtx->rows () * mtx->cols ());
	    rc = Value_P (shape_Z, LOC);
	    int p = 0;
	    for (int i = 0; i < mtx->rows (); i++) {
	      for (int j = 0; j < mtx->cols (); j++, p++) 
		(*rc).set_ravel_Complex (p,
					 res.val (i, j).real (),
					 res.val (i, j).imag ());
	    }
	    rc->check_value(LOC);
	    Shape shape_W;
	    shape_W.add_shape_item(mtx->rows ());
	    shape_W.add_shape_item(mtx->cols ());
	    (*rc).set_shape (shape_W);
	  }
	  break;
	case OP_EIGENVALUES:
	  {
	    vector<complex<double>> res = getEigenvalues (mtx);
	    Shape shape_Z;
	    shape_Z.add_shape_item(mtx->cols ());
	    rc = Value_P (shape_Z, LOC);
	    for (int i = 0; i < mtx->cols (); i++) 
	      (*rc).set_ravel_Complex (i, res[i].real (), res[i].imag ());
	    rc->check_value(LOC);
	  }
	  break;
	case OP_DETERMINANT:
	  {
	    complex<double>det = getDet (mtx);
	    rc = (det.imag () == 0.0) ?
	      FloatScalar(det.real (), LOC) :
	      ComplexScalar(det.real (), det.imag (), LOC);
	  }
	  break;
	case OP_CROSS_PRODUCT:
	  {
	    vector<complex<double>> cp = getCross (mtx);
	    Shape shape_Z;
	    shape_Z.add_shape_item(mtx->cols ());
	    rc = Value_P (shape_Z, LOC);
	    bool is_cpx = false;
	    for (int i = 0; i < mtx->cols (); i++) {
	      if (cp[i].imag () != 0.0) {
		is_cpx = true;
		break;
	      }
	    }
	    for (int i = 0; i < mtx->cols (); i++) {
	      if (is_cpx) 
		(*rc).set_ravel_Complex (i, cp[i].real (), cp[i].imag ());
	      else
		(*rc).set_ravel_Float (i, cp[i].real ());
	    }
	    rc->check_value(LOC);
	  }
	  break;
	}

	delete mtx;
      }
      break;
    case OP_EIGENVECTORS:
    case OP_EIGENVALUES:
      UERR << "Scalar argument.  No eigen ops posible." << endl;
      RANK_ERROR;
      break;
    default:		// can't deal with it
      UERR << "Invalid rank.." << endl;
      RANK_ERROR;
      break;
    }
  }
  else {				// not numeric
    UERR << "Non-numeric argument." << endl;
    DOMAIN_ERROR;
  }

  return Token(TOK_APL_VALUE1, rc);
}

static Token
eval_B(Value_P B, const NativeFunction * caller)
{
  Value_P X = IntScalar (OP_DETERMINANT, LOC);	// default to determinant
  return eval_XB (X, B, caller);
}

static
complex<double>
magnitude (vector<complex<double>> &v)
{
  complex<double> rc (0.0, 0.0);
  loop (i, v.size ()) rc += v[i] * v[i];
  rc = sqrt (rc);
  return rc;
}

static Token
eval_AXB(Value_P A, Value_P X, Value_P B,
	 const NativeFunction * caller)
{
  int op = OP_UNKNOWN;
  
  if (X->is_char_string ()) {
    const UCS_string  ustr = X->get_UCS_ravel();
    UTF8_string which (ustr);
    switch(*which.c_str ()) {
    case 'a':
    case 'A': op = OP_VECTOR_ANGLE; break;
    case 'c':
    case 'C': op = OP_CROSS_PRODUCT; break;
    }
    if (op == OP_UNKNOWN) {
      SYNTAX_ERROR;
    }
  }
  else if (X->is_numeric_scalar()) 
    op = X->get_sole_integer ();

  Value_P rc = Str0(LOC);
  
  const CellType A_celltype = A->deep_cell_types();
  const CellType B_celltype = B->deep_cell_types();

  if ((A_celltype & CT_NUMERIC) &&
      (A_celltype & CT_NUMERIC)) {
  }
  else {				// not numeric
    UERR << "Non-numeric argument." << endl;
    DOMAIN_ERROR;
  }
  
  const ShapeItem A_count   = A->element_count();
  const auto      A_rank    = A->get_rank();
  const ShapeItem B_count   = B->element_count();
  const auto      B_rank    = B->get_rank();

  switch(op) {
  case OP_CROSS_PRODUCT:
    if (A_rank == 1 && B_rank == 1 && A_count == B_count && A_count == 3) {
      Matrix *mtx = new Matrix (3, 3);
      loop (c, 3) {
	const Cell & Av = A->get_cravel (c);
	const Cell & Bv = B->get_cravel (c);
	APL_Float Avr = Av.get_real_value ();
	APL_Float Avi = Av.is_complex_cell () ? Av.get_imag_value () : 0.0;
	APL_Float Bvr = Bv.get_real_value ();
	APL_Float Bvi = Bv.is_complex_cell () ? Bv.get_imag_value () : 0.0;
	mtx->val (0, c, complex (1.0, 0.0));
	mtx->val (1, c, complex (Avr, Avi));
	mtx->val (2, c, complex (Bvr, Bvi));
      }
      {
	vector<complex<double>> cp = getCross (mtx);
	Shape shape_Z;
	shape_Z.add_shape_item(mtx->cols ());
	rc = Value_P (shape_Z, LOC);
	bool is_cpx = false;
	for (int i = 0; i < mtx->cols (); i++) {
	  if (cp[i].imag () != 0.0) {
	    is_cpx = true;
	    break;
	  }
	}
	for (int i = 0; i < mtx->cols (); i++) {
	  if (is_cpx) 
	    (*rc).set_ravel_Complex (i, cp[i].real (), cp[i].imag ());
	  else
	    (*rc).set_ravel_Float (i, cp[i].real ());
	}
	rc->check_value(LOC);
      }
    }
    else {
      UERR << "Invalid rank.." << endl;
      RANK_ERROR;
    }
    break;
  case OP_VECTOR_ANGLE:
    if (A_rank == 1 && B_rank == 1 && A_count == B_count) {
      vector<complex<double>> Av (A_count);
      vector<complex<double>> Bv (B_count);
      loop (c, A_count) {
	const Cell & Ac = A->get_cravel (c);
	const Cell & Bc = B->get_cravel (c);
	APL_Float Avr = Ac.get_real_value ();
	APL_Float Avi = Ac.is_complex_cell () ? Ac.get_imag_value () : 0.0;
	APL_Float Bvr = Bc.get_real_value ();
	APL_Float Bvi = Bc.is_complex_cell () ? Bc.get_imag_value () : 0.0;
	Av[c] = complex<double> (Avr, Avi);
	Bv[c] = complex<double> (Bvr, Bvi);
      }
      complex<double> Amag = magnitude (Av);
      complex<double> Bmag = magnitude (Bv);
      complex<double> mag = Amag * Bmag;
      if (mag != complex<double>(0.0, 0.0)) {
	complex<double> dp (0.0, 0.0);
	loop (i, Av.size ()) dp += Av[i] * Bv[i];
	complex<double> an = acos (dp/mag);
	rc = ComplexScalar((APL_Float)an.real (), an.imag (), LOC);
      }
      else {
	UERR << "Invalid vector(s)." << endl;
	DOMAIN_ERROR;
      }
    }
    break;
  default:
    UERR << "Not a dyadic operation.\n";
    DOMAIN_ERROR;
    break;
  }
  

  return Token(TOK_APL_VALUE1, rc);
}

static Token
eval_AB(Value_P A, Value_P B, const NativeFunction * caller)
{
  Value_P X = IntScalar (OP_DETERMINANT, LOC);	// default to determinant
  return eval_AXB (A, X, B, caller);
}



void *
get_function_mux(const char * function_name)
{
  // mandatory
  if (!strcmp(function_name, "get_signature"))
    return reinterpret_cast<void *>(&get_signature);
  if (!strcmp(function_name, "eval_fill_B"))
    return reinterpret_cast<void *>(&eval_fill_B);
  if (!strcmp(function_name, "eval_fill_AB"))
    return reinterpret_cast<void *>(&eval_fill_AB);
  if (!strcmp(function_name, "eval_ident_Bx"))
    return reinterpret_cast<void *>(&eval_ident_Bx);

  // ad hoc
  if (!strcmp(function_name, "eval_B"))
    return reinterpret_cast<void *>(&eval_B);
  if (!strcmp(function_name, "eval_AB"))
    return reinterpret_cast<void *>(&eval_AB);
  if (!strcmp(function_name, "eval_XB"))
    return reinterpret_cast<void *>(&eval_XB);
  if (!strcmp(function_name, "eval_AXB"))
    return reinterpret_cast<void *>(&eval_AXB);

  return 0;
}


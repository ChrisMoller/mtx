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

// https://www.microapl.com/apl_help/ch_020_030_030.htm
// https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions

#include "../mtx_config.h"

#include <stdio.h>

#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_URL
#undef PACKAGE_VERSION
#undef VERSION

#include<alloca.h>
#include<strings.h>
#include<cmath>
#include<complex>
#include <random>
#include<iostream>
#include<fstream>
#include<string>

#include "Native_interface.hh"
#include "APL_types.hh"
#include "Shape.hh"
#include "Value.hh"

#include "eigens.hh"

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

using namespace std;

class NativeFunction;

extern "C" void * get_function_mux(const char * function_name);

enum {
  OP_UNKNOWN,
  OP_DETERMINANT,
  OP_CROSS_PRODUCT,
  OP_VECTOR_ANGLE,
  OP_EIGENVECTORS,
  OP_EIGENVALUES,
  OP_IDENT,
  OP_ROTATION_MATRIX,
  OP_HOMOGENEOUS_MATRIX,
  OP_NORM,
  OP_GAUSSIAN,
  OP_PRINT
};

static bool
close_fun(Cause cause, const NativeFunction * caller)
{
   return true;
}

bool (*close_fun_is_unused)(Cause, const NativeFunction *) = &close_fun;

static complex<double>
genRand (double rsdev, double isdev)
{
  random_device rrd{};
  random_device ird{};
  mt19937 rgen{rrd()};
  mt19937 igen{ird()};
  normal_distribution rd{0.0, rsdev};
  normal_distribution id{0.0, isdev};

  complex<double> rc = complex (rd (rgen), id (igen));
  
  return rc;
}

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
eval_ident_Bx(Value_P B, sAxis x)
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

//https://misc.flogisoft.com/bash/tip_colors_and_formatting

static inline void
showHelpLine (const char first, const char * rest)
{
  char *str;
  if (first == 0)
    asprintf (&str, "\t\e[4m%s\e[0m", rest);
  else
    asprintf (&str, "\t\e[4m%c\e[0m%s", first, rest);
  cout << str << endl;
  free (str);
}

static Value_P
genRands (Value_P B)
{
  const ShapeItem count   = B->element_count();
  Shape shape_Z;
  shape_Z.add_shape_item(count);
  Value_P rc = Value_P (shape_Z, LOC);

  bool is_cpx = false;
  for (int i = 0; i < count; i++) {
    const Cell & Bv = B->get_cravel (i);
    if (Bv.is_complex_cell () &&
	Bv.get_imag_value () != 0.0) {
      is_cpx = true;
      break;
    }
  }
  
  for (int i = 0; i < count; i++) {
    const Cell & Bv = B->get_cravel (i);
    APL_Float xvr = Bv.get_real_value ();
    APL_Float xvi = Bv.is_complex_cell ()
      ? Bv.get_imag_value () : 0.0;
      complex<double> val = genRand (xvr, xvi);
      if (is_cpx) 
	(*rc).set_ravel_Complex (i, val.real (), val.imag ());
      else
	(*rc).set_ravel_Float (i, val.real ());
  }
  rc->check_value(LOC);
  
  const auto rank = B->get_rank();
  Shape shape_W;
  if (rank > 1) {
    for (sAxis r = 0; r < rank; r++) {
      ShapeItem s = B->get_shape_item (r);
      shape_W.add_shape_item (s);
    }
    (*rc).set_shape (shape_W);
  }
  
  return rc;
}

static Value_P
genRotation (int tp, Value_P A, Value_P B)
{
  Value_P rc;
  Shape shape_Z;
  int sdim = (tp == 9) ? 3 : 4;
  shape_Z.add_shape_item(tp);
  rc = Value_P (shape_Z, LOC);
  vector<complex<double>> angs (3);
  for (int i = 0; i < 3; i++) {
    const Cell & Bv = B->get_cravel (i);
    angs[i] =
      complex<double>(Bv.get_real_value (),
		      Bv.is_complex_cell () ?
		      Bv.get_imag_value () : 0.0);
  }
  complex<double>cosa = cos (angs[0]);
  complex<double>sina = sin (angs[0]);
  complex<double>cosb = cos (angs[1]);
  complex<double>sinb = sin (angs[1]);
  complex<double>cosg = cos (angs[2]);
  complex<double>sing = sin (angs[2]);

  complex<double>t00 = cosa * cosb;
  complex<double>t01 = cosa * sinb * sing - sina * cosg;;
  complex<double>t02 = cosa * sinb * cosg + sina * sing;;
	      
  complex<double>t10 = sina * cosb;
  complex<double>t11 = sina * sinb * sing + cosa * cosg;;
  complex<double>t12 = sina * sinb * cosg - cosa * sing;;
	      
  complex<double>t20 = -sinb;
  complex<double>t21 =  cosb * sing;
  complex<double>t22 =  cosb * cosg;

  if (tp == 9) {
    (*rc).set_ravel_Complex (0,  t00.real (),  t00.imag ());
    (*rc).set_ravel_Complex (1,  t01.real (),  t01.imag ());
    (*rc).set_ravel_Complex (2,  t02.real (),  t02.imag ());
    (*rc).set_ravel_Complex (3,  t10.real (),  t10.imag ());
    (*rc).set_ravel_Complex (4,  t11.real (),  t11.imag ());
    (*rc).set_ravel_Complex (5,  t12.real (),  t12.imag ());
    (*rc).set_ravel_Complex (6,  t20.real (),  t20.imag ());
    (*rc).set_ravel_Complex (7,  t21.real (),  t21.imag ());
    (*rc).set_ravel_Complex (8,  t22.real (),  t22.imag ());
  }
  else {
    vector<complex<double>> trans (3);
    for (int i = 0; i < 3; i++) {
      const Cell & Av = A->get_cravel (i);
      trans[i] =
	complex<double>(Av.get_real_value (),
			Av.is_complex_cell () ?
			Av.get_imag_value () : 0.0);
    }
    (*rc).set_ravel_Complex (0,   t00.real (),  t00.imag ());
    (*rc).set_ravel_Complex (1,   t01.real (),  t01.imag ());
    (*rc).set_ravel_Complex (2,   t02.real (),  t02.imag ());
    (*rc).set_ravel_Complex (3,   0.0, 0.0);
    
    (*rc).set_ravel_Complex (4,   t10.real (),  t10.imag ());
    (*rc).set_ravel_Complex (5,   t11.real (),  t11.imag ());
    (*rc).set_ravel_Complex (6,   t12.real (),  t12.imag ());
    (*rc).set_ravel_Complex (7,   0.0, 0.0);
    
    (*rc).set_ravel_Complex (8,   t20.real (),  t20.imag ());
    (*rc).set_ravel_Complex (9,   t21.real (),  t21.imag ());
    (*rc).set_ravel_Complex (10,  t22.real (),  t22.imag ());
    (*rc).set_ravel_Complex (11,  0.0, 0.0);
    
    (*rc).set_ravel_Complex (12,  trans[0].real (), trans[0].imag ());
    (*rc).set_ravel_Complex (13,  trans[1].real (), trans[1].imag ());
    (*rc).set_ravel_Complex (14,  trans[2].real (), trans[2].imag ());
    (*rc).set_ravel_Complex (15,  1.0, 0.0);
  }
  
  rc->check_value(LOC);
  Shape shape_W;
  shape_W.add_shape_item (sdim);
  shape_W.add_shape_item (sdim);
  (*rc).set_shape (shape_W);
  return rc;
}

static Token
eval_XB(Value_P X, Value_P B, const NativeFunction * caller)
{
  Value_P rc = Str0(LOC);
  
  if (B->is_char_string ()) {
    cout << "\n\tgeneral form: mtx['\e[3mindex\e[0m']\n\n";
    cout << "\tvalid '\e[3mindex\e[0m' values: (underscored minimal)\n\n";
    showHelpLine ('d', "eterminant");
    showHelpLine ('c', "ross_product");
    showHelpLine ('i', "ident");
    showHelpLine ('r', "otation");
    showHelpLine ('n', "orm");
    showHelpLine (0, "eigenvector");
    showHelpLine (0, "eigenvalue");
    cout << "\n\tdefaults to determinant if the index is omitted.\n\n";
    return Token(TOK_APL_VALUE1, rc);
  }
  
  int op = OP_UNKNOWN;
  
  if (X->is_char_string ()) {
    const UCS_string  ustr = X->get_UCS_ravel();
    UTF8_string which (ustr);
    switch(*which.c_str ()) {
    case 'd':						// working
    case 'D': op = OP_DETERMINANT; break;
    case 'c':						// working
    case 'C': op = OP_CROSS_PRODUCT; break;
    case 'i':						// working
    case 'I': op = OP_IDENT; break;
    case 'r':						// working
    case 'R': op = OP_ROTATION_MATRIX; break;
    case 'n':						// working
    case 'N': op = OP_NORM; break;
    case 'g':						// working
    case 'G': op = OP_GAUSSIAN; break;
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
	case OP_GAUSSIAN:
	  {
	    const Cell & Bv = B->get_cravel (0);
	    APL_Float xvr = Bv.get_real_value ();
	    APL_Float xvi = Bv.is_complex_cell ()
	      ? Bv.get_imag_value () : 0.0;
	    complex<double> val = genRand (xvr, xvi);
	    rc = (xvi == 0.0) ?
	      FloatScalar (val.real (), LOC) :
	      ComplexScalar (val.real (), val.imag (), LOC);
	  }	
	  break;
	case OP_NORM:
	  UERR << "Scalar argument.  Normalisation posible." << endl;
	  RANK_ERROR;
	  break;
	case OP_ROTATION_MATRIX:
	  {
	    Shape shape_Z;
	    shape_Z.add_shape_item(4);
	    rc = Value_P (shape_Z, LOC);
	    const Cell & Bv = B->get_cravel (0);
	    APL_Float xvr = Bv.get_real_value ();
	    APL_Float xvi = Bv.is_complex_cell ()
	      ? Bv.get_imag_value () : 0.0;
	    complex<double>theta (xvr, xvi);
	    complex<double>cosx = cos (theta);
	    complex<double>sinx = sin (theta);
	    (*rc).set_ravel_Complex (0,  cosx.real (),  cosx.imag ());
	    (*rc).set_ravel_Complex (1, -sinx.real (), -sinx.imag ());
	    (*rc).set_ravel_Complex (2,  sinx.real (),  sinx.imag ());
	    (*rc).set_ravel_Complex (3,  cosx.real (),  cosx.imag ());
	    rc->check_value(LOC);
	    Shape shape_W;
	    shape_W.add_shape_item (2);
	    shape_W.add_shape_item (2);
	    (*rc).set_shape (shape_W);
	  }
	  break;
	case OP_IDENT:
	  {
	    int dim = B->get_sole_integer();
	    Shape shape_Z;
	    shape_Z.add_shape_item(dim * dim);
	    rc = Value_P (shape_Z, LOC);

	    int p = 0;
	    for (int i = 0; i < dim; i++) {
	      for (int j = 0; j < dim; j++, p++) 
		(*rc).set_ravel_Complex (p,
					 ((i == j) ? 1.0 : 0.0),
					 0.0);
	    }
	    rc->check_value(LOC);
	    Shape shape_W;
	    shape_W.add_shape_item (dim);
	    shape_W.add_shape_item (dim);
	    (*rc).set_shape (shape_W);
	  }
	  break;
	default:
	  DOMAIN_ERROR;
	  break;
	}
      }
      break;
    case 1:		// vector --only count == 1 vectors will work for det
      {
	switch(op) {
	case OP_GAUSSIAN:
	  rc = genRands (B);
	  break;
	case OP_NORM:
	  {
	    complex<double> sum (0.0, 0.0);;
	    for (int i = 0; i < count; i++) {
	      const Cell & Bv = B->get_cravel (i);
	      complex<double> val =
		complex<double>(Bv.get_real_value (),
				Bv.is_complex_cell () ?
				Bv.get_imag_value () : 0.0);
	      sum += val * val;
	    }
	    sum = sqrt (sum);

	    Shape shape_Z;
	    shape_Z.add_shape_item(count);
	    rc = Value_P (shape_Z, LOC);

	    for (int i = 0; i < count; i++) {
	      const Cell & Bv = B->get_cravel (i);
	      complex<double> val =
		complex<double>(Bv.get_real_value (),
				Bv.is_complex_cell () ?
				Bv.get_imag_value () : 0.0);
	      val /= sum;
	      (*rc).set_ravel_Complex (i, val.real (), val.imag ());
	    }
	    rc->check_value(LOC);
	  }
	  break;
	case OP_ROTATION_MATRIX:
	  {
	    if (count == 3)
	      rc = genRotation (9, nullptr, B);
	    else
	      RANK_ERROR;
	  }
	  break;
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
	default:
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

	if (op != OP_CROSS_PRODUCT &&
	    op != OP_GAUSSIAN &&
	    op != OP_NORM &&
	    rows != cols) {
	  UERR << "Not a square matrix." << endl;
	  RANK_ERROR;
	}

	switch(op) {
	case OP_NORM:
	  {
	    complex<double> sum (0.0, 0.0);;
	    for (int i = 0; i < count; i++) {
	      const Cell & Bv = B->get_cravel (i);
	      complex<double> val =
		complex<double>(Bv.get_real_value (),
				Bv.is_complex_cell () ?
				Bv.get_imag_value () : 0.0);
	      sum += val * val;
	    }
	    sum = sqrt (sum);

	    Shape shape_Z;
	    shape_Z.add_shape_item(rows);
	    shape_Z.add_shape_item(cols);
	    rc = Value_P (shape_Z, LOC);

	    for (int i = 0; i < count; i++) {
	      const Cell & Bv = B->get_cravel (i);
	      complex<double> val =
		complex<double>(Bv.get_real_value (),
				Bv.is_complex_cell () ?
				Bv.get_imag_value () : 0.0);
	      val /= sum;
	      (*rc).set_ravel_Complex (i, val.real (), val.imag ());
	    }
	    rc->check_value(LOC);
	  }
	  break;
	case OP_GAUSSIAN:
	  rc = genRands (B);
	  break;
	case OP_CROSS_PRODUCT:
	  rows++;
	  if (rows != cols) {
	    UERR <<
	"For cross product, the shape of the argument must be [n-1 n]" << endl;
	    RANK_ERROR;
	  }
	}
	
	Matrix *mtx = new Matrix (rows, cols);

	if (op != OP_GAUSSIAN) {
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
	case OP_GAUSSIAN:	// do nothing
	  break;
	case OP_IDENT:
	  RANK_ERROR;
	  break;
	case OP_EIGENVECTORS:
	  {
#ifdef HAVE_EIGEN
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
#else
	  UERR << "Eigen operations not supported" << endl;
	  DOMAIN_ERROR;
#endif
	  }
	  break;
	case OP_EIGENVALUES:
	  {
#ifdef HAVE_EIGEN
	    vector<complex<double>> res = getEigenvalues (mtx);
	    Shape shape_Z;
	    shape_Z.add_shape_item(mtx->cols ());
	    rc = Value_P (shape_Z, LOC);
	    for (int i = 0; i < mtx->cols (); i++) 
	      (*rc).set_ravel_Complex (i, res[i].real (), res[i].imag ());
	    rc->check_value(LOC);
#else
	  UERR << "Eigen operations not supported" << endl;
	  DOMAIN_ERROR;
#endif
	  }
	  break;
	case OP_DETERMINANT:
	  {
	    complex<double>det = getDet (mtx);
	    rc = (det.imag () == 0.0) ?
	      FloatScalar(det.real (), LOC) :
	      ComplexScalar(det.real (), det.imag (), LOC);
	    rc->check_value(LOC);
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
    default:		// can't deal with it
      if (op == OP_GAUSSIAN)
	rc = genRands (B);
      else
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
    case 'a':						// working
    case 'A': op = OP_VECTOR_ANGLE; break;
    case 'c':						// working
    case 'C': op = OP_CROSS_PRODUCT; break;
    case 'p':						// working
    case 'P': op = OP_PRINT; break;
    case 'h':						// working
    case 'H': op = OP_HOMOGENEOUS_MATRIX; break;
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

  if (op == OP_PRINT) {
    if ((A_celltype & CT_NUMERIC) &&
	B->is_char_string ()) {
      const ShapeItem A_count   = A->element_count();
      const auto      A_rank    = A->get_rank();
      const UCS_string  ustr = B->get_UCS_ravel();
      UTF8_string fn (ustr);
      FILE *ofile = fopen (fn.c_str (), "w");
      if (!ofile) {
	UERR << "Open failure on " << ustr << endl;
	DOMAIN_ERROR;
      }
      if (A_rank <= 1) {
	for (int i = 0; i < A_count; i++) {
	  const Cell & Av = A->get_cravel (i);
	  APL_Float Avr = Av.get_real_value ();
	  if (Av.is_complex_cell ()) {
	    APL_Float Avi = Av.get_imag_value ();
	    fprintf (ofile, "%gj%g ", Avr, Avi);
	  }
	  else
	    fprintf (ofile, "%g ", Avr);
	}
	fprintf (ofile, "\n");
      }
      else {
	int end_line = (int)(A->get_shape_item (A_rank-1));
	int end_grid = end_line * (int)(A->get_shape_item (A_rank-2));
#define STR_LEN 256
	char str[STR_LEN];
	bool is_cpx = false;
	int max_len = -1;
	for (int i = 0; i < A_count; i++) {
	  int len;
	  const Cell & Av = A->get_cravel (i);
	  APL_Float Avr = Av.get_real_value ();
	  if (Av.is_complex_cell ()) {
	    APL_Float Avi = Av.get_imag_value ();
	    len = snprintf (str, STR_LEN, "%gj%g", Avr, Avi);
	    if (Avi != 0.0) is_cpx = true;
	  }
	  else
	    len = snprintf (str, STR_LEN, "%g", Avr);
	  if (max_len < len) max_len = len;
	}
	int *rho = (int *)alloca (A_rank * sizeof(int));
	bzero (rho, A_rank * sizeof(int));
	for (int i = 0; i < A_count; i++) {
	  if (A_rank > 2) {
	    if (0 == i%end_grid) {
	      fprintf (ofile, "\n[");
	      for (int j = 0; j < A_rank - 2; j++)
		fprintf (ofile, "%d ", rho[j]);
	      fprintf (ofile, "* *]:\n");
	      bool carry = 1;
	      for (int j = A_rank - 3; j >= 0; j--) {
		rho[j] += carry;
		if (rho[j] >= A->get_shape_item (j)) {
		  rho[j] = 0;
		  carry = 1;
		}
		else
		  carry = 0;
	      }
	    }
	  }
	  const Cell & Av = A->get_cravel (i);
	  APL_Float Avr = Av.get_real_value ();
	  char str[STR_LEN];
	  if (is_cpx && Av.is_complex_cell ()) {
	    APL_Float Avi = Av.get_imag_value ();
	    snprintf (str, STR_LEN, "%gj%g", Avr, Avi);
	  }
	  else
	    snprintf (str, STR_LEN, "%g", Avr);
	  fprintf (ofile, "%*s ", max_len, str);
	  if (0 == (i+1)%end_line) fprintf (ofile, "\n");
	}
      }

      fclose (ofile);
      COUT << "File " << fn << " printed.\n";
      return Token(TOK_APL_VALUE1, rc);
    }
    else {
      UERR << "Incompatible arguments.\n";
      DOMAIN_ERROR;
    }
  }

  if ((A_celltype & CT_NUMERIC) &&
      (B_celltype & CT_NUMERIC)) {
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
  case OP_HOMOGENEOUS_MATRIX:
    {
      if (A_count == 3 && B_count == 3) 
	rc = genRotation (16, A, B);
      else
	RANK_ERROR;
    }
    break;
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


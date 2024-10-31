
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

#include "Matrix.hh"


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

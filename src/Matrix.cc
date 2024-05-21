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

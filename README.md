# MTX

mtx is collection of matrix operations not natively supported by gnu APL.
including support for:
<ul>
<li>Matrix determinants</li>
<li>Matrix eigenvalues and eigenvectors</li>
<li>Identity matrices</li>
<li>Vector cross products</li>
<li>Vector interior angles</li>
<li>Vector or scalar rotation matrices</li>
<li>Homogeneous matrices</li>
<li>Vector/matrix normalisation</li>
<li>Gaussian complex random values</li>
</ul>

## Installation

See the INSTALL file on how to build mtx, but once it's built and installed
it must be fixed in the workspace:

&nbsp;&nbsp;&nbsp;&nbsp;     'libmtx.so' ⎕fx 'mtx'

using any function name that pleases you instead of 'mtx'.

The general form of the use of mtx is either monadic:

&nbsp;&nbsp;&nbsp;&nbsp;mtx[<em>string</em>] y  

or dyadic:

&nbsp;&nbsp;&nbsp;&nbsp;x mtx[<em>string</em>] y  

where <em>string</em> is one of:
<ul>
<li><em>d</em> -- determinant</li>
<li><em>eigenvalue</em></li>
<li><em>eigenvector</em></li>
<li><em>i</em> -- identity</li>
<li><em>c</em> -- cross</li>
<li><em>a</em> -- angle</li>
<li><em>r</em> -- rotation</li>
<li><em>g</em> -- gaussian</li>
<li><em>p</em> -- print</li>
<li><em>h</em> -- homogeneous</li>
<li><em>n</em> -- normalise</li>
</ul>

(The single-character strings can actually be any string starting with that
letter--spell out, if you like, <em>d</em> as <em>determinant</em>.  All
strings are case-insensitive.)

It's handy to create named lambdas for the mtx opertions.  The ones I use are:
<ul>
<li>det    ← {mtx['d'] ⍵}</li>
<li>eval   ← {mtx['eigenvalue'] ⍵}</li>
<li>evec   ← {mtx['eigenvector'] ⍵}</li>
<li>ident  ← {mtx['i'] ⍵}</li>
<li>cross  ← {⍺ mtx['c'] ⍵}</li>
<li>crossm ← {mtx['c'] ⍵}</li>
<li>angle  ← {⍺ mtx['a'] ⍵}</li>
<li>rotate ← {⍺ mtx['r'] ⍵}</li>
<li>grand  ← {mtx['g'] ⍵}</li>
<li>print  ← {⍺ mtx['p'] ⍵}</li>
<li>homgeneous  ← {⍺ mtx['h'] ⍵}</li>
<li>norm  ← {mtx['n'] ⍵}</li>
</ul>

## Details

### Monadic

#### Determinant

The argument for determinants must be a real or complex square matrix and the
function returns a real result if possible or a complex result if necessary.
E.g.,

>det 3 3⍴5 3 1 9 7 6 2 8 4

¯114

#### Eigenvalues, eigenvectors

The argument for these operations must be a real or complex square matrix.  The
eigenvalue function returns a vector of the values.  The eigenvector
function returns a matrix of the same shape of the argument where each row is
the eigenvector corresponding to that index of the eigenvalue vector.  I.e.,
eigenvalue[i] corresponds to eigenvector[i;].  E.g.

>t←4 4 ⍴ ¯1 1 ¯1 1 ¯8 4 ¯2 1 27 9 3 1 64 16 4 1

>evec t

<pre>
 0.0999           0.111         ¯0.293          ¯0.945
 0.0431J0.00969  ¯0.0709J0.139   0.517J¯0.016    0.84J0.0414
 0.0431J¯0.00969 ¯0.0709J¯0.139  0.517J0.016     0.84J¯0.0414
¯0.145            0.357          0.919           0.0812

</pre>

>eval t

¯2.34 3.23 15.1

(Eigenvalue ¯2.34 will correspond to eigenvector [¯0.158 0.637  ¯0.755]. and
so on.)


#### Identity

The argument for this operation must be a scalar integer and it returns a
complex square matrix of that dimension with the value 1.0j0.0 on the diagonal
and 0.0j0.0 elsewhere.

>ident 5

<pre>
1 0 0 0 0
0 1 0 0 0
0 0 1 0 0
0 0 0 1 0
0 0 0 0 1
</pre>


#### Rotate

If the argument is a scalar, the function returns a 2 × 2 2D rotation
transform matrix.  If the argument is a ⍴3 vector, the function returns
a 3 × 3  rotation transform matrix.

>rotate d2r 45

<pre>
0.707 ¯0.707
0.707  0.707
</pre>

>rotate d2r 30 45 60

<pre>
 0.612 0.28   0.739
 0.354 0.739 ¯0.573
¯0.707 0.612  0.354
</pre>


#### Normalise

For vector or matrix arguments, returns a value of the same shape such that
+/,(norm ⍵)*2 = 1.  (Not valid for scalar arguments.)

>norm 2 4 6

0.267 0.535 0.802

>+/,(norm 2 4 6)*2

1

>norm 2 3⍴⍳6

<pre>
0.105 0.21  0.314
0.419 0.524 0.629
</pre>

>+/,(norm 2 3⍴⍳6)*2

1

#### Gaussian

Gaussian returns normally distributed (i.e., gaussian) complex random numbers
of mean 0.0j0.0 and a standard deviation of the right argument.  I.e., repeated 
instances of  **mtx['g'] 2.0** would return normally distributed random
numbers.  If the argument is complex, the result will also be complex.  The
result will be of the same shape of the argument, e.g., **mtx['g'] 2 3 4**
will retirn a ⍴ 3 vector of randoms with standard veviations of 2.0, 3.0, and
4.0.  All results are of mean 0.0j0.0 but can be adjusted by addition,
subtraction, or whatever.

E.g., given 500 Gaussian-distributes

>t←grand (500 2⍴5.0 2.0)

>t print 's500.data'

rendered by gnuplot yields:

![A scatterplot of 500 normally distributed points in 2-space](s500.jpg
 "Normal distribution")

For 5000 points:

>grand (5000 2⍴5.0 2.0)) print 's5000.data'

![A scatterplot of 5000 normally distributed points in 2-space](s5000.jpg
 "Normal distribution")

(See Pretty-print (**mtx['p]**) below.)


### Dyadic

#### Angle

Both arguments must be rank 1 real or complex vectors and must be of the same
length.  The function returns a complex scalar value representing the angle
between the vectors.  (Frankly, I have no idea what it means if the result has
a non-zero imaginary component...)

>r2d (1 0) angle (1 1)

45.0

>r2d (1 0 0) angle 1 1 1

54.7


#### Print

Pretty-prints matrices to a file, converting all the APL ¯ "negative" symbols
to standard minus signs:

>*matrix* print '*filename*'

Actually, the "matrix" argument may be of any shape.

#### Homogeneous

Returns a homogeneous matrix of the given translation and rotation.

>1 2 3 homogeneous d2r 30 45 60

<pre>
 0.612 0.28   0.739 0
 0.354 0.739 ¯0.573 0
¯0.707 0.612  0.354 0
 1     2      3     1
</pre>

### Ambivalent

#### Cross product

In dynadic form, both arguments must be rank 1 real or complex vectors and
must be of length 3.  In monadic form, the argument must be a real or complex
matrix (rank 2) of shape <em>n</em>-1 × <em>n</em> and composed of <em>n</em>-1
real or complex vectors of length  <em>n</em>.

The function returns the real or complex cross product of those
vectors, i.e., a vector of the same length orthogonal to the arguments.
For reasons I'm not enough of a mathematician to understand, cross products
are only valid in 3-space and 7-space, so the only valid arguments are of
shapes [2 3] or [6 7].  mtx, however, doesn't check this and will happily
give you a result in any dimensionality and leave it your imagination what
it may mean.

>(1 2 3) cross 4 5 6

¯3 6 ¯3

>crossm (2 3⍴⍳6)

¯3 6 ¯3


I may add more functionality in later releases.  I'm open to suggestions.

Chris Moller  
henrik@henrikmoller.me

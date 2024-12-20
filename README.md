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
<li>Gaussian complex random values</li>
<li>Vector/matrix normalisation</li>
<li>Homogeneous matrices</li>
<li>Covariance</li>
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
<li><em>n</em> -- normalise</li>
<li><em>h</em> -- homogeneous</li>
<li><em>C</em> -- covariance</li>
<li><em>p</em> -- print</li>
</ul>

(The single-character strings can actually be any string starting with that
letter--spell out, if you like, <em>d</em> as <em>determinant</em>.)

It's handy to create named lambdas for the mtx opertions.  The ones I use are:
<ul>
<li>det    ← {mtx['d'] ⍵}</li>
<li>eval   ← {mtx['eigenvalue'] ⍵}</li>
<li>evec   ← {mtx['eigenvector'] ⍵}</li>
<li>ident  ← {mtx['i'] ⍵}</li>
<li>crossd ← {⍺ mtx['c'] ⍵}</li>
<li>crossm ← {mtx['c'] ⍵}</li>
<li>angle  ← {⍺ mtx['a'] ⍵}</li>
<li>rotate ← {⍺ mtx['r'] ⍵}</li>
<li>grand  ← {mtx['g'] ⍵}</li>
<li>norm   ← {mtx['n'] ⍵}</li>
<li>homgeneous  ← {⍺ mtx['h'] ⍵}</li>
<li>covd   ← {⍺ mtx['C'] ⍵}</li>
<li>covm   ← {mtx['C'] ⍵}</li>
<li>print  ← {⍺ mtx['p'] ⍵}</li>
</ul>

cov and cross are both ambivalent, but there's no apparent way to create
ambivalent lambdas, so I also use:

<pre>
z←l cross r
→(0≠⎕nc 'l')/Dyadic
z←mtx['c'] r
→0
Dyadic:
z←l mtx['c'] r
</pre>

<pre>
z←l cov r
→(0≠⎕nc 'l')/Dyadic
z←mtx['C'] r
→0
Dyadic:
z←l mtx['C'] r
</pre>

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

¯6.41 5.55J3.09 5.55J¯3.09 2.32

(Eigenvalue ¯6.41 will correspond to eigenvector [0.0999 0.111 ¯0.293 ¯0.945]
and so on.)


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
0     0.135 0.27
0.405 0.539 0.674
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

>r2d 1 0 angle 1 1

45.0

>r2d 1 0 0 angle 1 1 1

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

>1 2 3 crossd 4 5 6

¯3 6 ¯3

>crossm (2 3⍴⍳6)

¯3 6 ¯3

#### Covariance

Monadic covariance interprets each row of a matrix argument as samples of a
given variable, computes the covariance of each sample set with respect to
every other sample set, and presents the covariances as a matrix.  For example a
matrix:

<pre>
 6  3  5 15
 8  9 14 13
16 12  4  1
11 10  7  2
</pre>

would represent samples of four variables, [ 6  3  5 15 ] being sample of one
variable, [ 8  9 14 13] samples of another variable, and so on.  Labelling the
matrix M and these variables A, B, C, and D respectively,

>covm M

yields a matrix:

<pre>
cov(A, A)  cov(A, B)  cov(A, C)  cov(A, D)
cov(B, A)  cov(B, B)  cov(C, C)  cov(D, D)
cov(C, A)  cov(C, B)  cov(C, C)  cov(D, D)
cov(D, A)  cov(D, B)  cov(C, C)  cov(D, D)
</pre>

i.e., in this example:

<pre>
 28.3    7    ¯24.1 ¯18.8
  7      8.67 ¯19.3  ¯9.33
¯24.1  ¯19.3   48.3  26.2
¯18.8   ¯9.33  26.2  16.3
</pre>

The dyadic form of covariance returns the covariance of the left and right
arguments.  E.g.:

>2 1 3 covd 3 1 2

0.5

### Principal Component Analysis

The principal axis of sampled data is the eigensystem of the covariance of the
sample vectors:

<pre>
      m←4 4⍴16?16
      m
 9  7  3  2
 8 16  1 15
13 10 14  5
12  6  4 11

      evec covm m
¯0.00278 ¯0.862   0.445 ¯0.245
¯0.771   ¯0.0651 ¯0.401 ¯0.491
¯0.308   ¯0.395  ¯0.325  0.802
 0.558   ¯0.312  ¯0.732 ¯0.236

      eval covm m
63.1 16.6 11.1 0
</pre>

A typical graphical example using a two-variable system of 100 randomised
samples each results in eigenvalues of [1142.1j0 13.9427j0] and eigenvectors
of:

<pre>
   0.870314    0.492497 
  -0.492497    0.870314
</pre>

as represented in:

![PCA example data](pca.jpg "PCA example data")

Note that

>r2d (0.870314    0.492497) angle (¯0.492497    0.870314)

90

I.e., the eigenvectors are orthogonal, as expected.

Note also that the eigenvectors are origin-based.  In the image above, they've
been arbitrarily translated to where you seem them, at the centroid of the
sample data.

The APL for all this:

---
<pre>
m←c pca a;x;y;c;d;el;ec;xb;yb;xa;ya;co;av;s;⎕io;m0;m1;h
⎕io←0
e←tand a
x←4×grand c⍴1
y←4×grand c⍴1
xb←x+⍳c
yb←y+e×⍳c
xa←(+/xb)÷⍴xb
ya←(+/yb)÷⍴yb
av←xa,ya
d←2 c⍴xb,yb
⊣(⍉2 c⍴d) print 'pca.data'
co←covm d
el←eval co
ec←20×evec co

⊣'eigenvectors' print 'pcaeigensystem.txt'
⊣(ec÷20) print '>pcaeigensystem.txt'
⊣'eigenvalues' print '>pcaeigensystem.txt'
⊣el print '>pcaeigensystem.txt'

xa← 100×(1↑el)÷+/el
xb← 100×(¯1↑el)÷+/el

m0←atand ÷/⌽1 2↑ec
m1←atand ÷/⌽¯1 2↑ec

h←'w' ⎕fio[3] 'pcalabels.gp'   ⍝ open

s←⍕"set label \"eigenvalue λ0 = ", (1↑el), " (", xa, "%)\" at graph .1,.8\n"
⊣'%s' s ⎕fio[22] h  ⍝ fwrite

s←⍕"set label \"eigenvalue λ1 = ", (¯1↑el), " (", xb, "%)\" at graph .1,.75\n"
⊣'%s' s ⎕fio[22] h  ⍝ fwrite

s←⍕"set label \"slope λ0 = ", m0, " degress (nominally ", a, ")\" at graph .1,.7\n"
⊣'%s' s ⎕fio[22] h  ⍝ fwrite

s←⍕"set label \"slope λ1 = ", m1, " degrees (nominally ", (a-90), ")\" at graph .1,.65\n"
⊣'%s' s ⎕fio[22] h  ⍝ fwrite

⊣⎕fio[4] h


⊣(2 2⍴(av-(,1 2↑ec)÷2),av+(,1 2↑ec)÷2) print 'pcae1.data'
⊣(2 2⍴(av-(,¯1 2↑ec)÷2),av+(,¯1 2↑ec)÷2) print 'pcae2.data'

</pre>
---

The jpg was created by gnuplot using:

---
<pre>
set terminal jpeg enhanced size 512,360 font NimbusRoman 10

set output "pca.jpg"

set title "Principal component analysis of 2D biased random samples"

load "pcalabels.gp"

plot  [0:120][0:80] "pca.data" using 1:2 with points pt 7 ps .75 t 'Samples', \
   "pcae1.data" using 1:2 with lines lc "red" t "λ0 axis", \
   "pcae2.data" using 1:2 with lines lc "green" t "λ1 axis"
</pre>
---

The gnuplot data was created by:

>100 pca 30

<pre>



</pre>


### The End

I may add more functionality in later releases.  I'm open to suggestions.

Chris Moller  
henrik@henrikmoller.me

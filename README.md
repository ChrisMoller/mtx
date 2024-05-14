# MTX

mtx is collection of matrix operations not natively supported by gnu APL.
For this frst release it includes support for finding the determinant of
a matrix and the cross product of vectors.  

To use mtx, it must be fixed in the workspace:

&nbsp;&nbsp;&nbsp;&nbsp;     'libmtx.so' ⎕fx 'mtx'

using any function name that pleases you instead of 'mtx'.

The general form of the use of mtx is either monadic:

&nbsp;&nbsp;&nbsp;&nbsp;mtx y  
&nbsp;&nbsp;&nbsp;&nbsp;mtx['c'] y  
&nbsp;&nbsp;&nbsp;&nbsp;mtx['d'] y  

or dyadic:

&nbsp;&nbsp;&nbsp;&nbsp;x mtx y

In the monadic form, an "axis" of 'c' (or any string that starts with 'c' or
'C') will do the cross product of the arguement.  'd" or 'D' get the
determinant.  If neither is specified, the default is determinant.  The dyadic
form always yields the cross product.   

(It's
handy to create a couple of named lambdas for these:  

&nbsp;&nbsp;&nbsp;&nbsp;cross ← {⍺ mtx['c'] ⍵}   

and

&nbsp;&nbsp;&nbsp;&nbsp;det ← {mtx['d'] ⍵}   
)

The argument for determinants must be a real or complex square matrix. For
cross products in the dyadic form, both arguments must be of rank 1 and of
length 3.  In monadic form, the argument must have a shape of [n-1 n] though
for reasons I'm not enough of a mathematician to understand, cross products
are only valid in 3-space and 7-space, so the only valid arguments are of
shapes [2 3] or [6 7].  mtx, however, doesn't check this and will happily
give you a result in any dimensionality and leave it your imagination what
it may mean.

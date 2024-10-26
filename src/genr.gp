
set terminal jpeg enhanced size 512,360 font NimbusRoman 10

set output "genr.jpg"

set title "Principal component analysis of 2D biased random samples"

load "labels.gp"

plot  [0:120][0:80] "genr.data" using 1:2 with points pt 7 ps .75 t 'Samples', \
   "genre1.data" using 1:2 with lines lc "red" t "λ0 axis", \
   "genre2.data" using 1:2 with lines lc "green" t "λ1 axis"

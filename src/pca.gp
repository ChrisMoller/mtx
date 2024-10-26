
set terminal jpeg enhanced size 512,360 font NimbusRoman 10

set output "pca.jpg"

set title "Principal component analysis of 2D biased random samples"

load "pcalabels.gp"

plot  [0:120][0:80] "pca.data" using 1:2 with points pt 7 ps .75 t 'Samples', \
   "pcae1.data" using 1:2 with lines lc "red" t "λ0 axis", \
   "pcae2.data" using 1:2 with lines lc "green" t "λ1 axis"

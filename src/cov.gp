
set terminal jpeg enhanced size 512,480 font NimbusRoman 10

set output "cov.jpg"

set title "cov"

plot  "cov.data" using 1:2 with points pt 7 ps .75 not, \
      "evec.data" using 1:2 with lines not

